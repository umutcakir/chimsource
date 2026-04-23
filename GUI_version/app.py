"""
ChiMSource Web GUI — Flask application with real-time progress via SSE.

Launch with:  python app.py
Then open:    http://localhost:5000
"""

import os
import sys
import json
import time
import uuid
import re
import threading
import queue
from flask import (Flask, render_template, request, jsonify,
                   Response, send_file, stream_with_context)
import pandas as pd
import multiprocessing

from engine import run_analysis, load_codon_table

# ---------------------------------------------------------------------------
# App setup
# ---------------------------------------------------------------------------
app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 2 * 1024 * 1024 * 1024  # 2 GB upload limit

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
UPLOAD_DIR = os.path.join(BASE_DIR, 'uploads')
RESULTS_DIR = os.path.join(BASE_DIR, 'results')
os.makedirs(UPLOAD_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)

# In-memory store for running jobs  {job_id: {...}}
jobs = {}

# ---------------------------------------------------------------------------
# Routes — pages
# ---------------------------------------------------------------------------

@app.route('/')
def index():
    cpu_count = multiprocessing.cpu_count()
    return render_template('index.html', cpu_count=cpu_count)

# ---------------------------------------------------------------------------
# Routes — API
# ---------------------------------------------------------------------------

@app.route('/api/submit', methods=['POST'])
def submit():
    """Accept files + parameters, start background analysis, return job_id."""
    try:
        nuc_file = request.files.get('nucleotide_file')
        pep_file = request.files.get('peptide_file')
        if not nuc_file or not pep_file:
            return jsonify(error="Both nucleotide and peptide files are required."), 400

        job_id = uuid.uuid4().hex[:12]
        job_dir = os.path.join(RESULTS_DIR, job_id)
        os.makedirs(job_dir, exist_ok=True)

        nuc_path = os.path.join(job_dir, nuc_file.filename)
        pep_path = os.path.join(job_dir, pep_file.filename)
        nuc_file.save(nuc_path)
        pep_file.save(pep_path)

        # Collect parameters from the form
        params = {
            'nucleotide_file': nuc_path,
            'peptide_file': pep_path,
            'codon_table_id': int(request.form.get('codon_table', 1)),
            'max_gap': int(request.form.get('max_gap', 2)),
            'num_threads': int(request.form.get('num_threads', 1)),
            'max_transcript_length': int(request.form.get('max_transcript_length', 30000)),
            'max_flanking_seq': int(request.form.get('max_flanking_seq', 500)),
            'no_reverse_complement_check': request.form.get('no_reverse_complement_check') == 'true',
            'verbose': request.form.get('verbose') == 'true',
            'save_log': request.form.get('save_log') != 'false',  # default ON
            'output_base': os.path.join(job_dir,
                           re.sub(r'[^\w\-.]', '_', request.form.get('output_name', 'results')) or 'results'),
        }

        # Prepare event queue for SSE
        q = queue.Queue()
        jobs[job_id] = {
            'queue': q,
            'status': 'running',
            'params': params,
            'result_file': None,
            'error': None,
        }

        # Launch background thread
        t = threading.Thread(target=_run_job, args=(job_id,), daemon=True)
        t.start()

        return jsonify(job_id=job_id)

    except Exception as e:
        return jsonify(error=str(e)), 500


@app.route('/api/progress/<job_id>')
def progress(job_id):
    """SSE endpoint — streams log and progress events."""
    job = jobs.get(job_id)
    if not job:
        return jsonify(error="Job not found"), 404

    def generate():
        q = job['queue']
        while True:
            try:
                msg = q.get(timeout=60)
            except queue.Empty:
                # keep-alive
                yield "data: {}\n\n"
                continue
            yield f"data: {json.dumps(msg)}\n\n"
            if msg.get('type') in ('done', 'error'):
                break
    return Response(stream_with_context(generate()),
                    mimetype='text/event-stream',
                    headers={'Cache-Control': 'no-cache',
                             'X-Accel-Buffering': 'no'})


@app.route('/api/results/<job_id>')
def get_results(job_id):
    """Return a page of results as JSON (paginated)."""
    job = jobs.get(job_id)
    if not job or job['status'] != 'done':
        return jsonify(error="Results not ready"), 404
    tsv_path = job['result_file']
    if not tsv_path or not os.path.exists(tsv_path):
        return jsonify(error="No results file"), 404
    try:
        df = pd.read_csv(tsv_path, sep='\t')
        page = int(request.args.get('page', 1))
        per_page = int(request.args.get('per_page', 50))
        total_rows = len(df)
        total_pages = max(1, (total_rows + per_page - 1) // per_page)
        page = max(1, min(page, total_pages))
        start = (page - 1) * per_page
        end = min(start + per_page, total_rows)
        page_df = df.iloc[start:end]
        return jsonify(
            total_rows=total_rows,
            total_pages=total_pages,
            page=page,
            per_page=per_page,
            columns=list(df.columns),
            rows=page_df.values.tolist(),
        )
    except Exception as e:
        return jsonify(error=str(e)), 500


@app.route('/api/download/<job_id>/<filetype>')
def download(job_id, filetype):
    """Download result TSV or parameters TXT."""
    job = jobs.get(job_id)
    if not job:
        return jsonify(error="Job not found"), 404
    output_base = job['params']['output_base']
    if filetype == 'tsv':
        path = output_base + '.tsv'
    elif filetype == 'params':
        path = output_base + '_parameters.txt'
    elif filetype == 'log':
        path = output_base + '_log.txt'
    else:
        return jsonify(error="Unknown file type"), 400
    if not os.path.exists(path):
        return jsonify(error="File not found"), 404
    return send_file(path, as_attachment=True)

# ---------------------------------------------------------------------------
# Background worker
# ---------------------------------------------------------------------------

def _run_job(job_id):
    job = jobs[job_id]
    q = job['queue']
    params = job['params']
    save_log = params.pop('save_log', True)
    log_path = params['output_base'] + '_log.txt'
    log_file = None

    if save_log:
        try:
            log_file = open(log_path, 'w', encoding='utf-8')
        except Exception:
            log_file = None

    from datetime import datetime

    def log_cb(msg):
        q.put({'type': 'log', 'message': msg})
        if log_file:
            try:
                log_file.write(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}\n")
                log_file.flush()
            except Exception:
                pass

    def progress_cb(stage, done, total):
        q.put({'type': 'progress', 'stage': stage,
               'done': done, 'total': total})

    try:
        combined, params_text, output_tsv = run_analysis(
            **params,
            log_callback=log_cb,
            progress_callback=progress_cb,
        )
        job['status'] = 'done'
        job['result_file'] = output_tsv
        n_results = len(combined) if combined is not None else 0
        job['has_log'] = save_log and log_file is not None
        q.put({'type': 'done', 'total_results': n_results, 'has_log': job.get('has_log', False)})
    except Exception as e:
        job['status'] = 'error'
        job['error'] = str(e)
        if log_file:
            log_file.write(f"[{datetime.now().strftime('%H:%M:%S')}] ERROR: {str(e)}\n")
        q.put({'type': 'error', 'message': str(e)})
    finally:
        if log_file:
            try:
                log_file.close()
            except Exception:
                pass

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    multiprocessing.freeze_support()  # Required for Windows multiprocessing
    # Allow setting port via env or default to 5000
    port = int(os.environ.get('CHIMSOURCE_PORT', 5000))
    print(f"\n{'='*60}")
    print(f"  ChiMSource Web GUI")
    print(f"  Open your browser at:  http://localhost:{port}")
    print(f"{'='*60}\n")
    app.run(host='0.0.0.0', port=port, debug=False, threaded=True)
