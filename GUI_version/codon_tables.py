"""
All 27 NCBI genetic code tables, bundled directly.

Source: INSDC / NCBI Taxonomy
https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
https://www.insdc.org/submitting-standards/genetic-code-tables/

This module eliminates the need for biotite or biopython just to load
codon tables.  The tables are static data maintained by NCBI; IDs 7 and
8 were retired and are intentionally absent.
"""

NCBI_CODON_TABLE_NAMES = {
    1: "The Standard Code",
    2: "Vertebrate Mitochondrial",
    3: "Yeast Mitochondrial",
    4: "Mold, Protozoan, Coelenterate Mito.; Mycoplasma; Spiroplasma",
    5: "Invertebrate Mitochondrial",
    6: "Ciliate, Dasycladacean and Hexamita Nuclear",
    9: "Echinoderm and Flatworm Mitochondrial",
    10: "Euplotid Nuclear",
    11: "Bacterial, Archaeal and Plant Plastid",
    12: "Alternative Yeast Nuclear",
    13: "Ascidian Mitochondrial",
    14: "Alternative Flatworm Mitochondrial",
    15: "Blepharisma Nuclear",
    16: "Chlorophycean Mitochondrial",
    21: "Trematode Mitochondrial",
    22: "Scenedesmus obliquus Mitochondrial",
    23: "Thraustochytrium Mitochondrial",
    24: "Rhabdopleuridae Mitochondrial",
    25: "Candidate Division SR1 and Gracilibacteria",
    26: "Pachysolen tannophilus Nuclear",
    27: "Karyorelict Nuclear",
    28: "Condylostoma Nuclear",
    29: "Mesodinium Nuclear",
    30: "Peritrich Nuclear",
    31: "Blastocrithidia Nuclear",
    32: "Balanophoraceae Plastid",
    33: "Cephalodiscidae Mitochondrial UAA-Tyr",
}

# ---------------------------------------------------------------------------
# Helper: build a codon dict from the NCBI amino-acid string
# ---------------------------------------------------------------------------
_BASES1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
_BASES2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
_BASES3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
_CODONS = [_BASES1[i] + _BASES2[i] + _BASES3[i] for i in range(64)]

def _build(aa_string):
    return {_CODONS[i]: aa_string[i] for i in range(64)}

# ---------------------------------------------------------------------------
# All 27 tables, keyed by NCBI transl_table ID
# ---------------------------------------------------------------------------
NCBI_CODON_TABLES = {
    1:  _build("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
    2:  _build("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"),
    3:  _build("FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
    4:  _build("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
    5:  _build("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG"),
    6:  _build("FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
    9:  _build("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"),
    10: _build("FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
    11: _build("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
    12: _build("FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
    13: _build("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG"),
    14: _build("FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"),
    15: _build("FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
    16: _build("FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
    21: _build("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG"),
    22: _build("FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
    23: _build("FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
    24: _build("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG"),
    25: _build("FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
    26: _build("FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
    27: _build("FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
    28: _build("FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
    29: _build("FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
    30: _build("FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
    31: _build("FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
    32: _build("FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG"),
    33: _build("FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG"),
}


def load_codon_table(table_id):
    """Load a codon table by NCBI ID.  Returns (dict, name_string)."""
    if table_id == 0:
        table_id = 1  # 0 = default = standard
    if table_id not in NCBI_CODON_TABLES:
        available = sorted(NCBI_CODON_TABLES.keys())
        raise ValueError(
            f"Unknown codon table #{table_id}. "
            f"Available IDs: {available}"
        )
    return NCBI_CODON_TABLES[table_id], f"{NCBI_CODON_TABLE_NAMES[table_id]} (#{table_id})"
