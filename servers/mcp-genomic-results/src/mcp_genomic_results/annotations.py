"""Clinical annotation databases for somatic variants and copy number alterations.

Hardcoded lookup tables for ClinVar, COSMIC, and therapy associations relevant
to high-grade serous ovarian carcinoma (HGSOC).  These are intentionally
simplified for demonstration and are NOT a substitute for real clinical
annotation pipelines (e.g., OncoKB, CIViC, ACMG).
"""

from typing import Dict, Optional, Any

# ---------------------------------------------------------------------------
# Somatic variant annotations (keyed by VCF ID like "TP53_R175H")
# ---------------------------------------------------------------------------

VARIANT_ANNOTATIONS: Dict[str, Dict[str, Any]] = {
    "TP53_R175H": {
        "gene": "TP53",
        "protein_change": "p.R175H",
        "clinvar_id": "VCV000012347",
        "cosmic_id": "COSV57271936",
        "classification": "Pathogenic",
        "oncogenicity": "Oncogenic",
        "frequency_in_hgsoc": 0.96,
        "therapies": [
            {"drug": "APR-246 (eprenetapopt)", "evidence": "Phase III", "status": "Investigational"},
        ],
        "notes": "Most common TP53 hotspot in HGSOC. Loss of tumor suppressor function.",
    },
    "PIK3CA_E545K": {
        "gene": "PIK3CA",
        "protein_change": "p.E545K",
        "clinvar_id": "VCV000013652",
        "cosmic_id": "COSV55514880",
        "classification": "Pathogenic",
        "oncogenicity": "Oncogenic",
        "frequency_in_hgsoc": 0.04,
        "therapies": [
            {"drug": "Alpelisib (BYL719)", "evidence": "FDA-approved (breast)", "status": "Off-label"},
            {"drug": "Copanlisib", "evidence": "Phase II", "status": "Investigational"},
        ],
        "notes": "PI3K pathway activating mutation. Uncommon in HGSOC but actionable.",
    },
    "PTEN_LOH": {
        "gene": "PTEN",
        "protein_change": "LOH (splice acceptor)",
        "clinvar_id": "VCV000233105",
        "cosmic_id": None,
        "classification": "Pathogenic",
        "oncogenicity": "Oncogenic",
        "frequency_in_hgsoc": 0.07,
        "therapies": [
            {"drug": "AKT inhibitors (capivasertib)", "evidence": "Phase III", "status": "Investigational"},
        ],
        "notes": "PTEN loss activates PI3K/AKT pathway. Synergistic with PIK3CA mutation.",
    },
}

# ---------------------------------------------------------------------------
# Copy-number annotations (keyed by gene symbol)
# ---------------------------------------------------------------------------

CN_ANNOTATIONS: Dict[str, Dict[str, Any]] = {
    "MYC": {
        "cytoband": "8q24.21",
        "alteration": "Amplification",
        "clinical_significance": "Associated with aggressive disease and poor prognosis",
        "frequency_in_hgsoc": 0.30,
        "therapies": [
            {"drug": "BET inhibitors (JQ1)", "evidence": "Preclinical", "status": "Investigational"},
        ],
    },
    "CCNE1": {
        "cytoband": "19q12",
        "alteration": "Amplification",
        "clinical_significance": "Platinum resistance marker; excludes HRD benefit",
        "frequency_in_hgsoc": 0.20,
        "therapies": [
            {"drug": "CDK2 inhibitors", "evidence": "Phase I/II", "status": "Investigational"},
            {"drug": "Wee1 inhibitors (adavosertib)", "evidence": "Phase II", "status": "Investigational"},
        ],
    },
    "AKT2": {
        "cytoband": "19q13.2",
        "alteration": "Amplification",
        "clinical_significance": "PI3K/AKT pathway activation; synergistic with PTEN loss",
        "frequency_in_hgsoc": 0.12,
        "therapies": [
            {"drug": "Capivasertib (AKT inhibitor)", "evidence": "Phase III", "status": "Investigational"},
        ],
    },
    "RB1": {
        "cytoband": "13q14.2",
        "alteration": "Deletion",
        "clinical_significance": "Cell cycle deregulation; may predict CDK4/6 inhibitor resistance",
        "frequency_in_hgsoc": 0.08,
        "therapies": [],
    },
    "CDKN2A": {
        "cytoband": "9p21.3",
        "alteration": "Deletion",
        "clinical_significance": "Loss of p16 cell cycle checkpoint; aggressive phenotype",
        "frequency_in_hgsoc": 0.05,
        "therapies": [
            {"drug": "CDK4/6 inhibitors (palbociclib)", "evidence": "Phase II", "status": "Off-label"},
        ],
    },
}

# ---------------------------------------------------------------------------
# Gene lists
# ---------------------------------------------------------------------------

HRD_GENES = [
    "BRCA1", "BRCA2", "PALB2", "RAD51C", "RAD51D", "ATM", "ATR",
    "CHEK1", "CHEK2", "BARD1", "BRIP1", "FANCA", "FANCD2",
]

OVC_GENE_PANEL = [
    "TP53", "BRCA1", "BRCA2", "PIK3CA", "PTEN", "KRAS", "BRAF",
    "ARID1A", "MYC", "CCNE1", "AKT2", "RB1", "CDKN2A", "NF1",
    "CDK12", "EMSY", "RAD51C", "RAD51D",
]


# ---------------------------------------------------------------------------
# Lookup functions
# ---------------------------------------------------------------------------

def get_variant_annotation(variant_id: str) -> Optional[Dict[str, Any]]:
    """Look up clinical annotation for a somatic variant by its VCF ID.

    Args:
        variant_id: VCF ID field value (e.g., "TP53_R175H").

    Returns:
        Annotation dict or None if not in the database.
    """
    return VARIANT_ANNOTATIONS.get(variant_id)


def get_cn_annotation(gene: str) -> Optional[Dict[str, Any]]:
    """Look up clinical annotation for a copy-number altered gene.

    Args:
        gene: Gene symbol (e.g., "MYC").

    Returns:
        Annotation dict or None if not in the database.
    """
    return CN_ANNOTATIONS.get(gene)
