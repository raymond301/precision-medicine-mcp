# PatientOne Clinical Profile (Canonical Reference)

**Patient ID:** PAT001-OVC-2025
**Name:** Sarah Anderson (synthetic)
**Age:** 58 years
**Status:** 100% synthetic data — safe to share, no IRB required

---

## Diagnosis

- **Cancer Type:** Stage IV High-Grade Serous Ovarian Carcinoma (HGSOC)
- **Treatment Status:** Platinum-resistant
- **BRCA1:** Germline mutation (c.5266dupC, pathogenic)

## Key Genomic Findings

Source: [`data/patient-data/PAT001-OVC-2025/genomics/somatic_variants.vcf`](../../../data/patient-data/PAT001-OVC-2025/genomics/somatic_variants.vcf)

| Gene | Variant | Type | VAF | Clinical Significance |
|------|---------|------|-----|----------------------|
| TP53 | R175H | Missense | 73% | Pathogenic — prognostic |
| PIK3CA | E545K | Missense | 42% | Pathogenic — PI3K inhibitor responsive |
| PTEN | LOH | Splice acceptor | 85% | Pathogenic — PI3K/AKT pathway activation |
| BRCA1 | Germline | c.5266dupC | — | Pathogenic — PARP inhibitor responsive |

**Copy Number Variants:** MYC amplification, CCNE1 amplification, AKT2 amplification, RB1 deletion, CDKN2A deletion

**Tumor Mutational Burden (TMB):** 8.2 mutations/Mb (high)

## Lab Results

- **CA-125 trajectory:** 1456 → 22 → 389 → 289 U/mL (platinum-resistant pattern)

## Multi-Omics Data (15 PDX Samples)

- 7 resistant samples (high PIK3CA, AKT1, MTOR, ABCB1 activation)
- 8 sensitive samples (lower pathway activation)
- Modalities: RNA-seq, proteomics, phosphoproteomics

## Spatial Transcriptomics

- 900 spots x 31 genes (simplified Visium layout)
- 6 tissue regions: tumor core, tumor proliferative, tumor interface, stroma, stroma immune, necrotic/hypoxic
- Key markers: Ki67, CD3D, CD8A, HIF1A, COL1A1

## Treatment Recommendation

Based on integrated analysis (BRCA1 germline + PIK3CA E545K + spatial microenvironment):
- **Olaparib** (PARP inhibitor) — for BRCA1 germline mutation
- **Alpelisib** (PI3K inhibitor) — for PIK3CA E545K
- Checkpoint inhibitors — based on CD8+ T cell infiltration and TMB

## Data Files

**Total:** ~4.9 MB across 5 modalities (20 files)

| Modality | Files | Location |
|----------|-------|----------|
| Clinical | 2 JSON | `clinical/` |
| Genomic | 2 (VCF + CNS) | `genomics/` |
| Multi-omics | 4 CSV | `multiomics/` |
| Spatial | 4 CSV | `spatial/` |
| Imaging | 7 TIFF (placeholders) | `imaging/` |

**GCS path:** `gs://sample-inputs-patientone/patient-data/PAT001-OVC-2025/`
**Local path:** `data/patient-data/PAT001-OVC-2025/`

---

**Full dataset documentation:** [`data/patient-data/PAT001-OVC-2025/README.md`](../../../data/patient-data/PAT001-OVC-2025/README.md)
