TEST 1: Clinical Data and Genomic Analysis (SYNTHETIC_DATA Mode)
=================================================================

> **Data Mode:** This test uses **SYNTHETIC_DATA** — `*_DRY_RUN=false`. Servers parse the actual generated files in `/data/patient-data/PAT001-OVC-2025/`. No heavy bioinformatics tools required, but Python parsing dependencies must be installed. See [Data Modes Guide](../../data-modes-guide.md) for details.

Patient ID: PAT001-OVC-2025

## Prerequisites

| Requirement | Details |
|------------|---------|
| `MOCKEPIC_DRY_RUN` | `false` |
| `FGBIO_DRY_RUN` | `false` |
| `GENOMIC_RESULTS_DRY_RUN` | `false` |
| `MOCKTCGA_DRY_RUN` | `false` |
| Python deps | Standard libraries (json, csv parsing) |
| Data files | `data/patient-data/PAT001-OVC-2025/clinical/` and `genomics/` |

## PART 1: Clinical Data (use mcp-mockepic)

The patient data files are located at `data/patient-data/PAT001-OVC-2025/clinical/`.

1. For patient PAT001-OVC-2025, retrieve:
   - Patient demographics (name, age, family history)
   - Genetic mutations noted in family history
   - BRCA1 status

2. Retrieve lab results for this patient:
   - CA-125 tumor marker trends over time
   - What does the CA-125 trajectory indicate about treatment response?
   - Evidence of platinum resistance in the lab data?

Files parsed by the server:
- `data/patient-data/PAT001-OVC-2025/clinical/patient_demographics.json`
- `data/patient-data/PAT001-OVC-2025/clinical/lab_results.json`

## PART 2: Genomic Analysis (use mcp-genomic-results and mcp-mocktcga)

3. Parse somatic variants for patient PAT001-OVC-2025:
   Use genomic-results to read the VCF file:
   - File: `data/patient-data/PAT001-OVC-2025/genomics/somatic_variants.vcf`
   - Look for TP53, PIK3CA, PTEN mutations

   Parse copy number alterations:
   - File: `data/patient-data/PAT001-OVC-2025/genomics/copy_number_results.cns`
   - Look for MYC, CCNE1, AKT2 amplifications and RB1, CDKN2A deletions

4. Compare to TCGA-OV cohort (use mcp-mocktcga):
   For a patient with:
   - BRCA1 germline mutation
   - TP53 somatic mutation
   - Stage IV HGSOC
   - Platinum-resistant progression

   Questions:
   - What TCGA molecular subtype does this match?
   - What is the typical prognosis for this profile?
   - What pathways are commonly activated in platinum-resistant ovarian cancer?

## How This Differs from DRY_RUN

| Aspect | DRY_RUN | SYNTHETIC_DATA (this test) |
|--------|---------|---------------------------|
| Data source | Hardcoded inline mock values | Parsed from actual files on disk |
| VCF parsing | Returns fixed variant list | Reads and parses `somatic_variants.vcf` |
| CNS parsing | Returns fixed CNV list | Reads and parses `copy_number_results.cns` |
| JSON parsing | Returns fixed demographics | Reads and parses `patient_demographics.json` |
| File I/O | None | Real file reads from `data/` directory |

## Expected Results

Results should match the synthetic data files. Validate that the server correctly parsed:

**Clinical (from JSON files):**
- Patient name, age, diagnosis from `patient_demographics.json`
- CA-125 values and trends from `lab_results.json`
- BRCA1 status from clinical records

**Genomic (from VCF/CNS files):**
- Somatic variants parsed from `somatic_variants.vcf`
- Copy number calls parsed from `copy_number_results.cns`
- Variants should include chromosomal positions, ref/alt alleles, and annotations

**TCGA:**
- Subtype matching and prognosis (mocktcga uses same logic in both modes)

## Output Format

Please provide:
1. Patient summary (demographics, genetic risk) — from parsed JSON
2. CA-125 trend analysis — from parsed lab results
3. Key somatic mutations identified — from parsed VCF
4. Copy number alterations — from parsed CNS
5. TCGA subtype and pathway analysis
6. Note any discrepancies between parsed values and expected synthetic data

## Validation Checkpoints

- [ ] Server correctly reads files from `data/patient-data/PAT001-OVC-2025/`
- [ ] VCF parsing returns structured variant records (not hardcoded mock)
- [ ] CNS parsing returns copy number segments (not hardcoded mock)
- [ ] JSON parsing returns patient demographics (not hardcoded mock)
- [ ] Results are consistent with the synthetic data design
