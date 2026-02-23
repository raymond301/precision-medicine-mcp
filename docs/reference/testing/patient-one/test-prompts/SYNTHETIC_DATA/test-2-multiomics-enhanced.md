TEST 2: Multi-Omics Resistance Analysis (SYNTHETIC_DATA Mode)
==============================================================

> **Data Mode:** This test uses **SYNTHETIC_DATA** — `*_DRY_RUN=false`. The multiomics server parses the actual CSV files in `/data/patient-data/PAT001-OVC-2025/multiomics/`. No heavy bioinformatics tools required, but Python parsing dependencies (pandas, numpy, scipy) must be installed. See [Data Modes Guide](../../data-modes-guide.md) for details.

Patient ID: PAT001-OVC-2025

## Prerequisites

| Requirement | Details |
|------------|---------|
| `MULTIOMICS_DRY_RUN` | `false` |
| Python deps | pandas, numpy, scipy (for CSV parsing, stats) |
| Data files | `data/patient-data/PAT001-OVC-2025/multiomics/` |

## Data Files (Real CSVs)

Files parsed by the server from `data/patient-data/PAT001-OVC-2025/multiomics/`:
- `sample_metadata.csv` — Sample annotations (includes Batch column)
- `pdx_rna_seq.csv` — RNA-seq expression values
- `pdx_proteomics.csv` — TMT proteomics data
- `pdx_phosphoproteomics.csv` — Phosphoproteomics data

## Analysis Steps

### STEP 0: DATA PREPROCESSING

1. **Validate Data Quality:**
   ```
   Use mcp-multiomics tool: validate_multiomics_data

   Inputs:
   - rna_path: data/patient-data/PAT001-OVC-2025/multiomics/pdx_rna_seq.csv
   - protein_path: data/patient-data/PAT001-OVC-2025/multiomics/pdx_proteomics.csv
   - phospho_path: data/patient-data/PAT001-OVC-2025/multiomics/pdx_phosphoproteomics.csv
   - metadata_path: data/patient-data/PAT001-OVC-2025/multiomics/sample_metadata.csv
   ```

   Check for batch effects, missing value patterns, sample consistency, and outliers.

2. **Preprocess Data:**
   ```
   Use mcp-multiomics tool: preprocess_multiomics_data

   Apply: Batch correction (ComBat), KNN imputation, quantile normalization, outlier removal
   ```

3. **Visualize QC:**
   ```
   Use mcp-multiomics tool: visualize_data_quality

   Generate before/after PCA plots, correlation heatmaps, missing value patterns
   ```

### STEP 1: DATA INTEGRATION

4. **Integrate Preprocessed Data:**
   ```
   Use mcp-multiomics tool: integrate_omics_data

   Load preprocessed files from all 3 modalities
   ```

### STEP 2: ASSOCIATION TESTING & META-ANALYSIS

5. **Focus on KEY RESISTANCE GENES:**
   - **PI3K/AKT pathway:** PIK3CA, AKT1, MTOR, PTEN
   - **Drug resistance:** ABCB1 (MDR1)
   - **Anti-apoptotic:** BCL2L1
   - **Tumor suppressor:** TP53

6. **Run Stouffer's Meta-Analysis:**
   ```
   Use mcp-multiomics tool: calculate_stouffer_meta

   For each gene: combine p-values from RNA, Protein, Phospho
   Method: Stouffer's Z-score with FDR correction after combination
   ```

### STEP 3: UPSTREAM REGULATOR PREDICTION

7. **Predict Therapeutic Targets:**
   ```
   Use mcp-multiomics tool: predict_upstream_regulators

   Input: Significant genes from Stouffer's (q < 0.05)
   Analyze for kinases, transcription factors, drug responses
   ```

## How This Differs from DRY_RUN

| Aspect | DRY_RUN | SYNTHETIC_DATA (this test) |
|--------|---------|---------------------------|
| Data source | Hardcoded inline values | Parsed from actual CSVs on disk |
| CSV parsing | Returns fixed sample counts | Reads real CSV structure and values |
| Batch effects | Returns simulated correlation | Calculates actual PCA from data |
| Missing values | Returns fixed percentage | Counts actual NaN values in CSVs |
| Stouffer's | Returns predetermined Z-scores | Calculates from actual expression values |
| File I/O | None | Reads 4 CSV files from `data/` directory |

## Expected Results

Results come from parsing the actual synthetic CSV files. Values may differ slightly from DRY_RUN hardcoded values but should show the same biological patterns:

**Sample Summary:**
- Sample counts from `sample_metadata.csv`
- Resistant and sensitive group sizes from metadata

**Gene-Level Results:**
- Expression fold changes calculated from real CSV values
- Stouffer's Z-scores computed from actual p-values
- FDR-corrected q-values

**Pathway Analysis:**
- PI3K/AKT/mTOR pathway activation pattern should be preserved
- Resistance genes (PIK3CA, AKT1, ABCB1) should trend upward
- Tumor suppressors (PTEN, TP53) should trend downward

## Output Format

Please provide:

1. **Preprocessing Summary** — Validation results from actual CSV parsing
2. **Sample Summary** — Counts from real metadata file
3. **Gene-Level Results** — Stouffer's Z-scores computed from actual data
4. **Upstream Regulator Predictions** — Based on computed significant genes
5. **Pathway Interpretation** — PI3K/AKT pathway activation assessment

## Validation Checkpoints

- [ ] Server correctly reads CSVs from `data/patient-data/PAT001-OVC-2025/multiomics/`
- [ ] Sample metadata parsed with correct group assignments
- [ ] Expression values are real numbers (not hardcoded mock)
- [ ] Stouffer's Z-scores calculated from actual per-modality p-values
- [ ] Biological patterns (resistance gene upregulation) preserved
- [ ] FDR correction applied after Stouffer's combination
