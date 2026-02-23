TEST 3: Spatial Transcriptomics Analysis (SYNTHETIC_DATA Mode)
==============================================================

> **Data Mode:** This test uses **SYNTHETIC_DATA** — `*_DRY_RUN=false`. The spatialtools server parses the actual CSV files in `/data/patient-data/PAT001-OVC-2025/spatial/`. No heavy bioinformatics tools required, but Python parsing dependencies (pandas, numpy, scipy) must be installed. See [Data Modes Guide](../../data-modes-guide.md) for details.

Patient ID: PAT001-OVC-2025

**This test uses TABULAR DATA (CSV files), NOT images.**
Use ONLY the spatialtools server. DeepCell is NOT needed.

## Prerequisites

| Requirement | Details |
|------------|---------|
| `SPATIALTOOLS_DRY_RUN` | `false` |
| Python deps | pandas, numpy, scipy (for CSV parsing, spatial stats) |
| Data files | `data/patient-data/PAT001-OVC-2025/spatial/` |

## Data Files (Real CSVs)

Files parsed by the server from `data/patient-data/PAT001-OVC-2025/spatial/`:
- `visium_spatial_coordinates.csv` — Spot x,y coordinates
- `visium_gene_expression.csv` — Gene expression per spot
- `visium_region_annotations.csv` — Region labels per spot

## Analysis Steps

### Pre-Flight Validation

1. **File existence and integrity check:**
   - Confirm all 3 CSV files are accessible and parseable
   - Coordinates file: row count, column names
   - Expression file: row count, gene columns
   - Regions file: row count, unique region labels

2. **Data consistency check:**
   - Same number of spots across all files
   - Expected genes present in expression file
   - Expected regions present in annotations

### Spatial Analysis

1. **Load spatial data:**
   ```
   Use mcp-spatialtools to load from:
   - data/patient-data/PAT001-OVC-2025/spatial/visium_spatial_coordinates.csv
   - data/patient-data/PAT001-OVC-2025/spatial/visium_gene_expression.csv
   - data/patient-data/PAT001-OVC-2025/spatial/visium_region_annotations.csv
   ```
   - How many spatial spots?
   - What spatial regions are identified?
   - How many spots per region?

2. **Focus on KEY GENES (8 genes):**
   - **Proliferation:** MKI67, PCNA
   - **Resistance:** PIK3CA, AKT1, ABCB1
   - **Immune:** CD3D, CD8A, CD68

3. **Gene expression by region:**
   - Which regions have high proliferation?
   - Where are resistance markers expressed?
   - Where are immune cells located?

4. **Spatial patterns:**
   - Spatial heterogeneity in resistance markers?
   - Immune cell exclusion from tumor regions?
   - Tumor-stroma interface patterns?

5. **Generate visualizations:**
   - Spatial heatmap, region composition, gene expression heatmap, spatial autocorrelation

## How This Differs from DRY_RUN

| Aspect | DRY_RUN | SYNTHETIC_DATA (this test) |
|--------|---------|---------------------------|
| Data source | Hardcoded inline values | Parsed from actual CSVs on disk |
| Spot counts | Returns fixed 900 | Counts actual rows in CSV |
| Region labels | Returns fixed 6 regions | Reads unique values from annotations CSV |
| Expression values | Returns predetermined means | Calculates from actual expression matrix |
| Moran's I | Returns fixed statistics | Computes from real spatial coordinates + expression |
| File I/O | None | Reads 3 CSV files from `data/` directory |

## Expected Results

Results come from parsing the actual synthetic CSV files:

**Spatial Structure:**
- Total spots: parsed from CSV (designed as ~900)
- Regions: parsed from annotations (designed as 6 regions)
- Spot counts per region: calculated from actual data

**Expression Patterns:**
- Gene expression means calculated from real CSV values per region
- Spatial autocorrelation computed from real coordinates
- Patterns should match the synthetic design (proliferation in tumor, immune in stroma)

**Spatial Findings:**
- Proliferation pattern from real MKI67/PCNA values
- Resistance marker distribution from real PIK3CA/AKT1/ABCB1 values
- Immune localization from real CD3D/CD8A/CD68 values

## Output Format

Please provide:

1. **Spatial Structure** — Spot counts and regions from parsed CSV
2. **Gene Expression Heatmap** — 8 genes x regions from actual data
3. **Key Spatial Findings** — Proliferation, resistance, immune patterns
4. **Visualizations** — Spatial heatmaps, region charts, autocorrelation
5. **Clinical Interpretation** — Heterogeneity, immune exclusion, TME classification

## Validation Checkpoints

- [ ] Server correctly reads CSVs from `data/patient-data/PAT001-OVC-2025/spatial/`
- [ ] Spot count matches actual CSV row count (not hardcoded 900)
- [ ] Region labels read from annotations (not hardcoded list)
- [ ] Expression values are real numbers from CSV (not predetermined)
- [ ] Moran's I calculated from actual spatial data
- [ ] Biological patterns (immune exclusion, tumor proliferation) preserved
