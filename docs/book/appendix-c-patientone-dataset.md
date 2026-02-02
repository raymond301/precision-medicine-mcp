# Appendix C: PatientOne Complete Dataset

*Reference guide for the complete PAT001-OVC-2025 synthetic dataset*

---

## Dataset Overview

**Patient ID**: PAT001-OVC-2025
**Diagnosis**: Stage IV high-grade serous ovarian carcinoma
**Age**: 58 years
**Status**: 100% synthetic (safe to share, no IRB required)
**Purpose**: Testing, training, and demonstration of AI-orchestrated precision medicine workflows

**Location**: [`data/patient-data/PAT001-OVC-2025/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/data/patient-data/PAT001-OVC-2025)

**Public access**: `gs://precision-medicine-mcp-public/patient-data/PAT001-OVC-2025/`

---

## Complete File Manifest

### Clinical Data (2 files)

**Location**: `clinical/`

| File | Format | Size | Description |
|------|--------|------|-------------|
| `patient-001.json` | FHIR R4 | 4 KB | Patient demographics, conditions, medications |
| `genomic-test-order-001.json` | FHIR R4 | 2 KB | Genomic test orders, ServiceRequest resources |

**Key fields**:
- Patient: Sarah Anderson, 58yo, Female
- Diagnosis: Ovarian cancer (C56.9), Stage IV
- BRCA status: BRCA1+
- Medications: Carboplatin, paclitaxel, bevacizumab

---

### Genomics Data (1 file)

**Location**: `genomics/`

| File | Format | Size | Description |
|------|--------|------|-------------|
| `somatic-variants.vcf` | VCF 4.2 | 12 KB | 8 pathogenic somatic mutations |

**Pathogenic variants**:

| Gene | Variant | Type | VAF | ClinVar | Actionability |
|------|---------|------|-----|---------|---------------|
| TP53 | R175H | Missense | 85% | Pathogenic (★★★★★) | Prognostic |
| PIK3CA | E545K | Missense | 42% | Pathogenic (★★★★) | PI3K inhibitor (alpelisib) |
| KRAS | G12D | Missense | 38% | Pathogenic (★★★★★) | Prognostic (poor) |
| PTEN | R130* | Nonsense | 55% | Pathogenic (★★★★★) | PI3K/AKT pathway activation |
| ARID1A | Q456* | Nonsense | 48% | Pathogenic (★★★★) | SWI/SNF complex loss |
| CTNNB1 | S45F | Missense | 31% | Pathogenic (★★★★) | Wnt pathway activation |
| FBXW7 | R465C | Missense | 27% | Likely pathogenic (★★★) | Tumor suppressor loss |
| PPP2R1A | P179R | Missense | 33% | Pathogenic (★★★★) | Platinum resistance |

**Tumor mutational burden (TMB)**: 8.2 mutations/Mb (high)

---

### Multi-Omics Data (4 files)

**Location**: `multiomics/`

| File | Format | Size | Samples | Features | Description |
|------|--------|------|---------|----------|-------------|
| `rna_counts.tsv` | TSV | 280 KB | 15 PDX models | 2,000 genes | RNA-seq expression (TPM) |
| `protein_abundance.tsv` | TSV | 95 KB | 15 PDX models | 350 proteins | RPPA protein levels |
| `phospho_abundance.tsv` | TSV | 110 KB | 15 PDX models | 450 phospho-sites | Phosphoproteomics |
| `sample_metadata.tsv` | TSV | 3 KB | 15 PDX models | 8 annotations | Treatment response, histology |

**PDX model cohort**:
- **n = 15** ovarian cancer patient-derived xenografts
- **Treatment groups**: Carboplatin (n=5), Olaparib (n=5), Combination (n=5)
- **Response**: Complete response (n=4), Partial response (n=7), Progressive disease (n=4)

**Key genes/proteins**:
- **TP53**: Mutant expression confirmed (RNA, protein)
- **PIK3CA**: Elevated protein and phospho-AKT (S473)
- **MKI67**: 25% proliferative cells (protein)
- **CD8A/CD8B**: Moderate T cell infiltration

---

### Spatial Transcriptomics Data (3 files)

**Location**: `spatial/`

| File | Format | Size | Spots | Genes | Description |
|------|--------|------|-------|-------|-------------|
| `filtered_feature_bc_matrix.h5` | H5AD | 45 MB | 900 | 7,000 | 10X Visium gene expression |
| `spatial_coords.tsv` | TSV | 25 KB | 900 | 2 | Spot x,y coordinates |
| `tissue_positions.csv` | CSV | 30 KB | 900 | 5 | Tissue coordinates, selection |

**Spatial regions** (6 annotated):

| Region | Spots | % | Key Markers | Description |
|--------|-------|---|-------------|-------------|
| Tumor proliferative | 245 | 27% | MKI67+, PCNA+, TOP2A+ | Actively dividing tumor cells |
| Tumor quiescent | 198 | 22% | MKI67-, TP53+ | TP53-mutant non-dividing cells |
| Hypoxic niche | 112 | 12% | HIF1A+, CA9+, VEGFA+ | Oxygen-deprived regions |
| Immune infiltrate | 156 | 17% | CD8A+, CD8B+, GZMB+ | T cell rich areas |
| Stromal | 134 | 15% | COL1A1+, FAP+, ACTA2+ | Cancer-associated fibroblasts |
| Necrotic | 55 | 6% | Low UMI, low genes | Dead tissue |

**Moran's I results** (spatial autocorrelation):

| Gene | Moran's I | p-value | Interpretation |
|------|-----------|---------|----------------|
| MKI67 | 0.78 | < 0.001 | Highly clustered (proliferative zones) |
| CD8A | 0.65 | < 0.001 | Clustered (immune infiltrates) |
| HIF1A | 0.71 | < 0.001 | Clustered (hypoxic regions) |
| TP53 | 0.42 | < 0.001 | Moderately clustered |

---

### Imaging Data (7 files)

**Location**: `imaging/`

| File | Format | Size | Modality | Description |
|------|--------|------|----------|-------------|
| `he_section1.tif` | TIFF | 120 MB | H&E | Hematoxylin & eosin stain |
| `he_section2.tif` | TIFF | 118 MB | H&E | Second H&E section |
| `mxif_cd8_dapi.tif` | TIFF | 95 MB | MxIF | CD8 + DAPI (2-plex) |
| `mxif_ki67_tp53_dapi.tif` | TIFF | 98 MB | MxIF | Ki67 + TP53 + DAPI (3-plex) |
| `mxif_panck_vim_dapi.tif` | TIFF | 102 MB | MxIF | Pan-cytokeratin + Vimentin + DAPI |
| `nuclear_masks.png` | PNG | 12 MB | Segmentation | DeepCell nuclear masks |
| `cell_phenotypes.csv` | CSV | 85 KB | Annotation | Cell classifications |

**Imaging dimensions**:
- H&E: 2048 × 2048 pixels, 0.5 μm/pixel
- MxIF: 2048 × 2048 pixels, 0.5 μm/pixel
- Channels: DAPI (blue), FITC (green), TRITC (red), Cy5 (far-red)

**Cell counts**:
- Total cells segmented: 3,847
- Tumor cells (PanCK+): 2,234 (58%)
- Stromal cells (VIM+): 891 (23%)
- T cells (CD8+): 412 (11%)
- Proliferating cells (Ki67+): 961 (25%)
- TP53+ cells: 1,823 (47%)

---

## Data Access Methods

### Method 1: Clone Repository (Recommended for Local)

```bash
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp/data/patient-data/PAT001-OVC-2025
ls -lh
```

### Method 2: Download from GCS (Cloud Run)

```bash
# Public bucket (no authentication required)
gsutil -m cp -r gs://precision-medicine-mcp-public/patient-data/PAT001-OVC-2025/ ./data/
```

### Method 3: Use in Jupyter Notebooks

```python
import os
from pathlib import Path

# Set data directory
DATA_DIR = Path("../../../data/patient-data/PAT001-OVC-2025")

# Load clinical FHIR
import json
with open(DATA_DIR / "clinical/patient-001.json") as f:
    patient_fhir = json.load(f)

# Load genomic VCF
vcf_path = DATA_DIR / "genomics/somatic-variants.vcf"

# Load spatial h5ad
import scanpy as sc
adata = sc.read_10x_h5(DATA_DIR / "spatial/filtered_feature_bc_matrix.h5")
```

---

## Data Generation Details

**All data is 100% synthetic** and generated using:

1. **Clinical FHIR**: Generated with FHIR test data generator
2. **Genomic VCF**: Simulated with realistic VAF distributions (beta distribution)
3. **Multi-omics**: Sampled from TCGA ovarian cancer cohort + added noise
4. **Spatial**: Synthetic Visium data with spatial autocorrelation (geoR package)
5. **Imaging**: Generated with synthetic microscopy tools + DeepCell segmentation

**No real patient data was used**. Safe to share publicly, no IRB approval required.

---

## Usage in Book Chapters

| Chapter | Data Used | Purpose |
|---------|-----------|---------|
| 1 | All modalities | Complete PatientOne story |
| 3 | All modalities | Production validation |
| 4 | Clinical FHIR | FHIR integration |
| 5 | Genomic VCF | Variant calling, annotation |
| 6 | Multi-omics TSV | HAllA, Stouffer meta-analysis |
| 7 | Spatial H5AD | STAR, ComBat, Moran's I |
| 8 | Imaging TIFF | DeepCell segmentation |
| 9 | Multi-omics + Spatial | GEARS treatment prediction |
| 10 | Spatial H5AD | Quantum fidelity |
| 11 | Imaging TIFF | Histopathology |
| 15 | All modalities | Research workflows |
| 16 | All modalities | Teaching exercises |

---

## Companion Jupyter Notebooks

All 18 chapter notebooks use PatientOne data. See:
- **Notebook setup**: [`docs/book/companion-notebooks/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/book/companion-notebooks/README.md)
- **Data loading examples**: Each chapter's notebook

---

## Expected Analysis Results

**Complete analysis time**: 35 minutes (AI-orchestrated) vs 40 hours (traditional)

**Key findings**:
- TP53 R175H mutation (85% VAF) → Poor prognosis
- PIK3CA E545K (42% VAF) → PI3K inhibitor responsive
- TMB 8.2 mutations/Mb → High (immunotherapy candidate)
- 25% Ki67+ proliferating cells → Aggressive tumor
- CD8+ T cell infiltration (11%) → Moderate immune response
- Hypoxic regions (12%) → Potential checkpoint blockade target

**Treatment recommendation**: Olaparib (PARP inhibitor) + alpelisib (PI3K inhibitor) combination (82% predicted efficacy)

---

## Additional Resources

- **Full dataset README**: [`data/patient-data/PAT001-OVC-2025/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/data/patient-data/PAT001-OVC-2025/README.md)
- **Analysis scripts**: [`docs/test-docs/patient-one-scenario/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/docs/test-docs/patient-one-scenario)
- **Demo walkthrough**: [`docs/demos/FULL_PATIENTONE_DEMO.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/demos/FULL_PATIENTONE_DEMO.md)

---

**This appendix provides complete reference for the PatientOne synthetic dataset. All data is public domain (CC0 1.0 Universal).**

---
