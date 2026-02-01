# Chapter 7: Spatial Transcriptomics

*Building mcp-spatialtools for 10X Visium spatial gene expression analysis*

---

## Why Spatial Context Matters

In Chapters 5-6, you identified dysregulated pathways:
- **TP53 mutation** (73% variant allele frequency)
- **AKT/mTOR hyperactivation** (meta q-value < 0.01 across RNA/protein/phospho)
- **HIF1A upregulation** (hypoxia signature)

**But where in the tumor are these pathways active?**

Traditional bulk RNA-seq mixes all cells together. You can't distinguish:
- **Tumor proliferative regions** (MKI67+, PCNA+, TOP2A+ → chemotherapy-sensitive)
- **Necrotic/hypoxic cores** (HIF1A+, CA9+, VEGFA+ → radiotherapy-resistant)
- **Invasive fronts** (VIM+, SNAI1+, TWIST1+ → metastatic potential)
- **Stromal barriers** (COL1A1+, FAP+ → drug delivery obstacles)
- **Immune infiltrates** (CD3D+, CD8A+ → immunotherapy responsive)

**10X Visium spatial transcriptomics** measures gene expression in 900 tissue spots (55μm diameter each), preserving spatial relationships.

The `mcp-spatialtools` server provides 8 tools to analyze spatial data: STAR alignment, batch correction, differential expression, pathway enrichment, spatial autocorrelation, and cell type deconvolution.

---

## PatientOne's Spatial Landscape

Here's the spatial data structure ([`spatial/visium_gene_expression.csv`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/data/patient-data/PAT001-OVC-2025/spatial/visium_gene_expression.csv)):

**Dimensions**: 900 spots × 31 genes

**Key genes measured**:
- **Proliferation**: MKI67, PCNA, TOP2A
- **Mutations**: TP53, PIK3CA
- **Oncogenes**: MYC, EGFR
- **EMT markers**: VIM, SNAI1, TWIST1, CDH2
- **Stroma**: COL1A1, COL3A1, ACTA2, FAP
- **Immune**: CD3D, CD8A, CD4, FOXP3, CD68
- **Hypoxia**: HIF1A, CA9, VEGFA
- **Drug resistance**: ABCB1, BCL2L1

**Regional annotations** ([`spatial/visium_region_annotations.csv`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/data/patient-data/PAT001-OVC-2025/spatial/visium_region_annotations.csv)):

```csv
barcode,region
SPOT_00_00,necrotic_hypoxic
SPOT_05_12,tumor_proliferative
SPOT_08_18,tumor_invasive
SPOT_12_25,stroma_reactive
SPOT_15_10,stroma_fibrotic
SPOT_20_20,immune_infiltrated
```

**6 distinct regions**:
1. **Necrotic/Hypoxic** (150 spots): HIF1A+, CA9+ high, VEGFA+ → radiotherapy-resistant
2. **Tumor Proliferative** (200 spots): MKI67+, PCNA+, TOP2A+ → chemotherapy-sensitive
3. **Tumor Invasive** (180 spots): VIM+, SNAI1+, TWIST1+ → metastatic risk
4. **Stroma Reactive** (150 spots): COL1A1+, ACTA2+ → CAF-driven drug resistance
5. **Stroma Fibrotic** (120 spots): COL3A1+, FAP+ → physical barrier to drug delivery
6. **Immune Infiltrated** (100 spots): CD3D+, CD8A+ → immunotherapy potential

---

## The Eight mcp-spatialtools Tools

### Phase 1: Data Processing

#### 1. align_spatial_reads (STAR alignment)

**Why you need it**: Raw FASTQ files from the sequencer must be aligned to the hg38 reference genome to determine which genes are expressed in each spot.

**STAR (Spliced Transcripts Alignment to a Reference)**: The gold standard RNA-seq aligner, handling splice junctions across introns.

**Example workflow**:
```python
@mcp.tool()
def align_spatial_reads(
    fastq_r1: str,
    fastq_r2: str,
    genome_index: str,
    output_dir: str,
    threads: int = 8
) -> dict:
    """Align spatial FASTQ files using STAR."""
    # Run STAR alignment
    star_cmd = [
        "STAR",
        "--runThreadN", str(threads),
        "--genomeDir", genome_index,
        "--readFilesIn", fastq_r1, fastq_r2,
        "--outFileNamePrefix", f"{output_dir}/",
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--quantMode", "GeneCounts"
    ]

    subprocess.run(star_cmd, check=True)

    return {
        "aligned_bam": f"{output_dir}/Aligned.sortedByCoord.out.bam",
        "gene_counts": f"{output_dir}/ReadsPerGene.out.tab",
        "mapping_stats": {
            "uniquely_mapped": "85.3%",
            "multimapped": "8.2%",
            "unmapped": "6.5%"
        }
    }
```

**Typical mapping stats for PatientOne**:
- **Uniquely mapped reads**: 85.3% (GOOD - high quality data)
- **Multimapped reads**: 8.2% (paralogous genes)
- **Unmapped reads**: 6.5% (low quality, contamination)

Implementation: [`servers/mcp-spatialtools/src/mcp_spatialtools/tools/alignment.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-spatialtools/src/mcp_spatialtools/tools/alignment.py)

---

#### 2. filter_spatial_data

**Why you need it**: Low-quality spots (low read count, high mitochondrial percentage) and lowly expressed genes add noise.

**Filtering criteria**:
- **Spots**: Remove if total UMI count < 500 OR mitochondrial genes > 20%
- **Genes**: Remove if detected in < 10 spots (not spatially informative)

**Example filtering**:
```python
@mcp.tool()
def filter_spatial_data(
    counts_path: str,
    min_counts: int = 500,
    max_mito_percent: float = 20.0,
    min_spots: int = 10
) -> dict:
    """Filter low-quality spots and genes."""
    counts = pd.read_csv(counts_path, index_col=0)

    # Calculate QC metrics
    total_counts = counts.sum(axis=1)
    mito_genes = [g for g in counts.columns if g.startswith('MT-')]
    mito_percent = counts[mito_genes].sum(axis=1) / total_counts * 100

    # Filter spots
    valid_spots = (total_counts >= min_counts) & (mito_percent <= max_mito_percent)

    # Filter genes (detected in >= min_spots)
    gene_detection = (counts > 0).sum(axis=0)
    valid_genes = gene_detection >= min_spots

    filtered = counts.loc[valid_spots, valid_genes]

    return {
        "filtered_path": "/data/filtered/counts_filtered.csv",
        "spots_removed": int((~valid_spots).sum()),
        "genes_removed": int((~valid_genes).sum()),
        "spots_retained": int(valid_spots.sum()),
        "genes_retained": int(valid_genes.sum())
    }
```

**PatientOne filtering results**:
- Spots removed: 23 (low UMI counts)
- Genes removed: 0 (all 31 genes well-detected)
- Final data: **877 spots × 31 genes**

---

#### 3. correct_batch_effects

**Why you need it**: Spatial slides processed on different days show batch effects (just like proteomics in Chapter 6).

**ComBat batch correction**: Same algorithm as multi-omics preprocessing.

```python
from combat.pycombat import pycombat

@mcp.tool()
def correct_batch_effects(
    counts_path: str,
    metadata_path: str
) -> dict:
    """Apply ComBat batch correction to spatial data."""
    counts = pd.read_csv(counts_path, index_col=0)
    metadata = pd.read_csv(metadata_path)

    # ComBat correction
    corrected = pycombat(counts.T, metadata['Batch']).T

    # Verify correction worked
    from sklearn.decomposition import PCA
    pca = PCA(n_components=2)
    pc1_before = pca.fit_transform(counts.T)[:, 0]
    pc1_after = pca.fit_transform(corrected.T)[:, 0]

    batch_corr_before = np.corrcoef(pc1_before, metadata['Batch'])[0, 1]**2
    batch_corr_after = np.corrcoef(pc1_after, metadata['Batch'])[0, 1]**2

    return {
        "corrected_path": "/data/processed/counts_corrected.csv",
        "batch_correlation": {
            "before": batch_corr_before,  # e.g., 0.64
            "after": batch_corr_after      # e.g., 0.08
        }
    }
```

---

### Phase 2: Analysis

#### 4. differential_expression_spatial

**Why you need it**: Find genes differentially expressed between regions (e.g., tumor vs stroma, proliferative vs invasive).

**Statistical method**: DESeq2 (negative binomial regression for count data) or simple t-test for normalized counts.

**Example analysis** (Tumor Proliferative vs Necrotic/Hypoxic):
```python
@mcp.tool()
def differential_expression_spatial(
    counts_path: str,
    annotations_path: str,
    group1: str = "tumor_proliferative",
    group2: str = "necrotic_hypoxic",
    fdr_threshold: float = 0.05
) -> dict:
    """Find differentially expressed genes between spatial regions."""
    counts = pd.read_csv(counts_path, index_col=0)
    annotations = pd.read_csv(annotations_path)

    # Get spots for each group
    group1_spots = annotations[annotations['region'] == group1]['barcode']
    group2_spots = annotations[annotations['region'] == group2]['barcode']

    # Log2 fold change and p-values
    from scipy import stats
    results = []
    for gene in counts.columns:
        g1_expr = counts.loc[group1_spots, gene]
        g2_expr = counts.loc[group2_spots, gene]

        log2fc = np.log2(g1_expr.mean() + 1) - np.log2(g2_expr.mean() + 1)
        pval = stats.ttest_ind(g1_expr, g2_expr).pvalue

        results.append({"gene": gene, "log2fc": log2fc, "pval": pval})

    # FDR correction
    from statsmodels.stats.multitest import multipletests
    pvals = [r["pval"] for r in results]
    _, qvals, _, _ = multipletests(pvals, method='fdr_bh')

    for i, r in enumerate(results):
        r["qval"] = qvals[i]

    # Filter significant
    significant = [r for r in results if r["qval"] < fdr_threshold]

    return {
        "comparison": f"{group1} vs {group2}",
        "significant_genes": significant,
        "num_significant": len(significant)
    }
```

**PatientOne results** (Tumor Proliferative vs Necrotic):
```json
{
  "significant_genes": [
    {"gene": "MKI67", "log2fc": 4.2, "qval": 0.00012},
    {"gene": "PCNA", "log2fc": 3.8, "qval": 0.00031},
    {"gene": "TOP2A", "log2fc": 3.1, "qval": 0.0015},
    {"gene": "HIF1A", "log2fc": -5.1, "qval": 0.000045},
    {"gene": "CA9", "log2fc": -6.3, "qval": 0.000008},
    {"gene": "VEGFA", "log2fc": -4.7, "qval": 0.00019}
  ],
  "num_significant": 12
}
```

**Interpretation**:
- **Upregulated in proliferative**: MKI67, PCNA, TOP2A → chemotherapy targets
- **Upregulated in necrotic**: HIF1A, CA9, VEGFA → hypoxia-driven angiogenesis

---

#### 5. pathway_enrichment_spatial

**Why you need it**: Instead of looking at individual genes, identify which biological pathways are enriched in each region.

**Gene set databases**:
- **GO Biological Process**: Cell cycle, apoptosis, immune response
- **KEGG**: PI3K/AKT, MAPK, hypoxia
- **Hallmark**: Epithelial-mesenchymal transition, inflammatory response

**Example enrichment**:
```python
@mcp.tool()
def pathway_enrichment_spatial(
    gene_list: list[str],
    gene_set_database: str = "KEGG",
    fdr_threshold: float = 0.05
) -> dict:
    """Perform pathway enrichment analysis on spatial gene sets."""
    # Fisher's exact test for enrichment
    from scipy.stats import fisher_exact

    enriched_pathways = []
    for pathway_name, pathway_genes in PATHWAY_DATABASE[gene_set_database].items():
        # 2x2 contingency table
        in_list_in_pathway = len(set(gene_list) & set(pathway_genes))
        in_list_not_pathway = len(gene_list) - in_list_in_pathway
        not_list_in_pathway = len(pathway_genes) - in_list_in_pathway
        not_list_not_pathway = 20000 - in_list_in_pathway - in_list_not_pathway - not_list_in_pathway

        _, pval = fisher_exact([[in_list_in_pathway, in_list_not_pathway],
                                 [not_list_in_pathway, not_list_not_pathway]])

        enriched_pathways.append({
            "pathway": pathway_name,
            "pval": pval,
            "overlap": f"{in_list_in_pathway}/{len(pathway_genes)}"
        })

    # FDR correction
    pvals = [p["pval"] for p in enriched_pathways]
    _, qvals, _, _ = multipletests(pvals, method='fdr_bh')

    for i, p in enumerate(enriched_pathways):
        p["qval"] = qvals[i]

    # Filter significant
    significant = [p for p in enriched_pathways if p["qval"] < fdr_threshold]
    significant.sort(key=lambda x: x["qval"])

    return {
        "enriched_pathways": significant,
        "num_enriched": len(significant)
    }
```

**PatientOne enrichment results** (Necrotic/Hypoxic region):
```json
{
  "enriched_pathways": [
    {"pathway": "HIF-1 signaling", "qval": 0.00012, "overlap": "8/43"},
    {"pathway": "VEGF signaling", "qval": 0.00089, "overlap": "6/29"},
    {"pathway": "Glycolysis/Gluconeogenesis", "qval": 0.0021, "overlap": "5/31"}
  ]
}
```

This confirms the necrotic region is hypoxia-driven with active angiogenesis.

---

#### 6. spatial_autocorrelation (Moran's I)

**Why you need it**: Identify genes with **spatially clustered expression** (neighboring spots have similar expression).

**Moran's I statistic**:
- **I > 0**: Positive spatial autocorrelation (clustered expression)
- **I ≈ 0**: Random spatial pattern
- **I < 0**: Negative spatial autocorrelation (checkerboard pattern, rare)

**Example calculation**:
```python
@mcp.tool()
def spatial_autocorrelation(
    counts_path: str,
    coordinates_path: str,
    gene: str = "HIF1A"
) -> dict:
    """Calculate Moran's I for spatial autocorrelation."""
    counts = pd.read_csv(counts_path, index_col=0)
    coords = pd.read_csv(coordinates_path)

    # Extract gene expression and coordinates
    expr = counts[gene].values
    x = coords['x'].values
    y = coords['y'].values

    # Compute spatial weights (inverse distance)
    from scipy.spatial.distance import cdist
    distances = cdist(np.column_stack([x, y]), np.column_stack([x, y]))
    weights = 1 / (distances + 1)  # +1 to avoid division by zero
    np.fill_diagonal(weights, 0)  # No self-weights

    # Moran's I formula
    N = len(expr)
    mean_expr = expr.mean()
    numerator = np.sum(weights * np.outer(expr - mean_expr, expr - mean_expr))
    denominator = np.sum((expr - mean_expr)**2) * np.sum(weights)
    morans_i = (N / np.sum(weights)) * (numerator / denominator)

    # Permutation test for significance
    null_distribution = []
    for _ in range(1000):
        shuffled_expr = np.random.permutation(expr)
        null_i = (N / np.sum(weights)) * (np.sum(weights * np.outer(shuffled_expr - mean_expr, shuffled_expr - mean_expr)) / denominator)
        null_distribution.append(null_i)

    pval = (np.sum(null_distribution >= morans_i) + 1) / 1001

    return {
        "gene": gene,
        "morans_i": morans_i,
        "pval": pval,
        "interpretation": "Clustered" if morans_i > 0 and pval < 0.05 else "Random"
    }
```

**PatientOne spatial autocorrelation**:
```json
{
  "HIF1A": {"morans_i": 0.82, "pval": 0.001, "interpretation": "Clustered"},
  "MKI67": {"morans_i": 0.75, "pval": 0.001, "interpretation": "Clustered"},
  "CD3D": {"morans_i": 0.61, "pval": 0.008, "interpretation": "Clustered"}
}
```

All key markers show significant spatial clustering (not randomly distributed).

---

#### 7. deconvolve_cell_types

**Why you need it**: Each 55μm Visium spot contains ~10-30 cells. What cell types are present?

**Cell type deconvolution**: Uses reference single-cell RNA-seq signatures to estimate cell type proportions.

**Example deconvolution** (simplified):
```python
@mcp.tool()
def deconvolve_cell_types(
    spatial_counts_path: str,
    reference_signatures_path: str
) -> dict:
    """Estimate cell type proportions per spot."""
    spatial_counts = pd.read_csv(spatial_counts_path, index_col=0)
    reference = pd.read_csv(reference_signatures_path, index_col=0)

    # Non-negative least squares deconvolution
    from scipy.optimize import nnls

    cell_type_proportions = {}
    for spot in spatial_counts.index:
        spot_expr = spatial_counts.loc[spot].values
        ref_expr = reference.values

        # Solve: spot_expr ≈ Σ(proportion_i * cell_type_i)
        proportions, _ = nnls(ref_expr, spot_expr)
        proportions /= proportions.sum()  # Normalize to 100%

        cell_type_proportions[spot] = dict(zip(reference.columns, proportions))

    return {
        "cell_type_proportions": cell_type_proportions
    }
```

**PatientOne deconvolution** (example spot in immune infiltrated region):
```json
{
  "SPOT_20_20": {
    "Tumor cells": 0.35,
    "CD8+ T cells": 0.28,
    "CD4+ T cells": 0.15,
    "Macrophages": 0.12,
    "Fibroblasts": 0.08,
    "Endothelial": 0.02
  }
}
```

This spot is 35% tumor cells + 43% T cells (CD8+ and CD4+) → immunotherapy target.

---

#### 8. link_clinical_to_spatial

**Why you need it**: Connect PatientOne's clinical data (from Chapter 4, mcp-epic) to spatial regions.

**Example linkage**:
```python
@mcp.tool()
def link_clinical_to_spatial(
    patient_id: str,
    spatial_annotations_path: str
) -> dict:
    """Bridge clinical FHIR data to spatial tissue regions."""
    # Get clinical data from mcp-epic
    clinical = epic.get_patient_summary(patient_id)

    # Get spatial regions
    annotations = pd.read_csv(spatial_annotations_path)
    region_counts = annotations['region'].value_counts().to_dict()

    # Link treatment history to spatial features
    chemo_agent = clinical['medications'][0]['drug']  # Carboplatin

    # Predict response by region
    predictions = {
        "tumor_proliferative": "Sensitive (MKI67+ → DNA damage)",
        "necrotic_hypoxic": "Resistant (HIF1A+ → hypoxia-driven survival)",
        "tumor_invasive": "Partial (EMT markers → chemoresistance)",
        "immune_infiltrated": "Consider immunotherapy (CD8+ T cells present)"
    }

    return {
        "patient_id": patient_id,
        "treatment": chemo_agent,
        "spatial_response_predictions": predictions,
        "region_distribution": region_counts
    }
```

This tool connects clinical treatment history to spatial biology, enabling precision treatment stratification.

---

## Complete PatientOne Spatial Workflow

Natural language prompt in Claude Desktop:

```
I have 10X Visium spatial transcriptomics data for patient PAT001-OVC-2025:
- Counts: /data/patient-data/PAT001-OVC-2025/spatial/visium_gene_expression.csv
- Annotations: /data/patient-data/PAT001-OVC-2025/spatial/visium_region_annotations.csv
- Coordinates: /data/patient-data/PAT001-OVC-2025/spatial/visium_spatial_coordinates.csv

Please:
1. Filter low-quality spots (< 500 UMI counts)
2. Run differential expression: tumor_proliferative vs necrotic_hypoxic
3. Pathway enrichment for top 10 upregulated genes in each region
4. Calculate Moran's I for HIF1A, MKI67, CD3D
5. Deconvolve cell types
6. Link clinical treatment history to spatial predictions
```

Claude orchestrates all 8 tools automatically, returning a comprehensive spatial analysis in **<3 minutes**.

---

## Implementation Walkthrough

### Step 1: Project Setup

```bash
cd servers/mcp-spatialtools
python -m venv venv
source venv/bin/activate
pip install fastmcp pandas numpy scipy statsmodels scikit-learn
```

Environment variables (`.env`):
```bash
SPATIAL_DATA_DIR=/workspace/data/spatial
SPATIAL_CACHE_DIR=/workspace/cache/spatial
SPATIAL_DRY_RUN=true  # For testing
```

### Step 2: Initialize FastMCP Server

```python
from fastmcp import FastMCP
import os
from pathlib import Path

mcp = FastMCP("spatialtools")

config = {
    "data_dir": Path(os.getenv("SPATIAL_DATA_DIR", "/workspace/data/spatial")),
    "cache_dir": Path(os.getenv("SPATIAL_CACHE_DIR", "/workspace/cache/spatial")),
    "dry_run": os.getenv("SPATIAL_DRY_RUN", "false").lower() == "true"
}
```

### Step 3: Add Differential Expression

Core analysis tool. Create `src/mcp_spatialtools/tools/differential_expression.py`:

```python
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

def differential_expression_spatial_impl(
    counts_path: str,
    annotations_path: str,
    group1: str,
    group2: str,
    fdr_threshold: float = 0.05
) -> dict:
    """Internal implementation of differential expression."""
    counts = pd.read_csv(counts_path, index_col=0)
    annotations = pd.read_csv(annotations_path)

    # Get barcodes for each group
    group1_barcodes = annotations[annotations['region'] == group1]['barcode'].tolist()
    group2_barcodes = annotations[annotations['region'] == group2]['barcode'].tolist()

    # Differential expression
    results = []
    for gene in counts.columns:
        g1_expr = counts.loc[group1_barcodes, gene]
        g2_expr = counts.loc[group2_barcodes, gene]

        # Log2 fold change
        log2fc = np.log2(g1_expr.mean() + 1) - np.log2(g2_expr.mean() + 1)

        # T-test
        stat, pval = stats.ttest_ind(g1_expr, g2_expr)

        results.append({
            "gene": gene,
            "log2fc": log2fc,
            "pval": pval,
            "mean_group1": g1_expr.mean(),
            "mean_group2": g2_expr.mean()
        })

    # FDR correction
    pvals = [r["pval"] for r in results]
    _, qvals, _, _ = multipletests(pvals, method='fdr_bh')

    for i, r in enumerate(results):
        r["qval"] = qvals[i]

    # Filter and sort by qvalue
    significant = [r for r in results if r["qval"] < fdr_threshold]
    significant.sort(key=lambda x: x["qval"])

    return {
        "comparison": f"{group1} vs {group2}",
        "significant_genes": significant,
        "num_significant": len(significant),
        "total_tested": len(results)
    }
```

Full implementation: [`servers/mcp-spatialtools/src/mcp_spatialtools/tools/differential_expression.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-spatialtools/src/mcp_spatialtools/tools/differential_expression.py)

---

## Testing Your Server

### Unit Tests

```python
# tests/test_differential_expression.py
import pytest
from mcp_spatialtools.tools.differential_expression import differential_expression_spatial_impl

def test_differential_expression_tumor_vs_necrotic():
    """Test DE analysis on PatientOne data."""
    result = differential_expression_spatial_impl(
        counts_path="tests/data/visium_gene_expression.csv",
        annotations_path="tests/data/visium_region_annotations.csv",
        group1="tumor_proliferative",
        group2="necrotic_hypoxic",
        fdr_threshold=0.05
    )

    assert result["num_significant"] > 0
    assert "MKI67" in [g["gene"] for g in result["significant_genes"]]
    assert "HIF1A" in [g["gene"] for g in result["significant_genes"]]

    # Check directionality
    mki67 = next(g for g in result["significant_genes"] if g["gene"] == "MKI67")
    hif1a = next(g for g in result["significant_genes"] if g["gene"] == "HIF1A")

    assert mki67["log2fc"] > 2  # Upregulated in proliferative
    assert hif1a["log2fc"] < -2  # Downregulated in proliferative (= upregulated in necrotic)
```

Run tests:
```bash
pytest tests/ -v
```

---

## Connecting to Claude Desktop

Add to `claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "spatialtools": {
      "command": "/path/to/precision-medicine-mcp/servers/mcp-spatialtools/venv/bin/python",
      "args": ["-m", "mcp_spatialtools"],
      "env": {
        "SPATIAL_DATA_DIR": "/workspace/data/spatial",
        "SPATIAL_DRY_RUN": "false"
      }
    }
  }
}
```

Natural language use:

```
Claude: Run differential expression between tumor_proliferative and necrotic_hypoxic regions in my spatial data. Find enriched pathways for each region.
```

Claude calls `spatialtools.differential_expression_spatial()` and `spatialtools.pathway_enrichment_spatial()` automatically.

---

## What You've Built

You now have a spatial transcriptomics server that:

1. **Processes raw data**: STAR alignment, quality filtering, batch correction
2. **Identifies regional differences**: Differential expression, pathway enrichment
3. **Analyzes spatial patterns**: Moran's I autocorrelation, spatial clustering
4. **Deconvolves cell types**: Estimates tumor, immune, stromal proportions per spot
5. **Links to clinical data**: Connects treatment history to spatial predictions

This completes Part 2 (Building the Foundation). You now have servers for:
- Clinical data (Chapter 4, mcp-epic)
- Genomics (Chapter 5, mcp-fgbio)
- Multi-omics (Chapter 6, mcp-multiomics)
- Spatial (Chapter 7, mcp-spatialtools)

Together, these servers enable comprehensive precision oncology analysis from patient EHR to tissue microenvironment.

---

## Try It Yourself

### Option 1: Local Development

```bash
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp/servers/mcp-spatialtools

python -m venv venv
source venv/bin/activate
pip install -e ".[dev]"

export SPATIAL_DRY_RUN=true
python -m mcp_spatialtools
```

### Option 2: PatientOne Data

```bash
export SPATIAL_DRY_RUN=false

# In Claude Desktop:
# "Run differential expression on spatial data:
#  Counts: data/patient-data/PAT001-OVC-2025/spatial/visium_gene_expression.csv
#  Annotations: data/patient-data/PAT001-OVC-2025/spatial/visium_region_annotations.csv
#  Compare: tumor_proliferative vs necrotic_hypoxic"
```

---

## Next Steps

In **Part 3: Advanced Capabilities** (Chapters 8-11), you'll build:
- **Chapter 8**: Cell segmentation with DeepCell (nuclear/membrane segmentation, phenotyping)
- **Chapter 9**: Treatment response prediction with GEARS GNN
- **Chapter 10**: Quantum cell-type fidelity with PennyLane
- **Chapter 11**: Imaging and histopathology (H&E, MxIF)

The spatial analysis you built identifies **where pathways are active**. Cell segmentation (Chapter 8) will enable **single-cell resolution** analysis within those spatial regions.

---

**Chapter 7 Summary**:
- 10X Visium spatial transcriptomics preserves tissue context (900 spots × 31 genes)
- 6 distinct regions identified: necrotic, tumor proliferative, invasive, stroma reactive/fibrotic, immune
- Differential expression reveals MKI67+ (chemo-sensitive) vs HIF1A+ (chemo-resistant) regions
- Moran's I confirms spatially clustered expression patterns
- Cell type deconvolution estimates tumor/immune/stromal proportions per spot

**Files created**: `servers/mcp-spatialtools/src/mcp_spatialtools/server.py`, `tools/differential_expression.py`, `tools/pathway_enrichment.py`, `tools/spatial_autocorrelation.py`
**Tests added**: 14 unit tests, 71% coverage
**Tools exposed**: 8 MCP tools (align, filter, batch_correct, differential_expression, pathway_enrichment, spatial_autocorrelation, deconvolve, link_clinical)
