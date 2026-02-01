# Chapter 6: Multi-Omics Integration

*Building mcp-multiomics for RNA, protein, and phosphoproteomics analysis*

---

## The Multi-Omics Problem

You have three data types from PatientOne's PDX models (15 ovarian cancer xenografts):

1. **RNA-seq** (transcriptomics): What genes are transcribed?
2. **Proteomics**: What proteins are actually made?
3. **Phosphoproteomics**: What signaling pathways are active?

**The challenge**: A gene upregulated at the RNA level might not show elevated protein (post-transcriptional regulation), and high protein doesn't guarantee phosphorylation-driven activity (kinase regulation).

**Example from PatientOne**:
- **TP53 RNA**: Low expression (mutated copy)
- **TP53 protein**: Paradoxically high (misfolded protein accumulates)
- **AKT phosphorylation**: High (PI3K pathway activated despite TP53 dysfunction)

You can't make treatment decisions from RNA alone. You need **multi-omics meta-analysis** to find which pathways are consistently dysregulated across all three modalities.

The `mcp-multiomics` server provides 7 tools to integrate, validate, and analyze multi-omics data using HAllA association testing and Stouffer meta-analysis.

---

## The Seven Tools

### Phase 1: Data Validation (Run FIRST)

#### 1. validate_multiomics_data

**Why you need it**: Real-world proteomics data has batch effects (mass spec runs ~18 samples at a time), missing values (~20-30% typical), and sample naming inconsistencies.

**What it checks**:
- Sample overlap across modalities
- Missing value patterns (are they random or systematic?)
- Batch effect detection (PC1 correlated with batch?)
- Outlier samples (median absolute deviation > 3)

**Example validation result**:
```json
{
  "validation_status": "WARNING",
  "sample_overlap": {
    "rna_only": 0,
    "protein_only": 0,
    "common": 15
  },
  "missing_patterns": {
    "rna": "0.5%",
    "protein": "23.4%",
    "phospho": "31.2%"
  },
  "batch_effects": {
    "protein_batch_correlation": 0.68,
    "phospho_batch_correlation": 0.71
  },
  "warnings": [
    "Protein batch effect detected (PC1 ~ Batch, R²=0.68)",
    "Phospho missing >30% - consider imputation"
  ],
  "recommendations": [
    "Apply ComBat batch correction for protein and phospho",
    "Use KNN imputation for missing values"
  ]
}
```

**Natural language use**: *"Validate my multi-omics data quality: RNA at `/data/pdx_rna_seq.csv`, Protein at `/data/pdx_proteomics.csv`. Check for batch effects and sample naming issues."*

---

#### 2. preprocess_multiomics_data

**Why you need it**: Batch effects confound biological signal. Without correction, you'll detect "batch 1 vs batch 2" instead of "treatment response vs resistant".

**Preprocessing pipeline**:
1. **ComBat batch correction** (removes batch-to-batch variation)
2. **KNN imputation** (fills missing values using k-nearest neighbors)
3. **Quantile normalization** (makes sample distributions comparable)
4. **Outlier removal** (MAD > 3 threshold)

**Code example** (simplified):
```python
@mcp.tool()
def preprocess_multiomics_data(
    rna_path: str,
    protein_path: str,
    metadata_path: str,
    batch_correction: bool = True,
    normalize_method: str = "quantile",
    imputation_method: str = "knn"
) -> dict:
    """Preprocess multi-omics data with batch correction."""
    # Load data
    rna = pd.read_csv(rna_path, index_col=0)
    protein = pd.read_csv(protein_path, index_col=0)
    metadata = pd.read_csv(metadata_path)

    # Apply ComBat batch correction
    if batch_correction:
        protein_corrected = combat(protein, batch=metadata['Batch'])
    else:
        protein_corrected = protein

    # KNN imputation
    imputer = KNNImputer(n_neighbors=5)
    protein_imputed = imputer.fit_transform(protein_corrected.T).T

    # Quantile normalization
    protein_normalized = quantile_normalize(protein_imputed)

    return {
        "preprocessed_paths": {
            "protein": "/data/preprocessed/protein_processed.csv"
        },
        "batch_correction_results": {
            "before_pc1_batch_corr": 0.68,
            "after_pc1_batch_corr": 0.12
        },
        "imputation_stats": {
            "values_imputed": 2341
        }
    }
```

**Before/after batch correction**:
```
Before: PC1 explains 45% variance, R² with Batch = 0.68 (BAD)
After:  PC1 explains 38% variance, R² with Batch = 0.12 (GOOD)
```

Implementation: [`servers/mcp-multiomics/src/mcp_multiomics/tools/preprocessing.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-multiomics/src/mcp_multiomics/tools/preprocessing.py)

---

#### 3. visualize_data_quality

Generates QC plots to verify preprocessing worked:
- **PCA plots**: Samples colored by batch (should cluster by phenotype, not batch)
- **Missing value heatmaps**: Before/after imputation
- **Correlation matrices**: Inter-sample correlation

**Example plot interpretation**:
```
Before preprocessing: Samples cluster by MS batch (red/blue groups)
After preprocessing:  Samples cluster by treatment response (responders/resistant)
```

---

### Phase 2: Integration and Analysis

#### 4. integrate_omics_data

Aligns samples across modalities and normalizes within each modality.

**What it does**:
- Finds common samples (PatientOne: 15 samples present in all 3 modalities)
- Applies Z-score normalization per modality
- Filters features with >50% missing
- Returns aligned data matrices

**Example output**:
```json
{
  "common_samples": ["Sample_01", "Sample_02", ..., "Sample_15"],
  "feature_counts": {
    "rna": 1850,
    "protein": 1243,
    "phospho": 987
  },
  "qc_metrics": {
    "samples_aligned": 15,
    "features_filtered": {
      "rna": 150,
      "protein": 87,
      "phospho": 123
    }
  }
}
```

---

#### 5. run_halla_analysis

**HAllA (HMP All-Against-All Association)**: Tests all pairwise associations between features within a modality (e.g., which RNA-RNA pairs correlate?).

**The scalability problem**:
- 1850 RNA features × 1850 = **3.4 million pairwise tests**
- Naive implementation: **Days of compute time**

**Solution: Chunking strategy** (bioinformatician feedback, 2025):
```python
# Chunk features into groups of 1000
chunks = [features[i:i+1000] for i in range(0, len(features), 1000)]

# Run HAllA on each chunk (1000×1000 = 1M tests, ~5 minutes)
for chunk in chunks:
    halla_results = halla.run(chunk, chunk, method='spearman')
```

**Result**: 3.4 million tests in ~15 minutes instead of days.

**Natural language use**: *"Run HAllA on the RNA data at `/data/pdx_rna_seq.csv` to find co-expressed gene pairs"*

Implementation: [`servers/mcp-multiomics/src/mcp_multiomics/tools/halla.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-multiomics/src/mcp_multiomics/tools/halla.py)

---

#### 6. calculate_stouffer_meta

**Stouffer's Meta-Analysis**: Combines p-values across modalities to find consistently dysregulated features.

**The critical workflow** (per bioinformatician review):

```
┌─────────────────────────────────────────────────────────────────┐
│ CORRECT FDR TIMING FOR MULTI-OMICS META-ANALYSIS               │
├─────────────────────────────────────────────────────────────────┤
│ Step 1: HAllA Analysis                                          │
│         → Use NOMINAL p-values (NOT FDR-corrected)             │
│         → Returns: p_value_nominal for each association         │
│                                                                  │
│ Step 2: Stouffer's Meta-Analysis (THIS TOOL)                   │
│         → Input: NOMINAL p-values from each modality           │
│         → Combine p-values across modalities                    │
│         → Output: meta_p_values (still nominal)                 │
│                                                                  │
│ Step 3: FDR Correction (APPLIED AFTER COMBINATION)             │
│         → Input: meta_p_values from Step 2                      │
│         → Apply: Benjamini-Hochberg FDR                         │
│         → Output: q_values (FDR-corrected)                      │
└─────────────────────────────────────────────────────────────────┘

WHY THIS ORDER MATTERS:
- Applying FDR before Stouffer's → loses power (over-conservative)
- Applying FDR after Stouffer's → correct statistical framework
```

**Stouffer's Z-score method**:
1. Convert p-values to Z-scores: `Z = Φ⁻¹(1 - p/2)`
2. Combine Z-scores (weighted by sample size): `Z_meta = Σ(w_i × Z_i) / √(Σw_i²)`
3. Convert back to p-value: `p_meta = 2 × (1 - Φ(|Z_meta|))`

**Example implementation**:
```python
class StoufferMetaAnalysis:
    """Combine p-values from multiple omics modalities."""

    def p_to_z(self, p_values: np.ndarray, effect_sizes: np.ndarray = None) -> np.ndarray:
        """Convert p-values to Z-scores."""
        # Clip to avoid numerical issues
        p_values = np.clip(p_values, 1e-300, 1 - 1e-15)

        # Two-tailed Z-scores
        z_scores = stats.norm.ppf(1 - p_values / 2)

        # Apply directionality from effect sizes
        if effect_sizes is not None:
            signs = np.sign(effect_sizes)
            z_scores = z_scores * signs

        return z_scores

    def combine_z_scores(self, z_scores_list: List[np.ndarray], weights: List[float]) -> np.ndarray:
        """Combine Z-scores from multiple modalities."""
        # Weighted sum
        z_meta = sum(w * z for w, z in zip(weights, z_scores_list))

        # Normalize by sqrt of sum of squared weights
        z_meta /= np.sqrt(sum(w**2 for w in weights))

        return z_meta

    def z_to_p(self, z_scores: np.ndarray) -> np.ndarray:
        """Convert Z-scores back to two-tailed p-values."""
        return 2 * (1 - stats.norm.cdf(np.abs(z_scores)))

    def apply_fdr_correction(self, p_values: np.ndarray) -> np.ndarray:
        """Apply Benjamini-Hochberg FDR correction (AFTER meta-analysis)."""
        _, q_values, _, _ = multipletests(p_values, method='fdr_bh')
        return q_values
```

**Example result** (AKT pathway from PatientOne):
```json
{
  "feature": "AKT1",
  "modality_results": {
    "rna": {"p_value": 0.023, "log2fc": 1.2},
    "protein": {"p_value": 0.041, "log2fc": 0.8},
    "phospho": {"p_value": 0.0087, "log2fc": 2.1}
  },
  "meta_analysis": {
    "meta_p_value": 0.00031,
    "meta_q_value": 0.012,
    "combined_z_score": 3.58,
    "direction": "UP"
  }
}
```

AKT1 shows consistent upregulation across all 3 modalities (meta q-value = 0.012 < 0.05), making it a high-confidence therapeutic target.

Implementation: [`servers/mcp-multiomics/src/mcp_multiomics/tools/stouffer.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-multiomics/src/mcp_multiomics/tools/stouffer.py)

---

#### 7. predict_upstream_regulators

**Why you need it**: Finding dysregulated proteins is only half the story. You need to know **what controls them** to identify druggable targets.

**What it predicts**:
- **Kinases**: Which kinases phosphorylate the dysregulated proteins?
- **Transcription factors**: Which TFs control gene expression?
- **Drug responses**: Which drugs target these regulators?

**Example** (PatientOne AKT pathway):
```json
{
  "upstream_kinases": [
    {"kinase": "PDK1", "confidence": 0.92, "targets": ["AKT1", "AKT2"]},
    {"kinase": "mTORC2", "confidence": 0.87, "targets": ["AKT1"]}
  ],
  "drug_targets": [
    {"drug": "Capivasertib", "target": "AKT1/2/3", "clinical_phase": "Phase 3"},
    {"drug": "Ipatasertib", "target": "AKT1/2/3", "clinical_phase": "Phase 2"}
  ]
}
```

This tells you: AKT is hyperphosphorylated → likely by PDK1/mTORC2 → druggable with AKT inhibitors (Capivasertib).

---

## The Complete PatientOne Multi-Omics Workflow

Here's how you'd run the full analysis in Claude Desktop:

```
I have multi-omics data from 15 ovarian cancer PDX models:
- RNA: /data/patient-data/PAT001-OVC-2025/multiomics/pdx_rna_seq.csv
- Protein: /data/patient-data/PAT001-OVC-2025/multiomics/pdx_proteomics.csv
- Phospho: /data/patient-data/PAT001-OVC-2025/multiomics/pdx_phosphoproteomics.csv
- Metadata: /data/patient-data/PAT001-OVC-2025/multiomics/sample_metadata.csv

Please:
1. Validate data quality and check for batch effects
2. Preprocess with ComBat batch correction and KNN imputation
3. Integrate the three modalities
4. Run Stouffer meta-analysis to find consistently dysregulated pathways
5. Predict upstream regulators for the top 10 hits
```

Claude orchestrates all 7 tools automatically:

```python
# Step 1: Validate
validation = multiomics.validate_multiomics_data(
    rna_path="/data/.../pdx_rna_seq.csv",
    protein_path="/data/.../pdx_proteomics.csv",
    phospho_path="/data/.../pdx_phosphoproteomics.csv",
    metadata_path="/data/.../sample_metadata.csv"
)

# Step 2: Preprocess (if validation shows batch effects)
preprocessed = multiomics.preprocess_multiomics_data(
    rna_path="/data/.../pdx_rna_seq.csv",
    protein_path="/data/.../pdx_proteomics.csv",
    metadata_path="/data/.../sample_metadata.csv",
    batch_correction=True,
    normalize_method="quantile",
    imputation_method="knn"
)

# Step 3: Integrate
integrated = multiomics.integrate_omics_data(
    rna_path=preprocessed["preprocessed_paths"]["rna"],
    protein_path=preprocessed["preprocessed_paths"]["protein"],
    phospho_path=preprocessed["preprocessed_paths"]["phospho"]
)

# Step 4: Run HAllA on each modality (chunked)
halla_rna = multiomics.run_halla_analysis(
    data_path=integrated["integrated_data"]["rna"]
)
halla_protein = multiomics.run_halla_analysis(
    data_path=integrated["integrated_data"]["protein"]
)
halla_phospho = multiomics.run_halla_analysis(
    data_path=integrated["integrated_data"]["phospho"]
)

# Step 5: Stouffer meta-analysis
meta_results = multiomics.calculate_stouffer_meta(
    rna_results=halla_rna["associations"],
    protein_results=halla_protein["associations"],
    phospho_results=halla_phospho["associations"],
    fdr_threshold=0.05
)

# Step 6: Upstream regulators
regulators = multiomics.predict_upstream_regulators(
    top_features=meta_results["significant_features"][:10]
)
```

**Results for PatientOne**:
- **15 common samples** across all 3 modalities
- **1850 RNA features**, 1243 proteins, 987 phosphosites
- **23 pathways** with meta q-value < 0.05 (AKT/mTOR, PI3K, TP53, DNA damage response)
- **12 upstream kinases** predicted (PDK1, mTORC2, ATM, ATR, CHK1/2)
- **8 druggable targets** identified (AKT inhibitors, PI3K inhibitors, PARP inhibitors)

---

## Implementation Walkthrough

### Step 1: Project Setup

```bash
cd servers/mcp-multiomics
python -m venv venv
source venv/bin/activate
pip install fastmcp pandas numpy scipy statsmodels scikit-learn
```

Set environment variables (`.env`):
```bash
MULTIOMICS_DRY_RUN=true  # For testing
MULTIOMICS_DATA_DIR=/workspace/data
MULTIOMICS_CACHE_DIR=/workspace/cache
```

### Step 2: Initialize FastMCP Server

Create `src/mcp_multiomics/server.py`:

```python
from fastmcp import FastMCP
from pathlib import Path
import os

mcp = FastMCP("multiomics")

# Configuration
config = {
    "dry_run": os.getenv("MULTIOMICS_DRY_RUN", "false").lower() == "true",
    "data_dir": Path(os.getenv("MULTIOMICS_DATA_DIR", "/workspace/data")),
    "cache_dir": Path(os.getenv("MULTIOMICS_CACHE_DIR", "/workspace/cache"))
}
```

### Step 3: Add Stouffer Meta-Analysis

The core statistical method. Create `src/mcp_multiomics/tools/stouffer.py`:

```python
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

@mcp.tool()
def calculate_stouffer_meta(
    rna_results: dict,
    protein_results: dict,
    phospho_results: dict,
    fdr_threshold: float = 0.05
) -> dict:
    """Combine p-values across modalities using Stouffer's method."""

    # Extract p-values (NOMINAL, not FDR-corrected)
    p_rna = np.array([r["p_value"] for r in rna_results])
    p_protein = np.array([r["p_value"] for r in protein_results])
    p_phospho = np.array([r["p_value"] for r in phospho_results])

    # Convert to Z-scores
    z_rna = stats.norm.ppf(1 - p_rna / 2)
    z_protein = stats.norm.ppf(1 - p_protein / 2)
    z_phospho = stats.norm.ppf(1 - p_phospho / 2)

    # Combine with equal weights (could use sample size weights)
    weights = [1.0, 1.0, 1.0]
    z_meta = (weights[0] * z_rna + weights[1] * z_protein + weights[2] * z_phospho) / np.sqrt(sum(w**2 for w in weights))

    # Convert back to p-values
    p_meta = 2 * (1 - stats.norm.cdf(np.abs(z_meta)))

    # Apply FDR correction AFTER combining (CRITICAL)
    _, q_meta, _, _ = multipletests(p_meta, method='fdr_bh')

    # Find significant features
    significant = q_meta < fdr_threshold

    return {
        "meta_p_values": p_meta.tolist(),
        "meta_q_values": q_meta.tolist(),
        "significant_features": np.where(significant)[0].tolist(),
        "num_significant": int(significant.sum())
    }
```

Full implementation: [`servers/mcp-multiomics/src/mcp_multiomics/tools/stouffer.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-multiomics/src/mcp_multiomics/tools/stouffer.py)

---

## Dry-Run Mode for Safe Testing

Before processing real proteomics data:

```bash
export MULTIOMICS_DRY_RUN=true
python -m mcp_multiomics
```

Dry-run returns synthetic results with warning banners:

```
╔═══════════════════════════════════════════════════════════════════════════╗
║                    ⚠️  SYNTHETIC DATA WARNING ⚠️                          ║
╚═══════════════════════════════════════════════════════════════════════════╝

This result was generated in DRY_RUN mode and does NOT represent real analysis.

🔴 CRITICAL: Do NOT use this data for research decisions or publications.
🔴 All values are SYNTHETIC/MOCKED and have no scientific validity.
```

---

## Testing Your Server

### Unit Tests

```python
# tests/test_stouffer.py
import pytest
import numpy as np
from mcp_multiomics.tools.stouffer import StoufferMetaAnalysis

def test_stouffer_combine_strong_signals():
    """Test combining strong signals from all 3 modalities."""
    stouffer = StoufferMetaAnalysis()

    # Strong signals (p < 0.01) in all 3 modalities
    p_rna = np.array([0.005, 0.008, 0.003])
    p_protein = np.array([0.007, 0.004, 0.009])
    p_phospho = np.array([0.002, 0.006, 0.001])

    # Convert to Z-scores
    z_rna = stouffer.p_to_z(p_rna)
    z_protein = stouffer.p_to_z(p_protein)
    z_phospho = stouffer.p_to_z(p_phospho)

    # Combine
    z_meta = stouffer.combine_z_scores([z_rna, z_protein, z_phospho], [1, 1, 1])
    p_meta = stouffer.z_to_p(z_meta)

    # Meta p-values should be even stronger
    assert np.all(p_meta < 0.001)

def test_stouffer_weak_signal_filtering():
    """Test filtering weak signals (only 1 modality significant)."""
    stouffer = StoufferMetaAnalysis()

    # Strong in RNA, weak in protein/phospho
    p_rna = np.array([0.001])
    p_protein = np.array([0.45])
    p_phospho = np.array([0.67])

    z_meta = stouffer.combine_z_scores(
        [stouffer.p_to_z(p_rna), stouffer.p_to_z(p_protein), stouffer.p_to_z(p_phospho)],
        [1, 1, 1]
    )
    p_meta = stouffer.z_to_p(z_meta)

    # Meta p-value should be non-significant (1 out of 3 modalities)
    assert p_meta[0] > 0.05
```

Run tests:
```bash
pytest tests/test_stouffer.py -v
```

---

## Connecting to Claude Desktop

Add to `claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "multiomics": {
      "command": "/path/to/precision-medicine-mcp/servers/mcp-multiomics/venv/bin/python",
      "args": ["-m", "mcp_multiomics"],
      "env": {
        "MULTIOMICS_DATA_DIR": "/workspace/data",
        "MULTIOMICS_DRY_RUN": "false"
      }
    }
  }
}
```

Natural language use:

```
Claude: Run Stouffer meta-analysis on my RNA, protein, and phospho results to find consistently dysregulated pathways.
```

Claude calls `multiomics.calculate_stouffer_meta()` automatically.

---

## What You've Built

You now have a multi-omics integration server that:

1. **Validates data quality**: Batch effects, missing values, outliers
2. **Preprocesses data**: ComBat batch correction, KNN imputation, normalization
3. **Integrates modalities**: Aligns samples, filters features
4. **Finds associations**: HAllA with chunking (3.4M tests in 15 minutes)
5. **Combines evidence**: Stouffer meta-analysis with correct FDR timing
6. **Predicts regulators**: Kinases, TFs, drug targets
7. **Generates visualizations**: PCA, heatmaps, correlation matrices

This server bridges genomics (Chapter 5) to spatial transcriptomics (Chapter 7) by identifying which pathways to investigate in tissue context.

---

## Try It Yourself

### Option 1: Local Development

```bash
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp/servers/mcp-multiomics

python -m venv venv
source venv/bin/activate
pip install -e ".[dev]"

export MULTIOMICS_DRY_RUN=true
python -m mcp_multiomics
```

### Option 2: PatientOne Data

```bash
export MULTIOMICS_DRY_RUN=false

# In Claude Desktop:
# "Run Stouffer meta-analysis on:
#  - RNA: data/patient-data/PAT001-OVC-2025/multiomics/pdx_rna_seq.csv
#  - Protein: data/patient-data/PAT001-OVC-2025/multiomics/pdx_proteomics.csv
#  - Phospho: data/patient-data/PAT001-OVC-2025/multiomics/pdx_phosphoproteomics.csv"
```

---

## Next Steps

In **Chapter 7: Spatial Transcriptomics**, you'll build `mcp-spatialtools` to analyze 10X Visium data and map the pathways you found here (AKT/mTOR, PI3K) to specific tissue regions.

The multi-omics analysis you built identifies **which pathways** are dysregulated. Spatial transcriptomics reveals **where in the tumor** those pathways are active.

---

**Chapter 6 Summary**:
- Multi-omics integration requires preprocessing (batch correction, imputation)
- HAllA with chunking enables scalable all-vs-all association testing
- Stouffer meta-analysis combines p-values BEFORE FDR correction (critical timing)
- Upstream regulator prediction connects pathways to druggable targets

**Files created**: `servers/mcp-multiomics/src/mcp_multiomics/server.py`, `tools/stouffer.py`, `tools/preprocessing.py`, `tools/halla.py`
**Tests added**: 18 unit tests, 67% coverage
**Tools exposed**: 7 MCP tools (validate, preprocess, visualize, integrate, halla, stouffer, upstream_regulators)
