# Chapter 6: Multi-Omics Integration

*Building mcp-multiomics for RNA, protein, and phosphoproteomics analysis*

---

## The Multi-Omics Problem

PatientOne's PDX models (15 ovarian cancer xenografts) generate three data types:
1. **RNA-seq** (transcriptomics): What genes are transcribed?
2. **Proteomics**: What proteins are made?
3. **Phosphoproteomics**: What signaling pathways are active?

**Challenge**: A gene upregulated at RNA level might not show elevated protein (post-transcriptional regulation). High protein doesn't guarantee phosphorylation-driven activity (kinase regulation).

**PatientOne example**:
- **TP53 RNA**: Low expression (mutated)
- **TP53 protein**: High (misfolded protein accumulates)
- **AKT phosphorylation**: High (PI3K pathway activated despite TP53 dysfunction)

You need **multi-omics meta-analysis** to find consistently dysregulated pathways across all modalities.

The `mcp-multiomics` server provides 7 tools for HAllA association testing and Stouffer meta-analysis.

---

## The Seven Tools

### Phase 1: Data Validation

#### 1. validate_multiomics_data

Checks sample overlap, missing value patterns, batch effects, outlier samples.

```python
@mcp.tool()
def validate_multiomics_data(rna_path: str, protein_path: str, metadata_path: str) -> dict:
    """Validate multi-omics data quality: batch effects, missing values, sample naming."""
    # Load data, check sample overlap, detect batch effects (PC1 ~ Batch correlation)
    # Identify missing patterns, outlier samples (MAD > 3)
    # Full implementation: servers/mcp-multiomics/src/mcp_multiomics/tools/preprocessing.py
```

**Example result**:
```json
{
  "validation_status": "WARNING",
  "sample_overlap": {"common": 15},
  "missing_patterns": {"rna": "0.5%", "protein": "23.4%"},
  "batch_effects": {"protein_batch_correlation": 0.68},
  "warnings": ["Protein batch effect detected (PC1 ~ Batch, R²=0.68)"],
  "recommendations": ["Apply ComBat batch correction", "Use KNN imputation"]
}
```

---

#### 2. preprocess_multiomics_data

Applies ComBat batch correction, KNN imputation, quantile normalization, outlier removal.

```python
@mcp.tool()
def preprocess_multiomics_data(rna_path: str, protein_path: str, metadata_path: str, batch_correction: bool = True) -> dict:
    """Preprocess multi-omics data with batch correction and imputation."""
    # ComBat batch correction → KNN imputation → Quantile normalization → Outlier removal
    # Full implementation: servers/mcp-multiomics/src/mcp_multiomics/tools/preprocessing.py
```

**Before/after batch correction**:
```
Before: PC1 explains 45% variance, R² with Batch = 0.68 (BAD)
After:  PC1 explains 38% variance, R² with Batch = 0.12 (GOOD)
```

---

#### 3. visualize_data_quality

Generates PCA plots, missing value heatmaps, correlation matrices to verify preprocessing.

---

### Phase 2: Integration and Analysis

#### 4. integrate_omics_data

Aligns samples across modalities and normalizes within each.

```python
@mcp.tool()
def integrate_omics_data(rna_path: str, protein_path: str, phospho_path: str) -> dict:
    """Align samples across modalities, apply Z-score normalization per modality."""
    # Find common samples (15 in PatientOne), apply Z-score normalization
    # Filter features with >50% missing
    # Full implementation: servers/mcp-multiomics/src/mcp_multiomics/tools/integration.py
```

---

#### 5. run_halla_analysis

**HAllA (HMP All-Against-All Association)**: Tests pairwise associations between features.

**Scalability problem**: 1850 RNA features × 1850 = 3.4 million pairwise tests → days of compute.

**Solution: Chunking strategy**:
```python
# Chunk features into groups of 1000
chunks = [features[i:i+1000] for i in range(0, len(features), 1000)]
# Run HAllA on each chunk (1000×1000 = 1M tests, ~5 minutes)
for chunk in chunks:
    halla_results = halla.run(chunk, chunk, method='spearman')
```

**Result**: 3.4 million tests in ~15 minutes instead of days.

Implementation: [`servers/mcp-multiomics/src/mcp_multiomics/tools/halla.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-multiomics/src/mcp_multiomics/tools/halla.py)

---

#### 6. calculate_stouffer_meta

**Stouffer's Meta-Analysis**: Combines p-values across modalities to find consistently dysregulated features.

**Critical workflow** (per bioinformatician review):

```
CORRECT FDR TIMING FOR MULTI-OMICS META-ANALYSIS:
Step 1: HAllA Analysis → Use NOMINAL p-values (NOT FDR-corrected)
Step 2: Stouffer's Meta-Analysis → Combine NOMINAL p-values across modalities
Step 3: FDR Correction → Apply Benjamini-Hochberg AFTER combination

WHY: Applying FDR before Stouffer's loses power (over-conservative)
```

**Stouffer's Z-score method**:
```python
class StoufferMetaAnalysis:
    def p_to_z(self, p_values: np.ndarray, effect_sizes: np.ndarray = None) -> np.ndarray:
        """Convert p-values to Z-scores."""
        z_scores = stats.norm.ppf(1 - p_values / 2)
        if effect_sizes is not None: z_scores *= np.sign(effect_sizes)
        return z_scores

    def combine_z_scores(self, z_scores_list: List[np.ndarray], weights: List[float]) -> np.ndarray:
        """Combine Z-scores from multiple modalities."""
        z_meta = sum(w * z for w, z in zip(weights, z_scores_list))
        return z_meta / np.sqrt(sum(w**2 for w in weights))

    def apply_fdr_correction(self, p_values: np.ndarray) -> np.ndarray:
        """Apply Benjamini-Hochberg FDR correction AFTER meta-analysis."""
        _, q_values, _, _ = multipletests(p_values, method='fdr_bh')
        return q_values
```

**Example result** (AKT pathway):
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
    "combined_z_score": 3.58
  }
}
```

AKT1 shows consistent upregulation across all 3 modalities (meta q < 0.05), making it a high-confidence therapeutic target.

Implementation: [`servers/mcp-multiomics/src/mcp_multiomics/tools/stouffer.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-multiomics/src/mcp_multiomics/tools/stouffer.py)

---

#### 7. predict_upstream_regulators

Predicts kinases, transcription factors, and drug responses.

```python
@mcp.tool()
def predict_upstream_regulators(top_features: list) -> dict:
    """Predict upstream kinases, TFs, and druggable targets."""
    # Use kinase-substrate databases, TF binding motifs, drug target databases
    # Full implementation: servers/mcp-multiomics/src/mcp_multiomics/tools/upstream.py
```

**Example** (PatientOne AKT pathway):
```json
{
  "upstream_kinases": [
    {"kinase": "PDK1", "confidence": 0.92, "targets": ["AKT1", "AKT2"]},
    {"kinase": "mTORC2", "confidence": 0.87, "targets": ["AKT1"]}
  ],
  "drug_targets": [
    {"drug": "Capivasertib", "target": "AKT1/2/3", "clinical_phase": "Phase 3"}
  ]
}
```

This tells you: AKT hyperphosphorylated → likely by PDK1/mTORC2 → druggable with AKT inhibitors.

---

## The Complete PatientOne Workflow

Natural language prompt:
```
I have multi-omics data from 15 ovarian cancer PDX models. Please:
1. Validate data quality and check for batch effects
2. Preprocess with ComBat batch correction and KNN imputation
3. Integrate the three modalities
4. Run Stouffer meta-analysis to find consistently dysregulated pathways
5. Predict upstream regulators for the top 10 hits
```

**Results for PatientOne**:
- **15 common samples** across all 3 modalities
- **1850 RNA features**, 1243 proteins, 987 phosphosites
- **23 pathways** with meta q < 0.05 (AKT/mTOR, PI3K, TP53, DNA damage response)
- **12 upstream kinases** predicted
- **8 druggable targets** identified

---

## Implementation Walkthrough

### Project Setup

```bash
cd servers/mcp-multiomics
python -m venv venv && source venv/bin/activate
pip install fastmcp pandas numpy scipy statsmodels scikit-learn
```

Environment variables (`.env`):
```bash
MULTIOMICS_DRY_RUN=true  # For testing
MULTIOMICS_DATA_DIR=/workspace/data
```

### Stouffer Meta-Analysis Core

```python
@mcp.tool()
def calculate_stouffer_meta(rna_results: dict, protein_results: dict, phospho_results: dict, fdr_threshold: float = 0.05) -> dict:
    """Combine p-values across modalities using Stouffer's method."""
    # Extract NOMINAL p-values, convert to Z-scores, combine with weights
    # Apply FDR correction AFTER combining (CRITICAL)
    # Full implementation: servers/mcp-multiomics/src/mcp_multiomics/tools/stouffer.py
```

---

## Testing Your Server

```python
def test_stouffer_combine_strong_signals():
    """Test combining strong signals from all 3 modalities."""
    stouffer = StoufferMetaAnalysis()
    p_rna = np.array([0.005, 0.008, 0.003])
    p_protein = np.array([0.007, 0.004, 0.009])
    p_phospho = np.array([0.002, 0.006, 0.001])
    # Meta p-values should be even stronger
    assert np.all(p_meta < 0.001)
```

Test coverage: **67%**, 18 unit tests

---

## Connecting to Claude Desktop

```json
{
  "mcpServers": {
    "multiomics": {
      "command": "/path/to/venv/bin/python",
      "args": ["-m", "mcp_multiomics"],
      "env": {"MULTIOMICS_DRY_RUN": "false"}
    }
  }
}
```

---

## What You've Built

A multi-omics integration server providing:
1. **Data validation**: Batch effects, missing values, outliers
2. **Preprocessing**: ComBat, KNN imputation, normalization
3. **Integration**: Sample alignment, Z-score normalization
4. **Association testing**: HAllA with chunking (3.4M tests in 15 minutes)
5. **Meta-analysis**: Stouffer with correct FDR timing
6. **Regulator prediction**: Kinases, TFs, drug targets
7. **Visualization**: PCA, heatmaps, correlation matrices

This bridges genomics (Chapter 5) to spatial transcriptomics (Chapter 7) by identifying which pathways to investigate in tissue context.

---

## Try It Yourself

```bash
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp/servers/mcp-multiomics
python -m venv venv && source venv/bin/activate
pip install -e ".[dev]"
export MULTIOMICS_DRY_RUN=true && python -m mcp_multiomics
```

---

## Next Steps

In **Chapter 7: Spatial Transcriptomics**, you'll build `mcp-spatialtools` to analyze 10X Visium data and map the pathways you found here (AKT/mTOR, PI3K) to specific tissue regions.

---

**Chapter 6 Summary**:
- Multi-omics integration requires preprocessing (batch correction, imputation)
- HAllA with chunking enables scalable all-vs-all association testing
- Stouffer meta-analysis combines p-values BEFORE FDR correction (critical timing)
- Upstream regulator prediction connects pathways to druggable targets
- PatientOne: 23 pathways with meta q < 0.05, 8 druggable targets

**Files created**: `servers/mcp-multiomics/src/mcp_multiomics/server.py`, `tools/stouffer.py`, `tools/preprocessing.py`, `tools/halla.py`
**Tests added**: 18 unit tests, 67% coverage
**Tools exposed**: 7 MCP tools (validate, preprocess, visualize, integrate, halla, stouffer, upstream_regulators)
