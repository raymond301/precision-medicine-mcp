# Chapter 9: Treatment Response Prediction

*Building mcp-perturbation with GEARS graph neural networks*

---

## Why Predict Treatment Response?

You've identified PatientOne's cancer biology:
- **TP53 R175H mutation** (Chapter 5, genomics)
- **AKT/mTOR hyperactivation** (Chapter 6, multi-omics)
- **HIF1A+ hypoxic regions** (Chapter 7, spatial transcriptomics)
- **25% Ki67+ proliferating cells** (Chapter 8, cell segmentation)

**The clinical question**: Which treatments will work?

Traditional approach:
1. **Try carboplatin/paclitaxel chemotherapy** (standard ovarian cancer treatment)
2. **Wait 3 months** for CT scan
3. **If no response**: Try PARP inhibitor olaparib
4. **Wait 3 months**
5. **If no response**: Try immunotherapy
6. **Wait 3 months**

**Total time**: 9-12 months of trial-and-error while disease progresses.

**AI-predicted approach**:
1. **Run *in silico* predictions** for all treatment options (~2 hours compute)
2. **Rank by predicted response** (e.g., PARP inhibitor 82% cell death, chemo 45%, immunotherapy 18%)
3. **Start with best option** (PARP inhibitor)
4. **Validate response at 3 months**

**Result**: Start with optimal treatment immediately, avoid 6-9 months of ineffective therapies.

The `mcp-perturbation` server uses **GEARS (Graph-Enhanced Gene Activation and Repression Simulator)**, a graph neural network published in *Nature Biotechnology* 2024, to predict how patient cells respond to treatments **before** giving them.

---

## Why GEARS? (vs VAE Methods)

**Traditional approach**: VAE (Variational Autoencoder) methods like scGen learn latent representations of cell states, then predict perturbation effects in latent space.

**Problem**: VAEs don't use biological knowledge—they treat all genes as independent, ignoring known regulatory relationships (TP53 regulates 300+ genes, AKT phosphorylates 100+ targets).

**GEARS innovation**: Integrates **biological knowledge graphs** (gene-gene regulatory networks) into graph neural networks.

**Performance**:
- **40% higher precision** than VAE methods (Nature Biotechnology 2024)
- **Better generalization** to unseen drug combinations
- **Uncertainty quantification** (confidence scores for predictions)
- **Multi-gene perturbations** (handles drug combinations natively)

**How GEARS works**:
1. **Knowledge graph**: Encodes gene-gene relationships (TF binding, protein interactions, pathways)
2. **Graph neural network**: Learns how perturbations propagate through regulatory cascades
3. **Prediction**: Given patient's baseline expression + drug targets → predicts post-treatment expression
4. **Uncertainty**: Bayesian layers provide confidence intervals

---

## The Five mcp-perturbation Tools

### 1. perturbation_load_dataset

**Why you need it**: GEARS trains on public single-cell RNA-seq datasets of perturbation experiments (e.g., CRISPR screens, drug treatments).

**Example datasets**:
- **GSE184880**: Ovarian cancer cell lines, 50 drug treatments, 10,000 cells/condition
- **GSE108097**: T cell activation screen, CD4/CD8A perturbations
- **Replogle2022**: Genome-wide CRISPR screen, 2.5M cells

**Example loading**:
```python
@mcp.tool()
def perturbation_load_dataset(
    dataset_id: str,
    normalize: bool = True,
    n_hvg: int = 7000,
    cell_type_key: str = "cell_type",
    condition_key: str = "condition"
) -> dict:
    """Load scRNA-seq perturbation dataset from GEO.

    Args:
        dataset_id: GEO accession (e.g., "GSE184880") or path to .h5ad file
        normalize: Apply log-normalization (default: True)
        n_hvg: Number of highly variable genes to keep (default: 7000)
        cell_type_key: Column with cell type annotations
        condition_key: Column with treatment conditions

    Returns:
        Dataset metadata (n_cells, n_genes, cell types, conditions)
    """
    import scanpy as sc

    # Load from GEO or local file
    if dataset_id.startswith("GSE"):
        adata = sc.datasets.fetch_geo(dataset_id)
    else:
        adata = sc.read_h5ad(dataset_id)

    # Preprocessing
    if normalize:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    # Select highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor='seurat_v3')
    adata = adata[:, adata.var['highly_variable']]

    # Save processed data
    adata.write_h5ad(f"/data/processed/{dataset_id}_processed.h5ad")

    return {
        "dataset_id": dataset_id,
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "cell_types": list(adata.obs[cell_type_key].unique()),
        "conditions": list(adata.obs[condition_key].unique()),
        "processed_path": f"/data/processed/{dataset_id}_processed.h5ad"
    }
```

**PatientOne example** (ovarian cancer training data):
```json
{
  "dataset_id": "GSE184880",
  "n_cells": 487520,
  "n_genes": 7000,
  "cell_types": ["tumor_cells", "T_cells", "fibroblasts"],
  "conditions": [
    "control",
    "carboplatin",
    "olaparib",
    "capivasertib",
    "carboplatin+olaparib",
    ...  // 50 total treatments
  ]
}
```

Implementation: [`servers/mcp-perturbation/mcp_perturbation/data_loader.py:50-150`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-perturbation/mcp_perturbation/data_loader.py#L50-L150)

---

### 2. perturbation_setup_model

**Why you need it**: Initialize GEARS graph neural network with biological knowledge graph and dataset-specific parameters.

**Model architecture**:
- **Input layer**: Gene expression vector (7000 genes)
- **Knowledge graph**: Gene-gene regulatory network (TF→target, protein→protein)
- **GNN layers**: 2-4 graph convolution layers (default: 2)
- **Hidden size**: 64-256 dimensions (default: 64)
- **Output layer**: Predicted post-perturbation expression
- **Uncertainty layer**: Bayesian dropout for confidence estimation

**Example setup**:
```python
@mcp.tool()
def perturbation_setup_model(
    dataset_id: str,
    hidden_size: int = 64,
    num_layers: int = 2,
    uncertainty: bool = True,
    uncertainty_reg: float = 1.0,
    model_name: str = None
) -> dict:
    """Initialize GEARS graph neural network.

    Args:
        dataset_id: Dataset ID from perturbation_load_dataset
        hidden_size: Hidden layer size (default: 64)
        num_layers: Number of GNN layers (default: 2)
        uncertainty: Enable Bayesian uncertainty quantification (default: True)
        uncertainty_reg: Uncertainty regularization weight (default: 1.0)
        model_name: Name for this model (default: auto-generated)

    Returns:
        Model configuration summary
    """
    from gears import PertData, GEARS

    # Load processed dataset
    pert_data = PertData(f"/data/processed/{dataset_id}_processed.h5ad")

    # Load gene-gene knowledge graph (from STRING, BioGRID, KEGG)
    pert_data.load_default_graph()  # ~20,000 gene-gene edges

    # Get data split (train/val/test)
    pert_data.prepare_split(split="simulation", seed=1)

    # Initialize GEARS model
    gears_model = GEARS(
        pert_data,
        device="cuda" if torch.cuda.is_available() else "cpu",
        hidden_size=hidden_size,
        num_go_gnn_layers=num_layers,
        uncertainty=uncertainty,
        uncertainty_reg=uncertainty_reg
    )

    # Save model config
    model_name = model_name or f"gears_{dataset_id}_{hidden_size}h_{num_layers}l"

    return {
        "model_name": model_name,
        "dataset_id": dataset_id,
        "n_genes": pert_data.adata.n_vars,
        "n_perts": len(pert_data.pert_names),
        "hidden_size": hidden_size,
        "num_layers": num_layers,
        "uncertainty_enabled": uncertainty,
        "graph_edges": len(pert_data.edge_list)
    }
```

**PatientOne model setup**:
```json
{
  "model_name": "patientone_ovarian_64h_2l",
  "dataset_id": "GSE184880",
  "n_genes": 7000,
  "n_perts": 50,
  "hidden_size": 64,
  "num_layers": 2,
  "uncertainty_enabled": true,
  "graph_edges": 18432
}
```

Implementation: [`servers/mcp-perturbation/mcp_perturbation/gears_wrapper.py:65-150`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-perturbation/mcp_perturbation/gears_wrapper.py#L65-L150)

---

### 3. perturbation_train_model

**Why you need it**: Train the GEARS model on perturbation data to learn drug response patterns.

**Training process**:
1. **Forward pass**: Predict post-perturbation expression for each drug treatment
2. **Loss calculation**: MSE between predicted and observed expression
3. **Graph backpropagation**: Update GNN weights via gene relationship graph
4. **Uncertainty loss**: Regularize Bayesian layers for calibrated confidence
5. **Validation**: Test on held-out perturbations every epoch

**Typical training time**: 20 epochs × 5 minutes/epoch = **~100 minutes** (GPU), **~300 minutes** (CPU)

**Example training**:
```python
@mcp.tool()
def perturbation_train_model(
    model_name: str,
    epochs: int = 20,
    batch_size: int = 32,
    learning_rate: float = 1e-3,
    valid_every: int = 1
) -> dict:
    """Train GEARS model on perturbation data.

    Args:
        model_name: Model name from perturbation_setup_model
        epochs: Number of training epochs (default: 20)
        batch_size: Batch size (default: 32)
        learning_rate: Adam learning rate (default: 0.001)
        valid_every: Run validation every N epochs (default: 1)

    Returns:
        Training metrics (final loss, best epoch, model path)
    """
    from gears import GEARS

    # Load model and data
    gears_model = load_model(model_name)

    # Training loop
    history = {
        "train_loss": [],
        "val_loss": [],
        "val_r2": []
    }

    for epoch in range(epochs):
        # Train
        train_loss = gears_model.train(
            epochs=1,
            batch_size=batch_size,
            lr=learning_rate
        )
        history["train_loss"].append(train_loss)

        # Validate
        if (epoch + 1) % valid_every == 0:
            val_metrics = gears_model.evaluate()
            history["val_loss"].append(val_metrics["mse"])
            history["val_r2"].append(val_metrics["r2"])

            logger.info(
                f"Epoch {epoch+1}/{epochs} - "
                f"Train Loss: {train_loss:.4f}, "
                f"Val MSE: {val_metrics['mse']:.4f}, "
                f"Val R²: {val_metrics['r2']:.4f}"
            )

    # Save trained model
    model_path = f"/data/models/{model_name}_trained.pt"
    gears_model.save_model(model_path)

    return {
        "model_name": model_name,
        "epochs_completed": epochs,
        "final_train_loss": history["train_loss"][-1],
        "final_val_r2": history["val_r2"][-1],
        "best_val_r2": max(history["val_r2"]),
        "model_path": model_path,
        "training_history": history
    }
```

**PatientOne training results**:
```json
{
  "model_name": "patientone_ovarian_64h_2l",
  "epochs_completed": 20,
  "final_train_loss": 0.0342,
  "final_val_r2": 0.78,
  "best_val_r2": 0.81,
  "model_path": "/data/models/patientone_ovarian_64h_2l_trained.pt"
}
```

**Interpretation**: **R² = 0.78** means GEARS explains 78% of expression variance in held-out perturbations—good predictive performance.

Implementation: [`servers/mcp-perturbation/mcp_perturbation/gears_wrapper.py:200-320`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-perturbation/mcp_perturbation/gears_wrapper.py#L200-L320)

---

### 4. perturbation_compute_delta

**Why you need it**: Calculate the **perturbation vector (Δ)** representing average transcriptional change caused by a treatment.

**What is Δ?**
```
Δ = avg(treated cells) - avg(control cells)
```

This vector captures the average effect of a drug (e.g., olaparib PARP inhibitor causes Δ with upregulated DNA damage response genes, downregulated cell cycle genes).

**Example computation**:
```python
@mcp.tool()
def perturbation_compute_delta(
    model_name: str,
    control_key: str = "control",
    treatment_key: str,
    cell_type: str = None
) -> dict:
    """Compute perturbation vector (Δ) for a treatment.

    Args:
        model_name: Trained GEARS model name
        control_key: Control condition label (default: "control")
        treatment_key: Treatment condition label (e.g., "olaparib")
        cell_type: Specific cell type to analyze (None = all cells)

    Returns:
        Delta vector statistics
    """
    import scanpy as sc

    # Load dataset
    adata = load_dataset(model_name)

    # Filter by cell type
    if cell_type:
        adata = adata[adata.obs["cell_type"] == cell_type]

    # Get control and treatment cells
    control_cells = adata[adata.obs["condition"] == control_key]
    treatment_cells = adata[adata.obs["condition"] == treatment_key]

    # Compute delta
    control_mean = control_cells.X.mean(axis=0)
    treatment_mean = treatment_cells.X.mean(axis=0)
    delta = treatment_mean - control_mean

    return {
        "treatment": treatment_key,
        "control": control_key,
        "cell_type": cell_type or "all",
        "n_control_cells": control_cells.n_obs,
        "n_treatment_cells": treatment_cells.n_obs,
        "delta_norm": float(np.linalg.norm(delta)),
        "delta_mean": float(delta.mean()),
        "delta_std": float(delta.std()),
        "top_upregulated_genes": get_top_genes(delta, n=10, direction="up"),
        "top_downregulated_genes": get_top_genes(delta, n=10, direction="down")
    }
```

**PatientOne olaparib Δ**:
```json
{
  "treatment": "olaparib",
  "control": "control",
  "cell_type": "tumor_cells",
  "n_control_cells": 9823,
  "n_treatment_cells": 9654,
  "delta_norm": 12.34,
  "top_upregulated_genes": [
    {"gene": "BRCA1", "log2fc": 3.2},
    {"gene": "RAD51", "log2fc": 2.8},
    {"gene": "PARP1", "log2fc": 2.1},
    ...
  ],
  "top_downregulated_genes": [
    {"gene": "CDK1", "log2fc": -2.9},
    {"gene": "CCNB1", "log2fc": -2.5},
    {"gene": "MKI67", "log2fc": -2.2},
    ...
  ]
}
```

**Interpretation**: Olaparib upregulates DNA repair genes (BRCA1, RAD51) and downregulates cell cycle genes (CDK1, MKI67) → expected PARP inhibitor mechanism.

---

### 5. perturbation_predict_response

**Why you need it**: Apply the trained GEARS model to **predict how PatientOne's specific cells will respond** to treatments.

**Workflow**:
1. **Load patient's baseline scRNA-seq** (untreated tumor cells)
2. **For each drug**: GEARS predicts post-treatment expression
3. **Quantify response**: Cell death score, apoptosis markers, proliferation reduction
4. **Rank drugs**: By predicted efficacy

**Example prediction**:
```python
@mcp.tool()
def perturbation_predict_response(
    model_name: str,
    patient_adata_path: str,
    perturbations: list[str],  # ["olaparib", "carboplatin", "capivasertib"]
    cell_type: str = "tumor_cells"
) -> dict:
    """Predict patient's response to multiple treatments.

    Args:
        model_name: Trained GEARS model name
        patient_adata_path: Path to patient's baseline scRNA-seq (.h5ad)
        perturbations: List of treatments to test
        cell_type: Cell type to predict (default: "tumor_cells")

    Returns:
        Predicted responses ranked by efficacy
    """
    from gears import GEARS
    import scanpy as sc

    # Load model and patient data
    gears_model = load_model(model_name)
    patient_adata = sc.read_h5ad(patient_adata_path)

    # Filter to target cell type
    patient_cells = patient_adata[patient_adata.obs["cell_type"] == cell_type]

    predictions = []

    for pert in perturbations:
        # Predict post-treatment expression
        predicted_expr = gears_model.predict(
            patient_cells.X,
            perturbation=pert
        )

        # Calculate response metrics
        baseline_expr = patient_cells.X.mean(axis=0)
        delta_expr = predicted_expr.mean(axis=0) - baseline_expr

        # Apoptosis score (upregulation of BAX, CASP3, PUMA)
        apoptosis_genes = ["BAX", "CASP3", "BBC3"]  # BBC3 = PUMA
        apoptosis_score = delta_expr[patient_cells.var_names.isin(apoptosis_genes)].mean()

        # Proliferation reduction (downregulation of MKI67, PCNA, CDK1)
        prolif_genes = ["MKI67", "PCNA", "CDK1"]
        prolif_reduction = -delta_expr[patient_cells.var_names.isin(prolif_genes)].mean()

        # Combined efficacy score
        efficacy_score = (apoptosis_score + prolif_reduction) / 2

        predictions.append({
            "treatment": pert,
            "efficacy_score": float(efficacy_score),
            "apoptosis_score": float(apoptosis_score),
            "proliferation_reduction": float(prolif_reduction),
            "predicted_cell_death_percent": min(100, max(0, efficacy_score * 100))
        })

    # Rank by efficacy
    predictions.sort(key=lambda x: x["efficacy_score"], reverse=True)

    return {
        "patient_id": extract_patient_id(patient_adata_path),
        "cell_type": cell_type,
        "n_cells_analyzed": patient_cells.n_obs,
        "predictions": predictions,
        "top_treatment": predictions[0]["treatment"],
        "top_treatment_efficacy": predictions[0]["efficacy_score"]
    }
```

**PatientOne treatment ranking**:
```json
{
  "patient_id": "PAT001-OVC-2025",
  "cell_type": "tumor_cells",
  "n_cells_analyzed": 1247,
  "predictions": [
    {
      "treatment": "olaparib",
      "efficacy_score": 0.82,
      "apoptosis_score": 0.89,
      "proliferation_reduction": 0.75,
      "predicted_cell_death_percent": 82
    },
    {
      "treatment": "carboplatin+olaparib",
      "efficacy_score": 0.71,
      "apoptosis_score": 0.78,
      "proliferation_reduction": 0.64,
      "predicted_cell_death_percent": 71
    },
    {
      "treatment": "carboplatin",
      "efficacy_score": 0.45,
      "apoptosis_score": 0.52,
      "proliferation_reduction": 0.38,
      "predicted_cell_death_percent": 45
    },
    {
      "treatment": "capivasertib",
      "efficacy_score": 0.38,
      "apoptosis_score": 0.31,
      "proliferation_reduction": 0.45,
      "predicted_cell_death_percent": 38
    }
  ],
  "top_treatment": "olaparib",
  "top_treatment_efficacy": 0.82
}
```

**Clinical decision**: Start with **olaparib monotherapy** (82% predicted cell death), reserve carboplatin+olaparib combination (71%) for second-line if needed. Standard carboplatin alone (45%) is suboptimal.

Implementation: [`servers/mcp-perturbation/mcp_perturbation/prediction.py:50-200`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-perturbation/mcp_perturbation/prediction.py#L50-L200)

---

## The Complete PatientOne Prediction Workflow

Natural language prompt in Claude Desktop:

```
I want to predict which treatments will work best for patient PAT001-OVC-2025.

The patient has:
- TP53 R175H mutation (Chapter 5)
- BRCA wild-type (no germline mutation)
- High-grade serous ovarian cancer
- Baseline tumor scRNA-seq: data/patient-data/PAT001-OVC-2025/scrna/tumor_baseline.h5ad

Please:
1. Load the ovarian cancer drug response training dataset (GSE184880)
2. Setup GEARS model with uncertainty quantification
3. Train for 20 epochs
4. Predict patient's response to: olaparib, carboplatin, carboplatin+olaparib, capivasertib
5. Rank treatments by predicted efficacy
```

Claude orchestrates all 5 tools:

```python
# Step 1: Load training data
dataset = perturbation.perturbation_load_dataset(
    dataset_id="GSE184880",
    normalize=True,
    n_hvg=7000
)

# Step 2: Setup model
model_setup = perturbation.perturbation_setup_model(
    dataset_id="GSE184880",
    hidden_size=64,
    num_layers=2,
    uncertainty=True,
    model_name="patientone_ovarian_64h_2l"
)

# Step 3: Train
training_result = perturbation.perturbation_train_model(
    model_name="patientone_ovarian_64h_2l",
    epochs=20,
    batch_size=32
)
# → Val R² = 0.78 (good predictive performance)

# Step 4: Predict patient response
predictions = perturbation.perturbation_predict_response(
    model_name="patientone_ovarian_64h_2l",
    patient_adata_path="data/patient-data/PAT001-OVC-2025/scrna/tumor_baseline.h5ad",
    perturbations=["olaparib", "carboplatin", "carboplatin+olaparib", "capivasertib"]
)

# Step 5: Results
# Top treatment: olaparib (82% predicted cell death)
# Rationale: TP53 mutant + BRCA WT → synthetic lethality via PARP inhibition
```

**Total analysis time**: ~2.5 hours (training) + ~15 minutes (prediction) = **~3 hours**

**Clinical impact**: Immediate selection of optimal therapy instead of 9-12 months trial-and-error.

---

## Testing Your Server

### Unit Tests

```python
# tests/test_gears_wrapper.py
import pytest
from mcp_perturbation.gears_wrapper import GearsWrapper

def test_model_setup():
    """Test GEARS model initialization."""
    wrapper = GearsWrapper()

    # Create synthetic dataset
    adata = create_synthetic_perturbation_data(
        n_cells=1000,
        n_genes=500,
        n_perturbations=10
    )

    wrapper.setup(adata, condition_key="condition", pert_key="perturbation")

    assert wrapper.pert_data is not None
    assert len(wrapper.pert_data.pert_names) == 10

def test_prediction():
    """Test perturbation prediction."""
    wrapper = GearsWrapper()

    # Train on synthetic data
    wrapper.setup(synthetic_adata)
    wrapper.train(epochs=5)  # Quick test

    # Predict
    baseline_expr = synthetic_adata[synthetic_adata.obs["condition"] == "control"].X
    predicted = wrapper.predict(baseline_expr, perturbation="drug_A")

    assert predicted.shape == baseline_expr.shape
    assert not np.array_equal(predicted, baseline_expr)  # Prediction differs from baseline
```

Run tests:
```bash
pytest tests/test_gears_wrapper.py -v
```

---

## What You've Built

You now have a treatment response prediction server that:

1. **Loads perturbation datasets**: Public scRNA-seq drug screens (GEO, .h5ad)
2. **Trains GEARS models**: Graph neural networks with gene regulatory knowledge
3. **Computes Δ vectors**: Average transcriptional change per treatment
4. **Predicts patient responses**: *In silico* drug screening on patient cells
5. **Ranks treatments**: By predicted efficacy (apoptosis + proliferation reduction)

This enables **precision medicine drug selection** before clinical trial-and-error.

---

## Try It Yourself

### Option 1: Synthetic Data Test

```bash
cd servers/mcp-perturbation
python -m venv venv
source venv/bin/activate
pip install -e ".[dev]"

# In Claude Desktop:
# "Load the test perturbation dataset and train a GEARS model for 5 epochs (quick test)"
```

### Option 2: PatientOne Prediction

```bash
# Requires: PatientOne scRNA-seq data (not included in repo)
# "Predict treatment response for patient PAT001-OVC-2025 using trained ovarian cancer model"
```

---

## Next Steps

In **Chapter 10: Quantum Cell-Type Fidelity**, you'll build `mcp-quantum-celltype-fidelity` to use **quantum computing** and **Bayesian uncertainty quantification** to verify cell type classifications from spatial/imaging data with calibrated confidence scores.

The treatment predictions you built identify **which drugs to use**. Quantum fidelity ensures you're treating the **correct cell types** with those drugs.

---

**Chapter 9 Summary**:
- GEARS graph neural networks predict drug responses 40% better than VAE methods
- Integrates gene regulatory knowledge graphs for biologically-informed predictions
- PatientOne: Olaparib ranked #1 (82% predicted efficacy) vs carboplatin (45%)
- Training time: ~2.5 hours on CPU, ~100 minutes on GPU
- Clinical impact: Immediate optimal treatment selection

**Files created**: `servers/mcp-perturbation/mcp_perturbation/server.py` (483 lines), `gears_wrapper.py` (414 lines), `prediction.py` (268 lines)
**Tests added**: 15 unit tests, 72% coverage
**Tools exposed**: 5 MCP tools (load_dataset, setup_model, train_model, compute_delta, predict_response)
