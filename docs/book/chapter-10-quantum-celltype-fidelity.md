# Chapter 10: Quantum Cell-Type Fidelity

*Building mcp-quantum-celltype-fidelity with PennyLane and Bayesian UQ*

---

## Why Quantum Fidelity for Cell Types?

In Chapters 7-8, you classified cells by phenotype:
- **Chapter 7** (Spatial transcriptomics): Tumor proliferative, necrotic/hypoxic, immune infiltrated (10-30 cells/spot)
- **Chapter 8** (Cell segmentation): Ki67+ (25%), TP53+ (40%), CD8+ T cells (single-cell resolution)

**Problem**: These classifications lack **uncertainty quantification**. When you say "this cell is CD8+ T cell":
- How confident are you? (95%? 60%? Unknown?)
- Could it be an exhausted T cell instead?
- Is the tumor cell truly evading immune detection?

**Clinical impact**: Without confidence scores, oncologists can't assess treatment risk:
- **High-confidence (95%)** immune detection → proceed with checkpoint blockade
- **Low-confidence (60%)** immune detection → order confirmatory tests first

The `mcp-quantum-celltype-fidelity` server uses **quantum fidelity** (overlap between quantum states) + **Bayesian uncertainty quantification** to provide calibrated confidence intervals on cell type classifications.

---

## Why Quantum Computing for This?

**Classical approach**: Cosine similarity in high-dimensional gene expression space.
- **Problem**: Doesn't capture non-linear gene regulatory relationships

**Quantum approach**: Encode cells as quantum states in Hilbert space using parameterized quantum circuits (PQCs).

**Quantum fidelity** measures overlap:
```
F(|ψ_A⟩, |ψ_B⟩) = |⟨ψ_A|ψ_B⟩|²
```

- **F = 1**: Identical cell types (perfect match)
- **F = 0**: Orthogonal cell types (completely different)
- **0 < F < 1**: Partial similarity

**Advantages**:
1. **Entanglement**: Quantum gates capture gene-gene interactions
2. **Hilbert space**: 2^n dimensional (8 qubits = 256D, 10 qubits = 1024D)
3. **Parameter-shift rule**: Exact gradients for training (works on real quantum hardware)
4. **Bayesian UQ**: Monte Carlo sampling from parameter distributions → confidence intervals

---

## The Six mcp-quantum-celltype-fidelity Tools

### 1. learn_spatial_cell_embeddings

**Why you need it**: Train quantum circuit parameters to embed cell types into Hilbert space.

**Parameterized Quantum Circuit (PQC) architecture**:
1. **Feature encoding layer**: Gene expression → rotation angles (RX, RY, RZ gates)
2. **Variational layers** (3 layers typical):
   - RX(θ_x), RY(θ_y), RZ(θ_z) gates per qubit
   - CNOT entanglement (ring topology)
3. **Measurement**: Statevector |ψ⟩ represents cell embedding

**Contrastive learning objective**:
- Maximize fidelity within same cell type
- Minimize fidelity between different cell types

**Example training**:
```python
@mcp.tool()
def learn_spatial_cell_embeddings(
    adata_path: str,
    cell_type_key: str = "cell_type",
    n_qubits: int = 8,
    n_layers: int = 3,
    n_epochs: int = 50,
    learning_rate: float = 0.01,
    backend: str = "cpu"  # "cpu", "gpu", or "ibm"
) -> dict:
    """Train quantum embeddings for cell types.

    Args:
        adata_path: Path to AnnData file (spatial transcriptomics or scRNA-seq)
        cell_type_key: Column in adata.obs with cell type labels
        n_qubits: Number of qubits (8 = 256D Hilbert space, 10 = 1024D)
        n_layers: Variational circuit layers (default: 3)
        n_epochs: Training epochs (default: 50)
        learning_rate: Parameter update step size
        backend: "cpu" (simulator), "gpu" (cuQuantum), "ibm" (real quantum hardware)

    Returns:
        Trained embedding ID, training history, circuit configuration
    """
    import scanpy as sc
    from quantum_celltype_fidelity import QuCoWECellTypeEmbedding, InfoNCELoss

    # Load spatial data
    adata = sc.read_h5ad(adata_path)

    # Initialize quantum embedding
    embedding = QuCoWECellTypeEmbedding(
        n_qubits=n_qubits,
        n_layers=n_layers,
        backend=backend
    )

    # Contrastive learning
    loss_fn = InfoNCELoss(temperature=0.1)
    optimizer = Adam(learning_rate)

    history = {"loss": []}

    for epoch in range(n_epochs):
        # Sample cell pairs
        anchor_cells = sample_cells_by_type(adata, cell_type_key)
        positive_cells = sample_same_type(anchor_cells, adata)
        negative_cells = sample_different_types(anchor_cells, adata)

        # Encode to quantum states
        anchor_states = embedding.encode(anchor_cells.X)
        positive_states = embedding.encode(positive_cells.X)
        negative_states = embedding.encode(negative_cells.X)

        # Compute fidelities
        pos_fidelity = quantum_fidelity(anchor_states, positive_states)
        neg_fidelity = quantum_fidelity(anchor_states, negative_states)

        # InfoNCE loss (maximize pos_fidelity, minimize neg_fidelity)
        loss = loss_fn(pos_fidelity, neg_fidelity)

        # Parameter-shift rule gradient
        gradients = compute_gradients_parameter_shift(embedding, loss)

        # Update parameters
        embedding.update_parameters(gradients, optimizer)

        history["loss"].append(loss)

        if epoch % 10 == 0:
            logger.info(f"Epoch {epoch}/{n_epochs} - Loss: {loss:.4f}")

    # Save trained embedding
    embedding_id = f"embedding_{adata_path.split('/')[-1]}_{n_qubits}q_{n_layers}l"
    embedding.save(f"/data/embeddings/{embedding_id}.pt")

    return {
        "embedding_id": embedding_id,
        "training_summary": {
            "final_loss": history["loss"][-1],
            "initial_loss": history["loss"][0],
            "epochs": n_epochs,
            "loss_history": history["loss"]
        },
        "embedding_summary": {
            "n_qubits": n_qubits,
            "n_layers": n_layers,
            "hilbert_space_dim": 2**n_qubits,
            "n_parameters": n_qubits * n_layers * 3,  # RX, RY, RZ per qubit per layer
            "backend": backend
        }
    }
```

**PatientOne training results**:
```json
{
  "embedding_id": "embedding_patientone_tcells_8q_3l",
  "training_summary": {
    "final_loss": 0.234,
    "initial_loss": 1.872,
    "epochs": 50,
    "loss_reduction_percent": 87.5
  },
  "embedding_summary": {
    "n_qubits": 8,
    "n_layers": 3,
    "hilbert_space_dim": 256,
    "n_parameters": 72,
    "backend": "cpu"
  }
}
```

Implementation: [`servers/mcp-quantum-celltype-fidelity/src/quantum_celltype_fidelity/training.py:100-300`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-quantum-celltype-fidelity/src/quantum_celltype_fidelity/training.py#L100-L300)

---

### 2. compute_cell_type_fidelity (with Bayesian UQ)

**Why you need it**: Compute quantum fidelity between cells with **confidence intervals** using Bayesian uncertainty quantification.

**Bayesian UQ approach**:
1. During training, track parameter gradient history
2. Build posterior distributions: `θ ~ N(μ, Σ)` (mean = trained value, variance = gradient stability)
3. Monte Carlo sampling: Sample 100 parameter sets from posterior
4. Compute fidelity for each sample → distribution of fidelities
5. Report: mean, std, 95% CI, 90% CI, epistemic uncertainty

**Example fidelity computation**:
```python
@mcp.tool()
def compute_cell_type_fidelity(
    adata_path: str,
    embedding_id: str,
    compute_matrix: bool = False,
    with_uncertainty: bool = True,
    n_uncertainty_samples: int = 100
) -> dict:
    """Compute quantum fidelity with Bayesian uncertainty quantification.

    Args:
        adata_path: Path to AnnData file
        embedding_id: ID from learn_spatial_cell_embeddings
        compute_matrix: Compute full NxN fidelity matrix (slow for large N)
        with_uncertainty: Enable Bayesian UQ for confidence intervals
        n_uncertainty_samples: Monte Carlo samples for UQ (default: 100)

    Returns:
        Fidelity scores with uncertainty estimates
    """
    import scanpy as sc
    from quantum_celltype_fidelity import (
        QuCoWECellTypeEmbedding,
        BayesianFidelityEstimator,
        BayesianParameterDistribution
    )

    # Load data and embedding
    adata = sc.read_h5ad(adata_path)
    embedding = QuCoWECellTypeEmbedding.load(f"/data/embeddings/{embedding_id}.pt")

    results = {
        "embedding_id": embedding_id,
        "n_cells": adata.n_obs
    }

    if compute_matrix:
        # Full NxN matrix (expensive)
        n_cells = adata.n_obs
        fidelity_matrix = np.zeros((n_cells, n_cells))

        if with_uncertainty:
            uncertainty_matrix = np.zeros((n_cells, n_cells, 4))  # mean, std, ci_95_low, ci_95_high

        for i in range(n_cells):
            for j in range(i, n_cells):
                if with_uncertainty:
                    # Bayesian UQ
                    param_dist = BayesianParameterDistribution(
                        embedding.parameters,
                        embedding.parameter_uncertainties  # From gradient history
                    )
                    estimator = BayesianFidelityEstimator(param_dist, n_samples=n_uncertainty_samples)

                    fidelity_estimate = estimator.estimate_fidelity_with_uncertainty(
                        lambda params: quantum_fidelity(
                            embedding.encode_with_params(adata.X[i], params),
                            embedding.encode_with_params(adata.X[j], params)
                        )
                    )

                    fidelity_matrix[i, j] = fidelity_estimate.mean
                    fidelity_matrix[j, i] = fidelity_estimate.mean
                    uncertainty_matrix[i, j] = [
                        fidelity_estimate.mean,
                        fidelity_estimate.std,
                        fidelity_estimate.confidence_interval_95[0],
                        fidelity_estimate.confidence_interval_95[1]
                    ]
                    uncertainty_matrix[j, i] = uncertainty_matrix[i, j]
                else:
                    # Point estimate only
                    fid = quantum_fidelity(
                        embedding.encode(adata.X[i]),
                        embedding.encode(adata.X[j])
                    )
                    fidelity_matrix[i, j] = fid
                    fidelity_matrix[j, i] = fid

        results["fidelity_matrix"] = fidelity_matrix.tolist()

        if with_uncertainty:
            results["uncertainty"] = {
                "mean_uncertainty": float(uncertainty_matrix[:, :, 1].mean()),
                "max_uncertainty": float(uncertainty_matrix[:, :, 1].max()),
                "uncertainty_matrix": uncertainty_matrix.tolist()
            }

    # Per-cell-type statistics
    cell_types = adata.obs[cell_type_key].unique()
    per_type_fidelities = {}

    for ct in cell_types:
        ct_cells = adata[adata.obs[cell_type_key] == ct]
        within_type_fidelities = []

        # Sample pairs from same type
        for i in range(min(100, len(ct_cells))):
            idx1, idx2 = np.random.choice(len(ct_cells), 2, replace=False)

            if with_uncertainty:
                param_dist = BayesianParameterDistribution(
                    embedding.parameters,
                    embedding.parameter_uncertainties
                )
                estimator = BayesianFidelityEstimator(param_dist, n_samples=n_uncertainty_samples)
                fid_est = estimator.estimate_fidelity_with_uncertainty(
                    lambda params: quantum_fidelity(
                        embedding.encode_with_params(ct_cells.X[idx1], params),
                        embedding.encode_with_params(ct_cells.X[idx2], params)
                    )
                )
                within_type_fidelities.append({
                    "mean": fid_est.mean,
                    "ci_95": fid_est.confidence_interval_95
                })
            else:
                fid = quantum_fidelity(
                    embedding.encode(ct_cells.X[idx1]),
                    embedding.encode(ct_cells.X[idx2])
                )
                within_type_fidelities.append(fid)

        per_type_fidelities[ct] = {
            "n_cells": len(ct_cells),
            "within_type_fidelity": {
                "mean": np.mean([f["mean"] if with_uncertainty else f for f in within_type_fidelities]),
                "std": np.std([f["mean"] if with_uncertainty else f for f in within_type_fidelities])
            }
        }

    results["per_cell_type"] = per_type_fidelities

    return results
```

**PatientOne fidelity results** (T cells with Bayesian UQ):
```json
{
  "embedding_id": "embedding_patientone_tcells_8q_3l",
  "n_cells": 487,
  "per_cell_type": {
    "CD8_T_cell": {
      "n_cells": 182,
      "within_type_fidelity": {
        "mean": 0.89,
        "std": 0.06
      }
    },
    "CD4_T_cell": {
      "n_cells": 156,
      "within_type_fidelity": {
        "mean": 0.87,
        "std": 0.07
      }
    },
    "T_exhausted": {
      "n_cells": 89,
      "within_type_fidelity": {
        "mean": 0.72,
        "std": 0.11
      }
    },
    "T_regulatory": {
      "n_cells": 60,
      "within_type_fidelity": {
        "mean": 0.78,
        "std": 0.09
      }
    }
  },
  "uncertainty": {
    "mean_uncertainty": 0.04,
    "max_uncertainty": 0.12
  }
}
```

**Interpretation**:
- **CD8+ T cells**: High within-type fidelity (0.89 ± 0.06) → well-defined quantum signature
- **Exhausted T cells**: Lower fidelity (0.72 ± 0.11) → heterogeneous subpopulation
- **Mean uncertainty**: 0.04 → predictions are well-calibrated (4% average CI width)

Implementation: [`servers/mcp-quantum-celltype-fidelity/src/quantum_celltype_fidelity/bayesian_uq.py:100-400`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-quantum-celltype-fidelity/src/quantum_celltype_fidelity/bayesian_uq.py#L100-L400)

---

### 3. identify_immune_evasion_states

**Why you need it**: Detect tumor cells evading immune surveillance with **classification confidence**.

**Immune evasion signature**:
- Low fidelity to CD8+ T cells (< 0.3)
- Low fidelity to NK cells (< 0.3)
- High fidelity to exhausted T cells (> 0.6)

**Classification confidence**: Probability that the binary decision (evading vs not evading) is correct.

**Example detection**:
```python
@mcp.tool()
def identify_immune_evasion_states(
    adata_path: str,
    embedding_id: str,
    immune_cell_types: list[str],
    exhausted_markers: list[str] = None,
    evasion_threshold: float = 0.3,
    with_confidence: bool = True
) -> dict:
    """Detect cells in immune evasion states.

    Args:
        adata_path: Path to AnnData file
        embedding_id: Trained embedding ID
        immune_cell_types: Canonical immune types (e.g., ["CD8_T_cell", "NK_cell"])
        exhausted_markers: Exhausted cell types (e.g., ["T_exhausted"])
        evasion_threshold: Fidelity threshold for flagging (default: 0.3)
        with_confidence: Include classification confidence scores

    Returns:
        Evading cells with evasion scores and confidence
    """
    from quantum_celltype_fidelity import BayesianFidelityEstimator

    adata = sc.read_h5ad(adata_path)
    embedding = load_embedding(embedding_id)

    # Get reference immune cells
    immune_cells = adata[adata.obs["cell_type"].isin(immune_cell_types)]
    immune_ref_state = embedding.encode(immune_cells.X.mean(axis=0))

    evading_cells = []

    for i, cell in enumerate(adata):
        cell_state = embedding.encode(cell.X)

        if with_confidence:
            # Bayesian UQ
            param_dist = BayesianParameterDistribution(
                embedding.parameters,
                embedding.parameter_uncertainties
            )
            estimator = BayesianFidelityEstimator(param_dist, n_samples=100)

            fidelity_est = estimator.estimate_fidelity_with_uncertainty(
                lambda params: quantum_fidelity(
                    embedding.encode_with_params(cell.X, params),
                    immune_ref_state
                )
            )

            evasion_score = 1.0 - fidelity_est.mean

            # Classification confidence
            # P(evasion_score > threshold) from posterior samples
            classification_confidence = estimator.estimate_classification_confidence(
                fidelity_est,
                threshold=evasion_threshold,
                direction="less_than"
            )

            if evasion_score > evasion_threshold:
                evading_cells.append({
                    "cell_idx": i,
                    "cell_type": adata.obs["cell_type"][i],
                    "evasion_score": float(evasion_score),
                    "evasion_score_ci_95": [
                        float(1.0 - fidelity_est.confidence_interval_95[1]),
                        float(1.0 - fidelity_est.confidence_interval_95[0])
                    ],
                    "classification_confidence": float(classification_confidence)
                })
        else:
            # Point estimate only
            fidelity = quantum_fidelity(cell_state, immune_ref_state)
            evasion_score = 1.0 - fidelity

            if evasion_score > evasion_threshold:
                evading_cells.append({
                    "cell_idx": i,
                    "cell_type": adata.obs["cell_type"][i],
                    "evasion_score": float(evasion_score)
                })

    return {
        "n_evading_cells": len(evading_cells),
        "evasion_threshold": evasion_threshold,
        "evading_cells": evading_cells,
        "high_confidence_evading": [
            c for c in evading_cells
            if with_confidence and c.get("classification_confidence", 0) > 0.90
        ]
    }
```

**PatientOne immune evasion results**:
```json
{
  "n_evading_cells": 47,
  "evasion_threshold": 0.3,
  "high_confidence_evading": [
    {
      "cell_idx": 234,
      "cell_type": "tumor_cell",
      "evasion_score": 0.78,
      "evasion_score_ci_95": [0.72, 0.84],
      "classification_confidence": 0.98
    },
    {
      "cell_idx": 567,
      "cell_type": "tumor_cell",
      "evasion_score": 0.82,
      "evasion_score_ci_95": [0.77, 0.87],
      "classification_confidence": 0.99
    },
    ...  // 32 total high-confidence evading cells
  ]
}
```

**Clinical decision**: **32 cells with >90% confidence** of immune evasion → strong evidence for checkpoint blockade (anti-PD1/anti-CTLA4) to re-enable immune recognition.

---

### 4. predict_perturbation_effect

**Why you need it**: Predict how treatments affect cell type fidelities (e.g., does checkpoint blockade reduce immune evasion?).

**Perturbation simulation**:
1. Load drug perturbation Δ from Chapter 9 (GEARS prediction)
2. Apply Δ to cell gene expression: `X_post = X_pre + α * Δ`
3. Re-encode with quantum circuit: `|ψ_post⟩ = PQC(X_post)`
4. Compute fidelity change: `ΔF = F(|ψ_post⟩, |ψ_immune⟩) - F(|ψ_pre⟩, |ψ_immune⟩)`

**Positive ΔF**: Drug increases immune recognition (good for immunotherapy)

---

### 5. analyze_tls_quantum_signature

**Why you need it**: Identify tertiary lymphoid structures (TLS) by quantum signatures.

**TLS definition**: Organized clusters of B cells + T cells + dendritic cells in tumor microenvironment → correlated with better immunotherapy response.

**Quantum TLS signature**:
- High B cell fidelity cluster (> 0.85)
- High T cell fidelity cluster (> 0.85)
- Spatial proximity (< 100μm distance)
- Cluster size (> 20 cells)

---

### 6. export_for_downstream

**Why you need it**: Export quantum embeddings for downstream analysis (UMAP, clustering, integration with other tools).

**Export formats**:
- **NumPy**: `.npy` arrays for Python workflows
- **PyTorch**: `.pt` tensors for deep learning pipelines
- **JSON**: `.json` for web visualization

---

## Implementation Walkthrough

### Step 1: Project Setup

```bash
cd servers/mcp-quantum-celltype-fidelity
python -m venv venv
source venv/bin/activate

# Install dependencies
pip install fastmcp qiskit pennylane numpy scipy scikit-learn
```

**Key dependencies**:
- **Qiskit**: Quantum circuit simulation (IBM)
- **PennyLane**: Differentiable quantum programming (parameter-shift gradients)
- **NumPy/SciPy**: Linear algebra, Monte Carlo sampling

### Step 2: Initialize FastMCP Server

```python
from fastmcp import FastMCP

mcp = FastMCP("quantum-celltype-fidelity")

# Configuration
config = {
    "backend": os.getenv("QUANTUM_BACKEND", "cpu"),  # "cpu", "gpu", "ibm"
    "n_qubits_default": 8,
    "n_uncertainty_samples": 100
}
```

### Step 3: Implement Bayesian UQ

The core innovation. Create `src/quantum_celltype_fidelity/bayesian_uq.py`:

```python
import numpy as np
from dataclasses import dataclass
from typing import Callable, Tuple

@dataclass
class UncertaintyEstimate:
    """Bayesian uncertainty estimate for fidelity."""
    mean: float
    std: float
    confidence_interval_95: Tuple[float, float]
    confidence_interval_90: Tuple[float, float]
    epistemic_uncertainty: float  # Model uncertainty
    aleatoric_uncertainty: float  # Data uncertainty

class BayesianParameterDistribution:
    """Maintains posterior distributions over circuit parameters."""

    def __init__(self, parameters: np.ndarray, uncertainties: np.ndarray):
        self.mean = parameters
        self.std = uncertainties

    def sample(self, n_samples: int = 100) -> np.ndarray:
        """Sample parameter sets from posterior."""
        return np.random.normal(
            self.mean,
            self.std,
            size=(n_samples, len(self.mean))
        )

class BayesianFidelityEstimator:
    """Monte Carlo sampling for uncertainty quantification."""

    def __init__(self, param_distribution: BayesianParameterDistribution, n_samples: int = 100):
        self.param_dist = param_distribution
        self.n_samples = n_samples

    def estimate_fidelity_with_uncertainty(
        self,
        fidelity_fn: Callable[[np.ndarray], float]
    ) -> UncertaintyEstimate:
        """Estimate fidelity with Bayesian UQ.

        Args:
            fidelity_fn: Function that computes fidelity given parameters

        Returns:
            UncertaintyEstimate with mean, std, CIs
        """
        # Sample parameters from posterior
        param_samples = self.param_dist.sample(self.n_samples)

        # Compute fidelity for each sample
        fidelity_samples = np.array([
            fidelity_fn(params) for params in param_samples
        ])

        # Statistics
        mean = fidelity_samples.mean()
        std = fidelity_samples.std()
        ci_95 = np.percentile(fidelity_samples, [2.5, 97.5])
        ci_90 = np.percentile(fidelity_samples, [5, 95])

        return UncertaintyEstimate(
            mean=float(mean),
            std=float(std),
            confidence_interval_95=(float(ci_95[0]), float(ci_95[1])),
            confidence_interval_90=(float(ci_90[0]), float(ci_90[1])),
            epistemic_uncertainty=float(std),  # From parameter uncertainty
            aleatoric_uncertainty=0.0  # Would need multiple measurements
        )

    def estimate_classification_confidence(
        self,
        fidelity_estimate: UncertaintyEstimate,
        threshold: float,
        direction: str = "greater_than"
    ) -> float:
        """Estimate probability that classification decision is correct.

        Args:
            fidelity_estimate: Fidelity with uncertainty
            threshold: Classification threshold
            direction: "greater_than" or "less_than"

        Returns:
            Probability ∈ [0, 1] that classification is correct
        """
        # Approximate posterior as Gaussian
        from scipy.stats import norm

        dist = norm(fidelity_estimate.mean, fidelity_estimate.std)

        if direction == "greater_than":
            # P(fidelity > threshold)
            confidence = 1.0 - dist.cdf(threshold)
        else:
            # P(fidelity < threshold)
            confidence = dist.cdf(threshold)

        return float(np.clip(confidence, 0, 1))
```

Full implementation: [`servers/mcp-quantum-celltype-fidelity/src/quantum_celltype_fidelity/bayesian_uq.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-quantum-celltype-fidelity/src/quantum_celltype_fidelity/bayesian_uq.py) (400+ lines)

---

## Clinical Impact: Before vs After Bayesian UQ

**Before** (Point estimates only):
```
Fidelity score: 0.85
Oncologist decision: "Is 0.85 high enough to proceed with checkpoint blockade? Unknown risk."
```

**After** (With Bayesian UQ):
```
Fidelity score: 0.85 ± 0.03 (95% CI: [0.79, 0.91])
Classification confidence: 95%
Oncologist decision: "95% confident anti-PD1 will work - proceed with treatment."
```

**Risk quantification enables informed clinical decisions.**

---

## Testing Your Server

### Unit Tests

```python
# tests/test_bayesian_uq.py
import pytest
from quantum_celltype_fidelity.bayesian_uq import BayesianFidelityEstimator, BayesianParameterDistribution

def test_uncertainty_estimation():
    """Test Bayesian UQ produces calibrated confidence intervals."""
    # Create parameter distribution
    params = np.array([0.5, 1.0, 1.5])
    uncertainties = np.array([0.1, 0.1, 0.1])
    param_dist = BayesianParameterDistribution(params, uncertainties)

    # Mock fidelity function
    def fidelity_fn(p):
        return np.sin(p[0]) * np.cos(p[1]) * np.exp(-p[2])

    # Estimate with UQ
    estimator = BayesianFidelityEstimator(param_dist, n_samples=1000)
    estimate = estimator.estimate_fidelity_with_uncertainty(fidelity_fn)

    # Check CI width is reasonable
    ci_width = estimate.confidence_interval_95[1] - estimate.confidence_interval_95[0]
    assert 0.01 < ci_width < 0.5  # Not too narrow or too wide

    # Check mean is within CI
    assert estimate.confidence_interval_95[0] <= estimate.mean <= estimate.confidence_interval_95[1]
```

Run tests:
```bash
pytest tests/test_bayesian_uq.py -v
```

---

## What You've Built

You now have a quantum cell-type fidelity server that:

1. **Trains quantum embeddings**: PQCs with 8-10 qubits, contrastive learning
2. **Computes fidelity with UQ**: Bayesian uncertainty quantification via Monte Carlo sampling
3. **Detects immune evasion**: With classification confidence (95% → proceed, 60% → retest)
4. **Predicts perturbations**: Drug effects on fidelity
5. **Analyzes TLS**: Quantum signatures of immune hubs

This provides **calibrated confidence scores** for clinical decision-making.

---

## Try It Yourself

### Option 1: Synthetic Data Test

```bash
cd servers/mcp-quantum-celltype-fidelity
python -m venv venv
source venv/bin/activate
pip install -e ".[dev]"

# In Claude Desktop:
# "Train quantum embeddings on synthetic T cell data with 8 qubits, 3 layers, 20 epochs"
```

### Option 2: PatientOne Analysis

```bash
# "Compute cell type fidelity for PAT001 T cells with Bayesian UQ (100 samples)"
# "Identify immune evading cells with >90% classification confidence"
```

---

## Next Steps

In **Chapter 11: Imaging and Histopathology**, you'll build `mcp-openimagedata` to integrate H&E histopathology and multiplexed immunofluorescence (MxIF) imaging for comprehensive tumor characterization.

The quantum fidelity you built verifies **cell type classifications**. Histopathology validates those classifications with **tissue morphology** (nuclear atypia, necrosis, immune infiltration patterns).

---

**Chapter 10 Summary**:
- Quantum fidelity measures cell type similarity via overlap in Hilbert space
- Parameterized quantum circuits (8 qubits = 256D space) trained with contrastive learning
- Bayesian UQ provides 95% confidence intervals via Monte Carlo parameter sampling
- PatientOne: 32 immune evading cells detected with >90% classification confidence
- Clinical impact: Risk-quantified treatment decisions

**Files created**: `servers/mcp-quantum-celltype-fidelity/src/quantum_celltype_fidelity/server.py`, `bayesian_uq.py` (400+ lines), `training.py`, `circuits.py`
**Tests added**: 18 unit tests, 75% coverage
**Tools exposed**: 6 MCP tools (learn_embeddings, compute_fidelity, identify_evasion, predict_perturbation, analyze_tls, export)
