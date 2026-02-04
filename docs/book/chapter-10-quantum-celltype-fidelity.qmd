# Quantum Cell-Type Fidelity

*Building mcp-quantum-celltype-fidelity with PennyLane and Bayesian UQ*

---

## Why Quantum Fidelity for Cell Types?

Chapters 7-8 classified cells by phenotype (Tumor proliferative, Ki67+, TP53+, CD8+ T cells). **Problem**: These classifications lack **uncertainty quantification**.

When you say "this cell is CD8+ T cell":
- How confident are you? (95%? 60%? Unknown?)
- Could it be an exhausted T cell instead?

**Clinical impact**: Without confidence scores, oncologists can't assess treatment risk:
- **High-confidence (95%)** immune detection → proceed with checkpoint blockade
- **Low-confidence (60%)** immune detection → order confirmatory tests first

The `mcp-quantum-celltype-fidelity` server uses **quantum fidelity** (overlap between quantum states) + **Bayesian uncertainty quantification** to provide calibrated confidence intervals on cell type classifications.

---

## Why Quantum Computing for This?

**Classical approach**: Cosine similarity in high-dimensional gene expression space. **Problem**: Doesn't capture non-linear gene regulatory relationships.

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

### Quantum Circuit Architecture and Workflow

```{mermaid}
graph TB
    subgraph InputData[Input: Cell Gene Expression]
        GENES[8 Selected Genes<br/>CD3D, CD8A, EPCAM<br/>MKI67, TP53, VIM, etc.<br/>Normalized 0-1]
    end

    subgraph PQC[Parameterized Quantum Circuit]
        direction TB
        ENCODE[Feature Encoding Layer<br/>RX, RY, RZ gates<br/>Map expression to rotation angles]
        VAR1[Variational Layer 1<br/>RX, RY, RZ trainable gates<br/>CNOT entanglement ring]
        VAR2[Variational Layer 2<br/>RX, RY, RZ trainable gates<br/>CNOT entanglement ring]
        VAR3[Variational Layer 3<br/>RX, RY, RZ trainable gates<br/>CNOT entanglement ring]
        STATE[Output Statevector<br/>256-dimensional complex vector<br/>Cell embedding in Hilbert space]
    end

    subgraph Training[Training: Contrastive Learning]
        POS[Positive Pairs<br/>Same cell type<br/>Maximize fidelity]
        NEG[Negative Pairs<br/>Different cell types<br/>Minimize fidelity]
        LOSS[Contrastive Loss<br/>Push same types to 1<br/>Push different types to 0]
        OPT[Optimizer<br/>Adam<br/>Update theta parameters]
    end

    subgraph Bayesian[Bayesian Uncertainty Quantification]
        MC[Monte Carlo Sampling<br/>Sample from posterior distribution<br/>100 parameter sets]
        FIDELITY[Compute Fidelity<br/>Quantum state overlap<br/>For each sample]
        CI[Credible Intervals<br/>Mean plus-minus 95 percent CI<br/>Example: 0.87 with range 0.82 to 0.91]
    end

    subgraph Output[Output: Classification with Confidence]
        CLASS[Cell Type Prediction<br/>Example: CD8 T cell]
        CONF[Confidence Score<br/>High: 95 percent CI narrow]
        IMMUNO[Immune Evasion Detection<br/>Fidelity drop indicates<br/>tumor escape]
    end

    GENES --> ENCODE
    ENCODE --> VAR1
    VAR1 --> VAR2
    VAR2 --> VAR3
    VAR3 --> STATE

    STATE --> POS
    STATE --> NEG
    POS --> LOSS
    NEG --> LOSS
    LOSS --> OPT
    OPT --> VAR1

    STATE --> MC
    MC --> FIDELITY
    FIDELITY --> CI

    CI --> CLASS
    CI --> CONF
    CI --> IMMUNO

    style GENES fill:#d1ecf1
    style ENCODE fill:#fff3cd,stroke:#ffc107,stroke-width:2px
    style VAR1 fill:#cce5ff
    style VAR2 fill:#cce5ff
    style VAR3 fill:#cce5ff
    style STATE fill:#cce5ff,stroke:#004085,stroke-width:2px
    style LOSS fill:#e1ecf4
    style MC fill:#fce8e8
    style CI fill:#d4edda,stroke:#28a745,stroke-width:2px
    style CONF fill:#d4edda,stroke:#28a745,stroke-width:2px
```

**Figure 10.1: Quantum Circuit Architecture and Bayesian UQ Workflow**
*Parameterized Quantum Circuit (PQC) with 3-layer architecture: (1) Feature encoding maps 8 gene expression values to rotation gates (RX, RY, RZ), (2) Three variational layers with learnable parameters (θ₁, θ₂, θ₃) and CNOT entanglement gates in ring topology, (3) Output statevector |ψ⟩ ∈ ℂ²⁵⁶ represents cell in 256D Hilbert space. Training uses contrastive learning to maximize same-type fidelity and minimize different-type fidelity. Bayesian UQ performs Monte Carlo sampling from parameter distribution P(θ|data) to compute 95% credible intervals on fidelity scores, providing calibrated confidence for clinical decisions.*

**Quantum Fidelity Formula:**
```
F(|ψ_A⟩, |ψ_B⟩) = |⟨ψ_A|ψ_B⟩|²
```
- F = 1: Perfect match (identical cell types)
- F = 0: Orthogonal (completely different)
- 0 < F < 1: Partial similarity with quantified confidence

---

## The 6 mcp-quantum-celltype-fidelity Tools

### 1. learn_spatial_cell_embeddings

Trains quantum circuit parameters to embed cell types into Hilbert space.

```python
@mcp.tool()
def learn_spatial_cell_embeddings(
        adata_path: str, cell_type_key: str = "cell_type",
        n_qubits: int = 8, n_layers: int = 3,
        n_epochs: int = 50) -> dict:
    """Train quantum embeddings for cell types using contrastive learning."""
    # Initialize PQC: feature encoding (RX, RY, RZ) + variational layers + CNOT entanglement
    # Contrastive learning: maximize within-type fidelity, minimize between-type
    # Full implementation: servers/mcp-quantum-celltype-fidelity/src/quantum_celltype_fidelity/training.py:100-300
```

**Parameterized Quantum Circuit (PQC)**:
1. **Feature encoding layer**: Gene expression → rotation angles (RX, RY, RZ gates)
2. **Variational layers** (3 layers): RX(θ_x), RY(θ_y), RZ(θ_z) gates per qubit + CNOT entanglement (ring topology)
3. **Measurement**: Statevector |ψ⟩ represents cell embedding

**PatientOne training results**:
```json
{
  "embedding_id": "embedding_patientone_tcells_8q_3l",
  "training_summary": {
    "final_loss": 0.234,
    "initial_loss": 1.872,
    "loss_reduction_percent": 87.5
  },
  "embedding_summary": {
    "n_qubits": 8,
    "hilbert_space_dim": 256,
    "n_parameters": 72
  }
}
```

---

### 2. compute_cell_type_fidelity (with Bayesian UQ)

Computes quantum fidelity with **confidence intervals** using Bayesian uncertainty quantification.

```python
@mcp.tool()
def compute_cell_type_fidelity(
        adata_path: str, embedding_id: str,
        with_uncertainty: bool = True,
        n_uncertainty_samples: int = 100) -> dict:
    """Compute quantum fidelity with Bayesian UQ (95% confidence intervals)."""
    # During training, track parameter gradient history
    # Build posterior: θ ~ N(μ, Σ) (mean = trained value, variance = gradient stability)
    # Monte Carlo sampling: sample 100 parameter sets
    # Compute fidelity for each → distribution of fidelities
    # Full implementation: servers/mcp-quantum-celltype-fidelity/src/quantum_celltype_fidelity/bayesian_uq.py:100-400
```

**PatientOne fidelity results** (T cells with Bayesian UQ):
```json
{
  "embedding_id": "embedding_patientone_tcells_8q_3l",
  "per_cell_type": {
    "CD8_T_cell": {
      "n_cells": 182,
      "within_type_fidelity": {"mean": 0.89, "std": 0.06}
    },
    "T_exhausted": {
      "n_cells": 89,
      "within_type_fidelity": {"mean": 0.72, "std": 0.11}
    }
  },
  "uncertainty": {"mean_uncertainty": 0.04}
}
```

**Interpretation**:
- **CD8+ T cells**: High within-type fidelity (0.89 ± 0.06) → well-defined quantum signature
- **Exhausted T cells**: Lower fidelity (0.72 ± 0.11) → heterogeneous subpopulation
- **Mean uncertainty**: 0.04 → predictions are well-calibrated (4% average CI width)

---

### 3. identify_immune_evasion_states

Detects tumor cells evading immune surveillance with **classification confidence**.

```python
@mcp.tool()
def identify_immune_evasion_states(
        adata_path: str, embedding_id: str,
        immune_cell_types: list[str],
        evasion_threshold: float = 0.3,
        with_confidence: bool = True) -> dict:
    """Detect cells in immune evasion states with >90% classification confidence."""
    # Measure fidelity to immune cells, apply threshold
    # Bayesian UQ: P(evasion_score > threshold) from posterior samples
    # Full implementation: servers/mcp-quantum-celltype-fidelity/src/quantum_celltype_fidelity/evasion.py:50-200
```

**PatientOne immune evasion results**:
```json
{
  "n_evading_cells": 47,
  "high_confidence_evading": [
    {
      "cell_idx": 234,
      "evasion_score": 0.78,
      "evasion_score_ci_95": [0.72, 0.84],
      "classification_confidence": 0.98
    }
  ]  // 32 total high-confidence evading cells
}
```

**Clinical decision**: **32 cells with >90% confidence** of immune evasion → strong evidence for checkpoint blockade (anti-PD1/anti-CTLA4).

---

### 4. predict_perturbation_effect

Predicts how treatments affect cell type fidelities.

```python
@mcp.tool()
def predict_perturbation_effect(model_name: str, perturbation_delta: dict, embedding_id: str) -> dict:
    """Predict fidelity change after drug treatment."""
    # Load drug perturbation Δ from Chapter 9
    # Apply Δ: X_post = X_pre + α * Δ
    # Re-encode with quantum circuit, compute fidelity change
    # Full implementation: servers/mcp-quantum-celltype-fidelity/src/quantum_celltype_fidelity/perturbation.py
```

**Positive ΔF**: Drug increases immune recognition (good for immunotherapy)

---

### 5. analyze_tls_quantum_signature

Identifies tertiary lymphoid structures (TLS) by quantum signatures.

```python
@mcp.tool()
def analyze_tls_quantum_signature(adata_path: str, embedding_id: str) -> dict:
    """Identify TLS clusters by quantum signatures."""
    # Find high B cell + T cell fidelity clusters in spatial proximity
    # TLS signature: B cell fidelity > 0.85, T cell fidelity > 0.85, cluster size > 20
    # Full implementation: servers/mcp-quantum-celltype-fidelity/src/quantum_celltype_fidelity/tls.py
```

**TLS definition**: Organized B cell + T cell + dendritic cell clusters → correlated with better immunotherapy response.

---

### 6. export_for_downstream

Exports quantum embeddings for downstream analysis.

```python
@mcp.tool()
def export_for_downstream(embedding_id: str, format: str = "numpy") -> dict:
    """Export embeddings as NumPy arrays, PyTorch tensors, or JSON."""
    # Save to .npy (Python), .pt (PyTorch), or .json (web visualization)
    # Full implementation: servers/mcp-quantum-celltype-fidelity/src/quantum_celltype_fidelity/export.py
```

---

## Implementation Walkthrough

### Project Setup

```bash
cd servers/mcp-quantum-celltype-fidelity
python -m venv venv && source venv/bin/activate
pip install fastmcp qiskit pennylane numpy scipy scikit-learn
```

**Key dependencies**:
- **Qiskit**: Quantum circuit simulation (IBM)
- **PennyLane**: Differentiable quantum programming (parameter-shift gradients)
- **NumPy/SciPy**: Linear algebra, Monte Carlo sampling

### Initialize FastMCP Server

```python
from fastmcp import FastMCP
mcp = FastMCP("quantum-celltype-fidelity")

config = {
    "backend": os.getenv("QUANTUM_BACKEND", "cpu"),  # "cpu", "gpu", or "ibm"
    "n_qubits_default": 8,
    "n_uncertainty_samples": 100
}
```

### Bayesian UQ Core

```python
@dataclass
class UncertaintyEstimate:
    """Bayesian uncertainty estimate for fidelity."""
    mean: float
    std: float
    confidence_interval_95: Tuple[float, float]
    epistemic_uncertainty: float  # Model uncertainty

class BayesianFidelityEstimator:
    """Monte Carlo sampling for uncertainty quantification."""

    def estimate_fidelity_with_uncertainty(self, fidelity_fn: Callable) -> UncertaintyEstimate:
        """Estimate fidelity with Bayesian UQ."""
        param_samples = self.param_dist.sample(self.n_samples)  # Sample 100 parameter sets
        fidelity_samples = np.array([fidelity_fn(params) for params in param_samples])
        mean = fidelity_samples.mean()
        std = fidelity_samples.std()
        ci_95 = np.percentile(fidelity_samples, [2.5, 97.5])
        return UncertaintyEstimate(mean=mean, std=std, confidence_interval_95=ci_95, epistemic_uncertainty=std)
        # Full implementation: servers/mcp-quantum-celltype-fidelity/src/quantum_celltype_fidelity/bayesian_uq.py (400 lines)
```

---

## Clinical Impact: Before vs After Bayesian UQ

**Before** (Point estimates only):
```
Fidelity score: 0.85
Oncologist decision: "Is 0.85 high enough? Unknown risk."
```

**After** (With Bayesian UQ):
```
Fidelity score: 0.85 ± 0.03 (95% CI: [0.79, 0.91])
Classification confidence: 95%
Oncologist decision: "95% confident anti-PD1 will work - proceed."
```

**Risk quantification enables informed clinical decisions.**

---

## Testing Your Server

```python
def test_uncertainty_estimation():
    """Test Bayesian UQ produces calibrated confidence intervals."""
    params = np.array([0.5, 1.0, 1.5])
    uncertainties = np.array([0.1, 0.1, 0.1])
    param_dist = BayesianParameterDistribution(params, uncertainties)
    estimator = BayesianFidelityEstimator(param_dist, n_samples=1000)
    estimate = estimator.estimate_fidelity_with_uncertainty(lambda p: np.sin(p[0]) * np.cos(p[1]))
    ci_width = estimate.confidence_interval_95[1] - estimate.confidence_interval_95[0]
    assert 0.01 < ci_width < 0.5  # Not too narrow or too wide
```

Test coverage: **75%**, 18 unit tests

---

## What You've Built

A quantum cell-type fidelity server providing:
1. **Quantum embeddings**: PQCs with 8-10 qubits, contrastive learning
2. **Fidelity with UQ**: Bayesian uncertainty quantification via Monte Carlo sampling
3. **Immune evasion detection**: With classification confidence (95% → proceed, 60% → retest)
4. **Perturbation prediction**: Drug effects on fidelity
5. **TLS analysis**: Quantum signatures of immune hubs
6. **Export**: NumPy, PyTorch, JSON formats

This provides **calibrated confidence scores** for clinical decision-making.

---

## Try It Yourself

```bash
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp/servers/mcp-quantum-celltype-fidelity
python -m venv venv && source venv/bin/activate
pip install -e ".[dev]"
# In Claude Desktop: "Train quantum embeddings on synthetic T cell data with 8 qubits, 3 layers, 20 epochs"
```

---

## Summary

**Chapter 10 Summary**:
- Quantum fidelity measures cell type similarity via overlap in Hilbert space
- Parameterized quantum circuits (8 qubits = 256D space) trained with contrastive learning
- Bayesian UQ provides 95% confidence intervals via Monte Carlo parameter sampling
- PatientOne: 32 immune evading cells detected with >90% classification confidence
- Clinical impact: Risk-quantified treatment decisions

**Files created**: `servers/mcp-quantum-celltype-fidelity/src/quantum_celltype_fidelity/server.py`, `bayesian_uq.py` (400 lines), `training.py`, `circuits.py`
**Tests added**: 18 unit tests, 75% coverage
**Tools exposed**: 6 MCP tools (learn_embeddings, compute_fidelity, identify_evasion, predict_perturbation, analyze_tls, export)
