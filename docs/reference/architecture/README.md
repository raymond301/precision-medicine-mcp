# 🏗️ Precision Medicine MCP Server Architectures

Architecture documentation for the Precision Medicine MCP system, covering both modality-specific workflows and cross-cutting system design.

---

## System Overview

```mermaid
graph LR
    subgraph Input["📁 Data Sources"]
        EHR[Clinical<br/>Epic FHIR]
        SEQ[Genomics<br/>VCF/FASTQ]
        OMICS[Multiomics<br/>RNA/Protein]
        SPATIAL[Spatial<br/>Visium]
        IMG[Imaging<br/>H&E/MxIF]
        SCRNA[Single-cell<br/>scRNA-seq]
    end

    subgraph MCP["🔧 12 MCP Servers (69 Tools)"]
        direction TB
        S1[Clinical<br/>epic/mockepic]
        S2[Genomic<br/>fgbio/tcga]
        S3[Spatial<br/>spatialtools]
        S4[Multiomics<br/>multiomics]
        S5[Imaging<br/>openimagedata<br/>deepcell]
        S6[AI/ML<br/>huggingface<br/>seqera]
        S7[Perturbation<br/>perturbation<br/>GEARS]
        S8[Quantum<br/>quantum-celltype-fidelity]
    end

    subgraph Output["📊 Analysis Outputs"]
        TREAT[Treatment<br/>Recommendations]
        PREDICT[Response<br/>Predictions]
        VIZ[Visualizations<br/>& Reports]
        INSIGHTS[Therapeutic<br/>Targets]
    end

    EHR --> S1
    SEQ --> S2
    OMICS --> S4
    SPATIAL --> S3
    IMG --> S5
    SCRNA --> S7

    S1 --> MCP
    S2 --> MCP
    S3 --> MCP
    S4 --> MCP
    S5 --> MCP
    S6 --> MCP
    S7 --> MCP
    S8 --> MCP

    MCP --> TREAT
    MCP --> PREDICT
    MCP --> VIZ
    MCP --> INSIGHTS

    style Input fill:#e1f5ff,stroke:#0288d1,stroke-width:2px
    style MCP fill:#fff3e0,stroke:#f57c00,stroke-width:2px
    style Output fill:#f1f8e9,stroke:#689f38,stroke-width:2px
```

---

## 📊 Architecture by Analysis Modality

9 analysis modalities, 12 specialized servers, 69 tools:

| Modality | Servers | Tools | Status | Documentation |
|----------|---------|-------|--------|---------------|
| 🧬 **Clinical Data** | mcp-epic, mcp-mockepic | 7 (4+3) | ✅ Production/Mock | [clinical/README.md](clinical/README.md) |
| 🧪 **Genomic Cohorts** | mcp-tcga | 5 | ❌ Mocked (GDC-ready) | [genomic/README.md](genomic/README.md) |
| 🖼️ **Imaging** | mcp-openimagedata, mcp-deepcell | 9 (5+4) | ✅ Production (100%/100%) | [imaging/README.md](imaging/README.md) |
| 🔬 **Multiomics** | mcp-multiomics | 10 | ✅ Production (85%) | [multiomics/README.md](multiomics/README.md) |
| 📍 **Spatial Transcriptomics** | mcp-fgbio, mcp-spatialtools | 18 (4+14) | ✅ Production (95%) | [spatial-transcriptomics/README.md](spatial-transcriptomics/README.md) |
| 🎯 **Perturbation Prediction** | mcp-perturbation | 8 | ✅ Production (GEARS) | [perturbation/README.md](perturbation/README.md) |
| ⚛️ **Quantum Cell Type Fidelity** | mcp-quantum-celltype-fidelity | 6 | ✅ Production (Qiskit + Bayesian UQ) | [quantum/README.md](quantum/README.md) |
| 🤖 **AI/ML Inference** | mcp-huggingface | 3 | ❌ Mocked (HF-ready) | [ai-ml/README.md](ai-ml/README.md) |
| ⚙️ **Workflow Orchestration** | mcp-seqera | 3 | ❌ Mocked (Seqera-ready) | [workflow/README.md](workflow/README.md) |


---

## 🧬 1. Clinical Data Retrieval

**EHR integration for patient clinical context**

**Servers:** mcp-epic (real FHIR, local only) • mcp-mockepic (synthetic, GCP deployed)

**Key Features:**
- FHIR R4 data retrieval (demographics, diagnoses, labs, medications)
- HIPAA Safe Harbor de-identification
- Clinical-spatial outcome linkage

**Workflow:** `Patient EHR → FHIR API → De-identification → Clinical Data`

📖 **[Detailed Architecture →](clinical/README.md)**

---

## 🧪 2. Genomic Cohort Analysis

**TCGA cohort comparison for population-level genomic context**

**Server:** mcp-tcga (33 cancer types, 11,000+ samples)

**Key Features:**
- Gene expression comparison (z-scores, percentiles, p-values)
- Somatic mutation frequency queries
- Survival stratification (Kaplan-Meier, hazard ratios)

**Workflow:** `TCGA Database → Statistical Comparison → Survival Analysis → Integration`

📖 **[Detailed Architecture →](genomic/README.md)**

---

## 🖼️ 3. Imaging Analysis

**Histology and multiplexed immunofluorescence (MxIF) image processing**

**Servers:** mcp-openimagedata (100% real) • mcp-deepcell (100% real)

**Key Workflows:**
- **H&E (Brightfield):** Morphology assessment, necrosis identification
- **MxIF (Fluorescence):** Cell segmentation and quantification

📖 **[Detailed Architecture →](imaging/README.md)**

---

## 🔬 4. Multiomics Integration

**PDX multi-omics data integration with preprocessing and therapeutic target prediction**

**Server:** mcp-multiomics (21 tools, 85% real)

**Key Features:**
- Preprocessing pipeline (batch correction, KNN imputation, QC visualization)
- Association testing (HAllA with chunking)
- Meta-analysis (Stouffer's method)
- Therapeutic targets (kinase/TF/drug prediction)

**Workflow:** `RNA/Protein/Phospho → Validate → Preprocess → Integrate → HAllA → Meta-Analysis → Upstream Regulators`

📖 **[Detailed Architecture →](multiomics/README.md)**

---

## 📍 5. Spatial Transcriptomics

**Spatial gene expression analysis with tissue context**

**Servers:** mcp-fgbio (reference genomes, FASTQ QC) • mcp-spatialtools (spatial analysis, DE, pathway enrichment)

**Key Features:**
- **Analysis Tools (10):** Spatial autocorrelation (Moran's I), differential expression, batch correction, pathway enrichment, cell type deconvolution
- **Visualization Tools (4):** Spatial heatmaps, gene expression heatmaps, region composition charts
- **Bridge Tool:** Integrates with mcp-multiomics

**Workflows:** CSV (current) • FASTQ with STAR alignment (implemented)

📖 **[Detailed Architecture →](spatial-transcriptomics/README.md)**

---

## 🎯 6. Perturbation Prediction

**GEARS-based treatment response prediction using graph neural networks**

**Server:** mcp-perturbation (8 tools, production)

**Key Features:**
- **Model Training:** Setup and train GEARS GNN models on single-cell perturbation datasets
- **Response Prediction:** Predict cellular responses to genetic/pharmacological perturbations
- **Differential Expression:** Identify genes most affected by perturbations
- **Treatment Screening:** Test multiple therapies to find optimal responses

**Workflow:** `scRNA-seq Data → Load Dataset → Setup GEARS Model → Train → Predict Response → Differential Expression → Treatment Recommendations`

**Use Cases:**
- Predict T-cell response to checkpoint inhibitors (PD1/CTLA4)
- Screen PARP inhibitors vs platinum therapy for ovarian cancer
- Identify biomarkers of treatment sensitivity/resistance

**Technology:** GEARS (Graph-Enhanced Gene Activation Modeling) - Nature Biotechnology 2024

📖 **[Detailed Architecture →](perturbation/README.md)**

---

## 🤖 7. AI/ML Model Inference

**Genomic foundation model inference for cell type prediction and sequence embedding**

**Server:** mcp-huggingface (DNABERT-2, Geneformer, Nucleotide-Transformer)

**Key Features:**
- Cell type annotation from single-cell expression data
- DNA/RNA sequence embeddings (768-dimensional vectors)
- Variant effect scoring via embedding distance

**Workflow:** `Single-cell Data / DNA Sequences → Load Model → Generate Predictions → Integration`

📖 **[Detailed Architecture →](ai-ml/README.md)**

---

## ⚛️ 7. Quantum Cell Type Fidelity

**Quantum computing for cell type validation and immune evasion detection**

**Server:** mcp-quantum-celltype-fidelity (6 tools, production)

**Key Features:**
- **Quantum Embeddings:** Parameterized quantum circuits (PQCs) with 8-10 qubits
- **Fidelity Computation:** Quantum state overlap F = |⟨ψ_a|ψ_b⟩|² for cell similarity
- **Bayesian Uncertainty Quantification:** 95%/90% confidence intervals for clinical decisions (Phase 1, Jan 2026)
- **Immune Evasion Detection:** Identify tumor cells evading immune surveillance
- **TLS Analysis:** Characterize tertiary lymphoid structures with quantum signatures
- **Perturbation Prediction:** Simulate drug effects on quantum cell states

**Workflow:** `Spatial Data → Feature Encoding → Quantum Circuits → Contrastive Training → Fidelity Analysis → Immune Evasion Detection`

**Use Cases:**
- Detect tumor cells mimicking immune cells (low fidelity to canonical types)
- Characterize TLS immune organization via quantum coherence
- Validate GEARS perturbation predictions with quantum state changes
- Spatial mapping of cell type fidelity across tissue

**Technology:** Qiskit 1.0+ with parameter-shift rule for exact gradients

**Integration:** Works with mcp-perturbation for dual quantum+GEARS validation

📖 **[Detailed Architecture →](quantum/README.md)**

---

## ⚙️ 8. Workflow Orchestration

**Nextflow pipeline execution and monitoring via Seqera Platform**

**Server:** mcp-seqera (90+ nf-core workflows)

**Key Features:**
- nf-core pipeline submission (rnaseq, sarek, spatial, etc.)
- Multi-cloud execution (AWS, Azure, GCP, local HPC)
- Real-time progress monitoring and resource tracking

**Workflow:** `Select Pipeline → Configure → Choose Compute → Submit → Monitor → Retrieve Results`

📖 **[Detailed Architecture →](workflow/README.md)**

---

## 🏥 End-to-End Example: PatientOne

**Complete precision medicine workflow combining all 12 MCP servers**

**Use Case:** Stage IV High-Grade Serous Ovarian Cancer (HGSOC), platinum-resistant
**Patient:** PAT001-OVC-2025 (synthetic test case)
**Data Modalities:** Clinical (FHIR) • Genomic (VCF) • Multiomics (RNA/Protein/Phospho) • Spatial (Visium) • Imaging (H&E, MxIF) • Perturbation (scRNA-seq)

**Tests:**
- 🧬 TEST_1: Clinical data retrieval (mcp-epic)
- 🔬 TEST_2: Multiomics integration (mcp-multiomics)
- 📍 TEST_3: Spatial transcriptomics (mcp-spatialtools)
- 🖼️ TEST_4: Imaging analysis (mcp-openimagedata, mcp-deepcell)
- 🎯 TEST_5: Perturbation prediction (mcp-perturbation)
- ⚛️ TEST_6: Quantum cell type fidelity (mcp-quantum-celltype-fidelity)
- 🔄 TEST_7: Complete end-to-end workflow

📖 **[PatientOne Workflow →](../test-docs/patient-one-scenario/README.md)**
📖 **[PatientOne Architecture →](../test-docs/patient-one-scenario/architecture/overview.md)**

___

```mermaid

sequenceDiagram
    autonumber
    actor User as Clinician / Researcher
    participant AI as AI Orchestrator (Claude/Gemini)
    box Silver MCP Server Layer
        participant Clinical as Clinical (FHIR)
        participant Bio as Omics/Imaging/Spatial
        participant Quantum as Quantum/Perturbation
    end
    participant Output as Orchestrated Outputs

    User->>AI: Submit complex query (e.g., "Predict target for Patient X")
    
    activate AI
    AI->>AI: Analyze query & route to tools
    
    par Parallel Data Retrieval
        AI->>Clinical: Request EHR/FHIR data
        Clinical-->>AI: Patient history
        AI->>Bio: Fetch VCF/Imaging/Spatial data
        Bio-->>AI: Genomic & Spatial profiles
    end

    opt Advanced Modeling
        AI->>Quantum: Run Qiskit simulation / GEARS GNN
        Quantum-->>AI: Perturbation results
    end

    AI->>AI: Synthesize & Aggregate Data
    
    AI->>Output: Generate Targets & Visualizations
    Output-->>User: Present Final Insights
    deactivate AI

```

---

**Organization Principle:**
- `docs/architecture/` = High-level design & workflows by modality + cross-cutting architecture
- `servers/` = Detailed tool specifications & implementation
- `docs/` = Operational guides & deployment
- `tests/` = End-to-end use cases & validation
