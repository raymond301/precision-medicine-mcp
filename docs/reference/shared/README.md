# Platform Overview (Canonical Reference)

The Precision Medicine MCP platform uses multiple MCP servers providing specialized tools for AI-orchestrated precision oncology analysis. For current counts, see the [Server Registry](server-registry.md).

---

## Server Status Matrix

| Server | Tools | Status | Description |
|--------|-------|--------|-------------|
| mcp-epic | 4 | ✅ Production (local only) | Epic FHIR with HIPAA de-identification |
| mcp-mockepic | 3 | 🎭 Mock by design | Synthetic EHR for demos |
| mcp-fgbio | 4 | ✅ Production | Reference genomes, FASTQ QC, gene annotations |
| mcp-multiomics | 10 | ✅ Production | RNA/Protein/Phospho integration, Stouffer meta-analysis |
| mcp-spatialtools | 14 | ✅ Production | Spatial transcriptomics, cell deconvolution |
| mcp-perturbation | 8 | ✅ Production | Perturbation prediction (GEARS, Nature Biotech 2024) |
| mcp-quantum-celltype-fidelity | 6 | ✅ Production | Quantum cell type fidelity (Qiskit) |
| mcp-openimagedata | 5 | ✅ Production | Histology image processing, MxIF compositing |
| mcp-deepcell | 4 | ✅ Production | DeepCell-TF cell segmentation (Cloud Run) |
| mcp-cell-classify | 3 | ✅ Production | Cell phenotype classification |
| mcp-tcga | 5 | ❌ Mocked | TCGA cohort comparison (GDC-ready) |
| mcp-huggingface | 3 | ❌ Mocked | ML model inference (HF-ready) |
| mcp-seqera | 3 | ❌ Mocked | Nextflow workflows (Seqera-ready) |
| mcp-patient-report | 5 | ✅ Production | PDF report generation |
| mcp-genomic-results | 4 | ✅ Production | Somatic variant/CNV parsing, HRD scoring |

**Summary:** 11 production-ready, 1 mock by design, 3 framework/utility (mocked)

## Architecture

```
┌─────────────────────────────────────────────────┐
│              User Interface Layer                │
│   Streamlit App │ Jupyter Notebook │ Dashboard   │
├─────────────────────────────────────────────────┤
│           AI Orchestration (Claude/Gemini)        │
│   Natural language → tool calls → integration    │
├─────────────────────────────────────────────────┤
│          MCP Servers (see server-registry.md)      │
│   Clinical │ Genomic │ Spatial │ Imaging │ ...   │
├─────────────────────────────────────────────────┤
│                 Data Layer                        │
│   GCS │ Epic FHIR │ Local Files │ TCGA/GDC      │
└─────────────────────────────────────────────────┘
```

## Data Modalities

1. **Clinical** — Epic FHIR with automatic de-identification
2. **Genomic** — WES/WGS somatic variants, CNV, germline risk, HRD scoring
3. **Multi-omics** — RNA-seq + proteomics + phosphoproteomics integration
4. **Spatial** — 10x Visium spatial transcriptomics, cell deconvolution
5. **Imaging** — H&E histopathology, multiplex immunofluorescence

---

