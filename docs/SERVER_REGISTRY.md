# MCP Server Registry - Quick Reference

**Total Servers:** 12 | **Production Ready:** 7 (58%) | **Total Tools:** 124

📁 **[Individual Server Documentation →](../servers/README.md)**

---

## Production Servers ✅

| Server | Tools | Status | Key Capabilities | Documentation |
|--------|-------|--------|------------------|---------------|
| **mcp-fgbio** | 9 | 95% Real | FASTQ/VCF QC, genome refs, variant calling | [README](../servers/mcp-fgbio/README.md) |
| **mcp-multiomics** | 21 | 85% Real | HAllA integration, Stouffer meta-analysis, upstream regulators | [README](../servers/mcp-multiomics/README.md) |
| **mcp-spatialtools** | 23 | 95% Real | Spatial DE, STAR alignment, ComBat, pathway enrichment | [README](../servers/mcp-spatialtools/README.md) |
| **mcp-perturbation** | 8 | 100% Real | GEARS GNN treatment response, perturbation prediction | [README](../servers/mcp-perturbation/README.md) |
| **mcp-quantum-celltype-fidelity** | 6 | 100% Real | Quantum PQCs, fidelity analysis, Bayesian UQ, immune evasion | [README](../servers/mcp-quantum-celltype-fidelity/README.md) |
| **mcp-deepcell** | 4 | 100% Real | DeepCell-TF segmentation, nuclear/membrane models, MxIF phenotyping | [README](../servers/mcp-deepcell/README.md) |
| **mcp-epic** | 9 | 100% Real | FHIR R4 API, real EHR integration (local deployment only) | [README](../servers/mcp-epic/README.md) |

---

## Partial Implementation ⚠️

| Server | Tools | Status | Notes | Documentation |
|--------|-------|--------|-------|---------------|
| **mcp-openimagedata** | 7 | 60% Real | Basic histology working, registration/advanced features mocked | [README](../servers/mcp-openimagedata/README.md) |

---

## Mock Servers (For Workflow Testing) ❌

| Server | Tools | Purpose | Documentation |
|--------|-------|---------|---------------|
| **mcp-tcga** | 11 | TCGA cohort queries, survival analysis (synthetic) | [README](../servers/mcp-tcga/README.md) |
| **mcp-huggingface** | 7 | ML model inference (API ready, awaiting models) | [README](../servers/mcp-huggingface/README.md) |
| **mcp-seqera** | 7 | Nextflow workflow orchestration (demo) | [README](../servers/mcp-seqera/README.md) |
| **mcp-mockepic** | 8 | Synthetic FHIR data for testing (by design) | [README](../servers/mcp-mockepic/README.md) |

---

## Status Legend

- ✅ **Production Ready**: Real APIs, extensively tested, validated outputs
- ⚠️ **Partial**: Core features real, some components mocked
- ❌ **Mocked**: Demonstration/workflow testing only - **DO NOT USE FOR RESEARCH**

---

## Quick Find

### By Analysis Type
- **Clinical Data**: mcp-epic (real EHR), mcp-mockepic (synthetic)
- **Genomics**: mcp-fgbio (QC/variants), mcp-tcga (cohort comparison - mocked)
- **Multi-omics**: mcp-multiomics (integration/meta-analysis)
- **Spatial**: mcp-spatialtools (spatial transcriptomics)
- **Imaging**: mcp-deepcell (cell segmentation), mcp-openimagedata (histology - partial)
- **Treatment**: mcp-perturbation (GEARS prediction), mcp-quantum-celltype-fidelity (quantum fidelity)
- **AI/ML**: mcp-huggingface (model inference - mocked)
- **Workflows**: mcp-seqera (Nextflow - mocked)

### By Production Readiness
- **Ready for Research**: mcp-fgbio, mcp-multiomics, mcp-spatialtools, mcp-perturbation, mcp-quantum-celltype-fidelity, mcp-deepcell, mcp-epic
- **Not Ready**: mcp-tcga, mcp-huggingface, mcp-seqera (synthetic data)
- **Mock by Design**: mcp-mockepic (testing only)
- **Partial**: mcp-openimagedata (verify per tool)

---

## GCP Deployment Status

**Deployed to Cloud Run:** 11/12 servers (mcp-epic is local-only)  
**Region:** us-central1  
**Transport:** SSE (Server-Sent Events) over HTTPS  
**Last Deployment:** 2026-01-31  

See [Hospital Deployment Guide](for-hospitals/README.md) for HIPAA-compliant production setup.

---

**Last Updated:** 2026-01-31
**Maintained By:** Precision Medicine MCP Team
