# MCP Server Registry - Quick Reference

**Total Servers:** 15 | **Production Ready:** 11 (73%) | **Total Tools:** 80

📁 **[Individual Server Documentation →](../../servers/README.md)**

---

## Production Servers ✅

| Server | Tools | Status | Key Capabilities | Documentation |
|--------|-------|--------|------------------|---------------|
| **mcp-fgbio** | 4 | 95% Real | FASTQ/VCF QC, genome refs, variant calling | [README](../../servers/mcp-fgbio/README.md) |
| **mcp-multiomics** | 10 | 95% Real | HAllA integration, Stouffer meta-analysis, upstream regulators, heatmap, PCA | [README](../../servers/mcp-multiomics/README.md) |
| **mcp-spatialtools** | 14 | 95% Real | Spatial DE, STAR alignment, ComBat, pathway enrichment | [README](../../servers/mcp-spatialtools/README.md) |
| **mcp-perturbation** | 8 | 100% Real | GEARS GNN treatment response, perturbation prediction | [README](../../servers/mcp-perturbation/README.md) |
| **mcp-quantum-celltype-fidelity** | 6 | 100% Real | Quantum PQCs, fidelity analysis, Bayesian UQ, immune evasion | [README](../../servers/mcp-quantum-celltype-fidelity/README.md) |
| **mcp-deepcell** | 3 | 100% Real | DeepCell-TF segmentation, nuclear/membrane models, per-cell marker quantification | [README](../../servers/mcp-deepcell/README.md) |
| **mcp-cell-classify** | 3 | 100% Real | Cell phenotype classification, multi-marker phenotyping, phenotype visualization | [README](../../servers/mcp-cell-classify/README.md) |
| **mcp-epic** | 4 | 100% Real | FHIR R4 API, real EHR integration (local deployment only) | [Source](../../servers/mcp-epic/) |
| **mcp-openimagedata** | 5 | 100% Real | PIL image loading, scikit-image registration + feature extraction, MxIF compositing, H&E annotation | [README](../../servers/mcp-openimagedata/README.md) |
| **mcp-patient-report** | 5 | 100% Real | Patient-facing PDF reports, plain-language summaries, clinician review gate | [README](../../servers/mcp-patient-report/README.md) |
| **mcp-genomic-results** | 4 | 100% Real | Somatic variant/CNV parsing, clinical annotations, HRD scoring | [README](../../servers/mcp-genomic-results/README.md) |

---

## Mock Servers (For Workflow Testing) ❌

| Server | Tools | Purpose | Documentation |
|--------|-------|---------|---------------|
| **mcp-tcga** | 5 | TCGA cohort queries, survival analysis (synthetic) | [README](../../servers/mcp-tcga/README.md) |
| **mcp-huggingface** | 3 | ML model inference (API ready, awaiting models) | [Source](../../servers/mcp-huggingface/) |
| **mcp-seqera** | 3 | Nextflow workflow orchestration (demo) | [Source](../../servers/mcp-seqera/) |
| **mcp-mockepic** | 3 | Synthetic FHIR data for testing (by design) | [Source](../../servers/mcp-mockepic/) |

---

## Status Legend

- ✅ **Production Ready**: Real APIs, extensively tested, validated outputs
- ⚠️ **Partial**: Core features real, some components mocked
- ❌ **Mocked**: Demonstration/workflow testing only - **DO NOT USE FOR RESEARCH**

---

## Quick Find

### By Analysis Type
- **Clinical Data**: mcp-epic (real EHR), mcp-mockepic (synthetic)
- **Genomics**: mcp-fgbio (QC/variants), mcp-genomic-results (somatic/CNV/HRD), mcp-tcga (cohort comparison - mocked)
- **Multi-omics**: mcp-multiomics (integration/meta-analysis)
- **Spatial**: mcp-spatialtools (spatial transcriptomics)
- **Imaging**: mcp-deepcell (cell segmentation + quantification), mcp-cell-classify (phenotype classification), mcp-openimagedata (histology + registration + features)
- **Treatment**: mcp-perturbation (GEARS prediction), mcp-quantum-celltype-fidelity (quantum fidelity)
- **Reports**: mcp-patient-report (patient-facing summaries)
- **AI/ML**: mcp-huggingface (model inference - mocked)
- **Workflows**: mcp-seqera (Nextflow - mocked)

### By Production Readiness
- **Ready for Research**: mcp-fgbio, mcp-multiomics, mcp-spatialtools, mcp-perturbation, mcp-quantum-celltype-fidelity, mcp-deepcell, mcp-cell-classify, mcp-epic, mcp-openimagedata, mcp-patient-report, mcp-genomic-results
- **Not Ready**: mcp-tcga, mcp-huggingface, mcp-seqera (synthetic data)
- **Mock by Design**: mcp-mockepic (testing only)

---

**Last Updated:** 2026-02-09
**Maintained By:** Precision Medicine MCP Team
