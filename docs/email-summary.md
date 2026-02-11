# Precision Medicine MCP System - Quick Summary

**For sharing via email or quick reference**

---

## System Overview

**A production-ready AI-orchestrated platform for precision oncology research**, integrating clinical (FHIR), genomic, multi-omics, spatial transcriptomics, and imaging data through 15 specialized MCP servers (136 bioinformatics tools) deployed on GCP Cloud Run. The system uses Claude or Gemini 3 to orchestrate complex multi-step analyses via natural language, reducing analysis time from 40 hours of manual bioinformatics work to 35 minutes of AI-orchestrated execution. Includes HIPAA compliance (de-identification, 10-year audit logs, bias detection), Bayesian uncertainty quantification for confident clinical decisions, and achieves 95% cost savings (~$1-2 per analysis vs. $3,200 traditional). **PatientOne** demonstrates the complete workflow: Stage IV ovarian cancer analysis integrating clinical data, somatic variants, multi-omics (RNA/protein/phospho), spatial transcriptomics, and H&E imaging—with treatment response prediction via GEARS perturbation modeling and quantum fidelity analysis.

**Technical implementation**: Multi-provider Streamlit UI supporting Claude (Anthropic's native MCP integration) and Gemini 3 (custom SSE-based MCP client with manual tool orchestration), with 11 production servers (fgbio, multiomics, spatialtools, perturbation, quantum-celltype-fidelity, deepcell, cell-classify, openimagedata, patient-report, genomic-results, epic), 1 local-only (epic), and 4 mock servers. Includes live monitoring dashboard for health tracking and token usage. The system handles GCS folder URIs (loads all files automatically), supports up to 30 tool-calling iterations for complex workflows, and includes comprehensive documentation for hospitals, researchers, developers, educators, patients, and funders. Validated against actual GCP deployment with 2026 pricing: ~$0.02-0.21 per Cloud Run analysis, free tier covers ~83 hours/month of testing.

---

## MCP Servers (15 Total)

### Production Servers (11 deployed)
1. **mcp-fgbio** - Genomic reference validation
2. **mcp-multiomics** - Multi-omics data integration
3. **mcp-spatialtools** - Spatial transcriptomics analysis
4. **mcp-perturbation** - Treatment response prediction (GEARS)
5. **mcp-quantum-celltype-fidelity** - Quantum computing cell analysis (Qiskit)
6. **mcp-deepcell** - Cell segmentation (DeepCell-TF)
7. **mcp-cell-classify** - Cell phenotype classification
8. **mcp-openimagedata** - Histology image processing (registration, feature extraction, MxIF compositing)
9. **mcp-patient-report** - Patient-facing PDF report generation
10. **mcp-genomic-results** - Somatic variant/CNV parsing with clinical annotations and HRD scoring

### Local Only (1)
11. **mcp-epic** - Clinical FHIR data (Epic integration)

### Mock by Design (1)
12. **mcp-mockepic** - Mock FHIR data (intentionally synthetic for demos)

### Mock Servers (3)
13. **mcp-tcga** - Cancer genomics data
14. **mcp-seqera** - Workflow automation platform
15. **mcp-huggingface** - AI model inference

**Total:** 136 bioinformatics tools across all servers (including Bayesian uncertainty quantification for quantum predictions)

---

## Key Metrics

- **Time Reduction:** 40 hours → 35 minutes (95% faster)
- **Cost Reduction:** $3,200 → $1-2 per analysis (95% savings)
- **Data Integration:** 5 modalities (clinical, genomic, multi-omics, spatial, imaging)
- **HIPAA Compliant:** De-identification, 10-year audit logs, VPC isolation
- **Deployment:** GCP Cloud Run (us-central1)
- **UI:** Streamlit with Claude & Gemini 3 support
- **Monitoring:** Live dashboard with token usage tracking

---

## Live Demo

- **Streamlit UI:** https://streamlit-mcp-chat-ondu7mwjpa-uc.a.run.app
- **PatientOne Scenario:** Stage IV ovarian cancer complete analysis
- **Test Data:** Available in GCS at `gs://sample-inputs-patientone/`

---

## Documentation

- **For Hospitals:** [docs/for-hospitals/](for-hospitals/)
- **For Researchers:** [docs/for-researchers/](for-researchers/)
- **For Developers:** [docs/for-developers/](for-developers/)
- **Executive Summary:** [docs/EXECUTIVE_SUMMARY.md](for-funders/EXECUTIVE_SUMMARY.md)
- **Architecture:** [docs/architecture/](reference/architecture)

---

**Repository:** https://github.com/lynnlangit/precision-medicine-mcp

**Status:** Production-ready, validated deployment (February 2026)

**License:** Apache 2.0
