# Precision Medicine MCP System - Quick Summary

**For sharing via email or quick reference**

---

```mermaid
graph LR
    U[Streamlit UI<br>Claude / Gemini] --> O[AI Orchestrator]
    O --> S1[Clinical<br>mcp-epic]
    O --> S2[Genomic<br>fgbio · genomic-results]
    O --> S3[Multi-Omics<br>multiomics]
    O --> S4[Spatial<br>spatialtools]
    O --> S5[Imaging<br>deepcell · openimagedata]
    O --> S6[Prediction<br>perturbation · quantum]
    O --> S7[Reporting<br>patient-report · cell-classify]
    S1 & S2 & S3 & S4 & S5 & S6 & S7 --> GCP[GCP Cloud Run]
```

## System Overview

A production-ready AI-orchestrated platform for precision oncology that integrates clinical, genomic, multi-omics, spatial transcriptomics, and imaging data through **15 MCP servers (80 tools)** on GCP Cloud Run — reducing production analysis from 40 hours to 2-5 hours at $324-702 per patient vs. $6,000-9,000 traditional ($3,137 avg savings). DRY_RUN demos complete in 25-35 minutes at ~$1. The system is HIPAA-compliant, supports Claude and Gemini via a multi-provider Streamlit UI, and includes a live monitoring dashboard, comprehensive audience-specific documentation, and a validated end-to-end demo (**PatientOne**: Stage IV ovarian cancer).

---

## MCP Servers (15 Total, 80 Tools)

11 production-ready, 1 local-only (Epic FHIR), 1 mock by design, 3 framework/utility. Covers genomics, multi-omics, spatial transcriptomics, imaging, cell segmentation, perturbation prediction, quantum computing, and patient reporting.

> **Full details:** See [Server Registry](reference/shared/server-registry.md)

---

## Key Metrics

- **Time Reduction:** 40 hours → 2-5 hours production (8-20x faster with AI orchestration)
- **Cost per Patient:** $324-702 production (incl. infrastructure) vs. $6,000-9,000 traditional
- **At Scale:** $16/patient at 500 patients/month (~$8,000/month total)
- **Data Integration:** 5 modalities (clinical, genomic, multi-omics, spatial, imaging)
- **HIPAA Compliant:** De-identification, 10-year audit logs, VPC isolation
- **Deployment:** 15 servers on GCP Cloud Run (us-central1), scale-to-zero
- **UI:** Streamlit with Claude & Gemini 3 support + live monitoring dashboard

---

## Live Demo

- **Streamlit UI (Developer):** streamlit app — Multi-provider (Claude + Gemini), 14 example prompts, 13 of 15 servers selectable (excludes local-only mcp-epic and mcp-genomic-results)
- **Streamlit UI (Student):** student streamlit app -  for 6-week bioinformatics course (Gemini-only, 6 tested prompts, 3 default servers). URL shared with enrolled students only.
- **PatientOne Scenario:** Stage IV ovarian cancer complete analysis
- **Test Data:** Available in GCS at `gs://sample-inputs-patientone/`
- **6-Week Course:** [Study group curriculum](../ui/streamlit-app-students/for-instructors/6-wk-course.md) for bioinformatics educators

---

## Documentation

- **For Hospitals:** [docs/for-hospitals/](for-hospitals/)
- **For Researchers:** [docs/for-researchers/](for-researchers/)
- **For Developers:** [docs/for-developers/](for-developers/)
- **Executive Summary:** [docs/EXECUTIVE_SUMMARY.md](for-funders/EXECUTIVE_SUMMARY.md)
- **Architecture:** [docs/architecture/](reference/architecture)

---

**Repository:** https://github.com/lynnlangit/precision-medicine-mcp

**Status:** Production-ready, validated deployment (February 2026). Includes dedicated student app for 6-week bioinformatics course.

**License:** Apache 2.0
