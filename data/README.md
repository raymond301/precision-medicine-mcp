# Data Directory

This directory contains sample data for testing and demonstrating the Precision Medicine MCP servers.

**⚠️ All data is 100% synthetic** - Created for demonstration and testing purposes only.

---

## 🏥 Synthetic Patient Datasets

| Patient ID | Cancer Type | Stage | Key Features | Documentation |
|------------|-------------|-------|--------------|---------------|
| **PAT001-OVC-2025** | Ovarian Cancer | Stage IV | BRCA1+, platinum-resistant | [📖 Details →](patient-data/PAT001-OVC-2025/README.md) |
| **PAT002-BC-2026** | Breast Cancer | Stage IIA | BRCA2+, ER+/PR+/HER2- | [📖 Details →](patient-data/PAT002-BC-2026/README.md) |

---

## PAT001: Stage IV Ovarian Cancer (PatientOne)

**Diagnosis:** High-Grade Serous Ovarian Carcinoma, platinum-resistant

**Data modalities:**
- Clinical data (demographics, CA-125 timeline)
- Genomic variants (VCF with BRCA1, TP53, PIK3CA mutations)
- Multi-omics data (RNA-seq, proteomics, phosphoproteomics)
- Spatial transcriptomics (10x Visium, 900 spots, 31 genes)
- Imaging data (H&E, immunofluorescence)

---

## PAT002: Stage IIA Breast Cancer

**Diagnosis:** ER+/PR+/HER2- Invasive Ductal Carcinoma, BRCA2 germline mutation

**Data modalities:**
- Clinical data (FHIR resources, CEA/CA 15-3 markers)
- Genomic variants (VCF with BRCA2, PIK3CA mutations)
- Multi-omics data (RNA-seq, proteomics, phosphoproteomics - pre/post treatment)
- Spatial transcriptomics (10x Visium, 900 spots, 35 genes)
- Imaging data (H&E, ER/PR/HER2/Ki67 immunofluorescence)
- Perturbation data (PD-1 knockout CRISPR screen)

---


