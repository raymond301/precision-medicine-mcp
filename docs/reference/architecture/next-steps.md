# Next Steps & Enhancements by MCP Server

**Last Updated:** 2026-02-15

---

## Summary

| # | Server | Status | Top Enhancement | Priority |
|---|--------|--------|-----------------|----------|
| 1 | mcp-epic | Production | FHIR paging + retry/backoff | Medium |
| 2 | mcp-mockepic | Mock (by design) | FHIR R4-compliant output format | Low |
| 3 | mcp-cell-classify | Production | FlowSOM/Leiden clustering methods | Medium |
| 4 | mcp-deepcell | Production (Phase 1) | Morphology, tracking, export tools (Phase 2) | High |
| 5 | mcp-fgbio | Production (95% Real) | Finish `fetch_reference_genome`; add BAM tools | Medium |
| 6 | mcp-genomic-results | Production (100% Real) | TMB calculation; MSI detection | High |
| 8 | mcp-multiomics | ~95% real | Replace hardcoded upstream regulator DBs with live API | Medium |
| 9 | mcp-openimagedata | Production | OpenSlide WSI support (.svs, .ndpi) | Medium |
| 10 | mcp-patient-report | Production (100% Real) | Implement `approve_patient_report` workflow | High |
| 11 | mcp-perturbation | Production (GEARS) | GEO dataset download; GPU acceleration | Medium |
| 12 | mcp-quantum-celltype-fidelity | Production (CPU) | GPU backend (cuQuantum); IBM Quantum hardware | Medium |
| 13 | mcp-seqera | Mocked | Wire up real Seqera Platform API | High |
| 14 | mcp-spatialtools | 95% real | Harmony/Scanorama batch correction | Low |
| 15 | mcp-tcga | Mocked | Wire up real GDC API | High |

---

## Detailed Enhancements

### 1. mcp-epic (Production)

| Enhancement | Type | Effort |
|-------------|------|--------|
| Add retry/backoff for Epic FHIR API calls | Reliability | Small |
| Add FHIR Bundle paging for large result sets | Feature | Small |
| Add `get_patient_procedures` tool | Feature | Small |
| Add `get_diagnostic_reports` tool | Feature | Small |

### 2. mcp-mockepic (Mock by design)

| Enhancement | Type | Effort |
|-------------|------|--------|
| Return valid FHIR R4 JSON (not flat dicts) | Quality | Medium |
| Add configurable synthetic patient profiles | Feature | Small |
| Implement real code paths for `link_spatial_to_clinical` | Feature | Medium |

### 3. mcp-cell-classify (Production)

| Enhancement | Type | Effort |
|-------------|------|--------|
| Add FlowSOM / Leiden unsupervised clustering | Feature | Medium |
| Add continuous intensity heatmap visualization | Feature | Small |
| Add batch processing tool (directory of masks) | Feature | Small |
| Add probabilistic confidence calibration | Quality | Medium |

### 4. mcp-deepcell (Production, Phase 1)

| Enhancement | Type | Effort |
|-------------|------|--------|
| `measure_cell_morphology` — area, eccentricity, solidity | Feature | Medium |
| `track_cells_over_time` — cell tracking across time series | Feature | Large |
| `export_cell_data` — CSV, AnnData, CellProfiler formats | Feature | Medium |
| Add comprehensive test suite (unit + integration) | Quality | Medium |

### 5. mcp-fgbio (Production, 95% Real)

| Enhancement | Type | Effort |
|-------------|------|--------|
| Finish `fetch_reference_genome` real download path | Fix | Small |
| Add `index_reference_genome` (BWA/STAR/bowtie2) | Feature | Medium |
| Add BAM/CRAM consensus calling tools (fgbio Java) | Feature | Large |
| Add retry/backoff for NCBI FTP downloads | Reliability | Small |

### 6. mcp-genomic-results (Production, 100% Real)

| Enhancement | Type | Effort |
|-------------|------|--------|
| Add TMB (Tumor Mutational Burden) calculation | Feature | Medium |
| Add MSI (Microsatellite Instability) detection | Feature | Medium |
| Upgrade HRD scoring from simplified POC to clinical-grade | Quality | Large |
| Replace hardcoded annotations with ClinVar/VEP API | Feature | Large |
| Add FHIR DiagnosticReport export | Feature | Medium |

### 7. mcp-multiomics (~95% real)

| Enhancement | Type | Effort |
|-------------|------|--------|
| ~~Implement `create_multiomics_heatmap`~~ | ~~Fix~~ | Done |
| ~~Implement `run_multiomics_pca`~~ | ~~Fix~~ | Done |
| Replace hardcoded upstream regulator DBs with live API | Feature | Large |
| Test and document R-based HAllA integration path | Quality | Small |

### 9. mcp-openimagedata (Production)

| Enhancement | Type | Effort |
|-------------|------|--------|
| Add OpenSlide integration for WSI formats (.svs, .ndpi, .czi) | Feature | Large |
| Implement true non-rigid deformable registration | Feature | Large |
| Add automated necrosis/cellularity detection | Feature | Medium |
| Add spatial correlation with gene expression data | Feature | Medium |

### 10. mcp-patient-report (Production, 100% Real)

| Enhancement | Type | Effort |
|-------------|------|--------|
| Implement `approve_patient_report` full workflow | Feature | Medium |
| Add FHIR DocumentReference integration | Feature | Medium |
| Add automated readability scoring (Flesch-Kincaid) | Quality | Small |
| Improve multi-modal consistency checks (placeholder) | Quality | Medium |
| Create distinct one-page report template | Feature | Small |

### 11. mcp-perturbation (Production, GEARS)

| Enhancement | Type | Effort |
|-------------|------|--------|
| Implement GEO dataset download (`_download_from_geo`) | Feature | Medium |
| Add authentication to Cloud Run endpoint | Security | Small |
| GPU upgrade for faster GEARS training | Performance | Medium |
| Test with real PatientOne data | Quality | Small |

### 12. mcp-quantum-celltype-fidelity (Production, CPU)

| Enhancement | Type | Effort |
|-------------|------|--------|
| GPU backend via cuQuantum (5-10x speedup) | Performance | Large |
| IBM Quantum hardware execution (currently NotImplementedError) | Feature | Large |
| Implement aleatoric uncertainty quantification | Feature | Medium |
| Add 3D spatial coordinate support | Feature | Medium |
| Amplitude amplification for rare cell types | Feature | Medium |
| Wire up GEARS integration for combined analysis | Feature | Medium |

### 13. mcp-seqera (Mocked)

| Enhancement | Type | Effort |
|-------------|------|--------|
| Wire up real Seqera Platform API (launch, monitor) | Feature | Large |
| Implement retry/circuit-breaker for API calls | Reliability | Small |
| Add workflow monitoring / log retrieval tool | Feature | Medium |
| Add compute environment management tool | Feature | Medium |

### 14. mcp-spatialtools (95% real)

| Enhancement | Type | Effort |
|-------------|------|--------|
| Harmony batch correction | Feature | Medium |
| Scanorama batch correction | Feature | Medium |
| ComBat-seq for count data | Feature | Small |
| Reactome / WikiPathways databases | Feature | Medium |
| Local Moran's I (LISA) | Feature | Medium |
| Geary's C, Getis-Ord Gi* statistics | Feature | Medium |
| Gene ID conversion (ENSEMBL/Entrez to Symbol) | Feature | Small |
| Performance optimization for 10K+ spots | Performance | Medium |

### 15. mcp-tcga (Mocked)

| Enhancement | Type | Effort |
|-------------|------|--------|
| Wire up real GDC Data Portal API | Feature | Large |
| Implement retry/circuit-breaker for GDC API | Reliability | Small |
| Add `get_methylation_data` tool | Feature | Medium |
| Add `get_copy_number_data` tool (GISTIC2) | Feature | Medium |
| Add mutation co-occurrence analysis | Feature | Medium |
| Add CPTAC proteomics integration | Feature | Large |

---

## Cross-Cutting Enhancements

| Enhancement | Servers Affected | Effort |
|-------------|-----------------|--------|
| Retry/backoff + circuit breaker for all external APIs | epic, fgbio, seqera, tcga | Medium |
| FHIR R4 compliance for clinical data exchange | epic, mockepic, patient-report, genomic-results | Large |
| GPU acceleration | quantum, perturbation, deepcell | Large |
| Authentication on Cloud Run endpoints | perturbation, seqera, quantum | Small |
| DICOM image format support | openimagedata, deepcell | Large |

---

**See also:** [Architecture Overview](README.md) | [Server Registry](../shared/server-registry.md)
