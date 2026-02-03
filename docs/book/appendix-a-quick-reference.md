# Appendix A: Quick Reference Guides

*Essential reference for MCP servers, tools, and prompts*

---

## MCP Server Quick Reference

### Production Servers (7/12 - 58%)

| Server | Tools | Status | Key Capabilities |
|--------|-------|--------|------------------|
| **mcp-spatialtools** | 23 | ✅ Production | STAR alignment, ComBat, Moran's I, pathway enrichment |
| **mcp-multiomics** | 21 | ✅ Production | HAllA, Stouffer meta-analysis, multi-omics integration |
| **mcp-fgbio** | 9 | ✅ Production | VCF parsing, variant annotation, quality metrics |
| **mcp-deepcell** | 4 | ✅ Production | Nuclear/membrane segmentation, phenotype classification |
| **mcp-perturbation** | 8 | ✅ Production | GEARS GNN, drug response prediction |
| **mcp-quantum-celltype-fidelity** | 6 | ✅ Production | Quantum PQCs, Bayesian UQ, cell-type fidelity |
| **mcp-mockepic** | 12 | ✅ Production | Synthetic FHIR data for testing |
| **mcp-epic** | 12 | ⚠️ Partial | FHIR R4 integration (requires Epic credentials) |
| **mcp-openimagedata** | 9 | ❌ Mocked | Histopathology imaging (framework) |
| **mcp-tcga** | 7 | ❌ Mocked | TCGA cohort comparisons (framework) |
| **mcp-huggingface** | 6 | ❌ Mocked | ML model inference (framework) |
| **mcp-seqera** | 7 | ❌ Mocked | Nextflow orchestration (framework) |

**Total**: 12 servers, 124 tools

**Full registry**: [`SERVER_REGISTRY.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/SERVER_REGISTRY.md)

---

## Top 20 Most-Used Tools

### Clinical Data
1. **mockepic_get_patient_demographics** - Retrieve patient info (age, sex, diagnosis)
2. **mockepic_get_genomic_test_order** - Get ordered genomic tests
3. **mockepic_get_medications** - Current medication list

### Genomics
4. **fgbio_parse_vcf** - Parse VCF files, extract variants
5. **fgbio_annotate_variants** - Add gene annotations (ExAC, ClinVar)
6. **fgbio_quality_metrics** - Calculate QC metrics

### Multi-Omics
7. **multiomics_run_halla** - HAllA association discovery (RNA-protein)
8. **multiomics_stouffer_meta_analysis** - Cross-omics meta-analysis
9. **multiomics_load_multiomics_data** - Load RNA/protein/phospho data

### Spatial Transcriptomics
10. **spatialtools_load_visium** - Load 10X Visium data
11. **spatialtools_morans_i** - Spatial autocorrelation
12. **spatialtools_combat_correction** - Batch effect removal
13. **spatialtools_differential_expression** - Spatial DE analysis
14. **spatialtools_pathway_enrichment** - Pathway analysis (Reactome, KEGG)

### Cell Segmentation
15. **deepcell_segment_cells** - Nuclear/membrane segmentation
16. **deepcell_classify_cell_states** - Proliferating vs quiescent

### Treatment Prediction
17. **perturbation_predict_response** - Drug response prediction (GEARS)
18. **perturbation_train_model** - Train GEARS model

### Quantum Fidelity
19. **learn_spatial_cell_embeddings** - Train quantum embeddings
20. **compute_cell_type_fidelity** - Compute fidelity with Bayesian UQ

**Full tool documentation**: Each server's README.md in [`servers/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/servers)

---

## Prompt Template Library

### Template 1: Complete Patient Analysis

```
Analyze patient {PATIENT_ID} using all available modalities:

1. Clinical: Retrieve demographics, diagnoses, medications, genomic test orders
2. Genomics: Parse VCF, annotate pathogenic variants (ClinVar)
3. Multi-omics: Load RNA/protein data, run HAllA associations, Stouffer meta-analysis
4. Spatial: Load Visium data, calculate Moran's I for {GENE_LIST}, identify tumor regions
5. Segmentation: Segment cells, classify proliferating vs quiescent
6. Treatment: Predict response to {DRUG_LIST}

Provide integrated clinical summary with treatment recommendations.
```

**Variables**: `{PATIENT_ID}`, `{GENE_LIST}`, `{DRUG_LIST}`

### Template 2: Genomic Variant Analysis

```
Analyze genomic variants for patient {PATIENT_ID}:

1. Parse VCF file: {VCF_PATH}
2. Annotate with ClinVar, ExAC, COSMIC
3. Filter for pathogenic variants (ClinVar significance 4+)
4. Identify actionable mutations with FDA-approved therapies
5. Calculate tumor mutational burden (TMB)

Summarize key mutations and therapeutic implications.
```

**Variables**: `{PATIENT_ID}`, `{VCF_PATH}`

### Template 3: Spatial Transcriptomics Analysis

```
Analyze spatial transcriptomics for patient {PATIENT_ID}:

1. Load 10X Visium data: {H5_PATH}
2. Quality filter: min 200 UMIs, 100 genes per spot
3. ComBat batch correction (if multiple sections)
4. Calculate Moran's I for: {MARKER_GENES}
5. Identify spatial domains via clustering
6. Run pathway enrichment (Reactome) on tumor regions

Provide spatial map with annotated regions.
```

**Variables**: `{PATIENT_ID}`, `{H5_PATH}`, `{MARKER_GENES}`

### Template 4: Drug Response Prediction

```
Predict drug response for patient {PATIENT_ID}:

1. Load patient baseline scRNA-seq: {ADATA_PATH}
2. Train GEARS model on {DATASET_ID} (20 epochs)
3. Predict response to: {DRUG_LIST}
4. Rank treatments by efficacy score
5. Provide confidence intervals (Bayesian UQ)

Recommend top treatment with rationale.
```

**Variables**: `{PATIENT_ID}`, `{ADATA_PATH}`, `{DATASET_ID}`, `{DRUG_LIST}`

### Template 5: Quantum Cell-Type Verification

```
Verify cell-type classifications with quantum fidelity:

1. Load spatial data: {ADATA_PATH}
2. Train quantum embeddings (8 qubits, 3 layers, 20 epochs)
3. Compute fidelity with Bayesian UQ (100 samples)
4. Identify immune evasion states (fidelity threshold 0.3)
5. Analyze TLS quantum signatures

Report high-confidence immune evading cells.
```

**Variables**: `{ADATA_PATH}`

**Full prompt library**: [`docs/prompt-library/PROMPT_INVENTORY.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/prompt-library/PROMPT_INVENTORY.md)

---

## Common Error Solutions

### Error: "ModuleNotFoundError: No module named 'mcp_spatialtools'"

**Cause**: Server dependencies not installed or venv not activated

**Solution**:
```bash
cd servers/mcp-spatialtools
python -m venv venv
source venv/bin/activate
pip install -e ".[dev]"
```

### Error: "FileNotFoundError: PatientOne data not found"

**Cause**: PatientOne dataset not downloaded

**Solution**:
```bash
# Download from GCS
gsutil -m cp -r gs://precision-medicine-mcp-public/patient-data/PAT001-OVC-2025/ ./data/

# Or clone repository with data
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
```

### Error: "TensorFlow incompatible with Python 3.11+"

**Cause**: DeepCell requires TensorFlow 2.8.x (Python 3.10 only)

**Solution**:
```bash
# Install Python 3.10 for mcp-deepcell
brew install python@3.10  # macOS
cd servers/mcp-deepcell
python3.10 -m venv venv
source venv/bin/activate
pip install -e ".[dev]"
```

Full story: [`servers/mcp-deepcell/DEPENDENCY_ISSUES.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/DEPENDENCY_ISSUES.md)

### Error: "STAR index not found"

**Cause**: STAR aligner genome index not built

**Solution**:
```bash
# Download reference genome (hg38)
wget http://ftp.ensembl.org/pub/release-104/fasta/\
  homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Build STAR index
STAR --runMode genomeGenerate \
     --genomeDir /data/star_index \
     --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa \
     --sjdbGTFfile Homo_sapiens.GRCh38.104.gtf \
     --runThreadN 8
```

See: [`servers/mcp-spatialtools/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-spatialtools/README.md#star-alignment-setup)

### Error: "Cloud Run deployment timeout"

**Cause**: Build takes longer than default 10-minute timeout

**Solution**:
```bash
# Edit cloudbuild.yaml, set timeout: 1200s
gcloud builds submit --timeout=20m
```

Full troubleshooting: **Chapter 12** (pages 175-178)

---

## API Endpoint Reference

### SSE Transport (Cloud Run)

**Base URL format**: `https://mcp-{SERVER_NAME}-{PROJECT_ID}.run.app/sse`

**Example endpoints**:
- Spatial: `https://mcp-spatialtools-my-project.run.app/sse`
- Multi-omics: `https://mcp-multiomics-my-project.run.app/sse`
- DeepCell: `https://mcp-deepcell-my-project.run.app/sse`

**Usage with Claude API**:
```python
import anthropic

client = anthropic.Anthropic()

response = client.beta.messages.create(
    model="claude-sonnet-4-5",
    max_tokens=2048,
    messages=[{"role": "user", "content": "Analyze PatientOne spatial data"}],
    mcp_servers=[{
        "type": "url",
        "url": "https://mcp-spatialtools-PROJECT_ID.run.app/sse",
        "name": "spatialtools"
    }],
    tools=[{"type": "mcp_toolset", "mcp_server_name": "spatialtools"}],
    betas=["mcp-client-2025-11-20"]
)
```

Full API documentation: **Chapter 12** (pages 167-179)

---

## PatientOne Quick Facts

**Patient ID**: PAT001-OVC-2025
**Diagnosis**: Stage IV high-grade serous ovarian carcinoma
**Age**: 58 years
**Status**: 100% synthetic (safe to share, no IRB required)

**Key mutations**:
- TP53 R175H (pathogenic, 85% VAF)
- BRCA1 wild-type
- PIK3CA E545K (pathogenic, 42% VAF)

**Data files** (17 total):
- Clinical: 2 FHIR R4 resources
- Genomics: 1 VCF (8 pathogenic variants)
- Multi-omics: 4 files (RNA/protein/phospho, 15 PDX models)
- Spatial: 3 files (10X Visium, 900 spots, 6 regions)
- Imaging: 7 files (H&E, MxIF)

**Treatment recommendations**:
1. Olaparib (PARP inhibitor) - 82% predicted efficacy
2. Carboplatin + olaparib - 71% predicted efficacy
3. Carboplatin alone - 45% predicted efficacy

Full dataset: [`data/patient-data/PAT001-OVC-2025/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/data/patient-data/PAT001-OVC-2025/README.md)

---

## Performance Benchmarks

### Analysis Time

| Workflow | Traditional | AI-Orchestrated | Speedup |
|----------|-------------|-----------------|---------|
| Complete patient analysis | 40 hours | 35 minutes | 68x |
| Genomic QC + annotation | 10 hours | 8 minutes | 75x |
| Multi-omics integration | 8 hours | 5 minutes | 96x |
| Spatial transcriptomics | 12 hours | 12 minutes | 60x |
| Cell segmentation | 4 hours | 3 minutes | 80x |
| Treatment prediction | 6 hours | 7 minutes | 51x |

### Cost per Analysis

| Deployment | Cost per Patient | Components |
|------------|------------------|------------|
| **AI-orchestrated (Cloud Run)** | $1.20-2.00 | Cloud Run: $0.02-0.21<br/>Claude API: $0.50-1.00<br/>Gemini API: $0.30-0.80 |
| **Traditional manual** | $3,200 | Personnel time: $200/hr × 16 hrs |
| **Savings** | **95% reduction** | $3,198 per patient |

### Server Resource Usage

| Server | CPU | RAM | Storage | GPU |
|--------|-----|-----|---------|-----|
| mcp-spatialtools | 2 vCPU | 4Gi | 10GB | No |
| mcp-multiomics | 2 vCPU | 4Gi | 5GB | No |
| mcp-deepcell | 2 vCPU | 4Gi | 15GB | Optional |
| mcp-perturbation | 2 vCPU | 4Gi | 10GB | Optional |
| mcp-quantum-celltype-fidelity | 1 vCPU | 2Gi | 2GB | No |

Full benchmarks: **Chapter 3** (pages 30-44)

---

## Additional Resources

### Key Documentation

- **Architecture overview**: [`docs/architecture/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/architecture/README.md)
- **Server registry**: [`docs/SERVER_REGISTRY.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/SERVER_REGISTRY.md)
- **Prompt library**: [`docs/prompt-library/PROMPT_INVENTORY.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/prompt-library/PROMPT_INVENTORY.md)
- **Quick reference**: [`docs/for-developers/QUICK_REFERENCE.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-developers/QUICK_REFERENCE.md)

### External Resources

- **MCP Specification**: [modelcontextprotocol.io](https://modelcontextprotocol.io)
- **FastMCP Framework**: [github.com/jlowin/fastmcp](https://github.com/jlowin/fastmcp)
- **Claude Desktop**: [claude.com/claude-desktop](https://claude.com/claude-desktop)
- **Claude API Docs**: [docs.anthropic.com](https://docs.anthropic.com)

---

**This appendix provides quick reference for common tasks. For detailed implementation guides, see the main chapters and linked repository documentation.**

---
