# For Researchers

*Exploratory analysis, prompt engineering, and research workflows*

---

## Why Researchers Need This

Chapters 1-14 focused on clinical deployment. **Researchers have different needs**:

**Clinical workflow**: Standardized PatientOne analysis → treatment recommendations
**Research workflow**: Exploratory analysis → hypothesis generation → validation → publication

**Researcher requirements**:
- **Exploratory analysis**: "Show me all spatially variable genes" (no predefined list)
- **Hypothesis testing**: "Is HIF1A expression correlated with platinum resistance?"
- **Reproducibility**: Methods sections for publications
- **Cost-effectiveness**: $25-104 per patient (vs $6,000 traditional)
- **Data sharing**: Export results for collaborators

**PatientOne for research**: Synthetic dataset (100% safe to publish, share, teach).

---

## Research Use Cases

### 1. Tumor Microenvironment Characterization

**Research question**: How do spatial patterns of immune cells relate to treatment resistance?

**Workflow**:
```
1. Load spatial transcriptomics (Visium 900 spots × 31 genes)
2. Cell type deconvolution (tumor, fibroblasts, immune, hypoxic)
3. Spatial neighborhood analysis (immune-excluded vs infiltrated regions)
4. Pathway enrichment by region (tumor vs stroma vs immune)
5. Correlation with clinical outcomes
```

**Example prompt**:
```
Perform cell type deconvolution on PatientOne spatial data.
Identify spatially variable genes using Moran's I (p < 0.05).
Correlate immune infiltration with spatial pathway enrichment scores.
Create visualization showing immune exclusion zones.
```

**Expected analysis time**: 10-15 minutes
**Cost**: ~$0.50 (Claude API) + $0.02 (Cloud Run) = $0.52

**Publications enabled**: Spatial heterogeneity studies, immune infiltration patterns, treatment response prediction.

Full workflow: [`docs/for-researchers/README.md#tumor-microenvironment-characterization`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-researchers/README.md#tumor-microenvironment-characterization)

### 2. Drug Resistance Mechanisms

**Research question**: Which pathways are activated in platinum-resistant tumors?

**Workflow**:
```
1. Load multi-omics (RNA + protein + phospho from 15 samples)
2. Stratify by response (responders vs non-responders)
3. Differential expression across all modalities
4. Stouffer meta-analysis (combine p-values across modalities)
5. Pathway enrichment on concordant hits
6. Validate with genomics (variant-pathway mapping)
```

**Example prompt**:
```
Integrate PatientOne RNA, protein, and phospho data using Stouffer's method.
Identify pathways activated concordantly across all 3 modalities (FDR < 0.05).
Map to drug targets with FDA-approved therapies.
```

**Expected analysis time**: 15-25 minutes
**Cost**: ~$0.75

**Publications enabled**: Resistance biomarker discovery, mechanism-of-action studies, combination therapy rationale.

Full workflow: [`docs/for-researchers/README.md#drug-resistance-mechanisms`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-researchers/README.md#drug-resistance-mechanisms)

### 3. Biomarker Discovery & Validation

**Research question**: Can we identify a prognostic signature for ovarian cancer?

**Workflow**:
```
Discovery cohort (PatientOne + synthetic):
1. Feature selection (differential expression, log2FC > 1)
2. Pathway enrichment (biological relevance)
3. Candidate biomarkers (top genes/pathways)

Validation cohort (TCGA):
4. Load TCGA ovarian cancer cohort
5. Test biomarkers in independent dataset
6. Clinical correlation (link to survival, response)
```

**Example prompt**:
```
Identify top 20 differentially expressed genes in PatientOne tumor vs normal
(Mann-Whitney U test, FDR < 0.05, log2FC > 1).

For each gene, check expression in TCGA ovarian cancer cohort.
Correlate with overall survival and platinum response status.
```

**Expected analysis time**: 20-30 minutes (if TCGA server available)
**Cost**: ~$1.00

**Publications enabled**: Biomarker validation studies, prognostic signature development, clinical utility assessment.

Full workflow: [`docs/for-researchers/README.md#biomarker-discovery--validation`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-researchers/README.md#biomarker-discovery--validation)

### Research Interface

![Streamlit UI Preview](streamlit-ui-preview.png){width=100%}

**Figure 15.1: Streamlit Research Interface**
*Web-based chat interface for exploratory analysis. Researchers can interact with all 12 MCP servers through natural language, visualize results, and export data for publications. Supports both Claude and Gemini AI models.*

---

## Prompt Engineering Patterns

**Effective prompts** follow these patterns:

### Pattern 1: Be Specific

**❌ Vague**:
```
Analyze PatientOne data
```

**✅ Specific**:
```
Perform spatial pathway enrichment on PatientOne tumor regions,
focusing on cancer-related KEGG pathways with FDR < 0.05.
Return top 10 pathways with p-values and gene lists.
```

### Pattern 2: Include Parameters

**❌ Missing parameters**:
```
Run differential expression analysis
```

**✅ With parameters**:
```
Run differential expression analysis comparing PatientOne tumor vs normal samples,
using Mann-Whitney U test with FDR < 0.05 threshold and log2FC > 1.
```

### Pattern 3: Chain Multiple Steps

**❌ Single step**:
```
Load spatial data
```

**✅ Multi-step workflow**:
```
For PatientOne spatial data:
1. Load Visium dataset and summarize (spots, genes, regions)
2. Run spatial differential expression (tumor vs normal, FDR < 0.05)
3. Perform pathway enrichment on upregulated genes (KEGG pathways)
4. Create spatial visualization showing top pathway expression
```

### Pattern 4: Specify Expected Output

**❌ No output guidance**:
```
Find upregulated genes
```

**✅ With output guidance**:
```
Find upregulated genes in PatientOne tumor regions (log2FC > 1, FDR < 0.05).
Return top 20 genes ranked by fold change.
Create volcano plot showing all genes with significant genes highlighted.
```

Full prompt engineering guide: [`docs/prompt-library/README.md#prompt-best-practices`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/prompt-library/README.md#prompt-best-practices)

---

## Prompt Library (90+ Prompts)

**Curated prompts** organized by audience and complexity:

### Clinical & Genomic Prompts

```
"Load PatientOne VCF file and identify pathogenic variants
(TP53, BRCA1, PIK3CA, PTEN). Interpret clinical significance."

"Get patient demographics and treatment history from Epic FHIR.
Extract CA-125 trend over time. Summarize response to prior therapy."
```

**Servers**: mcp-epic, mcp-fgbio
**Time**: 5-10 minutes

### Multi-Omics Prompts

```
"Integrate PatientOne RNA, protein, and phospho data using Stouffer meta-analysis.
Identify concordant pathway activations (FDR < 0.05 across all modalities).
Return top 5 pathways with combined p-values."

"Run HAllA association analysis on PatientOne multi-omics data.
Identify cross-modal associations between RNA and protein (FDR < 0.05).
```

**Servers**: mcp-multiomics
**Time**: 15-25 minutes

### Spatial Transcriptomics Prompts

```
"Load PatientOne Visium data (900 spots × 31 genes).
Perform spatial pathway enrichment focusing on tumor regions.
Identify spatially variable genes using Moran's I (p < 0.05).
Create visualization showing top pathway spatial distribution."

"Apply ComBat batch correction to PatientOne spatial data.
Create PCA plots before and after correction.
Verify batch effect removal."
```

**Servers**: mcp-spatialtools
**Time**: 10-15 minutes

### Complete Workflows

```
"Perform comprehensive multi-modal analysis for PatientOne:
1. Load clinical data (demographics, diagnoses, medications)
2. Analyze genomic variants (TP53, BRCA1 status)
3. Integrate multi-omics (RNA, protein, phospho using Stouffer)
4. Analyze spatial transcriptomics (pathway enrichment)
5. Synthesize results into treatment recommendations"
```

**Servers**: All 12 servers (124 tools)
**Time**: 35 minutes

Full prompt library: [`docs/prompt-library/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/docs/prompt-library) (90+ prompts)

Prompt index: [`docs/prompt-library/INDEX.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/prompt-library/INDEX.md)

---

## Statistical Methods (Reproducibility)

**All analyses use peer-reviewed methods**:

| Analysis | Method | Multiple Testing | Reference |
|----------|--------|------------------|-----------|
| Differential expression | Mann-Whitney U | Benjamini-Hochberg FDR | Non-parametric test |
| Pathway enrichment | Fisher's exact test | FDR < 0.05 | Hypergeometric distribution |
| Spatial autocorrelation | Moran's I | FDR < 0.05 | Spatial statistics |
| Batch correction | ComBat | N/A | Empirical Bayes |
| Meta-analysis | Stouffer's Z-score | FDR after combination | P-value combination |

**Methods section template** (for publications):

```markdown
Spatial pathway enrichment was performed using mcp-spatialtools
(version 0.3.0) with Fisher's exact test on 44 curated pathways
(KEGG, Hallmark, GO_BP, Drug_Resistance). FDR correction was
applied using the Benjamini-Hochberg method with α = 0.05.
Spatial graphs were constructed using k=6 nearest neighbors.
```

**Reproducibility features**:
- Tool versions logged (server commits, library versions)
- Parameters tracked (thresholds, methods, corrections)
- Data provenance (file paths, checksums)
- Random seeds (where applicable)

Full statistical methods: [`docs/for-researchers/README.md#statistical-methods`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-researchers/README.md#statistical-methods)

---

## Cost Analysis for Researchers

### Per-Patient Analysis Costs

| Analysis Type | Compute | API Tokens | Total | Traditional |
|--------------|---------|------------|-------|-------------|
| **Demo (DRY_RUN)** | ~$0 | ~$0.32 | **$0.32** | N/A |
| **Small files (PatientOne)** | $24-48 | $1-2 | **$25-50** | $6,000 |
| **Production (large files)** | $90-100 | $2-4 | **$92-104** | $6,000-9,000 |

**Savings**: $5,896-8,908 per patient (98% cost reduction)

### Cohort Analysis Costs (100 patients)

| Component | Cost | Notes |
|-----------|------|-------|
| Per-patient analysis | $25-104 × 100 | $2,500-10,400 total |
| Infrastructure | $1,000/month | GCP Cloud Run, storage |
| **Total (annual)** | **$14,500-22,400** | vs $600,000 traditional |

**Annual savings**: $577,600-585,500 (96% cost reduction)

### Grant Budget Example (NIH R01)

**Budget line items**:
- Infrastructure: $12,000/year (Cloud Run + storage)
- Claude API: $2,400/year (100 patients × $24 avg)
- Personnel: $50,000/year (1 bioinformatician @ 0.5 FTE)
- **Total**: $64,400/year

**Traditional alternative**:
- Personnel: $120,000/year (1.5 FTE @ 100% manual analysis)
- Infrastructure: $10,000/year (on-premise compute)
- **Total**: $130,000/year

**Savings**: $65,600/year (50% reduction)

Full cost analysis: [`docs/for-researchers/README.md#estimated-cost-analysis`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-researchers/README.md#estimated-cost-analysis)

---

## PatientOne Synthetic Dataset

**100% synthetic** = safe to publish, share, teach.

**Data availability**:
- **Location**: [`data/patient-data/PAT001-OVC-2025/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/data/patient-data/PAT001-OVC-2025)
- **Format**: FHIR JSON, VCF, CSV matrices, 10X Visium, TIFF images
- **License**: CC0 1.0 Universal (Public Domain)
- **DOI**: [To be assigned upon publication]

**Data modalities**:

| Modality | Demonstration | Production |
|----------|--------------|------------|
| Clinical | FHIR resources (demographics, CA-125) | Real Epic FHIR (HIPAA de-identified) |
| Genomics | VCF: TP53, PIK3CA, PTEN, BRCA1 variants | Whole exome sequencing (WES) |
| Multi-omics | 15 samples, 38 KB matrices | 15 samples, 2.7 GB raw |
| Spatial | 900 spots × 31 genes (315 KB) | 3,000-5,000 spots × 18,000 genes (100-500 MB) |
| Imaging | H&E, MxIF placeholders (4.1 MB) | Full resolution slides (500 MB - 2 GB) |

**Use in publications**:
- No privacy concerns (100% synthetic)
- No IRB approval needed (not human subjects research)
- Safe to share in supplementary materials
- Ideal for methods papers, software validation

Full data guide: [`data/patient-data/PAT001-OVC-2025/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/data/patient-data/PAT001-OVC-2025)

---

## Server Implementation Status

**4/12 servers production-ready** for research:

### Production (Real Analysis)
1. **mcp-fgbio** (95% real): Genomic QC, VCF parsing, variant annotation
2. **mcp-multiomics** (85% real): HAllA integration, Stouffer meta-analysis
3. **mcp-spatialtools** (95% real): STAR alignment, pathway enrichment, Moran's I
4. **mcp-epic** (90% real): FHIR integration, HIPAA de-identification

### Partial (30-60% real)
5. **mcp-openimagedata** (60% real): Image loading, registration (feature extraction mocked)
6. **mcp-quantum-celltype-fidelity** (40% real): Bayesian UQ (PQC mocked)

### Mocked (Framework Only)
7. **mcp-deepcell**: Segmentation framework (DeepCell integration planned)
8. **mcp-perturbation**: Treatment prediction framework (GEARS integration planned)
9. **mcp-tcga**: TCGA query framework (GDC API integration planned)
10. **mcp-huggingface**: ML framework (model integration planned)
11. **mcp-seqera**: Workflow framework (Nextflow integration planned)
12. **mcp-mockepic**: Synthetic FHIR only (for testing)

**Roadmap**: Chapters 8-11 documented advanced servers (DeepCell, GEARS, Quantum, Imaging). Production integration: 6-12 months.

Full server status: [`docs/architecture/servers.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/architecture/servers.md)

---

## Example Research Workflow

**Research question**: Is immune exclusion associated with platinum resistance in PatientOne?

**Hypothesis**: Spatially excluded immune cells correlate with activated PI3K/AKT pathway (resistance mechanism).

**Analysis**:

```
1. Spatial analysis:
"Load PatientOne Visium data. Perform cell type deconvolution
(tumor, immune, fibroblasts). Calculate immune cell fraction per spot.
Identify immune-excluded regions (immune fraction < 10%)."

Result: 35% of tumor spots are immune-excluded.

2. Pathway analysis:
"Perform spatial pathway enrichment on immune-excluded regions.
Focus on PI3K/AKT and drug resistance pathways. Compare to
immune-infiltrated regions."

Result: PI3K/AKT 3.2× higher in immune-excluded vs infiltrated (p = 0.003).

3. Multi-omics validation:
"Integrate PatientOne multi-omics data using Stouffer meta-analysis.
Test PI3K/AKT activation across RNA, protein, phospho modalities."

Result: PI3K/AKT activated concordantly across all 3 modalities (combined p = 2.1e-5).

4. Clinical correlation:
"Extract CA-125 trend from Epic FHIR. Correlate with PI3K/AKT activation score."

Result: High PI3K/AKT correlates with rising CA-125 (ρ = 0.78, p = 0.01).
```

**Total analysis time**: 25 minutes
**Total cost**: $0.85

**Conclusion**: Immune exclusion + PI3K/AKT activation explains platinum resistance. **Treatment**: PI3K inhibitor (alpelisib) + immune checkpoint blockade (pembrolizumab).

**Publication-ready**: Methods section, reproducible workflow, synthetic data (safe to share).

---

## Frequently Asked Questions

**Q: Can I use this for my research publication?**
**A**: Yes! Platform designed for research use. Synthetic PatientOne data is 100% safe to publish. For real patient data, ensure IRB approval.

**Q: How do I cite this platform?**
**A**: Citation information will be provided upon publication. For now, reference GitHub repository and specific tool versions used.

**Q: Can I add custom pathways or gene signatures?**
**A**: Yes! See [`docs/for-developers/ADD_NEW_MODALITY_SERVER.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-developers/ADD_NEW_MODALITY_SERVER.md) for how to add custom analysis tools.

**Q: What data formats are supported?**
**A**:
- Clinical: FHIR JSON
- Genomics: VCF, BAM, FASTQ
- Multi-omics: CSV matrices (samples × features)
- Spatial: 10X Visium format, Seurat objects
- Imaging: TIFF, PNG, DICOM

**Q: How do I ensure reproducibility?**
**A**: Platform automatically tracks tool versions, parameters, data provenance, and random seeds. Export logs for methods sections.

**Q: Can I share my results with collaborators?**
**A**: Yes! Export results as CSV, JSON, or images. PatientOne synthetic data is safe to share without restrictions.

Full FAQ: [`docs/for-researchers/README.md#frequently-asked-questions`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-researchers/README.md#frequently-asked-questions)

---

## What You've Learned

**Research workflows**:
1. **Exploratory analysis**: Hypothesis generation from multi-modal data
2. **Prompt engineering**: Specific, parameterized, multi-step prompts
3. **Statistical rigor**: Peer-reviewed methods, multiple testing correction
4. **Cost-effectiveness**: $25-104 per patient vs $6,000 traditional (98% savings)
5. **Reproducibility**: Tool versions, parameters, data provenance tracked
6. **Publishing**: Synthetic PatientOne data safe to share, methods sections provided

**Research use cases**:
- Tumor microenvironment characterization ($0.52, 10-15 min)
- Drug resistance mechanisms ($0.75, 15-25 min)
- Biomarker discovery & validation ($1.00, 20-30 min)

**Resources**:
- 90+ curated prompts in prompt library
- 100% synthetic PatientOne dataset (safe to publish)
- 4/12 servers production-ready (real analysis)
- Statistical methods documentation

---

## Summary

**Chapter 15 Summary**:
- Researchers need exploratory workflows (not just standardized clinical)
- Prompt engineering: Be specific, include parameters, chain steps, specify output
- 90+ prompts in library (clinical, multiomics, spatial, workflows)
- Cost: $25-104 per patient (98% savings), $14.5K-22.4K/year for 100-patient cohort
- Reproducibility: Tool versions, parameters, provenance tracked
- PatientOne: 100% synthetic, safe to publish, CC0 license
- 4/12 servers production-ready for real research

**Files**: [`docs/for-researchers/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/docs/for-researchers), [`docs/prompt-library/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/docs/prompt-library)
**Research use cases**: Tumor microenvironment, drug resistance, biomarker discovery
**Cost**: $0.32 (demo) to $104 (production) per patient
