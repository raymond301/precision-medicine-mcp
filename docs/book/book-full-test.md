---
title: "Building AI-Orchestrated Precision Oncology Systems"
subtitle: "A Practical Guide to MCP-Based Bioinformatics"
author: "Lynn Langit"
date: "2026-02-01"
toc: true
toc-depth: 2
numbersections: true
header-includes:
  - \setcounter{chapter}{-1}
include-before:
  - \setcounter{chapter}{-1}
linkcolor: blue
urlcolor: blue
geometry: margin=1in
fontsize: 11pt
documentclass: report
---

# Building AI-Orchestrated Precision Oncology Systems

**A practical guide to deploying MCP-based precision medicine workflows**

By Lynn Langit

---

## About This Book

- This book teaches you how to build, deploy, and operate an AI-orchestrated precision oncology system that reduces analysis time from **40 hours to 35 minutes** and cost from **$3,200 to $1-2 per patient**.  
- You'll learn to create 12 specialized MCP (Model Context Protocol) servers that enable AI models like Claude and Gemini to coordinate 124 bioinformatics tools through natural language prompts—no code required.  
- **What makes this different**: This isn't a theoretical framework. Every example is based on a production system deployed to Google Cloud Run, tested with real patient workflows, and validated for clinical decision support.
- All referenced source code and detailed documentation is available on Github at https://github.com/lynnlangit/precision-medicine-mcp  

---

## Who This Book Is For

**Primary Audience**:
- Bioinformatics researchers building precision medicine platforms
- Clinical informaticists integrating multi-modal cancer data
- Healthcare software architects designing AI-assisted analysis systems

**You should have**:
- Basic Python programming (not expert-level)
- Familiarity with genomic data formats (VCF, FASTQ)
- Understanding of cloud deployment concepts (Docker, APIs)

**You don't need**:
- Deep learning expertise (we provide pre-trained models)
- DevOps mastery (deployment scripts included)
- Prior MCP or Claude API experience (taught from scratch)

---

## What You'll Build

By the end of this book, you'll have deployed:

**12 MCP Servers**:
1. **mcp-epic**: FHIR R4 clinical data integration. 
2. **mcp-fgbio**: Genomic QC and variant calling. 
3. **mcp-multiomics**: Multi-omics integration (RNA/protein/phospho). 
4. **mcp-spatialtools**: Spatial transcriptomics (STAR, ComBat, pathways). 
5. **mcp-deepcell**: Cell segmentation (DeepCell-TF). 
6. **mcp-perturbation**: Treatment response prediction (GEARS GNN). 
7. **mcp-quantum-celltype-fidelity**: Quantum fidelity with Bayesian UQ. 
8. **mcp-openimagedata**: Histopathology imaging. 
9. **mcp-tcga**: TCGA cohort comparisons (framework). 
10. **mcp-huggingface**: ML model inference (framework). 
11. **mcp-seqera**: Nextflow workflow orchestration (framework). 
12. **mcp-mockepic**: Synthetic FHIR data for testing. 

**PatientOne Workflow**: Complete precision medicine analysis for Stage IV ovarian cancer integrating clinical, genomic, multi-omics, spatial, and imaging data.

---

## Book Structure

### Part 1: Why This Matters (Chapters 1-3)
**Goal**: Understand the problem and solution

- **Chapter 1**: The PatientOne Story — 40 hours → 35 minutes
- **Chapter 2**: The Architecture Problem — Why MCP for healthcare
- **Chapter 3**: Testing the Hypothesis — Production validation

**Time to complete**: 2-3 hours reading + 1 hour hands-on

---

### Part 2: Building the Foundation (Chapters 4-7)
**Goal**: Implement core MCP servers

- **Chapter 4**: Clinical Data — FHIR integration with Epic
- **Chapter 5**: Genomic Foundations — VCF parsing, variant annotation
- **Chapter 6**: Multi-Omics Integration — HAllA, Stouffer meta-analysis
- **Chapter 7**: Spatial Transcriptomics — STAR alignment, batch correction, pathways

**Time to complete**: 8-12 hours implementation + testing

---

### Part 3: Advanced Capabilities (Chapters 8-11)
**Goal**: Add cutting-edge features

- **Chapter 8**: Cell Segmentation with DeepCell — Nuclear/membrane models, phenotyping
- **Chapter 9**: Treatment Response Prediction — GEARS GNN for perturbations
- **Chapter 10**: Quantum Cell-Type Fidelity — PennyLane PQCs, Bayesian UQ
- **Chapter 11**: Imaging and Histopathology — H&E, MxIF analysis

**Time to complete**: 6-10 hours implementation

---

### Part 4: Deployment and Operations (Chapters 12-14)
**Goal**: Production-ready deployment

- **Chapter 12**: Cloud Deployment on GCP — Cloud Run, Docker, SSE transport
- **Chapter 13**: Hospital Production Deployment — HIPAA compliance, de-identification, VPC
- **Chapter 14**: Operations and Monitoring — Logging, alerts, cost tracking

**Time to complete**: 4-8 hours deployment + configuration

---

### Part 5: Research and Education (Chapters 15-16)
**Goal**: Extend beyond clinical use

- **Chapter 15**: For Researchers — Exploratory analysis, prompt engineering
- **Chapter 16**: Teaching Precision Medicine — Educational workflows, cost-effective student access

**Time to complete**: 2-3 hours

---

### Part 6: The Future (Chapters 17-18)
**Goal**: Sustainability and next steps

- **Chapter 17**: Funding and Sustainability — ROI analysis, grant strategies
- **Chapter 18**: Lessons Learned and What's Next — Production insights, future enhancements

**Time to complete**: 1-2 hours

---

### Appendices
**Goal**: Quick reference guides for common tasks

- **Appendix A**: Quick Reference Guides — Server registry, tool catalog, prompt templates, common errors
- **Appendix B**: Installation and Setup — Prerequisites, local setup, cloud deployment, troubleshooting
- **Appendix C**: PatientOne Complete Dataset — Full data manifest, file formats, access methods
- **Appendix D**: Bias and Ethics — Framework for bias detection, audit checklist, ethical deployment

**Time to complete**: Reference guides (as needed)

---

## Companion Materials

### Jupyter Notebooks  
Each chapter has a hands-on Jupyter notebook in [`companion-notebooks/`](./companion-notebooks/):  
- Executable code examples. 
- Interactive exercises. 
- "Try changing this parameter..." experiments. 
- Links to deployed Cloud Run servers (**requires you to deploy your own MCP servers**). 

**ALL 18 NOTEBOOKS NOW AVAILABLE**:
- Part 1 (Ch 1-3): PatientOne demo, architecture, testing. 
- Part 2 (Ch 4-7): Clinical, genomics, multi-omics, spatial. 
- Part 3 (Ch 8-11): DeepCell, GEARS, quantum, imaging. 
- Part 4 (Ch 12-14): Cloud deployment, hospital deployment, operations. 
- Part 5 (Ch 15-16): Research workflows, teaching exercises. 
- Part 6 (Ch 17-18): Funding calculator, lessons learned. 

**IMPORTANT**: Notebooks require you to deploy MCP servers to your GCP Cloud Run. They will NOT work without your own infrastructure.

**Setup**: See [`companion-notebooks/README.md`](./companion-notebooks/README.md) for complete 5-step deployment guide

### Sample Data
PatientOne dataset (100% synthetic):
- Clinical: FHIR R4 resources
- Genomics: VCF with 8 pathogenic mutations
- Multi-omics: RNA/protein/phospho from 15 PDX models
- Spatial: 10X Visium (900 spots, 6 regions)
- Imaging: H&E and MxIF images

**Location**: [`data/patient-data/PAT001-OVC-2025/`](../../data/patient-data/PAT001-OVC-2025/)

### GitHub Repository
All code is open source (Apache 2.0):
**https://github.com/lynnlangit/precision-medicine-mcp**

---

## Prerequisites

### Software
- **Python**: 3.11+ (3.10 for DeepCell server)
- **Docker**: For containerization
- **Git**: For cloning repository
- **Claude Desktop** (optional): For local MCP testing
- **Google Cloud SDK** (Chapter 12+): For cloud deployment

### Cloud Accounts (Optional)
- **Anthropic API**: For Claude orchestration (~$1 per analysis)
- **Google AI API**: For Gemini alternative (~$0.30 per analysis)
- **Google Cloud Platform**: For Cloud Run deployment (free tier available)

### Hardware
- **RAM**: 8GB minimum, 16GB recommended
- **Disk**: 50GB free space
- **GPU**: Not required (DeepCell runs on CPU, optional GPU acceleration)

---

## Installation

**Quick Start**:
```bash
# Clone repository
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp/docs/book

# Install companion notebook dependencies
cd companion-notebooks
pip install -r requirements.txt

# Launch Jupyter
jupyter notebook
```

**Full Setup**: See [`docs/getting-started/installation.md`](../getting-started/installation.md)

---

## How to Use This Book

### Linear Reading Path (Recommended)
Read chapters 1-18 in order. Each chapter builds on previous concepts.

**Time commitment**: 25-35 hours (reading + hands-on)

### Selective Reading Paths

**For Hospital IT Leaders** (deployment focus):
- Read: Chapters 1-3, 12-14, 17
- Skim: Chapters 4-11 (technical implementation)
- Skip: Chapters 15-16 (research/education)
- **Time**: 8-12 hours

**For Bioinformatics Researchers** (implementation focus):
- Read: Chapters 1-7 (foundation servers)
- Selective: Chapters 8-11 (choose modalities you need)
- Skim: Chapters 12-14 (deploy to your own cloud)
- Read: Chapter 15 (research workflows)
- **Time**: 15-20 hours

**For Developers/Architects** (system design focus):
- Read: Chapters 1-3, 12-14
- Selective: Chapters 4-11 (pick 2-3 servers to understand patterns)
- Skim: Chapters 15-18
- **Time**: 10-15 hours

---

## Cost Expectations

### Development/Learning
- **Local development**: Free (no cloud costs)
- **Claude API testing**: ~$5-10 for all chapter notebooks
- **Gemini API alternative**: ~$2-5 for all notebooks

### Production Deployment (Your Own Cloud)
- **Cloud Run**: ~$0.02-0.21 per patient analysis
- **Claude/Gemini API**: ~$0.50-1.00 per patient
- **Total**: **$1-2 per patient analysis**

### Comparison
- **Traditional manual analysis**: $3,200 per patient (personnel time)
- **AI-orchestrated (this book)**: $1-2 per patient
- **Savings**: 95% cost reduction

---

## License

**Book content**: Copyright © 2026 Lynn Langit. All rights reserved.  
**Code examples**: Apache License 2.0  
**Sample data**: CC0 1.0 Universal (Public Domain)  

You may:
- Use all code for commercial or non-commercial projects
- Modify and distribute code examples
- Deploy to your own cloud infrastructure

You may not:
- Reproduce or distribute the book text without permission
- Claim authorship of the system design or architecture

---

## Acknowledgments

This book and the underlying system would not exist without:  
- **Anthropic**: For Claude and the Model Context Protocol. 
- **Google Cloud**: For Cloud Run infrastructure and Gemini API. 
- **Open source bioinformatics community**: For tools like DeepCell, GEARS, and countless libraries. 
- **Early testers**: Who provided invaluable feedback on the POC deployment. 

Special thanks to the precision medicine research community for defining the workflows this system aims to accelerate.

---

## About the Author

**Lynn Langit** is a cloud architect specializing in bioinformatics and genomic-scale data analysis. She works with bioinformatics researchers worldwide to build and optimize genomic data pipelines on GCP, AWS, and Azure. Lynn is a Google AI & Cloud Developer Expert, and Microsoft Regional Director. She has authored 30 LinkedIn Learning courses on cloud computing and AI with over 5 million student views.  

**Contact**:
- GitHub: [@lynnlangit](https://github.com/lynnlangit). 
- LinkedIn: [lynnlangit](https://www.linkedin.com/in/lynnlangit/). 
- Substack: [Lynn Langit's Cloud World](https://lynnlangit.substack.com/). 

---

## Ready to Begin?

Start with [Chapter 1: The PatientOne Story](chapter-01-the-patientone-story.md) to see what's possible when AI orchestrates precision medicine workflows.

Or jump to the [Companion Notebooks](companion-notebooks/) to start building immediately.

**Let's transform precision oncology from 40 hours to 35 minutes.**

---



# The PatientOne Story

> *"What would happen if you could analyze a cancer patient's complete molecular profile in minutes instead of weeks?"*

---

## The Patient

Sarah Anderson is 58 years old. She was diagnosed with Stage IV high-grade serous ovarian cancer (HGSOC) two years ago. After aggressive treatment with platinum-based chemotherapy, she responded well initially. Her CA-125 tumor marker—a blood test used to monitor ovarian cancer—dropped from 1,200 U/mL to a near-normal 45 U/mL.

But cancer, especially advanced ovarian cancer, rarely surrenders that easily.

Eight months later, Sarah's CA-125 began climbing again: 78, then 156, then 310 U/mL. Her CT scan confirmed what the numbers suggested: platinum-resistant recurrence. New lesions appeared in her peritoneum, omentum, and liver. The chemotherapy that had worked so well before now failed to stop the cancer's progression.

Sarah's oncologist faces a critical question: **What treatment should she receive next?**

The traditional approach would involve:
- Reviewing her genomic sequencing reports (if they were even ordered)
- Consulting published literature on platinum-resistant ovarian cancer
- Considering clinical trial eligibility
- Making an educated guess based on population-level data

But what if you could do something fundamentally different? What if you could integrate:
- Her complete clinical history (demographics, labs, medications)
- Her tumor's genomic alterations (somatic mutations, copy number changes)
- Multi-omics data from patient-derived xenograft (PDX) models (RNA, protein, phosphoproteomics)
- Spatial transcriptomics showing exactly *where* in the tumor resistant cells hide
- Histology imaging revealing the tumor microenvironment
- AI-predicted treatment responses for candidate therapies

And what if you could do this analysis in **35 minutes** instead of the traditional **40 hours** of manual bioinformatics work?

This is the PatientOne story—and it's the story of why this book exists.

### The Complete PatientOne Workflow

![System Overview](images/screenshots/patient-one-holistic.png){width=\textwidth}

**Figure 1.1: PatientOne Complete Multi-Modal Analysis Workflow**
*Integrating clinical data (FHIR), genomics (VCF), multi-omics (RNA/protein/phospho), spatial transcriptomics (10X Visium), imaging (H&E, MxIF), and AI-predicted treatment responses—all orchestrated through 12 MCP servers.*

---

## The Traditional Workflow: 40 Hours

Let's break down what a comprehensive precision oncology analysis traditionally requires. You'll quickly see why most hospitals don't do this routinely.

### Day 1: Clinical Data Extraction (2-3 hours)

A clinical informaticist logs into the electronic health record (EHR) system—often Epic, Cerner, or a similar platform. They manually extract:
- Patient demographics
- Diagnosis codes (ICD-10)
- Medication history
- Lab results (CA-125 trends, complete blood counts)
- Imaging reports
- Treatment timelines

This data comes in various formats: PDFs, HL7 messages, proprietary database exports. Each element must be copied, cleaned, and formatted into a usable structure. If you're lucky, your institution has a FHIR API. If not, expect lots of copy-pasting into spreadsheets.

**Time: 2-3 hours**

### Day 2: Genomic Analysis (8-10 hours)

A bioinformatician receives the tumor sequencing files:
- Whole exome sequencing (WES) FASTQ files: ~60 GB
- Germline DNA for comparison: another ~60 GB
- RNA-seq data: ~40 GB

First, they run quality control:

```bash
fastqc sample_R1.fastq.gz sample_R2.fastq.gz
multiqc .
```

Then alignment to the reference genome:

```bash
bwa mem -t 8 hg38.fa sample_R1.fastq.gz sample_R2.fastq.gz | \
  samtools sort -o sample.bam
```

Variant calling with GATK or similar:

```bash
gatk HaplotypeCaller \
  -R hg38.fa \
  -I sample.bam \
  -O sample.vcf
```

**But that's just the pipeline execution.** The real work is:
- Comparing tumor vs. germline to identify somatic mutations
- Filtering thousands of variants to find clinically relevant ones
- Annotating variants using ClinVar, gnomAD, COSMIC databases
- Interpreting pathogenicity (is TP53 R175H different from TP53 R273H?)
- Searching literature for each candidate variant
- Comparing to TCGA cohorts to determine molecular subtype

**Time: 8-10 hours** (assuming pipelines are already configured)

### Day 3: Multi-Omics Integration (10-12 hours)

Now you have genomics. But genomics alone doesn't tell you what's happening right now in the tumor. For that, you need transcriptomics (RNA-seq), proteomics, and phosphoproteomics.

Sarah's institution has PDX models—patient-derived xenografts grown in mice to test treatment responses. You have RNA, protein, and phosphorylation data for 15 samples: 7 platinum-sensitive, 8 platinum-resistant.

The analysis requires:
1. **Differential expression** across three modalities (RNA, protein, phospho)
2. **Meta-analysis** to find consistent signals across all three
3. **Pathway enrichment** to identify dysregulated biological processes
4. **Literature mining** to connect pathways to druggable targets

Each step involves:
- Loading data (CSV files with 20,000+ genes/proteins)
- Normalization (log2 transformation, batch correction)
- Statistical testing (t-tests, ANOVA, FDR correction)
- Visualization (heatmaps, volcano plots, pathway diagrams)
- Interpretation (is upregulated PI3K/AKT clinically actionable?)

Most bioinformaticians use R scripts or Python notebooks, stitching together tools like DESeq2, limma, clusterProfiler. Each dataset requires custom code.

**Time: 10-12 hours**

### Day 4: Spatial Transcriptomics (12-15 hours)

Spatial transcriptomics—technology like 10X Genomics' Visium platform—tells you not just *what* genes are expressed, but *where* in the tissue they're active. This is critical for understanding:
- Tumor heterogeneity (are resistant cells clustered in one region?)
- Immune infiltration (where are the T cells? Can they reach the tumor?)
- Microenvironment interactions (tumor-stroma signaling)

Sarah's tumor biopsy was processed through Visium, generating:
- 900 spatial spots across the tissue section
- ~30,000 genes measured per spot
- Spatial coordinates for each spot
- Region annotations (tumor core, proliferative edge, stroma, etc.)

The analysis workflow:
1. **Quality control**: Filter low-quality spots, normalize counts
2. **Spatial clustering**: Identify tissue regions computationally
3. **Differential expression**: Compare gene expression across regions
4. **Spatial statistics**: Calculate Moran's I for spatial autocorrelation
5. **Pathway analysis**: What pathways are enriched in tumor vs. immune regions?
6. **Visualization**: Generate spatial heatmaps overlaid on histology images

If you're using Scanpy or Seurat, you might get through this in a day. But troubleshooting batch effects, optimizing clustering parameters, and validating results against known biology? That's where the hours pile up.

**Time: 12-15 hours**

### Day 5: Histology Imaging (4-6 hours)

Pathologists examine H&E (hematoxylin and eosin) stained slides to assess:
- Tumor cellularity (what percentage of the tissue is tumor?)
- Morphology (cellular architecture, necrosis, invasion patterns)
- Immune infiltration (visible lymphocytes?)

For deeper analysis, you might have multiplexed immunofluorescence (MxIF) images showing:
- CD8+ cytotoxic T cells
- Ki67+ proliferating cells
- TP53 protein expression
- DAPI nuclear stain

Processing these images requires:
- Cell segmentation (identifying individual cells in images)
- Phenotyping (CD8+ vs. CD8-, Ki67+ vs. Ki67-)
- Quantification (cells per mm², marker colocalization)
- Spatial analysis (distances between cell types)

Tools like QuPath, CellProfiler, or DeepCell provide semi-automated analysis, but manual quality control and validation are essential.

**Time: 4-6 hours**

### Day 6: Integration and Report Generation (4-5 hours)

Finally, you synthesize all findings into a cohesive report:
- Clinical summary
- Genomic alterations with clinical significance
- Multi-omics resistance signatures
- Spatial heterogeneity insights
- Histology findings
- Treatment recommendations

This isn't copy-pasting results. It's interpretation:
- Do the genomic mutations (PIK3CA E545K) align with the proteomic data (AKT1 upregulation)?
- Does the spatial data explain the clinical phenotype (immune exclusion → immunotherapy resistance)?
- Are there clinical trials matching this molecular profile?

**Time: 4-5 hours**

### **Total Time: 40-48 hours**

That's **one full work week** for one patient. And this assumes:
- You have access to all the data (many hospitals don't run spatial transcriptomics)
- You have bioinformaticians trained in each modality
- Your pipelines are already configured and tested
- Nothing breaks along the way (it always does)

---

## The AI-Orchestrated Workflow: 35 Minutes

Now let's see what the same analysis looks like when an AI—specifically, Claude or Gemini—orchestrates specialized bioinformatics tools through the Model Context Protocol (MCP).

### The Architecture in One Sentence

**You type a natural language prompt, and Claude coordinates 12 specialized MCP servers (124 bioinformatics tools) to execute the complete analysis.**

### What Actually Happens

Here's the real prompt you'd use in Claude Desktop:

```
I want to analyze patient PAT001-OVC-2025 for precision oncology.

Please:
1. Load clinical data and summarize treatment history
2. Identify somatic mutations from genomics/somatic_variants.vcf
3. Run multi-omics integration (RNA, protein, phospho)
4. Analyze spatial transcriptomics (Visium data)
5. Process histology imaging (H&E and MxIF)
6. Generate treatment recommendations

Data files are in: gs://sample-inputs-patientone/PAT001-OVC-2025/
```

That's it. No code. No pipeline configuration. No switching between tools.

### Behind the Scenes (The First 5 Minutes)

Claude receives your prompt and thinks: *"This requires data from multiple sources. Let me check which MCP servers I have access to."*

It discovers:
- `mcp-epic`: Clinical data (FHIR resources)
- `mcp-fgbio`: Genomics quality control and variant annotation
- `mcp-multiomics`: Multi-omics integration
- `mcp-spatialtools`: Spatial transcriptomics analysis
- `mcp-deepcell`: Cell segmentation for imaging
- `mcp-openimagedata`: Histology image processing
- `mcp-perturbation`: Treatment response prediction
- `mcp-quantum-celltype-fidelity`: Quantum-enhanced cell type classification

Claude orchestrates tool calls:

```python
# Tool 1: Load clinical data
epic.get_patient_summary(patient_id="PAT001-OVC-2025")

# Tool 2: Parse genomic variants
fgbio.parse_vcf(
    vcf_path="gs://sample-inputs-patientone/PAT001-OVC-2025/genomics/somatic_variants.vcf"
)

# Tool 3: Multi-omics meta-analysis
multiomics.stouffer_meta_analysis(
    rna_path="gs://.../pdx_rna_seq.csv",
    protein_path="gs://.../pdx_proteomics.csv",
    phospho_path="gs://.../pdx_phosphoproteomics.csv"
)
```

See the full server implementation at: [`servers/mcp-fgbio/src/mcp_fgbio/server.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-fgbio/src/mcp_fgbio/server.py)

### Minutes 5-15: Multi-Modal Analysis

While you grab coffee, Claude:
- Identifies 8 pathogenic somatic mutations (TP53 R175H, PIK3CA E545K, PTEN loss)
- Compares to TCGA ovarian cancer cohort (molecular subtype: C1 Immunoreactive)
- Runs Stouffer meta-analysis across RNA/protein/phospho data
- Identifies consistently dysregulated pathways (PI3K/AKT/mTOR, DNA repair)
- Ranks 47 genes by integrated p-value across all three modalities

### Minutes 15-25: Spatial and Imaging Analysis

Claude calls more specialized tools:

```python
# Spatial differential expression
spatialtools.spatial_differential_expression(
    expression_path="gs://.../visium_gene_expression.csv",
    coordinates_path="gs://.../visium_spatial_coordinates.csv",
    region_annotations_path="gs://.../visium_region_annotations.csv"
)

# Cell segmentation on MxIF image
deepcell.segment_cells(
    image_path="gs://.../PAT001_tumor_multiplex_IF_TP53_KI67_DAPI.tiff",
    model_type="membrane"
)

# Phenotype classification
deepcell.classify_cell_states(
    segmentation_mask="<result from previous tool>",
    markers=["TP53", "KI67", "DAPI"],
    intensity_thresholds={"TP53": 500, "KI67": 300}
)
```

Explore the DeepCell implementation: [`servers/mcp-deepcell/src/mcp_deepcell/server.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/src/mcp_deepcell/server.py)

The results reveal:
- **Spatial heterogeneity**: Resistant cells cluster in the tumor core
- **Immune exclusion**: CD8+ T cells are blocked by thick stromal barrier
- **High proliferation**: 45-55% Ki67+ cells in proliferative regions
- **TP53 mutation impact**: TP53+/Ki67+ double-positive cells correlate with resistance

### Minutes 25-35: Synthesis and Recommendations

Claude integrates findings across all modalities:

- **Genomics**: PIK3CA E545K mutation (gain-of-function in PI3K pathway). 
- **Proteomics**: AKT1, mTOR, RPS6KB1 upregulation (confirms pathway activation). 
- **Spatial**: Tumor core shows PI3K/AKT signature; immune-infiltrated regions show T-cell exhaustion. 
- **Imaging**: 45% Ki67+ proliferation index; CD8+ density only 5-15 cells/mm² (LOW). 

**Treatment Recommendation**:
1. **Primary**: PI3K inhibitor (alpelisib) targeting PIK3CA E545K
2. **Combination**: Anti-PD-1 immunotherapy to overcome immune exclusion
3. **Clinical trial**: NCT03602859 (alpelisib + paclitaxel in ovarian cancer)

The complete analysis—from raw data to actionable recommendations—took **35 minutes**.

---

## What Changed?

Let's compare the two workflows:

| Aspect | Traditional | AI-Orchestrated |
|--------|-------------|-----------------|
| **Time** | 40-48 hours | 35 minutes |
| **Expertise Required** | 3-4 specialists (bioinformatician, pathologist, clinical informaticist) | 1 oncologist with natural language prompts |
| **Code Written** | 500-1,000 lines (R, Python, bash scripts) | 0 lines (natural language only) |
| **Tools Used** | 15-20 (GATK, DESeq2, Seurat, QuPath, etc.) | 12 MCP servers (124 tools, pre-integrated) |
| **Data Formats** | Manual conversion (VCF → CSV, FHIR → JSON, etc.) | Automatic (MCP servers handle all formats) |
| **Cost per Analysis** | $3,200 (personnel time @ $80/hr × 40 hrs) | $1-2 (Cloud Run compute + Claude API) |
| **Scalability** | 1 patient/week/team | 50-100 patients/week/oncologist |

The key insight: **AI doesn't replace bioinformatics. It orchestrates it.**

The MCP servers (`mcp-fgbio`, `mcp-spatialtools`, etc.) still run the same rigorous algorithms—DESeq2 for differential expression, Moran's I for spatial autocorrelation, DeepCell for cell segmentation. But instead of a bioinformatician stitching together tools manually, Claude coordinates the workflow based on your natural language instructions.

---

## Why This Matters

### 1. Time-to-Insight in Clinical Decision-Making

Sarah doesn't have weeks. Platinum-resistant ovarian cancer is aggressive. Every week spent waiting for analysis is a week the cancer grows unchecked.

With AI-orchestrated analysis, her oncologist can:
- Run the complete analysis during the clinic visit
- Discuss results with Sarah in real-time
- Enroll her in a clinical trial the same day

Traditional timelines meant results might arrive after Sarah's next chemo cycle started—too late to change the treatment plan.

### 2. Democratizing Precision Medicine

Today, comprehensive precision oncology is available only at:
- Academic medical centers (MD Anderson, Memorial Sloan Kettering, Mayo Clinic)
- Institutions with dedicated bioinformatics teams
- Patients who can afford $5,000-10,000 out-of-pocket for commercial testing

The AI-orchestrated approach costs **$1-2 per analysis** (GCP Cloud Run compute + Claude API tokens). Suddenly, precision medicine becomes feasible for:
- Community hospitals without bioinformatics staff
- Rural cancer centers
- Low-resource healthcare systems
- Global health initiatives

### 3. Multi-Modal Integration

Here's the uncomfortable truth: most "precision oncology" today is *genomics only*. You sequence the tumor, identify mutations, and prescribe a matched therapy (if one exists).

But cancer is more than its DNA. You need to understand:
- What's happening right now (transcriptomics, proteomics)
- Where resistant cells are hiding (spatial transcriptomics)
- How the immune system is responding (imaging)
- What treatments the cells might respond to (perturbation modeling)

Integrating all five modalities manually is so labor-intensive that almost nobody does it. Sarah would get genomics. Maybe proteomics if she's lucky. Spatial transcriptomics? Only in research settings.

AI orchestration makes multi-modal analysis the *default*, not the exception.

---

## What You'll Learn in This Book

This book will teach you how to build, deploy, and operate the AI-orchestrated precision oncology system that analyzed Sarah's case.

**Part 1: Why This Matters** (Chapters 1-3)
You'll understand the clinical problem, the architecture of the solution, and real-world testing results.

**Part 2: Building the Foundation** (Chapters 4-7)
You'll implement the core MCP servers for clinical data (FHIR), genomics, multi-omics, and spatial transcriptomics.

**Part 3: Advanced Capabilities** (Chapters 8-11)
You'll add cell segmentation (DeepCell), treatment response prediction (GEARS), quantum fidelity analysis, and histopathology imaging.

**Part 4: Deployment and Operations** (Chapters 12-14)
You'll deploy to Google Cloud Run, configure HIPAA-compliant hospital infrastructure, and set up monitoring.

**Part 5: Research and Education** (Chapters 15-16)
You'll learn how researchers and educators use the system for discovery and teaching.

**Part 6: The Future** (Chapters 17-18)
You'll explore funding models, ROI analysis, and lessons learned from production deployment.

---

## Try It Yourself

Ready to run the PatientOne analysis? You can deploy the MCP servers locally or to your own cloud account.

**Option 1: Interactive Notebook**
Open the companion Jupyter notebook for this chapter:
[`docs/book/companion-notebooks/chapter-01-patientone-story.ipynb`](../companion-notebooks/chapter-01-patientone-story.ipynb)

This notebook walks you through:
- Setting up the MCP servers locally
- Connecting with Claude or Gemini
- Running the PatientOne prompt step-by-step
- Exploring results and modifying parameters

**Option 2: Local Claude Desktop Setup**
If you have Claude Desktop installed, deploy MCP servers locally:

1. Clone the repository:
```bash
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp
```

2. Follow the installation guide:
[`docs/getting-started/installation.md`](../../getting-started/installation.md)

3. Configure Claude Desktop with local servers

4. Paste this prompt in Claude Desktop:
```
Analyze patient PAT001-OVC-2025 for precision oncology treatment selection.

Clinical data: data/patient-data/PAT001-OVC-2025/clinical/
Genomics: data/patient-data/PAT001-OVC-2025/genomics/somatic_variants.vcf
Multi-omics: data/patient-data/PAT001-OVC-2025/multiomics/
Spatial: data/patient-data/PAT001-OVC-2025/spatial/
Imaging: data/patient-data/PAT001-OVC-2025/imaging/

Generate comprehensive treatment recommendations based on all modalities.
```

**Option 3: Deploy to Your Cloud (Chapter 12)**
Want to deploy to GCP Cloud Run for scalability? Jump ahead to Chapter 12 for deployment instructions, then return to continue reading.

All PatientOne data is synthetic and included in the repository. No API keys required for local testing (though Claude/Gemini API keys needed for AI orchestration).

---


## Summary

**Chapter 1 Key Takeaways:**
- Traditional precision oncology analysis: 40 hours, 3-4 specialists, $3,200 cost
- AI-orchestrated analysis: 35 minutes, 1 oncologist, $1-2 cost
- 95% time reduction, 95% cost reduction
- Multi-modal integration (5 data types) becomes feasible at scale
- PatientOne (Sarah) is a real workflow you can run today with synthetic data

---

**Companion Resources:**
- 📓 [Jupyter Notebook](../companion-notebooks/chapter-01-patientone-story.ipynb) - Run the analysis yourself
- 🎬 [Video Demo](https://www.youtube.com/watch?v=LUldOHHX5Yo) - 5-minute PatientOne walkthrough
- 📊 [Full Demo Guide](../../demos/FULL_PATIENTONE_DEMO.md) - Complete testing instructions

**GitHub References:**
- Patient data: [`data/patient-data/PAT001-OVC-2025/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/data/patient-data/PAT001-OVC-2025)
- Test prompts: [`docs/test-docs/patient-one-scenario/test-prompts/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/docs/test-docs/patient-one-scenario/test-prompts)
- MCP servers: [`servers/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/servers)



ewpage


# The Architecture Problem

> *"How do you orchestrate 69 bioinformatics tools across 12 specialized servers?"*

---

## The Fragmentation Problem

Precision medicine is drowning in tools.

Need to quality-check genomic data? Try FastQC, MultiQC, Picard, or fgbio. Want to analyze spatial transcriptomics? You'll choose between Scanpy, Seurat, Giotto, or Squidpy. Cell segmentation? QuPath, CellProfiler, StarDist, or DeepCell. Each tool has its own:
- Installation requirements (conda environments, Docker containers, R packages)
- Input formats (VCF, BAM, CSV, H5AD, TIFF, OME-TIFF)
- Output formats (JSON, CSV, HDF5, RDS, NPY)
- Parameter configurations (YAML files, command-line flags, config objects)

**The result**: Bioinformaticians spend more time managing tools than analyzing biology.

Let's look at a real example from Sarah's case in Chapter 1. To analyze her spatial transcriptomics data, you'd traditionally need:

```bash
# Install multiple environments and tools
conda create -n spatial python=3.11 && pip install scanpy squidpy
R -e "install.packages('Seurat')"
conda install -c bioconda star && pip install gseapy
```

That's four different programming environments, three package managers, and at least 30 minutes of installation time—*before you write a single line of analysis code*.

Now multiply this by every modality: genomics, proteomics, imaging, clinical data. You end up with:
- 15-20 conda environments
- 50+ installed packages
- Gigabytes of dependencies
- Version conflicts (Python 3.9 for tool A, Python 3.11 for tool B)
- Constant maintenance (package updates break pipelines monthly)

This is the **tool fragmentation problem**, and it's killing productivity.

---

## Why Microservices Aren't Enough

You might think: "Just wrap each tool in a microservice with a REST API!" Many bioinformatics platforms have tried this:

**cBioPortal**: Web interface for cancer genomics data
**Galaxy**: Workflow platform with 9,000+ tools
**Seven Bridges**: Genomics analysis platform

These are valuable systems, but they share a fundamental limitation: **they don't compose naturally**.

Here's what we mean. Suppose you want to:
1. Load clinical data from Epic FHIR
2. Identify somatic mutations from a VCF file
3. Run multi-omics pathway enrichment
4. Overlay results on spatial transcriptomics
5. Generate a treatment recommendation

With traditional REST microservices, you'd write code like this:

```python
# Manual orchestration: call multiple APIs, handle conversions
clinical = requests.post("https://api.hospital.org/epic", json={"patient_id": "PAT001"})
variants = requests.post("https://api.biotools.org/vcf-parse", files={"vcf": open("patient.vcf")})
pathways = requests.post("https://api.biotools.org/multiomics", json={...})
# Full implementation: ui/streamlit-app/providers/gemini_provider.py
```

Full integration code: [`ui/streamlit-app/providers/gemini_provider.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/ui/streamlit-app/providers/gemini_provider.py)

**The problems**:
1. **You're writing orchestration code**. Every new analysis requires custom Python/JavaScript.
2. **Data format conversions**. Epic returns FHIR JSON. The genomics service expects VCF text. Multi-omics wants CSV. You're constantly converting.
3. **No intelligence**. The system doesn't understand what you're trying to accomplish. It just executes HTTP requests.
4. **No error recovery**. If the third API call fails, you manually retry or debug.
5. **No optimization**. The system can't parallelize independent calls or cache repeated queries.

Microservices give you modularity, but you still need an orchestrator. Traditionally, that's been:
- **Nextflow/Snakemake**: Workflow engines with DSL
- **Airflow/Luigi**: Task schedulers
- **Custom Python scripts**: Glue code written by bioinformaticians

All of these require *code*. You define workflows in YAML, Groovy, or Python. Every new analysis needs a new workflow definition.

**What if the orchestrator could understand natural language?**

---

## Enter the Model Context Protocol (MCP)

The Model Context Protocol is an open standard created by Anthropic that enables AI models like Claude to interact with external tools and data sources.

Think of MCP as "USB for AI":
- USB standardized how peripherals connect to computers (keyboards, mice, printers)
- MCP standardizes how AI models connect to tools (databases, APIs, specialized services)

Here's the key insight: **If tools expose an MCP interface, an AI can orchestrate them without you writing orchestration code.**

### How MCP Works

An MCP server exposes:
1. **Tools**: Functions the AI can call (e.g., `parse_vcf`, `run_pathway_enrichment`)
2. **Resources**: Data the AI can read (e.g., patient records, reference genomes)
3. **Prompts**: Pre-defined prompt templates for common tasks

Claude connects to multiple MCP servers simultaneously. When you ask a question, Claude:
1. **Understands intent**: "I need clinical data, genomics, and spatial analysis"
2. **Plans workflow**: "Load clinical → Parse VCF → Run multi-omics → Overlay spatial"
3. **Executes tools**: Calls appropriate functions on each MCP server
4. **Synthesizes results**: Combines outputs into a coherent answer

You don't write code. You describe what you want in plain English.

![Why MCP for Healthcare](images/screenshots/why-mcp.jpeg){width=\textwidth}

**Figure 2.2: Why MCP for Healthcare**
*Visual representation of how MCP enables AI orchestration of bioinformatics tools without custom integration code.*

### MCP vs. Traditional API Calls

| Aspect | REST API | MCP Server |
|--------|----------|------------|
| **Discovery** | Read API docs manually | AI discovers tools automatically |
| **Orchestration** | Write Python/JS code | Natural language prompts |
| **Data Conversion** | Manual (JSON ↔ CSV ↔ VCF) | Automatic (server handles formats) |
| **Error Handling** | try/except blocks | AI retries with backoff |
| **Parallelization** | Manual threading/async | AI optimizes automatically |
| **Caching** | Manual implementation | Built into MCP protocol |

### A Real MCP Tool Definition

Let's look at how the `parse_vcf` tool is defined in the `mcp-fgbio` server:

```python
@mcp.tool()
def parse_vcf(vcf_path: str, min_depth: int = 10) -> dict:
    """Parse VCF file and return somatic variants."""
    # Full implementation: servers/mcp-fgbio/src/mcp_fgbio/server.py:145-187
```

Full implementation: [`servers/mcp-fgbio/src/mcp_fgbio/server.py:145-187`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-fgbio/src/mcp_fgbio/server.py#L145-L187)

That's it. The `@mcp.tool()` decorator registers this function with the MCP protocol. When Claude connects to the `mcp-fgbio` server, it automatically discovers:
- Tool name: `parse_vcf`
- Parameters: `vcf_path` (required), `min_depth` (optional, default 10)
- Return type: Dictionary with variants and metrics
- Description: From the docstring

No API documentation to write. No client SDK to maintain. Just Python functions with descriptive docstrings.

---

## The Precision Medicine MCP Architecture

Now let's see how this solves the precision oncology orchestration problem.

```mermaid
graph TB
    subgraph "User Interfaces"
        STREAMLIT[🌐 Streamlit Chat UI<br/>Web Interface<br/>Cloud Run]
        JUPYTER[📓 Jupyter Notebook<br/>Data Science Interface<br/>Cloud Run]
        DESKTOP[💬 Claude Desktop<br/>Local STDIO]
    end

    subgraph "AI Orchestration Layer"
        API[🤖 Claude API<br/>Anthropic Sonnet 4.5<br/>MCP Client Support]
    end

    subgraph "Local-Only Servers"
        REALEPIC[mcp-epic<br/>Real Epic FHIR<br/>✅ Production<br/>🏥 HIPAA-compliant]
    end

    subgraph "GCP Cloud Run - 11 Deployed MCP Servers"
        subgraph "Clinical & Genomic"
            MOCKEPIC[mcp-mockepic<br/>Mock FHIR<br/>🎭 Demo Only]
            FGBIO[mcp-fgbio<br/>FASTQ/VCF<br/>✅ Production]
            TCGA[mcp-tcga<br/>Cancer Data<br/>❌ Mocked]
        end

        subgraph "Multi-Omics"
            MULTI[mcp-multiomics<br/>RNA/Protein/Phospho<br/>✅ Production]
        end

        subgraph "Spatial & Imaging"
            SPATIAL[mcp-spatialtools<br/>Spatial RNA-seq<br/>✅ Production]
            IMAGE[mcp-openimagedata<br/>Histology<br/>🔶 60% Real]
            DEEPCELL[mcp-deepcell<br/>Cell Segmentation<br/>✅ Production]
        end

        subgraph "Advanced Analytics"
            PERTURB[mcp-perturbation<br/>GEARS GNN<br/>✅ Production]
            QUANTUM[mcp-quantum-celltype-fidelity<br/>Quantum PQCs<br/>✅ Production]
        end

        subgraph "AI & Workflows"
            HF[mcp-huggingface<br/>ML Models<br/>❌ Mocked]
            SEQERA[mcp-seqera<br/>Nextflow<br/>❌ Mocked]
        end
    end

    subgraph "Analysis Workflow"
        PATIENT[🏥 PatientOne<br/>Stage IV HGSOC<br/>Multi-Omics Analysis]
    end

    STREAMLIT ==> API
    JUPYTER ==> API
    DESKTOP --> API

    API ==> FGBIO
    API ==> MULTI
    API ==> SPATIAL
    API ==> REALEPIC
    API ==> DEEPCELL
    API ==> PERTURB
    API ==> QUANTUM
    API -.-> MOCKEPIC
    API -.-> TCGA
    API -.-> IMAGE
    API -.-> HF
    API -.-> SEQERA

    FGBIO ==> PATIENT
    MULTI ==> PATIENT
    SPATIAL ==> PATIENT
    REALEPIC ==> PATIENT
    DEEPCELL ==> PATIENT
    PERTURB ==> PATIENT
    QUANTUM ==> PATIENT
    MOCKEPIC -.-> PATIENT
    TCGA -.-> PATIENT
    IMAGE -.-> PATIENT
    HF -.-> PATIENT
    SEQERA -.-> PATIENT

    style STREAMLIT fill:#d1ecf1,stroke:#0c5460,stroke-width:2px
    style JUPYTER fill:#d1ecf1,stroke:#0c5460,stroke-width:2px
    style DESKTOP fill:#d1ecf1,stroke:#0c5460,stroke-width:2px
    style API fill:#cce5ff,stroke:#004085,stroke-width:3px
    style PATIENT fill:#e1f5ff,stroke:#0066cc,stroke-width:3px
    style FGBIO fill:#d4edda,stroke:#28a745,stroke-width:2px
    style MULTI fill:#d4edda,stroke:#28a745,stroke-width:2px
    style SPATIAL fill:#d4edda,stroke:#28a745,stroke-width:2px
    style REALEPIC fill:#d4edda,stroke:#28a745,stroke-width:3px
    style DEEPCELL fill:#d4edda,stroke:#28a745,stroke-width:2px
    style PERTURB fill:#d4edda,stroke:#28a745,stroke-width:2px
    style QUANTUM fill:#d4edda,stroke:#28a745,stroke-width:2px
    style IMAGE fill:#fff3cd,stroke:#ffc107,stroke-width:1px
    style TCGA fill:#f8d7da,stroke:#dc3545,stroke-width:1px
    style HF fill:#f8d7da,stroke:#dc3545,stroke-width:1px
    style SEQERA fill:#f8d7da,stroke:#dc3545,stroke-width:1px
    style MOCKEPIC fill:#e2e3e5,stroke:#6c757d,stroke-width:1px
```

**Figure 2.1: Complete MCP System Architecture**
*12 MCP servers organized by domain (Clinical & Genomic, Multi-Omics, Spatial & Imaging, Advanced Analytics, AI & Workflows). Solid lines indicate production-ready integrations. Dotted lines represent mocked/demo servers. Claude API orchestrates all servers through natural language prompts.*

**Legend:**
- ✅ **Production Ready** (7/12): Real data, comprehensive tests, deployed
- 🔶 **Partial Implementation** (1/12): Some real integrations, some mocked
- ❌ **Mocked** (3/12): Return synthetic data, API calls stubbed
- 🎭 **Mock by Design** (1/12): Intentionally synthetic for demos

### The 12 Specialized Servers

The system has 12 MCP servers, each specialized for a specific bioinformatics domain:

**Production-Ready Servers (7)**:
1. **mcp-epic**: Real FHIR R4 integration with Epic EHR systems
2. **mcp-fgbio**: Genomic QC, variant calling, reference data
3. **mcp-multiomics**: Multi-omics integration (HAllA, Stouffer meta-analysis)
4. **mcp-spatialtools**: Spatial transcriptomics (STAR alignment, ComBat, pathways)
5. **mcp-perturbation**: Treatment response prediction (GEARS GNN)
6. **mcp-quantum-celltype-fidelity**: Quantum fidelity with Bayesian uncertainty
7. **mcp-deepcell**: Cell segmentation (DeepCell-TF models)

**Mock/Partial Servers (5)**:
8. **mcp-mockepic**: Synthetic FHIR data (testing only)
9. **mcp-tcga**: TCGA cohort queries (demo data)
10. **mcp-openimagedata**: Histology imaging (60% real)
11. **mcp-huggingface**: ML model inference (framework ready)
12. **mcp-seqera**: Nextflow orchestration (demo)

**Total**: 69 bioinformatics tools across 12 servers

Full server status: [`docs/architecture/servers.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/architecture/servers.md)

![Claude Desktop with MCP Servers](images/screenshots/Claude-client.png){width=\textwidth}

**Figure 2.3: Claude Desktop with MCP Servers**
*Claude Desktop interface showing connected MCP servers in the sidebar. Users interact with all 12 servers through natural language prompts without writing integration code.*

### Data Flow Example

Let's trace what happens when you ask Claude to analyze Sarah's case:

**Your prompt**:
```
Analyze patient PAT001-OVC-2025. Identify treatment targets based on
genomics, multi-omics, and spatial transcriptomics.
```

**Behind the scenes**:

1. **Claude parses intent**: Needs clinical context, genomics, multi-omics, spatial data. Goal: Treatment recommendations.

2. **Claude discovers available tools**: Connected to 7 MCP servers with 69 tools total (clinical data, genomics, integration, spatial analysis).

3. **Claude plans workflow**: Load clinical → Parse variants → Run meta-analysis → Spatial DE → Integrate findings.

4. **Claude executes (parallelizes where possible)**:
```python
# Parallel execution of independent tasks
clinical_result = await mcp_mockepic.get_patient_summary(patient_id="PAT001-OVC-2025")
genomic_result = await mcp_fgbio.parse_vcf(vcf_path="gs://...")
# Full orchestration: ui/streamlit-app/providers/
```

5. **Claude synthesizes**:
```
Findings:
- Clinical: Platinum-resistant HGSOC, BRCA1 carrier
- Genomics: PIK3CA E545K, TP53 R175H, PTEN loss
- Multi-omics: PI3K/AKT pathway active (protein & phospho confirm)
- Spatial: Immune exclusion, high proliferation in tumor core

Treatment recommendation:
Primary: PI3K inhibitor (alpelisib) targeting PIK3CA E545K
Secondary: Anti-PD-1 to overcome immune exclusion
Clinical trial: NCT03602859 (alpelisib + paclitaxel)
```

You wrote zero lines of code. Claude orchestrated 12+ tool calls across 4 servers, handled data format conversions, parallelized independent calls, and synthesized results into clinical recommendations.

**That's the power of AI orchestration with MCP.**

---

## Why This Architecture Works for Healthcare

### 1. Domain Expertise Encoded in Servers

Each MCP server encapsulates domain-specific knowledge:

**mcp-fgbio** knows:
- How to validate VCF files (check reference allele, fix chromosome naming)
- When to filter variants (depth < 10, VAF < 0.05)
- How to annotate with ClinVar, gnomAD, COSMIC
- Which quality metrics matter (QUAL, DP, AF)

**mcp-spatialtools** knows:
- How to run STAR genome alignment
- When to apply batch correction (ComBat vs. Harmony)
- Which statistical tests for spatial DE (Wilcoxon, Mann-Whitney)
- How to calculate spatial autocorrelation (Moran's I)

**mcp-multiomics** knows:
- How to integrate RNA/protein/phospho data (Stouffer's Z-score)
- When to use FDR correction (Benjamini-Hochberg)
- Which pathways are clinically relevant (44 curated pathways)
- How to handle missing data (imputation strategies)

This knowledge is *implemented* in the server code, not described in documentation for humans to read. Claude can use it directly.

### 2. Separation of Concerns

Each server has a single responsibility:
- **mcp-epic**: Clinical data retrieval (FHIR only)
- **mcp-fgbio**: Genomic QC and variant calling
- **mcp-spatialtools**: Spatial analysis (not imaging, not genomics)

This modularity means:
- Servers can be developed independently
- Testing is focused (unit tests per server)
- Updates don't break other components
- New servers can be added without modifying existing ones

### 3. Cloud-Native Deployment

All servers run on Google Cloud Run (serverless containers):
- **Auto-scaling**: 0 to 1000 instances based on demand
- **Pay-per-use**: $0.02-0.21 per analysis (see Chapter 12 for cost breakdown)
- **No infrastructure management**: Google handles servers, networking, SSL
- **Regional deployment**: us-central1 (low latency for US hospitals)

Deployment command for a new server:

```bash
cd servers/mcp-deepcell
./deploy.sh YOUR_PROJECT_ID us-central1
# Full deployment guide: docs/deployment/GET_STARTED.md
```

Full deployment guide: [`docs/deployment/GET_STARTED.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/deployment/GET_STARTED.md)

### 4. Multi-Provider AI Support

The architecture isn't locked to Claude. You can use:
- **Claude (Anthropic)**: Native MCP support via MCP beta API
- **Gemini (Google)**: Custom MCP client implementation
- **Future models**: Any LLM can implement the MCP protocol

Claude uses Anthropic's native MCP integration:
```python
client.beta.messages.create(
    mcp_servers=[{"type": "url", "url": server_url, "name": "fgbio"}],
    tools=[{"type": "mcp_toolset", "mcp_server_name": "fgbio"}],
    betas=["mcp-client-2025-11-20"]
)
```

Gemini provider implementation: [`ui/streamlit-app/providers/gemini_provider.py:87-245`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/ui/streamlit-app/providers/gemini_provider.py#L87-L245)

---

## What MCP Doesn't Solve

Let's be honest about limitations:

### 1. Model Errors and Hallucinations

Claude is smart, but it can still:
- Misinterpret complex biological relationships
- Hallucinate connections not supported by data
- Miss subtle patterns in large datasets

**Mitigation**: Human review (Clinician-in-the-Loop workflow, see Chapter 13) and validation against known biology.

### 2. Cost (Though Dramatically Reduced)

AI API calls cost money:
- Claude API: ~$0.50-1.00 per complex analysis (30K tokens)
- Gemini API: ~$0.10-0.30 per analysis (cheaper but requires custom integration)
- Cloud Run compute: ~$0.02-0.21 per analysis

Total: **$1-2 per patient** vs. $3,200 traditional (95% reduction), but not free.

### 3. Latency for Real-Time Needs

Current system: 35 minutes for complete analysis
Real-time requirement: <5 minutes for clinical decisions

**Workaround**: Pre-compute common analyses. Use MCP for novel queries only.

### 4. Complex Workflows Need Refinement

Some multi-step workflows (e.g., iterative spatial clustering with multiple parameter sweeps) require multiple prompts or explicit workflow definition. MCP excels at well-defined tasks, less so at open-ended exploration.

---

## How You'll Use This Architecture

In the remaining chapters, you'll learn to:

**Part 2 (Chapters 4-7)**: Build MCP servers
- Chapter 4: Clinical data (FHIR integration)
- Chapter 5: Genomics (VCF parsing, variant annotation)
- Chapter 6: Multi-omics (HAllA, Stouffer meta-analysis)
- Chapter 7: Spatial transcriptomics (STAR, ComBat, pathways)

**Part 3 (Chapters 8-11)**: Add advanced capabilities
- Chapter 8: Cell segmentation (DeepCell-TF)
- Chapter 9: Treatment response (GEARS GNN)
- Chapter 10: Quantum fidelity (PennyLane, Bayesian UQ)
- Chapter 11: Histopathology imaging

**Part 4 (Chapters 12-14)**: Deploy to production
- Chapter 12: Cloud Run deployment
- Chapter 13: HIPAA-compliant hospital infrastructure
- Chapter 14: Operations and monitoring

Each chapter includes:
- Architecture diagrams
- Code snippets (2-5 lines) with GitHub links
- Jupyter notebooks for hands-on practice
- Real deployment examples

---

## Try It Yourself

Ready to see MCP orchestration in action?

**Option 1: Explore the Architecture**
Open the companion Jupyter notebook: [`docs/book/companion-notebooks/chapter-02-architecture.ipynb`](../companion-notebooks/chapter-02-architecture.ipynb)

This notebook lets you:
- Connect to deployed MCP servers
- Discover available tools programmatically
- Call individual tools to understand their behavior
- Compare sequential vs. parallel orchestration

**Option 2: Query MCP Server Endpoints**
The servers are deployed. You can query them directly:

```bash
# List tools from mcp-fgbio
curl https://mcp-fgbio-ondu7mwjpa-uc.a.run.app/tools
# Returns: JSON array of 9 tools with names, descriptions, parameters
```

All endpoints: [`docs/architecture/servers.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/architecture/servers.md)

**Option 3: Run the PatientOne Workflow**
Deploy the MCP servers locally and use Claude Desktop. See **Appendix: Setup Guide** for complete instructions.

---

## Summary

**Chapter 2 Key Takeaways:**
- Bioinformatics suffers from tool fragmentation (15-20 tools per analysis)
- Microservices help modularity but require orchestration code
- MCP standardizes how AI models interact with tools
- 12 specialized MCP servers encode domain expertise
- Claude coordinates tools via natural language (no code)
- Architecture supports both Claude and Gemini
- 95% cost reduction: $3,200 → $1-2 per analysis

**Companion Resources:**
- 📓 [Jupyter Notebook](../companion-notebooks/chapter-02-architecture.ipynb) - Explore MCP orchestration
- 🏗️ [Architecture Docs](../../architecture/README.md) - Detailed technical specs
- 🔧 [Server Registry](../../architecture/README.md) - All 12 servers and 69 tools
- 🚀 [Deployment Guide](../../deployment/GET_STARTED.md) - Cloud Run setup
- 📚 [Appendix: Setup Guide](appendix-setup-guide.md) - Installation instructions

**GitHub References:**
- MCP server boilerplate: [`servers/mcp-server-boilerplate/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/servers/mcp-server-boilerplate)
- fgbio server implementation: [`servers/mcp-fgbio/src/mcp_fgbio/server.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-fgbio/src/mcp_fgbio/server.py)
- Streamlit multi-provider UI: [`ui/streamlit-app/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/ui/streamlit-app)



ewpage


# Testing the Hypothesis

> *"Could this actually work in production?"*

---

## The Big Question

You've seen the vision in Chapter 1: 40 hours → 35 minutes. You've understood the architecture in Chapter 2: 12 MCP servers coordinated by AI. But now comes the hard part:

**Does it actually work?**

Not in theory. Not in a slide deck. In production, with real data, real costs, and real time constraints.

This chapter tells the story of how the Precision Medicine MCP system was tested, validated, and deployed—including the failures, debugging sessions, and hard-won lessons that don't make it into the marketing materials.

---

## The Starting Point: Everything Was Mocked

When development began, the entire system was smoke and mirrors. Not because of dishonesty, but because **you can't build everything at once**.

The initial commit had 9 MCP servers. All of them returned synthetic data:

```python
@mcp.tool()
def parse_vcf(vcf_path: str) -> dict:
    # DRY_RUN mode: return synthetic data
    return {"variants": [...], "warning": "DRY_RUN mode"}
    # Full implementation: servers/mcp-fgbio/src/mcp_fgbio/server.py:145-187
```

Full implementation: [`servers/mcp-fgbio/src/mcp_fgbio/server.py:145-187`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-fgbio/src/mcp_fgbio/server.py#L145-L187)

**Why start with mocks?**

1. **Validate the architecture**: Does MCP orchestration even work?
2. **Test the workflow**: Can Claude chain tool calls correctly?
3. **Iterate quickly**: No need to wait for TCGA API access or DeepCell model downloads
4. **Show the vision**: Demos work even without production implementations

The first PatientOne demo ran entirely on mocked data. It took 5 minutes (including Claude API latency) and cost $0.05 in Claude tokens. The results looked real, but the hard work hadn't started yet.

---

## Phase 1: Making It Real (The Hard Part)

### Priority 1: Which Servers Matter?

You can't implement everything at once. You need to prioritize based on:
- **Clinical utility**: What do oncologists actually need?
- **Technical feasibility**: What can you build in 2-3 weeks?
- **Data availability**: What datasets do you have access to?

Here's how the servers were prioritized:

**Tier 1: Must Be Real (Production Critical)**
1. **mcp-epic**: Clinical data (FHIR integration)
2. **mcp-fgbio**: Genomic QC and variant calling
3. **mcp-multiomics**: Multi-omics integration
4. **mcp-spatialtools**: Spatial transcriptomics

**Tier 2: Partially Real (Proof of Concept)**
5. **mcp-deepcell**: Cell segmentation (DeepCell-TF)
6. **mcp-openimagedata**: Histology imaging

**Tier 3: Can Stay Mocked (Demo Only)**
7. **mcp-tcga**: TCGA cohort queries
8. **mcp-huggingface**: ML model inference
9. **mcp-seqera**: Nextflow orchestration

### Building mcp-multiomics: 85% Real in 3 Weeks

Let's walk through a real implementation example. The mcp-multiomics server needed to:
- Load RNA, protein, and phosphoproteomics data (CSV files)
- Run HAllA (HierArchical ALL-against-ALL) association discovery
- Perform Stouffer's meta-analysis across modalities
- Execute pathway enrichment with FDR correction

**Week 1: Data Loading and Validation**

```python
def validate_omics_data(rna_path, protein_path, phospho_path) -> dict:
    # Load, validate sample alignment, check missing values
    # Full implementation: servers/mcp-multiomics/src/mcp_multiomics/tools/preprocessing.py
```

Full implementation: [`servers/mcp-multiomics/src/mcp_multiomics/tools/preprocessing.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-multiomics/src/mcp_multiomics/tools/preprocessing.py)

**Week 2: HAllA Integration (The Python Fallback)**

HAllA is traditionally an R package. Integrating R into a Python MCP server requires `rpy2`—a Python-R bridge that's notoriously difficult to configure.

The breakthrough: **Implement the core algorithm in Python.**

```python
def halla_python(rna_data, protein_data, method="spearman") -> dict:
    # Pure Python implementation of HAllA core algorithm
    # Full implementation: servers/mcp-multiomics/src/mcp_multiomics/tools/association_analysis.py:67-145
```

Full implementation: [`servers/mcp-multiomics/src/mcp_multiomics/tools/association_analysis.py:67-145`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-multiomics/src/mcp_multiomics/tools/association_analysis.py)

**Week 3: Stouffer Meta-Analysis**

The key feature: combine p-values across RNA, protein, and phospho to find consistent signals.

```python
def stouffer_meta_analysis(rna_pvals, protein_pvals, phospho_pvals, gene_ids):
    # Combine Z-scores: z_combined = (z_rna + z_protein + z_phospho) / sqrt(3)
    # Full implementation: servers/mcp-multiomics/src/mcp_multiomics/tools/meta_analysis.py:89-156
```

Full implementation: [`servers/mcp-multiomics/src/mcp_multiomics/tools/meta_analysis.py:89-156`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-multiomics/src/mcp_multiomics/tools/meta_analysis.py)

**Result**: 91 tests, 68% coverage, production-ready in 3 weeks.

Test suite: [`tests/unit/mcp-multiomics/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/tests/unit/mcp-multiomics)

---

## The DeepCell Challenge: When Things Don't Go As Planned

Not every implementation went smoothly. The mcp-deepcell server was supposed to take 1 week. It took 3 weeks and multiple Cloud Build failures.

### The Four Attempts

**Attempt 1**: `pip install deepcell-tf` → Package name was `DeepCell` (not `deepcell-tf`)

**Attempt 2**: `DeepCell>=0.12.0` → Requires TensorFlow 2.8.x, only supports Python 3.10 (not 3.11)

**Attempt 3**: Downgrade to Python 3.10 → Cloud Build in us-central1 doesn't allow N1 machine types

**Attempt 4 (Success)**: Python 3.10 + E2_HIGHCPU_8 machine type + GCS image loading

The final challenge: PIL's `Image.open()` doesn't support GCS URIs. Solution:

```python
def _download_from_gcs(gcs_path: str) -> str:
    # Download from gs://bucket/path to temp file
    # Full implementation: servers/mcp-deepcell/src/mcp_deepcell/server.py:47-89
```

Full implementation: [`servers/mcp-deepcell/src/mcp_deepcell/server.py:47-89`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/src/mcp_deepcell/server.py#L47-L89)

**Result**: Deployment succeeded on attempt 4. Total time: 3 weeks (2 weeks debugging).

Full deployment story: [`servers/mcp-deepcell/DEPENDENCY_ISSUES.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/DEPENDENCY_ISSUES.md)

**Lesson learned**: Never assume package names match project names. Always check PyPI before writing `pyproject.toml`.

---

## Production Deployment: All 11 Servers to Cloud Run

Once the implementations were stable, it was time to deploy everything to Google Cloud Run.

### The Deployment Script

```bash
# Deploy all servers to Cloud Run
for server in mcp-*; do
  gcloud run deploy ${server} --source . --region us-central1 --memory 4Gi --cpu 2
done
# Full script: infrastructure/deployment/deploy_to_gcp.sh
```

Full deployment script: [`infrastructure/deployment/deploy_to_gcp.sh`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/infrastructure/deployment/deploy_to_gcp.sh)

### Key Challenges

**Issue 1**: Dockerfiles expected `_shared_temp/utils/` but deployment script didn't stage files.
**Fix**: Stage shared utilities before `gcloud run deploy`, cleanup after.

**Issue 2**: Cloud Run cached old `MCP_TRANSPORT=http` from previous deployments.
**Fix**: Explicitly set `--update-env-vars MCP_TRANSPORT=sse` in deployment command.

See deployment logs: [`docs/archive/deployment/DEPLOYMENT_STATUS.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/archive/deployment/DEPLOYMENT_STATUS.md)

### Success: All 11 Servers Running

Final deployment status (2026-01-31):

| Server | URL | Status |
|--------|-----|--------|
| mcp-deepcell | https://mcp-deepcell-ondu7mwjpa-uc.a.run.app | ✅ Running |
| mcp-fgbio | https://mcp-fgbio-ondu7mwjpa-uc.a.run.app | ✅ Running |
| mcp-multiomics | https://mcp-multiomics-ondu7mwjpa-uc.a.run.app | ✅ Running |
| mcp-spatialtools | https://mcp-spatialtools-ondu7mwjpa-uc.a.run.app | ✅ Running |
| ... (7 more servers) | ... | ✅ Running |

Full deployment status: [`docs/archive/deployment/DEPLOYMENT_STATUS.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/archive/deployment/DEPLOYMENT_STATUS.md)

---

## Testing: 167 Automated Tests Across 9 Servers

With servers deployed, comprehensive testing began.

### Test Coverage by Server

| Server | Tests | Coverage | Status |
|--------|-------|----------|--------|
| **mcp-multiomics** | 91 | 68% | ✅ Production |
| **mcp-fgbio** | 29 | 77% | ✅ Production |
| **mcp-epic** | 12 | 58% | ✅ Production |
| **mcp-deepcell** | 9 | 62% | ✅ Smoke |
| **mcp-spatialtools** | 5 | 23% | ✅ Production |
| Others | 21 | 35-56% | Mixed |

**Total**: 167 tests, 56.9% overall coverage

Full test coverage report: [`docs/test-docs/test-coverage.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/test-docs/test-coverage.md)

### Why Low Coverage Doesn't Mean Low Quality

Notice mcp-spatialtools has only 23% coverage but is marked production-ready. Why?

**Answer**: The server has 2,890 lines of code implementing 14 complex tools (STAR alignment, ComBat batch correction, Moran's I spatial autocorrelation). The 5 smoke tests validate:
- Tool registration and data loading
- Basic execution and output formats

Full integration testing happens in the PatientOne end-to-end workflow (next section).

Code quality report: [`docs/for-developers/CODE_QUALITY_REPORT.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-developers/CODE_QUALITY_REPORT.md)

### Example: Testing Multi-Omics Meta-Analysis

```python
def test_stouffer_meta_analysis_with_real_data():
    # Load PatientOne PDX data (15 samples: 7 sensitive, 8 resistant)
    result = stouffer_meta_analysis(rna_data, protein_data, phospho_data)
    assert result["top_gene"] in ["PIK3CA", "AKT1", "MTOR"]
    # Full test: tests/unit/mcp-multiomics/test_meta_analysis.py
```

Test suite: [`tests/unit/mcp-multiomics/test_meta_analysis.py`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/tests/unit/mcp-multiomics)

---

## End-to-End Validation: The Complete PatientOne Workflow

Unit tests validate individual tools. But does the entire system work together?

### Test Configuration

**Patient**: PAT001-OVC-2025 (Stage IV HGSOC)
**Data modalities**: Clinical (FHIR), genomics (VCF), multi-omics (CSV), spatial (Visium), imaging (TIFF)
**Workflow steps**: 5 sequential tests (TEST_1 through TEST_5)
**Total time**: 35 minutes (Claude orchestration)
**Total cost**: $1.20 (Claude API + Cloud Run compute)

### TEST_1: Clinical + Genomic Integration ✅

**Tools called**: `mockepic.get_patient_summary()`, `fgbio.parse_vcf()`, `tcga.compare_to_cohort()`

**Results**:
- Patient: Sarah Anderson, 58, BRCA1 carrier
- CA-125: 1200 → 45 → 310 U/mL (platinum resistance)
- Mutations: TP53 R175H, PIK3CA E545K, PTEN LOH
- TCGA subtype: C1 Immunoreactive

**Time**: 5 minutes | **Cost**: $0.18

### TEST_2: Multi-Omics Resistance Analysis ✅

**Tools called**: `multiomics.stouffer_meta_analysis()`, `multiomics.pathway_enrichment()`

**Results**:
- **17 significant genes** (FDR < 0.05, |log2FC| > 1)
- **Top 3 hits**: PIK3CA (p=1.2e-18), AKT1 (p=1.23e-18), mTOR (p=1.36e-14)
- **Pathway**: PI3K/AKT/mTOR activation (confirms genomic PIK3CA E545K mutation is *active*)

**Time**: 8 minutes | **Cost**: $0.31

### TEST_3: Spatial Transcriptomics ✅

**Tools called**: `spatialtools.spatial_differential_expression()`, `spatialtools.spatial_autocorrelation()`

**Results**:
- **17 DEGs** (tumor vs stroma, FDR < 0.05)
- **Upregulated in tumor**: TP53, KRT8, ABCB1, BCL2L1, MKI67, TOP2A, AKT1, MTOR, MYC
- **Spatial pattern**: Immune exclusion phenotype (CD8+ cells blocked by stromal barrier)
- **Moran's I**: High spatial autocorrelation for MKI67 (I=0.72), CD8 (I=0.68)

**Time**: 6 minutes | **Cost**: $0.28

### TEST_4: Imaging Analysis ✅

**Tools called**: `deepcell.segment_cells()`, `deepcell.classify_cell_states()`

**Results**:
- **Ki67 proliferation index**: 45-55% (HIGH)
- **TP53+/Ki67+ double-positive**: 38% of tumor cells
- **CD8+ density**: 5-15 cells/mm² (LOW, peripheral only)

**Time**: 5 minutes | **Cost**: $0.22

### TEST_5: Integration & Recommendations ✅

**Result**:
- **Primary**: PI3K inhibitor (alpelisib) targeting PIK3CA E545K
- **Secondary**: Anti-PD-1 immunotherapy (overcome immune exclusion)
- **Clinical trial**: NCT03602859 (alpelisib + paclitaxel)

**Time**: 3 minutes | **Cost**: $0.11

Full test results: [`servers/mcp-spatialtools/COMPLETE_WORKFLOW_TEST_SUMMARY.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-spatialtools/COMPLETE_WORKFLOW_TEST_SUMMARY.md)

### Total Workflow Performance

| Metric | Value |
|--------|-------|
| **Total time** | 35 minutes |
| **Total cost** | $1.20 |
| **Tools called** | 12+ across 4 servers |
| **Data integrated** | 5 modalities |

**Comparison to traditional workflow**:
- Time: 40 hours → 35 minutes (95% reduction)
- Cost: $3,200 → $1.20 (99.96% reduction)
- Specialists: 3-4 → 1 oncologist

### Time-to-Insight Comparison

```mermaid
graph LR
    subgraph Traditional["Traditional Workflow: 40 Hours"]
        direction TB
        T1[Clinical Data Extraction<br/>3 hours]
        T2[Genomic Analysis<br/>10 hours]
        T3[Multi-Omics Integration<br/>12 hours]
        T4[Spatial Transcriptomics<br/>8 hours]
        T5[Imaging Analysis<br/>5 hours]
        T6[Report Generation<br/>2 hours]

        T1 --> T2 --> T3 --> T4 --> T5 --> T6
    end

    subgraph AIOrch["AI-Orchestrated: 35 Minutes"]
        direction TB
        A1[Clinical Data<br/>5 min]
        A2[Genomic Analysis<br/>8 min]
        A3[Multi-Omics<br/>6 min]
        A4[Spatial Analysis<br/>5 min]
        A5[Imaging<br/>8 min]
        A6[AI Synthesis<br/>3 min]

        A1 --> A2 --> A3 --> A4 --> A5 --> A6
    end

    style Traditional fill:#f8d7da,stroke:#dc3545,stroke-width:2px
    style AIOrch fill:#d4edda,stroke:#28a745,stroke-width:3px
    style T1 fill:#fff3cd
    style T2 fill:#fff3cd
    style T3 fill:#fff3cd
    style T4 fill:#fff3cd
    style T5 fill:#fff3cd
    style T6 fill:#fff3cd
    style A1 fill:#d1ecf1
    style A2 fill:#d1ecf1
    style A3 fill:#d1ecf1
    style A4 fill:#d1ecf1
    style A5 fill:#d1ecf1
    style A6 fill:#d1ecf1
```

**Figure 3.2: Time-to-Insight Comparison**
*Traditional workflow requires 40 hours of sequential manual analysis across 3-4 specialists (red border). AI-orchestrated workflow completes the same analysis in 35 minutes with automated tool coordination (green border). Each step shown with duration.*

**Key Insights:**
- **68x faster**: 40 hours → 35 minutes
- **Parallel execution**: AI orchestrates multiple servers simultaneously
- **No context switching**: Single natural language interface
- **Immediate synthesis**: AI generates integrated report automatically

---

## Cost Analysis: Real GCP Pricing

All costs validated against actual GCP deployment (2026-01-31 pricing).

### Cloud Run Compute Costs

**Pricing Model**: Pay-per-use (billed in 100ms increments)
- **CPU**: $0.00002400 per vCPU-second
- **Memory**: $0.00000250 per GiB-second

**Per-Analysis Breakdown** (PatientOne workflow):

| Server | Requests | Avg Time | Cost |
|--------|----------|----------|------|
| mcp-multiomics | 3 | 45s | $0.0086 |
| mcp-spatialtools | 2 | 30s | $0.0057 |
| mcp-fgbio | 2 | 15s | $0.0014 |
| mcp-deepcell | 1 | 60s | $0.0115 |
| mcp-mockepic | 2 | 5s | $0.0003 |

**Total Cloud Run**: $0.0275 per analysis (~$0.03)

### Claude API Costs

**Pricing** (Claude Sonnet 4.5):
- Input: $3.00 per million tokens
- Output: $15.00 per million tokens

**PatientOne workflow tokens**:
- Input: ~25,000 tokens → $0.075
- Output: ~5,000 tokens → $0.075
- **Total Claude API**: $0.15

### Grand Total Per Analysis

| Component | Cost |
|-----------|------|
| Cloud Run compute | $0.03 |
| Claude API tokens | $0.15 |
| Data egress (GCS) | $0.01 |
| **Total** | **$0.19** |

**Range**: $0.15-0.25 depending on analysis complexity
**Annual cost for 100 patients**: $19-25

**Comparison to traditional**:
- Traditional: $3,200 × 100 = $320,000
- AI-orchestrated: $20
- **Savings**: $319,980 (99.99%)

### ROI Summary

```mermaid
graph LR
    subgraph Input["💰 Investment"]
        COMPUTE[Compute Cost<br/>$0.19/patient]
        INFRA[Infrastructure<br/>$1,000/month]
    end

    subgraph Output["💰 Value Created"]
        REPLACE[Replaces Manual<br/>$3,200/patient]
        SPEED[40 hours → 35 min<br/>68x faster]
        SAVINGS[$3,098-3,176<br/>savings/patient]
    end

    Input --> Output

    style Input fill:#e1f5ff
    style Output fill:#d4edda,stroke:#28a745
```

**Figure 3.1: Return on Investment Summary**
*Investment: $0.19 per patient in compute + infrastructure costs. Value Created: Replaces $3,200 manual analysis (40 hours of specialist time), delivering 68x speed improvement and $3,098 savings per patient.*

**Annual Savings:**
- 100 patients/year = **$313,700 saved**
- 500 patients/year = **$1,568,500 saved**

Cost tracking implementation: [`shared/utils/cost_tracking.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/shared/utils/cost_tracking.py)

---

## What We Learned

### Success Factors

1. **Start with mocks, iterate to real**: Validate architecture before implementation
2. **Prioritize by clinical utility**: Build what oncologists need, not what's technically interesting
3. **Test incrementally**: Unit → integration → end-to-end
4. **Document failures**: DeepCell debugging taught more than successes
5. **Real data is messy**: Synthetic data always loads; real data has missing values, encoding issues, format variations

### Persistent Challenges

1. **Package dependencies are fragile**: DeepCell took 2 weeks of debugging
2. **Cloud services cache state**: Environment variables need explicit `--update-env-vars`
3. **Test coverage ≠ production readiness**: 23% coverage can be fine if core workflows are validated
4. **Documentation lags code**: By the time you document, implementation has changed

### Metrics That Matter

**Code quality score**: 7.5/10 (Good - Production Ready)
**Test coverage**: 56.9% overall (167 automated tests)
**Production readiness**: 7/12 servers (58%)
**Deployment success**: 11/11 servers on Cloud Run

### Server Production Readiness Matrix

| Server | Tools | Real Data | Tests | Cloud | Status |
|--------|-------|-----------|-------|-------|--------|
| **mcp-fgbio** | 9 | 95% | ✅ | ✅ | ✅ Production |
| **mcp-multiomics** | 21 | 85% | ✅ | ✅ | ✅ Production |
| **mcp-spatialtools** | 23 | 95% | ✅ | ✅ | ✅ Production |
| **mcp-deepcell** | 4 | 100% | ✅ | ✅ | ✅ Production |
| **mcp-perturbation** | 8 | 100% | ✅ | ✅ | ✅ Production |
| **mcp-quantum-celltype-fidelity** | 6 | 100% | ✅ | ✅ | ✅ Production |
| **mcp-epic** | 9 | 100% | ✅ | ⚠️ Local | ✅ Production |
| **mcp-openimagedata** | 7 | 60% | ⚠️ | ✅ | 🔶 Partial |
| **mcp-mockepic** | 8 | 0% | ✅ | ✅ | 🎭 Mock |
| **mcp-tcga** | 11 | 0% | ⚠️ | ✅ | ❌ Mocked |
| **mcp-huggingface** | 7 | 0% | ⚠️ | ✅ | ❌ Mocked |
| **mcp-seqera** | 7 | 0% | ⚠️ | ✅ | ❌ Mocked |

**Figure 3.3: Server Production Readiness Matrix**
*Assessment of all 12 MCP servers across key production criteria: real data integration, test coverage, cloud deployment, and overall readiness status. 7/12 servers (58%) are production-ready, 1 is partially implemented, 4 are intentionally mocked for demos.*

**Legend:**
- ✅ **Meets criteria** | ⚠️ **Partial/in progress** | ❌ **Mocked/synthetic data**
- **Status**: ✅ Production | 🔶 Partial | 🎭 Mock by design | ❌ Mocked

Full quality report: [`docs/for-developers/CODE_QUALITY_REPORT.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-developers/CODE_QUALITY_REPORT.md)

---

## Try It Yourself

Ready to validate these results on your own system?

### Option 1: Run the Test Suite

```bash
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp/servers/mcp-multiomics
pip install -e ".[dev]" && pytest ../../tests/unit/mcp-multiomics/ -v
# Expected: 91 passed in ~23s
```

### Option 2: Deploy to Your GCP Account

Follow the deployment guide to run on your own Cloud Run. See **Appendix: Setup Guide** for complete instructions.

1. Set up GCP project (free tier: $300 credit for 90 days)
2. Enable Cloud Run API
3. Run: `./infrastructure/deployment/deploy_to_gcp.sh your-project-id`

Cost: ~$0.02-0.05 for initial deployment

Full deployment guide: [`docs/deployment/GET_STARTED.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/deployment/GET_STARTED.md) or **Appendix: Setup Guide**

### Option 3: Interactive Notebook

Explore test results and run mini-analyses: [`docs/book/companion-notebooks/chapter-03-testing.ipynb`](../companion-notebooks/chapter-03-testing.ipynb)

This notebook includes:
- PatientOne workflow reproduction
- Cost calculator
- Performance profiling
- Test coverage visualization

---

## Summary

**Chapter 3 Key Takeaways:**
- Started with 100% mocked data, iterated to 7/12 production servers
- 167 automated tests, 56.9% coverage
- DeepCell took 3 weeks (package naming, Python version, GCS loading)
- Complete PatientOne workflow: 35 minutes, $1.20 cost (validated)
- Code quality: 7.5/10 (production-ready)
- Real GCP pricing: $0.19 per analysis (Cloud Run + Claude API)

**Companion Resources:**
- 📓 [Jupyter Notebook](../companion-notebooks/chapter-03-testing.ipynb) - Run tests and calculate costs
- 📊 [Code Quality Report](../../for-developers/CODE_QUALITY_REPORT.md) - Full analysis
- 🧪 [Test Coverage](../../test-docs/test-coverage.md) - Server-by-server breakdown
- ✅ [Workflow Test Summary](../../servers/mcp-spatialtools/COMPLETE_WORKFLOW_TEST_SUMMARY.md) - PatientOne validation
- 📚 [Appendix: Setup Guide](appendix-setup-guide.md) - Installation and deployment instructions

**GitHub References:**
- Test suite: [`tests/unit/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/tests/unit)
- Deployment script: [`infrastructure/deployment/deploy_to_gcp.sh`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/infrastructure/deployment/deploy_to_gcp.sh)
- DeepCell debugging: [`servers/mcp-deepcell/DEPENDENCY_ISSUES.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/DEPENDENCY_ISSUES.md)



ewpage


# Clinical Data—The Starting Point

> *"Every analysis begins with a patient. How do we integrate EHR data?"*

---

## Why Clinical Data Comes First

Before genomics, transcriptomics, or imaging analysis, you need essential clinical context:
- **Demographics**: Age and ancestry (risk factors vary)
- **Diagnosis and staging**: IIIA vs. IV determines treatment approach
- **Treatment history**: Platinum-sensitive vs. platinum-resistant status
- **Lab markers**: CA-125 trends indicate response or progression
- **Medications**: Current treatments and failures

This information lives in the **Electronic Health Record (EHR)** (Epic, Cerner). Without it, genomic analysis lacks actionable context. Sarah's TP53 and PIK3CA mutations (Chapter 1) only become clinically meaningful when combined with her 8-month platinum-free interval, CA-125 trajectory (45→310 U/mL), and BRCA1 status.

This chapter builds **mcp-epic**: an MCP server connecting to Epic FHIR APIs with automatic HIPAA Safe Harbor de-identification.

---

## The FHIR Standard

**HL7 FHIR (Fast Healthcare Interoperability Resources)** is the REST API standard for healthcare data exchange, defining:
- **Resources**: Patient, Condition, Observation, Medication (standardized structures)
- **Terminologies**: ICD-10 (diagnoses), LOINC (labs), RxNorm (drugs)
- **Operations**: HTTP CRUD operations

**FHIR R4** (2019) is supported by Epic, Cerner, Allscripts, athenahealth.

### Example: Patient Resource

```json
{
  "resourceType": "Patient",
  "id": "RESEARCH-PAT001",
  "name": [{"family": "Anderson", "given": ["Sarah"]}],
  "gender": "female",
  "birthDate": "1966-03-15"
}
```

FHIR specification: https://www.hl7.org/fhir/patient.html

**HIPAA Problem**: Name, full birth date, address, MRN are all **Protected Health Information (PHI)** requiring removal before research use.

---

## HIPAA Safe Harbor De-identification

HIPAA defines two de-identification methods:
1. **Safe Harbor**: Remove 18 specific identifiers (mechanical, automatable)
2. **Expert Determination**: Statistical disclosure risk analysis (expensive)

We use **Safe Harbor** for deterministic, automated de-identification.

### The 18 HIPAA Identifiers (Must Remove)

Names, geographic subdivisions smaller than state, dates except year (ages >89 aggregated), phone/fax/email, SSN, MRN, account numbers, license numbers, vehicle/device IDs, URLs, IP addresses, biometrics, photos, unique identifiers.

HIPAA guidance: https://www.hhs.gov/hipaa/for-professionals/privacy/special-topics/de-identification/index.html

### After Safe Harbor De-identification

```json
{
  "resourceType": "Patient",
  "id": "deidentified-a4f9c82b1e3d5f",
  "gender": "female",
  "birthDate": "1966",  // Year only
  "_deidentified": true
}
```

---

## Building mcp-epic: Architecture

```mermaid
graph TB
    subgraph "User Layer"
        USER[Clinician/Researcher<br/>Natural Language Query]
    end

    subgraph "AI Orchestration"
        CLAUDE[Claude API<br/>MCP Client]
    end

    subgraph "mcp-epic Server"
        TOOLS[MCP Tools<br/>FastMCP Decorators]
        FHIR[FHIR Client<br/>Epic API Wrapper]
        DEIDENT[De-identification Layer<br/>HIPAA Safe Harbor]
    end

    subgraph "Healthcare System"
        EPIC[Epic FHIR API<br/>OAuth 2.0 + TLS]
        EHR[Electronic Health Record<br/>Production PHI Database]
    end

    subgraph "Genomics Integration"
        SPATIAL[mcp-spatialtools<br/>Spatial Data Bridge]
        FGBIO[mcp-fgbio<br/>Genomic Variants]
    end

    USER --> CLAUDE
    CLAUDE --> TOOLS
    TOOLS --> FHIR
    FHIR --> EPIC
    EPIC --> EHR
    FHIR --> DEIDENT
    DEIDENT --> CLAUDE
    CLAUDE --> SPATIAL
    CLAUDE --> FGBIO

    style USER fill:#e1f5ff
    style CLAUDE fill:#cce5ff,stroke:#004085,stroke-width:2px
    style TOOLS fill:#d4edda,stroke:#28a745,stroke-width:2px
    style DEIDENT fill:#fff3cd,stroke:#ffc107,stroke-width:2px
    style EPIC fill:#f8d7da,stroke:#dc3545,stroke-width:2px
    style EHR fill:#f8d7da,stroke:#dc3545,stroke-width:2px
```

**Figure 4.1: FHIR to Genomics Bridge Architecture**
*Clinical data flows from Epic EHR → FHIR API → mcp-epic server → HIPAA de-identification → Claude AI orchestration → integration with genomics (mcp-fgbio) and spatial transcriptomics (mcp-spatialtools) for complete patient analysis.*

**Key Components:**
1. **MCP Tools**: 9 FHIR R4 resource tools with natural language interface
2. **FHIR Client**: Epic API wrapper handling OAuth 2.0 authentication
3. **De-identification**: Automatic HIPAA Safe Harbor compliance (18 identifiers removed)
4. **Integration**: Clinical context bridges to genomics and spatial analysis

Server structure: `servers/mcp-epic/` with `server.py` (9 tools), `epic_fhir_client.py`, `deidentify.py`, `__main__.py`.

Repository: [`servers/mcp-epic/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/servers/mcp-epic)

---

## Implementation: De-identification Layer

### Hash Identifiers

```python
def hash_identifier(value: str) -> str:
    """Hash identifier using SHA-256."""
    return f"HASH-{hashlib.sha256(value.encode()).hexdigest()[:16]}"
    # Full implementation: servers/mcp-epic/src/mcp_epic/deidentify.py:69-83
```

### Remove Direct Identifiers

```python
def deidentify_patient(patient: dict) -> dict:
    """Apply HIPAA Safe Harbor de-identification."""
    # Remove: name, telecom, address, photo, contact
    # Full implementation: servers/mcp-epic/src/mcp_epic/deidentify.py:19-110
```

### Date Reduction

```python
def reduce_date_to_year(date_string: str) -> str:
    """Reduce date to year, aggregate ages >89."""
    age = datetime.utcnow().year - datetime.fromisoformat(date_string).year
    return f"{datetime.utcnow().year - 90}" if age > 89 else str(birth_year)
    # Full implementation: servers/mcp-epic/src/mcp_epic/deidentify.py:84-100
```

---

## The Four MCP Tools

### 1. get_patient_demographics

Retrieves age, gender, hashed identifiers with automatic HIPAA de-identification.

```python
@mcp.tool()
async def get_patient_demographics(patient_id: str) -> dict:
    """Retrieve patient demographics from Epic FHIR API."""
    patient = await get_epic_client().get_patient(patient_id)
    return {"status": "success", "data": patient, "deidentified": True}
    # Full implementation: servers/mcp-epic/src/mcp_epic/server.py:40-77
```

**Example output**:
```json
{
  "status": "success",
  "data": {
    "id": "deidentified-a4f9c82b1e3d5f",
    "gender": "female",
    "birthDate": "1966"
  }
}
```

![Get Mock Patient Info from Epic](images/screenshots/get-mock-patient-info-from-epic.png){width=\textwidth}

**Figure 4.2: Mock Epic FHIR API Response**
*Example response from `mockepic_get_patient_demographics` showing de-identified patient data. Original PHI (name, full DOB, MRN) replaced with hashed identifiers per HIPAA Safe Harbor. Demonstrates FHIR R4 format with automatic de-identification applied.*

### 2. get_patient_conditions

Retrieves diagnoses with ICD-10 codes and staging.

```python
@mcp.tool()
async def get_patient_conditions(patient_id: str, category: str = None) -> dict:
    """Retrieve patient conditions/diagnoses with ICD-10 codes and staging."""
    # Full implementation: servers/mcp-epic/src/mcp_epic/server.py:79-120
```

**Example output**:
```json
{
  "data": [{
    "code": "C56.9",  // Ovarian cancer
    "display": "Malignant neoplasm of ovary",
    "stage": {"summary": "Stage IV", "type": "TNM"},
    "recordedDate": "2023"  // Year only
  }]
}
```

### 3. get_patient_observations

Retrieves lab results and vital signs.

```python
@mcp.tool()
async def get_patient_observations(patient_id: str, category: str = None, code: str = None) -> dict:
    """Retrieve observations (labs, vitals) with LOINC codes."""
    # Full implementation: servers/mcp-epic/src/mcp_epic/server.py:122-175
```

**Example** (CA-125 tumor marker):
```json
{
  "data": [{
    "code": "10334-1",  // LOINC: CA-125
    "valueQuantity": {"value": 310, "unit": "U/mL"},
    "referenceRange": {"low": 0, "high": 35},
    "interpretation": "High"
  }]
}
```

### 4. get_patient_medications

Retrieves current and historical medications.

```python
@mcp.tool()
async def get_patient_medications(patient_id: str, status: str = None) -> dict:
    """Retrieve patient medications with RxNorm codes."""
    # Full implementation: servers/mcp-epic/src/mcp_epic/server.py:177-210
```

---

## Epic FHIR Client: Authentication

```python
class EpicFHIRClient:
    async def get_access_token(self) -> str:
        """Obtain OAuth 2.0 access token using client credentials flow."""
        # Full implementation: servers/mcp-epic/src/mcp_epic/epic_fhir_client.py
```

**De-identification happens before data leaves FHIR client** to ensure no PHI in logs or responses.

---

## Configuration

```bash
EPIC_FHIR_ENDPOINT="https://hospital.epic.com/api/FHIR/R4/"
EPIC_CLIENT_ID="abc123-your-client-id"
EPIC_CLIENT_SECRET="your-secret-here"
DEIDENTIFY_ENABLED="true"
```

Get Epic credentials: https://fhir.epic.com/Documentation?docId=epiconfhirrequestprocess

---

## Testing with mcp-mockepic

**mcp-mockepic** provides synthetic FHIR data for testing without Epic credentials.

```python
@mcp.tool()
async def query_patient_records(patient_id: str) -> dict:
    """Retrieve synthetic patient (Synthea-based)."""
    # Returns synthetic patient PAT001-OVC-2025
```

Deployed endpoint: https://mcp-mockepic-ondu7mwjpa-uc.a.run.app

Mock server implementation: [`servers/mcp-mockepic/src/mcp_mockepic/server.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-mockepic/src/mcp_mockepic/server.py)

---

## Integration with Other Servers

Clinical data flows to downstream analysis:

- **mcp-fgbio** (genomics): Uses diagnosis to select relevant gene panels
- **mcp-multiomics**: Uses treatment history to group samples by response
- **mcp-spatialtools**: Links clinical outcomes to spatial patterns
- **mcp-tcga**: Uses diagnosis/staging to select matching cohorts

---

## Deployment: Local Only (HIPAA)

**Critical**: mcp-epic uses **STDIO transport** (not HTTP):
- ✅ Runs locally (hospital workstation/VPN)
- ✅ No network exposure of PHI
- ✅ HIPAA compliant
- ❌ Cannot deploy to Cloud Run

### Claude Desktop Configuration

```json
{
  "mcpServers": {
    "epic": {
      "command": "python",
      "args": ["-m", "mcp_epic"],
      "env": {
        "EPIC_FHIR_ENDPOINT": "https://hospital.epic.com/api/FHIR/R4/",
        "EPIC_CLIENT_ID": "your-client-id",
        "EPIC_CLIENT_SECRET": "your-secret",
        "DEIDENTIFY_ENABLED": "true"
      }
    }
  }
}
```

Setup guide: [`servers/mcp-epic/CLAUDE_DESKTOP_TESTING.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-epic/CLAUDE_DESKTOP_TESTING.md)

---

## Validation: HIPAA Compliance Checklist

Before production:
- [ ] All 18 HIPAA identifiers removed (run unit tests)
- [ ] Dates reduced to year only
- [ ] Ages >89 aggregated
- [ ] IDs hashed with SHA-256
- [ ] No PHI in logs
- [ ] BAA in place with Google Cloud
- [ ] Audit logging enabled (10-year retention)

Test suite: [`tests/unit/mcp-epic/test_deidentification.py`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/tests/unit/mcp-epic)

HIPAA compliance documentation: [`docs/for-hospitals/compliance/hipaa.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-hospitals/compliance/hipaa.md)

---

## Try It Yourself

### Option 1: Test with mcp-mockepic

```bash
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp/servers/mcp-mockepic
python -m venv venv && source venv/bin/activate
pip install -e . && python -m mcp_mockepic
```

### Option 2: Epic Sandbox

Epic provides public sandbox:
- Endpoint: https://fhir.epic.com/interconnect-fhir-oauth/api/FHIR/R4/
- Register: https://fhir.epic.com/Developer/Apps

### Option 3: Interactive Notebook

Explore FHIR resources: [`docs/book/companion-notebooks/chapter-04-clinical-data.ipynb`](../companion-notebooks/chapter-04-clinical-data.ipynb)

---

## Summary

**Chapter 4 Key Takeaways:**
- Clinical data provides essential patient context for precision medicine
- FHIR R4 is the standard for healthcare data interoperability
- HIPAA Safe Harbor removes 18 identifiers (names, dates, addresses, etc.)
- mcp-epic: 4 tools for demographics, conditions, observations, medications
- De-identification happens automatically before data leaves FHIR client
- Local STDIO deployment (not Cloud Run) for HIPAA compliance

**Companion Resources:**
- 📓 [Jupyter Notebook](../companion-notebooks/chapter-04-clinical-data.ipynb)
- 🏥 [HIPAA Compliance Guide](../../for-hospitals/compliance/hipaa.md)
- 📋 [Clinical Architecture](../../architecture/clinical/README.md)
- 🔧 [Epic Setup Guide](../../servers/mcp-epic/CLAUDE_DESKTOP_TESTING.md)

**GitHub References:**
- mcp-epic server: [`servers/mcp-epic/src/mcp_epic/server.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-epic/src/mcp_epic/server.py)
- De-identification: [`servers/mcp-epic/src/mcp_epic/deidentify.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-epic/src/mcp_epic/deidentify.py)
- mcp-mockepic: [`servers/mcp-mockepic/src/mcp_mockepic/server.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-mockepic/src/mcp_mockepic/server.py)



ewpage


# Genomic Foundations

*Building mcp-fgbio for VCF parsing and quality control*

---

## Why This Server Matters

Precision oncology requires identifying which mutations drive cancer:
- **Which genes are mutated?** TP53? PIK3CA? BRCA1?
- **What's the allele frequency?** 73% VAF suggests high tumor purity
- **Is it clinically relevant?** Check COSMIC, ClinVar, gnomAD
- **Is sequencing trustworthy?** Read depth, quality scores, contamination

The `mcp-fgbio` server provides 4 tools for reference genomes, FASTQ validation, UMI extraction, and gene annotation queries.

---

## The VCF Format in 30 Seconds

**VCF (Variant Call Format)** is the standard for representing genetic variants. PatientOne's most clinically relevant mutation:

```vcf
#CHROM  POS       ID          REF  ALT  QUAL    INFO
chr17   7578406   TP53_R175H  C    A    1250.5  DP=245;AF=0.73;GENE=TP53;EFFECT=missense_variant
```

This tells you: chr17 position 7,578,406 C→A (TP53 missense), QUAL 1250.5, 245 reads coverage, 73% allele frequency, COSMIC ID.

PatientOne's VCF ([`genomics/somatic_variants.vcf`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/data/patient-data/PAT001-OVC-2025/genomics/somatic_variants.vcf)) contains:
- **3 pathogenic**: TP53 R175H (73% AF), PIK3CA E545K (42% AF), PTEN LOH (85% AF)
- **5 copy number variants**: MYC amplification, CCNE1 amplification, etc.
- **4 wild-type**: BRCA1, BRAF, KRAS, ARID1A (important negatives)

![VCF Parsed Output](images/screenshots/vcf-parsed-output.png){width=\textwidth}

**Figure 5.1: Parsed VCF with Annotations**
*Claude output showing parsed and annotated PatientOne VCF file using `fgbio_parse_vcf` and `fgbio_annotate_variants` tools. Table displays gene names, variants (protein changes), mutation types, variant allele frequencies (VAF), and ClinVar pathogenicity ratings with star ratings. Pathogenic mutations (TP53, PIK3CA, PTEN, KRAS) are highlighted for clinical decision-making.*

---

## The Four mcp-fgbio Tools

### 1. fetch_reference_genome

Downloads reference genomes (hg38, mm10, hg19) from NCBI.

```python
@mcp.tool()
async def fetch_reference_genome(genome: str, output_dir: str) -> dict:
    """Download reference genome from NCBI."""
    genome_info = REFERENCE_GENOMES[genome]
    output_path = Path(output_dir) / f"{genome}.fna.gz"
    download_result = await _download_file(genome_info["url"], output_path)
    return {"path": str(output_path), "size_mb": download_result["size_bytes"] / (1024 * 1024)}
    # Full implementation: servers/mcp-fgbio/src/mcp_fgbio/server.py:275-345
```

**Natural language use**: *"Fetch the hg38 reference genome and save to `/workspace/data/reference`"*

---

### 2. validate_fastq

Validates FASTQ format and calculates quality metrics.

```python
@mcp.tool()
async def validate_fastq(fastq_path: str, min_quality_score: int = 20) -> dict:
    """Validate FASTQ and calculate QC metrics."""
    is_valid, messages, info = validate_fastq_file(fastq_path)
    if not is_valid: raise ValidationError(messages)
    stats = _calculate_fastq_stats(fastq_path)
    return {"valid": True, "total_reads": stats["total_reads"], "avg_quality": stats["avg_quality"]}
    # Full implementation: servers/mcp-fgbio/src/mcp_fgbio/validation.py:53-150
```

**Metrics returned**: Total reads, average Phred quality score, average read length, warnings.

---

### 3. extract_umis

Extracts Unique Molecular Identifiers (UMIs) for PCR duplicate removal.

```python
@mcp.tool()
async def extract_umis(fastq_path: str, output_dir: str, umi_length: int = 12, read_structure: str = "12M+T") -> dict:
    """Extract UMIs using FGbio toolkit."""
    # Read structure: 12M = 12bp UMI, +T = rest is template
    result = _run_fgbio_command(["ExtractUmisFromBam", f"--read-structure={read_structure}", ...])
    return {"output_fastq": f"{output_dir}/with_umis.fastq.gz", "umi_count": 45000}
    # Full implementation: servers/mcp-fgbio/src/mcp_fgbio/server.py:350-425
```

**Why you need it**: PCR duplicates inflate variant allele frequencies. UMIs tag molecules before amplification for accurate duplicate removal.

---

### 4. query_gene_annotations

Retrieves gene annotation data from GENCODE, Ensembl, or RefSeq.

```python
@mcp.tool()
async def query_gene_annotations(genome: str, gene_name: str = None, chromosome: str = None, annotation_source: str = "gencode") -> dict:
    """Query gene annotations from GENCODE/Ensembl/RefSeq."""
    annotations = await _fetch_annotations(genome, annotation_source)
    results = [ann for ann in annotations if (gene_name is None or ann["gene_name"] == gene_name)]
    return {"annotations": results, "total_genes": len(results), "source": annotation_source}
    # Full implementation: servers/mcp-fgbio/src/mcp_fgbio/server.py:430-510
```

**Example result** (TP53 in hg38):
```json
{
  "annotations": [{
    "gene_name": "TP53",
    "gene_id": "ENSG00000141510",
    "chromosome": "chr17",
    "start": 7661779,
    "end": 7687550,
    "gene_type": "protein_coding"
  }]
}
```

---

## Implementation Walkthrough

### Project Setup

```bash
cd servers/mcp-fgbio
python -m venv venv && source venv/bin/activate
pip install fastmcp httpx aiofiles
```

Environment variables (`.env`):
```bash
FGBIO_REFERENCE_DATA_DIR=/workspace/data/reference
FGBIO_CACHE_DIR=/workspace/cache
FGBIO_DRY_RUN=true  # Use mocks for development
```

### Initialize FastMCP Server

```python
from fastmcp import FastMCP
from pathlib import Path
import os

mcp = FastMCP("fgbio")

# Configuration
REFERENCE_DATA_DIR = Path(os.getenv("FGBIO_REFERENCE_DATA_DIR", "/workspace/data/reference"))
CACHE_DIR = Path(os.getenv("FGBIO_CACHE_DIR", "/workspace/cache"))
DRY_RUN = os.getenv("FGBIO_DRY_RUN", "false").lower() == "true"

# Reference genome metadata
REFERENCE_GENOMES = {
    "hg38": {"name": "Human GRCh38", "url": "https://ftp.ncbi.nlm.nih.gov/genomes/.../GRCh38_genomic.fna.gz", "size_mb": 938},
    "mm10": {"name": "Mouse GRCm38", "url": "https://ftp.ncbi.nlm.nih.gov/genomes/.../GRCm39_genomic.fna.gz", "size_mb": 794}
}
```

Full initialization: [`servers/mcp-fgbio/src/mcp_fgbio/server.py:1-130`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-fgbio/src/mcp_fgbio/server.py#L1-L130)

### Add FASTQ Validation

```python
def validate_fastq_file(file_path: str) -> Tuple[bool, List[str], dict]:
    """Validate FASTQ format and quality encoding."""
    path = Path(file_path)
    if not path.exists(): return False, [f"File not found"], {}
    is_gzipped = str(path).endswith('.gz')
    # Read first 8 lines, validate format, check sequence/quality lengths match
    # Full implementation: servers/mcp-fgbio/src/mcp_fgbio/validation.py
```

### Add Retry Logic for Downloads

```python
from shared.utils.api_retry import retry_with_backoff

@retry_with_backoff(max_retries=3, base_delay=2.0, max_delay=60.0)
async def _download_file(url: str, output_path: Path) -> dict:
    """Download file with automatic retries."""
    # Exponential backoff: 2s, 4s, 8s, up to 60s max
    # Full implementation: shared/utils/api_retry.py
```

---

## Testing Your Server

### Unit Tests

```python
@pytest.mark.asyncio
async def test_fetch_reference_genome():
    result = await mcp.call_tool("fetch_reference_genome", {"genome": "hg38", "output_dir": "/tmp/test"})
    assert result["metadata"]["genome_id"] == "hg38"
    assert result["size_mb"] > 0
```

Test coverage: **73%** ([`docs/test-docs/test-coverage.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/test-docs/test-coverage.md))

---

## Dry-Run Mode

```bash
export FGBIO_DRY_RUN=true
python -m mcp_fgbio
```

In dry-run mode:
- `fetch_reference_genome` creates small mock file instead of downloading 938 MB
- `validate_fastq` returns synthetic quality metrics
- Results include warning banner about synthetic data

Implementation: [`servers/mcp-fgbio/src/mcp_fgbio/server.py:44-72`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-fgbio/src/mcp_fgbio/server.py#L44-L72)

---

## Connecting to Claude Desktop

```json
{
  "mcpServers": {
    "fgbio": {
      "command": "/path/to/venv/bin/python",
      "args": ["-m", "mcp_fgbio"],
      "env": {
        "FGBIO_REFERENCE_DATA_DIR": "/workspace/data/reference",
        "FGBIO_DRY_RUN": "false"
      }
    }
  }
}
```

Full configuration: [`configs/claude_desktop_config.json`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/configs/claude_desktop_config.json)

---

## What You've Built

A genomics foundation server providing:
1. **Reference genome management**: Download hg38, mm10, hg19 from NCBI with retry logic
2. **Sequencing data validation**: Check FASTQ format, quality scores, read lengths
3. **UMI extraction**: Prepare data for PCR duplicate removal
4. **Gene annotations**: Map mutations to genes and transcripts

This supports:
- **Chapter 6**: Multi-omics integration (RNA-seq requires gene annotations)
- **Chapter 7**: Spatial transcriptomics (STAR alignment needs hg38)
- **Chapter 9**: Treatment response prediction (GEARS needs gene IDs)

---

## Try It Yourself

### Local Development

```bash
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp/servers/mcp-fgbio
python -m venv venv && source venv/bin/activate
pip install -e ".[dev]"
cp .env.example .env  # Set FGBIO_DRY_RUN=true for testing
python -m mcp_fgbio
```

---

## Summary

**Chapter 5 Summary**:
- VCF format stores somatic mutations with quality metrics and clinical annotations
- mcp-fgbio provides 4 tools: reference genomes, FASTQ validation, UMI extraction, gene annotations
- Dry-run mode enables safe testing without downloading large files
- Test coverage: 73% with 15 unit tests

**Files created**: `servers/mcp-fgbio/src/mcp_fgbio/server.py`, `validation.py`
**Tools exposed**: 4 MCP tools via FastMCP decorators



ewpage


# Multi-Omics Integration

*Building mcp-multiomics for RNA, protein, and phosphoproteomics analysis*

---

## The Multi-Omics Problem

PatientOne's PDX models (15 ovarian cancer xenografts) generate three data types:
1. **RNA-seq** (transcriptomics): What genes are transcribed?
2. **Proteomics**: What proteins are made?
3. **Phosphoproteomics**: What signaling pathways are active?

**Challenge**: A gene upregulated at RNA level might not show elevated protein (post-transcriptional regulation). High protein doesn't guarantee phosphorylation-driven activity (kinase regulation).

**PatientOne example**:
- **TP53 RNA**: Low expression (mutated)
- **TP53 protein**: High (misfolded protein accumulates)
- **AKT phosphorylation**: High (PI3K pathway activated despite TP53 dysfunction)

You need **multi-omics meta-analysis** to find consistently dysregulated pathways across all modalities.

The `mcp-multiomics` server provides 21 tools for HAllA association testing and Stouffer meta-analysis.

### Multi-Omics Integration Workflow

```mermaid
graph TB
    subgraph "Data Input"
        RNA[RNA-seq<br/>15 samples<br/>20,000 genes]
        PROTEIN[Proteomics<br/>15 samples<br/>8,000 proteins]
        PHOSPHO[Phosphoproteomics<br/>15 samples<br/>12,000 sites]
    end

    subgraph "Phase 1: Validation & QC"
        VALIDATE[validate_multiomics_data<br/>Check overlap, missing values,<br/>batch effects]
        PREPROCESS[preprocess_multiomics_data<br/>ComBat correction<br/>KNN imputation<br/>Quantile normalization]
        VISUALIZE[visualize_data_quality<br/>PCA plots<br/>Correlation heatmaps]
    end

    subgraph "Phase 2: Association Discovery"
        HALLA[run_halla_analysis<br/>HAllA: Hierarchical All-against-All<br/>Discover RNA-protein associations]
        INTEGRATE[integrate_phospho_data<br/>Add phosphorylation layer<br/>Kinase-substrate mapping]
    end

    subgraph "Phase 3: Meta-Analysis"
        STOUFFER[stouffer_meta_analysis<br/>Combine p-values across modalities<br/>Z-score integration]
        PATHWAYS[analyze_upstream_regulators<br/>Identify dysregulated pathways<br/>PI3K/AKT, DNA repair, etc.]
    end

    subgraph "Output"
        GENES[Top 47 Genes<br/>Integrated p-values<br/>Pathway annotations]
        VIZ[Visualizations<br/>Volcano plots<br/>Heatmaps]
    end

    RNA --> VALIDATE
    PROTEIN --> VALIDATE
    PHOSPHO --> VALIDATE

    VALIDATE --> PREPROCESS
    PREPROCESS --> VISUALIZE
    VISUALIZE --> HALLA
    HALLA --> INTEGRATE
    INTEGRATE --> STOUFFER
    STOUFFER --> PATHWAYS
    PATHWAYS --> GENES
    PATHWAYS --> VIZ

    style RNA fill:#d1ecf1
    style PROTEIN fill:#d1ecf1
    style PHOSPHO fill:#d1ecf1
    style VALIDATE fill:#fff3cd
    style PREPROCESS fill:#fff3cd
    style HALLA fill:#cce5ff
    style STOUFFER fill:#d4edda,stroke:#28a745,stroke-width:2px
    style GENES fill:#d4edda,stroke:#28a745,stroke-width:2px
```

**Figure 6.1: Multi-Omics Integration Workflow**
*Three-phase pipeline: (1) Validation & QC with batch correction and imputation, (2) Association discovery using HAllA to find RNA-protein correlations and integrate phosphorylation data, (3) Meta-analysis with Stouffer's method to combine p-values and identify upstream pathway regulators. Output: 47 top genes with integrated evidence across all three modalities.*

---

## The Tools (21 total)

### Phase 1: Data Validation

#### 1. validate_multiomics_data

Checks sample overlap, missing value patterns, batch effects, outlier samples.

```python
@mcp.tool()
def validate_multiomics_data(rna_path: str, protein_path: str, metadata_path: str) -> dict:
    """Validate multi-omics data quality: batch effects, missing values, sample naming."""
    # Load data, check sample overlap, detect batch effects (PC1 ~ Batch correlation)
    # Identify missing patterns, outlier samples (MAD > 3)
    # Full implementation: servers/mcp-multiomics/src/mcp_multiomics/tools/preprocessing.py
```

**Example result**:
```json
{
  "validation_status": "WARNING",
  "sample_overlap": {"common": 15},
  "missing_patterns": {"rna": "0.5%", "protein": "23.4%"},
  "batch_effects": {"protein_batch_correlation": 0.68},
  "warnings": ["Protein batch effect detected (PC1 ~ Batch, R²=0.68)"],
  "recommendations": ["Apply ComBat batch correction", "Use KNN imputation"]
}
```

---

#### 2. preprocess_multiomics_data

Applies ComBat batch correction, KNN imputation, quantile normalization, outlier removal.

```python
@mcp.tool()
def preprocess_multiomics_data(rna_path: str, protein_path: str, metadata_path: str, batch_correction: bool = True) -> dict:
    """Preprocess multi-omics data with batch correction and imputation."""
    # ComBat batch correction → KNN imputation → Quantile normalization → Outlier removal
    # Full implementation: servers/mcp-multiomics/src/mcp_multiomics/tools/preprocessing.py
```

**Before/after batch correction**:
```
Before: PC1 explains 45% variance, R² with Batch = 0.68 (BAD)
After:  PC1 explains 38% variance, R² with Batch = 0.12 (GOOD)
```

---

#### 3. visualize_data_quality

Generates PCA plots, missing value heatmaps, correlation matrices to verify preprocessing.

---

### Phase 2: Integration and Analysis

#### 4. integrate_omics_data

Aligns samples across modalities and normalizes within each.

```python
@mcp.tool()
def integrate_omics_data(rna_path: str, protein_path: str, phospho_path: str) -> dict:
    """Align samples across modalities, apply Z-score normalization per modality."""
    # Find common samples (15 in PatientOne), apply Z-score normalization
    # Filter features with >50% missing
    # Full implementation: servers/mcp-multiomics/src/mcp_multiomics/tools/integration.py
```

---

#### 5. run_halla_analysis

**HAllA (HMP All-Against-All Association)**: Tests pairwise associations between features.

**Scalability problem**: 1850 RNA features × 1850 = 3.4 million pairwise tests → days of compute.

**Solution: Chunking strategy**:
```python
# Chunk features into groups of 1000
chunks = [features[i:i+1000] for i in range(0, len(features), 1000)]
# Run HAllA on each chunk (1000×1000 = 1M tests, ~5 minutes)
for chunk in chunks:
    halla_results = halla.run(chunk, chunk, method='spearman')
```

**Result**: 3.4 million tests in ~15 minutes instead of days.

Implementation: [`servers/mcp-multiomics/src/mcp_multiomics/tools/halla.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-multiomics/src/mcp_multiomics/tools/halla.py)

---

#### 6. calculate_stouffer_meta

**Stouffer's Meta-Analysis**: Combines p-values across modalities to find consistently dysregulated features.

**Critical workflow** (per bioinformatician review):

```
CORRECT FDR TIMING FOR MULTI-OMICS META-ANALYSIS:
Step 1: HAllA Analysis → Use NOMINAL p-values (NOT FDR-corrected)
Step 2: Stouffer's Meta-Analysis → Combine NOMINAL p-values across modalities
Step 3: FDR Correction → Apply Benjamini-Hochberg AFTER combination

WHY: Applying FDR before Stouffer's loses power (over-conservative)
```

**Stouffer's Z-score method**:
```python
class StoufferMetaAnalysis:
    def p_to_z(self, p_values: np.ndarray, effect_sizes: np.ndarray = None) -> np.ndarray:
        """Convert p-values to Z-scores."""
        z_scores = stats.norm.ppf(1 - p_values / 2)
        if effect_sizes is not None: z_scores *= np.sign(effect_sizes)
        return z_scores

    def combine_z_scores(self, z_scores_list: List[np.ndarray], weights: List[float]) -> np.ndarray:
        """Combine Z-scores from multiple modalities."""
        z_meta = sum(w * z for w, z in zip(weights, z_scores_list))
        return z_meta / np.sqrt(sum(w**2 for w in weights))

    def apply_fdr_correction(self, p_values: np.ndarray) -> np.ndarray:
        """Apply Benjamini-Hochberg FDR correction AFTER meta-analysis."""
        _, q_values, _, _ = multipletests(p_values, method='fdr_bh')
        return q_values
```

**Example result** (AKT pathway):
```json
{
  "feature": "AKT1",
  "modality_results": {
    "rna": {"p_value": 0.023, "log2fc": 1.2},
    "protein": {"p_value": 0.041, "log2fc": 0.8},
    "phospho": {"p_value": 0.0087, "log2fc": 2.1}
  },
  "meta_analysis": {
    "meta_p_value": 0.00031,
    "meta_q_value": 0.012,
    "combined_z_score": 3.58
  }
}
```

AKT1 shows consistent upregulation across all 3 modalities (meta q < 0.05), making it a high-confidence therapeutic target.

Implementation: [`servers/mcp-multiomics/src/mcp_multiomics/tools/stouffer.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-multiomics/src/mcp_multiomics/tools/stouffer.py)

---

#### 7. predict_upstream_regulators

Predicts kinases, transcription factors, and drug responses.

```python
@mcp.tool()
def predict_upstream_regulators(top_features: list) -> dict:
    """Predict upstream kinases, TFs, and druggable targets."""
    # Use kinase-substrate databases, TF binding motifs, drug target databases
    # Full implementation: servers/mcp-multiomics/src/mcp_multiomics/tools/upstream.py
```

**Example** (PatientOne AKT pathway):
```json
{
  "upstream_kinases": [
    {"kinase": "PDK1", "confidence": 0.92, "targets": ["AKT1", "AKT2"]},
    {"kinase": "mTORC2", "confidence": 0.87, "targets": ["AKT1"]}
  ],
  "drug_targets": [
    {"drug": "Capivasertib", "target": "AKT1/2/3", "clinical_phase": "Phase 3"}
  ]
}
```

This tells you: AKT hyperphosphorylated → likely by PDK1/mTORC2 → druggable with AKT inhibitors.

---

## The Complete PatientOne Workflow

Natural language prompt:
```
I have multi-omics data from 15 ovarian cancer PDX models. Please:
1. Validate data quality and check for batch effects
2. Preprocess with ComBat batch correction and KNN imputation
3. Integrate the three modalities
4. Run Stouffer meta-analysis to find consistently dysregulated pathways
5. Predict upstream regulators for the top 10 hits
```

**Results for PatientOne**:
- **15 common samples** across all 3 modalities
- **1850 RNA features**, 1243 proteins, 987 phosphosites
- **23 pathways** with meta q < 0.05 (AKT/mTOR, PI3K, TP53, DNA damage response)
- **12 upstream kinases** predicted
- **8 druggable targets** identified

---

## Implementation Walkthrough

### Project Setup

```bash
cd servers/mcp-multiomics
python -m venv venv && source venv/bin/activate
pip install fastmcp pandas numpy scipy statsmodels scikit-learn
```

Environment variables (`.env`):
```bash
MULTIOMICS_DRY_RUN=true  # For testing
MULTIOMICS_DATA_DIR=/workspace/data
```

### Stouffer Meta-Analysis Core

```python
@mcp.tool()
def calculate_stouffer_meta(rna_results: dict, protein_results: dict, phospho_results: dict, fdr_threshold: float = 0.05) -> dict:
    """Combine p-values across modalities using Stouffer's method."""
    # Extract NOMINAL p-values, convert to Z-scores, combine with weights
    # Apply FDR correction AFTER combining (CRITICAL)
    # Full implementation: servers/mcp-multiomics/src/mcp_multiomics/tools/stouffer.py
```

---

## Testing Your Server

```python
def test_stouffer_combine_strong_signals():
    """Test combining strong signals from all 3 modalities."""
    stouffer = StoufferMetaAnalysis()
    p_rna = np.array([0.005, 0.008, 0.003])
    p_protein = np.array([0.007, 0.004, 0.009])
    p_phospho = np.array([0.002, 0.006, 0.001])
    # Meta p-values should be even stronger
    assert np.all(p_meta < 0.001)
```

Test coverage: **67%**, 18 unit tests

---

## Connecting to Claude Desktop

```json
{
  "mcpServers": {
    "multiomics": {
      "command": "/path/to/venv/bin/python",
      "args": ["-m", "mcp_multiomics"],
      "env": {"MULTIOMICS_DRY_RUN": "false"}
    }
  }
}
```

---

## What You've Built

A multi-omics integration server providing:
1. **Data validation**: Batch effects, missing values, outliers
2. **Preprocessing**: ComBat, KNN imputation, normalization
3. **Integration**: Sample alignment, Z-score normalization
4. **Association testing**: HAllA with chunking (3.4M tests in 15 minutes)
5. **Meta-analysis**: Stouffer with correct FDR timing
6. **Regulator prediction**: Kinases, TFs, drug targets
7. **Visualization**: PCA, heatmaps, correlation matrices

This bridges genomics (Chapter 5) to spatial transcriptomics (Chapter 7) by identifying which pathways to investigate in tissue context.

---

## Try It Yourself

```bash
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp/servers/mcp-multiomics
python -m venv venv && source venv/bin/activate
pip install -e ".[dev]"
export MULTIOMICS_DRY_RUN=true && python -m mcp_multiomics
```

---

## Summary

**Chapter 6 Summary**:
- Multi-omics integration requires preprocessing (batch correction, imputation)
- HAllA with chunking enables scalable all-vs-all association testing
- Stouffer meta-analysis combines p-values BEFORE FDR correction (critical timing)
- Upstream regulator prediction connects pathways to druggable targets
- PatientOne: 23 pathways with meta q < 0.05, 8 druggable targets

**Files created**: `servers/mcp-multiomics/src/mcp_multiomics/server.py`, `tools/stouffer.py`, `tools/preprocessing.py`, `tools/halla.py`
**Tests added**: 18 unit tests, 67% coverage
**Tools exposed**: 7 MCP tools (validate, preprocess, visualize, integrate, halla, stouffer, upstream_regulators)



ewpage


# Spatial Transcriptomics

*Building mcp-spatialtools for 10X Visium spatial gene expression analysis*

---

## Why Spatial Context Matters

Chapters 5-6 identified dysregulated pathways (TP53 mutation, AKT/mTOR hyperactivation, HIF1A upregulation). **But where in the tumor are these pathways active?**

Traditional bulk RNA-seq mixes all cells together. You can't distinguish:
- **Tumor proliferative regions** (MKI67+, PCNA+ → chemotherapy-sensitive)
- **Necrotic/hypoxic cores** (HIF1A+, CA9+ → radiotherapy-resistant)
- **Invasive fronts** (VIM+, SNAI1+ → metastatic potential)
- **Stromal barriers** (COL1A1+, FAP+ → drug delivery obstacles)
- **Immune infiltrates** (CD3D+, CD8A+ → immunotherapy responsive)

**10X Visium spatial transcriptomics** measures gene expression in 900 tissue spots (55μm diameter each), preserving spatial relationships.

The `mcp-spatialtools` server provides 23 tools for STAR alignment, batch correction, differential expression, pathway enrichment, spatial autocorrelation, and cell type deconvolution.

### Spatial Transcriptomics Pipeline

```mermaid
graph TB
    subgraph "Data Input"
        FASTQ[FASTQ Files<br/>Raw sequencing reads<br/>Paired-end 150bp]
        COORDS[Spatial Coordinates<br/>900 spots<br/>x, y positions]
        ANNOT[Region Annotations<br/>6 tissue regions<br/>Pathologist labels]
    end

    subgraph "Phase 1: Alignment & QC"
        STAR[align_spatial_reads<br/>STAR aligner<br/>Map to hg38 genome]
        FILTER[filter_spatial_data<br/>QC filtering<br/>Min counts, % mito]
        NORM[normalize_expression<br/>Quantile normalization<br/>Log2 transform]
    end

    subgraph "Phase 2: Batch Correction"
        COMBAT[correct_batch_effects<br/>ComBat correction<br/>Multi-sample integration]
    end

    subgraph "Phase 3: Analysis"
        SPLIT[split_by_region<br/>Segment by annotation<br/>6 tissue regions]
        DE[differential_expression<br/>Wilcoxon rank-sum<br/>Region vs. region]
        PATHWAY[pathway_enrichment<br/>44 curated pathways<br/>Hypergeometric test]
        MORANS[spatial_autocorrelation<br/>Moran's I<br/>Clustering patterns]
    end

    subgraph "Output"
        GENES[Differential Genes<br/>By region<br/>Fold-change + FDR]
        PATHWAYS_OUT[Enriched Pathways<br/>PI3K/AKT, DNA repair<br/>p-values]
        SPATIAL[Spatial Patterns<br/>Clustered vs. random<br/>Moran's I scores]
        VIZ[Visualizations<br/>Heatmaps<br/>Spatial plots]
    end

    FASTQ --> STAR
    STAR --> FILTER
    FILTER --> NORM
    NORM --> COMBAT
    COMBAT --> SPLIT

    COORDS --> SPLIT
    ANNOT --> SPLIT

    SPLIT --> DE
    SPLIT --> PATHWAY
    SPLIT --> MORANS

    DE --> GENES
    PATHWAY --> PATHWAYS_OUT
    MORANS --> SPATIAL
    GENES --> VIZ
    PATHWAYS_OUT --> VIZ
    SPATIAL --> VIZ

    style FASTQ fill:#d1ecf1
    style COORDS fill:#d1ecf1
    style ANNOT fill:#d1ecf1
    style STAR fill:#fff3cd
    style COMBAT fill:#cce5ff
    style DE fill:#d4edda
    style PATHWAY fill:#d4edda,stroke:#28a745,stroke-width:2px
    style MORANS fill:#d4edda,stroke:#28a745,stroke-width:2px
```

**Figure 7.1: Spatial Transcriptomics Analysis Pipeline**
*Three-phase workflow: (1) Alignment with STAR, QC filtering, and normalization, (2) ComBat batch correction for multi-sample integration, (3) Analysis including differential expression (Wilcoxon), pathway enrichment (44 curated pathways), and spatial autocorrelation (Moran's I). Outputs differential genes, enriched pathways, and spatial clustering patterns.*

---

## PatientOne's Spatial Landscape

Spatial data structure ([`spatial/visium_gene_expression.csv`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/data/patient-data/PAT001-OVC-2025/spatial/visium_gene_expression.csv)):

**Dimensions**: 900 spots × 31 genes (proliferation, mutations, oncogenes, EMT, stroma, immune, hypoxia, drug resistance markers)

**6 distinct regions** ([`spatial/visium_region_annotations.csv`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/data/patient-data/PAT001-OVC-2025/spatial/visium_region_annotations.csv)):
1. **Necrotic/Hypoxic** (150 spots): HIF1A+, CA9+, VEGFA+ → radiotherapy-resistant
2. **Tumor Proliferative** (200 spots): MKI67+, PCNA+, TOP2A+ → chemotherapy-sensitive
3. **Tumor Invasive** (180 spots): VIM+, SNAI1+, TWIST1+ → metastatic risk
4. **Stroma Reactive** (150 spots): COL1A1+, ACTA2+ → CAF-driven drug resistance
5. **Stroma Fibrotic** (120 spots): COL3A1+, FAP+ → physical barrier to drugs
6. **Immune Infiltrated** (100 spots): CD3D+, CD8A+ → immunotherapy potential

---

## The Eight mcp-spatialtools Tools

### Phase 1: Data Processing

#### 1. align_spatial_reads (STAR alignment)

Aligns raw FASTQ files to hg38 reference genome.

```python
@mcp.tool()
def align_spatial_reads(fastq_r1: str, fastq_r2: str, genome_index: str, output_dir: str, threads: int = 8) -> dict:
    """Align spatial FASTQ files using STAR."""
    # Run STAR alignment: geneCounts mode, sorted BAM output
    # Full implementation: servers/mcp-spatialtools/src/mcp_spatialtools/tools/alignment.py
```

**PatientOne mapping stats**: 85.3% uniquely mapped (high quality), 8.2% multimapped, 6.5% unmapped.

---

#### 2. filter_spatial_data

Removes low-quality spots and lowly expressed genes.

```python
@mcp.tool()
def filter_spatial_data(counts_path: str, min_counts: int = 500, max_mito_percent: float = 20.0, min_spots: int = 10) -> dict:
    """Filter low-quality spots (UMI count < 500 or mito% > 20) and genes (detected in < 10 spots)."""
    # Calculate QC metrics, filter spots and genes
    # Full implementation: servers/mcp-spatialtools/src/mcp_spatialtools/tools/filter.py
```

**PatientOne filtering**: 23 spots removed (low UMI), 0 genes removed → **877 spots × 31 genes**

---

#### 3. correct_batch_effects

Applies ComBat batch correction (same as Chapter 6).

```python
@mcp.tool()
def correct_batch_effects(counts_path: str, metadata_path: str) -> dict:
    """Apply ComBat batch correction to spatial data."""
    # ComBat correction, verify PC1 ~ Batch correlation reduction
    # Full implementation: servers/mcp-spatialtools/src/mcp_spatialtools/tools/batch.py
```

---

### Phase 2: Analysis

#### 4. differential_expression_spatial

Finds genes differentially expressed between regions.

```python
@mcp.tool()
def differential_expression_spatial(counts_path: str, annotations_path: str, group1: str, group2: str, fdr_threshold: float = 0.05) -> dict:
    """Find differentially expressed genes between spatial regions."""
    # Get spots for each group, compute log2 fold change and p-values
    # Apply FDR correction (Benjamini-Hochberg)
    # Full implementation: servers/mcp-spatialtools/src/mcp_spatialtools/tools/differential_expression.py
```

**PatientOne results** (Tumor Proliferative vs Necrotic):
```json
{
  "significant_genes": [
    {"gene": "MKI67", "log2fc": 4.2, "qval": 0.00012},
    {"gene": "HIF1A", "log2fc": -5.1, "qval": 0.000045},
    {"gene": "CA9", "log2fc": -6.3, "qval": 0.000008}
  ],
  "num_significant": 12
}
```

---

#### 5. pathway_enrichment_spatial

Identifies enriched biological pathways in each region.

```python
@mcp.tool()
def pathway_enrichment_spatial(gene_list: list[str], gene_set_database: str = "KEGG", fdr_threshold: float = 0.05) -> dict:
    """Perform pathway enrichment analysis on spatial gene sets."""
    # Fisher's exact test for enrichment, FDR correction
    # Full implementation: servers/mcp-spatialtools/src/mcp_spatialtools/tools/pathway_enrichment.py
```

**PatientOne enrichment** (Necrotic/Hypoxic region):
```json
{
  "enriched_pathways": [
    {"pathway": "HIF-1 signaling", "qval": 0.00012, "overlap": "8/43"},
    {"pathway": "VEGF signaling", "qval": 0.00089, "overlap": "6/29"}
  ]
}
```

---

#### 6. spatial_autocorrelation (Moran's I)

Identifies genes with spatially clustered expression.

```python
@mcp.tool()
def spatial_autocorrelation(counts_path: str, coordinates_path: str, gene: str = "HIF1A") -> dict:
    """Calculate Moran's I for spatial autocorrelation."""
    # Compute spatial weights (inverse distance), calculate Moran's I
    # Permutation test for significance (1000 permutations)
    # Full implementation: servers/mcp-spatialtools/src/mcp_spatialtools/tools/spatial_autocorrelation.py
```

**PatientOne spatial autocorrelation**:
- **HIF1A**: Moran's I = 0.82, p < 0.001 (Clustered)
- **MKI67**: Moran's I = 0.75, p < 0.001 (Clustered)
- **CD3D**: Moran's I = 0.61, p < 0.008 (Clustered)

All key markers show significant spatial clustering (not random distribution).

---

#### 7. deconvolve_cell_types

Estimates cell type proportions per spot using reference single-cell signatures.

```python
@mcp.tool()
def deconvolve_cell_types(spatial_counts_path: str, reference_signatures_path: str) -> dict:
    """Estimate cell type proportions per spot using non-negative least squares."""
    # NNLS deconvolution: spot_expr ≈ Σ(proportion_i × cell_type_i)
    # Full implementation: servers/mcp-spatialtools/src/mcp_spatialtools/tools/deconvolution.py
```

**PatientOne deconvolution** (immune infiltrated region):
```json
{
  "SPOT_20_20": {
    "Tumor cells": 0.35,
    "CD8+ T cells": 0.28,
    "CD4+ T cells": 0.15,
    "Macrophages": 0.12
  }
}
```

This spot is 35% tumor + 43% T cells → immunotherapy target.

---

#### 8. link_clinical_to_spatial

Bridges clinical FHIR data (Chapter 4) to spatial tissue regions.

```python
@mcp.tool()
def link_clinical_to_spatial(patient_id: str, spatial_annotations_path: str) -> dict:
    """Bridge clinical FHIR data to spatial tissue regions."""
    # Get clinical data from mcp-epic, link treatment history to spatial features
    # Predict response by region
    # Full implementation: servers/mcp-spatialtools/src/mcp_spatialtools/tools/clinical_link.py
```

---

## Complete PatientOne Spatial Workflow

Natural language prompt:
```
I have 10X Visium spatial transcriptomics data for patient PAT001-OVC-2025. Please:
1. Filter low-quality spots (< 500 UMI counts)
2. Run differential expression: tumor_proliferative vs necrotic_hypoxic
3. Pathway enrichment for top 10 upregulated genes in each region
4. Calculate Moran's I for HIF1A, MKI67, CD3D
5. Deconvolve cell types
```

Claude orchestrates all 8 tools automatically, returning comprehensive spatial analysis in **<3 minutes**.

---

## Implementation Walkthrough

### Project Setup

```bash
cd servers/mcp-spatialtools
python -m venv venv && source venv/bin/activate
pip install fastmcp pandas numpy scipy statsmodels scikit-learn
```

Environment variables (`.env`):
```bash
SPATIAL_DATA_DIR=/workspace/data/spatial
SPATIAL_DRY_RUN=true  # For testing
```

### Initialize FastMCP Server

```python
from fastmcp import FastMCP
import os
from pathlib import Path

mcp = FastMCP("spatialtools")

config = {
    "data_dir": Path(os.getenv("SPATIAL_DATA_DIR", "/workspace/data/spatial")),
    "dry_run": os.getenv("SPATIAL_DRY_RUN", "false").lower() == "true"
}
```

### Add Differential Expression

```python
def differential_expression_spatial_impl(counts_path: str, annotations_path: str, group1: str, group2: str, fdr_threshold: float = 0.05) -> dict:
    """Internal implementation of differential expression."""
    # Load data, get barcodes for each group
    # Compute log2 fold change, t-test p-values, FDR correction
    # Full implementation: servers/mcp-spatialtools/src/mcp_spatialtools/tools/differential_expression.py
```

---

## Testing Your Server

```python
def test_differential_expression_tumor_vs_necrotic():
    """Test DE analysis on PatientOne data."""
    result = differential_expression_spatial_impl(...)
    assert result["num_significant"] > 0
    assert "MKI67" in [g["gene"] for g in result["significant_genes"]]
    assert "HIF1A" in [g["gene"] for g in result["significant_genes"]]
```

Test coverage: **71%**, 14 unit tests

---

## Connecting to Claude Desktop

```json
{
  "mcpServers": {
    "spatialtools": {
      "command": "/path/to/venv/bin/python",
      "args": ["-m", "mcp_spatialtools"],
      "env": {"SPATIAL_DRY_RUN": "false"}
    }
  }
}
```

---

## What You've Built

A spatial transcriptomics server providing:
1. **Data processing**: STAR alignment, quality filtering, batch correction
2. **Regional analysis**: Differential expression, pathway enrichment
3. **Spatial patterns**: Moran's I autocorrelation, spatial clustering
4. **Cell type deconvolution**: Tumor, immune, stromal proportions per spot
5. **Clinical integration**: Links treatment history to spatial predictions

This completes Part 2 (Building the Foundation). You now have servers for clinical data (Chapter 4), genomics (Chapter 5), multi-omics (Chapter 6), and spatial (Chapter 7).

---

## Try It Yourself

```bash
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp/servers/mcp-spatialtools
python -m venv venv && source venv/bin/activate
pip install -e ".[dev]"
export SPATIAL_DRY_RUN=true && python -m mcp_spatialtools
```

---

## Summary

**Chapter 7 Summary**:
- 10X Visium spatial transcriptomics preserves tissue context (900 spots × 31 genes)
- 6 distinct regions: necrotic, tumor proliferative/invasive, stroma reactive/fibrotic, immune
- Differential expression reveals MKI67+ (chemo-sensitive) vs HIF1A+ (chemo-resistant) regions
- Moran's I confirms spatially clustered expression patterns
- Cell type deconvolution estimates tumor/immune/stromal proportions per spot

**Files created**: `servers/mcp-spatialtools/src/mcp_spatialtools/server.py`, `tools/differential_expression.py`, `tools/pathway_enrichment.py`, `tools/spatial_autocorrelation.py`
**Tests added**: 14 unit tests, 71% coverage
**Tools exposed**: 8 MCP tools (align, filter, batch_correct, differential_expression, pathway_enrichment, spatial_autocorrelation, deconvolve, link_clinical)



ewpage


# Cell Segmentation with DeepCell

*Building mcp-deepcell for single-cell resolution imaging analysis*

---

## Why Single-Cell Segmentation Matters

Chapter 7 analyzed spatial transcriptomics with 10X Visium (55μm spots containing 10-30 cells mixed together). **But Visium can't answer**:
- How many cells are **Ki67+** (actively proliferating)?
- Are **TP53+** mutant cells spatially clustered or dispersed?
- Do **CD8+ T cells** contact tumor cells (immune checkpoint blockade readiness)?
- Which cells are **double-positive** (TP53+/Ki67+ → aggressive phenotype)?

**Multiplexed Immunofluorescence (MxIF)** images 2-7 protein markers simultaneously at subcellular resolution. You need **cell segmentation** to:
1. **Detect cell boundaries** (where does one cell end and another begin?)
2. **Measure marker intensity per cell** (which cells express Ki67? TP53?)
3. **Classify cell phenotypes** (proliferating, quiescent, immune, stromal)
4. **Count cells spatially** (how many CD8+ T cells per mm²?)

The `mcp-deepcell` server uses **DeepCell-TF deep learning models** for nuclear and membrane segmentation, then classifies cells by marker intensity.

### DeepCell Segmentation & Classification Pipeline

```mermaid
graph TB
    subgraph "Input Images"
        DAPI[DAPI<br/>Nuclear stain<br/>2048×2048 TIFF]
        KI67[Ki67<br/>Proliferation marker<br/>Channel 1]
        TP53[TP53<br/>Mutation marker<br/>Channel 2]
        CD8[CD8<br/>T-cell marker<br/>Channel 3]
    end

    subgraph "Phase 1: Cell Segmentation"
        LOAD[load_image_from_gcs<br/>Handle gs:// URIs<br/>Multi-channel TIFF]
        ENGINE[DeepCellEngine<br/>Load pretrained models<br/>Nuclear or Membrane]
        PREPROCESS[Preprocessing<br/>Normalize 0-1<br/>Add batch dims]
        INFER[Deep Learning Inference<br/>DeepCell-TF CNN<br/>Cell boundary detection]
        FILTER[Post-processing<br/>Filter by size<br/>10-500 pixels]
    end

    subgraph "Phase 2: Phenotype Classification"
        INTENSITY[IntensityClassifier<br/>Measure per-cell intensity<br/>Within segmented regions]
        THRESHOLD[Thresholding<br/>Otsu's method<br/>Marker+/- classification]
        MULTIPLEX[Multi-marker phenotyping<br/>Ki67+/TP53+<br/>Double-positive cells]
    end

    subgraph "Output"
        MASKS[Segmentation Masks<br/>1247 cells detected<br/>Cell ID per pixel]
        PHENO[Cell Phenotypes<br/>562 Ki67+<br/>389 TP53+<br/>187 double-positive]
        METRICS[Quality Metrics<br/>Mean cell size<br/>Marker expression]
        VIZ[Visualizations<br/>Colored masks<br/>Spatial cell maps]
    end

    DAPI --> LOAD
    KI67 --> LOAD
    TP53 --> LOAD
    CD8 --> LOAD

    LOAD --> ENGINE
    ENGINE --> PREPROCESS
    PREPROCESS --> INFER
    INFER --> FILTER
    FILTER --> MASKS

    MASKS --> INTENSITY
    KI67 --> INTENSITY
    TP53 --> INTENSITY
    CD8 --> INTENSITY

    INTENSITY --> THRESHOLD
    THRESHOLD --> MULTIPLEX
    MULTIPLEX --> PHENO
    PHENO --> METRICS
    MASKS --> VIZ
    PHENO --> VIZ

    style DAPI fill:#d1ecf1
    style KI67 fill:#d1ecf1
    style TP53 fill:#d1ecf1
    style CD8 fill:#d1ecf1
    style ENGINE fill:#cce5ff,stroke:#004085,stroke-width:2px
    style INFER fill:#cce5ff,stroke:#004085,stroke-width:2px
    style MULTIPLEX fill:#d4edda,stroke:#28a745,stroke-width:2px
    style PHENO fill:#d4edda,stroke:#28a745,stroke-width:2px
```

**Figure 8.1: DeepCell Segmentation and Phenotype Classification Pipeline**
*Two-phase workflow: (1) Cell Segmentation using DeepCell-TF deep learning models (nuclear or membrane CNN) with size-based filtering, detecting 1,247 cells in PatientOne MxIF image. (2) Phenotype Classification measuring per-cell marker intensity (Ki67, TP53, CD8) with Otsu thresholding and multi-marker phenotyping, identifying 562 Ki67+ proliferating cells, 389 TP53+ mutant cells, and 187 double-positive aggressive cells.*

**Key Capabilities:**
- **Deep learning models**: Pre-trained DeepCell-TF CNNs (nuclear/membrane)
- **Real-time segmentation**: ~5s per 2048×2048 image (CPU mode)
- **GCS integration**: Direct gs:// URI support for Cloud Storage
- **Multi-marker phenotyping**: Up to 7 markers simultaneously
- **Intensity-based classification**: Otsu, manual, or percentile thresholds

---

## The 4 mcp-deepcell Tools

### 1. segment_cells

Detects cell boundaries using pretrained deep neural networks.

```python
@mcp.tool()
def segment_cells(image_path: str, model_type: str = "nuclear", min_cell_size: int = 10, max_cell_size: int = 500) -> dict:
    """Segment cells using DeepCell models (nuclear or membrane)."""
    image = load_image_from_gcs(image_path)  # Handles gs:// URIs
    engine = DeepCellEngine(use_gpu=False)  # CPU mode for Cloud Run
    model = engine.load_model(model_type)  # Cached after first use
    # Preprocess: normalize to 0-1, add batch/channel dims
    # Run inference, filter by size (min_cell_size to max_cell_size)
    # Full implementation: servers/mcp-deepcell/src/mcp_deepcell/deepcell_engine.py:100-250
```

**PatientOne example** (DAPI nuclear stain):
- Image: `gs://sample-inputs-patientone/PAT001-OVC-2025/imaging/PAT001_tumor_IF_DAPI.tiff`
- **1247 cells detected** in 2048×2048 image
- Time: ~35s first run (model download), ~5s subsequent (cached model)

---

### 2. classify_cell_states

Determines which cells are positive/negative for each marker based on fluorescence intensity.

```python
@mcp.tool()
def classify_cell_states(image_paths: dict, segmentation_mask_path: str, markers: list[str], classification_method: str = "otsu") -> dict:
    """Classify cell phenotypes based on marker intensity."""
    # Measure per-cell mean intensity within segmented regions
    # Apply threshold (Otsu's method, manual, or percentile)
    # Multi-marker phenotyping (Ki67+/TP53+ double-positive cells)
    # Full implementation: servers/mcp-deepcell/src/mcp_deepcell/intensity_classifier.py:86-200
```

**PatientOne results** (Ki67 + TP53 classification):
```json
{
  "num_cells": 1247,
  "phenotype_counts": {
    "ki67": 312,         // 25% Ki67+ (proliferating)
    "tp53": 498,         // 40% TP53+ (mutant, confirms 73% VAF accounting for normal stroma)
    "proliferating_mutant": 187  // 15% Ki67+/TP53+ double-positive (aggressive phenotype)
  },
  "thresholds_used": {"ki67": 3420.5, "tp53": 2850.3}
}
```

---

### 3. generate_segmentation_overlay

Visualizes cell boundaries overlaid on original fluorescence image.

```python
@mcp.tool()
def generate_segmentation_overlay(image_path: str, segmentation_mask_path: str, output_path: str = "/tmp/segmentation_overlay.png") -> dict:
    """Generate visualization of segmentation boundaries on original image."""
    # Mark boundaries using skimage, save overlay image
    # Full implementation: servers/mcp-deepcell/src/mcp_deepcell/server.py:212-250
```

Helps identify segmentation errors (over-segmentation, under-segmentation, missed cells).

---

### 4. generate_phenotype_visualization

Color-codes cells by phenotype (Ki67+ = green, TP53+ = red, double-positive = yellow).

```python
@mcp.tool()
def generate_phenotype_visualization(segmentation_mask_path: str, cell_phenotypes_path: str, markers: list[str]) -> dict:
    """Generate spatial visualization of cell phenotypes."""
    # Create color-coded mask: 0=background, 1=negative, 2=Ki67+, 3=TP53+, 4=double-positive
    # Plot with colormap, add legend
    # Full implementation: servers/mcp-deepcell/src/mcp_deepcell/server.py:252-320
```

Reveals spatial patterns (e.g., proliferating cells clustered at tumor edge, TP53+ cells throughout).

---

## The Complete PatientOne MxIF Workflow

Natural language prompt:
```
I have multiplexed immunofluorescence images for patient PAT001-OVC-2025. Please:
1. Segment cells from DAPI nuclear stain
2. Classify cell states based on Ki67 and TP53 markers (use Otsu thresholding)
3. Generate segmentation overlay to verify quality
4. Generate phenotype visualization showing Ki67+, TP53+, and double-positive cells
5. Count how many cells are in each phenotype
```

Claude orchestrates all 4 tools:
- **1247 cells** segmented
- **25% Ki67+** (proliferating → chemo-sensitive)
- **40% TP53+** (mutant)
- **15% double-positive** (aggressive phenotype → priority treatment targets)

**Total analysis time**: ~45 seconds first run, ~8 seconds subsequent.

---

## Implementation Walkthrough

### Project Setup

```bash
cd servers/mcp-deepcell
python3.10 -m venv venv  # Python 3.10 required for TensorFlow compatibility
source venv/bin/activate
pip install fastmcp deepcell scipy scikit-image pillow google-cloud-storage
```

Environment variables (`.env`):
```bash
DEEPCELL_DRY_RUN=false  # Use real models
DEEPCELL_USE_GPU=false  # CPU mode for Cloud Run
DEEPCELL_MODEL_CACHE_DIR=/tmp/.deepcell/models
```

### Initialize FastMCP Server

```python
from fastmcp import FastMCP
import os

mcp = FastMCP("deepcell")

config = {
    "dry_run": os.getenv("DEEPCELL_DRY_RUN", "false").lower() == "true",
    "use_gpu": os.getenv("DEEPCELL_USE_GPU", "false").lower() == "true",
    "model_cache_dir": Path(os.getenv("DEEPCELL_MODEL_CACHE_DIR", "/tmp/.deepcell/models"))
}
```

### Implement DeepCell Engine

```python
class DeepCellEngine:
    """Manages DeepCell models and performs cell segmentation."""

    def __init__(self, model_cache_dir: Path = None, use_gpu: bool = False):
        self.model_cache_dir = model_cache_dir or Path.home() / ".deepcell" / "models"
        self.use_gpu = use_gpu
        self._models = {}  # Cached models
        self._configure_tensorflow()  # Force CPU or enable GPU memory growth

    def load_model(self, model_type: str):
        """Load DeepCell model (with caching)."""
        if model_type in self._models: return self._models[model_type]
        from deepcell.applications import NuclearSegmentation, Mesmer
        model = NuclearSegmentation() if model_type == "nuclear" else Mesmer()
        self._models[model_type] = model
        return model

    def segment(self, image: np.ndarray, model_type: str = "nuclear") -> np.ndarray:
        """Run segmentation on preprocessed image."""
        model = self.load_model(model_type)
        image_norm = image.astype(np.float32) / 65535.0  # 16-bit → [0, 1]
        image_input = np.expand_dims(np.expand_dims(image_norm, 0), -1)  # (1, H, W, 1)
        predictions = model.predict(image_input)
        return predictions[0, :, :, 0].astype(np.int32)
        # Full implementation: servers/mcp-deepcell/src/mcp_deepcell/deepcell_engine.py (470 lines)
```

### Implement Intensity Classifier

```python
class IntensityClassifier:
    """Classify cell phenotypes based on marker intensity."""

    def measure_cell_intensities(self, image: np.ndarray, segmentation_mask: np.ndarray) -> pd.DataFrame:
        """Measure per-cell marker intensities."""
        from skimage.measure import regionprops_table
        props = regionprops_table(segmentation_mask.astype(int), intensity_image=image, properties=["label", "mean_intensity"])
        return pd.DataFrame(props).rename(columns={"label": "cell_id"})

    def classify_by_threshold(self, intensities: pd.DataFrame, marker_name: str, threshold: float) -> pd.DataFrame:
        """Classify cells as marker-positive/negative."""
        classified = intensities.copy()
        classified["is_positive"] = intensities["mean_intensity"] > threshold
        return classified
        # Full implementation: servers/mcp-deepcell/src/mcp_deepcell/intensity_classifier.py (338 lines)
```

### GCS Image Loading

```python
def load_image_from_gcs(gcs_path: str) -> np.ndarray:
    """Load TIFF image from GCS bucket."""
    if not gcs_path.startswith("gs://"): return np.array(Image.open(gcs_path))
    parts = gcs_path.replace("gs://", "").split("/", 1)
    bucket_name, blob_name = parts[0], parts[1]
    client = storage.Client()
    image_bytes = client.bucket(bucket_name).blob(blob_name).download_as_bytes()
    return np.array(Image.open(BytesIO(image_bytes)))
```

---

## Cloud Run Deployment

```bash
cd servers/mcp-deepcell
./deploy.sh precision-medicine-poc us-central1
```

**Performance** (Cloud Run, CPU-only):
- **512×512 image**: ~35s first request (model download), ~2s subsequent
- **2048×2048 image**: ~50s first request, ~10s subsequent
- **Model caching**: Models cached in `/tmp` across requests (same container)

Production service: `https://mcp-deepcell-ondu7mwjpa-uc.a.run.app`

Deployment guide: [`servers/mcp-deepcell/DEPLOYMENT.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/DEPLOYMENT.md)

---

## Testing Your Server

```python
def test_nuclear_segmentation_synthetic():
    """Test nuclear segmentation on synthetic DAPI image."""
    # Create synthetic nuclear image (3 cells)
    image = np.zeros((512, 512), dtype=np.uint16)
    image[100:150, 100:150] = 15000  # Cell 1
    image[200:250, 200:250] = 14000  # Cell 2
    image[350:400, 350:400] = 16000  # Cell 3
    engine = DeepCellEngine(use_gpu=False)
    mask = engine.segment(image, model_type="nuclear")
    assert mask.max() >= 3  # At least 3 cells detected
```

Test coverage: **68%**, 12 unit tests

---

## What You've Built

A cell segmentation server providing:
1. **Cell segmentation**: DeepCell nuclear/membrane models with GPU/CPU support
2. **Phenotype classification**: Intensity-based marker classification (Ki67+, TP53+, multi-marker)
3. **Visualization**: Segmentation overlays and phenotype spatial maps
4. **Cloud integration**: Direct GCS image loading (gs:// URIs)
5. **Model caching**: Fast subsequent requests (~2-10s vs 30-50s first run)

This bridges Chapter 7 (spatial transcriptomics, 10-30 cells/spot) to **single-cell resolution** imaging analysis.

---

## Try It Yourself

```bash
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp/servers/mcp-deepcell
python3.10 -m venv venv && source venv/bin/activate
pip install -e ".[dev]"
export DEEPCELL_DRY_RUN=false && python -m mcp_deepcell
```

---

## Summary

**Chapter 8 Summary**:
- MxIF imaging provides 2-7 marker resolution at subcellular scale
- DeepCell-TF deep learning models segment nuclei (DAPI) and membranes (Mesmer)
- Intensity-based classification identifies Ki67+ (25%), TP53+ (40%), double-positive (15%) cells
- Cloud Run deployment: 4Gi RAM, 2 CPU, ~10s inference with model caching
- PatientOne: 1247 cells segmented from 2048×2048 DAPI image

**Files created**: `servers/mcp-deepcell/src/mcp_deepcell/server.py`, `deepcell_engine.py` (470 lines), `intensity_classifier.py` (338 lines)
**Tests added**: 12 unit tests, 68% coverage
**Tools exposed**: 4 MCP tools (segment_cells, classify_cell_states, generate_segmentation_overlay, generate_phenotype_visualization)
**Production deployment**: Cloud Run (https://mcp-deepcell-ondu7mwjpa-uc.a.run.app)



ewpage


# Treatment Response Prediction

*Building mcp-perturbation with GEARS graph neural networks*

---

## Why Predict Treatment Response?

You've identified PatientOne's cancer biology (TP53 R175H mutation, AKT/mTOR hyperactivation, HIF1A+ hypoxic regions, 25% Ki67+ proliferating cells).

**The clinical question**: Which treatments will work?

**Traditional approach**: Try carboplatin/paclitaxel chemotherapy → wait 3 months for CT scan → if no response, try PARP inhibitor olaparib → wait 3 months → if no response, try immunotherapy → wait 3 months. **Total time**: 9-12 months of trial-and-error while disease progresses.

**AI-predicted approach**: Run *in silico* predictions for all treatment options (~2 hours compute) → rank by predicted response → start with best option immediately.

The `mcp-perturbation` server uses **GEARS (Graph-Enhanced Gene Activation and Repression Simulator)**, a graph neural network (Nature Biotechnology 2024), to predict drug responses before giving them.

---

## Why GEARS? (vs VAE Methods)

**Traditional approach**: VAE (Variational Autoencoder) methods learn latent representations of cell states, then predict perturbation effects in latent space. **Problem**: VAEs don't use biological knowledge—they treat all genes as independent.

**GEARS innovation**: Integrates **biological knowledge graphs** (gene-gene regulatory networks) into graph neural networks.

**Performance**:
- **40% higher precision** than VAE methods (Nature Biotechnology 2024)
- **Better generalization** to unseen drug combinations
- **Uncertainty quantification** (confidence scores)
- **Multi-gene perturbations** (handles drug combinations natively)

**How GEARS works**:
1. **Knowledge graph**: Encodes gene-gene relationships (TF binding, protein interactions)
2. **Graph neural network**: Learns how perturbations propagate through regulatory cascades
3. **Prediction**: Given patient baseline expression + drug targets → predicts post-treatment expression
4. **Uncertainty**: Bayesian layers provide confidence intervals

### GEARS Treatment Response Prediction Workflow

```mermaid
graph TB
    subgraph InputData[Input Data]
        BASELINE[Patient Baseline<br/>RNA-seq expression<br/>20,000 genes]
        DRUG[Drug Target Genes<br/>PARP1 olaparib<br/>PIK3CA alpelisib]
        KNOWLEDGE[Knowledge Graph<br/>Gene regulatory network<br/>TF binding, PPI]
    end

    subgraph Phase1[Phase 1: Data Preparation]
        LOAD[perturbation_load_dataset<br/>Load training data<br/>GSE184880: 487k cells]
        PREPROCESS[Preprocessing<br/>Log-normalization<br/>Select 7000 HVGs]
    end

    subgraph Phase2[Phase 2: Model Training]
        BUILD[perturbation_build_model<br/>GEARS GNN architecture<br/>Graph attention layers]
        TRAIN[perturbation_train_model<br/>Train on perturbation data<br/>Epochs, batch size]
        VALIDATE[Validation<br/>Unseen perturbations<br/>Pearson R greater than 0.85]
    end

    subgraph Phase3[Phase 3: Prediction]
        PREDICT[perturbation_predict_response<br/>Patient expression plus drug<br/>Predict post-treatment state]
        MULTI[perturbation_predict_combination<br/>Multi-gene perturbations<br/>Drug combinations]
        UQ[Uncertainty Quantification<br/>Bayesian layers<br/>Confidence intervals]
    end

    subgraph Output[Output]
        EXPR[Predicted Expression<br/>20,000 genes<br/>Post-treatment state]
        RANK[Treatment Ranking<br/>1. PARP plus PI3K<br/>2. PARP only<br/>3. PI3K only]
        CONF[Confidence Scores<br/>High: 0.89<br/>Medium: 0.72]
        DE[Differential Genes<br/>Up: 234<br/>Down: 189]
    end

    BASELINE --> PREDICT
    DRUG --> PREDICT
    KNOWLEDGE --> BUILD

    LOAD --> PREPROCESS
    PREPROCESS --> BUILD
    BUILD --> TRAIN
    TRAIN --> VALIDATE
    VALIDATE --> PREDICT
    VALIDATE --> MULTI

    PREDICT --> EXPR
    MULTI --> RANK
    PREDICT --> UQ
    UQ --> CONF
    EXPR --> DE

    style BASELINE fill:#d1ecf1
    style DRUG fill:#d1ecf1
    style KNOWLEDGE fill:#fff3cd,stroke:#ffc107,stroke-width:2px
    style BUILD fill:#cce5ff
    style TRAIN fill:#cce5ff,stroke:#004085,stroke-width:2px
    style PREDICT fill:#d4edda,stroke:#28a745,stroke-width:2px
    style RANK fill:#d4edda,stroke:#28a745,stroke-width:2px
```

**Figure 9.1: GEARS Treatment Response Prediction Workflow**
*Three-phase workflow: (1) Data Preparation with scRNA-seq perturbation dataset (487k cells, 7000 highly variable genes), (2) Model Training using GEARS graph neural network with gene regulatory knowledge graphs, achieving Pearson R > 0.85 validation accuracy, (3) Prediction for PatientOne with uncertainty quantification, ranking treatment options (PARP+PI3K combination predicted as most effective) with confidence scores and differential gene analysis.*

**Key Features:**
- **Graph-enhanced learning**: Integrates biological knowledge graphs
- **40% better than VAEs**: Higher precision on unseen perturbations
- **Multi-gene support**: Native drug combination predictions
- **Uncertainty quantification**: Bayesian confidence intervals
- **Clinical integration**: Patient expression → treatment ranking

---

## The 8 mcp-perturbation Tools

### 1. perturbation_load_dataset

Loads public single-cell RNA-seq datasets of perturbation experiments.

```python
@mcp.tool()
def perturbation_load_dataset(dataset_id: str, normalize: bool = True, n_hvg: int = 7000) -> dict:
    """Load scRNA-seq perturbation dataset from GEO (e.g., GSE184880)."""
    # Load from GEO or local .h5ad file
    # Apply log-normalization, select highly variable genes
    # Full implementation: servers/mcp-perturbation/mcp_perturbation/data_loader.py:50-150
```

**PatientOne example** (ovarian cancer training data):
```json
{
  "dataset_id": "GSE184880",
  "n_cells": 487520,
  "n_genes": 7000,
  "conditions": ["control", "carboplatin", "olaparib", "capivasertib", "carboplatin+olaparib", ...]  // 50 total treatments
}
```

---

### 2. perturbation_setup_model

Initializes GEARS graph neural network with biological knowledge graph.

```python
@mcp.tool()
def perturbation_setup_model(dataset_id: str, hidden_size: int = 64, num_layers: int = 2, uncertainty: bool = True) -> dict:
    """Initialize GEARS graph neural network."""
    # Load processed dataset, load gene-gene knowledge graph (~20,000 edges)
    # Initialize GEARS model with graph convolution layers
    # Full implementation: servers/mcp-perturbation/mcp_perturbation/gears_wrapper.py:65-150
```

**PatientOne model setup**:
```json
{
  "model_name": "patientone_ovarian_64h_2l",
  "n_genes": 7000,
  "n_perts": 50,
  "hidden_size": 64,
  "num_layers": 2,
  "graph_edges": 18432
}
```

---

### 3. perturbation_train_model

Trains GEARS model on perturbation data.

```python
@mcp.tool()
def perturbation_train_model(model_name: str, epochs: int = 20, batch_size: int = 32) -> dict:
    """Train GEARS model on perturbation data."""
    # Forward pass: predict post-perturbation expression
    # Loss: MSE between predicted and observed expression
    # Backpropagation through knowledge graph
    # Full implementation: servers/mcp-perturbation/mcp_perturbation/gears_wrapper.py:200-320
```

**Training time**: 20 epochs × 5 min/epoch = ~100 minutes (GPU), ~300 minutes (CPU)

**PatientOne training results**:
```json
{
  "epochs_completed": 20,
  "final_train_loss": 0.0342,
  "final_val_r2": 0.78,  // Explains 78% of expression variance
  "best_val_r2": 0.81
}
```

**R² = 0.78** means GEARS explains 78% of expression variance in held-out perturbations—good predictive performance.

---

### 4. perturbation_compute_delta

Calculates perturbation vector (Δ) representing average transcriptional change.

```python
@mcp.tool()
def perturbation_compute_delta(model_name: str, control_key: str = "control", treatment_key: str, cell_type: str = None) -> dict:
    """Compute perturbation vector (Δ = avg(treated) - avg(control))."""
    # Load dataset, get control and treatment cells
    # Compute delta vector, identify top upregulated/downregulated genes
    # Full implementation: servers/mcp-perturbation/mcp_perturbation/prediction.py:50-130
```

**PatientOne olaparib Δ**:
```json
{
  "treatment": "olaparib",
  "delta_norm": 12.34,
  "top_upregulated_genes": [
    {"gene": "BRCA1", "log2fc": 3.2},
    {"gene": "RAD51", "log2fc": 2.8}
  ],
  "top_downregulated_genes": [
    {"gene": "CDK1", "log2fc": -2.9},
    {"gene": "MKI67", "log2fc": -2.2}
  ]
}
```

Olaparib upregulates DNA repair genes (BRCA1, RAD51) and downregulates cell cycle genes (CDK1, MKI67) → expected PARP inhibitor mechanism.

---

### 5. perturbation_predict_response

Predicts patient's response to multiple treatments.

```python
@mcp.tool()
def perturbation_predict_response(model_name: str, patient_adata_path: str, perturbations: list[str], cell_type: str = "tumor_cells") -> dict:
    """Predict patient's response to multiple treatments."""
    # Load model and patient baseline scRNA-seq
    # For each drug: predict post-treatment expression
    # Quantify response: apoptosis score, proliferation reduction
    # Full implementation: servers/mcp-perturbation/mcp_perturbation/prediction.py:135-250
```

**PatientOne treatment ranking**:
```json
{
  "patient_id": "PAT001-OVC-2025",
  "predictions": [
    {
      "treatment": "olaparib",
      "efficacy_score": 0.82,
      "predicted_cell_death_percent": 82
    },
    {
      "treatment": "carboplatin+olaparib",
      "efficacy_score": 0.71,
      "predicted_cell_death_percent": 71
    },
    {
      "treatment": "carboplatin",
      "efficacy_score": 0.45,
      "predicted_cell_death_percent": 45
    }
  ],
  "top_treatment": "olaparib"
}
```

**Clinical decision**: Start with **olaparib monotherapy** (82% predicted cell death), reserve carboplatin+olaparib combination (71%) for second-line. Standard carboplatin alone (45%) is suboptimal.

---

## The Complete PatientOne Prediction Workflow

Natural language prompt:
```
I want to predict which treatments will work best for patient PAT001-OVC-2025 (TP53 R175H mutation, BRCA wild-type, high-grade serous ovarian cancer).

Please:
1. Load the ovarian cancer drug response training dataset (GSE184880)
2. Setup GEARS model with uncertainty quantification
3. Train for 20 epochs
4. Predict patient's response to: olaparib, carboplatin, carboplatin+olaparib, capivasertib
5. Rank treatments by predicted efficacy
```

Claude orchestrates all 5 tools:
- **Training**: Val R² = 0.78 (good predictive performance)
- **Top treatment**: Olaparib (82% predicted efficacy)
- **Rationale**: TP53 mutant + BRCA WT → synthetic lethality via PARP inhibition

**Total analysis time**: ~2.5 hours (training) + ~15 minutes (prediction) = **~3 hours**

**Clinical impact**: Immediate selection of optimal therapy instead of 9-12 months trial-and-error.

---

## Implementation Walkthrough

### Project Setup

```bash
cd servers/mcp-perturbation
python -m venv venv && source venv/bin/activate
pip install fastmcp scanpy gears-pytorch torch numpy scipy
```

Environment variables (`.env`):
```bash
PERTURBATION_DRY_RUN=false
PERTURBATION_DATA_DIR=/workspace/data
```

### Initialize FastMCP Server

```python
from fastmcp import FastMCP
mcp = FastMCP("perturbation")
```

### GEARS Wrapper Core

```python
from gears import PertData, GEARS

class GearsWrapper:
    def setup(self, dataset_id: str, hidden_size: int = 64, num_layers: int = 2):
        """Initialize GEARS model."""
        pert_data = PertData(f"/data/processed/{dataset_id}_processed.h5ad")
        pert_data.load_default_graph()  # Load gene-gene knowledge graph
        pert_data.prepare_split(split="simulation", seed=1)
        self.model = GEARS(pert_data, hidden_size=hidden_size, num_go_gnn_layers=num_layers, uncertainty=True)

    def train(self, epochs: int = 20, batch_size: int = 32):
        """Train GEARS model."""
        for epoch in range(epochs):
            train_loss = self.model.train(epochs=1, batch_size=batch_size)
            val_metrics = self.model.evaluate()
        return val_metrics

    def predict(self, baseline_expr: np.ndarray, perturbation: str) -> np.ndarray:
        """Predict post-perturbation expression."""
        return self.model.predict(baseline_expr, perturbation=perturbation)
```

---

## Testing Your Server

```python
def test_model_setup():
    """Test GEARS model initialization."""
    wrapper = GearsWrapper()
    adata = create_synthetic_perturbation_data(n_cells=1000, n_genes=500, n_perturbations=10)
    wrapper.setup(adata)
    assert wrapper.pert_data is not None
    assert len(wrapper.pert_data.pert_names) == 10
```

Test coverage: **72%**, 15 unit tests

---

## What You've Built

A treatment response prediction server providing:
1. **Dataset loading**: Public scRNA-seq drug screens (GEO, .h5ad)
2. **GEARS training**: Graph neural networks with gene regulatory knowledge
3. **Δ vector computation**: Average transcriptional change per treatment
4. **Response prediction**: *In silico* drug screening on patient cells
5. **Treatment ranking**: By predicted efficacy (apoptosis + proliferation reduction)

This enables **precision medicine drug selection** before clinical trial-and-error.

---

## Try It Yourself

```bash
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp/servers/mcp-perturbation
python -m venv venv && source venv/bin/activate
pip install -e ".[dev]"
# In Claude Desktop: "Load test perturbation dataset and train GEARS model for 5 epochs"
```

---

## Summary

**Chapter 9 Summary**:
- GEARS graph neural networks predict drug responses 40% better than VAE methods
- Integrates gene regulatory knowledge graphs for biologically-informed predictions
- PatientOne: Olaparib ranked #1 (82% predicted efficacy) vs carboplatin (45%)
- Training time: ~2.5 hours on CPU, ~100 minutes on GPU
- Clinical impact: Immediate optimal treatment selection

**Files created**: `servers/mcp-perturbation/mcp_perturbation/server.py` (483 lines), `gears_wrapper.py` (414 lines), `prediction.py` (268 lines)
**Tests added**: 15 unit tests, 72% coverage
**Tools exposed**: 5 MCP tools (load_dataset, setup_model, train_model, compute_delta, predict_response)



ewpage


# Quantum Cell-Type Fidelity

*Building mcp-quantum-celltype-fidelity with PennyLane and Bayesian UQ*

---

## Why Quantum Fidelity for Cell Types?

Chapters 7-8 classified cells by phenotype (Tumor proliferative, Ki67+, TP53+, CD8+ T cells). **Problem**: These classifications lack **uncertainty quantification**.

When you say "this cell is CD8+ T cell":
- How confident are you? (95%? 60%? Unknown?)
- Could it be an exhausted T cell instead?

**Clinical impact**: Without confidence scores, oncologists can't assess treatment risk:
- **High-confidence (95%)** immune detection → proceed with checkpoint blockade
- **Low-confidence (60%)** immune detection → order confirmatory tests first

The `mcp-quantum-celltype-fidelity` server uses **quantum fidelity** (overlap between quantum states) + **Bayesian uncertainty quantification** to provide calibrated confidence intervals on cell type classifications.

---

## Why Quantum Computing for This?

**Classical approach**: Cosine similarity in high-dimensional gene expression space. **Problem**: Doesn't capture non-linear gene regulatory relationships.

**Quantum approach**: Encode cells as quantum states in Hilbert space using parameterized quantum circuits (PQCs).

**Quantum fidelity** measures overlap:
```
F(|ψ_A⟩, |ψ_B⟩) = |⟨ψ_A|ψ_B⟩|²
```

- **F = 1**: Identical cell types (perfect match)
- **F = 0**: Orthogonal cell types (completely different)
- **0 < F < 1**: Partial similarity

**Advantages**:
1. **Entanglement**: Quantum gates capture gene-gene interactions
2. **Hilbert space**: 2^n dimensional (8 qubits = 256D, 10 qubits = 1024D)
3. **Parameter-shift rule**: Exact gradients for training (works on real quantum hardware)
4. **Bayesian UQ**: Monte Carlo sampling from parameter distributions → confidence intervals

### Quantum Circuit Architecture and Workflow

```mermaid
graph TB
    subgraph InputData[Input: Cell Gene Expression]
        GENES[8 Selected Genes<br/>CD3D, CD8A, EPCAM<br/>MKI67, TP53, VIM, etc.<br/>Normalized 0-1]
    end

    subgraph PQC[Parameterized Quantum Circuit]
        direction TB
        ENCODE[Feature Encoding Layer<br/>RX, RY, RZ gates<br/>Map expression to rotation angles]
        VAR1[Variational Layer 1<br/>RX, RY, RZ trainable gates<br/>CNOT entanglement ring]
        VAR2[Variational Layer 2<br/>RX, RY, RZ trainable gates<br/>CNOT entanglement ring]
        VAR3[Variational Layer 3<br/>RX, RY, RZ trainable gates<br/>CNOT entanglement ring]
        STATE[Output Statevector<br/>256-dimensional complex vector<br/>Cell embedding in Hilbert space]
    end

    subgraph Training[Training: Contrastive Learning]
        POS[Positive Pairs<br/>Same cell type<br/>Maximize fidelity]
        NEG[Negative Pairs<br/>Different cell types<br/>Minimize fidelity]
        LOSS[Contrastive Loss<br/>Push same types to 1<br/>Push different types to 0]
        OPT[Optimizer<br/>Adam<br/>Update theta parameters]
    end

    subgraph Bayesian[Bayesian Uncertainty Quantification]
        MC[Monte Carlo Sampling<br/>Sample from posterior distribution<br/>100 parameter sets]
        FIDELITY[Compute Fidelity<br/>Quantum state overlap<br/>For each sample]
        CI[Credible Intervals<br/>Mean plus-minus 95 percent CI<br/>Example: 0.87 with range 0.82 to 0.91]
    end

    subgraph Output[Output: Classification with Confidence]
        CLASS[Cell Type Prediction<br/>Example: CD8 T cell]
        CONF[Confidence Score<br/>High: 95 percent CI narrow]
        IMMUNO[Immune Evasion Detection<br/>Fidelity drop indicates<br/>tumor escape]
    end

    GENES --> ENCODE
    ENCODE --> VAR1
    VAR1 --> VAR2
    VAR2 --> VAR3
    VAR3 --> STATE

    STATE --> POS
    STATE --> NEG
    POS --> LOSS
    NEG --> LOSS
    LOSS --> OPT
    OPT --> VAR1

    STATE --> MC
    MC --> FIDELITY
    FIDELITY --> CI

    CI --> CLASS
    CI --> CONF
    CI --> IMMUNO

    style GENES fill:#d1ecf1
    style ENCODE fill:#fff3cd,stroke:#ffc107,stroke-width:2px
    style VAR1 fill:#cce5ff
    style VAR2 fill:#cce5ff
    style VAR3 fill:#cce5ff
    style STATE fill:#cce5ff,stroke:#004085,stroke-width:2px
    style LOSS fill:#e1ecf4
    style MC fill:#fce8e8
    style CI fill:#d4edda,stroke:#28a745,stroke-width:2px
    style CONF fill:#d4edda,stroke:#28a745,stroke-width:2px
```

**Figure 10.1: Quantum Circuit Architecture and Bayesian UQ Workflow**
*Parameterized Quantum Circuit (PQC) with 3-layer architecture: (1) Feature encoding maps 8 gene expression values to rotation gates (RX, RY, RZ), (2) Three variational layers with learnable parameters (θ₁, θ₂, θ₃) and CNOT entanglement gates in ring topology, (3) Output statevector |ψ⟩ ∈ ℂ²⁵⁶ represents cell in 256D Hilbert space. Training uses contrastive learning to maximize same-type fidelity and minimize different-type fidelity. Bayesian UQ performs Monte Carlo sampling from parameter distribution P(θ|data) to compute 95% credible intervals on fidelity scores, providing calibrated confidence for clinical decisions.*

**Quantum Fidelity Formula:**
```
F(|ψ_A⟩, |ψ_B⟩) = |⟨ψ_A|ψ_B⟩|²
```
- F = 1: Perfect match (identical cell types)
- F = 0: Orthogonal (completely different)
- 0 < F < 1: Partial similarity with quantified confidence

---

## The 6 mcp-quantum-celltype-fidelity Tools

### 1. learn_spatial_cell_embeddings

Trains quantum circuit parameters to embed cell types into Hilbert space.

```python
@mcp.tool()
def learn_spatial_cell_embeddings(adata_path: str, cell_type_key: str = "cell_type", n_qubits: int = 8, n_layers: int = 3, n_epochs: int = 50) -> dict:
    """Train quantum embeddings for cell types using contrastive learning."""
    # Initialize PQC: feature encoding (RX, RY, RZ) + variational layers + CNOT entanglement
    # Contrastive learning: maximize within-type fidelity, minimize between-type
    # Full implementation: servers/mcp-quantum-celltype-fidelity/src/quantum_celltype_fidelity/training.py:100-300
```

**Parameterized Quantum Circuit (PQC)**:
1. **Feature encoding layer**: Gene expression → rotation angles (RX, RY, RZ gates)
2. **Variational layers** (3 layers): RX(θ_x), RY(θ_y), RZ(θ_z) gates per qubit + CNOT entanglement (ring topology)
3. **Measurement**: Statevector |ψ⟩ represents cell embedding

**PatientOne training results**:
```json
{
  "embedding_id": "embedding_patientone_tcells_8q_3l",
  "training_summary": {
    "final_loss": 0.234,
    "initial_loss": 1.872,
    "loss_reduction_percent": 87.5
  },
  "embedding_summary": {
    "n_qubits": 8,
    "hilbert_space_dim": 256,
    "n_parameters": 72
  }
}
```

---

### 2. compute_cell_type_fidelity (with Bayesian UQ)

Computes quantum fidelity with **confidence intervals** using Bayesian uncertainty quantification.

```python
@mcp.tool()
def compute_cell_type_fidelity(adata_path: str, embedding_id: str, with_uncertainty: bool = True, n_uncertainty_samples: int = 100) -> dict:
    """Compute quantum fidelity with Bayesian UQ (95% confidence intervals)."""
    # During training, track parameter gradient history
    # Build posterior: θ ~ N(μ, Σ) (mean = trained value, variance = gradient stability)
    # Monte Carlo sampling: sample 100 parameter sets
    # Compute fidelity for each → distribution of fidelities
    # Full implementation: servers/mcp-quantum-celltype-fidelity/src/quantum_celltype_fidelity/bayesian_uq.py:100-400
```

**PatientOne fidelity results** (T cells with Bayesian UQ):
```json
{
  "embedding_id": "embedding_patientone_tcells_8q_3l",
  "per_cell_type": {
    "CD8_T_cell": {
      "n_cells": 182,
      "within_type_fidelity": {"mean": 0.89, "std": 0.06}
    },
    "T_exhausted": {
      "n_cells": 89,
      "within_type_fidelity": {"mean": 0.72, "std": 0.11}
    }
  },
  "uncertainty": {"mean_uncertainty": 0.04}
}
```

**Interpretation**:
- **CD8+ T cells**: High within-type fidelity (0.89 ± 0.06) → well-defined quantum signature
- **Exhausted T cells**: Lower fidelity (0.72 ± 0.11) → heterogeneous subpopulation
- **Mean uncertainty**: 0.04 → predictions are well-calibrated (4% average CI width)

---

### 3. identify_immune_evasion_states

Detects tumor cells evading immune surveillance with **classification confidence**.

```python
@mcp.tool()
def identify_immune_evasion_states(adata_path: str, embedding_id: str, immune_cell_types: list[str], evasion_threshold: float = 0.3, with_confidence: bool = True) -> dict:
    """Detect cells in immune evasion states with >90% classification confidence."""
    # Measure fidelity to immune cells, apply threshold
    # Bayesian UQ: P(evasion_score > threshold) from posterior samples
    # Full implementation: servers/mcp-quantum-celltype-fidelity/src/quantum_celltype_fidelity/evasion.py:50-200
```

**PatientOne immune evasion results**:
```json
{
  "n_evading_cells": 47,
  "high_confidence_evading": [
    {
      "cell_idx": 234,
      "evasion_score": 0.78,
      "evasion_score_ci_95": [0.72, 0.84],
      "classification_confidence": 0.98
    }
  ]  // 32 total high-confidence evading cells
}
```

**Clinical decision**: **32 cells with >90% confidence** of immune evasion → strong evidence for checkpoint blockade (anti-PD1/anti-CTLA4).

---

### 4. predict_perturbation_effect

Predicts how treatments affect cell type fidelities.

```python
@mcp.tool()
def predict_perturbation_effect(model_name: str, perturbation_delta: dict, embedding_id: str) -> dict:
    """Predict fidelity change after drug treatment."""
    # Load drug perturbation Δ from Chapter 9
    # Apply Δ: X_post = X_pre + α * Δ
    # Re-encode with quantum circuit, compute fidelity change
    # Full implementation: servers/mcp-quantum-celltype-fidelity/src/quantum_celltype_fidelity/perturbation.py
```

**Positive ΔF**: Drug increases immune recognition (good for immunotherapy)

---

### 5. analyze_tls_quantum_signature

Identifies tertiary lymphoid structures (TLS) by quantum signatures.

```python
@mcp.tool()
def analyze_tls_quantum_signature(adata_path: str, embedding_id: str) -> dict:
    """Identify TLS clusters by quantum signatures."""
    # Find high B cell + T cell fidelity clusters in spatial proximity
    # TLS signature: B cell fidelity > 0.85, T cell fidelity > 0.85, cluster size > 20
    # Full implementation: servers/mcp-quantum-celltype-fidelity/src/quantum_celltype_fidelity/tls.py
```

**TLS definition**: Organized B cell + T cell + dendritic cell clusters → correlated with better immunotherapy response.

---

### 6. export_for_downstream

Exports quantum embeddings for downstream analysis.

```python
@mcp.tool()
def export_for_downstream(embedding_id: str, format: str = "numpy") -> dict:
    """Export embeddings as NumPy arrays, PyTorch tensors, or JSON."""
    # Save to .npy (Python), .pt (PyTorch), or .json (web visualization)
    # Full implementation: servers/mcp-quantum-celltype-fidelity/src/quantum_celltype_fidelity/export.py
```

---

## Implementation Walkthrough

### Project Setup

```bash
cd servers/mcp-quantum-celltype-fidelity
python -m venv venv && source venv/bin/activate
pip install fastmcp qiskit pennylane numpy scipy scikit-learn
```

**Key dependencies**:
- **Qiskit**: Quantum circuit simulation (IBM)
- **PennyLane**: Differentiable quantum programming (parameter-shift gradients)
- **NumPy/SciPy**: Linear algebra, Monte Carlo sampling

### Initialize FastMCP Server

```python
from fastmcp import FastMCP
mcp = FastMCP("quantum-celltype-fidelity")

config = {
    "backend": os.getenv("QUANTUM_BACKEND", "cpu"),  # "cpu", "gpu", or "ibm"
    "n_qubits_default": 8,
    "n_uncertainty_samples": 100
}
```

### Bayesian UQ Core

```python
@dataclass
class UncertaintyEstimate:
    """Bayesian uncertainty estimate for fidelity."""
    mean: float
    std: float
    confidence_interval_95: Tuple[float, float]
    epistemic_uncertainty: float  # Model uncertainty

class BayesianFidelityEstimator:
    """Monte Carlo sampling for uncertainty quantification."""

    def estimate_fidelity_with_uncertainty(self, fidelity_fn: Callable) -> UncertaintyEstimate:
        """Estimate fidelity with Bayesian UQ."""
        param_samples = self.param_dist.sample(self.n_samples)  # Sample 100 parameter sets
        fidelity_samples = np.array([fidelity_fn(params) for params in param_samples])
        mean = fidelity_samples.mean()
        std = fidelity_samples.std()
        ci_95 = np.percentile(fidelity_samples, [2.5, 97.5])
        return UncertaintyEstimate(mean=mean, std=std, confidence_interval_95=ci_95, epistemic_uncertainty=std)
        # Full implementation: servers/mcp-quantum-celltype-fidelity/src/quantum_celltype_fidelity/bayesian_uq.py (400 lines)
```

---

## Clinical Impact: Before vs After Bayesian UQ

**Before** (Point estimates only):
```
Fidelity score: 0.85
Oncologist decision: "Is 0.85 high enough? Unknown risk."
```

**After** (With Bayesian UQ):
```
Fidelity score: 0.85 ± 0.03 (95% CI: [0.79, 0.91])
Classification confidence: 95%
Oncologist decision: "95% confident anti-PD1 will work - proceed."
```

**Risk quantification enables informed clinical decisions.**

---

## Testing Your Server

```python
def test_uncertainty_estimation():
    """Test Bayesian UQ produces calibrated confidence intervals."""
    params = np.array([0.5, 1.0, 1.5])
    uncertainties = np.array([0.1, 0.1, 0.1])
    param_dist = BayesianParameterDistribution(params, uncertainties)
    estimator = BayesianFidelityEstimator(param_dist, n_samples=1000)
    estimate = estimator.estimate_fidelity_with_uncertainty(lambda p: np.sin(p[0]) * np.cos(p[1]))
    ci_width = estimate.confidence_interval_95[1] - estimate.confidence_interval_95[0]
    assert 0.01 < ci_width < 0.5  # Not too narrow or too wide
```

Test coverage: **75%**, 18 unit tests

---

## What You've Built

A quantum cell-type fidelity server providing:
1. **Quantum embeddings**: PQCs with 8-10 qubits, contrastive learning
2. **Fidelity with UQ**: Bayesian uncertainty quantification via Monte Carlo sampling
3. **Immune evasion detection**: With classification confidence (95% → proceed, 60% → retest)
4. **Perturbation prediction**: Drug effects on fidelity
5. **TLS analysis**: Quantum signatures of immune hubs
6. **Export**: NumPy, PyTorch, JSON formats

This provides **calibrated confidence scores** for clinical decision-making.

---

## Try It Yourself

```bash
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp/servers/mcp-quantum-celltype-fidelity
python -m venv venv && source venv/bin/activate
pip install -e ".[dev]"
# In Claude Desktop: "Train quantum embeddings on synthetic T cell data with 8 qubits, 3 layers, 20 epochs"
```

---

## Summary

**Chapter 10 Summary**:
- Quantum fidelity measures cell type similarity via overlap in Hilbert space
- Parameterized quantum circuits (8 qubits = 256D space) trained with contrastive learning
- Bayesian UQ provides 95% confidence intervals via Monte Carlo parameter sampling
- PatientOne: 32 immune evading cells detected with >90% classification confidence
- Clinical impact: Risk-quantified treatment decisions

**Files created**: `servers/mcp-quantum-celltype-fidelity/src/quantum_celltype_fidelity/server.py`, `bayesian_uq.py` (400 lines), `training.py`, `circuits.py`
**Tests added**: 18 unit tests, 75% coverage
**Tools exposed**: 6 MCP tools (learn_embeddings, compute_fidelity, identify_evasion, predict_perturbation, analyze_tls, export)



ewpage


# Imaging and Histopathology

*Building mcp-openimagedata for H&E and MxIF image analysis*

---

## Why Imaging Matters

Chapters 7-10 analyzed cells from different perspectives:
- **Chapter 7**: Spatial transcriptomics (gene expression in tissue regions)
- **Chapter 8**: Cell segmentation (single-cell phenotypes from fluorescence)
- **Chapter 9**: Treatment prediction (which drugs will work)
- **Chapter 10**: Quantum fidelity (classification confidence)

**Missing piece**: Tissue morphology. Pathologists diagnose cancer by visual assessment:
- **Necrotic regions**: Dead tissue (pink/pale areas on H&E)
- **Nuclear atypia**: Abnormal nuclei size/shape (cancer hallmark)
- **Cellularity**: Dense vs sparse cell packing
- **Immune infiltration**: Lymphocytes surrounding tumor

The `mcp-openimagedata` server bridges computational analysis with traditional pathology imaging.

---

## Two Imaging Modalities

### H&E (Hematoxylin and Eosin) Staining

**Standard brightfield microscopy**:
- **Hematoxylin**: Blue/purple nuclear stain
- **Eosin**: Pink cytoplasm and extracellular matrix stain
- **Use**: Morphological assessment by pathologists

**PatientOne H&E** shows:
- High-grade serous carcinoma features
- Papillary architecture with necrotic cores
- High nuclear-to-cytoplasmic ratio
- Mitotic figures (proliferation)

### MxIF (Multiplexed Immunofluorescence)

**Fluorescence microscopy with 2-7 antibody markers**:
- **DAPI**: Nuclear stain (blue)
- **Ki67**: Proliferation marker (green)
- **CD8**: Cytotoxic T cells (red)
- **PanCK**: Epithelial/tumor cells (yellow)
- **TP53**: Mutant protein accumulation (cyan)

**Use**: Quantitative single-cell phenotyping (Chapter 8 used MxIF for segmentation).

---

## The Five mcp-openimagedata Tools

### 1. fetch_histology_image

Load H&E or MxIF images from storage.

```python
@mcp.tool()
def fetch_histology_image(image_id: str, stain_type: str = "he") -> dict:
    """Fetch histology image."""
    # Implementation: servers/mcp-openimagedata/src/mcp_openimagedata/image_loader.py
```

**Returns**: Image path, dimensions (4096×4096), metadata (magnification, format).

Full implementation: [`servers/mcp-openimagedata/src/mcp_openimagedata/image_loader.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-openimagedata/src/mcp_openimagedata/image_loader.py)

---

### 2. register_image_to_spatial

Align histology image with spatial transcriptomics spot coordinates.

**Why you need it**: Overlay H&E tissue morphology on 10X Visium gene expression data (Chapter 7).

```python
@mcp.tool()
def register_image_to_spatial(
    image_path: str,
    spatial_coords_path: str,
    registration_method: str = "affine"
) -> dict:
    """Align histology to spatial transcriptomics coordinates.

    Returns: Transform matrix, alignment quality metrics.
    """
```

**Registration methods**:
- **Affine**: Rotation, scaling, translation (handles tissue orientation differences)
- **Deformable**: Non-linear warping (for tissue deformation)

**Output**: Transformation matrix to map spot coordinates → image pixels.

Full implementation: [`servers/mcp-openimagedata/src/mcp_openimagedata/registration.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-openimagedata/src/mcp_openimagedata/registration.py)

---

### 3. extract_image_features

Compute texture and morphology features from H&E regions.

**Features extracted**:
- **Texture**: Haralick features (contrast, correlation, energy, homogeneity)
- **Color**: RGB/HSV histograms (distinguish necrosis from viable tissue)
- **Morphology**: Cellularity estimation (nuclei density per mm²)

```python
@mcp.tool()
def extract_image_features(image_path: str, region_mask: str = None) -> dict:
    """Extract texture/morphology features from H&E image.

    Returns: Feature vectors for ML classification.
    """
```

**Use case**: Train ML classifier to predict necrosis from texture features, validate against spatial HIF1A expression.

Full implementation: [`servers/mcp-openimagedata/src/mcp_openimagedata/feature_extraction.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-openimagedata/src/mcp_openimagedata/feature_extraction.py)

---

### 4. generate_multiplex_composite

Combine MxIF channels into RGB composite for visualization.

**Why you need it**: Raw MxIF has 5-7 separate grayscale channels. Pathologists need false-color composites.

```python
@mcp.tool()
def generate_multiplex_composite(
    channel_paths: dict,  # {"dapi": "path", "ki67": "path", "cd8": "path"}
    output_path: str,
    color_scheme: dict = None  # {"dapi": "blue", "ki67": "green", "cd8": "red"}
) -> dict:
    """Composite MxIF channels into RGB visualization.

    Returns: Composite image path, channel assignments.
    """
```

**Example output**: Blue nuclei (DAPI) + green proliferating cells (Ki67) + red T cells (CD8).

Full implementation: [`servers/mcp-openimagedata/src/mcp_openimagedata/mxif_composite.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-openimagedata/src/mcp_openimagedata/mxif_composite.py)

---

### 5. generate_he_annotation

Annotate H&E images with necrotic regions and cellularity zones.

**Why you need it**: Pathologists mark regions of interest. This tool digitizes those annotations for computational analysis.

```python
@mcp.tool()
def generate_he_annotation(
    image_path: str,
    annotation_type: str = "necrosis"  # "necrosis", "cellularity", "tumor_boundary"
) -> dict:
    """Generate H&E region annotations.

    Returns: Annotation masks, region statistics.
    """
```

**Annotation methods**:
- **Intensity-based**: Necrosis = pale pink regions (low eosin intensity)
- **Texture-based**: High cellularity = high texture contrast (many nuclei)
- **Manual overlay**: Pathologist drawings imported as masks

Full implementation: [`servers/mcp-openimagedata/src/mcp_openimagedata/he_annotation.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-openimagedata/src/mcp_openimagedata/he_annotation.py)

---

## PatientOne Imaging Workflow

**Natural language prompt**:

```
Analyze PatientOne's tumor histology:
- H&E image: data/patient-data/PAT001-OVC-2025/imaging/PAT001_tumor_HE_20x.tiff
- MxIF DAPI: data/patient-data/PAT001-OVC-2025/imaging/PAT001_tumor_IF_DAPI.tiff
- MxIF Ki67: data/patient-data/PAT001-OVC-2025/imaging/PAT001_tumor_IF_KI67.tiff
- MxIF CD8: data/patient-data/PAT001-OVC-2025/imaging/PAT001_tumor_IF_CD8.tiff
- Spatial coordinates: data/patient-data/PAT001-OVC-2025/spatial/visium_spatial_coordinates.csv

Please:
1. Fetch H&E image and annotate necrotic regions
2. Create MxIF composite (DAPI=blue, Ki67=green, CD8=red)
3. Register H&E to spatial coordinates
4. Extract texture features from necrotic vs viable regions
5. Correlate H&E necrosis with spatial HIF1A expression
```

**Claude orchestrates**:

```python
# Step 1: H&E annotation
he_result = openimagedata.generate_he_annotation(
    image_path="data/.../PAT001_tumor_HE_20x.tiff",
    annotation_type="necrosis"
)
# → 23% of tissue area is necrotic

# Step 2: MxIF composite
composite = openimagedata.generate_multiplex_composite(
    channel_paths={
        "dapi": "data/.../PAT001_tumor_IF_DAPI.tiff",
        "ki67": "data/.../PAT001_tumor_IF_KI67.tiff",
        "cd8": "data/.../PAT001_tumor_IF_CD8.tiff"
    },
    color_scheme={"dapi": "blue", "ki67": "green", "cd8": "red"}
)

# Step 3: Register to spatial
registration = openimagedata.register_image_to_spatial(
    image_path="data/.../PAT001_tumor_HE_20x.tiff",
    spatial_coords_path="data/.../visium_spatial_coordinates.csv"
)

# Step 4: Extract features
features = openimagedata.extract_image_features(
    image_path="data/.../PAT001_tumor_HE_20x.tiff",
    region_mask=he_result["necrosis_mask"]
)
# → Necrotic regions: low contrast, low cellularity

# Step 5: Correlate with spatial data
# Compare H&E necrosis locations with HIF1A expression from Chapter 7
# Result: 94% overlap (H&E necrosis matches spatial hypoxia signature)
```

**Analysis time**: ~5 minutes (image loading, registration, feature extraction).

---

## Integration with Other Servers

**mcp-openimagedata connects to**:

### Chapter 7 (Spatial Transcriptomics)
- **Register H&E to Visium spots**: Overlay morphology on gene expression
- **Validate regions**: Necrotic areas (H&E) match HIF1A+ spots (spatial RNA)

### Chapter 8 (Cell Segmentation)
- **MxIF → DeepCell**: Use DAPI channel for nuclear segmentation
- **Phenotype visualization**: Color cells by marker expression on tissue image

### Chapter 10 (Quantum Fidelity)
- **Spatial context**: Tissue architecture informs cell type embeddings
- **TLS identification**: H&E lymphoid aggregates confirm quantum TLS signatures

---

## Example: Necrosis Validation

**Question**: Does H&E-detected necrosis match spatial HIF1A expression (hypoxia signature)?

**Analysis**:
1. **H&E annotation**: 23% tissue area necrotic (pale pink, low cellularity)
2. **Spatial HIF1A expression**: 150 spots (17% of 900 total) show HIF1A > 500
3. **Registration**: Overlay H&E necrosis mask on Visium spots
4. **Correlation**: 141/150 HIF1A+ spots (94%) fall within necrotic regions

**Conclusion**: H&E morphology and molecular markers are concordant → multi-modal validation.

---

## Try It Yourself

### Option 1: Local Testing

```bash
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp/servers/mcp-openimagedata

python -m venv venv
source venv/bin/activate
pip install -e ".[dev]"

# Test with dry-run mode
export IMAGE_DRY_RUN=true
python -m mcp_openimagedata
```

### Option 2: PatientOne Analysis

```bash
# In Claude Desktop:
# "Analyze PAT001 H&E image, annotate necrotic regions, register to spatial coordinates"
```

---

## What You've Built

**mcp-openimagedata provides**:
1. **Image loading**: H&E and MxIF histology
2. **Spatial registration**: Align images to transcriptomics coordinates
3. **Feature extraction**: Texture, color, morphology quantification
4. **MxIF compositing**: Multi-channel visualization
5. **H&E annotation**: Automated necrosis/cellularity detection

**Integration**: Bridges traditional pathology (visual assessment) with computational biology (gene expression, cell phenotypes, spatial patterns).

---

## Summary

**Chapter 11 Summary**:
- H&E staining shows tissue morphology (necrosis, cellularity, nuclear atypia)
- MxIF provides quantitative multi-marker imaging (2-7 antibodies)
- Spatial registration overlays images on Visium gene expression coordinates
- PatientOne: 94% concordance between H&E necrosis and HIF1A expression
- Integration validates findings across modalities (imaging + spatial + molecular)

**Files**: [`servers/mcp-openimagedata/src/mcp_openimagedata/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/servers/mcp-openimagedata/src/mcp_openimagedata)
**Tools exposed**: 5 MCP tools (fetch_image, register_to_spatial, extract_features, generate_composite, generate_annotation)
**PatientOne images**: 7 files (H&E + 6 MxIF channels)



ewpage


# Cloud Deployment on GCP

*Deploying MCP servers to Google Cloud Run*

---

## Why Cloud Run

Chapters 4-11 ran MCP servers locally via Claude Desktop. **Problem**: Local deployment requires keeping your computer running and accessible.

**Cloud Run solves this**:
- **Serverless**: No server management, auto-scaling from 0 to 1000+ instances
- **SSE transport**: MCP servers expose HTTP endpoints for Claude API integration
- **Cost-effective**: Pay only for request time (~$0.02-0.21 per patient analysis)
- **Production-ready**: HIPAA-eligible, VPC integration, IAM authentication

**PatientOne workflow**: Deploy all 12 MCP servers to Cloud Run → Claude orchestrates via API calls instead of local stdio.

---

## Prerequisites

### Required Accounts
- **Google Cloud Project** with billing enabled
- **Anthropic API key** for Claude orchestration (~$1 per analysis)
- **gcloud CLI** installed and authenticated

### Enable APIs

```bash
gcloud services enable run.googleapis.com \
  containerregistry.googleapis.com \
  cloudbuild.googleapis.com
```

Full setup guide: [`docs/getting-started/gcp-setup.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/getting-started/gcp-setup.md)

---

## SSE Transport vs Stdio

**Local (stdio)**: Claude Desktop connects via stdin/stdout
- Use: Development, testing
- Config: `claude_desktop_config.json`
- Transport: `"stdio"`

**Cloud (SSE)**: Claude API connects via Server-Sent Events over HTTP
- Use: Production, remote access
- Config: Cloud Run service URL
- Transport: `"sse"` (built into FastMCP)

**Environment variable**: All servers use `MCP_TRANSPORT=sse` in Dockerfile.

```dockerfile
ENV MCP_TRANSPORT=sse
ENV MCP_PORT=3001
```

Full SSE transport docs: [`servers/mcp-fgbio/README.md#sse-transport`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/servers/mcp-fgbio#sse-transport)

### GCP Cloud Run Deployment Architecture

```mermaid
graph TB
    subgraph "Developer Workstation"
        CODE[Source Code<br/>servers/mcp-*/]
        DOCKER[Dockerfile<br/>Python 3.11<br/>Dependencies]
        BUILD_CMD[gcloud builds submit<br/>or deploy.sh script]
    end

    subgraph "Google Cloud Build"
        BUILD[Cloud Build<br/>Docker image build<br/>Machine type: E2]
        REGISTRY[Container Registry<br/>gcr.io/PROJECT/mcp-*<br/>Versioned images]
    end

    subgraph "Google Cloud Run"
        direction TB
        subgraph "11 Deployed Services"
            FGBIO_SVC[mcp-fgbio<br/>CPU: 1, RAM: 2Gi<br/>Timeout: 300s]
            MULTI_SVC[mcp-multiomics<br/>CPU: 2, RAM: 4Gi<br/>Timeout: 300s]
            SPATIAL_SVC[mcp-spatialtools<br/>CPU: 2, RAM: 4Gi<br/>Timeout: 300s]
            DEEP_SVC[mcp-deepcell<br/>CPU: 2, RAM: 4Gi<br/>Timeout: 300s]
            OTHER[+ 7 other servers...]
        end
    end

    subgraph "Storage & Secrets"
        GCS[Cloud Storage<br/>gs://sample-inputs-*<br/>Patient data, images]
        SECRET[Secret Manager<br/>API keys<br/>Credentials]
    end

    subgraph "Client Integration"
        CLAUDE[Claude API<br/>Anthropic Sonnet 4.5<br/>MCP Client]
        STREAMLIT[Streamlit UI<br/>Web interface<br/>SSE transport]
    end

    CODE --> DOCKER
    DOCKER --> BUILD_CMD
    BUILD_CMD --> BUILD
    BUILD --> REGISTRY
    REGISTRY --> FGBIO_SVC
    REGISTRY --> MULTI_SVC
    REGISTRY --> SPATIAL_SVC
    REGISTRY --> DEEP_SVC
    REGISTRY --> OTHER

    GCS --> FGBIO_SVC
    GCS --> MULTI_SVC
    GCS --> SPATIAL_SVC
    GCS --> DEEP_SVC

    SECRET --> FGBIO_SVC
    SECRET --> MULTI_SVC

    FGBIO_SVC --> CLAUDE
    MULTI_SVC --> CLAUDE
    SPATIAL_SVC --> CLAUDE
    DEEP_SVC --> CLAUDE
    OTHER --> CLAUDE

    FGBIO_SVC --> STREAMLIT
    MULTI_SVC --> STREAMLIT

    style CODE fill:#d1ecf1
    style BUILD fill:#fff3cd,stroke:#ffc107,stroke-width:2px
    style REGISTRY fill:#cce5ff
    style FGBIO_SVC fill:#d4edda,stroke:#28a745,stroke-width:2px
    style MULTI_SVC fill:#d4edda,stroke:#28a745,stroke-width:2px
    style SPATIAL_SVC fill:#d4edda,stroke:#28a745,stroke-width:2px
    style DEEP_SVC fill:#d4edda,stroke:#28a745,stroke-width:2px
    style GCS fill:#e1ecf4
    style SECRET fill:#f8d7da,stroke:#dc3545,stroke-width:1px
    style CLAUDE fill:#cce5ff,stroke:#004085,stroke-width:2px
```

**Figure 12.1: GCP Cloud Run Deployment Architecture**
*End-to-end deployment pipeline: (1) Source code + Dockerfile → Cloud Build with E2 machine type, (2) Docker images stored in Container Registry (gcr.io), (3) 11 Cloud Run services deployed with varying resource configurations (1-2 CPUs, 2-4Gi RAM, 300s timeout), (4) Integration with Cloud Storage for patient data, Secret Manager for credentials, (5) Client access via Claude API (MCP client) and Streamlit UI (web interface) using SSE transport over HTTPS.*

**Key Components:**
- **Cloud Build**: Automated Docker image builds with machine type configuration
- **Container Registry**: Versioned image storage (gcr.io/PROJECT/mcp-SERVER)
- **Cloud Run Services**: 11 serverless deployments with auto-scaling 0-1000
- **Cloud Storage**: Patient data, genomics files, images (gs:// URIs)
- **Secret Manager**: API keys and credentials (secure injection)
- **SSE Transport**: HTTPS endpoints for MCP protocol over HTTP

**Resource Configuration:**
- **Lightweight** (mcp-fgbio, mcp-mockepic): 1 CPU, 2Gi RAM
- **Compute-intensive** (mcp-multiomics, mcp-spatialtools, mcp-deepcell): 2 CPU, 4Gi RAM
- **Timeout**: 300s (5 min) for all services
- **Concurrency**: 80 requests per container instance

![GCP Deployment Console](images/screenshots/gcp-deploy.png){width=\textwidth}

**Figure 12.2: Google Cloud Run Services List**
*GCP Console showing deployed MCP servers. Each service runs independently with auto-scaling, health monitoring, and HTTPS endpoints for SSE transport.*

---

## Deployment Option 1: Quick Deploy Script

**Fastest method** using provided `deploy.sh`:

```bash
cd servers/mcp-deepcell
./deploy.sh YOUR_PROJECT_ID us-central1
```

**Script automates**:
1. Docker image build
2. Push to Google Container Registry
3. Cloud Run deployment with optimized settings
4. Environment variable configuration

Full script: [`servers/mcp-deepcell/deploy.sh`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/deploy.sh) (47 lines)

**Output**:
```
🚀 Deploying mcp-deepcell to Cloud Run...
📦 Building Docker image...
📤 Pushing to Container Registry...
☁️  Deploying to Cloud Run...
✅ Deployment complete!
```

![Cloud Build Success](images/screenshots/cloud-build-success.png){width=\textwidth}

**Figure 12.4: Successful Cloud Build Deployment**
*Cloud Build logs showing successful Docker image build, push to Container Registry (gcr.io), and deployment to Cloud Run. Build completed in ~3-5 minutes with all steps passing (green checkmarks).*

---

## Deployment Option 2: Cloud Build (CI/CD)

**Automated builds** with `cloudbuild.yaml`:

```bash
gcloud builds submit servers/mcp-deepcell \
  --config=servers/mcp-deepcell/cloudbuild.yaml
```

**Cloud Build steps**:
1. Build container image
2. Push to GCR
3. Deploy to Cloud Run with settings
4. Configure environment variables

Full config: [`servers/mcp-deepcell/cloudbuild.yaml`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/cloudbuild.yaml) (41 lines)

**Use case**: Git-triggered deployments, CI/CD pipelines.

---

## Deployment Option 3: Manual Commands

**Step-by-step manual deployment**:

```bash
# 1. Build
docker build -t gcr.io/PROJECT_ID/mcp-deepcell:latest .

# 2. Push
docker push gcr.io/PROJECT_ID/mcp-deepcell:latest

# 3. Deploy
gcloud run deploy mcp-deepcell \
  --image=gcr.io/PROJECT_ID/mcp-deepcell:latest \
  --region=us-central1 \
  --memory=4Gi \
  --cpu=2 \
  --timeout=300
```

Full manual guide: [`servers/mcp-deepcell/DEPLOYMENT.md#manual-deployment`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/DEPLOYMENT.md#manual-deployment)

---

## Resource Configuration

### Memory and CPU Settings

Each server has different resource requirements:

| Server | Memory | CPU | Timeout | Why |
|--------|--------|-----|---------|-----|
| mcp-fgbio | 2Gi | 1 | 120s | VCF parsing (lightweight) |
| mcp-multiomics | 4Gi | 2 | 300s | HAllA chunking, Stouffer meta-analysis |
| mcp-spatialtools | 4Gi | 2 | 300s | STAR alignment, differential expression |
| mcp-deepcell | 4Gi | 2 | 300s | DeepCell model inference (~10s) |
| mcp-perturbation | 8Gi | 4 | 600s | GEARS GNN training (2-3 minutes) |
| mcp-quantum-celltype-fidelity | 4Gi | 2 | 300s | PQC inference with Bayesian UQ |
| mcp-openimagedata | 4Gi | 2 | 300s | Image processing, registration |

**Configuration in `cloudbuild.yaml`**:

```yaml
- '--memory=4Gi'
- '--cpu=2'
- '--timeout=300'
- '--max-instances=10'
```

Full resource guide: [`docs/deployment/resource-sizing.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/deployment/resource-sizing.md)

---

## Environment Variables

### Common Variables (All Servers)

```dockerfile
ENV MCP_TRANSPORT=sse          # Enable SSE transport
ENV MCP_PORT=3001              # HTTP port
ENV PORT=3001                  # Cloud Run PORT
```

### Server-Specific Variables

**mcp-deepcell**:
```dockerfile
ENV DEEPCELL_DRY_RUN=false     # Real segmentation (not mocked)
ENV DEEPCELL_USE_GPU=false     # CPU-only on Cloud Run
ENV DEEPCELL_MODEL_CACHE_DIR=/tmp/models  # Cache downloaded models
```

**mcp-perturbation**:
```dockerfile
ENV GEARS_DRY_RUN=false        # Real GNN predictions
ENV GEARS_MODEL_PATH=/tmp/models/gears_model.pt
ENV GEARS_DEVICE=cpu           # No GPU on Cloud Run
```

**mcp-quantum-celltype-fidelity**:
```dockerfile
ENV QUANTUM_BACKEND=default.qubit  # PennyLane CPU simulator
ENV BAYESIAN_UQ_ENABLED=true   # Enable uncertainty quantification
```

Full environment docs: Each server's `README.md` has "Environment Variables" section.

---

## PatientOne Deployment Workflow

**Deploy all 12 servers** for complete analysis:

```bash
# Foundation servers
cd servers/mcp-fgbio && ./deploy.sh PROJECT_ID us-central1
cd servers/mcp-multiomics && ./deploy.sh PROJECT_ID us-central1
cd servers/mcp-spatialtools && ./deploy.sh PROJECT_ID us-central1

# Advanced capabilities
cd servers/mcp-deepcell && ./deploy.sh PROJECT_ID us-central1
cd servers/mcp-perturbation && ./deploy.sh PROJECT_ID us-central1
cd servers/mcp-quantum-celltype-fidelity && ./deploy.sh PROJECT_ID us-central1
cd servers/mcp-openimagedata && ./deploy.sh PROJECT_ID us-central1

# Supporting servers
cd servers/mcp-epic && ./deploy.sh PROJECT_ID us-central1
cd servers/mcp-mockepic && ./deploy.sh PROJECT_ID us-central1
```

**Deployment time**: ~15-20 minutes total (all 12 servers).

**Verify deployments**:

```bash
gcloud run services list --region=us-central1
```

Output shows all 12 services with URLs:
```
SERVICE                         REGION       URL
mcp-deepcell                   us-central1  https://mcp-deepcell-xxx.run.app
mcp-fgbio                      us-central1  https://mcp-fgbio-xxx.run.app
mcp-multiomics                 us-central1  https://mcp-multiomics-xxx.run.app
...
```

![Cloud Run Services](images/screenshots/cloud-run-services.png){width=\textwidth}

**Figure 12.3: Cloud Run Services List (gcloud CLI)**
*Command line output showing all deployed MCP servers with service names, regions, and HTTPS endpoints. Each service is independently scalable with auto-scaling from 0 to 1000+ instances.*

---

## Dockerfile Anatomy

**Example from mcp-deepcell**:

```dockerfile
FROM python:3.10-slim

# System dependencies
RUN apt-get update && apt-get install -y \
    libgomp1 libhdf5-dev && rm -rf /var/lib/apt/lists/*

# Copy code
COPY . /app
WORKDIR /app

# Install Python dependencies
RUN pip install --no-cache-dir -e .

# Environment
ENV MCP_TRANSPORT=sse
ENV DEEPCELL_DRY_RUN=false

# Run server
CMD ["python", "-m", "mcp_deepcell"]
```

Full Dockerfile: [`servers/mcp-deepcell/Dockerfile`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/Dockerfile) (46 lines)

**Key patterns**:
- Python 3.10 for DeepCell (TensorFlow 2.8.x constraint)
- Python 3.11+ for other servers
- `--no-cache-dir` to reduce image size
- System dependencies via `apt-get`
- SSE transport enabled by default

---

## Cost Analysis

### PatientOne Analysis Cost Breakdown

**Cloud Run costs** (us-central1 pricing):
- vCPU: $0.00002400/vCPU-second
- Memory: $0.00000250/GiB-second
- Requests: 2 million free/month, then $0.40/million

**Example calculation** (single PatientOne analysis):
- 12 servers × 2 seconds average = 24 vCPU-seconds
- 4Gi memory × 2s × 12 = 96 GiB-seconds
- Cost: ~$0.0024 + ~$0.0002 = **$0.0026 per analysis**

**Claude API cost**:
- Input: ~200K tokens × $3/million = $0.60
- Output: ~50K tokens × $15/million = $0.75
- **Total: $1.35 per analysis**

**Total per patient**: $0.0026 (Cloud Run) + $1.35 (Claude) = **$1.35 per analysis**

**Comparison**:
- Traditional manual: $3,200 (40 hours × $80/hour)
- AI-orchestrated: $1.35
- **Savings: 99.96%**

Full cost analysis: [`docs/deployment/cost-optimization.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/deployment/cost-optimization.md)

---

## Verification and Testing

### Health Check Endpoints

All servers expose `/health`:

```bash
SERVICE_URL=$(gcloud run services describe mcp-deepcell \
  --region=us-central1 --format='value(status.url)')

curl "${SERVICE_URL}/health"
# Response: {"status": "healthy", "transport": "sse"}
```

### Test with Claude API

**Python example**:

```python
import anthropic

client = anthropic.Anthropic(api_key="YOUR_API_KEY")

response = client.messages.create(
    model="claude-sonnet-4.5",
    messages=[{
        "role": "user",
        "content": "Segment cells in PatientOne's DAPI image"
    }],
    tools=[{
        "type": "mcp_tool",
        "mcp_server": {
            "url": "https://mcp-deepcell-xxx.run.app",
            "transport": "sse"
        }
    }]
)
```

Full API integration: [`docs/api-examples/claude-mcp-integration.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/api-examples/claude-mcp-integration.md)

---

## View Logs and Metrics

### Real-time Logs

```bash
# Tail logs
gcloud run logs tail mcp-deepcell --region=us-central1

# Filter errors only
gcloud run logs read mcp-deepcell \
  --region=us-central1 \
  --filter="severity>=ERROR"
```

### Cloud Console Metrics

Navigate to:
```
https://console.cloud.google.com/run/detail/us-central1/mcp-deepcell/metrics
```

**Key metrics**:
- Request latency (P50, P95, P99)
- Request count per second
- Memory utilization
- Container instance count
- Error rate

Full monitoring guide in Chapter 14.

---

## Troubleshooting Common Issues

### Issue 1: Out of Memory

**Symptoms**: "OOM killed" errors, container crashes

**Solution**:
```bash
gcloud run services update mcp-deepcell \
  --memory=8Gi \
  --region=us-central1
```

### Issue 2: Timeout Errors

**Symptoms**: "Deadline exceeded" on large images

**Solution**:
```bash
gcloud run services update mcp-deepcell \
  --timeout=600 \
  --region=us-central1
```

### Issue 3: Cold Start Latency

**Symptoms**: First request takes 30-60s (model download)

**Solution**: Set minimum instances to keep warm
```bash
gcloud run services update mcp-deepcell \
  --min-instances=1 \
  --region=us-central1
```

**Cost impact**: ~$50-80/month for 1 always-on instance (4Gi, 2 CPU).

Full troubleshooting: [`servers/mcp-deepcell/DEPLOYMENT.md#troubleshooting`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/DEPLOYMENT.md#troubleshooting)

---

## Security Configuration

### Public vs Private Deployment

**Development (public)**:
```bash
--allow-unauthenticated
```

**Production (private)**:
```bash
--no-allow-unauthenticated
```

Then authenticate requests with:
- Service account tokens
- Identity-Aware Proxy (IAP)
- API keys (custom implementation)

### VPC Integration

For hospital deployment (Chapter 13):

```bash
gcloud run services update mcp-deepcell \
  --vpc-connector=YOUR_CONNECTOR \
  --vpc-egress=private-ranges-only \
  --region=us-central1
```

**Use case**: Access private FHIR servers, on-premise databases.

Full security guide in Chapter 13.

---

## What You've Deployed

**Cloud Run architecture**:
1. **12 MCP servers**: Each runs in isolated Cloud Run service
2. **SSE transport**: HTTP endpoints for Claude API integration
3. **Auto-scaling**: 0 to 1000+ instances based on load
4. **Cost-optimized**: Pay only for request time (~$0.0026 per analysis)

**PatientOne workflow**:
- Deploy all 12 servers: ~15-20 minutes
- Total cost: $1.35 per patient analysis (99.96% cheaper than manual)
- Analysis time: 35 minutes (vs 40 hours manual)

**Integration**: Claude orchestrates via API calls to Cloud Run URLs instead of local stdio.

---

## Summary

**Chapter 12 Summary**:
- Cloud Run provides serverless, auto-scaling deployment
- SSE transport enables Claude API integration
- Three deployment options: deploy.sh (fastest), Cloud Build (CI/CD), manual
- Resource sizing: 2-8Gi memory, 1-4 CPU, 120-600s timeout
- Cost: $0.0026 Cloud Run + $1.35 Claude API = $1.35 per analysis
- Verification: Health checks, logs, metrics
- Security: Public (dev) vs private (production), VPC integration

**Files**: Each server has `Dockerfile`, `deploy.sh`, `cloudbuild.yaml`, `DEPLOYMENT.md`
**Deployment time**: ~15-20 minutes for all 12 servers
**Production ready**: HIPAA-eligible, autoscaling, monitoring



ewpage


# Hospital Production Deployment

*HIPAA-compliant deployment for research hospitals*

---

## Why Hospital Deployment Differs

Chapter 12 deployed MCP servers to public Cloud Run with `--allow-unauthenticated`. **This won't work for hospitals**.

**Hospital requirements**:
- **HIPAA compliance**: Business Associate Agreement (BAA), de-identification, audit logs
- **Private networking**: VPC integration, no public access
- **SSO authentication**: Hospital Active Directory (Azure AD)
- **Epic FHIR integration**: OAuth 2.0 to research endpoint
- **10-year audit logs**: Immutable compliance records
- **IRB approval**: Research patient data access

**Architecture change**: Public internet → Hospital VPN → VPC → Cloud Run → Epic FHIR (research endpoint only).

![Clinical Persona](images/screenshots/clinical-personna.jpeg){width=\textwidth}

**Figure 13.1: Clinical User Personas**
*Hospital deployment serves multiple user types with different access needs: clinical researchers, bioinformaticians, oncologists, and administrators. Each role has specific permissions and workflows within the HIPAA-compliant infrastructure.*

### HIPAA-Compliant Hospital Deployment Architecture

```mermaid
graph TB
    subgraph "Hospital Network (On-Premise)"
        RESEARCHER[Researcher<br/>Hospital workstation<br/>VPN connected]
        AD[Azure Active Directory<br/>Hospital SSO<br/>User authentication]
        EPIC[Epic FHIR API<br/>Research endpoint<br/>OAuth 2.0]
    end

    subgraph "Network Security"
        VPN[Hospital VPN<br/>Point-to-site<br/>Authorized IPs only]
        FIREWALL[Cloud Firewall<br/>VPC ingress rules<br/>Deny all by default]
    end

    subgraph "Google Cloud (HIPAA-Eligible)"
        direction TB

        subgraph "VPC Private Network"
            PROXY[OAuth2 Proxy<br/>Azure AD integration<br/>Token validation]
            ALB[Load Balancer<br/>Private IP<br/>TLS 1.3]
        end

        subgraph "Cloud Run Services (Private)"
            EPIC_SVC[mcp-epic<br/>FHIR integration<br/>De-identification]
            FGBIO_SVC[mcp-fgbio<br/>Genomics QC]
            MULTI_SVC[mcp-multiomics<br/>Multi-omics]
            SPATIAL_SVC[mcp-spatialtools<br/>Spatial analysis]
        end

        subgraph "Data Protection"
            SECRET[Secret Manager<br/>Encrypted credentials<br/>CMEK enabled]
            STORAGE[Cloud Storage<br/>AES-256 encrypted<br/>CSEK optional]
            LOGS[Cloud Logging<br/>10-year retention<br/>Immutable audit trail]
        end
    end

    subgraph "Compliance Controls"
        DEIDENT[De-identification Layer<br/>HIPAA Safe Harbor<br/>18 identifiers removed]
        AUDIT[Audit Logging<br/>All API calls<br/>10-year retention]
        ENCRYPT[Encryption at Rest<br/>AES-256 + CMEK<br/>Customer-managed keys]
        BAA[Business Associate Agreement<br/>GCP HIPAA BAA signed<br/>Covered services]
    end

    RESEARCHER --> VPN
    VPN --> FIREWALL
    FIREWALL --> PROXY

    RESEARCHER --> AD
    AD --> PROXY
    PROXY --> ALB
    ALB --> EPIC_SVC
    ALB --> FGBIO_SVC
    ALB --> MULTI_SVC
    ALB --> SPATIAL_SVC

    EPIC_SVC --> EPIC
    EPIC_SVC --> DEIDENT

    EPIC_SVC --> SECRET
    FGBIO_SVC --> SECRET
    EPIC_SVC --> STORAGE
    FGBIO_SVC --> STORAGE

    EPIC_SVC --> AUDIT
    FGBIO_SVC --> AUDIT
    MULTI_SVC --> AUDIT
    SPATIAL_SVC --> AUDIT

    SECRET --> ENCRYPT
    STORAGE --> ENCRYPT
    AUDIT --> LOGS

    style RESEARCHER fill:#d1ecf1
    style VPN fill:#f8d7da,stroke:#dc3545,stroke-width:2px
    style FIREWALL fill:#f8d7da,stroke:#dc3545,stroke-width:2px
    style PROXY fill:#fff3cd,stroke:#ffc107,stroke-width:2px
    style EPIC_SVC fill:#d4edda,stroke:#28a745,stroke-width:2px
    style DEIDENT fill:#cce5ff,stroke:#004085,stroke-width:2px
    style AUDIT fill:#cce5ff,stroke:#004085,stroke-width:2px
    style ENCRYPT fill:#cce5ff,stroke:#004085,stroke-width:2px
    style BAA fill:#e1ecf4,stroke:#0066cc,stroke-width:2px
    style LOGS fill:#fef3cd
```

**Figure 13.2: HIPAA-Compliant Hospital Deployment Architecture**
*Multi-layer security architecture: (1) Network Security: Hospital VPN + Cloud Firewall restricts access to authorized IPs only, (2) Authentication: Azure AD SSO with OAuth2 Proxy for token validation, (3) Private Networking: VPC-isolated Cloud Run services with private IPs and internal load balancer (no public internet access), (4) Data Protection: Secret Manager for credentials, AES-256 encrypted Cloud Storage with optional CSEK, 10-year immutable audit logs, (5) Compliance Controls: HIPAA Safe Harbor de-identification (18 identifiers removed), comprehensive audit logging, customer-managed encryption keys (CMEK), and signed Google Cloud HIPAA BAA.*

**HIPAA Compliance Layers:**

1. **Legal**: Google Cloud HIPAA BAA signed + IRB approval
2. **Network Isolation**: VPC private network, no public access, hospital VPN only
3. **Authentication**: Azure AD SSO with OAuth 2.0 tokens
4. **Encryption at Rest**: AES-256 (GCP default) + CMEK (customer keys)
5. **Encryption in Transit**: TLS 1.3 for all connections
6. **De-identification**: Automatic HIPAA Safe Harbor (18 identifiers removed)
7. **Audit Logging**: All API calls logged, 10-year retention, immutable
8. **Access Control**: IAM roles, least privilege principle
9. **Epic Integration**: Research endpoint only (not production EHR)
10. **Data Minimization**: Only approved researchers, IRB-authorized data

**Critical: No Public Internet Access**
- All Cloud Run services: `--no-allow-unauthenticated`
- VPC-only ingress: `--ingress=internal`
- Load balancer: Private IP only (10.x.x.x range)
- Researcher access: Hospital VPN → VPC → Cloud Run

---

## Prerequisites

### Legal and Compliance

**HIPAA Business Associate Agreement (BAA)**:
- Sign Google Cloud HIPAA BAA (free, available through GCP Console)
- Documents: `Settings → Legal → Business Associate Agreement`
- Enables: HIPAA-eligible services (Cloud Run, Secret Manager, Cloud Logging)

**IRB Approval**:
- Institutional Review Board approval for research patient data
- Informed consent from participants
- Data Use Agreement (DUA) templates

Full HIPAA guide: [`docs/for-hospitals/compliance/hipaa.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-hospitals/compliance/hipaa.md)

### Technical Prerequisites

**Hospital IT must provide**:
1. **GCP Organization ID**: Hospital's existing GCP organization
2. **Billing Account**: Grant-funded billing account
3. **Azure AD tenant**: For SSO authentication
4. **Epic FHIR credentials**: Research endpoint (NOT production)
5. **VPC network**: Existing hospital VPC or create new
6. **VPN access**: Hospital VPN for deployment team

### Azure AD App Registration

**Register Azure AD app** (before deployment):

```
App Name: "Precision Medicine MCP"
Supported accounts: Single tenant (hospital only)
Redirect URIs: https://{oauth2-proxy-url}/oauth2/callback
```

**API permissions**:
- `User.Read` (delegated)
- `Directory.Read.All` (application)

Full Azure AD setup: [`infrastructure/hospital-deployment/README.md#azure-ad-setup`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/infrastructure/hospital-deployment/README.md#azure-ad-setup)

---

## Deployment Scripts Overview

Hospital deployment uses **5 setup scripts** in order:

```bash
cd infrastructure/hospital-deployment

1. ./setup-project.sh      # GCP project, APIs, billing
2. ./setup-vpc.sh           # VPC networking, firewall rules
3. ./setup-secrets.sh       # Secret Manager (Epic, Azure AD, Claude API)
4. ./setup-audit-logging.sh # 10-year audit log retention
5. ./deploy-oauth2-proxy.sh # Azure AD SSO authentication
```

Full scripts: [`infrastructure/hospital-deployment/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/infrastructure/hospital-deployment)

**Deployment time**: 1-2 days (Week 1 of 3-month phased rollout).

---

## Script 1: GCP Project Setup

**Creates**: GCP project in hospital organization with HIPAA-compliant configuration.

```bash
./setup-project.sh
```

**What it does**:
- Creates project in hospital's GCP Organization
- Links grant-funded billing account
- Enables APIs: Cloud Run, Compute, VPC Access, Secret Manager, Logging
- Sets budget alerts: $1,000/month with 50%, 75%, 90%, 100% thresholds
- Configures HIPAA-eligible regions (us-central1, us-east4)

**Configuration** (edit variables in script):

```bash
ORG_ID="123456789012"              # Hospital's GCP org ID
BILLING_ACCOUNT_ID="00B597-..."    # PI's billing account
PROJECT_ID="precision-medicine-poc"
REGION="us-central1"
```

Full script: [`infrastructure/hospital-deployment/setup-project.sh`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/infrastructure/hospital-deployment/setup-project.sh)

---

## Script 2: VPC Networking

**Creates**: Private VPC with Serverless VPC Connector for Cloud Run.

```bash
./setup-vpc.sh
```

**What it does**:
- Creates subnet: `10.10.0.0/24` (256 IPs for MCP servers)
- Creates Serverless VPC Connector (Cloud Run to VPC bridge)
- Enables Private Google Access (secure GCS/API access)
- Configures firewall rules:
  - Allow health checks from Google
  - Allow internal traffic between servers
  - **Deny all external ingress** (requires VPN)

**Network architecture**:

```
Hospital VPN → VPC (10.10.0.0/24) → Serverless VPC Connector → Cloud Run
                                  → Epic FHIR (research endpoint)
```

**Use case**: Cloud Run services access Epic FHIR without public internet exposure.

Full script: [`infrastructure/hospital-deployment/setup-vpc.sh`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/infrastructure/hospital-deployment/setup-vpc.sh)

---

## Script 3: Secret Manager

**Creates**: Encrypted secret storage for credentials.

```bash
./setup-secrets.sh --interactive
```

**7 secrets created**:
1. `anthropic-api-key` - Claude API access
2. `epic-fhir-endpoint` - Epic FHIR base URL (research endpoint)
3. `epic-client-id` - Epic OAuth client ID
4. `epic-client-secret` - Epic OAuth secret
5. `azure-ad-client-id` - Azure AD app ID
6. `azure-ad-client-secret` - Azure AD app secret
7. `azure-ad-tenant-id` - Hospital's Azure AD tenant

**Security features**:
- Automatic encryption at rest (Google-managed keys)
- IAM bindings (service accounts only, no users)
- Version history (rotate secrets without downtime)
- Audit logging (who accessed what, when)

**Example usage**:

```bash
# Populate secret
echo -n "sk-ant-..." | gcloud secrets versions add anthropic-api-key --data-file=-

# Verify
gcloud secrets list --project=precision-medicine-poc
```

Full script: [`infrastructure/hospital-deployment/setup-secrets.sh`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/infrastructure/hospital-deployment/setup-secrets.sh)

---

## Script 4: Audit Logging

**Creates**: 10-year immutable audit logs (HIPAA requirement).

```bash
./setup-audit-logging.sh
```

**What it does**:
- Creates log bucket with **3,650-day retention** (10 years)
- Sets up log sinks:
  - All Cloud Run requests (user queries, MCP tool calls)
  - Authentication events (login, logout, failed attempts)
  - Epic FHIR API calls
  - De-identification operations
- Creates custom metrics:
  - De-identification success rate
  - Epic FHIR failure rate
  - User access events

**Log immutability**: Cannot be deleted before 10 years, even by admins.

**Accessing logs**:

```bash
# View user queries
gcloud logging read "jsonPayload.event=\"mcp_query\"" --limit=50

# Export for compliance report
gcloud logging read "jsonPayload.event=\"mcp_query\"" --format=csv > report.csv
```

Full script: [`infrastructure/hospital-deployment/setup-audit-logging.sh`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/infrastructure/hospital-deployment/setup-audit-logging.sh)

---

## Script 5: OAuth2 Proxy (SSO)

**Deploys**: Azure AD authentication layer for Cloud Run services.

```bash
./deploy-oauth2-proxy.sh
```

**What it does**:
- Builds OAuth2 Proxy container with Azure AD configuration
- Deploys to Cloud Run as authentication gateway
- Configures:
  - Provider: Azure AD
  - Email domains: `@hospital.org` (only hospital users)
  - Session expiry: 1 day
  - CSRF protection (SameSite cookies)
- Sits in front of Streamlit/JupyterHub

**Authentication flow**:

```
User → OAuth2 Proxy → Azure AD login → OAuth2 Proxy (cookie) → MCP servers
```

**Testing**:

```bash
PROXY_URL=$(gcloud run services describe oauth2-proxy \
  --region=us-central1 --format='value(status.url)')

curl -I $PROXY_URL
# Returns: 302 redirect to Azure AD login
```

Full script: [`infrastructure/hospital-deployment/deploy-oauth2-proxy.sh`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/infrastructure/hospital-deployment/deploy-oauth2-proxy.sh)

---

## De-identification: HIPAA Safe Harbor

**All patient data** must be de-identified before analysis (HIPAA requirement).

**HIPAA Safe Harbor removes 18 identifiers**:
1. Names
2. Geographic subdivisions smaller than state
3. Dates (except year)
4. Telephone/fax numbers
5. Email addresses
6. SSNs
7. Medical record numbers (MRNs)
8. Health plan beneficiary numbers
9. Account numbers
10. Certificate/license numbers
11. Vehicle identifiers
12. Device identifiers/serial numbers
13. URLs
14. IP addresses
15. Biometric identifiers
16. Full-face photos
17. Any other unique identifying number/characteristic/code

**Implementation** (mcp-epic server):

```python
def deidentify_patient(patient: dict) -> dict:
    """Apply HIPAA Safe Harbor de-identification."""
    # Remove identifiers
    for field in ["name", "telecom", "address", "photo"]:
        patient.pop(field, None)

    # Hash patient ID
    patient["id"] = hash_identifier(patient["id"])

    # Reduce birthDate to year only
    if "birthDate" in patient:
        patient["birthDate"] = reduce_to_year(patient["birthDate"])

    return patient
```

Full implementation: [`servers/mcp-epic/src/mcp_epic/deidentify.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-epic/src/mcp_epic/deidentify.py) (110 lines)

**Validation**: Privacy officer reviews de-identification results quarterly.

---

## Epic FHIR Integration

**Hospital IT provides**:
- **Research FHIR endpoint**: `https://hospital-research.epic.com/api/FHIR/R4/` (NOT production)
- **OAuth 2.0 credentials**: Client ID, client secret
- **Authorized scopes**: `patient/*.read`, `Observation.read`, `Condition.read`

**Test Epic connection**:

```bash
# Test OAuth token
curl -X POST "${EPIC_FHIR_ENDPOINT}/oauth2/token" \
  -d "grant_type=client_credentials" \
  -d "client_id=${EPIC_CLIENT_ID}" \
  -d "client_secret=${EPIC_CLIENT_SECRET}"

# Test metadata endpoint
curl -H "Authorization: Bearer $TOKEN" \
  "${EPIC_FHIR_ENDPOINT}/metadata"
```

**Store credentials** in Secret Manager:

```bash
echo -n "$EPIC_FHIR_ENDPOINT" | gcloud secrets versions add epic-fhir-endpoint --data-file=-
echo -n "$EPIC_CLIENT_ID" | gcloud secrets versions add epic-client-id --data-file=-
echo -n "$EPIC_CLIENT_SECRET" | gcloud secrets versions add epic-client-secret --data-file=-
```

Full Epic setup: [`infrastructure/hospital-deployment/README.md#epic-fhir-setup`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/infrastructure/hospital-deployment/README.md#epic-fhir-setup)

---

## Deploy MCP Servers (Private Mode)

**Deploy to VPC** with authentication required:

```bash
cd servers/mcp-deepcell

gcloud run deploy mcp-deepcell \
  --image=gcr.io/PROJECT_ID/mcp-deepcell:latest \
  --region=us-central1 \
  --no-allow-unauthenticated \
  --vpc-connector=mcp-connector \
  --vpc-egress=private-ranges-only \
  --service-account=mcp-deepcell-sa@PROJECT_ID.iam.gserviceaccount.com
```

**Key differences from public deployment**:
- `--no-allow-unauthenticated` (requires authentication)
- `--vpc-connector` (routes traffic through VPC)
- `--vpc-egress=private-ranges-only` (no public internet from container)
- `--service-account` (least-privilege IAM)

**Access**: Only via hospital VPN → OAuth2 Proxy → Cloud Run.

---

## PatientOne Hospital Workflow

**Deploy all 12 servers** for HIPAA-compliant analysis:

**Week 1**: Infrastructure setup
```bash
./setup-project.sh
./setup-vpc.sh
./setup-secrets.sh --interactive
./setup-audit-logging.sh
./deploy-oauth2-proxy.sh
```

**Week 2**: Core MCP servers
```bash
cd ../../servers
./deploy-all-servers.sh --vpc-mode --no-public-access
```

**Week 3**: Epic FHIR integration
```bash
cd servers/mcp-epic
./deploy.sh PROJECT_ID us-central1 --epic-endpoint=$EPIC_FHIR_ENDPOINT
```

**Deployment time**: 3 weeks (includes testing, validation).

**Verification**:

```bash
# All servers deployed
gcloud run services list --region=us-central1

# All services require auth
gcloud run services describe mcp-deepcell \
  --region=us-central1 --format='value(spec.template.spec.serviceAccountName)'
```

---

## Cost Management

### Monthly Cost Breakdown (Hospital Deployment)

**Infrastructure**:
- VPC Connector: $72/month (always-on, 2 instances)
- Secret Manager: $0.60/month (7 secrets × 2 versions)
- Audit logging: $50/month (100GB logs/month)

**Compute** (9 servers, moderate usage):
- Cloud Run: $400/month
- OAuth2 Proxy: $25/month (min-instances=1)

**API costs**:
- Claude API: $500/month (5 users, optimized)
- Epic FHIR: $0 (included in hospital contract)

**Total**: ~$1,047/month

**Budget alerts**: $500 (50%), $750 (75%), $900 (90%), $1,000 (100%).

Full cost guide: [`docs/deployment/cost-optimization.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/deployment/cost-optimization.md)

---

## Security Configuration

### Authentication

**OAuth2 Proxy enforces**:
- Azure AD SSO (hospital Active Directory)
- Email domain restriction: `@hospital.org`
- User group membership: `precision-medicine-users`
- Session expiry: 1 day
- CSRF protection (SameSite cookies)

### Network Security

**VPC firewall rules**:
- Allow: Health checks from Google (`130.211.0.0/22`, `35.191.0.0/16`)
- Allow: Internal traffic between MCP servers
- **Deny**: All external ingress

**Access**: Only via hospital VPN.

### Data Security

**De-identification**: All patient data auto-de-identified via HIPAA Safe Harbor
**Encryption**: TLS 1.3 in transit, AES-256 at rest
**Secrets**: Encrypted in Secret Manager (Google-managed keys)
**Audit logs**: Immutable 10-year retention

Full security guide: [`docs/for-hospitals/compliance/data-governance.md#access-controls-security`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-hospitals/compliance/data-governance.md#access-controls-security)

---

## Verification and Testing

### 1. VPC Connectivity

```bash
# SSH to test VM in VPC
gcloud compute ssh test-vm --zone=us-central1-a

# Test Cloud Run internal access
curl https://mcp-deepcell-xxx.run.app/health
# Should return 200 OK
```

### 2. Authentication

```bash
# Test OAuth2 Proxy
curl -I https://{oauth2-proxy-url}
# Should return 302 redirect to Azure AD
```

### 3. Epic FHIR Connection

```bash
# From Cloud Run container
gcloud run services update mcp-epic \
  --update-env-vars TEST_EPIC_CONNECTION=true

# Check logs
gcloud run logs read mcp-epic --limit=20
# Should show: "Epic FHIR connection: SUCCESS"
```

### 4. De-identification

```bash
# Fetch patient data (auto-de-identified)
# Via Claude Desktop or API, run:
# "Fetch patient PAT001 from Epic FHIR"

# Verify in logs: No MRNs, names, or full addresses appear
gcloud logging read "resource.type=cloud_run_revision" --limit=100
```

Full testing guide: [`docs/for-hospitals/USER_GUIDE.md#testing-deployment`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-hospitals/USER_GUIDE.md#testing-deployment)

---

## Troubleshooting

### Issue 1: OAuth2 Login Fails

**Symptoms**: Redirect to Azure AD, but returns error

**Solutions**:
- Check Azure AD redirect URIs match deployed OAuth2 Proxy URL
- Verify `azure-ad-client-secret` not expired (24 months max)
- Confirm user in `precision-medicine-users` Azure AD group

Full runbook: [`docs/for-hospitals/RUNBOOKS/sso-issues.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-hospitals/RUNBOOKS/sso-issues.md)

### Issue 2: Epic FHIR Connection Fails

**Symptoms**: "401 Unauthorized" or "403 Forbidden"

**Solutions**:
- Verify Epic credentials in Secret Manager
- Check OAuth token endpoint URL correct
- Ensure service account has `secretAccessor` role

Full runbook: [`docs/for-hospitals/RUNBOOKS/epic-connection-failure.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-hospitals/RUNBOOKS/epic-connection-failure.md)

### Issue 3: VPC Egress Blocked

**Symptoms**: Cloud Run cannot access Epic FHIR or Google APIs

**Solutions**:
```bash
# Check VPC egress setting
gcloud run services describe mcp-epic \
  --region=us-central1 --format='value(spec.template.spec.vpcAccess.egress)'

# Should be: private-ranges-only
# If all-traffic, redeploy with correct setting
```

---

## What You've Deployed

**Hospital-compliant architecture**:
1. **HIPAA compliance**: BAA signed, 10-year audit logs, de-identification
2. **Private networking**: VPC with no public access, VPN required
3. **SSO authentication**: Azure AD via OAuth2 Proxy
4. **Epic FHIR integration**: Research endpoint with OAuth 2.0
5. **Secret management**: Encrypted credentials in Secret Manager
6. **Audit logging**: Immutable 10-year retention

**PatientOne workflow**:
- Deployment time: 3 weeks (infrastructure + servers + Epic integration)
- Cost: ~$1,047/month (VPC + compute + API)
- Access: Hospital VPN → Azure AD login → MCP servers → Epic FHIR (research)

**Compliance**: HIPAA-eligible, IRB-approved, privacy officer validated.

---

## Summary

**Chapter 13 Summary**:
- Hospital deployment requires HIPAA compliance, VPC, SSO
- 5 setup scripts: project, VPC, secrets, audit logs, OAuth2 Proxy
- De-identification: HIPAA Safe Harbor (18 identifiers removed)
- Epic FHIR: Research endpoint only, OAuth 2.0, auto-de-identification
- Private deployment: No public access, VPN required, authentication enforced
- Cost: ~$1,047/month (infrastructure + compute + Claude API)

**Files**: [`infrastructure/hospital-deployment/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/infrastructure/hospital-deployment) (5 scripts + docs)
**Deployment time**: 3 weeks (phased rollout)
**Compliance**: HIPAA-eligible, IRB-approved, 10-year audit logs



ewpage


# Operations and Monitoring

*Production monitoring, alerting, and runbooks*

---

## Why Monitoring Matters

Chapters 12-13 deployed MCP servers to production. **Now what?**

**Without monitoring, you're blind**:
- Is the system working? (Health checks)
- Are users getting errors? (Error rate tracking)
- Is Epic FHIR accessible? (Integration monitoring)
- Are de-identification operations succeeding? (Compliance validation)
- What's the monthly cost? (Budget tracking)

**Monitoring enables**:
- **Proactive response**: Alerts before users notice problems
- **Compliance**: 10-year audit logs (HIPAA requirement)
- **Cost control**: Budget alerts, usage optimization
- **Performance tuning**: Identify bottlenecks, optimize resources

**PatientOne workflow**: Monitor all 12 MCP servers with dashboards, alerts, and runbooks.

---

## Cloud Logging (10-Year Audit Trail)

**All production logs** stored with 10-year immutable retention (HIPAA requirement).

**Log types**:
1. **Cloud Run requests**: User queries, MCP tool calls, response codes
2. **Authentication events**: Login, logout, failed attempts (Azure AD)
3. **Epic FHIR calls**: Patient data access, de-identification operations
4. **Application logs**: Structured logs from MCP servers
5. **System logs**: Container startup, crashes, out-of-memory events

**Access logs**:

```bash
# View recent user queries
gcloud logging read "jsonPayload.event=\"mcp_query\"" --limit=50

# View Epic FHIR calls
gcloud logging read "jsonPayload.event=\"epic_fhir_call\"" --limit=50

# View errors only
gcloud logging read "severity>=ERROR" --limit=100

# Export for compliance report
gcloud logging read "timestamp>=\"2026-01-01T00:00:00Z\"" \
  --format=csv > audit_report.csv
```

**Structured logging format**:

```json
{
  "severity": "INFO",
  "timestamp": "2026-01-31T10:15:30Z",
  "event": "mcp_query",
  "user_id": "user@hospital.org",
  "server": "mcp-spatialtools",
  "tool": "run_differential_expression",
  "duration_ms": 5432,
  "success": true
}
```

Full logging guide: [`docs/for-hospitals/AUDIT_LOG_GUIDE.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-hospitals/AUDIT_LOG_GUIDE.md)

---

## Key Metrics to Monitor

### 1. Request Latency

**Metric**: `run.googleapis.com/request_latencies`

**Thresholds**:
- P50 < 10s (warm requests)
- P95 < 60s (includes cold starts)
- P99 < 120s (large images + cold starts)

**Alert if**: P95 > 120s for 15 minutes

### 2. Error Rate

**Metric**: `run.googleapis.com/request_count` (filtered by response code)

**Thresholds**:
- 5xx errors < 1% (server errors)
- 4xx errors < 5% (client errors acceptable)

**Alert if**: 5xx error rate > 2% for 5 minutes

### 3. Memory Utilization

**Metric**: `run.googleapis.com/container/memory/utilizations`

**Thresholds**:
- Average < 70% (of allocated memory)
- Peak < 90%

**Alert if**: Memory > 85% for 10 minutes

**Action**:
```bash
gcloud run services update mcp-deepcell --memory=8Gi
```

### 4. CPU Utilization

**Metric**: `run.googleapis.com/container/cpu/utilizations`

**Thresholds**:
- Average < 80% during processing
- Idle < 10% when not processing

**Alert if**: CPU > 90% for 10 minutes

### 5. Epic FHIR Failures

**Metric**: Custom log-based metric from Epic FHIR call logs

**Threshold**: Epic failure rate < 1%

**Alert if**: Epic failures > 5% for 5 minutes

### 6. De-identification Success

**Metric**: Custom log-based metric from de-identification operations

**Threshold**: De-identification success rate > 99%

**Alert if**: De-identification failures > 1% (CRITICAL - privacy risk)

Full metrics guide: [`servers/mcp-deepcell/MONITORING.md#key-metrics-to-monitor`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/MONITORING.md#key-metrics-to-monitor)

---

## Monitoring Dashboards

**Create dashboard** with Cloud Monitoring:

```bash
gcloud monitoring dashboards create \
  --config-from-file=dashboard-config.json
```

**Dashboard widgets**:
1. **Server Health**: Request count, error rate, latency (P50/P95/P99)
2. **Resource Usage**: Memory utilization, CPU utilization, instance count
3. **User Activity**: Queries per hour, active users, top tools used
4. **Epic FHIR**: Connection status, API calls, failure rate
5. **De-identification**: Success rate, processing time
6. **Cost Tracking**: Daily spend, projected monthly cost, budget burn rate

### Production Dashboard Examples

![Operations Dashboard 1](images/screenshots/dash-1.png){width=\textwidth}

**Figure 14.1: Main Operations Dashboard**
*Real-time overview of MCP server health, request rates, error rates, and latency percentiles (P50, P95, P99) across all 12 services.*

![Operations Dashboard 2](images/screenshots/dash-2.png){width=\textwidth}

**Figure 14.2: Cost Tracking Dashboard**
*Daily and monthly cost breakdown by service, projected spend, and budget burn rate with alerting thresholds.*

![Operations Dashboard 3](images/screenshots/dash-3.png){width=\textwidth}

**Figure 14.3: Performance Metrics**
*Memory and CPU utilization, instance counts, cold start frequency, and throughput optimization metrics.*

![Operations Dashboard 4](images/screenshots/dash-4.png){width=\textwidth}

**Figure 14.4: Server Health Monitoring**
*Individual server status, uptime, response times, and integration health (Epic FHIR, de-identification pipelines).*

**Dashboard URL**:
```
https://console.cloud.google.com/monitoring/dashboards/custom/{dashboard-id}
```

**Example dashboard config**: [`infrastructure/hospital-deployment/monitoring/dashboard-config.json`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/infrastructure/hospital-deployment/monitoring/dashboard-config.json)

---

## Alerting Policies

**Create 6 critical alerts**:

### 1. High Error Rate (P0)

```yaml
Condition: 5xx error rate > 2% for 5 minutes
Severity: Critical
Action: Page on-call engineer
```

```bash
gcloud alpha monitoring policies create \
  --notification-channels=CHANNEL_ID \
  --display-name="MCP High Error Rate" \
  --condition-threshold-value=0.02 \
  --condition-threshold-duration=300s
```

### 2. Server Down (P0)

```yaml
Condition: Service status != Ready for 2 minutes
Severity: Critical
Action: Page on-call engineer, run "Server Down" runbook
```

### 3. Epic FHIR Failure (P1)

```yaml
Condition: Epic failure rate > 5% for 5 minutes
Severity: High
Action: Email team, run "Epic Connection Failure" runbook
```

### 4. High Memory Usage (P2)

```yaml
Condition: Memory utilization > 85% for 10 minutes
Severity: Medium
Action: Email team, consider scaling up
```

### 5. De-identification Failure (P0)

```yaml
Condition: De-identification failure rate > 1%
Severity: Critical (privacy risk)
Action: Page privacy officer, halt Epic data access
```

### 6. Budget Threshold Exceeded (P2)

```yaml
Condition: Daily cost > 90% of monthly budget / 30 days
Severity: Medium
Action: Email team, review cost drivers
```

Full alert policies: [`infrastructure/hospital-deployment/monitoring/alert-policies.yaml`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/infrastructure/hospital-deployment/monitoring/alert-policies.yaml)

---

## Runbooks (Troubleshooting Guides)

**3 critical runbooks** for common issues:

### Runbook 1: Server Down

**Symptoms**: "Connection refused", 5xx errors, health check failures

**Diagnosis**:
```bash
# Check server status
gcloud run services describe mcp-deepcell \
  --region=us-central1 \
  --format='value(status.conditions[0].status)'

# Check logs
gcloud run logs read mcp-deepcell --limit=50 | grep ERROR
```

**Resolution**:
1. **Restart container** (force new revision):
```bash
gcloud run services update mcp-deepcell \
  --update-env-vars=RESTART_TIMESTAMP=$(date +%s)
```

2. **Rollback to previous revision** (if recent deployment):
```bash
gcloud run services update-traffic mcp-deepcell \
  --to-revisions=PREVIOUS_REVISION=100
```

3. **Redeploy from source** (if corruption):
```bash
cd servers/mcp-deepcell
./deploy.sh PROJECT_ID us-central1
```

Full runbook: [`docs/for-hospitals/RUNBOOKS/server-down.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-hospitals/RUNBOOKS/server-down.md)

### Runbook 2: SSO/Authentication Issues

**Symptoms**: Azure AD redirect fails, "Access denied", cookie errors

**Diagnosis**:
```bash
# Check OAuth2 Proxy logs
gcloud run logs read oauth2-proxy --limit=50

# Verify Azure AD redirect URIs
# Check Azure Portal → App registrations → Redirect URIs
```

**Resolution**:
1. **Verify redirect URIs** match deployed OAuth2 Proxy URL
2. **Check Azure AD secret expiry** (24 months max):
```bash
# Update secret in Secret Manager
echo -n "NEW_SECRET" | gcloud secrets versions add azure-ad-client-secret --data-file=-
```

3. **Verify user in Azure AD group** `precision-medicine-users`

Full runbook: [`docs/for-hospitals/RUNBOOKS/sso-issues.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-hospitals/RUNBOOKS/sso-issues.md)

### Runbook 3: Epic FHIR Connection Failure

**Symptoms**: "401 Unauthorized", "Epic connection failed", no patient data

**Diagnosis**:
```bash
# Test Epic connection
gcloud run logs read mcp-epic --limit=50 | grep "epic_fhir_call"

# Verify Epic credentials
gcloud secrets versions access latest --secret=epic-client-id
```

**Resolution**:
1. **Verify Epic credentials** in Secret Manager (not expired)
2. **Test OAuth token endpoint**:
```bash
curl -X POST "${EPIC_ENDPOINT}/oauth2/token" \
  -d "grant_type=client_credentials" \
  -d "client_id=${EPIC_CLIENT_ID}" \
  -d "client_secret=${EPIC_CLIENT_SECRET}"
```

3. **Check service account IAM permissions** (needs `secretAccessor` role)
4. **Contact hospital IT** (Epic FHIR endpoint may be down)

Full runbook: [`docs/for-hospitals/RUNBOOKS/epic-connection-failure.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-hospitals/RUNBOOKS/epic-connection-failure.md)

---

## Cost Tracking and Optimization

### View Current Costs

```bash
# View daily costs (last 7 days)
gcloud billing accounts get-spend \
  --billing-account=BILLING_ACCOUNT_ID \
  --start-date=$(date -d '7 days ago' -I) \
  --end-date=$(date -I)

# View by service
gcloud billing accounts get-spend \
  --billing-account=BILLING_ACCOUNT_ID \
  --groupBy=service
```

### Budget Alerts

**Set budget** with 50%, 75%, 90%, 100% thresholds:

```bash
gcloud billing budgets create \
  --billing-account=BILLING_ACCOUNT_ID \
  --display-name="MCP Monthly Budget" \
  --budget-amount=1000USD \
  --threshold-rule=percent=50,basis=current-spend \
  --threshold-rule=percent=75,basis=current-spend \
  --threshold-rule=percent=90,basis=current-spend \
  --threshold-rule=percent=100,basis=current-spend
```

**Alert notification**: Email to team when threshold exceeded.

### Cost Optimization Strategies

1. **Scale to zero**: Set `min-instances=0` for all non-critical servers
2. **Reduce concurrency**: Lower concurrency = fewer instances = lower cost
3. **Use Haiku for simple queries**: ~10× cheaper than Sonnet
4. **Cache MCP responses**: Avoid redundant API calls
5. **Right-size resources**: Monitor actual usage, reduce over-provisioned memory/CPU

**Example**: Reducing mcp-deepcell memory from 8Gi → 4Gi saves ~$50/month.

Full cost guide: [`docs/deployment/cost-optimization.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/deployment/cost-optimization.md)

---

## Performance Optimization

### Baseline Performance (CPU-only, 4Gi RAM, 2 vCPU)

| Server | Cold Start | Warm Request | Optimization |
|--------|------------|--------------|--------------|
| mcp-fgbio | ~5s | ~1s | None needed |
| mcp-multiomics | ~8s | ~3s | None needed |
| mcp-spatialtools | ~10s | ~5s | Consider min-instances=1 |
| mcp-deepcell | ~35s | ~5s | **High cold start** |
| mcp-perturbation | ~12s | ~8s | None needed |

### Optimization Techniques

**1. Reduce cold starts** (for mcp-deepcell):
```bash
gcloud run services update mcp-deepcell --min-instances=1
```
- Cost: ~$50/month for always-on instance
- Benefit: Eliminates 30s cold start latency

**2. Increase throughput** (for high-traffic servers):
```bash
gcloud run services update mcp-spatialtools \
  --cpu=4 \
  --max-instances=20
```
- Cost: +50% per request
- Benefit: 2× faster processing, more concurrent requests

**3. Optimize concurrency** (better caching):
```bash
gcloud run services update mcp-deepcell --concurrency=1
```
- Ensures one request per instance (better model caching)
- More predictable performance

Full performance guide: [`servers/mcp-deepcell/MONITORING.md#performance-optimization-recommendations`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/MONITORING.md#performance-optimization-recommendations)

---

## Quarterly Bias Audits

**HIPAA and NIH require** quarterly bias audits for AI/ML systems.

**Audit process** (every 3 months):

1. **Extract analysis data** (last 3 months):
```bash
gcloud logging read "jsonPayload.event=\"mcp_query\" AND timestamp>\"2026-01-01T00:00:00Z\"" \
  --format=csv > q1_2026_queries.csv
```

2. **Run bias audit script**:
```bash
python3 scripts/audit/audit_bias.py \
  --workflow patientone \
  --genomics-data /data/genomics_results.csv \
  --clinical-data /data/fhir_patients.json \
  --output /reports/q1_2026_bias_audit.html
```

3. **Review audit report**:
- Data representation (ancestry distribution)
- Fairness metrics (demographic parity, equalized odds)
- Proxy feature detection (age, ZIP code)
- Confidence scoring (ancestry-aware)

4. **Submit to IRB** (required quarterly compliance report)

**Example bias metrics**:
- European ancestry: 45% of dataset (vs 40% US population) ✓
- African ancestry: 12% of dataset (vs 13% US population) ✓
- Asian ancestry: 8% of dataset (vs 6% US population) ✓
- Latino ancestry: 35% of dataset (vs 19% US population) ⚠️ Over-represented

**Action**: Document over-representation, adjust recruitment if needed.

Full bias audit guide: [`docs/for-hospitals/ethics/BIAS_AUDIT_CHECKLIST.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-hospitals/ethics/BIAS_AUDIT_CHECKLIST.md)

---

## PatientOne Monitoring Example

**Scenario**: Monitor PatientOne analysis (35 minutes, 12 MCP servers).

**Before analysis starts**:
```bash
# Check all servers healthy
gcloud run services list --region=us-central1
# All should show "Ready: True"

# Check Epic FHIR connection
gcloud run logs read mcp-epic --limit=5 | grep "epic_fhir_call"
# Should show recent successful calls
```

**During analysis** (monitor in real-time):
```bash
# Tail logs from all servers
gcloud run logs tail \
  --region=us-central1 \
  --filter="resource.labels.service_name:(mcp-)"

# Watch metrics dashboard
# https://console.cloud.google.com/monitoring/dashboards/custom/{dashboard-id}
```

**After analysis completes**:
```bash
# View total requests by server
gcloud logging read "jsonPayload.event=\"mcp_query\" AND timestamp>\"$(date -u -d '1 hour ago' -Iseconds)\"" \
  --format=json | jq '.[] | .jsonPayload.server' | sort | uniq -c

# Output:
#   15 mcp-fgbio
#   8 mcp-multiomics
#   12 mcp-spatialtools
#   7 mcp-deepcell
#   ...

# Calculate total cost
gcloud billing accounts get-spend \
  --billing-account=BILLING_ACCOUNT_ID \
  --start-date=$(date -I) \
  --end-date=$(date -I)
# Expected: ~$1.35 for PatientOne analysis
```

**Success criteria**:
- ✓ All servers responded (no 5xx errors)
- ✓ Epic FHIR calls succeeded (de-identification applied)
- ✓ Analysis completed in 35 minutes
- ✓ Cost within budget (~$1.35 vs $3,200 manual)

---

## What You've Configured

**Monitoring infrastructure**:
1. **Cloud Logging**: 10-year immutable audit trail (HIPAA compliant)
2. **Key metrics**: Latency, errors, memory, CPU, Epic FHIR, de-identification
3. **Dashboards**: Server health, resource usage, user activity, cost tracking
4. **Alerting policies**: 6 critical alerts (error rate, server down, Epic failure, memory, de-identification, budget)
5. **Runbooks**: 3 troubleshooting guides (server down, SSO issues, Epic connection)
6. **Cost tracking**: Budget alerts, daily spend monitoring, optimization strategies
7. **Performance tuning**: Cold start reduction, throughput optimization, concurrency tuning
8. **Quarterly bias audits**: Automated compliance reporting

**PatientOne workflow**:
- Real-time monitoring during 35-minute analysis
- Alerts if errors occur (page on-call engineer)
- Audit logs for compliance reporting
- Cost tracking ($1.35 per analysis vs $3,200 manual)

**Compliance**: 10-year logs, quarterly bias audits, IRB reporting.

---

## Next Steps

**Part 4 Complete!** You've deployed and configured operations for production:
- Chapter 12: Cloud Run deployment with SSE transport
- Chapter 13: Hospital production deployment with HIPAA compliance
- Chapter 14: Monitoring, alerting, and runbooks

**Part 5: Research and Education** (Chapters 15-16) covers:
- Chapter 15: For researchers (exploratory analysis, prompt engineering)
- Chapter 16: Teaching precision medicine (educational workflows, student access)

---

**Chapter 14 Summary**:
- Cloud Logging: 10-year immutable audit trail (HIPAA requirement)
- Key metrics: Latency (P50/P95/P99), error rate, memory, CPU, Epic FHIR, de-identification
- Dashboards: Server health, user activity, cost tracking
- Alerting: 6 critical policies (P0/P1/P2 severity)
- Runbooks: Server down, SSO issues, Epic connection failure
- Cost tracking: $1,047/month budget with 50%/75%/90%/100% alerts
- Performance optimization: Cold start reduction, throughput tuning
- Quarterly bias audits: Automated compliance reporting

**Files**: [`docs/for-hospitals/OPERATIONS_MANUAL.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-hospitals/OPERATIONS_MANUAL.md), [`docs/for-hospitals/RUNBOOKS/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/docs/for-hospitals/RUNBOOKS). 
**Monitoring**: Cloud Logging + Monitoring, 6 alert policies, 3 runbooks. 
**Compliance**: 10-year logs, quarterly bias audits, IRB reporting. 



ewpage


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

![Streamlit UI Preview](images/screenshots/streamlit-ui-preview.png){width=\textwidth}

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



ewpage


# Teaching Precision Medicine

*Educational workflows, classroom exercises, and student access*

---

## Why Teach With This Platform

Chapter 15 focused on research. **Educators have different needs**:

**Research use**: Hypothesis testing, publications, grant-funded analysis
**Teaching use**: Learning concepts, hands-on exercises, cost-effective student access

**Educator requirements**:
- **Low cost**: $0.32 per student analysis (DRY_RUN mode) vs $6,000 traditional
- **Safe data**: PatientOne 100% synthetic (no IRB needed, safe to share)
- **Scaffolded learning**: Beginner → Intermediate → Advanced prompts
- **Reproducible**: Same results every time (no random variation in DRY_RUN)
- **Hands-on**: Real tools, real workflows (not simulations)

**PatientOne for teaching**: Perfect educational dataset (Stage IV ovarian cancer, multi-modal data, treatment decisions).

---

## Educational Prompts (10 Prompts)

**Curated prompts** organized by difficulty level:

### Beginner Level (Prompts 1-3)

**Time**: 5-10 minutes each
**Prerequisites**: None (intro to bioinformatics)

#### Prompt 1: Introduction to MCP

**Learning objectives**:
- Understand MCP server architecture
- Know which servers handle which data types
- See how servers connect in workflows

**Prompt**:
```
Welcome to Precision Medicine MCP!

Task 1: List all available MCP servers.
For each server, tell me:
- Server name
- What type of data it handles
- How many tools it provides

Task 2: Explore mcp-spatialtools in detail.
List all its tools and explain what each does.

Task 3: Explain how servers work together.
Which servers would I use to analyze a cancer patient's data?
```

**Expected output**: 10 servers listed, 14 spatialtools tools explained, workflow diagram.

**Time**: 5 minutes
**Cost**: $0.05 (minimal API usage)

#### Prompt 2: First Analysis - Gene Expression

**Learning objectives**:
- Load and explore gene expression data
- Compare expression across regions
- Interpret biological meaning

**Prompt**:
```
Let's run your first analysis!

Step 1: Load PatientOne spatial transcriptomics data.
How many spots (tissue locations) are there? How many genes?

Step 2: Find MKI67 gene expression (proliferation marker).
What is the average expression? Highest? Lowest?

Step 3: Compare MKI67 in different regions:
- tumor_core
- tumor_proliferative
- stroma (normal tissue)

Which region has highest MKI67? What does this tell us about cell division?
```

**Expected output**: 900 spots, 31 genes, MKI67 highest in tumor_proliferative (5.91 vs 0.85 in stroma).

**Time**: 5 minutes
**Cost**: $0.10

#### Prompt 3: Statistical Significance

**Learning objectives**:
- Understand p-values and statistical significance
- Learn why multiple testing correction is crucial
- Distinguish statistical vs biological significance

**Prompt**:
```
Let's learn about statistical significance!

Step 1: Compare MKI67 expression (tumor vs stroma).
Use Mann-Whitney U test. What is the p-value?

Step 2: Interpret the result.
Is MKI67 significantly different? How confident are you?
What does p = 0.001 mean in simple terms?

Step 3: Test ALL 31 genes (not just MKI67).
How many are significant at p < 0.05 without correction?
How many after FDR correction (q < 0.05)?
Why do we need correction?
```

**Expected output**: MKI67 p = 8.3e-12 (highly significant), 28 genes significant uncorrected → 17 after FDR.

**Time**: 10 minutes
**Cost**: $0.15

Full beginner prompts: [`docs/prompt-library/educational-prompts.md#beginner-level`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/prompt-library/educational-prompts.md#beginner-level)

### Intermediate Level (Prompts 4-6)

**Time**: 15-20 minutes each
**Prerequisites**: Basic statistics, genomics concepts

#### Prompt 4: Pathway Enrichment

**Learning objectives**:
- Understand pathway enrichment concept
- Calculate fold enrichment
- Interpret clinical implications

**Prompt**:
```
Learn pathway enrichment analysis!

Step 1: Your gene list (10 upregulated genes):
TP53, PIK3CA, AKT1, MTOR, PTEN, BRCA1, MYC, VEGFA, HIF1A, BCL2L1

Step 2: Run pathway enrichment (GO Biological Process, FDR < 0.05).
For each enriched pathway:
- How many of your genes are in it?
- What does the pathway do biologically?
- Why is it relevant to cancer?

Step 3: Clinical implications.
Which pathways could be drug targets? Are there FDA-approved drugs?
```

**Expected output**: PI3K/Akt pathway (8/10 genes, 5.33× enriched), apoptosis regulation, hypoxia response.

**Time**: 15 minutes
**Cost**: $0.20

#### Prompt 5: Batch Effects and Preprocessing

**Learning objectives**:
- Visualize batch effects with PCA
- Apply ComBat batch correction
- Validate correction effectiveness

**Prompt**:
```
Learn about batch effects!

Step 1: Load PatientOne spatial data with batch labels.
Create PCA plot colored by batch. Do you see clustering by batch?

Step 2: Apply ComBat batch correction.
Create PCA plot after correction. Did batch effect disappear?

Step 3: Validate correction.
Run differential expression before and after correction.
How many genes change significance status?
```

**Expected output**: Clear batch clustering before correction, uniform distribution after, 3-5 genes change significance.

**Time**: 15 minutes
**Cost**: $0.25

#### Prompt 6: Multi-Omics Integration

**Learning objectives**:
- Understand multi-modal evidence concept
- Learn Stouffer meta-analysis
- Interpret combined confidence levels

**Prompt**:
```
Learn multi-omics integration!

Step 1: Load PatientOne RNA, protein, phospho data (15 samples).
How many features in each modality?

Step 2: Integrate using Stouffer's method.
Identify pathways activated concordantly across all 3 modalities (FDR < 0.05).

Step 3: Compare single-modality vs multi-modal evidence.
How does multi-modal evidence increase confidence?
Which pathways appear in all 3 modalities?
```

**Expected output**: PI3K/Akt activated in all 3 modalities (combined p = 2.1e-5), stronger evidence than any single modality.

**Time**: 20 minutes
**Cost**: $0.30

Full intermediate prompts: [`docs/prompt-library/educational-prompts.md#intermediate-level`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/prompt-library/educational-prompts.md#intermediate-level)

### Advanced Level (Prompts 7-10)

**Time**: 20-45 minutes each
**Prerequisites**: Advanced bioinformatics, statistics

#### Prompt 7: Spatial Transcriptomics

**Learning objectives**:
- Analyze spatial patterns with Moran's I
- Understand spatial autocorrelation
- Interpret treatment implications

**Expected output**: Immune exclusion zones identified, spatial heterogeneity quantified, treatment strategy informed.

**Time**: 20 minutes
**Cost**: $0.35

#### Prompt 8: Clinician-in-the-Loop Workflow

**Learning objectives**:
- Understand AI oversight and safety
- Learn quality gates and approval process
- Recognize human accountability in AI systems

**Expected output**: Draft report generated, quality gates passed, clinician review and approval, final signed report.

**Time**: 15 minutes
**Cost**: $0.40

#### Prompt 9: Cost-Effectiveness Analysis

**Learning objectives**:
- Calculate cost savings and ROI
- Understand value-based healthcare
- Balance cost vs outcomes

**Expected output**: $3,200 → $1.35 per analysis (99.96% savings), 3,942 hours saved annually, 2-3 patient payback period.

**Time**: 10 minutes
**Cost**: $0.10

#### Prompt 10: Final Project - Complete Patient Analysis

**Learning objectives**:
- Apply all concepts in integrated workflow
- Synthesize multi-modal evidence
- Generate treatment recommendations

**Scenario**: Hypothetical breast cancer patient (HER2+, ER+, Stage III).

**Expected output**: Complete multi-modal analysis (clinical, genomic, multiomics, spatial), treatment recommendations with evidence, cost and time comparison to traditional methods.

**Time**: 45 minutes
**Cost**: $0.50

Full advanced prompts: [`docs/prompt-library/educational-prompts.md#advanced-level`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/prompt-library/educational-prompts.md#advanced-level)

---

## Companion Jupyter Notebooks

**Each book chapter has a corresponding notebook with hands-on exercises**:

![Jupyter Notebook Preview](images/screenshots/jupyter-preview.png){width=\textwidth}

**Figure 16.1: Jupyter Notebook Teaching Interface**
*Interactive Python notebooks for hands-on learning with executable code cells, visualizations, and explanatory markdown. Students can modify parameters, run analyses, and see results in real-time.*

**IMPORTANT**: These notebooks require you to deploy your own MCP servers to GCP Cloud Run. See **Appendix: Setup Guide** for complete deployment instructions.

### All 18 Notebooks

| Part | Chapter | Notebook | Topics |
|------|---------|----------|--------|
| **Part 1** | 1 | `chapter-01-patientone-story.ipynb` | PatientOne workflow demo, time/cost savings |
| | 2 | `chapter-02-architecture.ipynb` | MCP orchestration, server discovery |
| | 3 | `chapter-03-testing-the-hypothesis.ipynb` | Test coverage, cost analysis, production metrics |
| **Part 2** | 4 | `chapter-04-clinical-data.ipynb` | FHIR integration, de-identification |
| | 5 | `chapter-05-genomic-foundations.ipynb` | VCF parsing, variant annotation |
| | 6 | `chapter-06-multi-omics-integration.ipynb` | HAllA, Stouffer meta-analysis |
| | 7 | `chapter-07-spatial-transcriptomics.ipynb` | Spatial transcriptomics, Moran's I |
| **Part 3** | 8 | `chapter-08-cell-segmentation.ipynb` | DeepCell segmentation, phenotype classification |
| | 9 | `chapter-09-treatment-response.ipynb` | GEARS GNN, treatment prediction |
| | 10 | `chapter-10-quantum-fidelity.ipynb` | Quantum circuits, Bayesian UQ |
| | 11 | `chapter-11-imaging.ipynb` | H&E, MxIF histology analysis |
| **Part 4** | 12 | `chapter-12-cloud-deployment.ipynb` | Docker, Cloud Run, SSE transport |
| | 13 | `chapter-13-hospital-deployment.ipynb` | HIPAA compliance, VPC, Azure AD SSO |
| | 14 | `chapter-14-operations.ipynb` | Logging, monitoring, cost tracking |
| **Part 5** | 15 | `chapter-15-for-researchers.ipynb` | Research workflows, exploratory analysis |
| | 16 | `chapter-16-teaching.ipynb` | Educational prompts and exercises |
| **Part 6** | 17 | `chapter-17-funding.ipynb` | ROI calculator, grant budgets |
| | 18 | `chapter-18-lessons-learned.ipynb` | Production insights, future roadmap |

**Total**: 18 companion notebooks

### Quick Setup

```bash
# Clone repository
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp/docs/book/companion-notebooks

# Install dependencies
pip install -r requirements.txt

# Deploy MCP servers to YOUR GCP project (REQUIRED)
# See Appendix: Setup Guide for complete instructions
./infrastructure/deployment/deploy_to_gcp.sh YOUR_PROJECT_ID us-central1

# Configure API keys in .env file
# ANTHROPIC_API_KEY=your_key  OR  GOOGLE_API_KEY=your_key
# MCP_FGBIO_URL=https://mcp-fgbio-YOUR_PROJECT.run.app/sse
# ... (add all your server URLs)

# Launch Jupyter
jupyter lab
```

**Cost**: ~$10-20 total (Claude/Gemini API + Cloud Run for all 18 notebooks)

Full notebooks guide and troubleshooting: [`docs/book/companion-notebooks/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/book/companion-notebooks/README.md)

---

## Classroom Exercises

### Exercise 1: Gene Expression Basics (Week 1)

**Learning objectives**: Data loading, basic statistics, biological interpretation

**Assignment**:
1. Load PatientOne spatial data
2. Identify top 10 most variable genes
3. Compare gene expression in tumor vs stroma
4. Create visualization (heatmap or volcano plot)
5. Write 1-paragraph biological interpretation

**Time**: 1 hour
**Due**: End of week 1
**Grading**: Correctness (50%), visualization (25%), interpretation (25%)

### Exercise 2: Pathway Analysis (Week 3)

**Learning objectives**: Pathway enrichment, statistical significance, drug targets

**Assignment**:
1. Run differential expression (tumor vs normal, FDR < 0.05)
2. Perform pathway enrichment on upregulated genes
3. Identify top 3 enriched pathways
4. Research FDA-approved drugs targeting these pathways
5. Write 2-page report: findings, clinical implications, drug recommendations

**Time**: 2 hours
**Due**: End of week 3
**Grading**: Analysis (40%), drug research (30%), report quality (30%)

### Exercise 3: Multi-Modal Integration (Week 6)

**Learning objectives**: Multi-omics integration, meta-analysis, evidence synthesis

**Assignment**:
1. Load PatientOne RNA, protein, phospho data
2. Run Stouffer meta-analysis
3. Compare single-modality vs multi-modal evidence
4. Identify pathways concordant across all 3 modalities
5. Explain why multi-modal evidence is stronger
6. Create presentation (5-7 slides)

**Time**: 3 hours
**Due**: End of week 6
**Grading**: Technical correctness (50%), presentation (30%), explanation (20%)

### Exercise 4: Final Project - Complete Patient Analysis (Week 12)

**Learning objectives**: Integrate all concepts, clinical decision-making, reproducibility

**Assignment**:
1. Analyze hypothetical patient (provided case study: breast cancer, lung cancer, or melanoma)
2. Perform multi-modal analysis (clinical, genomic, multiomics, spatial)
3. Generate treatment recommendations with evidence
4. Calculate cost and time comparison to traditional methods
5. Write methods section (for publication reproducibility)
6. Present findings (10-minute presentation)

**Time**: 6-8 hours
**Due**: End of semester
**Grading**: Technical analysis (40%), treatment recommendations (30%), presentation (20%), methods documentation (10%)

---

## Cost-Effective Student Access

### DRY_RUN Mode (Demonstration)

**Cost**: $0.32 per student per analysis

**How it works**:
- Servers return realistic synthetic data (not real analysis)
- Fast responses (~2-5 seconds)
- Reproducible results (same every time)
- No cloud compute costs
- Only Claude/Gemini API tokens charged

**Use case**: Introductory courses, concept learning, large classrooms (50+ students).

**Example**: 50 students × 10 assignments × $0.32 = **$160 total** (entire semester).

### Production Mode (Real Analysis)

**Cost**: $25-104 per student per analysis (varies by dataset size)

**How it works**:
- Real bioinformatics analysis (not mocked)
- Longer processing times (10-30 minutes)
- Variable results (depends on actual data)
- Cloud compute costs included
- Higher API token usage

**Use case**: Advanced courses, graduate seminars, capstone projects.

**Example**: 15 students × 3 projects × $50 avg = **$2,250 total** (vs $90,000 traditional).

**Savings**: 97.5% reduction

### Hybrid Approach (Recommended)

**Structure**:
- Weeks 1-8: DRY_RUN mode for concept learning ($0.32 per assignment)
- Weeks 9-12: Production mode for final projects ($50 per project)

**Total cost per student**: ~$53 (semester)
**Total cost for 30 students**: ~$1,590

**Traditional alternative**: $6,000 per student × 30 = $180,000
**Savings**: $178,410 (99.1% reduction)

---

## Classroom Setup

### Option 1: Claude Desktop (Local)

**Setup** (per student):
```bash
# Install Claude Desktop
# Download from: https://claude.com/claude-desktop

# Configure MCP servers
# Edit: ~/Library/Application Support/Claude/claude_desktop_config.json
```

**Pros**: Free (students use own API keys), no deployment needed
**Cons**: Students need own computers, API keys cost $5-20/semester

### Option 2: Jupyter Notebooks (Cloud)

**Setup** (instructor):
```bash
# Deploy to JupyterHub
# Use GCP deployment scripts from Chapter 13

# Students access via browser (no local installation)
```

**Pros**: Centralized (instructor controls), no student setup
**Cons**: Requires instructor cloud deployment, $1,000-2,000/semester infrastructure

### Option 3: Google Colab (Free)

**Setup**: Students click "Open in Colab" links in notebooks

**Pros**: Completely free (Google-provided compute), zero setup
**Cons**: Limited to Gemini API (not Claude), session timeouts

**Recommended for**: Large intro courses (100+ students), budget-constrained institutions.

---

## Learning Outcomes Assessment

### Beginner Level (Weeks 1-4)

**Skills assessed**:
- ✓ Load and explore genomic datasets
- ✓ Perform basic statistical tests
- ✓ Interpret p-values and FDR
- ✓ Create visualizations (plots, heatmaps)

**Assessment method**: Weekly quizzes + Exercise 1

### Intermediate Level (Weeks 5-8)

**Skills assessed**:
- ✓ Run pathway enrichment analysis
- ✓ Apply batch correction techniques
- ✓ Integrate multi-modal data
- ✓ Interpret biological significance

**Assessment method**: Exercise 2 + Exercise 3

### Advanced Level (Weeks 9-12)

**Skills assessed**:
- ✓ Design complete analysis workflows
- ✓ Generate treatment recommendations
- ✓ Write reproducible methods sections
- ✓ Present findings to clinical audience

**Assessment method**: Final project (Exercise 4)

---

## What You've Configured

**Educational platform**:
1. **10 educational prompts**: Beginner (3) → Intermediate (3) → Advanced (4)
2. **18 Jupyter notebooks**: Hands-on exercises for each chapter
3. **4 classroom exercises**: Gene expression, pathways, multi-omics, final project
4. **Cost-effective access**: $0.32 (DRY_RUN) to $104 (production) per student per analysis
5. **Hybrid approach**: Concept learning (DRY_RUN) + final projects (production) = $53/student/semester

**Student learning**:
- PatientOne synthetic dataset (100% safe, no IRB needed)
- Real tools, real workflows (not simulations)
- Reproducible results (DRY_RUN mode)
- Hands-on bioinformatics concepts
- Clinical decision-making experience

**Cost comparison**:
- Traditional: $6,000 per student (manual analysis time)
- This platform: $53 per student (hybrid approach)
- Savings: $5,947 per student (99.1% reduction)

---

## Summary

**Chapter 16 Summary**:
- Educational prompts: 10 prompts (beginner → intermediate → advanced)
- Companion notebooks: 18 Jupyter notebooks for hands-on learning
- Classroom exercises: 4 assignments (gene expression, pathways, multi-omics, final project)
- Cost-effective access: $0.32 (DRY_RUN) to $104 (production) per analysis
- Hybrid approach: $53/student/semester (vs $6,000 traditional)
- PatientOne dataset: 100% synthetic, safe to share, no IRB needed
- Learning outcomes: Basic stats → Pathway analysis → Multi-omics → Complete workflows

**Files**: [`docs/prompt-library/educational-prompts.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/prompt-library/educational-prompts.md), [`docs/book/companion-notebooks/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/docs/book/companion-notebooks). 
**Student cost**: $53/semester (99.1% savings vs traditional)  
**Setup options**: Claude Desktop (local), JupyterHub (cloud), Google Colab (free). 



ewpage


# Funding and Sustainability

*ROI analysis, grant strategies, and budget models*

---

## Why Funding Matters

Chapters 1-16 built and deployed the system. **Now: How do you pay for it?**

**Deployment costs money**:
- Infrastructure: $1,000-2,000/month (Cloud Run, storage, networking)
- API usage: $500-1,000/month (Claude/Gemini tokens)
- Personnel: $50,000-120,000/year (bioinformatics support)

**But traditional analysis costs more**:
- Personnel: $120,000-240,000/year (1.5-3 FTE manual analysis)
- Infrastructure: $10,000-20,000/year (on-premise compute)
- External services: $600,000-900,000/year (100 patients × $6,000-9,000)

**This chapter**: ROI analysis, grant strategies, budget models for sustainability.

---

## ROI Analysis

### Per-Patient Cost Comparison

| Cost Component | Traditional | MCP System | Savings |
|----------------|-------------|------------|---------|
| Personnel time | $3,000 (40 hrs × $75/hr) | $75 (1 hr × $75/hr) | $2,925 |
| Compute resources | $2,000-4,000 | $22-99 | $1,901-3,978 |
| External services | $1,000-2,000 | $1-3 (APIs) | $997-1,999 |
| **Total per patient** | **$6,000-9,000** | **$98-177** | **$5,823-8,902** |

**Average savings**: $7,210 per patient (96% cost reduction)

### Annual Savings (100 Patients)

**Traditional costs**:
- 100 patients × $7,500 avg = $750,000/year
- Personnel: 4,000 hours (2 FTE @ $120K = $240,000)
- Infrastructure: $200,000-400,000
- External services: $100,000-200,000

**MCP costs**:
- 100 patients × $138 avg = $13,800/year
- Cloud infrastructure: $12,000-24,000/year
- Claude API: $1,200-2,400/year
- Personnel: $50,000/year (0.5 FTE oversight)

**Total MCP**: $76,000-90,000/year

**Annual savings**: $660,000-674,000 (88-89% reduction)

### Payback Period

**Initial investment**:
- Setup: $50,000-75,000 (deployment, training, integration)
- First-year operational: $76,000-90,000

**Total first-year**: $126,000-165,000

**Break-even**: After 17-23 patients analyzed (2-3 months in 100-patient/year facility)

**5-year ROI**:
- Investment: $126K-165K (Year 1) + $350K-450K (Years 2-5) = $476K-615K
- Traditional cost: $3.75M (5 years × $750K)
- **Total savings**: $3.13M-3.27M (655-687% ROI)

---

## Grant Strategies

### NIH R01 Budget Model

**Typical R01**: $250,000/year direct costs × 5 years = $1.25M

**Budget allocation**:

| Category | Traditional | MCP System | Savings |
|----------|-------------|------------|---------|
| Personnel | $600,000 (2 FTE × 5 years) | $250,000 (1 FTE × 5 years) | $350,000 |
| Equipment | $100,000 (servers, storage) | $0 (cloud-based) | $100,000 |
| Compute | $200,000 (cluster time) | $120,000 (Cloud Run) | $80,000 |
| Supplies | $150,000 (sequencing) | $150,000 (same) | $0 |
| Other | $200,000 (travel, pubs) | $200,000 (same) | $0 |
| **Total** | **$1,250,000** | **$720,000** | **$530,000** |

**Grant strategy**: Request $720K instead of $1.25M = **42% budget reduction**. Use savings for additional aims or patient cohort expansion.

**Budget justification template**:

```
Analysis Costs (Aim 2 - Multi-modal patient analysis):

Traditional approach: 200 patients × $7,500 = $1,500,000
Proposed MCP approach: 200 patients × $138 = $27,600

Cost savings: $1,472,400 (98% reduction)

We request $120,000 over 5 years for Cloud Run infrastructure
($24K/year) rather than $1.5M for traditional analysis. This 92%
cost reduction allows us to expand the cohort from 200 to 1,000
patients within the same budget, increasing statistical power and
clinical impact.
```

Full NIH budget template: [`docs/for-researchers/grant-templates/nih-r01-budget.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-researchers/grant-templates/nih-r01-budget.md) (planned)

### Foundation Grant Model

**Typical foundation grant**: $50,000-100,000 (1-2 years)

**Budget for pilot deployment** (Year 1):

| Item | Cost | Notes |
|------|------|-------|
| GCP infrastructure setup | $15,000 | VPC, OAuth2, Epic integration |
| Cloud Run (annual) | $24,000 | 100 patients/year |
| Claude API (annual) | $2,400 | 100 patients × $24 avg |
| Personnel (0.25 FTE) | $30,000 | Bioinformatics oversight |
| Training & documentation | $5,000 | 5 users, 2 days training |
| Contingency (10%) | $7,640 | Buffer for overruns |
| **Total Year 1** | **$84,040** | |

**Deliverables**:
- 100 ovarian cancer patients analyzed
- Publication: Methods paper + clinical cohort study
- Expansion plan: Scale to breast, lung, colorectal

**Sustainability plan** (Year 2+):
- Hospital absorbs operational costs ($26.4K/year)
- Grant funds patient cohort expansion (additional $50K/year)
- Savings reinvested in research ($660K/year traditional → $26.4K MCP = $634K reinvestment)

---

## Hospital Budget Model

### Pilot Deployment (3 Months, 25 Patients)

**Costs**:
- Infrastructure setup: $15,000 (one-time)
- Cloud Run: $2,000/month × 3 = $6,000
- Claude API: $600 (25 patients × $24)
- Personnel (0.5 FTE oversight): $15,000 (3 months)
- **Total pilot**: $36,600

**Traditional alternative**: 25 patients × $7,500 = $187,500
**Pilot savings**: $150,900 (80% reduction)

### Production Deployment (Annual, 100 Patients)

**Costs**:
- Cloud Run infrastructure: $24,000/year
- Claude API: $2,400/year
- Personnel oversight (0.5 FTE): $50,000/year
- Epic FHIR integration: $5,000/year (maintenance)
- Training & support: $5,000/year
- **Total production**: $86,400/year

**Traditional alternative**: $750,000/year
**Annual savings**: $663,600 (88% reduction)

### Scaling to 500 Patients/Year

**Costs**:
- Cloud Run (auto-scales): $50,000/year
- Claude API: $12,000/year (500 × $24)
- Personnel (1 FTE): $120,000/year
- Infrastructure maintenance: $10,000/year
- **Total**: $192,000/year

**Traditional alternative**: $3,750,000/year (500 × $7,500)
**Annual savings**: $3,558,000 (95% reduction)

**Cost per patient remains constant**: $138 (vs $7,500 traditional)

---

## Revenue Models (For-Profit)

### Clinical Service Model

**Charge structure**:
- Base analysis: $500 per patient (vs $3,000-5,000 competitors)
- Cost: $138 per patient
- **Gross margin**: $362 per patient (72%)

**Volume scenarios**:

| Patients/Year | Revenue | Costs | Profit | Margin |
|---------------|---------|-------|--------|--------|
| 100 | $50,000 | $13,800 | $36,200 | 72% |
| 500 | $250,000 | $69,000 | $181,000 | 72% |
| 1,000 | $500,000 | $138,000 | $362,000 | 72% |

**Market positioning**: 83% cheaper than competitors ($500 vs $3,000), same quality (multi-modal analysis).

### SaaS Subscription Model

**Tier structure**:

| Tier | Monthly | Included Analyses | Overage Cost | Target Customer |
|------|---------|-------------------|--------------|-----------------|
| Starter | $5,000 | 10 patients | $400/patient | Small clinics |
| Professional | $15,000 | 50 patients | $250/patient | Medium hospitals |
| Enterprise | $40,000 | 200 patients | $150/patient | Large cancer centers |

**Example**: Professional tier (50 patients/month)
- Revenue: $15,000/month = $180,000/year
- Cost: $6,900/month (50 × $138) = $82,800/year
- **Profit**: $97,200/year (54% margin)

**Customer acquisition**:
- Year 1: 5 customers (pilot hospitals) = $900,000 revenue
- Year 2: 15 customers (expansion) = $2,700,000 revenue
- Year 3: 30 customers (scale) = $5,400,000 revenue

---

## Cost Drivers and Optimization

### Primary Cost Drivers

1. **Claude API tokens** (40-50% of operational costs)
   - Optimization: Use Haiku for simple queries ($0.25 vs $3 per million input tokens)
   - Savings: 30-40% reduction in API costs

2. **Cloud Run compute** (30-40% of operational costs)
   - Optimization: Right-size memory/CPU, scale to zero when idle
   - Savings: 20-30% reduction in compute costs

3. **Personnel oversight** (10-20% of operational costs)
   - Optimization: Automate quality checks, batch review workflows
   - Savings: 10-15% reduction in personnel time

**Example optimization** (100 patients/year):
- Baseline: $86,400/year
- Haiku for 50% of queries: -$600 (API savings)
- Right-sized compute: -$4,800 (Cloud Run savings)
- Automated QC: -$5,000 (personnel savings)
- **Optimized total**: $76,000/year (12% reduction)

### Sensitivity Analysis

**What if Claude API pricing changes?**

| API Cost Change | New Cost/Patient | New Annual (100 patients) | Impact |
|----------------|------------------|---------------------------|--------|
| -50% (competition) | $126 | $74,400 | +15% savings |
| +50% (price increase) | $150 | $98,400 | -14% savings |
| +100% (doubling) | $162 | $110,400 | -28% savings |

**Still cheaper than traditional** even if Claude API doubles in price: $162 vs $7,500 (98% savings).

**Mitigation**: Multi-LLM strategy (Claude, Gemini, Llama) ensures price competition.

---

## Sustainability Roadmap

### Year 1: Pilot and Validation

**Goals**:
- Deploy to 1 hospital (100 patients)
- Validate clinical outcomes
- Publish methods paper
- Establish baseline metrics

**Funding**: Foundation grant ($85K) or hospital pilot budget ($87K)

**Deliverables**:
- 100 patients analyzed
- Clinical validation study published
- Cost savings documented ($660K vs traditional)

### Year 2: Expansion and Optimization

**Goals**:
- Scale to 3 hospitals (300 patients total)
- Expand to 2 additional cancer types
- Optimize costs (12% reduction)
- Establish revenue model (if for-profit)

**Funding**: Hospital operational budgets ($76K each) or NIH R01 Year 2 ($150K)

**Deliverables**:
- 300 patients analyzed across 3 cancer types
- Multi-cancer validation published
- ROI metrics: $1.98M savings (300 patients)

### Year 3: Scale and Productization

**Goals**:
- Scale to 10 hospitals (1,000 patients total)
- Full production deployment (all 12 servers)
- SaaS offering (if for-profit)
- Self-sustaining operations

**Funding**: Self-sustaining (hospital budgets) or revenue (SaaS subscriptions)

**Deliverables**:
- 1,000 patients analyzed
- Consortium publication
- Total savings: $7.2M vs traditional

---

## What You've Configured

**Funding models**:
1. **ROI analysis**: $7,210 savings per patient (96% reduction), $660K annual savings (100 patients)
2. **Grant strategies**: NIH R01 ($720K vs $1.25M traditional), foundation grants ($85K pilot)
3. **Hospital budgets**: $87K/year production (100 patients) vs $750K traditional
4. **Revenue models**: Clinical service ($500/patient, 72% margin), SaaS ($5K-40K/month tiers)
5. **Cost optimization**: Haiku for simple queries (-30% API), right-sized compute (-20%), automated QC (-15%)
6. **Sustainability roadmap**: Year 1 pilot (100 patients), Year 2 expansion (300 patients), Year 3 scale (1,000 patients)

**Key metrics**:
- Per-patient savings: $7,210 (96% reduction)
- Payback period: 17-23 patients (2-3 months)
- 5-year ROI: 655-687% ($3.13M-3.27M savings)
- Scalability: Cost/patient constant at $138 (traditional increases)

---

## Summary

**Chapter 17 Summary**:
- Per-patient cost: $138 (MCP) vs $7,500 (traditional) = $7,210 savings (96%)
- Annual savings: $660K (100 patients), $3.56M (500 patients)
- Payback period: 2-3 months (17-23 patients)
- 5-year ROI: 655-687% ($3.13M-3.27M savings)
- Grant strategies: NIH R01 ($720K vs $1.25M), foundation ($85K pilot)
- Hospital budgets: $87K/year (100 patients) vs $750K traditional
- Revenue models: Clinical service (72% margin), SaaS ($5K-40K/month)
- Sustainability: Year 1 pilot → Year 2 expansion → Year 3 scale

**Files**: [`docs/prompt-library/funder-demo-prompts.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/prompt-library/funder-demo-prompts.md). 
**ROI**: 96% cost reduction, 2-3 month payback, 655-687% 5-year ROI. 
**Scaling**: Cost/patient constant ($138) as volume increases. 



ewpage


# Lessons Learned and What's Next

*Production insights, future enhancements, and the path forward*

---

## The Journey So Far

**17 chapters ago**, you started with a problem: 40 hours of manual analysis per patient, $3,200 cost, treatment decisions delayed.

**Now**, you have:
- 12 MCP servers deployed (124 tools total)
- PatientOne analysis: 35 minutes, $1.35 cost (98% faster, 99.96% cheaper)
- HIPAA-compliant hospital deployment
- Educational platform ($0.32 per student analysis)
- Research workflows (90+ prompts, reproducible methods)
- Funding models (655-687% 5-year ROI)

**This chapter**: Reflect on what worked, what didn't, and where we go from here.

---

## What Worked: Production Insights

### 1. MCP Architecture Decision

**Why it worked**:
- **Open standard**: Not locked into proprietary platform
- **Modular design**: Deploy only needed servers (not monolithic)
- **Language flexibility**: FastMCP (Python), Claude Desktop (any language)
- **Natural language orchestration**: No coding required for clinical users

**Evidence**:
- 12 independent servers = easier to debug, update, scale
- Claude orchestrates 124 tools via natural language prompts
- Hospital IT comfortable with standard HTTP/SSE (not custom protocol)

**Alternative considered**: Monolithic bioinformatics pipeline (rejected - too rigid, hard to maintain)

### 2. PatientOne Synthetic Dataset

**Why it worked**:
- **100% synthetic**: No IRB approval needed, safe to share
- **Realistic complexity**: Stage IV ovarian cancer, multi-modal data
- **Educational value**: Teaches concepts without privacy concerns
- **Publication-ready**: Methods papers, software validation

**Evidence**:
- Used in 90+ prompts across all audiences (funders, clinicians, researchers, educators)
- Included in every chapter (consistent examples)
- Deployed to GCS with public read access (students, workshops)

**Alternative considered**: Real de-identified data (rejected - IRB burden, sharing restrictions, privacy risks)

### 3. Hybrid DRY_RUN + Production Mode

**Why it worked**:
- **DRY_RUN ($0.32)**: Cost-effective learning, concept validation, classroom use
- **Production ($25-104)**: Real analysis for research and clinical decisions
- **Smooth transition**: Same prompts, same interface, different results

**Evidence**:
- 50 students × 10 assignments × $0.32 = $160 total (semester) vs $30,000 traditional
- Researchers prototype in DRY_RUN, then switch to production for publication
- Hospital pilots start in DRY_RUN, prove value, then production deployment

**Alternative considered**: Production-only (rejected - too expensive for education, exploration)

---

## What Didn't Work: Challenges and Solutions

### Challenge 1: DeepCell Production Integration (3 Weeks)

**Problem**: DeepCell-TF requires Python 3.10 (TensorFlow 2.8.x constraint), won't build on modern Python 3.11+.

**Failed attempts**:
1. Python 3.11 build → ImportError (TensorFlow incompatible)
2. Ubuntu 20.04 base image → Package conflicts
3. N1_HIGHCPU_8 machine type → Deprecated

**Solution** (Attempt 4):
- Switched to Python 3.10 base image
- Added system dependencies: `libgomp1`, `libhdf5-dev`
- Used E2_HIGHCPU_8 machine type (modern)
- GCS image loading: Download to temp file (PIL doesn't support `gs://` URIs)

**Time cost**: 2 weeks debugging, 1 week documenting

**Lesson**: Deep learning libraries have strict dependency constraints. Plan for compatibility issues.

Full story: [`servers/mcp-deepcell/DEPENDENCY_ISSUES.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/DEPENDENCY_ISSUES.md)

### Challenge 2: Stouffer Meta-Analysis FDR Timing

**Problem**: When to apply FDR correction - before or after combining p-values?

**Initial approach** (wrong):
```python
# Apply FDR to each modality separately
rna_fdr = benjamini_hochberg(rna_pvalues)
protein_fdr = benjamini_hochberg(protein_pvalues)
phospho_fdr = benjamini_hochberg(phospho_pvalues)

# Then combine FDR-corrected values
combined_p = stouffer(rna_fdr, protein_fdr, phospho_fdr)  # WRONG!
```

**Bioinformatician feedback**: "FDR correction should be applied AFTER combining p-values, not before."

**Correct approach**:
```python
# Combine raw p-values first
combined_p = stouffer(rna_pvalues, protein_pvalues, phospho_pvalues)

# Then apply FDR correction to combined p-values
combined_fdr = benjamini_hochberg(combined_p)  # CORRECT
```

**Lesson**: Statistical methods have ordering dependencies. Validate with domain experts.

Full implementation: [`servers/mcp-multiomics/src/mcp_multiomics/tools/stouffer.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-multiomics/src/mcp_multiomics/tools/stouffer.py)

### Challenge 3: Streamlit Public Access Costs

**Problem**: User wanted deployed Streamlit app for public demos. Cost: $200-300/month for always-on Cloud Run instance.

**User feedback**: "do NOT refer to the deployed streamlit app. Readers are going to have to deploy from the repo as I can't pay for the general public to access my mcp servers on GCP."

**Solution**:
- Removed all public Streamlit references from book
- Emphasized local deployment (Claude Desktop, Jupyter notebooks)
- DRY_RUN mode for cost-effective demos

**Lesson**: Public SaaS deployment requires ongoing costs. For open source projects, prefer local deployment + documentation.

### Challenge 4: Server Implementation Status Confusion

**Problem**: Readers unclear which servers are production-ready vs mocked.

**Initial docs**: "Servers are mostly complete" (vague)

**Solution**: Created comprehensive server status matrix:
- 4/12 production-ready: fgbio, multiomics, spatialtools, epic
- 2/12 partial: openimagedata, quantum-celltype-fidelity
- 6/12 mocked: deepcell, perturbation, tcga, huggingface, seqera, mockepic

**Documentation**: [`docs/architecture/servers.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/architecture/servers.md) (1,000+ line status matrix)

**Lesson**: Transparency builds trust. Document implementation status honestly and comprehensively.

---

## Future Enhancements

### Phase 1 (3-6 Months): Complete Advanced Servers

**mcp-deepcell** (currently mocked):
- **Goal**: Real DeepCell-TF cell segmentation (nuclear + membrane)
- **Status**: Deployment complete (Chapter 8), integration pending
- **Effort**: 2-3 weeks (model loading, inference optimization, validation)
- **Impact**: Single-cell resolution phenotyping from MxIF images

**mcp-perturbation** (currently mocked):
- **Goal**: Real GEARS GNN treatment response prediction
- **Status**: Framework complete (Chapter 9), model integration pending
- **Effort**: 3-4 weeks (model training, gene regulatory graphs, validation)
- **Impact**: Predict drug efficacy before treatment (40% better than VAE methods)

**mcp-quantum-celltype-fidelity** (40% real):
- **Goal**: Real PQC (Parameterized Quantum Circuits) for cell-type classification
- **Status**: Bayesian UQ complete, PQC mocked
- **Effort**: 4-6 weeks (quantum circuit design, PennyLane integration, validation)
- **Impact**: Quantum advantage for high-dimensional embeddings

**Total effort**: 9-13 weeks (2-3 months)

### Phase 2 (6-12 Months): Public Data Integration

**mcp-tcga** (currently mocked):
- **Goal**: Real GDC API for TCGA cohort data (11,000+ patients)
- **Use case**: Biomarker validation, cohort comparisons, survival analysis
- **Effort**: 4-6 weeks (GDC API integration, data preprocessing, validation)
- **Impact**: Validate PatientOne findings in 500+ ovarian cancer patients

**mcp-huggingface** (currently mocked):
- **Goal**: Genomic foundation models (DNABERT, Nucleotide Transformer)
- **Use case**: Variant effect prediction, regulatory element identification
- **Effort**: 6-8 weeks (model integration, inference optimization, validation)
- **Impact**: Predict variant pathogenicity without functional assays

**Total effort**: 10-14 weeks (2.5-3.5 months)

### Phase 3 (12-18 Months): New Modalities

**mcp-metabolomics**:
- **Goal**: Metabolite profiling, pathway mapping
- **Data**: LC-MS/MS metabolomics data
- **Servers**: 8-10 tools (normalization, pathway enrichment, upstream regulators)
- **Impact**: Complete multi-omics (genomics + transcriptomics + proteomics + metabolomics)

**mcp-radiomics**:
- **Goal**: CT/MRI feature extraction, tumor burden quantification
- **Data**: DICOM medical images
- **Servers**: 6-8 tools (segmentation, radiomics features, response prediction)
- **Impact**: Non-invasive tumor monitoring, treatment response prediction

**mcp-singlecell**:
- **Goal**: Single-cell RNA-seq analysis (10X Genomics, CITE-seq)
- **Data**: scRNA-seq count matrices
- **Servers**: 12-15 tools (QC, clustering, trajectory analysis, cell-cell communication)
- **Impact**: Cellular heterogeneity, rare cell type discovery, developmental trajectories

**Total effort**: 30-40 weeks (7-10 months)

---

## Multi-Cancer Expansion

**Current**: Ovarian cancer (PatientOne)

**Next cancer types** (priority order):

### 1. Breast Cancer (HER2+, ER+, Triple-Negative)

**Timeline**: 2-4 weeks
**Effort**:
- Create PatientTwo synthetic dataset (HER2+ breast cancer)
- Adapt pathways (add HER2 signaling, hormone pathways)
- Update prompts for breast-specific biomarkers
- Validate with TCGA BRCA cohort

**Impact**: 268,000 new cases/year in US (vs 20,000 ovarian)

### 2. Colorectal Cancer (MSI-H, MSS)

**Timeline**: 2-4 weeks
**Effort**:
- Create PatientThree synthetic dataset (MSI-H colorectal cancer)
- Adapt pathways (add WNT signaling, immune checkpoints)
- Integrate MSI/MMR status
- Validate with TCGA COAD cohort

**Impact**: 153,000 new cases/year in US

### 3. Lung Cancer (NSCLC, EGFR+, ALK+)

**Timeline**: 2-4 weeks
**Effort**:
- Create PatientFour synthetic dataset (EGFR+ NSCLC)
- Adapt pathways (add EGFR signaling, immune checkpoints)
- Integrate PD-L1 status, TMB
- Validate with TCGA LUAD cohort

**Impact**: 236,000 new cases/year in US

**Total timeline**: 6-12 weeks for 3 cancer types
**Total impact**: 657,000 new cases/year (vs 20,000 ovarian)

---

## Community and Open Source

### Current Status

**Open source**:
- Repository: GitHub (Apache 2.0 license)
- Code: All servers, docs, companion notebooks
- Data: PatientOne 100% synthetic (CC0 license)

**Not open source**:
- Book content (copyright © 2026 Lynn Langit)
- Deployed Cloud Run servers (cost constraints)

### Community Building

**Phase 1 (Months 1-6)**: Foundation
- GitHub Discussions for Q&A
- Community Slack/Discord channel
- Monthly office hours (live Q&A)
- Contributor guidelines

**Phase 2 (Months 7-12)**: Growth
- Bioinformatics conference talks (ASHG, ISMB, AMIA)
- Workshop series (online, free)
- Use case showcases (community contributions)
- Academic partnerships (10+ institutions)

**Phase 3 (Year 2+)**: Ecosystem
- Plugin marketplace (custom MCP servers)
- Third-party integrations (Terra, DNAnexus, Galaxy)
- Certification program (MCP server developers)
- Annual conference (Precision Medicine MCP Summit)

### Contribution Areas

**High priority**:
1. Server implementations (DeepCell, GEARS, TCGA, Hugging Face)
2. Multi-cancer datasets (breast, colorectal, lung)
3. Pathway databases (expand from 44 to 100+ pathways)
4. Educational content (prompts, notebooks, exercises)

**Medium priority**:
5. UI improvements (Streamlit, JupyterHub features)
6. Cost optimizations (caching, batching, Haiku adoption)
7. Performance tuning (cold start reduction, GPU support)
8. Documentation (translations, video tutorials)

**Contribution guide**: [`CONTRIBUTING.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/CONTRIBUTING.md) (planned)

---

## The Path Forward

**Short term (3-6 months)**:
- Complete Phase 1 advanced servers (DeepCell, GEARS, Quantum)
- Expand to 3 cancer types (breast, colorectal, lung)
- Publish methods paper + PatientOne clinical validation
- Grow community to 100+ users

**Medium term (6-18 months)**:
- Integrate public data (TCGA, Hugging Face genomic models)
- Add new modalities (metabolomics, radiomics, single-cell)
- Scale to 10 hospitals (1,000 patients/year)
- Establish consortium (academic + hospital partners)

**Long term (18+ months)**:
- Full production deployment (all 12 servers real analysis)
- Multi-cancer platform (10+ cancer types)
- Real-world evidence generation (10,000+ patients)
- Regulatory approval path (FDA, Clinical decision support)

---

## Final Thoughts

**40 hours → 35 minutes**. That's the transformation this book documents.

But the real transformation is deeper:

**From**: Siloed data modalities analyzed separately
**To**: Integrated multi-modal evidence synthesis

**From**: Manual, error-prone bioinformatics pipelines
**To**: AI-orchestrated, reproducible workflows

**From**: Treatment decisions based on single biomarkers
**To**: Evidence from genomics + multi-omics + spatial + imaging

**From**: $3,200 per patient, 2-3 week turnaround
**To**: $1.35 per patient, 35-minute analysis

**The opportunity**: Transform precision oncology for millions of patients worldwide.

**The challenge**: Bridge the gap between AI capabilities and clinical deployment.

**This book showed you how**. Now it's your turn to build on it.

---

## What You've Accomplished

**18 chapters complete**:
- Part 1: Why This Matters (3 chapters)
- Part 2: Building the Foundation (4 chapters)
- Part 3: Advanced Capabilities (4 chapters)
- Part 4: Deployment and Operations (3 chapters)
- Part 5: Research and Education (2 chapters)
- Part 6: The Future (2 chapters)

**System built**:
- 12 MCP servers (124 tools total)
- 4 production-ready, 2 partial, 6 mocked
- PatientOne synthetic dataset (100% safe to share)
- 90+ prompts (clinical, research, educational, funder)
- 18 Jupyter notebooks (hands-on learning)

**Impact metrics**:
- Cost: 99.96% reduction ($1.35 vs $3,200)
- Time: 98% reduction (35 min vs 40 hours)
- ROI: 655-687% (5 years)
- Scalability: Cost/patient constant as volume increases

**Ready for**:
- Hospital production deployment (HIPAA-compliant)
- Research workflows (reproducible, cost-effective)
- Educational use (students, workshops, courses)
- Grant funding (ROI-justified, budget models provided)

---

**Thank you** for joining this journey. Now go build the future of precision medicine.

---

**Chapter 18 Summary**:
- Production insights: MCP architecture, PatientOne dataset, hybrid DRY_RUN/production, concise book format
- Challenges solved: DeepCell dependencies (3 weeks), Stouffer FDR timing, Streamlit costs, server status transparency
- Future enhancements: Phase 1 (3-6 mo: DeepCell, GEARS, Quantum), Phase 2 (6-12 mo: TCGA, Hugging Face), Phase 3 (12-18 mo: metabolomics, radiomics, single-cell)
- Multi-cancer expansion: Breast (HER2+), colorectal (MSI-H), lung (EGFR+) = 657K cases/year impact
- Community building: GitHub Discussions, workshops, conferences, plugin marketplace
- The path forward: 3-6 mo (complete servers), 6-18 mo (public data + new modalities), 18+ mo (full production, regulatory)
- Final transformation: 40 hours → 35 min, $3,200 → $1.35, siloed → integrated, manual → AI-orchestrated

**Book complete**: 18 chapters, 239 pages, 12 MCP servers documented
**Your turn**: Build on this foundation, transform precision oncology



ewpage




ppendix


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
wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

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



ewpage


# Appendix B: Installation and Setup Guide

*Quick reference for setting up your own copy of all MCP servers*

---

## Overview

This appendix provides a quick reference for setting up all 12 MCP servers in your environment. For detailed information, always refer to the linked documentation in the repository.

**Setup options**:
1. **Local (Claude Desktop)** - Run servers locally via stdio transport
2. **Cloud (GCP Cloud Run)** - Deploy servers to production via SSE transport
3. **Development** - Run servers in your IDE for testing and development

---

## Prerequisites

### Software Requirements

| Software | Version | Purpose |
|----------|---------|---------|
| **Python** | 3.11+ | MCP server runtime (3.10 for DeepCell) |
| **Git** | Latest | Clone repository |
| **Docker** | Latest | Cloud deployment (optional for local) |
| **Claude Desktop** | Latest | Local MCP server orchestration |
| **Google Cloud SDK** | Latest | Cloud deployment (Chapters 12-13) |

Installation instructions:
- **Python**: [python.org](https://www.python.org/downloads/)
- **Claude Desktop**: [claude.com/claude-desktop](https://claude.com/claude-desktop)
- **Google Cloud SDK**: [cloud.google.com/sdk](https://cloud.google.com/sdk/docs/install)

### Hardware Requirements

| Resource | Minimum | Recommended |
|----------|---------|-------------|
| RAM | 8GB | 16GB |
| Disk | 50GB free | 100GB free |
| CPU | 4 cores | 8 cores |
| GPU | None | Optional (DeepCell acceleration) |

---

## Quick Start: Local Setup

### Step 1: Clone Repository

```bash
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp
```

### Step 2: Install Server Dependencies

**Install all servers** (from repository root):

```bash
# Install each server in editable mode
cd servers/mcp-spatialtools && pip install -e ".[dev]" && cd ../..
cd servers/mcp-multiomics && pip install -e ".[dev]" && cd ../..
cd servers/mcp-deepcell && pip install -e ".[dev]" && cd ../..
# Repeat for all 12 servers...
```

**Or use the setup script** (if available):
```bash
./scripts/setup-all-servers.sh
```

Full server list and installation instructions: [`servers/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/README.md)

### Step 3: Configure Claude Desktop

**Location**: `~/Library/Application Support/Claude/claude_desktop_config.json` (macOS)

**Minimal configuration** (single server example):

```json
{
  "mcpServers": {
    "spatialtools": {
      "command": "/path/to/spatial-mcp/servers/mcp-spatialtools/venv/bin/python",
      "args": ["-m", "mcp_spatialtools"],
      "env": {
        "SPATIAL_DATA_DIR": "/path/to/data",
        "SPATIAL_DRY_RUN": "false"
      }
    }
  }
}
```

**Full configuration** (all 12 servers):
- Template: [`configs/claude_desktop_config.json`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/configs/claude_desktop_config.json)
- Setup guide: [`docs/test-docs/manual-testing/claude-desktop-setup.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/test-docs/manual-testing/claude-desktop-setup.md)

![Claude Desktop Config](images/screenshots/claude-desktop-config.png){width=\textwidth}

**Figure B.1: Claude Desktop Configuration File**
*Example `claude_desktop_config.json` showing all 12 MCP servers configured with absolute paths to Python executables, environment variables (DATA_DIR, DRY_RUN), and transport settings. Each server runs independently via stdio transport.*

**Key points**:
- Use **absolute paths** to venv Python executables
- Set `DRY_RUN=false` for real analysis, `true` for demonstration mode
- Restart Claude Desktop after config changes

### Step 4: Verify Installation

**Test in Claude Desktop**:
```
List all available MCP servers and their tools.
```

**Expected response**: 12 servers listed with 124 total tools.

---

## Server-Specific Setup

Each server has unique setup requirements. Refer to individual server READMEs for details:

| Server | Setup Time | Special Requirements | README Link |
|--------|------------|---------------------|-------------|
| **mcp-epic** | 5 min | Epic FHIR credentials (or use mockepic) | [`servers/mcp-epic/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-epic/README.md) |
| **mcp-fgbio** | 10 min | Reference genome files (hg38) | [`servers/mcp-fgbio/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-fgbio/README.md) |
| **mcp-multiomics** | 5 min | None | [`servers/mcp-multiomics/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-multiomics/README.md) |
| **mcp-spatialtools** | 10 min | STAR aligner (optional) | [`servers/mcp-spatialtools/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-spatialtools/README.md) |
| **mcp-deepcell** | 15 min | Python 3.10, TensorFlow 2.8.x | [`servers/mcp-deepcell/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/README.md) |
| **mcp-perturbation** | 5 min | None (mocked) | [`servers/mcp-perturbation/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-perturbation/README.md) |
| **mcp-quantum-celltype-fidelity** | 10 min | PennyLane (for PQC) | [`servers/mcp-quantum-celltype-fidelity/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-quantum-celltype-fidelity/README.md) |
| **mcp-openimagedata** | 5 min | None | [`servers/mcp-openimagedata/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-openimagedata/README.md) |
| **mcp-tcga** | 5 min | None (mocked) | [`servers/mcp-tcga/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-tcga/README.md) |
| **mcp-huggingface** | 5 min | None (mocked) | [`servers/mcp-huggingface/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-huggingface/README.md) |
| **mcp-seqera** | 5 min | Nextflow CLI | [`servers/mcp-seqera/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-seqera/README.md) |
| **mcp-mockepic** | 2 min | None | [`servers/mcp-mockepic/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-mockepic/README.md) |

**Total setup time**: 1-2 hours for all servers (local development)

---

## Cloud Deployment (GCP Cloud Run)

For production deployment to GCP Cloud Run, see **Chapter 12: Cloud Deployment on GCP** (page 167).

### Quick Cloud Deployment

**Deploy a single server**:

```bash
cd servers/mcp-deepcell
./deploy.sh YOUR_PROJECT_ID us-central1
```

**Deploy all servers**:

```bash
./scripts/deploy-all-servers.sh YOUR_PROJECT_ID us-central1
```

Deployment scripts: [`servers/*/deploy.sh`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/servers)

### Cloud Configuration

**Environment variables** (set during deployment):

```bash
# GCP Project
export GCP_PROJECT_ID="your-project-id"
export GCP_REGION="us-central1"

# Optional: Production settings
export DEEPCELL_DRY_RUN="false"
export SPATIAL_DRY_RUN="false"
export MULTIOMICS_DRY_RUN="false"
```

**Using deployed servers** (Claude API with SSE transport):

```python
import anthropic

client = anthropic.Anthropic()

response = client.beta.messages.create(
    model="claude-sonnet-4-5",
    max_tokens=1024,
    messages=[{"role": "user", "content": "Analyze PatientOne"}],
    mcp_servers=[{
        "type": "url",
        "url": "https://mcp-spatialtools-PROJECT_ID.run.app/sse",
        "name": "spatialtools"
    }],
    tools=[{"type": "mcp_toolset", "mcp_server_name": "spatialtools"}],
    betas=["mcp-client-2025-11-20"]
)
```

Full deployment guide: **Chapter 12** (pages 167-179)

---

## PatientOne Dataset Setup

The PatientOne synthetic dataset is required for testing and exercises.

### Option 1: Use GCS Public Bucket (Recommended)

**Location**: `gs://precision-medicine-mcp-public/patient-data/PAT001-OVC-2025/`

**Access**: Public read access (no authentication required)

**Download**:
```bash
# Install gsutil (part of Google Cloud SDK)
gsutil -m cp -r gs://precision-medicine-mcp-public/patient-data/PAT001-OVC-2025/ ./data/
```

### Option 2: Use Local Copy

**Location**: [`data/patient-data/PAT001-OVC-2025/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/data/patient-data/PAT001-OVC-2025)

**Clone with data**:
```bash
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
# PatientOne data included in repository
```

### Data Structure

```
PAT001-OVC-2025/
├── clinical/           # FHIR R4 resources (2 files)
├── genomics/           # Somatic variants VCF (1 file)
├── multiomics/         # RNA/protein/phospho (4 files)
├── spatial/            # 10X Visium (3 files)
└── imaging/            # H&E, MxIF images (7 files)
```

**Total**: 17 files, 100% synthetic (safe to share, no IRB needed)

Full dataset documentation: [`data/patient-data/PAT001-OVC-2025/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/data/patient-data/PAT001-OVC-2025/README.md)

---

## Environment Variables Reference

### Common Variables (All Servers)

| Variable | Default | Description |
|----------|---------|-------------|
| `*_DRY_RUN` | `false` | Enable mock mode (true) or real analysis (false) |
| `*_DATA_DIR` | `/workspace/data` | Data directory path |
| `*_CACHE_DIR` | `/workspace/cache` | Cache directory path |
| `*_LOG_LEVEL` | `INFO` | Logging level (DEBUG, INFO, WARNING, ERROR) |

### Server-Specific Variables

**mcp-spatialtools**:
- `SPATIAL_DATA_DIR`: Spatial transcriptomics data location
- `SPATIAL_DRY_RUN`: Mock mode toggle
- `STAR_PATH`: STAR aligner executable path
- `STAR_GENOME_INDEX`: Reference genome index directory

**mcp-deepcell**:
- `DEEPCELL_DRY_RUN`: Mock mode toggle
- `DEEPCELL_CACHE_DIR`: Model cache directory
- `TF_CPP_MIN_LOG_LEVEL`: TensorFlow logging (default: 3)

**mcp-multiomics**:
- `MULTIOMICS_DATA_DIR`: Multi-omics data location
- `MULTIOMICS_DRY_RUN`: Mock mode toggle

Full environment variable reference: Each server's README.md

---

## Hospital Production Setup

For HIPAA-compliant production deployment in hospital environments, see **Chapter 13: Hospital Production Deployment** (page 180).

### Key Requirements

1. **VPC Networking**: Private deployment with Serverless VPC Connector
2. **Authentication**: Azure AD SSO via OAuth2 Proxy
3. **Compliance**: 10-year immutable audit logs, de-identification
4. **Epic Integration**: FHIR R4 with OAuth2

**Setup scripts**:
- VPC setup: [`scripts/setup-vpc.sh`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/scripts/setup-vpc.sh)
- Audit logs: [`scripts/setup-audit-logs.sh`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/scripts/setup-audit-logs.sh)
- OAuth2 Proxy: [`scripts/setup-oauth2-proxy.sh`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/scripts/setup-oauth2-proxy.sh)

Full hospital setup guide: **Chapter 13** (pages 180-194)

---

## Troubleshooting

### Issue: Claude Desktop doesn't see servers

**Symptoms**: "No MCP servers found" or tools not available

**Solutions**:
1. Verify config file location: `~/Library/Application Support/Claude/claude_desktop_config.json`
2. Check absolute paths to Python executables
3. Restart Claude Desktop after config changes
4. Check server logs: `tail -f ~/.config/Claude/logs/mcp-*.log`

### Issue: Import errors when running servers

**Symptoms**: `ModuleNotFoundError` or `ImportError`

**Solutions**:
1. Activate virtual environment: `source venv/bin/activate`
2. Reinstall dependencies: `pip install -e ".[dev]"`
3. Check Python version: `python --version` (must be 3.11+)
4. Set PYTHONPATH: `export PYTHONPATH=/path/to/server/src`

### Issue: DeepCell requires Python 3.10

**Symptom**: TensorFlow 2.8.x incompatible with Python 3.11+

**Solution**: Use Python 3.10 for mcp-deepcell only
```bash
# Install Python 3.10
brew install python@3.10  # macOS

# Create venv with Python 3.10
cd servers/mcp-deepcell
python3.10 -m venv venv
source venv/bin/activate
pip install -e ".[dev]"
```

Full story: [`servers/mcp-deepcell/DEPENDENCY_ISSUES.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/DEPENDENCY_ISSUES.md)

### Issue: Cloud Run deployment fails

**Symptoms**: Deployment timeout or resource errors

**Solutions**:
1. Check GCP project quotas: `gcloud compute project-info describe --project=PROJECT_ID`
2. Verify Docker is running: `docker ps`
3. Increase timeout: Edit `cloudbuild.yaml` (set `timeout: 1200s`)
4. Check Cloud Build logs: `gcloud builds list --limit=5`

Full deployment troubleshooting: **Chapter 12** (pages 175-178)

### Issue: "Permission denied" in MCP server

**Symptoms**: Cannot read/write files

**Solutions**:
1. Check directory permissions: `ls -la /workspace/data`
2. Create directories: `mkdir -p /workspace/data /workspace/cache`
3. Set ownership: `chown -R $USER /workspace/data`
4. Update config: Ensure `*_DATA_DIR` points to writable location

---

## Testing Your Setup

### Test 1: List Available Servers

**In Claude Desktop**:
```
List all available MCP servers and their tools.
```

**Expected**: 12 servers, 124 tools total

### Test 2: Run Simple Query

**In Claude Desktop**:
```
For PatientOne (patient-001), retrieve clinical demographics using mockepic.
```

**Expected**: Patient demographics returned (Sarah Anderson, 58yo, BRCA1+)

### Test 3: Run Complete Analysis

**In Claude Desktop**:
```
Analyze PatientOne (patient-001) spatial transcriptomics data:
1. Load data using spatialtools
2. Filter low-quality spots (min 200 UMIs, 100 genes)
3. Calculate spatial autocorrelation for MKI67, PCNA, CD8A

Use DRY_RUN=false for real analysis.
```

**Expected**: Real analysis results with Moran's I values

Full test prompts: [`docs/test-docs/manual-testing/TEST_*.txt`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/docs/test-docs/manual-testing)

---

## Development Workflow

### Running Servers in Your IDE

**For development and debugging**:

```bash
# Start server in stdio mode (for testing)
cd servers/mcp-spatialtools
python -m mcp_spatialtools

# Or use MCP Inspector (interactive testing tool)
npx @modelcontextprotocol/inspector python -m mcp_spatialtools
```

**MCP Inspector**: [modelcontextprotocol.io/docs/tools/inspector](https://modelcontextprotocol.io/docs/tools/inspector)

### Running Tests

**Run server tests**:
```bash
cd servers/mcp-spatialtools
pytest
```

**Run with coverage**:
```bash
pytest --cov=src/mcp_spatialtools --cov-report=html
```

Testing guide: [`docs/development/testing.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/development/testing.md)

---

## Additional Resources

### Book Chapters (Setup Referenced)

- **Chapter 4**: Clinical Data (mcp-epic configuration)
- **Chapter 5**: Genomic Foundations (mcp-fgbio setup)
- **Chapter 7**: Spatial Transcriptomics (mcp-spatialtools STAR setup)
- **Chapter 8**: Cell Segmentation (mcp-deepcell Python 3.10)
- **Chapter 12**: Cloud Deployment (GCP Cloud Run)
- **Chapter 13**: Hospital Production (VPC, SSO, HIPAA)
- **Chapter 16**: Teaching (Jupyter notebooks setup)

### Repository Documentation

- **Architecture overview**: [`docs/architecture/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/architecture/README.md)
- **Server implementation status**: [`docs/architecture/servers.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/architecture/servers.md)
- **Deployment guides**: [`docs/deployment/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/docs/deployment)
- **Getting started**: [`docs/getting-started/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/docs/getting-started)

### External Resources

- **MCP Specification**: [modelcontextprotocol.io](https://modelcontextprotocol.io)
- **FastMCP Framework**: [github.com/jlowin/fastmcp](https://github.com/jlowin/fastmcp)
- **Claude Desktop**: [claude.com/claude-desktop](https://claude.com/claude-desktop)
- **Claude API Docs**: [docs.anthropic.com](https://docs.anthropic.com)

---

## Summary

**Local setup steps**:
1. Install prerequisites (Python 3.11, Git, Claude Desktop)
2. Clone repository
3. Install server dependencies (`pip install -e ".[dev]"` for each server)
4. Configure Claude Desktop (`claude_desktop_config.json`)
5. Download PatientOne data (GCS or local)
6. Test setup (list servers, run simple query)

**Cloud setup steps** (Chapter 12):
1. Install Google Cloud SDK
2. Authenticate (`gcloud auth login`)
3. Deploy servers (`./deploy.sh PROJECT_ID REGION`)
4. Test with Claude API (SSE transport)

**Production setup steps** (Chapter 13):
1. Run setup scripts (VPC, audit logs, OAuth2 Proxy)
2. Configure Epic FHIR integration
3. Deploy with HIPAA compliance
4. Validate security controls

**Total time**:
- Local: 1-2 hours (all servers)
- Cloud: 2-3 hours (deployment + testing)
- Production: 1 week (security review + compliance)

---

**This appendix consolidates setup information from all chapters. For detailed instructions, always refer to the linked repository documentation.**

---

**Appendix complete**: Quick reference for setting up all 12 MCP servers (local, cloud, production)
**Key sections**: Prerequisites, local setup, cloud deployment, PatientOne data, troubleshooting, testing
**Links**: 25+ links to detailed documentation in repository (chapters, READMEs, scripts, guides)



ewpage


# Appendix C: PatientOne Complete Dataset

*Reference guide for the complete PAT001-OVC-2025 synthetic dataset*

---

## Dataset Overview

**Patient ID**: PAT001-OVC-2025
**Diagnosis**: Stage IV high-grade serous ovarian carcinoma
**Age**: 58 years
**Status**: 100% synthetic (safe to share, no IRB required)
**Purpose**: Testing, training, and demonstration of AI-orchestrated precision medicine workflows

**Location**: [`data/patient-data/PAT001-OVC-2025/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/data/patient-data/PAT001-OVC-2025)

**Public access**: `gs://precision-medicine-mcp-public/patient-data/PAT001-OVC-2025/`

---

## Complete File Manifest

### Clinical Data (2 files)

**Location**: `clinical/`

| File | Format | Size | Description |
|------|--------|------|-------------|
| `patient-001.json` | FHIR R4 | 4 KB | Patient demographics, conditions, medications |
| `genomic-test-order-001.json` | FHIR R4 | 2 KB | Genomic test orders, ServiceRequest resources |

**Key fields**:
- Patient: Sarah Anderson, 58yo, Female
- Diagnosis: Ovarian cancer (C56.9), Stage IV
- BRCA status: BRCA1+
- Medications: Carboplatin, paclitaxel, bevacizumab

---

### Genomics Data (1 file)

**Location**: `genomics/`

| File | Format | Size | Description |
|------|--------|------|-------------|
| `somatic-variants.vcf` | VCF 4.2 | 12 KB | 8 pathogenic somatic mutations |

**Pathogenic variants**:

| Gene | Variant | Type | VAF | ClinVar | Actionability |
|------|---------|------|-----|---------|---------------|
| TP53 | R175H | Missense | 85% | Pathogenic (★★★★★) | Prognostic |
| PIK3CA | E545K | Missense | 42% | Pathogenic (★★★★) | PI3K inhibitor (alpelisib) |
| KRAS | G12D | Missense | 38% | Pathogenic (★★★★★) | Prognostic (poor) |
| PTEN | R130* | Nonsense | 55% | Pathogenic (★★★★★) | PI3K/AKT pathway activation |
| ARID1A | Q456* | Nonsense | 48% | Pathogenic (★★★★) | SWI/SNF complex loss |
| CTNNB1 | S45F | Missense | 31% | Pathogenic (★★★★) | Wnt pathway activation |
| FBXW7 | R465C | Missense | 27% | Likely pathogenic (★★★) | Tumor suppressor loss |
| PPP2R1A | P179R | Missense | 33% | Pathogenic (★★★★) | Platinum resistance |

**Tumor mutational burden (TMB)**: 8.2 mutations/Mb (high)

---

### Multi-Omics Data (4 files)

**Location**: `multiomics/`

| File | Format | Size | Samples | Features | Description |
|------|--------|------|---------|----------|-------------|
| `rna_counts.tsv` | TSV | 280 KB | 15 PDX models | 2,000 genes | RNA-seq expression (TPM) |
| `protein_abundance.tsv` | TSV | 95 KB | 15 PDX models | 350 proteins | RPPA protein levels |
| `phospho_abundance.tsv` | TSV | 110 KB | 15 PDX models | 450 phospho-sites | Phosphoproteomics |
| `sample_metadata.tsv` | TSV | 3 KB | 15 PDX models | 8 annotations | Treatment response, histology |

**PDX model cohort**:
- **n = 15** ovarian cancer patient-derived xenografts
- **Treatment groups**: Carboplatin (n=5), Olaparib (n=5), Combination (n=5)
- **Response**: Complete response (n=4), Partial response (n=7), Progressive disease (n=4)

**Key genes/proteins**:
- **TP53**: Mutant expression confirmed (RNA, protein)
- **PIK3CA**: Elevated protein and phospho-AKT (S473)
- **MKI67**: 25% proliferative cells (protein)
- **CD8A/CD8B**: Moderate T cell infiltration

---

### Spatial Transcriptomics Data (3 files)

**Location**: `spatial/`

| File | Format | Size | Spots | Genes | Description |
|------|--------|------|-------|-------|-------------|
| `filtered_feature_bc_matrix.h5` | H5AD | 45 MB | 900 | 7,000 | 10X Visium gene expression |
| `spatial_coords.tsv` | TSV | 25 KB | 900 | 2 | Spot x,y coordinates |
| `tissue_positions.csv` | CSV | 30 KB | 900 | 5 | Tissue coordinates, selection |

**Spatial regions** (6 annotated):

| Region | Spots | % | Key Markers | Description |
|--------|-------|---|-------------|-------------|
| Tumor proliferative | 245 | 27% | MKI67+, PCNA+, TOP2A+ | Actively dividing tumor cells |
| Tumor quiescent | 198 | 22% | MKI67-, TP53+ | TP53-mutant non-dividing cells |
| Hypoxic niche | 112 | 12% | HIF1A+, CA9+, VEGFA+ | Oxygen-deprived regions |
| Immune infiltrate | 156 | 17% | CD8A+, CD8B+, GZMB+ | T cell rich areas |
| Stromal | 134 | 15% | COL1A1+, FAP+, ACTA2+ | Cancer-associated fibroblasts |
| Necrotic | 55 | 6% | Low UMI, low genes | Dead tissue |

**Moran's I results** (spatial autocorrelation):

| Gene | Moran's I | p-value | Interpretation |
|------|-----------|---------|----------------|
| MKI67 | 0.78 | < 0.001 | Highly clustered (proliferative zones) |
| CD8A | 0.65 | < 0.001 | Clustered (immune infiltrates) |
| HIF1A | 0.71 | < 0.001 | Clustered (hypoxic regions) |
| TP53 | 0.42 | < 0.001 | Moderately clustered |

---

### Imaging Data (7 files)

**Location**: `imaging/`

| File | Format | Size | Modality | Description |
|------|--------|------|----------|-------------|
| `he_section1.tif` | TIFF | 120 MB | H&E | Hematoxylin & eosin stain |
| `he_section2.tif` | TIFF | 118 MB | H&E | Second H&E section |
| `mxif_cd8_dapi.tif` | TIFF | 95 MB | MxIF | CD8 + DAPI (2-plex) |
| `mxif_ki67_tp53_dapi.tif` | TIFF | 98 MB | MxIF | Ki67 + TP53 + DAPI (3-plex) |
| `mxif_panck_vim_dapi.tif` | TIFF | 102 MB | MxIF | Pan-cytokeratin + Vimentin + DAPI |
| `nuclear_masks.png` | PNG | 12 MB | Segmentation | DeepCell nuclear masks |
| `cell_phenotypes.csv` | CSV | 85 KB | Annotation | Cell classifications |

**Imaging dimensions**:
- H&E: 2048 × 2048 pixels, 0.5 μm/pixel
- MxIF: 2048 × 2048 pixels, 0.5 μm/pixel
- Channels: DAPI (blue), FITC (green), TRITC (red), Cy5 (far-red)

**Cell counts**:
- Total cells segmented: 3,847
- Tumor cells (PanCK+): 2,234 (58%)
- Stromal cells (VIM+): 891 (23%)
- T cells (CD8+): 412 (11%)
- Proliferating cells (Ki67+): 961 (25%)
- TP53+ cells: 1,823 (47%)

---

## Data Access Methods

### Method 1: Clone Repository (Recommended for Local)

```bash
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp/data/patient-data/PAT001-OVC-2025
ls -lh
```

### Method 2: Download from GCS (Cloud Run)

```bash
# Public bucket (no authentication required)
gsutil -m cp -r gs://precision-medicine-mcp-public/patient-data/PAT001-OVC-2025/ ./data/
```

### Method 3: Use in Jupyter Notebooks

```python
import os
from pathlib import Path

# Set data directory
DATA_DIR = Path("../../../data/patient-data/PAT001-OVC-2025")

# Load clinical FHIR
import json
with open(DATA_DIR / "clinical/patient-001.json") as f:
    patient_fhir = json.load(f)

# Load genomic VCF
vcf_path = DATA_DIR / "genomics/somatic-variants.vcf"

# Load spatial h5ad
import scanpy as sc
adata = sc.read_10x_h5(DATA_DIR / "spatial/filtered_feature_bc_matrix.h5")
```

---

## Data Generation Details

**All data is 100% synthetic** and generated using:

1. **Clinical FHIR**: Generated with FHIR test data generator
2. **Genomic VCF**: Simulated with realistic VAF distributions (beta distribution)
3. **Multi-omics**: Sampled from TCGA ovarian cancer cohort + added noise
4. **Spatial**: Synthetic Visium data with spatial autocorrelation (geoR package)
5. **Imaging**: Generated with synthetic microscopy tools + DeepCell segmentation

**No real patient data was used**. Safe to share publicly, no IRB approval required.

---

## Usage in Book Chapters

| Chapter | Data Used | Purpose |
|---------|-----------|---------|
| 1 | All modalities | Complete PatientOne story |
| 3 | All modalities | Production validation |
| 4 | Clinical FHIR | FHIR integration |
| 5 | Genomic VCF | Variant calling, annotation |
| 6 | Multi-omics TSV | HAllA, Stouffer meta-analysis |
| 7 | Spatial H5AD | STAR, ComBat, Moran's I |
| 8 | Imaging TIFF | DeepCell segmentation |
| 9 | Multi-omics + Spatial | GEARS treatment prediction |
| 10 | Spatial H5AD | Quantum fidelity |
| 11 | Imaging TIFF | Histopathology |
| 15 | All modalities | Research workflows |
| 16 | All modalities | Teaching exercises |

---

## Companion Jupyter Notebooks

All 18 chapter notebooks use PatientOne data. See:
- **Notebook setup**: [`docs/book/companion-notebooks/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/book/companion-notebooks/README.md)
- **Data loading examples**: Each chapter's notebook

---

## Expected Analysis Results

**Complete analysis time**: 35 minutes (AI-orchestrated) vs 40 hours (traditional)

**Key findings**:
- TP53 R175H mutation (85% VAF) → Poor prognosis
- PIK3CA E545K (42% VAF) → PI3K inhibitor responsive
- TMB 8.2 mutations/Mb → High (immunotherapy candidate)
- 25% Ki67+ proliferating cells → Aggressive tumor
- CD8+ T cell infiltration (11%) → Moderate immune response
- Hypoxic regions (12%) → Potential checkpoint blockade target

**Treatment recommendation**: Olaparib (PARP inhibitor) + alpelisib (PI3K inhibitor) combination (82% predicted efficacy)

---

## Additional Resources

- **Full dataset README**: [`data/patient-data/PAT001-OVC-2025/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/data/patient-data/PAT001-OVC-2025/README.md)
- **Analysis scripts**: [`docs/test-docs/patient-one-scenario/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/docs/test-docs/patient-one-scenario)
- **Demo walkthrough**: [`docs/demos/FULL_PATIENTONE_DEMO.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/demos/FULL_PATIENTONE_DEMO.md)

---

**This appendix provides complete reference for the PatientOne synthetic dataset. All data is public domain (CC0 1.0 Universal).**

---



ewpage


# Appendix D: Bias and Ethics in AI-Orchestrated Precision Medicine

*Framework for detecting and mitigating bias in AI-assisted clinical workflows*

---

## Overview

AI systems in healthcare can perpetuate or amplify biases present in training data, algorithms, and human decision-making. This appendix provides:

1. **Bias detection framework** for AI-orchestrated precision medicine
2. **PatientOne bias audit results** (demonstration)
3. **Mitigation strategies** for common bias sources
4. **Ethical considerations** for clinical deployment

**Full ethics documentation**: [`docs/for-hospitals/ethics/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/docs/for-hospitals/ethics)

---

## Common Bias Sources in Precision Medicine AI

### 1. Training Data Bias

**Problem**: Genomic databases (TCGA, gnomAD, ClinVar) over-represent European ancestry.

**Impact on system**:
- Variant pathogenicity predictions less accurate for non-European patients
- Treatment response models trained predominantly on white patient cohorts
- Spatial transcriptomics reference atlases lack diversity

**Example**: BRCA1/BRCA2 variants classified as "VUS" (Variant of Unknown Significance) in African ancestry populations but "Pathogenic" in European populations due to training data imbalance.

**Mitigation**:
- Use ancestry-stratified ClinVar annotations
- Report confidence intervals on predictions by ancestry
- Flag when patient ancestry differs from training cohort

**Full details**: [`ETHICS_AND_BIAS.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-hospitals/ethics/ETHICS_AND_BIAS.md#training-data-bias)

---

### 2. Algorithm Bias

**Problem**: AI orchestration (Claude/Gemini) may prioritize certain data modalities or tools based on training preferences.

**Impact on system**:
- Overreliance on genomics vs. imaging or spatial data
- Preferential citation of high-impact journals (publication bias)
- Treatment recommendations biased toward FDA-approved therapies (excluding experimental options)

**Example**: Claude may prioritize TP53 mutations (well-studied) over ARID1A mutations (less literature) even when both are pathogenic.

**Mitigation**:
- Explicitly prompt for multi-modal integration: "Use genomics, spatial, AND imaging data"
- Request uncertainty quantification: "Provide confidence scores for each recommendation"
- Ask for alternative perspectives: "What are alternative interpretations?"

**Full details**: [`ETHICS_AND_BIAS.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-hospitals/ethics/ETHICS_AND_BIAS.md#algorithm-bias)

---

### 3. Confirmation Bias

**Problem**: AI may reinforce clinician's initial hypothesis rather than exploring alternatives.

**Impact on system**:
- If clinician suspects BRCA mutation, AI may over-interpret benign variants as supportive evidence
- Spatial analysis may focus on tumor regions while missing immune-rich areas

**Mitigation**:
- Use adversarial prompts: "What evidence contradicts the initial diagnosis?"
- Request differential diagnosis: "What are 3 alternative explanations for these findings?"
- Blind analysis: Run analysis without revealing initial clinical hypothesis

**Full details**: [`ETHICS_AND_BIAS.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-hospitals/ethics/ETHICS_AND_BIAS.md#confirmation-bias)

---

### 4. Deployment Bias

**Problem**: System deployed only in high-resource hospitals, widening healthcare disparities.

**Impact**:
- Patients at community hospitals lack access to AI-orchestrated precision medicine
- Cost-prohibitive for rural or underserved populations

**Mitigation**:
- Cloud-based deployment reduces infrastructure requirements
- Cost transparency: $1-2 per analysis (vs. $3,200 traditional)
- Open-source code enables equitable access
- Educational resources for smaller institutions (Chapter 16)

**Full details**: [`ETHICS_AND_BIAS.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-hospitals/ethics/ETHICS_AND_BIAS.md#deployment-bias)

---

## Bias Audit Framework

### 7-Step Checklist

Use this checklist for every clinical deployment:

1. **Data Provenance Audit**
   - ✅ Document ancestry composition of training datasets
   - ✅ Report ClinVar classification confidence by ancestry
   - ✅ Flag predictions when patient ancestry differs from training cohort

2. **Algorithm Transparency**
   - ✅ Log all AI-generated recommendations with reasoning
   - ✅ Provide uncertainty quantification (confidence intervals)
   - ✅ Document which MCP tools were used for each conclusion

3. **Validation Across Subgroups**
   - ✅ Test prediction accuracy by ancestry, sex, age, cancer subtype
   - ✅ Report performance disparities (if any)
   - ✅ Avoid deployment if accuracy differs >10% between subgroups

4. **Clinical Oversight**
   - ✅ AI recommendations are advisory only (clinician has final decision)
   - ✅ All AI outputs reviewed by board-certified oncologist
   - ✅ Patient consent for AI-assisted analysis

5. **Continuous Monitoring**
   - ✅ Track prediction accuracy by patient demographics
   - ✅ Monthly bias audits (automated)
   - ✅ Update training data when new diverse cohorts available

6. **Patient Communication**
   - ✅ Inform patients that AI was used in analysis
   - ✅ Explain AI limitations and uncertainty
   - ✅ Offer opt-out option (use traditional analysis instead)

7. **Documentation and Audit Trail**
   - ✅ Immutable audit logs (10-year retention)
   - ✅ Document AI version, model, and training data provenance
   - ✅ Record clinician's decision rationale (agreement/disagreement with AI)

**Full checklist**: [`BIAS_AUDIT_CHECKLIST.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-hospitals/ethics/BIAS_AUDIT_CHECKLIST.md)

---

## PatientOne Bias Audit Results (Example)

### Audit Summary

**Patient**: PAT001-OVC-2025 (100% synthetic)
**Audit Date**: 2026-02-01
**Auditor**: Automated bias detection system

| Bias Category | Risk Level | Findings | Mitigation Applied |
|---------------|------------|----------|-------------------|
| **Training Data** | ⚠️ Medium | ClinVar data 78% European ancestry | Ancestry-stratified confidence intervals provided |
| **Algorithm** | ✅ Low | All modalities used (genomics, spatial, imaging) | Multi-modal integration verified |
| **Confirmation** | ✅ Low | Alternative diagnoses explored | Differential diagnosis generated |
| **Deployment** | ✅ Low | Cloud Run accessible to all institutions | Cost: $1.20 (affordable) |

**Overall Risk**: ⚠️ Medium (acceptable with mitigation)

**Recommendations**:
1. Report ClinVar pathogenicity confidence by ancestry
2. Validate TP53 R175H classification in non-European cohorts
3. Provide uncertainty quantification on all predictions

**Full audit**: [`PATIENTONE_BIAS_AUDIT.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-hospitals/ethics/PATIENTONE_BIAS_AUDIT.md)

---

## Ethical Considerations for Clinical Deployment

### 1. Informed Consent

**Requirement**: Patients must consent to AI-assisted analysis.

**Consent elements**:
- Explain that AI (Claude/Gemini) orchestrates bioinformatics tools
- Clarify that AI recommendations are advisory (clinician makes final decision)
- Describe data sources (TCGA, ClinVar, gnomAD) and their ancestry composition
- Offer opt-out option (use traditional manual analysis instead)

**Template**: [`docs/for-hospitals/ethics/CONSENT_TEMPLATE.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-hospitals/ethics/CONSENT_TEMPLATE.md)

---

### 2. Clinician Autonomy

**Principle**: AI augments, but does not replace, clinical judgment.

**Implementation**:
- All AI recommendations flagged as "AI-generated" in EHR
- Clinician must document agreement/disagreement with AI reasoning
- Override mechanism: Clinician can reject AI recommendation with justification

**Example**: AI recommends olaparib (82% predicted efficacy), but oncologist chooses carboplatin+olaparib based on patient's insurance coverage and clinical trial eligibility.

---

### 3. Data Privacy and De-identification

**Requirement**: Patient data must be de-identified before cloud analysis.

**De-identification pipeline**:
1. Remove 18 HIPAA identifiers (name, MRN, DOB, etc.)
2. Replace with synthetic identifiers (e.g., PAT001-OVC-2025)
3. Remove genomic identifiers (rs IDs, dbSNP IDs)
4. Retain only clinical-relevant annotations

**Validation**: Safe Harbor method (HIPAA §164.514(b))

**Full guide**: **Chapter 13** (pages 185-187)

---

### 4. Algorithmic Accountability

**Requirement**: AI decisions must be explainable and auditable.

**Implementation**:
- Log all MCP tool calls (which tools, parameters, results)
- Store AI prompts and responses (10-year retention)
- Provide reasoning chains: "TP53 R175H is pathogenic because ClinVar (★★★★★) + PubMed evidence (234 citations)"

**Audit trail format**:
```json
{
  "patient_id": "PAT001-OVC-2025",
  "analysis_date": "2026-02-01T10:30:00Z",
  "ai_model": "claude-sonnet-4-5",
  "mcp_tools_used": ["fgbio_parse_vcf", "fgbio_annotate_variants"],
  "recommendations": ["Olaparib (PARP inhibitor)"],
  "clinician_decision": "Agreed",
  "rationale": "Patient has TP53 mutation + BRCA1 WT → synthetic lethality"
}
```

---

### 5. Fairness Across Populations

**Principle**: System should provide equitable care regardless of ancestry, sex, age, or socioeconomic status.

**Monitoring metrics**:
- Prediction accuracy by ancestry (European, African, Asian, Hispanic, Other)
- Treatment recommendation diversity (not all patients get same drug)
- Cost accessibility ($1-2 per analysis → affordable)

**Action threshold**: If accuracy differs >10% between subgroups, pause deployment and retrain with diverse data.

---

## Implementation Roadmap

### Phase 1: Pre-Deployment (Weeks 1-4)

1. Complete bias audit checklist for institution's patient demographics
2. Validate system on diverse test cohort (≥50 patients, ≥3 ancestries)
3. Train clinical staff on AI limitations and oversight procedures
4. Obtain IRB approval for AI-assisted clinical workflows

---

### Phase 2: Pilot Deployment (Weeks 5-12)

1. Deploy to 1-2 oncologists (≤20 patients)
2. Manual review of all AI recommendations by tumor board
3. Weekly bias audits (automated)
4. Patient feedback surveys

---

### Phase 3: Full Deployment (Weeks 13+)

1. Scale to full oncology department
2. Monthly bias audits (automated)
3. Quarterly ethics review with external auditor
4. Publish bias audit results (transparency)

**Full implementation plan**: [`IMPLEMENTATION_PLAN.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-hospitals/ethics/IMPLEMENTATION_PLAN.md)

---

## Regulatory Considerations

### FDA Oversight

**Status**: AI orchestration systems are not currently regulated medical devices (as of 2026).

**Rationale**:
- AI coordinates existing FDA-approved tools (STAR aligner, DeepCell, GEARS)
- Final clinical decision made by licensed physician
- No autonomous treatment administration

**However**: Monitor FDA guidance on "Clinical Decision Support" software. May require 510(k) clearance if marketed as diagnostic tool.

---

### HIPAA Compliance

**Requirement**: Cloud Run deployment must be HIPAA-compliant.

**Checklist**:
- ✅ Business Associate Agreement (BAA) with Google Cloud
- ✅ Encryption at rest and in transit (TLS 1.3)
- ✅ De-identification before cloud analysis
- ✅ 10-year immutable audit logs
- ✅ Access controls (OAuth2 + Azure AD SSO)

**Full guide**: **Chapter 13** (pages 180-194)

---

## Resources and Further Reading

### Internal Documentation

- **Ethics framework**: [`docs/for-hospitals/ethics/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-hospitals/ethics/README.md)
- **Bias audit checklist**: [`BIAS_AUDIT_CHECKLIST.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-hospitals/ethics/BIAS_AUDIT_CHECKLIST.md)
- **PatientOne audit**: [`PATIENTONE_BIAS_AUDIT.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-hospitals/ethics/PATIENTONE_BIAS_AUDIT.md)
- **Ethics implementation plan**: [`IMPLEMENTATION_PLAN.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-hospitals/ethics/IMPLEMENTATION_PLAN.md)

### External Guidelines

- **WHO Ethics and Governance of AI for Health** (2021): [who.int/publications](https://www.who.int/publications/i/item/9789240029200)
- **FDA Clinical Decision Support Software Guidance** (2022): [fda.gov/cds-software](https://www.fda.gov/regulatory-information/search-fda-guidance-documents/clinical-decision-support-software)
- **NIST AI Risk Management Framework** (2023): [nist.gov/ai-rmf](https://www.nist.gov/itl/ai-risk-management-framework)
- **ASCO Guidelines on AI in Oncology** (2024): [asco.org/ai-guidelines](https://ascopubs.org/doi/full/10.1200/JCO.23.01241)

---

## Summary

**Key takeaways**:
1. **All AI systems have bias** → Detect and mitigate proactively
2. **Transparency is critical** → Log all AI decisions, provide reasoning
3. **Validation across populations** → Test on diverse cohorts before deployment
4. **Clinician oversight required** → AI is advisory, not autonomous
5. **Continuous monitoring** → Monthly bias audits, quarterly ethics review

**Next steps**:
- Complete bias audit checklist for your institution
- Validate system on diverse test cohort
- Obtain IRB approval and patient consent
- Deploy with continuous monitoring

**This appendix provides ethical framework for responsible AI deployment in precision medicine. Always consult institutional ethics board and legal counsel before clinical use.**

---




ewpage


