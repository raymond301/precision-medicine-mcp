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

![System Overview](images/screenshots/patient-one-holistic.png){width=100%}

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
