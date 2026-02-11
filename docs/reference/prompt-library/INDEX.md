# Prompt Library Index

**Complete index of all prompts in the Precision Medicine MCP library**

---

## Quick Navigation by Audience

| Audience | File | Prompts | Time Range |
|----------|------|---------|------------|
| **Funders & Decision-Makers** | [funder-demo-prompts.md](funder-demo-prompts.md) | 7 prompts | 5-10 min each |
| **Clinicians & Researchers** | [clinical-workflow-prompts.md](clinical-workflow-prompts.md) | 15 prompts | 5-45 min each |
| **Developers & QA** | [developer-test-prompts.md](developer-test-prompts.md) | 11 prompts | 1-15 min each |
| **Students & Educators** | [educational-prompts.md](educational-prompts.md) | 10 prompts | 5-20 min each |

**Total Prompts:** 43 curated prompts + 50+ discovered prompts = 90+ prompts available

---

## Funder Demo Prompts (7 Prompts)

**Audience:** Grant reviewers, hospital CFOs, funding agencies
**Focus:** ROI, cost savings, clinical impact, scalability

### High-Impact Demonstrations

1. **Executive Summary Demo** (5-10 min)
   - Complete PatientOne workflow overview
   - **Output:** Clinical impact, cost savings ($3,137/patient), time savings (98% reduction)
   - **Metrics:** $313K-1.6M annual savings, 2-3 patient payback period

2. **ROI Analysis** (5-10 min)
   - Financial justification for institutional adoption
   - **Output:** Cost comparison tables, break-even analysis, 5-year projections
   - **Metrics:** 96% cost reduction, 2-3 patient payback, 3,942 hours saved annually

3. **Clinical Impact Demo** (5-10 min)
   - Multi-modal evidence improves treatment decisions
   - **Output:** Confidence improvement (65% → 92%), response rate increase (30% → 50%)
   - **Metrics:** +4 months PFS, 20-30% fewer treatment failures

4. **Scalability Demo** (5-10 min)
   - Pilot to institutional deployment roadmap
   - **Output:** Scaling plan (5 → 100 users), cost per patient remains constant
   - **Metrics:** Sub-linear staffing, auto-scaling infrastructure

5. **HIPAA Compliance Demo** (5-10 min)
   - Security architecture and regulatory compliance
   - **Output:** 18 identifiers removed, 10-year audit logs, encryption specs
   - **Metrics:** 100% de-identification accuracy, 15/15 HIPAA controls

6. **Multi-Cancer Extensibility** (5-10 min)
   - Platform extension to breast, colorectal, other cancers
   - **Output:** 2-4 week deployment per cancer type, $7K-15K one-time cost
   - **Metrics:** $313K annual ROI per cancer type

7. **Real-World Evidence Generation** (5-10 min)
   - Research value for publications and clinical trials
   - **Output:** 73% faster time to publication, $2M-4M trial savings
   - **Metrics:** 3× publication rate, $120K-160K savings per study

---

## Clinical Workflow Prompts (15 Prompts)

**Audience:** Oncologists, clinical researchers, bioinformaticians
**Focus:** Patient analysis, treatment recommendations, scientific rigor

### Complete PatientOne Workflow (6 Prompts)

> **Sample Data:** PatientOne data (PAT001-OVC-2025) is in GCS at `gs://sample-inputs-patientone/patient-data/PAT001-OVC-2025/` with sub-directories for `multiomics/`, `spatial/`, and `imaging/`.

1. **Clinical Data & Genomic Profiling** (5-10 min)
   - Retrieve demographics, CA-125 trends, somatic mutations
   - **Output:** Patient summary, resistance evidence, TCGA subtype

2. **Multi-Omics Resistance Analysis** (15-25 min)
   - Preprocessing, integration, meta-analysis, upstream regulators
   - **Output:** PI3K/AKT activation, drug targets (alpelisib, capivasertib)

3. **Spatial Tumor Microenvironment** (10-15 min)
   - 900 spots across 6 regions, spatial patterns, visualizations
   - **Output:** Immune exclusion confirmed, heterogeneous resistance markers

4. **Histopathology & Immunofluorescence** (20-40 min)
   - H&E morphology, CD8 infiltration, Ki67 index, multiplex IF
   - **Output:** Cold tumor phenotype, high proliferation, TP53+/Ki67+ cells

5. **Integrated Multi-Modal Analysis** (5-10 min)
   - Synthesize findings, rank resistance mechanisms, treatment recommendations
   - **Output:** 3 treatment recommendations, monitoring strategy, integrated figure

6. **Clinician-in-the-Loop Validation** (30-45 min)
   - Draft report, quality gates, clinician review, final approval
   - **Output:** APPROVED clinical report with digital signature, audit trail

### Additional Clinical Use Cases (9 Prompts)

7. **BRCA1/2 Testing & PARP Inhibitor Eligibility** (5 min)
8. **Bevacizumab Continuation Assessment** (5 min)
9. **Immunotherapy Eligibility** (10 min) - PD-L1, TMB, MSI assessment
10. **Platinum Re-Challenge vs Alternative Chemotherapy** (5-10 min)
11. **Clinical Trial Matching** (10 min)
12. **Circulating Tumor DNA Monitoring** (5 min)
13. **Comprehensive Molecular Tumor Board Presentation** (15 min)
14. **Adverse Event Management** (3 min) - Hypersensitivity reaction
15. **End-of-Treatment Summary & Survivorship Plan** (10 min)

---

## Developer Test Prompts (11 Prompts)

**Audience:** Developers, DevOps engineers, QA testers
**Focus:** Server validation, integration testing, error handling

### Server Validation (9 Prompts)

1. **List All Available Tools** (<1 min)
   - Quick smoke test for any server
   - **Output:** Tool inventory, descriptions, parameter lists

2. **Validate mcp-multiomics Configuration** (2 min)
   - Server health check, tool categories, dry-run test
   - **Output:** 9 tools verified, documentation complete

3. **Validate mcp-spatialtools with PatientOne** (5 min)
   - Data loading, processing, output validation, error handling
   - **Output:** 900 spots loaded, Moran's I calculated, runtime <10s

4. **Validate mcp-fgbio Reference Genome Tool** (3 min)
   - Resource queries, gene annotations, FASTQ validation
   - **Output:** hg38 metadata, BRCA1 coordinates, tool functionality

5. **Test mcp-epic FHIR De-identification** (5 min)
   - HIPAA Safe Harbor validation, 18 identifiers removed
   - **Output:** 100% de-identification, audit log created

6. **Load Test - Concurrent Analysis Requests** (5-10 min)
   - 5 concurrent users, response time, throughput, error rate
   - **Output:** p95 <30s, 0% error rate, auto-scaling validated

7. **Integration Test - Multi-Server Workflow** (10-15 min)
   - Clinical + Multi-omics + Spatial end-to-end
   - **Output:** All 3 servers contribute, findings concordant

8. **Error Handling and Recovery Test** (5 min)
   - Invalid inputs, missing files, malformed data, timeouts
   - **Output:** Graceful errors, helpful messages, server recovers

9. **Data Quality and Integrity Check** (10 min)
   - VCF format, multi-omics completeness, spatial validation
   - **Output:** Quality report with preprocessing recommendations

### Additional Developer Prompts (2 Prompts)

10. **Cost and Token Usage Monitoring** (2 min)
    - Track Claude API usage across TEST_1-6 workflow
    - **Output:** ~58K tokens, $0.37 total cost

11. **Audit Log Verification** (3 min)
    - HIPAA compliance, 10-year retention, log immutability
    - **Output:** All access logged, retention configured

---

## Educational Prompts (10 Prompts)

**Audience:** Students, educators, new users
**Focus:** Learning bioinformatics concepts, hands-on tutorials

### Beginner Level (3 Prompts)

1. **Introduction to MCP** (5 min)
   - Explore servers and tools, understand architecture
   - **Learning:** MCP server structure, data types, server relationships

2. **First Analysis - Simple Gene Expression** (5 min)
   - Query MKI67 expression in PatientOne data
   - **Learning:** Load data, compare regions, interpret biology

3. **Understanding p-values and Statistical Significance** (10 min)
   - Mann-Whitney U test, multiple testing correction
   - **Learning:** p-values, FDR, statistical vs biological significance

### Intermediate Level (3 Prompts)

4. **Pathway Enrichment Analysis Tutorial** (15 min)
   - Find biological pathways from gene lists
   - **Learning:** Enrichment concept, fold enrichment, clinical implications

5. **Batch Effects and Data Preprocessing** (15 min)
   - Visualize batch effects with PCA, apply ComBat correction
   - **Learning:** Batch effects, PCA, quality control, validation

6. **Multi-Omics Integration** (20 min)
   - Combine RNA, protein, phospho evidence with Stouffer's method
   - **Learning:** Multi-modal evidence, meta-analysis, confidence levels

### Advanced Level (4 Prompts)

7. **Spatial Transcriptomics - Tumor Microenvironment** (20 min)
   - Analyze spatial patterns, Moran's I, immune exclusion
   - **Learning:** Spatial biology, autocorrelation, treatment implications

8. **Clinician-in-the-Loop Workflow** (15 min)
   - Understand clinical validation, quality gates, approval process
   - **Learning:** AI oversight, safety, human accountability

9. **Cost-Effectiveness Analysis** (10 min)
   - Calculate cost savings, ROI, value-based healthcare
   - **Learning:** Healthcare economics, ROI, balancing cost vs outcomes

10. **Final Project - Complete Patient Analysis** (45 min)
    - Capstone exercise: analyze hypothetical breast cancer patient
    - **Learning:** Apply all concepts in integrated workflow

---

## Prompt Discovery Inventory

**Location:** [PROMPT_INVENTORY.md](PROMPT_INVENTORY.md)

Comprehensive inventory of all 50+ prompts discovered in the repository:
- PatientOne clinical workflow prompts (6 TEST prompts)
- MCP server-specific examples (multiomics, spatialtools, fgbio)
- UI-specific prompts (Streamlit, Jupyter)
- Quickstart workflow examples
- Summary statistics by audience, server, complexity

---

## How to Use This Library

### 1. Choose by Audience

**I am a...**
- **Funder/Decision-Maker:** Start with [funder-demo-prompts.md](funder-demo-prompts.md)
- **Clinician/Researcher:** Start with [clinical-workflow-prompts.md](clinical-workflow-prompts.md)
- **Developer/QA:** Start with [developer-test-prompts.md](developer-test-prompts.md)
- **Student/Educator:** Start with [educational-prompts.md](educational-prompts.md)

### 2. Choose by Time Available

**Quick (<5 min):**
- Funder: Executive Summary Demo
- Clinical: BRCA Testing Prompt
- Developer: List Tools Test
- Educational: First Analysis

**Medium (5-15 min):**
- Funder: ROI Analysis
- Clinical: Spatial Tumor Analysis
- Developer: Integration Test
- Educational: Pathway Enrichment Tutorial

**Comprehensive (15-45 min):**
- Funder: Multi-Cancer Extensibility
- Clinical: Complete PatientOne Workflow (6 prompts)
- Developer: Data Quality Check
- Educational: Final Project

### 3. Choose by Learning Goal

**Understand Platform:**
- Educational Prompt 1 (MCP Architecture)
- Developer Prompt 1 (List Tools)

**Learn Analysis Methods:**
- Educational Prompts 3-6 (Statistics, Pathways, Batch Correction, Multi-Omics)
- Clinical Prompts 1-5 (Complete Workflow)

**Validate System:**
- Developer Prompts 1-9 (Server Validation)
- Clinical Prompt 6 (CitL Validation)

**Show ROI:**
- Funder Prompts 1-7 (All demonstrate value)

---

## Prompt Statistics

### By Audience
- **Funders:** 7 prompts
- **Clinicians:** 15 prompts
- **Developers:** 11 prompts
- **Educators:** 10 prompts
- **Total Curated:** 43 prompts

### By Complexity
- **Low (<5 min):** 15 prompts
- **Medium (5-15 min):** 20 prompts
- **High (15-45 min):** 8 prompts

### By Expected Output
- **Business Metrics:** 7 prompts (funders)
- **Clinical Reports:** 15 prompts (clinicians)
- **Technical Validation:** 11 prompts (developers)
- **Tutorial Results:** 10 prompts (educators)

### By MCP Servers Used
- **mcp-multiomics:** 18 prompts
- **mcp-spatialtools:** 22 prompts
- **mcp-fgbio:** 8 prompts
- **mcp-epic/mockepic:** 12 prompts
- **mcp-openimagedata:** 6 prompts
- **mcp-deepcell:** 6 prompts
- **Multiple servers:** 25 prompts

---

## Related Resources

- **[Main README](../../README.md)** - Repository overview
- **[PROMPT_INVENTORY.md](PROMPT_INVENTORY.md)** - Discovery of existing prompts
- **[For Developers Guide](../../for-developers/README.md)** - Developer documentation
- **[For Researchers Guide](../../for-researchers/README.md)** - Research workflows
- **[PatientOne Guide](../test-docs/patient-one-scenario/README.md)** - Complete case study
- **[Server Documentation](../../../servers/README.md)** - MCP server reference

---

## Contributing New Prompts

Have a useful prompt? Add it to the library!

**Guidelines:**
1. Choose appropriate file (funder/clinical/developer/educational)
2. Follow existing format (Time, Learning Objective, Expected Output)
3. Include both prompt text and expected output
4. Test with PatientOne data
5. Submit PR with documentation

**Prompt Template:**
```
### Prompt X: [Clear Title]

**Time:** X minutes | **Complexity:** Low/Medium/High | **Output:** [Output type]

[Prompt text with clear instructions]

**Expected Output:**
[Detailed expected output with metrics/examples]

**Learning Objectives:** (for educational)
- ✓ Objective 1
- ✓ Objective 2
```

---

**Document Version:** 1.0
**Date:** 2026-01-16
**Total Prompts:** 43 curated + 50+ discovered = 90+ prompts
**Status:** Complete prompt library ready for use
