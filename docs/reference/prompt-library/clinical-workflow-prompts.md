# Clinical Workflow Prompts

**Audience:** Oncologists, Clinical Researchers, Bioinformaticians
**Purpose:** Complete patient analysis workflows for precision oncology
**Time per prompt:** 5-45 minutes
**Expected output:** Clinical reports, treatment recommendations, scientific analysis

---

## Complete PatientOne Workflow (6 Prompts)

These prompts follow the TEST_1 through TEST_6 workflow for comprehensive multi-modal analysis.

### Prompt 1: Clinical Data and Genomic Profiling

**Time:** 5-10 minutes | **Complexity:** Medium | **Output:** Clinical summary + genomic profile

```
For patient PAT001-OVC-2025 with Stage IV high-grade serous ovarian carcinoma:

**PART 1: Clinical Data (use mcp-mockepic)**
1. Retrieve patient demographics:
   - Age, family history of cancer
   - BRCA1 germline mutation status
   - Prior treatment history

2. Analyze CA-125 tumor marker trends:
   - Baseline at diagnosis
   - Response to initial platinum therapy
   - Current levels (evidence of resistance?)

**PART 2: Genomic Analysis (use mcp-fgbio, mcp-tcga)**
3. Parse somatic variants from VCF:
   - TP53 mutations (expect R175H hotspot)
   - PIK3CA mutations (expect E545K)
   - PTEN alterations (expect LOH)
   - Copy number: MYC, CCNE1, AKT2 amplifications

4. Compare to TCGA-OV cohort:
   - What molecular subtype? (C1 immunoreactive vs C2 differentiated)
   - Prognosis for BRCA1+/TP53+ Stage IV HGSOC?
   - Common resistance pathways?

**Output Format:**
- Patient summary (demographics, genetic risk factors)
- CA-125 trajectory with resistance interpretation
- Key somatic mutations with clinical significance
- TCGA subtype and prognosis
- Initial treatment considerations
```

**Expected Output:**
- Patient: Sarah Anderson, 58yo, BRCA1 germline mutation
- CA-125: 1456 → 22 → 389 U/mL (platinum resistance confirmed)
- Mutations: TP53 R175H, PIK3CA E545K, PTEN LOH
- TCGA subtype: C1/C2 (poor prognosis)
- Key insight: PI3K pathway activated → potential PI3K inhibitor target

---

### Prompt 2: Multi-Omics Resistance Mechanisms

**Time:** 15-25 minutes | **Complexity:** High | **Output:** Resistance pathway analysis

```
For patient PAT001-OVC-2025, analyze platinum resistance using multi-omics data:

**STEP 0: Data Preprocessing (CRITICAL)**
1. Validate data quality (use mcp-multiomics):
   - Check for batch effects in proteomics
   - Identify missing value patterns
   - Detect outlier samples

2. Preprocess data:
   - Apply ComBat batch correction
   - KNN imputation for missing values
   - Quantile normalization

**STEP 1: Multi-Omics Integration**
3. Integrate preprocessed data:
   - RNA-seq (14 samples: 7 resistant, 7 sensitive PDX)
   - Proteomics (TMT data)
   - Phosphoproteomics

**STEP 2: Resistance Gene Identification**
4. Focus on key resistance pathways:
   - PI3K/AKT: PIK3CA, AKT1, MTOR, PTEN
   - Drug efflux: ABCB1 (MDR1)
   - Anti-apoptotic: BCL2L1
   - Tumor suppressor: TP53

5. Run Stouffer's meta-analysis:
   - Combine NOMINAL p-values across modalities
   - Apply FDR correction AFTER combination
   - Identify significant genes (q < 0.05)

**STEP 3: Upstream Regulator Prediction**
6. Predict therapeutic targets:
   - Activated kinases (AKT1, MTOR, PI3K)
   - Inhibited tumor suppressors (TP53)
   - Drug recommendations with evidence

**Output Format:**
- Preprocessing summary (batch correction metrics)
- Gene-level results table (7 genes: FC, Z-score, q-value)
- Upstream regulators (kinases, TFs, drug targets)
- Pathway interpretation (PI3K/AKT activation status)
- Therapeutic recommendations with evidence levels
```

**Expected Output:**
- Batch effects REDUCED: PC1-batch 0.82 → 0.15
- Resistance genes: PIK3CA (+2.3), AKT1 (+2.1), MTOR (+1.9), ABCB1 (+2.5) all q<0.001
- Tumor suppressors: PTEN (-2.1), TP53 (-1.5) both q<0.005
- Upstream regulators: AKT1/MTOR/PI3K activated (Z>2.5)
- Drug targets: Alpelisib (PI3K inhibitor), Capivasertib (AKT inhibitor)

---

### Prompt 3: Spatial Tumor Microenvironment Analysis

**Time:** 10-15 minutes | **Complexity:** Medium | **Output:** Spatial analysis report

```
For patient PAT001-OVC-2025, analyze spatial transcriptomics (10x Visium):

**Data:** 900 spots across 6 tissue regions
- tumor_core, tumor_proliferative, tumor_interface
- stroma_immune, stroma, necrotic_hypoxic

**Analysis Steps:**

1. **Spatial Structure Assessment:**
   - How many spots per region?
   - Spatial distribution of regions?

2. **Key Gene Expression by Region (focus on 8 genes):**
   - Proliferation: MKI67, PCNA
   - Resistance: PIK3CA, AKT1, ABCB1
   - Immune: CD3D, CD8A, CD68

3. **Spatial Patterns:**
   - Where is proliferation highest?
   - Are resistance markers heterogeneous?
   - Are immune cells excluded from tumor?

4. **Generate Visualizations:**
   - Spatial heatmap (6 key genes)
   - Region composition bar chart
   - Gene×region expression heatmap (8×6)
   - Spatial autocorrelation plot

**Output Format:**
- Spatial structure (900 spots, 6 regions with counts)
- Gene expression patterns by region
- Spatial findings (proliferation, resistance, immune exclusion)
- Visualizations (4 figures)
- Clinical interpretation (tumor microenvironment classification)
```

**Expected Output:**
- Structure: 900 spots (tumor_core: 69, tumor_proliferative: 124, stroma_immune: 212, etc.)
- Proliferation: HIGH in tumor_proliferative (MKI67/PCNA)
- Resistance: Concentrated in tumor regions (heterogeneous)
- Immune cells: EXCLUDED from tumor (located in stroma_immune)
- TME classification: "COLD" tumor (poor immunotherapy candidate)

---

### Prompt 4: Histopathology and Immunofluorescence

**Time:** 20-40 minutes | **Complexity:** High | **Output:** Imaging analysis report

```
For patient PAT001-OVC-2025, analyze histology and immunofluorescence images:

**Image Files:**
1. H&E histology (PAT001_tumor_HE_20x.tiff)
2. CD8 immunofluorescence (PAT001_tumor_IF_CD8.tiff)
3. Ki67 immunofluorescence (PAT001_tumor_IF_KI67.tiff)
4. Multiplex IF (PAT001_tumor_multiplex_IF_TP53_KI67_DAPI.tiff)

**Analysis Steps:**

1. **H&E Morphology (use mcp-openimagedata):**
   - Estimate tumor cellularity (%)
   - Identify necrotic regions (%)
   - Describe tissue architecture
   - Consistent with HGSOC diagnosis?

2. **CD8 T Cell Infiltration (use mcp-openimagedata + mcp-deepcell + mcp-cell-classify):**
   - Segment and quantify CD8+ cells (cells/mm²)
   - Spatial distribution (peripheral vs infiltrating)
   - Immune phenotype (hot/warm/cold)

3. **Ki67 Proliferation (use mcp-openimagedata + mcp-deepcell + mcp-cell-classify):**
   - Calculate Ki67 index (% positive cells)
   - Distribution pattern (uniform vs heterogeneous)
   - Proliferation level classification

4. **Multiplex IF Phenotyping (use mcp-deepcell for segmentation, mcp-cell-classify for classification):**
   - Segment cells based on DAPI (mcp-deepcell)
   - Quantify per-cell marker intensities (mcp-deepcell)
   - Classify: TP53+, Ki67+, TP53+/Ki67+ double-positive (mcp-cell-classify)
   - Are TP53-mutant cells proliferating?

5. **Generate Visualizations:**
   - H&E annotation (necrosis, cellularity)
   - CD8 segmentation overlay + density heatmap
   - Ki67 segmentation overlay + proliferation map
   - Multiplex composite (TP53/Ki67/DAPI channels)
   - Cell phenotype segmentation

**Output Format:**
- H&E summary (cellularity, necrosis, architecture)
- CD8 analysis (density, distribution, immune phenotype)
- Ki67 index (%, pattern, level)
- Multiplex results (cell counts, phenotypes, correlation)
- Overall imaging interpretation
```

**Expected Output:**
- H&E: 70-80% cellularity, 15-20% necrosis, HGSOC morphology confirmed
- CD8: LOW infiltration (5-15 cells/mm²), peripheral clustering, COLD phenotype
- Ki67: HIGH index (45-55%), heterogeneous distribution
- Multiplex: 65-75% TP53+, 45-55% Ki67+, 40-50% double-positive
- Interpretation: Aggressive, immune-excluded tumor

---

### Prompt 5: Integrated Multi-Modal Analysis

**Time:** 5-10 minutes | **Complexity:** Medium | **Output:** Clinical report + recommendations

```
Synthesize findings from Tests 1-4 for patient PAT001-OVC-2025:

**Reference Results:**
- Clinical/Genomic: BRCA1+, TP53 R175H, PIK3CA E545K, CA-125 resistance pattern
- Multi-omics: PI3K/AKT pathway activated, ABCB1 upregulated
- Spatial: Immune exclusion, resistance markers in tumor regions
- Imaging: High Ki67, low CD8 infiltration

**Integration Tasks:**

1. **Identify Primary Resistance Mechanisms (rank by evidence):**
   - Which mechanisms appear across ≥2 modalities?
   - Evidence strength (High/Medium/Low)?
   - Therapeutic implications?

2. **Multi-Modal Consistency Assessment:**
   - Create cross-reference table:
     * TP53 loss: Genomics + Imaging?
     * PI3K/AKT: Genomics + Multi-omics + Spatial?
     * High proliferation: Multi-omics + Spatial + Imaging?
     * Immune exclusion: Spatial + Imaging?

3. **Treatment Recommendations (Top 3):**
   - Targeted therapies with molecular evidence
   - Immunotherapy consideration (yes/no, why?)
   - Clinical trial opportunities
   - Expected efficacy for each recommendation

4. **Biomarkers for Monitoring:**
   - Molecular: Which genes/proteins to track?
   - Imaging: Ki67, CD8 infiltration trends?
   - Clinical: CA-125 trajectory, RECIST criteria

5. **Multi-Modal Visualization Synthesis:**
   - Generate integrated 4-panel figure:
     * Panel A: Spatial heatmap (resistance genes)
     * Panel B: H&E annotation (morphology)
     * Panel C: Multiplex IF (TP53/Ki67/CD8)
     * Panel D: Gene×region expression heatmap

**Output Format:**
- Executive summary (3-4 sentences)
- Resistance mechanisms (ranked 1-3 with evidence)
- Multi-modal consistency table
- Treatment recommendations (targeted, immuno, trials)
- Monitoring strategy (molecular, imaging, clinical)
- Integrated multi-modal figure with caption
```

**Expected Output:**
- Top mechanisms: (1) PI3K/AKT activation (4 modalities), (2) Drug efflux/ABCB1 (2 modalities), (3) Immune exclusion (2 modalities)
- Consistency: TP53 (genomics+imaging ✓), PI3K (genomics+multiomics+spatial ✓), Proliferation (multiomics+spatial+imaging ✓)
- Recommendations: (1) PI3K inhibitor+PARP inhibitor (high evidence), (2) MDR1 reversal (medium), (3) Checkpoint inhibitor (low - cold tumor)
- Monitoring: CA-125, circulating tumor DNA, Ki67 index

---

### Prompt 6: Clinician-in-the-Loop Validation

**Time:** 30-45 minutes (includes 20-30 min manual review) | **Complexity:** High | **Output:** Approved clinical report

```
Implement formal clinician review workflow for patient PAT001-OVC-2025:

**STEP 1: Generate Draft Report (automated, ~30 sec)**
- Consolidate findings from Tests 1-5
- Run automated quality checks:
  * Sample sizes adequate (≥30 spots/region)
  * FDR thresholds met (q<0.05)
  * Data completeness (>95%)
  * Cross-modal consistency verified
- Flag any quality issues

**STEP 2: Clinician Review (manual, 20-30 min)**
Oncologist reviews and validates:

1. **High-Level Decision:**
   - APPROVE / REVISE / REJECT
   - Rationale (2-3 sentences)

2. **Per-Finding Validation (10 top findings):**
   For each finding, mark CONFIRMED / UNCERTAIN / INCORRECT:
   - TP53 R175H mutation
   - PIK3CA E545K activation
   - BRCA1 germline pathogenic
   - ABCB1 upregulation (4.3× log2FC)
   - BCL2L1 anti-apoptotic signaling
   - Immune exclusion phenotype
   - VEGFA upregulation (hypoxia)
   - TP53+/Ki67+ proliferative cells
   - PI3K/AKT pathway activation
   - MDR1 drug efflux

3. **Guideline Compliance Check:**
   - NCCN alignment: ALIGNED / PARTIAL / NOT_ALIGNED
   - Institutional alignment: ALIGNED / PARTIAL / NOT_ALIGNED

4. **Treatment Recommendations Review:**
   For each recommendation, mark AGREE / DISAGREE:
   - PI3K inhibitor + PARP inhibitor
   - Anti-VEGF continuation
   - MDR1 reversal consideration
   - Immunotherapy (conditional)

5. **Digital Attestation:**
   - Reviewed all findings: ✓
   - Assessed compliance: ✓
   - Clinical judgment applied: ✓
   - Medical record acknowledgment: ✓

**STEP 3: Submit Review (automated, ~5 sec)**
- Validate review against JSON schema
- Generate digital signature (SHA-256)
- Log to audit trail (10-year retention)

**STEP 4: Finalize Report (automated, ~10 sec)**
- If APPROVED: Generate final clinical report
- Mark as "clinically_approved"
- Include attestation and signature
- Ready for tumor board presentation

**Output Format:**
- Draft report with quality checks
- Completed review form (6 sections)
- Signed review with digital signature
- Final approved clinical report (if APPROVE)
- Audit trail record
```

**Expected Output:**
- Quality checks: ALL PASS (sample sizes OK, FDR<0.05, completeness 99%)
- Decision: APPROVE (expected for PatientOne)
- Findings validated: 10/10 CONFIRMED
- Guideline compliance: NCCN ALIGNED, Institutional ALIGNED
- Treatments: 3 AGREE, 1 conditional
- Final report: CLINICALLY_APPROVED status
- Audit trail: Complete with 10-year retention

---

## Additional Clinical Use Cases (9 Prompts)

### Prompt 7: BRCA1/2 Testing and PARP Inhibitor Eligibility

**Time:** 5 minutes | **Complexity:** Low | **Output:** BRCA status + treatment eligibility

```
For patient [ID] with ovarian cancer:

1. Check BRCA1/2 germline mutation status (mcp-mockepic clinical history)
2. Check for somatic BRCA1/2 mutations (mcp-fgbio VCF analysis)
3. Assess homologous recombination deficiency (HRD) score
4. Determine PARP inhibitor eligibility:
   - BRCA1/2 germline: First-line maintenance (Category 1)
   - BRCA1/2 somatic: First-line maintenance (Category 1)
   - HRD+/BRCA-: Consider maintenance (Category 2A)
   - HRD-/BRCA-: Not recommended
5. List approved PARP inhibitors: Olaparib, Niraparib, Rucaparib
6. Provide NCCN guideline reference

Output: BRCA status, HRD score, PARP inhibitor eligibility, specific drug recommendations
```

---

### Prompt 8: Bevacizumab Continuation Assessment

**Time:** 5 minutes | **Complexity:** Low | **Output:** Anti-VEGF recommendation

```
For patient [ID] currently on bevacizumab:

1. Assess VEGFA expression (mcp-spatialtools spatial analysis)
2. Check for contraindications:
   - History of bowel perforation/fistula?
   - Recent surgery (<28 days)?
   - Uncontrolled hypertension?
   - Proteinuria (>2g/24hr)?
3. Evaluate response:
   - CA-125 trend (decreasing?)
   - Imaging response (RECIST)?
4. Recommendation:
   - Continue if VEGFA high + no contraindications + responding
   - Discontinue if progressing or contraindications present
5. Provide NCCN guideline reference

Output: VEGFA level, contraindications check, response status, continue/discontinue recommendation
```

---

### Prompt 9: Immunotherapy Eligibility (PD-L1, TMB, MSI)

**Time:** 10 minutes | **Complexity:** Medium | **Output:** Immunotherapy assessment

```
For patient [ID] with ovarian cancer:

1. **PD-L1 Assessment:**
   - Check PD-L1 expression (mcp-openimagedata IHC image or mcp-spatialtools)
   - Tumor proportion score (TPS): <1%, 1-49%, ≥50%

2. **Tumor Mutational Burden (TMB):**
   - Calculate from VCF (mcp-fgbio)
   - TMB-H threshold: ≥10 mutations/Mb

3. **Microsatellite Instability (MSI):**
   - Check MSI status (mcp-fgbio)
   - MSI-H or dMMR?

4. **Immune Infiltration:**
   - Quantify CD8+ T cells (mcp-deepcell)
   - Tumor microenvironment: Hot/Warm/Cold?

5. **Recommendation:**
   - Pembrolizumab if: PD-L1 CPS≥10 OR TMB-H OR MSI-H
   - Dostarlimab if: dMMR/MSI-H
   - NOT recommended if: PD-L1 negative + TMB low + MSS + cold TME

Output: PD-L1 score, TMB value, MSI status, CD8 count, immune phenotype, checkpoint inhibitor recommendation
```

---

### Prompt 10: Platinum Re-Challenge vs Alternative Chemotherapy

**Time:** 5-10 minutes | **Complexity:** Medium | **Output:** Chemotherapy recommendation

```
For patient [ID] with platinum-resistant ovarian cancer:

1. **Define Platinum Resistance:**
   - Progression-free interval: <6 months = resistant, 6-12 months = partially sensitive, >12 months = sensitive
   - Calculate PFI from clinical data (mcp-mockepic)

2. **Assess Resistance Mechanisms:**
   - ABCB1/MDR1 expression (mcp-multiomics or mcp-spatialtools)
   - Anti-apoptotic signaling (BCL2 family)
   - DNA repair pathway activation

3. **Treatment Options:**
   - If resistant (<6 months): Non-platinum agents (paclitaxel, pegylated liposomal doxorubicin, topotecan, gemcitabine)
   - If partially sensitive (6-12 months): Consider platinum re-challenge or non-platinum
   - If sensitive (>12 months): Platinum re-challenge recommended

4. **Combination Strategies:**
   - PARP inhibitor + chemotherapy (if BRCA+)
   - Bevacizumab + chemotherapy
   - PI3K inhibitor + chemotherapy (if PIK3CA mutant)

Output: PFI, resistance status, resistance mechanisms, recommended chemotherapy regimen
```

---

### Prompt 11: Clinical Trial Matching

**Time:** 10 minutes | **Complexity:** Medium | **Output:** Trial recommendations

```
For patient [ID] with Stage IV platinum-resistant HGSOC:

**Patient Molecular Profile:**
- BRCA1 germline mutation
- TP53 R175H somatic mutation
- PI3K/AKT pathway activation
- Immune exclusion phenotype

**Search Clinical Trials:**
1. **PI3K/AKT Inhibitor Trials:**
   - Alpelisib + Olaparib (NCT03740165 or similar)
   - Capivasertib + Olaparib
   - Eligibility: PIK3CA mutation or pathway activation + BRCA1/2 mutation

2. **PARP Inhibitor Trials:**
   - Novel PARP inhibitors
   - PARP + immune checkpoint combinations
   - Eligibility: BRCA1/2 or HRD+

3. **Immunotherapy Combination Trials:**
   - Checkpoint inhibitor + VEGF inhibitor
   - Checkpoint inhibitor + PARP inhibitor
   - Eligibility: Usually allows cold tumors if combination

4. **Novel Targeted Therapy:**
   - ATR inhibitors
   - WEE1 inhibitors
   - Eligibility: TP53 mutant

Provide:
- Top 3 matching trials with NCT numbers
- Eligibility criteria match assessment
- Trial phase and location
- Contact information
```

---

### Prompt 12: Circulating Tumor DNA (ctDNA) Monitoring

**Time:** 5 minutes | **Complexity:** Low | **Output:** ctDNA monitoring plan

```
For patient [ID] starting new treatment:

1. **Baseline ctDNA Panel:**
   - Track mutations found in tumor: TP53, PIK3CA, BRCA1, PTEN
   - Establish baseline ctDNA level (copies/mL)

2. **Monitoring Schedule:**
   - Cycle 2 (week 6): First on-treatment assessment
   - Cycle 4 (week 12): Confirm molecular response
   - Every 3 months during maintenance
   - At progression or symptoms

3. **Interpretation:**
   - Decreasing ctDNA: Molecular response (positive prognostic)
   - Stable ctDNA: Stable disease
   - Increasing ctDNA: Molecular progression (may precede imaging by 2-4 months)

4. **Action Plan:**
   - If ctDNA rises: Consider imaging, biopsy for resistance mechanisms
   - If new mutations appear: Update treatment plan

Output: Mutations to track, monitoring schedule, interpretation guide, action thresholds
```

---

### Prompt 13: Comprehensive Molecular Tumor Board Presentation

**Time:** 15 minutes | **Complexity:** High | **Output:** Tumor board presentation

```
Prepare molecular tumor board presentation for patient [ID]:

**Slide 1: Case Summary**
- Patient demographics
- Cancer type, stage, prior treatments
- Current status (progression, response, stable)

**Slide 2: Clinical Timeline**
- Diagnosis date
- Treatment history with responses
- CA-125 trajectory graph
- Current problem (resistance, progression)

**Slide 3: Genomic Landscape**
- Germline mutations (BRCA1/2, Lynch)
- Somatic mutations (TP53, PIK3CA, etc.)
- Copy number alterations
- Comparison to TCGA cohort

**Slide 4: Multi-Omics Resistance Profile**
- Activated pathways (PI3K/AKT, WNT, etc.)
- Drug resistance genes (ABCB1, BCL2L1)
- Upstream regulator prediction
- Pathway visualization

**Slide 5: Spatial Tumor Microenvironment**
- Spatial heatmap showing gene expression
- Immune infiltration patterns
- Proliferation zones
- Tumor heterogeneity assessment

**Slide 6: Histopathology & Imaging**
- H&E morphology
- Ki67 proliferation index
- CD8 infiltration quantification
- Multiplex IF phenotyping

**Slide 7: Integrated Findings**
- Multi-modal consistency table
- Primary resistance mechanisms (ranked)
- Druggable targets identified

**Slide 8: Treatment Recommendations**
- First-line: Drug + evidence + expected response
- Second-line: Alternative + rationale
- Clinical trial opportunities
- NCCN guideline compliance

**Slide 9: Monitoring Plan**
- Biomarkers to track
- Imaging schedule
- ctDNA monitoring
- Response assessment criteria

Output: 9-slide presentation with visualizations, suitable for oncology tumor board
```

---

### Prompt 14: Adverse Event Management - Hypersensitivity Reaction

**Time:** 3 minutes | **Complexity:** Low | **Output:** AE management plan

```
Patient [ID] developed hypersensitivity reaction during paclitaxel infusion:

1. **Grade Assessment:**
   - Grade 1: Mild transient rash, drug fever <38°C
   - Grade 2: Moderate rash, fever ≥38°C, dyspnea
   - Grade 3: Severe bronchospasm, hypotension, generalized urticaria
   - Grade 4: Anaphylaxis

2. **Immediate Management (if Grade 3-4):**
   - Stop infusion immediately
   - Epinephrine 0.3-0.5 mg IM
   - Diphenhydramine 50 mg IV
   - Hydrocortisone 100-200 mg IV
   - O2, IV fluids, monitor vitals

3. **Future Treatment:**
   - Grade 1-2: Consider re-challenge with aggressive premedication
   - Grade 3-4: Avoid paclitaxel, use alternative (docetaxel, carboplatin, PLD)

4. **Premedication Protocol (if re-challenge):**
   - Dexamethasone 20 mg PO 12h + 6h before
   - Diphenhydramine 50 mg IV 30 min before
   - H2 blocker (ranitidine 50 mg IV or famotidine 20 mg IV) 30 min before
   - Slow infusion rate (3-6 hours)

Output: AE grade, immediate management, re-challenge decision, premedication protocol
```

---

### Prompt 15: End-of-Treatment Summary and Survivorship Plan

**Time:** 10 minutes | **Complexity:** Low | **Output:** Survivorship care plan

```
For patient [ID] completing front-line treatment for ovarian cancer:

**Treatment Summary:**
1. Diagnosis: Stage, grade, histology
2. Treatment received:
   - Surgery: Type, date, outcome
   - Chemotherapy: Regimen, cycles, completion date
   - Targeted therapy: Drug, duration
3. Response: Complete/partial/stable/progressive
4. Residual disease: None / measurable

**Surveillance Plan:**
1. **Clinical Exams:**
   - Every 3 months for 2 years
   - Every 6 months for years 3-5
   - Annually thereafter

2. **CA-125:**
   - Every 3 months for 2 years (if elevated at diagnosis)
   - Note: Rising CA-125 alone not indication for treatment (GCIG criteria)

3. **Imaging:**
   - CT chest/abdomen/pelvis every 6 months for 2 years
   - Then annually or as clinically indicated
   - PET not routinely recommended

4. **Genetic Counseling:**
   - If BRCA1/2+: Family testing, cascade screening
   - Risk-reducing salpingo-oophorectomy for relatives

5. **Quality of Life:**
   - Neuropathy management (if present)
   - Fatigue assessment
   - Psychosocial support
   - Fertility preservation discussion (if applicable)

6. **Recurrence Signs:**
   - New/worsening abdominal pain or bloating
   - GI symptoms (obstruction, nausea)
   - Unexplained weight loss
   - **Instruct patient to report symptoms promptly**

Output: Treatment summary, surveillance schedule, genetic counseling referral, QOL plan
```

---

## Summary of Clinical Workflow Prompts

### By Clinical Phase:
- **Diagnosis & Molecular Profiling:** Prompts 1, 7
- **Treatment Planning:** Prompts 2, 3, 4, 5, 8, 9, 10, 11
- **Treatment Validation:** Prompt 6
- **Monitoring:** Prompt 12
- **Tumor Board:** Prompt 13
- **Adverse Events:** Prompt 14
- **Survivorship:** Prompt 15

### By Complexity:
- **Low (3-5 min):** Prompts 7, 8, 10, 12, 14, 15
- **Medium (5-15 min):** Prompts 1, 3, 5, 9, 11
- **High (15-45 min):** Prompts 2, 4, 6, 13

### By Output Type:
- **Clinical Reports:** Prompts 1, 2, 3, 4, 5, 6, 13
- **Treatment Recommendations:** Prompts 7, 8, 9, 10, 11
- **Monitoring Plans:** Prompts 12, 15
- **Acute Management:** Prompt 14

---

**Document Version:** 1.0
**Date:** 2026-01-16
**Target Audience:** Oncologists, Clinical Researchers, Bioinformaticians
