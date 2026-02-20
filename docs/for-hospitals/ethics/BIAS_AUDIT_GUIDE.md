# Bias Audit Checklist

**Purpose:** Practical checklist for conducting algorithmic bias audits in precision medicine workflows

**Target Audience:** Clinicians, Researchers, Quality Assurance Teams

**Time Required:** 3-5 days per audit

**Last Updated:** 2026-01-12

---

## Quick Start

This checklist guides you through three phases of bias auditing:
1. **Pre-Analysis** (Day 1): Prepare data and documentation
2. **During Analysis** (Days 2-3): Run audit script and analyze results
3. **Post-Analysis** (Days 4-5): Document findings and implement mitigations

**Prerequisites:**
- Access to de-identified patient data (≥100 patients per ancestry group)
- Python environment with audit script installed
- Reference dataset documentation
- Clinical domain expertise

---

## Pre-Analysis Checklist

### Data Preparation

- [ ] **Collect de-identified patient data**
  - Genomic data (VCF files, gene expression matrices)
  - Clinical data (FHIR resources, de-identified)
  - Spatial transcriptomics (expression matrices, coordinates)
  - Demographics (self-reported ancestry, NOT race-as-biology)

- [ ] **Verify sample size**
  - Minimum 100 patients per ancestry group for statistical power
  - If <100 for any group, note as limitation in report

- [ ] **Check de-identification compliance**
  - All 18 HIPAA identifiers removed
  - No dates more specific than year
  - Geographic subdivisions larger than state
  - Use HIPAA-compliant de-identification tool if needed

---

### Reference Dataset Documentation

- [ ] **Inventory all reference datasets used**
  - Genomics: ClinVar, COSMIC, gnomAD, GTEx, etc.
  - Spatial: SingleR, Human Cell Atlas, custom references
  - Clinical: Local EHR norms, MIMIC-IV, etc.

- [ ] **Document ancestry distribution for each dataset**
  - Create table: Dataset | European % | African % | Latino % | Asian % | Other %
  - Flag datasets with <10% representation for any major group
  - Example:
    ```
    Dataset    | European | African | Latino | Asian | Other | Flag
    -----------|----------|---------|--------|-------|-------|------
    ClinVar    | 70%      | 10%     | 10%    | 5%    | 5%    | ⚠️
    gnomAD     | 43%      | 21%     | 14%    | 9%    | 8%    | ✅
    GTEx       | 85%      | 10%     | 5%     | 0%    | 0%    | ⚠️⚠️
    All of Us  | 40%      | 25%     | 20%    | 8%    | 7%    | ✅
    ```

- [ ] **Check reference dataset versions**
  - gnomAD: v4.1 (latest as of Jan 2026)
  - ClinVar: Monthly updates (use current month)
  - GTEx: v8 (check if v9 released)
  - All of Us: Quarterly releases

- [ ] **Identify datasets requiring updates**
  - If reference is >1 year old, check for updates
  - If newer diverse dataset available (e.g., All of Us), consider switching

---

### Ethical Review

- [ ] **Confirm appropriate use of ancestry data**
  - ✅ Appropriate: Ancestry for genomic variant interpretation
  - ❌ Inappropriate: Race as biological proxy (e.g., kidney function)
  - ✅ Appropriate: Socioeconomic ZIP code for resource access
  - ❌ Inappropriate: Insurance status influencing clinical recommendations

- [ ] **Verify IRB approval (if applicable)**
  - Bias audit on existing de-identified data: May not require new IRB
  - Check institutional policy
  - Document IRB determination

- [ ] **Confirm no PHI in analysis**
  - Audit uses aggregate statistics only
  - No individual patient re-identification
  - Stratification by demographic groups is permitted under HIPAA

---

### Workflow Documentation

- [ ] **Document current workflow**
  - Which tools are used (mcp-spatialtools, mcp-multiomics, etc.)
  - Which reference datasets feed each tool
  - Which features/variables are inputs to algorithms
  - Where ancestry/demographic data is used

- [ ] **Identify potential bias entry points**
  - Data collection: Are certain populations undersampled?
  - Feature selection: Are demographic proxies used (ZIP code, insurance)?
  - Model training: Was training data diverse?
  - Clinical interpretation: Are clinicians aware of limitations?

- [ ] **Review past bias incidents (if any)**
  - Have clinicians reported suspected bias?
  - Were certain ancestries disproportionately affected?
  - What mitigations were implemented?

---

### Tool Setup

- [ ] **Install Python dependencies**
  ```bash
  cd /path/to/spatial-mcp
  pip install -r requirements.txt
  ```

- [ ] **Verify audit script is runnable**
  ```bash
  python infrastructure/audit/audit_bias.py --help
  ```

- [ ] **Test on sample data**
  ```bash
  python infrastructure/audit/audit_bias.py \
    --workflow patientone \
    --genomics-data data/genomics/sample.vcf \
    --output /tmp/test_audit.html
  ```

- [ ] **Prepare output directory**
  ```bash
  mkdir -p reports/bias_audits/$(date +%Y)
  ```

---

### Stakeholder Communication

- [ ] **Notify clinical team**
  - Inform clinicians that bias audit is being conducted
  - Request input on suspected biases or edge cases
  - Schedule review meeting for audit results

- [ ] **Notify compliance team**
  - Hospital quality assurance
  - Privacy officer (confirm de-identification)
  - Legal (if needed for regulatory submission)

- [ ] **Set expectations**
  - Audit will take 3-5 days
  - Results may require workflow changes
  - Mitigations will be implemented before production deployment

---

## During Analysis Checklist

### Data Representation Analysis

- [ ] **Calculate ancestry distribution of reference datasets**
  - Use documented percentages from Pre-Analysis phase
  - Compare to US population baseline:
    * European: 60%
    * Hispanic/Latino: 18%
    * African/Black: 13%
    * Asian: 6%
    * Native American: 1%
    * Multiracial: 2%

- [ ] **Flag underrepresented groups**
  - <10% representation = HIGH RISK
  - <5% representation = CRITICAL RISK
  - 0% representation = UNACCEPTABLE

- [ ] **Document findings**
  - Create "Reference Dataset Diversity Report" table
  - For each dataset, note risk level
  - Prioritize datasets for replacement (critical risk first)

---

### Algorithm Fairness Testing

- [ ] **Stratify test data by ancestry**
  - European/White
  - African/Black/African American
  - Hispanic/Latino
  - Asian (specify: East Asian, South Asian, etc.)
  - Native American
  - Multiracial/Other

- [ ] **Calculate performance metrics per group**
  - Accuracy, sensitivity, specificity
  - Positive predictive value (PPV), negative predictive value (NPV)
  - For continuous outcomes: RMSE, MAE
  - For rankings: AUC-ROC

- [ ] **Test for statistical significance**
  - Use chi-square test for categorical outcomes
  - Use t-test for continuous outcomes
  - Flag disparities >10% as requiring review
  - Flag disparities >20% as critical

- [ ] **Calculate fairness metrics**
  - **Demographic parity:** P(ŷ=1|A=a) equal across groups?
  - **Equalized odds:** Equal TPR and FPR across groups?
  - **Calibration:** Predicted probabilities match observed rates?
  - **Clinical utility parity:** Equal benefit (QALYs, lives saved)?

- [ ] **Document fairness violations**
  - Which metric is violated?
  - Which groups are affected?
  - What is the magnitude of disparity?
  - Example:
    ```
    Metric: Equalized Odds (Sensitivity)
    European: 92%
    African: 80%
    Disparity: 12 percentage points ⚠️
    Conclusion: Requires mitigation
    ```

---

### Output Stratification Analysis

- [ ] **Stratify treatment recommendations by ancestry**
  - Count aggressive vs. conservative treatments by group
  - Run chi-square test for independence
  - Flag if p < 0.05 (statistically significant disparity)

- [ ] **Analyze confidence scores by group**
  - Calculate mean confidence per ancestry
  - Flag if >10% difference between groups
  - Lower confidence for underrepresented groups is EXPECTED (indicates uncertainty)
  - Example:
    ```
    Group      | Mean Confidence | Interpretation
    -----------|-----------------|----------------
    European   | 0.92            | High confidence
    African    | 0.78            | Medium confidence ⚠️ (14% lower due to limited reference data)
    Asian      | 0.70            | Medium confidence ⚠️⚠️ (22% lower)
    ```

- [ ] **Check for demographic proxy influence**
  - Identify proxy features: ZIP code, insurance type, language preference
  - Check feature importance scores
  - Run model with/without proxies
  - If performance changes >5%, proxies are influential (BAD)
  - Example:
    ```
    Feature          | Importance | Flag
    -----------------|------------|------
    BRCA1 status     | 0.45       | ✅ (clinical)
    CA-125 level     | 0.32       | ✅ (clinical)
    ZIP code         | 0.08       | ⚠️⚠️ (proxy - should be 0)
    Insurance type   | 0.05       | ⚠️ (proxy - should be 0)
    ```

- [ ] **Stratify by socioeconomic indicators (if appropriate)**
  - Medicaid vs. commercial insurance (should NOT affect clinical recommendations)
  - Urban vs. rural (should only affect resource access, NOT treatment)
  - Education level (NOT collected in most clinical workflows)

---

### Run Audit Script

- [ ] **Execute bias audit script**
  ```bash
  python infrastructure/audit/audit_bias.py \
    --workflow patientone \
    --genomics-data data/genomics/patient_variants.vcf \
    --clinical-data data/fhir/patients_deidentified.json \
    --spatial-data data/spatial/spatial_expression.csv \
    --demographics data/demographics/ancestry_labels.csv \
    --output reports/bias_audit_$(date +%Y-%m-%d).html \
    --min-representation 0.10 \
    --max-disparity 0.10
  ```

- [ ] **Verify script completes without errors**
  - Check for Python exceptions
  - Verify all input files loaded correctly
  - Confirm output files generated

- [ ] **Review automated findings**
  - Open HTML report in browser
  - Read executive summary (Pass/Fail)
  - Review flagged issues
  - Check recommendations

---

### Manual Validation

- [ ] **Spot-check algorithmic results**
  - Manually review 10-20 cases per ancestry group
  - Do predictions make clinical sense?
  - Are confidence scores appropriate?
  - Are disclaimers present for underrepresented groups?

- [ ] **Validate stratified metrics**
  - Re-calculate fairness metrics manually for 1 group
  - Verify automated script is correct
  - Cross-check with clinical domain knowledge

- [ ] **Review edge cases**
  - Patients with rare variants
  - Patients from underrepresented ancestries
  - Patients with conflicting data (e.g., discordant ClinVar classifications)

- [ ] **Check for data quality issues**
  - Are there missing values?
  - Are ancestry labels accurate (self-reported preferred)?
  - Are there sample size imbalances?

---

### Generate Visualizations

- [ ] **Create fairness metrics bar chart**
  - X-axis: Ancestry groups
  - Y-axis: Sensitivity (or other metric)
  - Highlight disparities >10%

- [ ] **Create reference dataset diversity heatmap**
  - Rows: Datasets (ClinVar, gnomAD, GTEx, etc.)
  - Columns: Ancestries (European, African, Latino, Asian, Other)
  - Color: Green (≥20%), Yellow (10-20%), Red (<10%)

- [ ] **Create confidence score distribution plot**
  - Separate distributions per ancestry group
  - Overlay to show differences

- [ ] **Include visualizations in HTML report**

---

## Post-Analysis Checklist

### Document Findings

- [ ] **Write executive summary**
  - Overall assessment: PASS / CONDITIONAL PASS / FAIL
  - Risk level: LOW / MEDIUM / HIGH / CRITICAL
  - Number of findings: X warnings, Y critical issues
  - Recommended action: Deploy / Deploy with mitigations / Do not deploy

- [ ] **Categorize findings by severity**
  - **Critical:** >20% disparity, 0% representation, proxy features used
  - **High:** 10-20% disparity, <5% representation
  - **Medium:** 5-10% disparity, <10% representation
  - **Low:** <5% disparity, disclaimers missing

- [ ] **Provide specific examples**
  - "African patients have 12% lower sensitivity (80% vs. 92% for European)"
  - "GTEx has 0% Asian representation → flagged in all reports"
  - "ZIP code has 8% feature importance → removed from model"

- [ ] **Map findings to mitigation strategies**
  - For each finding, specify mitigation (see Mitigation section below)

---

### Risk Assessment

- [ ] **Evaluate overall bias risk**
  - Genomics risk: LOW / MEDIUM / HIGH
  - Gene expression risk: LOW / MEDIUM / HIGH
  - Clinical risk: LOW / MEDIUM / HIGH
  - Spatial transcriptomics risk: LOW / MEDIUM / HIGH
  - Multiomics risk: LOW / MEDIUM / HIGH

- [ ] **Determine deployment readiness**
  - If all risks LOW → Deploy
  - If any risk MEDIUM → Deploy with mitigations + quarterly monitoring
  - If any risk HIGH → Implement mitigations first, re-audit
  - If any risk CRITICAL → Do not deploy until resolved

- [ ] **Estimate mitigation effort**
  - Days to implement mitigations
  - Resources required (personnel, compute, data)
  - Timeline for re-audit

- [ ] **Document limitations**
  - Known biases that cannot be fully mitigated
  - Reference datasets that cannot be replaced (no diverse alternative)
  - Acknowledge in reports and patient communications

---

### Implement Mitigations

**For Underrepresented Reference Data:**

- [ ] **Add diverse reference datasets**
  - Replace ClinVar → ClinVar + gnomAD + All of Us
  - Replace GTEx → GTEx + TOPMed + Human Cell Atlas
  - Document in workflow

- [ ] **Add ancestry caveats to reports**
  - Example: "GTEx has 0% Asian representation. Results may not accurately reflect normal variation in Asian ancestries."
  - Display prominently in clinician reports

- [ ] **Adjust thresholds for underrepresented groups**
  - Increase fold-change cutoff (log2FC > 2.0 instead of 1.5)
  - Widen confidence intervals
  - Require manual review

---

**For Algorithmic Fairness Violations:**

- [ ] **Retrain with fairness constraints**
  - Use fairlearn library for fairness-aware training
  - Constrain demographic parity or equalized odds
  - Validate improved fairness metrics

- [ ] **Post-hoc calibration by group**
  - Train separate calibration models per ancestry
  - Validate calibration curves

- [ ] **Adjust decision thresholds by group**
  - If sensitivity is lower for African patients, lower threshold for that group
  - Document rationale

---

**For Demographic Proxy Features:**

- [ ] **Remove proxy features**
  - Drop ZIP code, insurance type, language preference from clinical models
  - Retain only for resource allocation (not clinical decisions)

- [ ] **Retrain model without proxies**
  - Verify performance does not degrade significantly (<5% accuracy loss acceptable)
  - If performance degrades >5%, investigate why (may indicate confounding)

- [ ] **Add audit log for proxy usage**
  - If proxies must be used (e.g., resource allocation), log every use
  - Enable review of inappropriate proxy usage

---

**For Confidence Score Disparities:**

- [ ] **Implement ancestry-aware confidence scoring**
  - Reduce confidence by X% if patient ancestry has <10% representation
  - Flag for manual review if confidence <0.70

- [ ] **Add uncertainty quantification**
  - Display confidence intervals, not just point estimates
  - Example: "Risk: 65% (95% CI: 55%-75%)"

---

### Clinical Validation

- [ ] **Schedule review meeting with clinical team**
  - Present audit findings
  - Discuss clinical significance of disparities
  - Get input on mitigations

- [ ] **Obtain sign-off from experts**
  - Medical geneticist (genomics findings)
  - Oncologist (treatment recommendations)
  - Bioinformatician (technical validation)
  - Bioethicist (fairness assessment)

- [ ] **Document expert recommendations**
  - Do experts agree with audit findings?
  - Are mitigations sufficient?
  - Are there additional concerns?

- [ ] **Decide on deployment**
  - Go / No-Go decision
  - If Go: Implement mitigations first or concurrent monitoring?
  - If No-Go: Timeline for re-audit

---

### Update Documentation

- [ ] **Update HIPAA compliance docs**
  - Add reference to bias audit in §8 (new section)
  - Cross-link to ETHICS_AND_BIAS.md

- [ ] **Update Operations Manual**
  - Add bias audit to monthly compliance checklist
  - Include in incident response procedures

- [ ] **Update Admin Guide**
  - Add "Bias Audit" section with instructions
  - Document how to run audit script

- [ ] **Update workflow documentation**
  - Note which diverse reference datasets are used
  - Document ancestry caveats
  - Update tool descriptions with bias mitigations

---

### Archive Audit Report

- [ ] **Store audit report**
  ```bash
  mkdir -p reports/bias_audits/$(date +%Y)
  mv reports/bias_audit_$(date +%Y-%m-%d).html \
     reports/bias_audits/$(date +%Y)/$(date +%Y-%m-%d)_initial_audit/
  ```

- [ ] **Generate compliance summary**
  - One-page PDF for compliance officers
  - Key findings, risk level, mitigations

- [ ] **Set retention policy**
  - Keep for 10 years (HIPAA alignment)
  - Include in annual compliance reviews

- [ ] **Update audit log**
  ```json
  {
    "audit_date": "2026-01-12",
    "auditor": "Jane Doe, PhD",
    "workflow": "PatientOne Complete Analysis",
    "risk_level": "MEDIUM",
    "findings": 2,
    "mitigations_implemented": 3,
    "next_audit": "2026-04-12"
  }
  ```

---

### Schedule Next Audit

- [ ] **Determine next audit date**
  - Quarterly audits: 3 months from now
  - Triggered audits: After workflow changes
  - Annual comprehensive audit: 12 months from now

- [ ] **Set calendar reminders**
  - 2 weeks before: Prepare data
  - 1 week before: Review checklist
  - Audit date: Run script and analyze

- [ ] **Assign responsibilities**
  - Who will run the audit?
  - Who will review results?
  - Who will implement mitigations?

- [ ] **Update audit schedule**
  ```markdown
  ## Bias Audit Schedule
  - 2026-01-12: Initial Audit (COMPLETE)
  - 2026-04-12: Quarterly Audit Q2 (SCHEDULED)
  - 2026-07-12: Quarterly Audit Q3 (SCHEDULED)
  - 2026-10-12: Quarterly Audit Q4 (SCHEDULED)
  - 2027-01-12: Annual Comprehensive Audit (SCHEDULED)
  ```

---

### Communicate Results

- [ ] **Present to clinical team**
  - Summarize findings (non-technical language)
  - Explain mitigations
  - Answer questions

- [ ] **Notify compliance officers**
  - Share audit report
  - Confirm compliance with ethical AI standards

- [ ] **Update stakeholders**
  - Hospital administration
  - Quality assurance
  - Privacy office

- [ ] **External communication (if required)**
  - FDA (if part of pre-market submission)
  - IRB (if research use)
  - Publications (if publishing methodology)

---

## Common Pitfalls to Avoid

### Data Issues

- ❌ **Small sample sizes:** <100 patients per group → Insufficient statistical power
  - ✅ Solution: Pool data across institutions or extend data collection period

- ❌ **Mislabeled ancestry:** Administrative race vs. self-reported ancestry
  - ✅ Solution: Use self-reported ancestry from FHIR Patient.extension (US Core)

- ❌ **Survivorship bias:** Only analyzing patients who completed treatment
  - ✅ Solution: Include all patients, stratify by completion status

---

### Analysis Issues

- ❌ **Cherry-picking metrics:** Only reporting fairness metrics that look good
  - ✅ Solution: Report all standard metrics (demographic parity, equalized odds, calibration)

- ❌ **Ignoring small disparities:** "5% difference isn't clinically meaningful"
  - ✅ Solution: Aggregate small disparities across population can have large impact

- ❌ **Post-hoc rationalization:** Finding excuses for bias instead of mitigating
  - ✅ Solution: Implement mitigations first, then re-assess

---

### Mitigation Issues

- ❌ **Band-aid solutions:** Adding disclaimers without addressing root cause
  - ✅ Solution: Replace biased reference datasets, retrain models

- ❌ **Over-correction:** Enforcing demographic parity when base rates differ
  - ✅ Solution: Use appropriate fairness metric (calibration for risk scores)

- ❌ **Ignoring trade-offs:** Fairness improvements may reduce overall accuracy
  - ✅ Solution: Document trade-offs, get clinical approval

---

## Quick Reference: Thresholds

| Metric | Threshold | Action |
|--------|-----------|--------|
| **Representation** | <10% | HIGH RISK - Find diverse alternative |
| **Representation** | <5% | CRITICAL RISK - Do not use |
| **Fairness Disparity** | >10% | Implement mitigation |
| **Fairness Disparity** | >20% | Critical - Do not deploy |
| **Proxy Importance** | >5% | Remove feature, retrain |
| **Confidence Difference** | >15% | Add ancestry-specific disclaimers |
| **Sample Size** | <100 per group | Extend data collection |

---

## Checklist Summary

**Pre-Analysis (Day 1):**
- ☑ Data collected and de-identified
- ☑ Reference datasets documented
- ☑ Ethical review confirmed
- ☑ Tools set up and tested
- ☑ Stakeholders notified

**During Analysis (Days 2-3):**
- ☑ Data representation analyzed
- ☑ Fairness metrics calculated
- ☑ Output stratification completed
- ☑ Audit script executed
- ☑ Manual validation performed

**Post-Analysis (Days 4-5):**
- ☑ Findings documented
- ☑ Risk assessment completed
- ☑ Mitigations implemented
- ☑ Clinical validation obtained
- ☑ Reports archived and communicated

---

**For detailed technical guidance, see:**
- [ETHICS_AND_BIAS.md](ETHICS_AND_BIAS.md) - Comprehensive framework

**Questions?** File an issue on GitHub or consult your institutional bioethics committee.

---

**Document Status:** ✅ Complete (Week 1 Deliverable)
**Last Updated:** 2026-01-12
**Version:** 1.0


---

# PatientOne Bias Audit Demonstration

**Purpose:** Concrete demonstration of bias audit methodology on PatientOne precision medicine workflow

**Patient:** PAT001-OVC-2025 (58-year-old woman, Stage IV HGSOC, platinum-resistant)

**Audit Date:** 2026-01-12

**Auditor:** Precision Medicine MCP Team

**Status:** DEMONSTRATION (Educational Example)

---

## Executive Summary

**Overall Assessment:** ✅ CONDITIONAL PASS

**Risk Level:** MEDIUM

**Key Findings:**
- 3 biases detected (2 MEDIUM risk, 1 LOW risk)
- 0 critical biases
- 5 biases checked and not detected (insurance, geographic, race coding, spatial algorithms, PDX models)

**Recommended Action:** Deploy with implemented mitigations + quarterly monitoring

**Mitigations Required:**
1. Add ancestry diversity warnings to all genomic reports
2. Flag BRCA variants with <5 studies in patient ancestry
3. Document GTEx 85% European composition in expression reports
4. Cross-validate with TOPMed/All of Us for non-European patients
5. Use ovarian-specific cell type references (GSE146026)

---

## Table of Contents

1. [Workflow Overview](#workflow-overview)
2. [Audit Objectives](#audit-objectives)
3. [Genomics Bias Analysis](#genomics-bias-analysis)
4. [Clinical Bias Analysis](#clinical-bias-analysis)
5. [Spatial Transcriptomics Bias Analysis](#spatial-transcriptomics-bias-analysis)
6. [Multiomics Bias Analysis](#multiomics-bias-analysis)
7. [Summary & Recommendations](#summary--recommendations)

---

## Workflow Overview

### Patient Profile

**Demographics:**
- **Age:** 58 years
- **Sex:** Female
- **Ancestry:** Self-reported European (White)
- **Diagnosis:** Stage IV High-Grade Serous Ovarian Carcinoma (HGSOC)
- **Clinical Status:** Platinum-resistant, post-chemotherapy

**Clinical History:**
- Diagnosed 2023 (age 61)
- Treatment: Carboplatin + Paclitaxel + Bevacizumab (2023-2024)
- Progression: Platinum-resistant disease after 6 cycles
- Current: Exploring precision medicine options

---

### Workflow Scope

**Data Modalities Analyzed:**
1. **Genomics:** BRCA1/BRCA2 sequencing, germline variants
2. **Spatial Transcriptomics:** Tumor gene expression, spatial patterns
3. **Multiomics:** RNA-Protein-Phospho integration from PDX models
4. **Clinical:** FHIR data (conditions, medications, observations)

**MCP Servers Used:**
- mcp-epic (clinical EHR data)
- mcp-fgbio (genomic reference data)
- mcp-spatialtools (spatial analysis)
- mcp-multiomics (multi-omics integration)
- mcp-openimagedata (H&E and MxIF histology)
- mcp-tcga (cohort comparison)

**Total Tool Calls:** 27 across 6 servers

**Analysis Duration:** 127 seconds

**Cost:** $1.45 (Claude Sonnet 4.5)

---

## Audit Objectives

### Primary Objectives

1. **Identify potential biases** across all data modalities
2. **Quantify bias severity** (LOW, MEDIUM, HIGH, CRITICAL)
3. **Document findings** with specific examples
4. **Propose mitigations** using diverse reference datasets

### Audit Questions

**Genomics:**
- Are BRCA variant databases Euro-centric?
- Do gene expression references represent diverse populations?
- Are pathway databases universal or population-specific?

**Clinical:**
- Does insurance status influence treatment recommendations?
- Is geographic location used inappropriately?
- Is race/ethnicity used correctly (ancestry for genomics, NOT biology)?

**Spatial:**
- Are cell type references from diverse populations?
- Do spatial algorithms work equally across tissue types?

**Multiomics:**
- Are PDX models representative of patient populations?
- Do drug sensitivity predictions generalize?

---

## Genomics Bias Analysis

### 1. BRCA1/BRCA2 Variant Interpretation

#### Current Data Sources

**Primary References:**
- **ClinVar:** Variant pathogenicity classifications
- **COSMIC:** Cancer-specific mutations
- **gnomAD:** Population allele frequencies

**Ancestry Distribution:**
```
Database  | European | African | Latino | Asian | Other | Source
----------|----------|---------|--------|-------|-------|--------
ClinVar   | ~70%     | ~10%    | ~10%   | ~5%   | ~5%   | Estimated
COSMIC    | ~75%     | ~8%     | ~8%    | ~6%   | ~3%   | Estimated
gnomAD    | 43%      | 21%     | 14%    | 9%    | 8%    | Published
```

---

#### PatientOne Finding

**Variant:** BRCA1 c.5266dupC (p.Gln1756Profs*74)

**ClinVar Classification:** Pathogenic

**Evidence:**
- 50+ studies in European populations
- Founder mutation in Ashkenazi Jewish population (1.2% allele frequency)
- Well-characterized: frameshift leading to premature stop codon
- Strong association with breast/ovarian cancer risk

**gnomAD Frequencies:**
```
Population              | Allele Frequency | Interpretation
------------------------|------------------|----------------
European (non-Finnish)  | 0.0001          | Rare
Ashkenazi Jewish        | 0.012           | Common (founder)
African/African American| 0.00003         | Very rare
Latino/Admixed American | 0.00008         | Very rare
East Asian              | 0.00001         | Extremely rare
South Asian             | 0.00002         | Extremely rare
```

---

#### Bias Assessment

**For PatientOne (European Ancestry):**
- ✅ **APPROPRIATE:** Variant is well-studied in patient's ancestry
- ✅ **CONFIDENCE:** HIGH (50+ European studies)
- ✅ **CLASSIFICATION:** Pathogenic (supported by gnomAD data)

**For Non-European Patients:**
- ⚠️ **LIMITED DATA:** <5 studies in African, Latino, Asian ancestries
- ⚠️ **LOWER CONFIDENCE:** Pathogenicity extrapolated from European data
- ⚠️ **POTENTIAL IMPACT:** May be classified as VUS (variant of uncertain significance) in other ancestries

---

#### Bias Severity

**Risk Level:** 🟡 MEDIUM

**Justification:**
- ClinVar is 70% European → Adequate for European patients
- gnomAD has 43% European, 21% African, 14% Latino → Good diversity
- However, variant-specific studies are Euro-centric (50+ European vs. <5 others)

---

#### Recommended Mitigations

**✅ Implemented in PatientOne Report:**

1. **Cross-reference with gnomAD:**
   - Always check population-specific allele frequencies
   - Flag if frequency >1% in any population (likely benign)

2. **Document study count by ancestry:**
   ```
   Variant: BRCA1 c.5266dupC
   Studies by Ancestry:
     - European: 50+ studies (HIGH confidence)
     - African: 2 studies (LOW confidence)
     - Latino: 1 study (LOW confidence)
     - Asian: 1 study (LOW confidence)
   ```

3. **Add ancestry-specific disclaimer:**
   ```
   ═══════════════════════════════════════════════════════════════
   ⚠️ ANCESTRY CONSIDERATION

   This variant has been extensively studied in European populations
   (50+ studies) but has limited data in non-European ancestries
   (<5 studies each).

   For patients of non-European ancestry, pathogenicity classification
   may have greater uncertainty. Genetic counseling is recommended.

   Recommended Validation:
   - Check gnomAD for population-specific frequencies
   - Reference All of Us (80% underrepresented groups)
   - Consider ClinGen expert panel review
   ═══════════════════════════════════════════════════════════════
   ```

**📊 Recommended Diverse References:**
- **gnomAD v4.1:** Use for population-specific frequencies (21% African, 14% Latino)
- **ClinGen:** Expert-curated variant assessments with ancestry context
- **All of Us:** Cross-validate with 80% underrepresented groups (when available)

---

### 2. Gene Expression Reference Ranges

#### Current Data Source

**Primary Reference:**
- **GTEx (Genotype-Tissue Expression) v8**
  - 948 deceased donors
  - 54 tissue types (including ovary)
  - **Ancestry Distribution:** 85% European, 10% African American, 5% Other, 0% Asian

---

#### PatientOne Finding

**Analysis:** Differential gene expression in ovarian tumor vs. GTEx normal ovary

**Top Upregulated Genes:**
```
Gene   | Patient FC | GTEx Normal | Interpretation
-------|------------|-------------|----------------
MUC16  | +8.2       | Low         | CA-125, tumor marker
TP53   | +6.5       | Medium      | Tumor suppressor (mutated)
BRCA1  | -3.8       | High        | DNA repair (germline mutation)
```

**Reference Range Issue:**
- GTEx normal ovary expression is **85% European**
- For PatientOne (European), this is APPROPRIATE
- For Asian patients, GTEx has **0% representation**

---

#### Bias Assessment

**For PatientOne (European Ancestry):**
- ✅ **APPROPRIATE:** Patient ancestry matches GTEx majority
- ✅ **CONFIDENCE:** HIGH for differential expression calls

**For Non-European Patients:**
- ⚠️ **GTEx BIAS:** 85% European, 0% Asian
- ⚠️ **POTENTIAL IMPACT:** Normal expression ranges may not generalize
- ⚠️ **EXAMPLE:** If baseline MUC16 expression varies by ancestry, fold-change interpretation may be incorrect

---

#### Bias Severity

**Risk Level:** 🟡 MEDIUM

**Justification:**
- GTEx is the gold-standard tissue expression atlas
- 85% European composition is known limitation
- No comparably comprehensive diverse alternative currently available
- Risk mitigated by using multiple references and larger fold-change thresholds

---

#### Recommended Mitigations

**✅ Implemented in PatientOne Report:**

1. **Document GTEx ancestry distribution:**
   ```
   ═══════════════════════════════════════════════════════════════
   ⚠️ REFERENCE DATA LIMITATION

   Gene expression reference ranges are from GTEx v8:
     - European (non-Finnish): 85%
     - African American: 10%
     - Other: 5%
     - Asian: 0%

   Patient Ancestry: European → APPROPRIATE match

   For Asian patients: Results may not accurately reflect normal
   variation. Validation with TOPMed or All of Us recommended.
   ═══════════════════════════════════════════════════════════════
   ```

2. **Apply larger fold-change thresholds for non-European patients:**
   - Standard: log2FC > 1.5 (2.8-fold change)
   - Non-European: log2FC > 2.0 (4-fold change) for conservatism

3. **Cross-validate with diverse references when patient ancestry differs:**
   - **TOPMed:** 180,000+ genomes, diverse US population
   - **Human Cell Atlas:** Single-cell references, global diversity
   - **All of Us:** Validate expression ranges in underrepresented groups

**📊 Recommended Diverse References:**
- **GTEx v8:** Use with explicit 85% European caveat
- **TOPMed RNA-seq:** Validate expression ranges in diverse populations
- **Human Cell Atlas:** Ovarian single-cell data for cell-type-specific expression
- **All of Us Genomic Data:** Cross-validate expression patterns (when available)

---

### 3. Pathway Enrichment Databases

#### Current Data Sources

**Primary References:**
- **KEGG (Kyoto Encyclopedia of Genes and Genomes)**
- **Reactome**
- **GO (Gene Ontology)**

**Ancestry Consideration:**
- Pathway definitions based on **aggregate human data**
- Not stratified by ancestry
- Assumption: Pathways are universal across populations

---

#### PatientOne Finding

**Enriched Pathways (from differential expression):**
```
Pathway                          | p-value  | Genes | Database
---------------------------------|----------|-------|----------
DNA repair (BRCA pathway)        | 1.2e-08  | 12    | KEGG
p53 signaling                    | 3.4e-07  | 15    | Reactome
Apoptosis                        | 8.9e-06  | 18    | GO
Cell cycle checkpoint            | 2.1e-05  | 14    | KEGG
```

---

#### Bias Assessment

**Current Understanding:**
- ✅ Core biological pathways (DNA repair, apoptosis, cell cycle) are **conserved across ancestries**
- ⚠️ Pathway **activity levels** may vary by ancestry due to:
  - Genetic modifiers (ancestry-specific variants affecting pathway genes)
  - Gene expression differences (baseline expression varies by population)
  - Environmental factors (diet, environment interacting with genetics)

**For PatientOne:**
- ✅ **APPROPRIATE:** DNA repair and p53 pathways are universal
- ✅ **CONFIDENCE:** HIGH for pathway enrichment findings

**Potential Concern:**
- ⚠️ **DRUG RESPONSE PATHWAYS:** Pharmacogenomic pathways may show ancestry-specific variation
- ⚠️ **EXAMPLE:** CYP450 metabolism varies significantly by ancestry → Drug dosing implications

---

#### Bias Severity

**Risk Level:** 🟢 LOW

**Justification:**
- Core cancer pathways (DNA repair, p53, apoptosis) are biologically conserved
- Pathway enrichment is statistical test, not absolute truth
- Cross-validation across multiple databases (KEGG + Reactome + GO) mitigates single-source bias

---

#### Recommended Mitigations

**✅ Implemented in PatientOne Report:**

1. **Cross-validate across multiple pathway databases:**
   - Use KEGG **AND** Reactome **AND** GO
   - Flag pathways found in only 1 database (lower confidence)
   - Prioritize pathways found in all 3 (higher confidence)

2. **Flag pathways with known ancestry-specific variation:**
   - Pharmacogenomic pathways (CYP450, drug metabolism)
   - Immune response pathways (HLA variation)
   - Nutrient metabolism (lactose tolerance, alcohol metabolism)

3. **Validate pathway associations in All of Us (when available):**
   - Check if pathway-phenotype associations hold in diverse cohorts
   - Example: "BRCA pathway enrichment → ovarian cancer risk" validated across ancestries?

**📊 Recommended Diverse References:**
- **KEGG + Reactome + GO:** Cross-validate across all three
- **All of Us Workbench:** Validate pathway associations in diverse cohorts (when available)
- **PharmGKB:** Ancestry-specific drug-gene interaction data

---

## Clinical Bias Analysis

### 1. Insurance Status & Treatment Recommendations

#### Audit Objective

**Question:** Does insurance status influence treatment recommendations?

**Ethical Requirement:** Insurance status should NEVER affect clinical recommendations (only billing and access to resources).

---

#### Data Source

**FHIR Resources:**
- `Patient` - Demographics
- `Condition` - Diagnosis (Stage IV HGSOC)
- `Observation` - Lab values (CA-125, BRCA status)
- `MedicationStatement` - Treatment history
- `Coverage` - Insurance information

---

#### PatientOne Finding

**Treatment Recommendation:** Based on:
1. **BRCA1 mutation:** Indicates PARP inhibitor sensitivity
2. **Platinum-resistant:** Rules out platinum re-challenge
3. **CA-125 elevated:** Confirms active disease

**Insurance Status:** NOT used in recommendation algorithm

**Verification:**
```python
# Reviewed mcp-epic tool calls
tool_calls = [
    "get_patient_demographics",      # ✅ Name, age, ancestry
    "get_patient_conditions",         # ✅ Stage IV HGSOC diagnosis
    "get_patient_observations",       # ✅ CA-125, BRCA status
    "get_patient_medications"         # ✅ Treatment history
]

# Coverage resource NOT queried
assert "get_patient_coverage" not in tool_calls  # ✅ PASS
```

**Treatment Recommendation:**
```
Recommended Next-Line Therapy:
1. PARP inhibitor (Olaparib or Niraparib) - BRCA1 mutation indicates sensitivity
2. Clinical trial enrollment - Platinum-resistant HGSOC
3. Consider immunotherapy - PD-L1 status pending

Insurance Status: NOT CONSIDERED in recommendation
```

---

#### Bias Assessment

**Result:** ✅ **NO BIAS DETECTED**

**Evidence:**
- Insurance/coverage data NOT queried during analysis
- Treatment recommendations based solely on:
  - Molecular markers (BRCA1, CA-125)
  - Clinical status (platinum-resistant, Stage IV)
  - Evidence-based guidelines (NCCN ovarian cancer guidelines)

---

#### Bias Severity

**Risk Level:** 🟢 NONE

**Conclusion:** Insurance bias is **NOT present** in PatientOne workflow.

---

### 2. Geographic & Socioeconomic Factors

#### Audit Objective

**Question:** Is geographic location (ZIP code) used as proxy for socioeconomic status to influence treatment eligibility?

**Ethical Requirement:** ZIP code should ONLY be used for:
- ✅ Resource access (finding nearest clinic, pharmacy)
- ✅ Logistics (scheduling, transportation assistance)
- ❌ NOT for treatment eligibility or clinical recommendations

---

#### Data Source

**FHIR Patient.address:**
```json
{
  "address": [{
    "use": "home",
    "line": ["456 Oak Street"],
    "city": "Boston",
    "state": "MA",
    "postalCode": "02108",
    "country": "US"
  }]
}
```

---

#### PatientOne Finding

**Address Usage:**
- ✅ Used for healthcare provider coordination (finding local oncologist)
- ✅ Used for pharmacy selection (nearest specialty pharmacy for PARP inhibitors)
- ❌ **NOT** used in treatment recommendation algorithm

**Verification:**
```python
# Check if ZIP code is in feature importance
feature_importance = {
    "BRCA1_status": 0.45,      # ✅ Clinical
    "CA125_level": 0.32,       # ✅ Clinical
    "platinum_resistant": 0.18, # ✅ Clinical
    "stage": 0.05,             # ✅ Clinical
    "postal_code": 0.0         # ✅ NOT USED
}

assert feature_importance["postal_code"] == 0.0  # ✅ PASS
```

---

#### Bias Assessment

**Result:** ✅ **NO BIAS DETECTED**

**Evidence:**
- ZIP code NOT in algorithm feature set
- No correlation between geographic location and treatment recommendations
- Address used appropriately (resource coordination only)

---

#### Bias Severity

**Risk Level:** 🟢 NONE

**Conclusion:** Geographic/socioeconomic bias is **NOT present** in PatientOne workflow.

---

### 3. Race/Ethnicity Coding

#### Audit Objective

**Question:** Is race/ethnicity used appropriately?

**Ethical Requirements:**
- ✅ **Appropriate:** Ancestry for genomic variant interpretation (genetic context)
- ❌ **Inappropriate:** Race as biological proxy (e.g., assuming kidney function differs by race)

---

#### Data Source

**FHIR Patient.extension (US Core Race/Ethnicity):**
```json
{
  "extension": [{
    "url": "http://hl7.org/fhir/us/core/StructureDefinition/us-core-race",
    "extension": [{
      "url": "ombCategory",
      "valueCoding": {
        "system": "urn:oid:2.16.840.1.113883.6.238",
        "code": "2106-3",
        "display": "White"
      }
    }, {
      "url": "text",
      "valueString": "European ancestry"
    }]
  }]
}
```

---

#### PatientOne Finding

**Ancestry Usage:**
1. **Genomic Variant Interpretation:**
   - ✅ APPROPRIATE: "European ancestry" used to select appropriate gnomAD population frequencies
   - ✅ APPROPRIATE: Checked if BRCA1 variant has sufficient studies in European populations (50+)

2. **Treatment Recommendation:**
   - ✅ APPROPRIATE: Ancestry NOT used to determine treatment eligibility
   - ✅ APPROPRIATE: No race-based treatment algorithms (e.g., no "use Drug X for White patients, Drug Y for Black patients")

3. **Clinical Calculations:**
   - ✅ APPROPRIATE: No race-based eGFR calculation (race coefficient removed)
   - ✅ APPROPRIATE: No race-based dosing adjustments (unless pharmacogenomic indication)

---

#### Bias Assessment

**Result:** ✅ **APPROPRIATE USE**

**Evidence:**
- Ancestry used ONLY for genomic context (matching patient to reference populations)
- Race NOT used as biological proxy
- Consistent with best practices: "Ancestry is genetic, race is social construct"

---

#### Best Practice Validation

**✅ Ancestry for Genomics (CORRECT):**
```
Patient: European ancestry
Variant: BRCA1 c.5266dupC
Reference: gnomAD European population (43% of database)
Interpretation: Pathogenic (consistent with European studies)
```

**❌ Race as Biology (WOULD BE INCORRECT):**
```
# EXAMPLE OF INAPPROPRIATE USE (NOT in PatientOne)
Patient: White race
eGFR Calculation: Multiply creatinine clearance by 1.0 (White) vs. 1.2 (Black)
Conclusion: Race-based medical algorithm (INAPPROPRIATE)
```

**PatientOne does NOT use race-based clinical algorithms.** ✅

---

#### Bias Severity

**Risk Level:** 🟢 NONE

**Conclusion:** Race/ethnicity coding is **APPROPRIATE** in PatientOne workflow.

---

## Spatial Transcriptomics Bias Analysis

### 1. Cell Type Reference Signatures

#### Current Data Source

**Primary Reference:**
- **SingleR default references:** Generic immune cell signatures
  - Source: Bulk RNA-seq from multiple tissue types
  - Ancestry: Not explicitly documented (assumed diverse but unknown)
  - Tissue specificity: Pan-tissue (not ovarian-specific)

---

#### PatientOne Finding

**Cell Type Deconvolution Results:**
```
Cell Type                | Proportion | Confidence
-------------------------|------------|------------
Tumor cells (epithelial) | 68%        | HIGH
T cells (CD8+)           | 12%        | MEDIUM
Fibroblasts (CAFs)       | 15%        | HIGH
Macrophages              | 3%         | MEDIUM
B cells                  | 2%         | LOW
```

**Concern:**
- Generic immune references may NOT capture ovarian cancer-specific cell states
- Example: Tumor-associated macrophages (TAMs) in ovary may differ from TAMs in lung
- Ovarian cancer has unique tumor microenvironment (ascites, peritoneal seeding)

---

#### Bias Assessment

**For PatientOne (European Ancestry):**
- ⚠️ **MODERATE CONCERN:** Using generic (not ovarian-specific) references
- ⚠️ **POTENTIAL IMPACT:** May misclassify ovarian-specific cell states
- ✅ **MITIGATION:** Cross-validate with H&E pathology review

**For Non-European Patients:**
- ⚠️ **UNKNOWN ANCESTRY:** SingleR reference ancestry not documented
- ⚠️ **POTENTIAL IMPACT:** If references are Euro-centric, deconvolution accuracy may vary

---

#### Bias Severity

**Risk Level:** 🟡 LOW-MEDIUM

**Justification:**
- Cell types (epithelial, immune, fibroblast) are biologically conserved
- However, cell **states** (e.g., exhausted T cells, M1 vs. M2 macrophages) may be tissue- and ancestry-specific
- Risk partially mitigated by pathology review

---

#### Recommended Mitigations

**✅ Implemented in PatientOne Workflow:**

1. **Use ovarian-specific reference when available:**
   - **GSE146026:** Single-cell RNA-seq of human ovarian tumors
   - Contains: Tumor cells, T cells, CAFs, macrophages, B cells, dendritic cells
   - Ovarian-specific cell states (e.g., ascites-associated macrophages)

2. **Cross-validate with H&E pathology review:**
   - Pathologist estimates: "70-80% tumor cells, moderate immune infiltration"
   - Deconvolution estimates: "68% tumor cells, 12% T cells"
   - ✅ **CONCORDANCE:** High agreement between deconvolution and pathology

3. **Document reference ancestry (when known):**
   - For Human Cell Atlas references, document donor demographics
   - Flag if patient ancestry is underrepresented in reference

**📊 Recommended Diverse References:**
- **GSE146026:** Ovarian cancer single-cell atlas (tissue-specific)
- **Human Cell Atlas:** Ovarian tissue data (global diversity)
- **Custom references:** Validate with local patient population when possible

---

### 2. Spatial Autocorrelation Algorithms

#### Audit Objective

**Question:** Do spatial algorithms (Moran's I, Geary's C) work equally well across tissue types and ancestries?

---

#### PatientOne Finding

**Spatial Autocorrelation Analysis:**
```
Gene    | Moran's I | p-value | Interpretation
--------|-----------|---------|----------------
MUC16   | 0.72      | <0.001  | Strong spatial clustering
TP53    | 0.58      | <0.001  | Moderate clustering
BRCA1   | 0.35      | 0.003   | Weak clustering
```

**Algorithm:** Moran's I (standard spatial statistics)

---

#### Bias Assessment

**Result:** ✅ **NO BIAS DETECTED**

**Justification:**
- Moran's I is a **mathematical algorithm** (not trained on data)
- Tissue-agnostic: Works on any spatial coordinate system
- Ancestry-agnostic: No population-specific parameters

**Formula (Ancestry-Independent):**
```
Moran's I = (N/W) * Σ(w_ij * (x_i - x̄) * (x_j - x̄)) / Σ(x_i - x̄)²

Where:
  N = number of spatial locations
  W = sum of spatial weights
  w_ij = spatial weight between locations i and j
  x_i = gene expression at location i
  x̄ = mean expression across all locations
```

---

#### Bias Severity

**Risk Level:** 🟢 NONE

**Conclusion:** Spatial autocorrelation algorithms are **NOT biased** by ancestry or tissue type.

---

## Multiomics Bias Analysis

### 1. PDX Model Representativeness

#### Current Data Source

**PDX (Patient-Derived Xenograft) Models:**
- Source: Ovarian cancer patients (diverse cohort)
- Modalities: RNA (gene expression), Protein (TMT proteomics), Phospho (phosphoproteomics)
- Use: Drug sensitivity prediction, upstream regulator analysis

---

#### PatientOne Finding

**Multiomics Integration:**
- Identified: 42 RNA-Protein-Phospho associations
- Top pathway: PI3K/AKT signaling (dysregulated)
- Predicted drug targets: PI3K inhibitors, AKT inhibitors

**PDX Cohort Diversity:**
- Ovarian cancer PDX models from multiple institutions
- Represents diverse patient populations (assumed, but not explicitly validated)

---

#### Bias Assessment

**Concern:**
- ⚠️ **PDX LIMITATION:** Xenografts lack human immune microenvironment
- ⚠️ **ANCESTRY UNKNOWN:** PDX model donor demographics not always documented
- ⚠️ **POTENTIAL IMPACT:** Drug response may differ in patients with intact immune system

**For PatientOne:**
- ✅ **MITIGATION:** Combined PDX predictions with patient's own spatial transcriptomics
  - PDX: Tumor-intrinsic biology (RNA, protein, phospho)
  - Patient spatial: Tumor microenvironment (immune infiltration)

---

#### Bias Severity

**Risk Level:** 🟡 LOW-MEDIUM

**Justification:**
- PDX models are valuable despite limitations
- Combining PDX data with patient-specific spatial data mitigates bias
- Drug predictions should be validated in clinical trials (diverse cohorts)

---

#### Recommended Mitigations

**✅ Implemented in PatientOne Workflow:**

1. **Combine PDX with patient-specific data:**
   - PDX: Tumor-intrinsic drug sensitivity
   - Patient spatial: Microenvironment context (immune infiltration, CAFs)
   - Integration: More comprehensive than PDX alone

2. **Validate drug predictions with clinical trial data:**
   - Check if PI3K inhibitor trials included diverse populations
   - Example: KEYNOTE-100 (pembrolizumab in ovarian cancer) → 23% non-White patients

3. **Document PDX limitations:**
   ```
   ═══════════════════════════════════════════════════════════════
   ⚠️ PDX MODEL LIMITATION

   Drug sensitivity predictions are based on Patient-Derived Xenograft
   (PDX) models, which have the following limitations:

   1. Lack of human immune microenvironment
   2. PDX donor demographics not always documented
   3. Drug response may differ in patients with intact immunity

   Mitigation: Combined PDX predictions with patient's spatial
   transcriptomics data to capture tumor microenvironment context.

   Recommendation: Validate predictions in clinical trials with
   diverse patient populations.
   ═══════════════════════════════════════════════════════════════
   ```

**📊 Recommended Diverse References:**
- **PDX cohorts:** Ensure diverse donor representation when available
- **Clinical trial data:** Cross-validate with diverse trial populations
- **Patient-specific spatial data:** Complement PDX with patient's own microenvironment

---

## Summary & Recommendations

### Biases Detected

#### 1. BRCA Variant Databases: Euro-Centric (MEDIUM Risk)

**Finding:**
- ClinVar/COSMIC are ~70% European
- BRCA1 c.5266dupC has 50+ European studies, <5 in other ancestries

**Impact:**
- ✅ PatientOne (European) → Well-studied, HIGH confidence
- ⚠️ Non-European patients → Lower confidence, may be VUS

**Mitigation:**
- Flag variants with <5 studies in patient ancestry
- Cross-reference gnomAD (43% European, 21% African, 14% Latino)
- Recommend genetic counseling for VUS in underrepresented ancestries

**Status:** ✅ MITIGATED (disclaimers added, gnomAD cross-reference implemented)

---

#### 2. GTEx Reference Ranges: 85% European (MEDIUM Risk)

**Finding:**
- GTEx is 85% European, 10% African, 0% Asian
- Gene expression reference ranges may not generalize

**Impact:**
- ✅ PatientOne (European) → Appropriate match
- ⚠️ Asian patients → 0% representation, HIGH uncertainty

**Mitigation:**
- Document GTEx 85% European composition in all reports
- Apply larger fold-change thresholds for non-European patients (log2FC > 2.0)
- Cross-validate with TOPMed/Human Cell Atlas for non-European patients

**Status:** ✅ MITIGATED (disclaimers added, validation strategy documented)

---

#### 3. Cell Type References: Generic, Not Cancer-Specific (LOW Risk)

**Finding:**
- Using generic immune cell references (not ovarian-specific)
- May miss ovarian cancer-specific cell states

**Impact:**
- ⚠️ Moderate concern for all patients (ancestry-independent)
- Potential misclassification of tissue-specific cell states

**Mitigation:**
- Use ovarian-specific references (GSE146026)
- Cross-validate with H&E pathology review

**Status:** ✅ MITIGATED (ovarian-specific references recommended, pathology concordance confirmed)

---

### Biases NOT Detected

#### 1. Insurance Status ✅

**Finding:** Insurance/coverage data NOT used in treatment recommendations

**Evidence:** Coverage resource NOT queried, feature importance = 0%

**Status:** ✅ NO BIAS

---

#### 2. Geographic Location ✅

**Finding:** ZIP code NOT used for treatment eligibility

**Evidence:** ZIP code feature importance = 0%, used only for resource coordination

**Status:** ✅ NO BIAS

---

#### 3. Race/Ethnicity Coding ✅

**Finding:** Ancestry used appropriately (genomics only, NOT as biological proxy)

**Evidence:** No race-based clinical algorithms (eGFR, dosing, treatment)

**Status:** ✅ APPROPRIATE USE

---

#### 4. Spatial Algorithms ✅

**Finding:** Moran's I is mathematical algorithm, ancestry- and tissue-agnostic

**Evidence:** Algorithm has no population-specific parameters

**Status:** ✅ NO BIAS

---

#### 5. PDX Models ✅ (with mitigation)

**Finding:** PDX models lack immune microenvironment, ancestry unknown

**Mitigation:** Combined with patient-specific spatial data

**Status:** ✅ MITIGATED

---

### Overall Risk Assessment

**Genomics:** 🟡 MEDIUM
- BRCA databases Euro-centric
- Mitigated by gnomAD cross-reference and disclaimers

**Gene Expression:** 🟡 MEDIUM
- GTEx 85% European, 0% Asian
- Mitigated by ancestry documentation and validation strategy

**Clinical:** 🟢 LOW
- No insurance or geographic bias
- Appropriate use of ancestry

**Spatial:** 🟢 LOW
- Generic cell references (tissue, not ancestry issue)
- Mitigated by ovarian-specific references

**Multiomics:** 🟢 LOW-MEDIUM
- PDX limitations acknowledged
- Mitigated by combining with patient spatial data

**OVERALL:** 🟡 MEDIUM (acceptable with implemented mitigations)

---

### Deployment Decision

**✅ RECOMMENDED: Deploy with Mitigations + Quarterly Monitoring**

**Conditions:**
1. All mitigations implemented (ancestry disclaimers, diverse references)
2. Quarterly bias audits to monitor for drift
3. Genetic counseling available for VUS cases
4. Continuous monitoring of fairness metrics

---

### Next Steps

**Immediate (Before Production):**
- [x] Add ancestry diversity warnings to all genomic reports
- [x] Document GTEx 85% European composition
- [x] Implement gnomAD cross-reference
- [x] Use ovarian-specific cell type references

**Quarterly (Every 3 Months):**
- [ ] Re-run fairness metrics on new patient data
- [ ] Check for updated diverse reference datasets
- [ ] Validate no bias drift over time

**Annual (Every 12 Months):**
- [ ] Comprehensive bias audit (all 4 steps)
- [ ] Update reference datasets (gnomAD, GTEx, All of Us)
- [ ] Review clinical validation with expert panel

---

**For methodology details, see:**
- [ETHICS_AND_BIAS.md](ETHICS_AND_BIAS.md) - Comprehensive framework

**Questions?** File an issue on GitHub or consult institutional bioethics committee.

---

**Document Status:** ✅ Complete (Week 1 Deliverable)
**Last Updated:** 2026-01-12
**Version:** 1.0
