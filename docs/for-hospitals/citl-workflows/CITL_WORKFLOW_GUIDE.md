# Clinician-in-the-Loop (CitL) Workflow Guide

**Target Audience:** Oncologists, Pathologists, Radiologists
**Purpose:** Step-by-step instructions for reviewing and validating AI-generated precision medicine reports
**Time Required:** 20-30 minutes per report
**Version:** 1.1
**Last Updated:** 2026-02-19

---

## Table of Contents

1. [Overview](#overview)
2. [Workflow Diagram](#workflow-diagram)
3. [Step-by-Step Instructions](#step-by-step-instructions)
4. [Section-by-Section Guidance](#section-by-section-guidance)
5. [Decision Scenarios](#decision-scenarios)
6. [Troubleshooting](#troubleshooting)
7. [Best Practices](#best-practices)
8. [FAQs](#faqs)
9. [Examples](#examples)
10. [Review Template](#review-template)

---

## Overview

### What is CitL Review?

The **Clinician-in-the-Loop (CitL) validation workflow** is a formal review and sign-off process where you, as the treating oncologist or specialist, validate AI-generated precision medicine analysis results before they are incorporated into clinical decision-making.

**Think of it as:** The AI system is your research assistant who has analyzed all the patient's multi-omics data (genomics, spatial transcriptomics, imaging) and prepared a draft report. Your job is to review the draft for accuracy, clinical validity, and guideline compliance before approving it for use.

### Why CitL Matters

**AI as Co-Pilot, Not Autopilot:**
- High-stakes oncology decisions require human expertise and judgment
- AI can miss clinical context, guideline nuances, or patient-specific factors
- Regulatory and ethical standards require human accountability for clinical decisions
- CitL creates an audit trail demonstrating appropriate oversight

**What You'll Validate:**
1. **Molecular findings accuracy** - Are the reported mutations, gene expression changes, and pathway alterations clinically valid?
2. **Clinical guideline compliance** - Do treatment recommendations follow NCCN and institutional protocols?
3. **Data quality** - Are sample sizes adequate? Are statistical thresholds met?
4. **Treatment appropriateness** - Are recommended therapies suitable for this specific patient?

### Your Role & Responsibility

**As the reviewing clinician, you attest that:**
- You have carefully reviewed all findings in the draft report
- You have assessed alignment with NCCN and institutional clinical guidelines
- You have evaluated automated quality flags and determined whether they are acceptable
- Your approval decision reflects your independent clinical judgment
- You understand this review becomes part of the permanent medical record (10-year retention per HIPAA)

**Decision Authority:**
- **APPROVE** - Report is accurate and ready for clinical use (proceed to tumor board, treatment planning)
- **REVISE** - Report requires modifications (specify issues, bioinformatician will re-analyze)
- **REJECT** - Report has critical errors requiring complete re-analysis (escalate to PI/bioinformatics team)

---

## Workflow Diagram

```
┌─────────────────────────────────────────────────────────────┐
│  AI Analysis Pipeline (Automated)                           │
│  ├─ TEST_1: Clinical data extraction (EHR/FHIR)            │
│  ├─ TEST_2: Multi-omics integration                        │
│  ├─ TEST_3: Spatial transcriptomics analysis               │
│  ├─ TEST_4: Imaging analysis                               │
│  └─ TEST_5: Cross-modal integration & synthesis            │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│  STEP 1: Generate Draft Report (Automated, ~30 sec)        │
│  ├─ Run quality checks (sample size, FDR, completeness)    │
│  ├─ Generate draft_report.json                             │
│  └─ Generate clinical_summary.txt                          │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│  STEP 2: Clinician Review (Manual, 20-30 min) ⭐ YOU       │
│  ├─ Open draft report files                                │
│  ├─ Complete CitL Review Template (7 sections)             │
│  ├─ Validate molecular findings (10 top findings)          │
│  ├─ Check NCCN guideline compliance                        │
│  ├─ Assess quality flags                                   │
│  └─ Make decision: APPROVE / REVISE / REJECT               │
└─────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────┐
│  STEP 3: Submit Review (Automated, ~5 sec)                 │
│  ├─ Validate review against JSON schema                    │
│  ├─ Generate digital signature (SHA-256)                   │
│  └─ Log to Cloud Logging audit trail                       │
└─────────────────────────────────────────────────────────────┘
                            ↓
                    ┌───────┴───────┐
                    │   Decision?   │
                    └───────┬───────┘
                            ↓
        ┌───────────────────┼───────────────────┐
        ↓                   ↓                   ↓
┌───────────────┐  ┌────────────────┐  ┌──────────────┐
│  APPROVE      │  │  REVISE        │  │  REJECT      │
├───────────────┤  ├────────────────┤  ├──────────────┤
│ Finalize      │  │ Re-analyze     │  │ Escalate to  │
│ report        │  │ with adjusted  │  │ PI/team      │
│               │  │ parameters     │  │              │
│ Ready for     │  │                │  │ Schedule     │
│ tumor board   │  │ Resubmit for   │  │ review       │
│               │  │ review         │  │ meeting      │
└───────────────┘  └────────────────┘  └──────────────┘
```

---

## Step-by-Step Instructions

### STEP 1: Access Draft Report (Automated - Already Complete)

**What Happened Before You Arrived:**
The bioinformatics team or automated pipeline has:
1. Run AI analysis pipeline (TEST_1 through TEST_5)
2. Generated draft report with quality checks
3. Notified you that a report is ready for review

**Files You'll Receive:**
```
results/PAT001-OVC-2025/
├── draft_report.json          # Structured findings (JSON format)
├── quality_checks.json         # Automated QC results
├── clinical_summary.txt        # Human-readable summary
└── citl_review_template.md     # Form to complete (copy provided)
```

**Your Action:** Open these files in your preferred editor or viewer.

---

### STEP 2: Review Clinical Summary (5-10 minutes)

**File to Open:** `clinical_summary.txt`

**What to Look For:**

1. **Patient Demographics & Context**
   - Verify patient ID matches your expectations
   - Review cancer type, stage, prior treatments
   - Note any contraindications or special considerations

2. **Quality Check Status**
   ```
   PASS: Sample sizes adequate (>=30 spots per region)
   PASS: FDR thresholds met (FDR < 0.05 for significant findings)
   PASS: Data completeness (>95% complete)
   WARNING: One region has minimum size (35 spots)
   ```

   **Interpretation:**
   - **All PASS:** Proceed with review, data quality is good
   - **Warnings present:** Note the warning, assess whether it affects clinical validity (you'll address this in Section 4)
   - **Critical failures:** Contact bioinformatics team before proceeding

3. **Top Molecular Findings**
   - Review list of top 10 genes/mutations/pathways identified
   - Note confidence scores (FDR values)
   - Identify any findings that seem inconsistent with clinical presentation

4. **Treatment Recommendations**
   - Review suggested therapies
   - Note molecular targets and expected efficacy
   - Consider patient-specific factors (comorbidities, prior treatments, performance status)

**Red Flags to Watch For:**
- Findings contradict known clinical presentation (e.g., ER+ markers in triple-negative breast cancer)
- Treatment recommendations contraindicated by patient history
- Suspicious FDR values (too many findings with FDR exactly = 0.05)
- Sample sizes below 30 spots in multiple regions

**If Red Flags Present:** Make note in your review form and consider REVISE or REJECT decision.

---

### STEP 3: Complete Review Form (15-20 minutes)

**File to Complete:** `citl_review_template.md` (copy to `citl_review_completed.json`)

**Overview of 7 Sections:**
1. **High-Level Decision** (REQUIRED) - APPROVE / REVISE / REJECT
2. **Per-Finding Validation** (REQUIRED) - Validate top 10 molecular findings
3. **Clinical Guideline Compliance** (REQUIRED) - NCCN + institutional protocols
4. **Quality Flags Assessment** - Assess automated QC flags
5. **Treatment Recommendations Review** - Evaluate treatment suggestions
6. **Attestation & Signature** (REQUIRED) - Digital signature
7. **Revision Instructions** - Complete if REVISE or REJECT

---

## Section-by-Section Guidance

### SECTION 1: High-Level Decision (REQUIRED)

**What to Do:**
1. Select ONE decision: `[x] APPROVE` or `[x] REVISE` or `[x] REJECT`
2. Write 2-3 sentence rationale explaining your decision

**Decision Criteria:**

**APPROVE when:**
- All findings are consistent with clinical presentation
- No critical quality flags
- Treatment recommendations follow NCCN guidelines
- Data quality is adequate for clinical decision-making
- You are confident presenting these results at tumor board

**REVISE when:**
- Findings are mostly valid but require minor corrections
- Quality flags need addressing (e.g., exclude small regions)
- Treatment recommendations need refinement
- Additional tests or analyses would strengthen conclusions
- You can provide specific instructions for re-analysis

**REJECT when:**
- Critical errors in findings (e.g., wrong patient data, major inconsistencies)
- Data quality is inadequate (multiple failed QC checks)
- Results are clinically implausible
- Complete re-analysis required (not just parameter adjustments)

**Example Rationales:**

**Good APPROVE rationale:**
> "All findings are consistent with clinical presentation and imaging. Molecular results align with observed platinum resistance. Treatment recommendations follow NCCN guidelines for recurrent HGSOC. Quality checks passed with one minor warning that does not affect clinical validity."

**Good REVISE rationale:**
> "Findings are valid but tumor_necrotic region (35 spots) should be excluded due to insufficient sample size. FDR values for DEGs in this region may be unreliable. Request re-analysis excluding this region and increasing minimum spots per region to 50."

**Good REJECT rationale:**
> "Critical inconsistency: Report identifies BRCA1 germline mutation, but patient's prior genetic testing (2024-03-15) showed BRCA wild-type. Suspect sample mix-up or data error. Recommend verifying sample identity before re-analysis."

**Poor rationales:**
- "Looks good" (too vague, no clinical reasoning)
- "I disagree with the AI" (lacks specifics)
- "Not sure about this" (uncertain - use REVISE instead)

---

### SECTION 2: Per-Finding Validation (REQUIRED)

**What to Do:**
For each of the 10 top molecular findings listed in `draft_report.json`, indicate validation status:
- `[x] CONFIRMED` - Finding is clinically valid and supported by data
- `[x] UNCERTAIN` - Requires additional validation or expert consultation
- `[x] INCORRECT` - Finding appears to be erroneous

**How to Validate Each Finding:**

**Finding Structure:**
```
Finding 1: TP53 (R175H) - Somatic Mutation
- Finding ID: MUT_1
- Confidence Score: FDR = 5.04e-20
- Clinical Significance: Tumor suppressor loss, TP53 pathway disrupted
- Validation Status: [CONFIRM / UNCERTAIN / INCORRECT]
- Comments: [Your assessment]
```

**Validation Questions:**

1. **Does this finding align with clinical presentation?**
   - Example: TP53 mutations are expected in 96% of HGSOC cases -> CONFIRMED
   - Example: HER2 amplification in histologically confirmed TNBC -> INCORRECT

2. **Is the confidence score appropriate?**
   - FDR < 0.01 -> High confidence
   - FDR 0.01-0.05 -> Moderate confidence (acceptable)
   - FDR > 0.05 -> Low confidence -> Mark UNCERTAIN, request re-analysis

3. **Does the clinical significance make sense?**
   - Check if gene function matches reported effect
   - Example: "BCL2L1 overexpression -> apoptosis resistance" -> Biologically plausible -> CONFIRMED

4. **Can you verify with other data?**
   - Cross-reference with imaging (e.g., high Ki67 expression should correlate with high proliferation on imaging)
   - Cross-reference with prior genomic reports if available
   - Check protein expression (IHC) if available

**Examples:**

CONFIRMED Example:
```
Finding 1: TP53 (R175H) - Somatic Mutation
- Validation Status: [x] CONFIRMED
- Comments: "TP53 R175H is a hotspot mutation in HGSOC (present in >95% of cases).
  Consistent with p53 staining pattern on IHC (diffuse nuclear positivity).
  Very high confidence (FDR = 5.04e-20). CONFIRMED."
```

UNCERTAIN Example:
```
Finding 5: NRAS (Q61K) - Somatic Mutation
- Validation Status: [x] UNCERTAIN
- Comments: "NRAS mutations uncommon in HGSOC (<5% prevalence). Confidence score
  is moderate (FDR = 0.03). Recommend orthogonal validation with targeted
  sequencing (Sanger or ddPCR) before incorporating into treatment decisions."
```

INCORRECT Example:
```
Finding 8: EGFR (Amplification) - Copy Number Variant
- Validation Status: [x] INCORRECT
- Comments: "EGFR amplification is rare in HGSOC (<2% prevalence) and not supported
  by IHC staining (EGFR 0/3+). Suspected false positive. Recommend excluding
  this finding from final report."
```

**Pro Tip:** If you have 2+ UNCERTAIN or 1+ INCORRECT findings, consider REVISE decision in Section 1.

---

### SECTION 3: Clinical Guideline Compliance (REQUIRED)

**What to Do:**
Assess whether treatment recommendations align with:
1. **NCCN Guidelines** (National Comprehensive Cancer Network)
2. **Institutional Protocols** (your hospital's formulary and treatment pathways)

**For Each Guideline Set:**

Select alignment level:
- `[x] ALIGNED` - All recommendations follow guidelines
- `[x] PARTIAL` - Some recommendations deviate (explain below)
- `[x] NOT_ALIGNED` - Significant deviations from guidelines (explain below)

**How to Assess:**

1. **Open NCCN Guidelines for Patient's Cancer Type**
   - Example: NCCN Guidelines for Ovarian Cancer (current version)
   - Navigate to relevant section (e.g., "Recurrent or Persistent Disease")

2. **Compare Recommended Therapies**
   - Check if suggested drugs are in NCCN-recommended regimens
   - Verify sequencing (e.g., platinum-based before PARP inhibitors)
   - Check Category ratings (Category 1 preferred, Category 2A acceptable, Category 2B conditional)

3. **Document Deviations**
   - If recommendation is Category 2B or off-label, explain justification
   - If recommendation deviates from institutional formulary, note alternative options

**Examples:**

ALIGNED Example:
```
NCCN Guidelines Alignment: [x] ALIGNED
NCCN Deviations: None

Institutional Protocol Alignment: [x] ALIGNED
Institutional Deviations: None

Comment: All recommended therapies (carboplatin/paclitaxel, olaparib, bevacizumab)
are NCCN Category 1 or 2A recommendations for platinum-sensitive recurrent HGSOC.
All drugs are on institutional formulary. No deviations.
```

PARTIAL Example:
```
NCCN Guidelines Alignment: [x] PARTIAL
NCCN Deviations:
- "Alpelisib (PI3K inhibitor) is Category 2B (limited data) for PIK3CA-mutant
  ovarian cancer. NCCN prefers clinical trial enrollment. However, patient has
  PIK3CA E545K mutation with strong rationale for PI3K targeting. Recommend
  discussing with patient: clinical trial (preferred) vs. off-label alpelisib."

Institutional Protocol Alignment: [x] ALIGNED
Institutional Deviations: None
```

NOT_ALIGNED Example (should trigger REVISE):
```
NCCN Guidelines Alignment: [x] NOT_ALIGNED
NCCN Deviations:
- "Recommendation for imatinib (PDGFR inhibitor) is not supported by NCCN
  guidelines for HGSOC. Imatinib has no established efficacy in this disease.
  Recommend removing this recommendation."

Institutional Protocol Alignment: [x] NOT_ALIGNED
Institutional Deviations:
- "Bevacizumab continuation beyond 15 months is not per institutional protocol
  (protocol limits to 12-15 months due to GI perforation risk). Recommend
  discontinuing bevacizumab at 15 months per protocol."
```

**Pro Tip:** PARTIAL alignment is common and acceptable if deviations are justified. NOT_ALIGNED typically requires REVISE decision.

---

### SECTION 4: Quality Flags Assessment

**What to Do:**
Review automated quality flags raised during analysis. For each flag, assess:
- `[x] ACCEPTABLE` - Flag noted but not concerning, does not affect clinical validity
- `[x] REQUIRES_ACTION` - Flag must be addressed before approval

**Common Quality Flags:**

**1. Sample Size Warning**
```
Flag: sample_size_warning
Severity: warning
Message: "Minimum region size: 35 spots (tumor_necrotic region)"
Recommendation: "Prefer >=50 spots for robust statistical power. Consider excluding
                small regions or interpreting results with caution."

Your Assessment: [ACCEPTABLE / REQUIRES_ACTION]
```

**When to mark ACCEPTABLE:**
- Region is 35-49 spots AND FDR values are very low (< 1e-10) -> Strong signal compensates for smaller sample
- Region is biologically unimportant (e.g., necrotic tissue) -> Can be excluded without affecting conclusions

**When to mark REQUIRES_ACTION:**
- Region is <35 spots -> Too small for reliable statistics
- FDR values are marginal (0.01-0.05) -> Insufficient power -> Request re-analysis excluding small regions

**2. Missing Data Warning**
```
Flag: missing_data_warning
Severity: warning
Message: "5% of spots have missing gene expression values (technical dropout)"
Recommendation: "Acceptable if <10%. Values were imputed using k-NN method."

Your Assessment: [ACCEPTABLE / REQUIRES_ACTION]
```

**When ACCEPTABLE:** <10% missing, imputation method is appropriate
**When REQUIRES_ACTION:** >10% missing, indicates poor data quality

**3. Batch Effect Warning**
```
Flag: batch_effect_detected
Severity: critical
Message: "Significant batch effect detected between sample batches (p < 0.001)"
Recommendation: "Apply batch correction (ComBat or Harmony) before analysis."

Your Assessment: [ACCEPTABLE / REQUIRES_ACTION]
```

**When ACCEPTABLE:** Batch correction was already applied, residual effect is minimal
**When REQUIRES_ACTION:** No batch correction applied -> REVISE and request correction

**Example Assessment:**
```
Flag 1: sample_size_warning (tumor_necrotic region, 35 spots)
Reviewer Assessment: [x] ACCEPTABLE
Comments: "Minimum region size of 35 spots is acceptable. While below ideal threshold
of 50, the strong effect sizes (log2FC > 3) and very low FDR values (< 1e-10) provide
sufficient confidence. Additionally, tumor_necrotic region findings are corroborated
by imaging (necrotic areas on CT). ACCEPTABLE."
```

**Pro Tip:** If ANY flag is marked REQUIRES_ACTION, you must select REVISE in Section 1 and provide re-analysis instructions in Section 7.

---

### SECTION 5: Treatment Recommendations Review

**What to Do:**
For each treatment recommendation, indicate:
- `[x] AGREE` - Recommendation is appropriate for this patient
- `[x] DISAGREE` - Recommendation is not appropriate (explain below)

**Evaluation Criteria:**

1. **Molecular Target Match**
   - Does the patient have the molecular alteration targeted by this therapy?
   - Example: PIK3CA inhibitor requires PIK3CA mutation -> Check if present

2. **Expected Efficacy**
   - Is the reported response rate realistic based on literature?
   - Is the patient likely to benefit given disease stage and prior treatments?

3. **Patient-Specific Factors**
   - Contraindications (e.g., bevacizumab contraindicated if recent GI perforation)
   - Performance status (ECOG 0-2 required for most trials)
   - Organ function (adequate renal/hepatic function)
   - Prior treatment history (e.g., platinum-free interval)

4. **Availability**
   - Is this therapy available at your institution?
   - Is it covered by insurance?
   - Are there clinical trials available?

**Examples:**

AGREE Example:
```
Recommendation 1: PI3K inhibitor (Alpelisib) + PARP inhibitor (Olaparib)
Molecular Target: PIK3CA E545K mutation + BRCA1 germline mutation
Expected Efficacy: 40-50% ORR in PIK3CA-mutant, PARP-sensitive tumors

Agreement: [x] AGREE
Comments: "Patient has both PIK3CA E545K (confirmed by NGS) and BRCA1 germline
mutation (confirmed by prior germline testing). Platinum-free interval is 8 months
(platinum-sensitive). Combination is rational based on synthetic lethality. Patient
meets eligibility criteria for clinical trial NCT12345678. AGREE - recommend
discussing clinical trial enrollment."
```

DISAGREE Example:
```
Recommendation 3: Bevacizumab continuation (cycle 20+)
Molecular Target: VEGFA overexpression
Expected Efficacy: Modest PFS benefit (2-3 months)

Agreement: [x] DISAGREE
Comments: "Patient has already received 19 cycles of bevacizumab (>18 months).
Institutional protocol limits bevacizumab to 12-15 months due to cumulative GI
perforation risk. Patient has history of diverticulosis (relative contraindication).
Risk exceeds benefit at this point. DISAGREE - recommend discontinuing bevacizumab
and considering alternative anti-angiogenic approaches if indicated."
```

**Pro Tip:** It's okay to DISAGREE with some recommendations - this demonstrates appropriate clinical oversight. Document your reasoning clearly.

---

### SECTION 6: Attestation & Signature (REQUIRED)

**What to Do:**
1. Check ALL attestation boxes (must be checked to submit):
   - `[x]` I have reviewed all findings in this report
   - `[x]` I have assessed clinical guideline compliance
   - `[x]` I have evaluated quality flags
   - `[x]` My decision reflects my clinical judgment
   - `[x]` I understand this review is part of the medical record

2. Fill in reviewer information:
   - Full Name: `Dr. Sarah Johnson`
   - Credentials: `MD, Gynecologic Oncology`
   - Email: `sarah.johnson@hospital.org`
   - Role: `oncologist`
   - Date: `2026-02-19` (YYYY-MM-DD)
   - Time: `14:30 PST` (HH:MM timezone)

3. Leave Digital Signature blank (auto-generated by submission script)

**Legal & Regulatory Notes:**

This attestation is legally binding:
- Your review becomes part of the permanent medical record
- 10-year retention per HIPAA requirements
- Audit trail includes your name, email, timestamp, and decision
- Digital signature (SHA-256 hash) ensures immutability

You are attesting that:
- You personally reviewed the report (not delegated to resident/fellow)
- Your decision reflects your independent clinical judgment (not influenced by administrative pressure)
- You are qualified to evaluate this cancer type and molecular findings

---

### SECTION 7: Revision Instructions (Required if REVISE or REJECT)

**When to Complete:** ONLY if you selected REVISE or REJECT in Section 1

**What to Do:**

1. **List Specific Issues to Address** (minimum 1 issue)
   ```
   Issue 1: Exclude tumor_necrotic region (35 spots) due to insufficient sample size.
            FDR values for DEGs in this region are unreliable.

   Issue 2: Re-run differential expression analysis with stricter FDR threshold
            (FDR < 0.01 instead of FDR < 0.05) to reduce false positives.

   Issue 3: Validate TP53 mutation with orthogonal sequencing method (Sanger or ddPCR)
            to confirm hotspot mutation.
   ```

2. **Specify Re-analysis Parameters**
   ```
   Regions to exclude:
   - tumor_necrotic (n=35 spots)

   FDR threshold adjustment:
   - Change from FDR < 0.05 to FDR < 0.01 for differential expression

   Minimum spots per region:
   - Increase from 30 to 50 spots minimum per region

   Additional tests required:
   - Sanger sequencing confirmation for TP53 R175H
   ```

3. **Set Resubmission Timeline**
   ```
   Expected resubmission date: 2026-02-26
   Estimated time for corrections: 3-5 days
   ```

**Pro Tip:** Be specific! Vague instructions like "improve data quality" are not actionable. Provide exact parameters to adjust.

---

## Decision Scenarios

### Scenario 1: Clean Approval (Most Common)

**Situation:**
- All quality checks passed
- Findings consistent with clinical presentation
- NCCN-aligned treatment recommendations
- No critical quality flags

**Your Actions:**
1. Section 1: Select **APPROVE**
2. Section 2: Mark 8-10 findings as **CONFIRMED**, 0-2 as **UNCERTAIN**
3. Section 3: Mark **ALIGNED** for both NCCN and institutional
4. Section 4: Mark all flags as **ACCEPTABLE** (if any)
5. Section 5: **AGREE** with most/all treatment recommendations
6. Section 6: Complete attestation
7. Section 7: Leave blank (not required for APPROVE)

**Estimated Time:** 20-25 minutes

**Next Steps:** Submit review -> Finalize report -> Present at tumor board

---

### Scenario 2: Revise Due to Quality Flag

**Situation:**
- Findings are valid but one region has inadequate sample size
- Quality flag raised: sample_size_warning
- You determine sample size is too small for reliable statistics

**Your Actions:**
1. Section 1: Select **REVISE**
   - Rationale: "Quality flag for tumor_necrotic region (35 spots) requires action. Sample size below threshold for reliable differential expression. Request re-analysis excluding this region."

2. Section 2: Mark findings from adequate regions as **CONFIRMED**, findings from small region as **UNCERTAIN**

3. Section 3: Mark **ALIGNED** (no guideline issues)

4. Section 4: Mark sample_size_warning as **REQUIRES_ACTION**
   - Comment: "35 spots insufficient for robust statistics. FDR values marginal (0.02-0.04). Exclude tumor_necrotic region and re-run analysis."

5. Section 5: **AGREE** with treatments (assuming they're not based on unreliable findings)

6. Section 6: Complete attestation

7. Section 7: **MUST COMPLETE**
   ```
   Issues to Address:
   1. Exclude tumor_necrotic region (35 spots) from differential expression analysis

   Re-analysis Parameters:
   - Regions to exclude: tumor_necrotic
   - Minimum spots per region: Increase to 50 spots

   Resubmission Timeline: 2026-02-22 (3 days)
   ```

**Estimated Time:** 25-35 minutes (additional time for Section 7)

**Next Steps:** Submit review -> Bioinformatician re-runs analysis -> You review version 2

---

### Scenario 3: Reject Due to Critical Error

**Situation:**
- Findings contradict known clinical data
- Example: Report shows BRCA1 mutation but prior genetic testing was wild-type
- Suspect sample mix-up or data error

**Your Actions:**
1. Section 1: Select **REJECT**
   - Rationale: "Critical inconsistency: Report identifies BRCA1 germline mutation, but patient's prior genetic testing (2024-03-15, GeneDx report #12345) showed BRCA wild-type. Suspect sample mix-up or data error. Cannot proceed with clinical decision-making until discrepancy is resolved. Recommend verifying sample identity (fingerprinting) and re-running from raw data."

2. Section 2: Mark contradictory findings as **INCORRECT**, others as **UNCERTAIN**

3. Section 3: Cannot assess (report validity in question)

4. Section 4: Mark critical flags as **REQUIRES_ACTION**

5. Section 5: Cannot assess (recommendations based on invalid data)

6. Section 6: Complete attestation

7. Section 7: **MUST COMPLETE**
   ```
   Issues to Address:
   1. Verify sample identity - confirm correct patient specimen used
   2. Reconcile BRCA1 discrepancy (report shows BRCA1 mutation, prior test wild-type)
   3. If sample mix-up confirmed, re-run entire analysis pipeline with correct sample
   4. If data processing error, re-run from raw sequencing data

   Re-analysis Parameters:
   - Verify sample barcode matches patient ID
   - Perform sample fingerprinting (SNP concordance check)
   - If sample correct: Investigate BRCA1 call (tumor vs. germline, technical artifact)

   Resubmission Timeline: 2026-03-05 (2 weeks for investigation + re-analysis)
   ```

**Estimated Time:** 30-40 minutes (requires investigation and documentation)

**Next Steps:** Submit review -> Escalate to PI/bioinformatics -> Schedule team meeting -> Corrective action -> Complete re-analysis

---

## Troubleshooting

### Problem: Draft report files missing or incomplete

**Symptoms:**
- `draft_report.json` not found
- Files are empty or corrupted

**Solutions:**
1. Contact bioinformatics team - analysis may still be running
2. Check file paths - ensure you're looking in correct patient directory
3. Re-generate draft report:
   ```bash
   python servers/mcp-patient-report/scripts/generate_patient_report.py --patient-id PAT001-OVC-2025 --output-dir ./results --generate-draft
   ```

---

### Problem: I don't understand a molecular finding

**Example:** "What does 'PIK3CA E545K' mean?"

**Solutions:**
1. **Look up gene function:**
   - GeneCards: https://www.genecards.org/
   - OMIM: https://www.omim.org/
   - Example: PIK3CA = phosphoinositide 3-kinase, E545K = hotspot mutation in helical domain

2. **Check clinical significance:**
   - ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/
   - OncoKB: https://www.oncokb.org/
   - Example: PIK3CA E545K = Oncogenic, targetable with PI3K inhibitors

3. **Mark as UNCERTAIN and request expert consultation:**
   ```
   Finding: PIK3CA E545K
   Validation Status: [x] UNCERTAIN
   Comments: "Requesting molecular pathology consultation to confirm clinical
   significance and treatment implications before finalizing approval."
   ```

---

### Problem: Treatment recommendation seems off-label or experimental

**Example:** Report suggests drug not in standard NCCN guidelines

**Solutions:**
1. **Check if supported by molecular rationale:**
   - Is there a strong biological justification?
   - Is there clinical trial data or case reports supporting use?

2. **Consult NCCN guidelines carefully:**
   - Drug may be Category 2B (acceptable with lower evidence)
   - May be recommended for specific molecular subtype

3. **Document in Section 3 (Guidelines) and Section 5 (Treatments):**
   ```
   NCCN Alignment: [x] PARTIAL
   Deviation: "Alpelisib is Category 2B for PIK3CA-mutant ovarian cancer.
   Prefer clinical trial enrollment per NCCN. However, patient has exhausted
   standard options and has strong molecular rationale (PIK3CA E545K).
   Acceptable as off-label use after discussing trial availability."
   ```

4. **Mark AGREE with caveat:**
   ```
   Recommendation: Alpelisib (PI3K inhibitor)
   Agreement: [x] AGREE
   Comments: "Off-label use acceptable given PIK3CA E545K mutation and exhausted
   standard therapies. Recommend prior authorization and informed consent
   discussion with patient."
   ```

---

### Problem: Quality flag is unclear or I'm unsure if it's acceptable

**Example:** "Batch effect detected - is this a problem?"

**Solutions:**
1. **Read the flag details carefully:**
   - What is the severity? (critical / warning / info)
   - What does the recommendation say?
   - Has correction been applied?

2. **Consult with bioinformatics team:**
   - Ask: "Was batch correction applied?"
   - Ask: "Does this affect findings validity?"

3. **If unsure, mark REQUIRES_ACTION and request clarification:**
   ```
   Flag: batch_effect_detected
   Reviewer Assessment: [x] REQUIRES_ACTION
   Comments: "Unclear whether batch correction was applied. Request confirmation
   from bioinformatics team. If not applied, request re-analysis with ComBat or
   Harmony batch correction before approval."
   ```

4. **Select REVISE in Section 1** and document in Section 7

---

### Problem: Schema validation errors when submitting review

**Symptoms:**
```
Schema validation failed:
   'status' is a required property
   Path: decision -> status
```

**Solutions:**
1. **Check required fields:**
   - Section 1: Decision status and rationale
   - Section 2: At least one finding validated
   - Section 3: NCCN and institutional alignment selected
   - Section 6: All attestation checkboxes checked

2. **Check data types:**
   - Dates must be YYYY-MM-DD format
   - Email must be valid email format
   - Enums must exactly match (e.g., "APPROVE" not "Approve")

3. **Check conditional requirements:**
   - If status = REVISE or REJECT -> Section 7 is required
   - If quality flags exist -> Section 4 is required

4. **Validate JSON syntax:**
   - Use online JSON validator: https://jsonlint.com/
   - Check for missing commas, brackets, quotes

---

## Best Practices

### Time Management

**Efficient Review (20-25 minutes):**
1. **Skim clinical summary (5 min)** - Get high-level overview
2. **Section 1 decision (2 min)** - Quick initial assessment (APPROVE/REVISE/REJECT)
3. **Section 2 findings (8 min)** - Validate top 10 findings (~1 min each)
4. **Section 3 guidelines (3 min)** - Check NCCN alignment
5. **Section 4-5 (5 min)** - Quick flag and treatment review
6. **Section 6 attestation (2 min)** - Complete signature
7. **Buffer (5 min)** - Double-check, finalize

**When to Take More Time:**
- Complex molecular findings you're unfamiliar with (consult references)
- Multiple quality flags requiring assessment
- Off-label or experimental treatment recommendations
- REVISE or REJECT decisions (Section 7 requires detailed instructions)

---

### Communication with Bioinformatics Team

**When to Reach Out:**
1. **Before starting review:** If files are missing or analysis incomplete
2. **During review:** If you encounter unclear findings or flags
3. **After REVISE decision:** Confirm they understand re-analysis instructions
4. **After REJECT decision:** Schedule meeting to discuss corrective action

**How to Communicate:**
- Be specific: "Exclude tumor_necrotic region" not "improve data quality"
- Provide rationale: Explain why change is needed
- Set timeline: Expected resubmission date
- Use patient ID consistently: Avoid confusion with multiple patients

---

### Documentation Standards

**Good Comments:**
- Specific, clinical reasoning
- References to supporting data (IHC, imaging, prior reports)
- Actionable instructions (for REVISE)

**Example:**
> "TP53 R175H confirmed. Hotspot mutation present in 96% of HGSOC cases (TCGA data). Consistent with p53 IHC staining pattern (diffuse nuclear positivity, consistent with missense mutation). Very high confidence (FDR = 5.04e-20). CONFIRMED."

**Poor Comments:**
- Vague, no reasoning
- No clinical context
- Non-actionable

**Example:**
> "Looks good"

---

### Quality Assurance

**Self-Check Before Submitting:**
- [ ] All REQUIRED sections completed (1, 2, 3, 6)
- [ ] Decision (APPROVE/REVISE/REJECT) matches findings
- [ ] If REVISE/REJECT: Section 7 completed with specific instructions
- [ ] If APPROVE: You're confident presenting at tumor board
- [ ] Attestation checkboxes all checked
- [ ] Reviewer information complete and accurate
- [ ] No typos in patient ID or critical fields

---

## FAQs

**Q: How long does a typical review take?**
A: 20-30 minutes for straightforward APPROVE, 25-35 minutes for REVISE, 30-40 minutes for REJECT.

**Q: Can I delegate review to a fellow or resident?**
A: No. The attestation requires that YOU personally reviewed the report. Fellows/residents can assist with data review, but final validation and signature must be yours.

**Q: What if I'm not confident in molecular findings?**
A: Mark findings as UNCERTAIN and request expert consultation (molecular pathologist, geneticist). You can still APPROVE with UNCERTAIN findings if they don't affect treatment decisions, or select REVISE and request additional validation.

**Q: What happens after I submit an APPROVE review?**
A: The finalize_patient_report.py script generates a final approved report with status "clinically_approved". This report is ready for tumor board presentation and clinical decision-making.

**Q: What happens after I submit a REVISE review?**
A: Bioinformatics team receives your revision instructions (Section 7), re-runs analysis with adjusted parameters, and generates a new draft report (version 2). You'll then review version 2 using the same CitL workflow.

**Q: What happens after I submit a REJECT review?**
A: Your REJECT decision triggers escalation to PI/bioinformatics team lead. A team meeting is scheduled to investigate the critical error, determine corrective action, and create a plan for re-analysis.

**Q: Can I change my decision after submitting?**
A: No. Once submitted, the review is immutable (digital signature ensures tamper-evidence). If you realize you made an error, contact bioinformatics team and PI to discuss next steps. A new review version may be required.

**Q: What if I disagree with NCCN guidelines?**
A: Document your reasoning in Section 3. If you have strong clinical justification for deviating from guidelines (e.g., patient-specific factors, recent literature), you can mark PARTIAL and explain the deviation. This is acceptable and demonstrates appropriate clinical judgment.

**Q: What if quality checks failed but findings still seem valid?**
A: Assess each flag individually in Section 4. Some flags are informational and don't affect validity (mark ACCEPTABLE). Others require action (mark REQUIRES_ACTION and select REVISE). Use your clinical judgment.

**Q: Can I approve a report with UNCERTAIN findings?**
A: Yes, IF the uncertain findings don't materially affect treatment decisions. Document your reasoning: "Finding X marked UNCERTAIN pending additional validation, but primary treatment recommendation (Y) is supported by other CONFIRMED findings." If uncertain findings are critical to treatment decisions, select REVISE and request additional validation.

**Q: How do I convert the markdown template to JSON?**
A: Currently manual (copy values from .md to .json following schema structure). A web-based form interface is planned for future release to automate this conversion.

**Q: Where can I find example completed reviews?**
A: See the [Examples](#examples) section below for three complete example reviews (APPROVE, REVISE, REJECT scenarios).

**Q: Who can I contact for help?**
A:
- Technical questions (files, scripts, schema): Bioinformatics team
- Clinical questions (molecular findings): Molecular pathology or genetics
- Workflow questions (process, timeline): Clinic administrator or PI

---

## Summary Checklist

Before submitting your review, confirm:

**Required Sections Complete:**
- [ ] Section 1: Decision (APPROVE/REVISE/REJECT) and rationale
- [ ] Section 2: All 10 findings validated (CONFIRMED/UNCERTAIN/INCORRECT)
- [ ] Section 3: NCCN and institutional guideline compliance assessed
- [ ] Section 6: All attestation boxes checked, reviewer info complete

**Conditional Sections:**
- [ ] Section 4: Quality flags assessed (if flags exist)
- [ ] Section 5: Treatment recommendations reviewed (if recommendations exist)
- [ ] Section 7: Revision instructions (if REVISE or REJECT)

**Quality Checks:**
- [ ] Decision is consistent with findings (e.g., APPROVE but all findings INCORRECT -> inconsistent)
- [ ] Comments are specific and clinically justified
- [ ] Patient ID is correct
- [ ] You're confident in your decision

**Ready to Submit:**
```bash
python servers/mcp-patient-report/scripts/citl_submit_review.py \
  --patient-id PAT001-OVC-2025 \
  --review-file ./results/PAT001-OVC-2025/citl_review_completed.json
```

---

## Examples

Complete example reviews demonstrating APPROVE, REVISE, and REJECT scenarios for clinicians learning the CitL workflow.

### Example 1: APPROVE - PatientOne Clean Approval

**Scenario:** High-grade serous ovarian carcinoma (HGSOC), Stage IV, platinum-resistant recurrence. All quality checks passed, findings consistent with clinical presentation, NCCN-aligned treatment recommendations.

**Review Time:** 25 minutes
**Expected Outcome:** Final approved report ready for tumor board

#### Complete Review Form

```json
{
  "patient_id": "PAT001-OVC-2025",
  "report_date": "2026-01-13T14:00:00Z",
  "reviewer": {
    "name": "Dr. Sarah Johnson",
    "email": "sarah.johnson@hospital.org",
    "credentials": "MD, Gynecologic Oncology",
    "role": "oncologist"
  },
  "review_date": "2026-01-13T14:30:00Z",

  "decision": {
    "status": "APPROVE",
    "rationale": "All findings are consistent with clinical presentation and imaging. Molecular results align with observed platinum resistance (ABCB1/MDR1 overexpression, BCL2L1 anti-apoptotic signaling). Treatment recommendations follow NCCN guidelines for recurrent HGSOC. Quality checks passed with one minor warning (tumor_necrotic region 35 spots) that does not affect clinical validity given strong effect sizes and corroborating imaging."
  },

  "per_finding_validation": [
    {
      "finding_id": "DEG_1",
      "gene": "TP53",
      "validation_status": "CONFIRMED",
      "comments": "TP53 R175H is a well-characterized hotspot mutation present in >95% of HGSOC cases. Consistent with p53 IHC staining pattern (diffuse nuclear positivity). Very high confidence (FDR = 5.04e-20). CONFIRMED."
    },
    {
      "finding_id": "DEG_2",
      "gene": "PIK3CA",
      "validation_status": "CONFIRMED",
      "comments": "PIK3CA E545K activation mutation (helical domain). Found in ~10-15% of HGSOC cases. Consistent with PI3K/AKT pathway upregulation observed in pathway analysis. High confidence (FDR = 2.31e-18). Actionable target for PI3K inhibitor therapy. CONFIRMED."
    },
    {
      "finding_id": "DEG_3",
      "gene": "BRCA1",
      "validation_status": "CONFIRMED",
      "comments": "BRCA1 germline mutation confirmed by prior genetic testing (GeneDx report 2024-03-15, variant c.5266dupC, p.Gln1756Profs). Consistent with homologous recombination deficiency signature. Rationale for PARP inhibitor therapy. CONFIRMED."
    },
    {
      "finding_id": "DEG_4",
      "gene": "ABCB1",
      "validation_status": "CONFIRMED",
      "comments": "ABCB1 (MDR1) overexpression consistent with platinum resistance. Patient has platinum-free interval of 4 months (platinum-resistant per GOG criteria). Log2FC = 3.8, FDR = 1.2e-15. Strong biological rationale. CONFIRMED."
    },
    {
      "finding_id": "DEG_5",
      "gene": "BCL2L1",
      "validation_status": "CONFIRMED",
      "comments": "BCL2L1 overexpression indicates anti-apoptotic signaling. Mechanism of therapy resistance. Log2FC = 3.2, FDR = 8.7e-14. Consistent with apoptosis resistance phenotype. Potential target for BCL2 inhibitors. CONFIRMED."
    },
    {
      "finding_id": "DEG_6",
      "gene": "CD8A",
      "validation_status": "CONFIRMED",
      "comments": "CD8A downregulation in tumor core regions consistent with immune exclusion phenotype observed on imaging (lack of T-cell infiltration on IHC). Log2FC = -2.9, FDR = 4.3e-12. Explains poor response to prior checkpoint inhibitor. CONFIRMED."
    },
    {
      "finding_id": "DEG_7",
      "gene": "VEGFA",
      "validation_status": "CONFIRMED",
      "comments": "VEGFA overexpression consistent with angiogenic phenotype. Patient responded to bevacizumab previously. Log2FC = 2.7, FDR = 1.8e-11. Rationale for continued anti-angiogenic therapy. CONFIRMED."
    },
    {
      "finding_id": "DEG_8",
      "gene": "MKI67",
      "validation_status": "CONFIRMED",
      "comments": "MKI67 overexpression indicates high proliferation rate. Consistent with Ki67 IHC (75% positivity). Log2FC = 3.5, FDR = 9.2e-13. Aggressive phenotype. CONFIRMED."
    },
    {
      "finding_id": "DEG_9",
      "gene": "AKT1",
      "validation_status": "CONFIRMED",
      "comments": "AKT1 upregulation consistent with PI3K/AKT pathway activation (PIK3CA E545K mutation). Log2FC = 2.3, FDR = 7.4e-10. Mechanistically linked to PIK3CA mutation. CONFIRMED."
    },
    {
      "finding_id": "DEG_10",
      "gene": "MTOR",
      "validation_status": "CONFIRMED",
      "comments": "MTOR upregulation downstream of PI3K/AKT activation. Log2FC = 1.9, FDR = 2.1e-9. Consistent with pathway analysis showing mTOR signaling activation. Potential combination target (PI3K + mTOR inhibitor). CONFIRMED."
    }
  ],

  "guideline_compliance": {
    "nccn_aligned": "ALIGNED",
    "nccn_deviations": [],
    "institutional_aligned": "ALIGNED",
    "institutional_deviations": []
  },

  "quality_flags_assessment": [
    {
      "flag_id": "sample_size_warning",
      "severity": "warning",
      "reviewer_assessment": "ACCEPTABLE",
      "comments": "Tumor_necrotic region has 35 spots (below ideal 50 threshold). However, strong effect sizes (log2FC > 3) and very low FDR values (< 1e-10) provide sufficient confidence. Additionally, necrotic regions corroborated by CT imaging (hypodense areas consistent with necrosis). Findings from this region are not critical to treatment decisions. ACCEPTABLE."
    }
  ],

  "treatment_recommendations_review": [
    {
      "recommendation_id": "REC_1",
      "therapy_name": "PI3K inhibitor (Alpelisib) + PARP inhibitor (Olaparib)",
      "agreement": "AGREE",
      "comments": "Combination is rational given PIK3CA E545K mutation + BRCA1 germline deficiency. Patient meets eligibility for clinical trial NCT04729387 (alpelisib + olaparib in PIK3CA-mutant, BRCA-deficient ovarian cancer). Platinum-free interval is 4 months (platinum-resistant), making PARP inhibitor monotherapy less effective. Combination addresses both PI3K activation and HR deficiency. AGREE - recommend discussing trial enrollment."
    },
    {
      "recommendation_id": "REC_2",
      "therapy_name": "Bevacizumab continuation",
      "agreement": "AGREE",
      "comments": "Patient responded to bevacizumab in prior line (PFS 9 months, above median). VEGFA overexpression provides biological rationale. No contraindications (no GI perforation history, no recent surgery). Duration is 8 cycles (within institutional protocol limit of 15 cycles). AGREE."
    },
    {
      "recommendation_id": "REC_3",
      "therapy_name": "Liposomal doxorubicin (alternative option)",
      "agreement": "AGREE",
      "comments": "NCCN Category 1 recommendation for platinum-resistant ovarian cancer. Appropriate as alternative if patient declines trial or is ineligible. Cardiac function adequate (LVEF 60%). AGREE as backup option."
    }
  ],

  "attestation": {
    "reviewed_all_findings": true,
    "assessed_compliance": true,
    "clinical_judgment": true,
    "medical_record_acknowledgment": true,
    "signature_hash": "a1b2c3d4e5f6789012345678901234567890123456789012345678901234abcd",
    "timestamp": "2026-01-13T14:55:00Z"
  },

  "revision_count": 0
}
```

#### Key Decision Points

**Why APPROVE?**
1. All 10 findings CONFIRMED with strong clinical rationale
2. Quality checks passed (one minor warning acceptable)
3. NCCN-aligned treatment recommendations
4. Findings consistent with clinical presentation (platinum resistance, imaging, IHC)
5. Confident presenting at tumor board

**Next Steps:**
1. Run finalization script:
   ```bash
   python servers/mcp-patient-report/scripts/finalize_patient_report.py --patient-id PAT001-OVC-2025
   ```
2. Present at tumor board (Tuesday)
3. Discuss clinical trial enrollment with patient
4. Document decision in EHR

---

### Example 2: REVISE - Quality Flag Requires Action

**Scenario:** Breast cancer, ER+/HER2-, Stage III. Quality flag raised for inadequate sample size in one region. Findings from that region have marginal FDR values. Clinician determines re-analysis is needed with region excluded.

**Review Time:** 30 minutes
**Expected Outcome:** Re-analysis with adjusted parameters, resubmission for review

#### Complete Review Form

```json
{
  "patient_id": "PAT002-BRC-2025",
  "report_date": "2026-01-12T10:00:00Z",
  "reviewer": {
    "name": "Dr. Michael Chen",
    "email": "michael.chen@hospital.org",
    "credentials": "MD, PhD, Medical Oncology",
    "role": "oncologist"
  },
  "review_date": "2026-01-12T15:45:00Z",

  "decision": {
    "status": "REVISE",
    "rationale": "Quality flag for tumor_edge region (28 spots) requires action. Sample size is below minimum threshold of 30 spots, and FDR values from this region are marginal (0.02-0.04), indicating insufficient statistical power. Findings from adequate regions (tumor_core: 420 spots, stroma: 380 spots) are valid. Request re-analysis excluding tumor_edge region and increasing minimum spots threshold to 50 to prevent similar issues."
  },

  "per_finding_validation": [
    {
      "finding_id": "DEG_1",
      "gene": "ESR1",
      "validation_status": "CONFIRMED",
      "comments": "ESR1 (estrogen receptor) high expression consistent with ER+ status (IHC: 95% positive). From tumor_core region (420 spots). FDR = 1.2e-25. CONFIRMED."
    },
    {
      "finding_id": "DEG_2",
      "gene": "PGR",
      "validation_status": "CONFIRMED",
      "comments": "PGR (progesterone receptor) moderate expression consistent with PR+ status (IHC: 60% positive). From tumor_core region. FDR = 3.4e-18. CONFIRMED."
    },
    {
      "finding_id": "DEG_3",
      "gene": "ERBB2",
      "validation_status": "CONFIRMED",
      "comments": "ERBB2 (HER2) low expression consistent with HER2- status (IHC 0). From tumor_core region. FDR = 5.1e-20. CONFIRMED."
    },
    {
      "finding_id": "DEG_4",
      "gene": "MKI67",
      "validation_status": "UNCERTAIN",
      "comments": "MKI67 upregulation from tumor_edge region (28 spots). FDR = 0.024 (marginal). Ki67 IHC shows heterogeneity (15-30% across tumor). Insufficient statistical power in small region. Mark UNCERTAIN pending re-analysis with larger sample."
    },
    {
      "finding_id": "DEG_5",
      "gene": "CCND1",
      "validation_status": "UNCERTAIN",
      "comments": "CCND1 overexpression from tumor_edge region (28 spots). FDR = 0.037 (marginal). Cyclin D1 IHC inconclusive. Insufficient power. Mark UNCERTAIN pending re-analysis."
    },
    {
      "finding_id": "DEG_6",
      "gene": "CDK4",
      "validation_status": "UNCERTAIN",
      "comments": "CDK4 upregulation from tumor_edge region (28 spots). FDR = 0.041 (marginal). Related to CCND1 finding. Both require validation with larger sample. Mark UNCERTAIN."
    },
    {
      "finding_id": "DEG_7",
      "gene": "FOXC1",
      "validation_status": "CONFIRMED",
      "comments": "FOXC1 downregulation in tumor_core (luminal phenotype). FDR = 2.8e-12. Consistent with ER+ luminal subtype. CONFIRMED."
    },
    {
      "finding_id": "DEG_8",
      "gene": "GATA3",
      "validation_status": "CONFIRMED",
      "comments": "GATA3 high expression in tumor_core. Luminal lineage marker. FDR = 1.5e-16. CONFIRMED."
    },
    {
      "finding_id": "DEG_9",
      "gene": "BCL2",
      "validation_status": "CONFIRMED",
      "comments": "BCL2 overexpression in tumor_core. Associated with ER+ tumors. FDR = 4.2e-14. CONFIRMED."
    },
    {
      "finding_id": "DEG_10",
      "gene": "PTEN",
      "validation_status": "CONFIRMED",
      "comments": "PTEN expression retained in tumor_core. No loss of function. FDR = 8.9e-11. PI3K pathway not constitutively activated. CONFIRMED."
    }
  ],

  "guideline_compliance": {
    "nccn_aligned": "ALIGNED",
    "nccn_deviations": [],
    "institutional_aligned": "ALIGNED",
    "institutional_deviations": []
  },

  "quality_flags_assessment": [
    {
      "flag_id": "sample_size_warning",
      "severity": "warning",
      "reviewer_assessment": "REQUIRES_ACTION",
      "comments": "Tumor_edge region has only 28 spots (below minimum 30). Three findings (MKI67, CCND1, CDK4) from this region have marginal FDR values (0.024-0.041), indicating insufficient statistical power. Cannot confidently use these findings for treatment decisions. REQUIRES_ACTION: Exclude tumor_edge region and re-run differential expression analysis."
    }
  ],

  "treatment_recommendations_review": [
    {
      "recommendation_id": "REC_1",
      "therapy_name": "Endocrine therapy (Aromatase inhibitor + CDK4/6 inhibitor)",
      "agreement": "AGREE",
      "comments": "Standard NCCN Category 1 recommendation for ER+/HER2- advanced breast cancer. ESR1, PGR expression confirmed. However, CDK4 overexpression (rationale for CDK4/6 inhibitor) comes from unreliable tumor_edge region. Still AGREE with recommendation as CDK4/6 inhibitors are standard regardless of CDK4 expression. Re-analysis will strengthen molecular rationale."
    },
    {
      "recommendation_id": "REC_2",
      "therapy_name": "Chemotherapy (alternative if endocrine-resistant)",
      "agreement": "AGREE",
      "comments": "Appropriate alternative if endocrine therapy fails. MKI67 status from tumor_edge region is uncertain, but clinical features (node-positive, large tumor) suggest intermediate-high risk. AGREE."
    }
  ],

  "attestation": {
    "reviewed_all_findings": true,
    "assessed_compliance": true,
    "clinical_judgment": true,
    "medical_record_acknowledgment": true,
    "signature_hash": "b2c3d4e5f678901234567890123456789012345678901234567890123456bcde",
    "timestamp": "2026-01-12T16:15:00Z"
  },

  "revision_instructions": {
    "issues_to_address": [
      "Exclude tumor_edge region (28 spots) due to insufficient sample size. FDR values from this region (0.024-0.041) are unreliable.",
      "Re-run differential expression analysis with minimum spot threshold increased to 50 to prevent similar issues in future analyses.",
      "Three findings (MKI67, CCND1, CDK4) from tumor_edge region should be re-evaluated after excluding this region. These may still be detected in adjacent regions if truly present."
    ],
    "reanalysis_parameters": {
      "regions_to_exclude": ["tumor_edge"],
      "min_spots_per_region": 50,
      "fdr_threshold": 0.05,
      "additional_tests_required": []
    },
    "resubmission_date": "2026-01-15"
  },

  "revision_count": 0
}
```

#### Key Decision Points

**Why REVISE (not APPROVE)?**
1. Quality flag requires action (sample size < 30)
2. Three findings UNCERTAIN due to marginal FDR values
3. Insufficient statistical power in one region
4. Other findings valid, no NCCN deviations
5. Cannot confidently use uncertain findings for treatment decisions

**Why REVISE (not REJECT)?**
- Majority of findings (7/10) are CONFIRMED from adequate regions
- Issue is localized to one small region
- Clear path to resolution (exclude problematic region)
- Treatment recommendations still appropriate
- No critical errors or data integrity issues

**Next Steps:**
1. Bioinformatics team receives revision instructions
2. Re-run analysis excluding tumor_edge region, min spots = 50
3. Generate new draft report (version 2)
4. Clinician reviews version 2 (expected time: 15-20 minutes, faster than initial)
5. Likely outcome: APPROVE on version 2

---

### Example 3: REJECT - Critical Data Error

**Scenario:** Lung adenocarcinoma, Stage IV. Report shows EGFR exon 19 deletion, but patient's prior molecular testing (6 months ago, FoundationOne) showed EGFR wild-type. Critical inconsistency suggests sample mix-up or data error.

**Review Time:** 35 minutes
**Expected Outcome:** Escalation to PI/bioinformatics team, investigation, complete re-analysis

#### Complete Review Form

```json
{
  "patient_id": "PAT003-NSCLC-2025",
  "report_date": "2026-01-11T09:00:00Z",
  "reviewer": {
    "name": "Dr. Jennifer Martinez",
    "email": "jennifer.martinez@hospital.org",
    "credentials": "MD, Thoracic Oncology",
    "role": "oncologist"
  },
  "review_date": "2026-01-11T14:20:00Z",

  "decision": {
    "status": "REJECT",
    "rationale": "Critical inconsistency: Report identifies EGFR exon 19 deletion (p.E746_A750del), but patient's prior molecular testing (FoundationOne, 2025-07-15, report #FMI-12345) showed EGFR wild-type (no mutations, amplifications, or deletions). Patient has not received EGFR TKI therapy that could select for resistant clones. Suspect sample mix-up (wrong patient specimen) or critical data processing error. Cannot proceed with clinical decision-making until discrepancy is resolved. Recommend verifying sample identity (SNP fingerprinting) and investigating EGFR call. If sample identity confirmed, re-run from raw sequencing data."
  },

  "per_finding_validation": [
    {
      "finding_id": "MUT_1",
      "gene": "EGFR",
      "validation_status": "INCORRECT",
      "comments": "EGFR exon 19 deletion (p.E746_A750del) contradicts prior FoundationOne testing (2025-07-15) which showed EGFR wild-type. Patient has not received EGFR TKI therapy. No biological explanation for acquiring EGFR mutation in 6 months without TKI selection pressure. Suspect sample mix-up or data error. Mark INCORRECT pending investigation."
    },
    {
      "finding_id": "MUT_2",
      "gene": "TP53",
      "validation_status": "UNCERTAIN",
      "comments": "TP53 R273H reported. Prior FoundationOne testing showed TP53 R248Q (different hotspot). Both are gain-of-function mutations, but different amino acid positions. If sample identity is correct, this suggests tumor evolution or sampling from different clone. If sample mix-up, this is from different patient. Mark UNCERTAIN pending sample verification."
    },
    {
      "finding_id": "MUT_3",
      "gene": "KRAS",
      "validation_status": "UNCERTAIN",
      "comments": "KRAS wild-type reported. Consistent with prior testing (KRAS WT). However, EGFR and KRAS mutations are typically mutually exclusive. If EGFR deletion is real, KRAS WT is expected. But EGFR deletion contradicts prior data. Overall validity UNCERTAIN pending investigation."
    },
    {
      "finding_id": "DEG_4",
      "gene": "ALK",
      "validation_status": "UNCERTAIN",
      "comments": "ALK not rearranged (consistent with prior testing). Cannot validate other findings until sample identity confirmed. Mark UNCERTAIN."
    },
    {
      "finding_id": "DEG_5",
      "gene": "ROS1",
      "validation_status": "UNCERTAIN",
      "comments": "ROS1 not rearranged (consistent with prior testing). Cannot validate until sample verified. Mark UNCERTAIN."
    },
    {
      "finding_id": "DEG_6",
      "gene": "BRAF",
      "validation_status": "UNCERTAIN",
      "comments": "BRAF wild-type. Consistent with prior testing. Mark UNCERTAIN pending resolution."
    },
    {
      "finding_id": "DEG_7",
      "gene": "MET",
      "validation_status": "UNCERTAIN",
      "comments": "MET amplification not detected. Consistent with prior testing. Mark UNCERTAIN."
    },
    {
      "finding_id": "DEG_8",
      "gene": "RET",
      "validation_status": "UNCERTAIN",
      "comments": "RET not rearranged. Consistent with prior. Mark UNCERTAIN."
    },
    {
      "finding_id": "DEG_9",
      "gene": "ERBB2",
      "validation_status": "UNCERTAIN",
      "comments": "ERBB2 (HER2) no mutations detected. Consistent with prior. Mark UNCERTAIN."
    },
    {
      "finding_id": "DEG_10",
      "gene": "NTRK1",
      "validation_status": "UNCERTAIN",
      "comments": "NTRK not rearranged. Consistent with prior. Mark UNCERTAIN."
    }
  ],

  "guideline_compliance": {
    "nccn_aligned": "NOT_ALIGNED",
    "nccn_deviations": [
      "Treatment recommendation for EGFR TKI (osimertinib) is based on EGFR exon 19 deletion that contradicts prior molecular testing. Cannot assess guideline compliance until EGFR discrepancy is resolved."
    ],
    "institutional_aligned": "NOT_ALIGNED",
    "institutional_deviations": [
      "EGFR TKI therapy contradicts institutional protocol requiring confirmed EGFR mutation. Prior testing showed wild-type."
    ]
  },

  "quality_flags_assessment": [
    {
      "flag_id": "mutation_concordance_error",
      "severity": "critical",
      "reviewer_assessment": "REQUIRES_ACTION",
      "comments": "CRITICAL: EGFR mutation call discordant with prior molecular testing. This is not a quality flag from automated checks, but a critical finding discovered during clinical review. REQUIRES_ACTION: Immediate investigation required."
    }
  ],

  "treatment_recommendations_review": [
    {
      "recommendation_id": "REC_1",
      "therapy_name": "EGFR TKI (Osimertinib)",
      "agreement": "DISAGREE",
      "comments": "Recommendation is based on EGFR exon 19 deletion that contradicts prior testing. Cannot administer EGFR TKI based on unverified molecular finding. DISAGREE - patient should NOT receive osimertinib until EGFR status is definitively confirmed."
    },
    {
      "recommendation_id": "REC_2",
      "therapy_name": "Platinum-based chemotherapy (alternative)",
      "agreement": "AGREE",
      "comments": "Standard NCCN recommendation for EGFR wild-type NSCLC (per prior FoundationOne testing). This remains appropriate treatment if EGFR mutation is ruled out. AGREE as appropriate alternative."
    }
  ],

  "attestation": {
    "reviewed_all_findings": true,
    "assessed_compliance": true,
    "clinical_judgment": true,
    "medical_record_acknowledgment": true,
    "signature_hash": "c3d4e5f6789012345678901234567890123456789012345678901234567cdef",
    "timestamp": "2026-01-11T14:55:00Z"
  },

  "revision_instructions": {
    "issues_to_address": [
      "CRITICAL: Verify sample identity. Confirm correct patient specimen was analyzed. Perform SNP fingerprinting to match sample to patient germline.",
      "Investigate EGFR exon 19 deletion call. Reconcile discrepancy with prior FoundationOne testing (2025-07-15, EGFR wild-type).",
      "If sample identity confirmed correct: Investigate technical reasons for EGFR discrepancy (sequencing error, variant caller settings, reference genome mismatch). Consider orthogonal validation (Sanger sequencing, ddPCR) for EGFR locus.",
      "If sample mix-up confirmed: Identify correct patient specimen, re-run complete analysis pipeline with correct sample.",
      "Review TP53 discrepancy (R273H vs R248Q from prior testing). If sample correct, this may indicate tumor evolution or multi-clonal disease.",
      "Do NOT proceed with EGFR TKI therapy until EGFR status definitively confirmed. Patient harm could result from inappropriate targeted therapy."
    ],
    "reanalysis_parameters": {
      "regions_to_exclude": [],
      "additional_tests_required": [
        "SNP fingerprinting for sample identity verification",
        "Orthogonal EGFR sequencing (Sanger or ddPCR) for exon 19 region",
        "Re-run from raw FASTQ files if data processing error suspected",
        "Review variant caller logs for EGFR locus",
        "Compare with prior FoundationOne VCF file if available"
      ]
    },
    "resubmission_date": "2026-01-25"
  },

  "revision_count": 0
}
```

#### Key Decision Points

**Why REJECT (not REVISE)?**
1. **Critical error:** EGFR mutation contradicts prior definitive testing
2. **Patient safety:** Recommending inappropriate targeted therapy (EGFR TKI)
3. **Data integrity concern:** Suspect sample mix-up (most serious error)
4. **Cannot use ANY findings:** Validity of entire report is questionable
5. **Requires investigation:** Not just parameter adjustment, need root cause analysis

**Next Steps:**
1. **Immediate:** Alert bioinformatics team and PI (critical error)
2. **Day 1:** Schedule emergency meeting (clinician, bioinformatics, lab director)
3. **Day 2-3:** Sample verification (SNP fingerprinting)
4. **Day 4-7:** Investigation and root cause analysis
   - If sample mix-up -> Identify correct sample, re-run analysis
   - If data error -> Debug pipeline, validate EGFR call, re-run from raw data
5. **Day 10-14:** Complete re-analysis and resubmit for review
6. **Patient care:** Continue with standard platinum chemotherapy (per prior FoundationOne results) while investigation ongoing. Do NOT start EGFR TKI.

---

### Examples Comparison Table

| Aspect | APPROVE | REVISE | REJECT |
|--------|---------|---------|---------|
| **Decision Status** | All findings valid, ready for clinical use | Findings mostly valid, needs corrections | Critical error, cannot use report |
| **Confidence Level** | High - comfortable presenting at tumor board | Medium - need adjustments before clinical use | Low - data validity in question |
| **Findings Validation** | 8-10 CONFIRMED, 0-2 UNCERTAIN, 0 INCORRECT | 5-8 CONFIRMED, 2-5 UNCERTAIN, 0-2 INCORRECT | 0-3 CONFIRMED, 5-10 UNCERTAIN, 1-3 INCORRECT |
| **Quality Flags** | All ACCEPTABLE | 1-2 REQUIRES_ACTION | Critical flags or data integrity concern |
| **Guideline Compliance** | ALIGNED or PARTIAL (justified) | ALIGNED (but molecular rationale needs strengthening) | NOT_ALIGNED (based on invalid data) |
| **Treatment Recs** | AGREE with all/most | AGREE with most (pending molecular confirmation) | DISAGREE with key recommendations |
| **Next Step** | Finalize report -> Tumor board | Re-analysis -> Resubmit for review | Investigation -> Root cause analysis -> Complete re-analysis |
| **Timeline** | Immediate (same day) | 3-5 days (re-analysis) | 1-2 weeks (investigation + re-analysis) |
| **Section 7 Required?** | No | Yes - specific re-analysis instructions | Yes - detailed investigation plan |
| **Review Time** | 20-25 minutes | 25-35 minutes | 30-40 minutes |
| **Typical Frequency** | 70-80% of reviews | 15-25% of reviews | 5-10% of reviews |

### Decision-Making Framework

**Ask yourself:**
1. **Would I present these findings at tumor board?**
   - Yes -> APPROVE
   - Not yet, but close -> REVISE
   - No, something is wrong -> REJECT

2. **Would I start treatment based on these findings?**
   - Yes -> APPROVE
   - Yes, but with additional validation -> REVISE
   - No, data may be incorrect -> REJECT

3. **What's the primary issue?**
   - No significant issues -> APPROVE
   - Statistical power, minor quality flags -> REVISE
   - Data validity, patient safety -> REJECT

4. **How confident am I in the findings?**
   - High (>90%) -> APPROVE
   - Medium (60-90%) -> REVISE
   - Low (<60%) -> REJECT

### Examples Best Practices

1. **Always review prior molecular testing reports** - Catch discrepancies early
2. **Don't ignore quality flags** - Even "minor" warnings can indicate real issues
3. **Be specific in revision instructions** - "Improve quality" is not actionable
4. **Document your reasoning** - Future reviewers and auditors will read your comments
5. **Err on the side of caution** - When in doubt, REVISE rather than APPROVE
6. **Patient safety first** - REJECT if any concern about data validity or patient harm

---

## Review Template

This form captures your clinical review and validation of AI-generated precision medicine analysis results. Complete all required sections below. This review becomes part of the permanent medical record and audit trail (10-year retention per HIPAA).

**Patient ID:** {{PATIENT_ID}}
**Report Date:** {{REPORT_DATE}}
**Reviewer:** {{REVIEWER_NAME}}
**Review Date:** {{REVIEW_DATE}}

**Time estimate:** 20-30 minutes
**Required:** All sections except comments fields
**Format:** Convert to JSON using provided template or automation tool

---

### TEMPLATE SECTION 1: HIGH-LEVEL DECISION (REQUIRED)

**Overall Assessment:**
- [ ] **APPROVE** - Report is accurate and ready for clinical use
- [ ] **REVISE** - Report requires modifications (complete Section 7)
- [ ] **REJECT** - Report has critical errors, requires re-analysis (complete Section 7)

**Rationale for Decision:** (2-3 sentences explaining your decision)

```
[Your rationale here - must be 10-1000 characters]
```

---

### TEMPLATE SECTION 2: PER-FINDING VALIDATION (REQUIRED)

For each major molecular finding in the draft report, indicate validation status. Review the top 10 findings listed in `draft_report.json` under `key_molecular_findings`.

#### Finding 1: {{GENE_NAME_1}} - {{FINDING_TYPE_1}}
- **Finding ID:** {{FINDING_ID_1}}
- **Confidence Score:** FDR = {{FDR_VALUE_1}}
- **Clinical Significance:** {{CLINICAL_SIGNIFICANCE_1}}
- **Validation Status:**
  - [ ] **CONFIRMED** - Finding is clinically valid and supported by data
  - [ ] **UNCERTAIN** - Requires additional validation or expert consultation
  - [ ] **INCORRECT** - Finding appears to be erroneous
- **Comments:** (Optional - explain your validation decision)
  ```
  [Your comments here]
  ```

#### Finding 2: {{GENE_NAME_2}} - {{FINDING_TYPE_2}}
- **Finding ID:** {{FINDING_ID_2}}
- **Confidence Score:** FDR = {{FDR_VALUE_2}}
- **Clinical Significance:** {{CLINICAL_SIGNIFICANCE_2}}
- **Validation Status:**
  - [ ] CONFIRMED
  - [ ] UNCERTAIN
  - [ ] INCORRECT
- **Comments:**
  ```
  [Your comments here]
  ```

#### Finding 3: {{GENE_NAME_3}} - {{FINDING_TYPE_3}}
- **Finding ID:** {{FINDING_ID_3}}
- **Confidence Score:** FDR = {{FDR_VALUE_3}}
- **Clinical Significance:** {{CLINICAL_SIGNIFICANCE_3}}
- **Validation Status:**
  - [ ] CONFIRMED
  - [ ] UNCERTAIN
  - [ ] INCORRECT
- **Comments:**
  ```
  [Your comments here]
  ```

#### Finding 4: {{GENE_NAME_4}} - {{FINDING_TYPE_4}}
- **Finding ID:** {{FINDING_ID_4}}
- **Confidence Score:** FDR = {{FDR_VALUE_4}}
- **Clinical Significance:** {{CLINICAL_SIGNIFICANCE_4}}
- **Validation Status:**
  - [ ] CONFIRMED
  - [ ] UNCERTAIN
  - [ ] INCORRECT
- **Comments:**
  ```
  [Your comments here]
  ```

#### Finding 5: {{GENE_NAME_5}} - {{FINDING_TYPE_5}}
- **Finding ID:** {{FINDING_ID_5}}
- **Confidence Score:** FDR = {{FDR_VALUE_5}}
- **Clinical Significance:** {{CLINICAL_SIGNIFICANCE_5}}
- **Validation Status:**
  - [ ] CONFIRMED
  - [ ] UNCERTAIN
  - [ ] INCORRECT
- **Comments:**
  ```
  [Your comments here]
  ```

#### Finding 6: {{GENE_NAME_6}} - {{FINDING_TYPE_6}}
- **Finding ID:** {{FINDING_ID_6}}
- **Confidence Score:** FDR = {{FDR_VALUE_6}}
- **Clinical Significance:** {{CLINICAL_SIGNIFICANCE_6}}
- **Validation Status:**
  - [ ] CONFIRMED
  - [ ] UNCERTAIN
  - [ ] INCORRECT
- **Comments:**
  ```
  [Your comments here]
  ```

#### Finding 7: {{GENE_NAME_7}} - {{FINDING_TYPE_7}}
- **Finding ID:** {{FINDING_ID_7}}
- **Confidence Score:** FDR = {{FDR_VALUE_7}}
- **Clinical Significance:** {{CLINICAL_SIGNIFICANCE_7}}
- **Validation Status:**
  - [ ] CONFIRMED
  - [ ] UNCERTAIN
  - [ ] INCORRECT
- **Comments:**
  ```
  [Your comments here]
  ```

#### Finding 8: {{GENE_NAME_8}} - {{FINDING_TYPE_8}}
- **Finding ID:** {{FINDING_ID_8}}
- **Confidence Score:** FDR = {{FDR_VALUE_8}}
- **Clinical Significance:** {{CLINICAL_SIGNIFICANCE_8}}
- **Validation Status:**
  - [ ] CONFIRMED
  - [ ] UNCERTAIN
  - [ ] INCORRECT
- **Comments:**
  ```
  [Your comments here]
  ```

#### Finding 9: {{GENE_NAME_9}} - {{FINDING_TYPE_9}}
- **Finding ID:** {{FINDING_ID_9}}
- **Confidence Score:** FDR = {{FDR_VALUE_9}}
- **Clinical Significance:** {{CLINICAL_SIGNIFICANCE_9}}
- **Validation Status:**
  - [ ] CONFIRMED
  - [ ] UNCERTAIN
  - [ ] INCORRECT
- **Comments:**
  ```
  [Your comments here]
  ```

#### Finding 10: {{GENE_NAME_10}} - {{FINDING_TYPE_10}}
- **Finding ID:** {{FINDING_ID_10}}
- **Confidence Score:** FDR = {{FDR_VALUE_10}}
- **Clinical Significance:** {{CLINICAL_SIGNIFICANCE_10}}
- **Validation Status:**
  - [ ] CONFIRMED
  - [ ] UNCERTAIN
  - [ ] INCORRECT
- **Comments:**
  ```
  [Your comments here]
  ```

---

### TEMPLATE SECTION 3: CLINICAL GUIDELINE COMPLIANCE (REQUIRED)

#### NCCN Guidelines Alignment
Review treatment recommendations against current NCCN (National Comprehensive Cancer Network) guidelines for this cancer type.

- [ ] **ALIGNED** - All recommendations follow NCCN guidelines
- [ ] **PARTIAL** - Some recommendations deviate from guidelines (explain below)
- [ ] **NOT_ALIGNED** - Significant deviations from NCCN guidelines (explain below)

**Deviations from NCCN Guidelines (if any):**
```
[List specific guideline deviations and justification]

Example:
- Recommendation for off-label alpelisib use: Justified by PIK3CA E545K mutation and platinum resistance despite Category 2B NCCN rating
```

#### Institutional Protocol Alignment
Review treatment recommendations against institutional formulary and clinical protocols.

- [ ] **ALIGNED** - All recommendations follow institutional protocols
- [ ] **PARTIAL** - Some recommendations deviate from protocols (explain below)
- [ ] **NOT_ALIGNED** - Significant deviations from institutional protocols (explain below)

**Deviations from Institutional Protocols (if any):**
```
[List specific protocol deviations and justification]

Example:
- Bevacizumab continuation: Not typically used in platinum-resistant setting per institutional protocol, but justified by absence of GI perforation risk factors
```

---

### TEMPLATE SECTION 4: QUALITY FLAGS ASSESSMENT

Review automated quality flags raised during analysis. Assess whether each flag is acceptable or requires action before approval.

#### Flag 1: {{FLAG_CHECK_NAME_1}}
- **Severity:** {{FLAG_SEVERITY_1}} (critical / warning / info)
- **Message:** {{FLAG_MESSAGE_1}}
- **Recommendation:** {{FLAG_RECOMMENDATION_1}}
- **Reviewer Assessment:**
  - [ ] **ACCEPTABLE** - Flag noted but not concerning, does not affect clinical validity
  - [ ] **REQUIRES_ACTION** - Flag must be addressed before approval
- **Comments:**
  ```
  [Explain your assessment of this quality flag]

  Example:
  Minimum region size of 35 spots is acceptable. While below ideal threshold of 50, the strong effect sizes (log2FC > 3) and very low FDR values (< 1e-10) provide sufficient confidence.
  ```

#### Flag 2: {{FLAG_CHECK_NAME_2}}
- **Severity:** {{FLAG_SEVERITY_2}}
- **Message:** {{FLAG_MESSAGE_2}}
- **Recommendation:** {{FLAG_RECOMMENDATION_2}}
- **Reviewer Assessment:**
  - [ ] ACCEPTABLE
  - [ ] REQUIRES_ACTION
- **Comments:**
  ```
  [Your assessment here]
  ```

#### Flag 3: {{FLAG_CHECK_NAME_3}}
- **Severity:** {{FLAG_SEVERITY_3}}
- **Message:** {{FLAG_MESSAGE_3}}
- **Recommendation:** {{FLAG_RECOMMENDATION_3}}
- **Reviewer Assessment:**
  - [ ] ACCEPTABLE
  - [ ] REQUIRES_ACTION
- **Comments:**
  ```
  [Your assessment here]
  ```

**If no quality flags were raised:**
No quality flags - all automated checks passed

---

### TEMPLATE SECTION 5: TREATMENT RECOMMENDATIONS REVIEW

Review each treatment recommendation for appropriateness given patient context, molecular profile, and clinical guidelines.

#### Recommendation 1: {{THERAPY_NAME_1}}
- **Recommendation ID:** {{REC_ID_1}}
- **Molecular Target:** {{MOLECULAR_TARGET_1}}
- **Expected Efficacy:** {{EXPECTED_EFFICACY_1}}
- **Agreement:**
  - [ ] **AGREE** - Recommendation is appropriate for this patient
  - [ ] **DISAGREE** - Recommendation is not appropriate (explain below)
- **Comments:**
  ```
  [Your assessment of this recommendation]

  Example:
  AGREE. PI3K inhibitor + PARP inhibitor combination is appropriate given PIK3CA E545K mutation and BRCA1 germline status. Patient meets eligibility criteria for clinical trial NCT12345678.
  ```

#### Recommendation 2: {{THERAPY_NAME_2}}
- **Recommendation ID:** {{REC_ID_2}}
- **Molecular Target:** {{MOLECULAR_TARGET_2}}
- **Expected Efficacy:** {{EXPECTED_EFFICACY_2}}
- **Agreement:**
  - [ ] AGREE
  - [ ] DISAGREE
- **Comments:**
  ```
  [Your assessment here]
  ```

#### Recommendation 3: {{THERAPY_NAME_3}}
- **Recommendation ID:** {{REC_ID_3}}
- **Molecular Target:** {{MOLECULAR_TARGET_3}}
- **Expected Efficacy:** {{EXPECTED_EFFICACY_3}}
- **Agreement:**
  - [ ] AGREE
  - [ ] DISAGREE
- **Comments:**
  ```
  [Your assessment here]
  ```

#### Recommendation 4: {{THERAPY_NAME_4}}
- **Recommendation ID:** {{REC_ID_4}}
- **Molecular Target:** {{MOLECULAR_TARGET_4}}
- **Expected Efficacy:** {{EXPECTED_EFFICACY_4}}
- **Agreement:**
  - [ ] AGREE
  - [ ] DISAGREE
- **Comments:**
  ```
  [Your assessment here]
  ```

---

### TEMPLATE SECTION 6: ATTESTATION & SIGNATURE (REQUIRED)

By completing this review, you attest to the following:

- [ ] **I have reviewed all findings in this report** (10 molecular findings validated)
- [ ] **I have assessed clinical guideline compliance** (NCCN + institutional protocols)
- [ ] **I have evaluated quality flags** (all automated QC flags assessed)
- [ ] **My decision reflects my clinical judgment** based on expertise and patient context
- [ ] **I understand this review is part of the medical record** (10-year retention, HIPAA-compliant)

**Reviewer Information:**

**Full Name:** {{REVIEWER_NAME}}
**Credentials:** {{CREDENTIALS}} (e.g., MD, Gynecologic Oncology)
**Email:** {{EMAIL}}
**Role:** {{ROLE}} (oncologist / pathologist / radiologist / bioinformatician)
**Date:** {{DATE}} (YYYY-MM-DD)
**Time:** {{TIME}} (HH:MM timezone)

**Digital Signature:**
*[Auto-generated SHA-256 hash by citl_submit_review.py - leave blank]*

---

### TEMPLATE SECTION 7: REVISION INSTRUCTIONS (Required if status = REVISE or REJECT)

**Complete this section ONLY if you selected REVISE or REJECT in Section 1.**

#### Issues to Address

List specific issues that must be corrected before resubmission:

1. **Issue 1:**
   ```
   [Specific issue description]

   Example:
   Exclude tumor_necrotic region (35 spots) due to insufficient sample size. FDR values for DEGs in this region are unreliable.
   ```

2. **Issue 2:**
   ```
   [Specific issue description]
   ```

3. **Issue 3:**
   ```
   [Specific issue description]
   ```

#### Re-analysis Parameters

Specify parameter adjustments for re-running the analysis:

**Regions to exclude:**
```
[List region IDs to exclude from analysis]

Example:
- tumor_necrotic (n=35 spots)
- stroma_artifact (n=28 spots)
```

**FDR threshold adjustment:**
```
[Specify new FDR threshold if needed]

Example:
- Change from FDR < 0.05 to FDR < 0.01 for more stringent criteria
```

**Minimum spots per region:**
```
[Specify minimum sample size threshold]

Example:
- Increase from 30 to 50 spots minimum per region
```

**Additional tests required:**
```
[List any additional analyses needed]

Example:
- Validate TP53 mutation with orthogonal sequencing method (Sanger)
- Request immunohistochemistry confirmation for Ki67 proliferation rate
```

#### Resubmission Timeline

**Expected resubmission date:** {{RESUBMISSION_DATE}} (YYYY-MM-DD)

**Estimated time for corrections:** {{ESTIMATED_TIME}} (hours/days/weeks)

---

### Template Notes for Reviewers

**This form structure:**
- Sections 1, 2, 3, 6: **REQUIRED** for all reviews
- Section 4: Required if quality flags exist, otherwise skip
- Section 5: Required if treatment recommendations exist
- Section 7: **REQUIRED** if status is REVISE or REJECT

**Converting to JSON:**
1. Complete all required sections in this markdown template
2. Use provided conversion script: `python scripts/convert_review_to_json.py`
3. OR manually structure as JSON following `shared/schemas/citl_review_schema.json`
4. Submit using: `python servers/mcp-patient-report/scripts/citl_submit_review.py --patient-id {{PATIENT_ID}} --review-file <path>`

**Typical completion time by status:**
- **APPROVE:** 20-25 minutes (straightforward validation)
- **REVISE:** 25-35 minutes (requires detailed revision instructions)
- **REJECT:** 30-40 minutes (requires detailed error documentation and escalation plan)

**Questions or issues?**
- Consult the guide sections above for detailed instructions
- Contact bioinformatics team for technical questions
- Contact clinic administrator for workflow questions

---

**Document Version:** 1.1
**Last Updated:** 2026-02-19
**Part of:** Precision Medicine MCP - Clinician-in-the-Loop Validation Workflow
**Contact:** bioinformatics@hospital.org
