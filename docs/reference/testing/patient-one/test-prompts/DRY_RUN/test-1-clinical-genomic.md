TEST 1: Clinical Data and Genomic Analysis
===========================================

> **Prompt library:** This is the canonical prompt for [Clinical Workflow Prompt 1](../../../../prompts/clinical-workflow.md#prompt-1-clinical-data-and-genomic-profiling). See that page for clinical interpretation context and expected output summary.

> **Data Mode:** This test uses **synthetic data** — DRY_RUN mode for mcp-fgbio, plus mock servers (mcp-mockepic and mcp-mocktcga are always synthetic). For real TCGA cohort data, use the [cBioPortal external connector](../../../../../for-researchers/CONNECT_EXTERNAL_MCP.md). For real EHR data, use mcp-epic (local, HIPAA-compliant). See [Data Modes Guide](../../data-modes-guide.md) for details.

Patient ID: PAT001-OVC-2025

⚠️ IMPORTANT: Use the MCP server tools to read files. The data has been copied to MCP-accessible locations.

**GCS Data Location:** `gs://sample-inputs-patientone/patient-data/PAT001-OVC-2025/`

## PART 1: Clinical Data (use mcp-mockepic)

The patient data files are located in GCS at `gs://sample-inputs-patientone/patient-data/PAT001-OVC-2025/`. Use the mockepic server to retrieve clinical information. MCP servers on Cloud Run can access GCS URIs directly.

1. For patient PAT001-OVC-2025, retrieve:
   - Patient demographics (name, age, family history)
   - Genetic mutations noted in family history
   - BRCA1 status

2. Retrieve lab results for this patient:
   - CA-125 tumor marker trends over time
   - What does the CA-125 trajectory indicate about treatment response?
   - Evidence of platinum resistance in the lab data?

Files to read:
- patient_demographics.json (contains Sarah Anderson, age 58, BRCA1 germline mutation)
- lab_results.json (contains CA-125 values: 1456→22→389→289 U/mL showing resistance)

## PART 2: Genomic Analysis (use mcp-fgbio and mcp-mocktcga)

3. Parse somatic variants for patient PAT001-OVC-2025:
   Use fgbio to read the VCF file with these expected mutations:
   - TP53 R175H (chr17:7,578,406 C>A) - should be present
   - PIK3CA E545K (chr3:178,936,091 G>A) - should be present
   - PTEN LOH (chr10:89,692,940 G>T) - should be present

   Copy number alterations to look for:
   - MYC, CCNE1, AKT2 amplifications
   - RB1, CDKN2A deletions

4. Compare to TCGA-OV cohort (use mcp-mocktcga):
   For a patient with:
   - BRCA1 germline mutation
   - TP53 R175H somatic mutation
   - Stage IV HGSOC
   - Platinum-resistant progression

   Questions:
   - What TCGA molecular subtype does this match? (Expect C1 immunoreactive or C2 differentiated)
   - What is the typical prognosis for BRCA1-mutant, TP53-mutant Stage IV HGSOC?
   - What pathways are commonly activated in platinum-resistant ovarian cancer?

## Expected Results to Validate:

**Clinical:**
- Patient: Sarah Elizabeth Anderson
- Age: 58 years old
- BRCA1: Pathogenic germline mutation
- CA-125 at diagnosis: 1456 U/mL
- CA-125 after treatment: 22 U/mL (normalized)
- CA-125 at progression: 389 U/mL (resistance)
- Current CA-125: 289 U/mL

**Genomic:**
- TP53 mutation: R175H (hotspot)
- PIK3CA mutation: E545K (activating)
- PTEN: Loss of heterozygosity
- Copy number: MYC, CCNE1, AKT2 amplified

**TCGA:**
- Subtype: C1 or C2
- Prognosis: Poor with Stage IV + platinum resistance
- Activated pathways: PI3K/AKT/mTOR

## Output Format:
Please provide:
1. Patient summary (demographics, genetic risk)
2. CA-125 trend analysis (with interpretation of resistance)
3. Key somatic mutations identified
4. TCGA subtype and pathway analysis
5. Clinical significance summary
