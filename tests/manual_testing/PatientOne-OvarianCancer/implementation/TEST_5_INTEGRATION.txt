TEST 5: Integrated Analysis & Clinical Recommendations
=======================================================

Patient ID: PAT001-OVC-2025

⚠️ NOTE: Run this test AFTER completing Tests 1-4

This test synthesizes findings from all previous tests - NO new data loading required.

## Integration & Clinical Recommendations

Based on the findings from Tests 1-4, please synthesize a comprehensive clinical report.

### Reference Results from Previous Tests:

**TEST 1 - Clinical & Genomic:**
- Patient: Sarah Anderson, 58yo, BRCA1 germline mutation
- CA-125 trajectory: 1456 → 22 → 389 → 289 U/mL (platinum resistance)
- Somatic mutations: TP53 R175H, PIK3CA E545K, PTEN LOH
- TCGA subtype: C1/C2 (poor prognosis with Stage IV + resistance)

**TEST 2 - Multi-Omics:**
- PI3K/AKT pathway: ACTIVATED across RNA/Protein/Phospho
- Key genes dysregulated: PIK3CA, AKT1, MTOR (up); PTEN (down)
- Drug resistance: ABCB1 upregulated (MDR1 efflux pump)
- Anti-apoptotic: BCL2L1 upregulated

**TEST 3 - Spatial Transcriptomics:**
- 900 spots across 6 regions
- Resistance markers: Concentrated in tumor regions (heterogeneous)
- Immune cells: EXCLUDED from tumor (located in stroma_immune region)
- Tumor microenvironment: Immunologically "COLD"

**TEST 4 - Imaging:**
- High proliferation (Ki67 ~45-55%)
- Low CD8+ infiltration (~5-15 cells/mm²)
- TP53-mutant cells are highly proliferative
- Immune exclusion phenotype confirmed

---

## Analysis Questions:

### 1. PRIMARY RESISTANCE MECHANISMS (Rank by Evidence Strength)

Based on ALL modalities (genomics, multi-omics, spatial, imaging), identify the 2-3 main mechanisms driving platinum resistance in this patient.

For each mechanism, provide:
- **Mechanism name**
- **Supporting evidence** (which tests/modalities show this?)
- **Strength of evidence** (High/Medium/Low)
- **Therapeutic implications**

Expected top mechanisms:
1. PI3K/AKT/mTOR pathway activation
2. Drug efflux (ABCB1/MDR1)
3. Anti-apoptotic signaling (BCL2L1)
4. TP53 loss of function

### 2. MULTI-MODAL CONSISTENCY

Which molecular alterations appear consistently across multiple data types?

Create a cross-reference table:

| Feature | Genomics | Multi-Omics | Spatial | Imaging | Consistent? |
|---------|----------|-------------|---------|---------|-------------|
| TP53 mutation/loss | TP53 R175H | ? | ? | TP53+ cells | Yes/No |
| PI3K/AKT activation | PIK3CA E545K | PIK3CA/AKT1 up | ? | ? | Yes/No |
| High proliferation | ? | ? | MKI67 high | Ki67 high | Yes/No |
| Immune exclusion | ? | ? | CD8 in stroma | CD8 at periphery | Yes/No |

### 3. THERAPEUTIC RECOMMENDATIONS

Based on the integrated data, provide:

**A. Targeted Therapy Recommendations (Top 3):**

For each recommendation:
- **Drug/class**
- **Molecular target**
- **Supporting evidence** (from which tests?)
- **Expected efficacy** (High/Medium/Low)
- **FDA approval status** for HGSOC

Expected recommendations should include:
- PI3K/AKT/mTOR inhibitors (e.g., alpelisib + olaparib)
- MDR1 inhibitors or chemotherapy modifications
- BCL2L1 inhibitors (if available)
- PARP inhibitors (BRCA1 mutation)

**B. Immunotherapy Consideration:**

- Should checkpoint inhibitors be considered? (Yes/No)
- Why or why not? (Cite spatial and imaging evidence)
- Expected response rate based on immune phenotype
- Could combination with other agents overcome exclusion?

**C. Clinical Trial Opportunities:**

For BRCA1-mutant, platinum-resistant, Stage IV HGSOC with:
- PI3K/AKT pathway activation
- Immune exclusion phenotype
- High proliferation

Suggest trial types:
- PARP + PI3K/AKT inhibitor combinations
- Novel immunotherapy combinations
- BRCA-targeted therapies

### 4. BIOMARKERS FOR MONITORING

**A. Molecular Biomarkers:**
- Which genes/proteins should be tracked for treatment response?
- How often should they be monitored?
- What change indicates resistance?

**B. Imaging Biomarkers:**
- Which imaging features predict resistance?
- Can spatial transcriptomics track treatment response?
- Should Ki67 or immune infiltration be monitored?

**C. Clinical Biomarkers:**
- CA-125 trends (what trajectory indicates response vs resistance?)
- RECIST criteria
- Circulating tumor DNA

---

## Output Format:

Please provide a **concise 1-2 page clinical report** with:

### Executive Summary (3-4 sentences)
Brief overview of patient case and key findings

### Section 1: Resistance Mechanisms
**Ranked by evidence strength:**
1. [Mechanism] - Evidence: [tests], Strength: [High/Medium/Low]
2. [Mechanism] - Evidence: [tests], Strength: [High/Medium/Low]
3. [Mechanism] - Evidence: [tests], Strength: [High/Medium/Low]

### Section 2: Multi-Modal Consistency
Table showing which findings are consistent across modalities

### Section 3: Treatment Recommendations
**A. Targeted Therapies (Ranked by expected efficacy):**
1. [Drug/class] - Target: [gene/pathway], Evidence: [tests]
2. [Drug/class] - Target: [gene/pathway], Evidence: [tests]
3. [Drug/class] - Target: [gene/pathway], Evidence: [tests]

**B. Immunotherapy:**
- Recommendation: Yes/No
- Rationale: [cite immune exclusion evidence]
- Combination strategies: [if applicable]

**C. Clinical Trials:**
- Trial type 1: [description]
- Trial type 2: [description]

### Section 4: Monitoring Strategy
**Molecular:** [genes/proteins to track]
**Imaging:** [features to monitor]
**Clinical:** [CA-125, RECIST, ctDNA]

### Section 5: Multi-Modal Visualization Synthesis ⭐ NEW
**Spatial Visualizations:** [List visualization files generated]
**Imaging Visualizations:** [List visualization files generated]
**Cross-Modal Analysis:** [Key findings from comparing visualizations]
**Integrated Figure:** [4-panel figure with caption]

---

## Detailed Visualization Instructions:

Generate integrated visualizations to synthesize findings across all modalities:

**A. Spatial Gene Expression Visualizations** (use mcp-spatialtools):

1. **Spatial Heatmaps for Key Resistance Genes**
   ```
   Use tool: generate_spatial_heatmap

   Input:
   - Expression file: /data/spatial/filtered_expression.csv
   - Coordinates file: /data/spatial/spatial_coordinates.csv
   - Genes: ["PIK3CA", "AKT1", "MTOR", "PTEN", "TP53", "MKI67"]
   - Colormap: "RdYlBu_r" (red=high, blue=low)

   Output: Multi-panel spatial plot showing:
   - PIK3CA/AKT1/MTOR: HIGH expression in tumor regions (resistance)
   - PTEN: LOW expression in tumor (loss of tumor suppressor)
   - TP53: Mutant/loss pattern
   - MKI67: High proliferation index
   ```

2. **Gene Expression Heatmap Matrix**
   ```
   Use tool: generate_gene_expression_heatmap

   Input:
   - Expression file: /data/spatial/filtered_expression.csv
   - Regions file: /data/spatial/region_annotations.csv
   - Genes: Resistance signature (PIK3CA, AKT1, MTOR, ABCB1, BCL2L1, PTEN, TP53)

   Output: Heatmap showing gene expression differences across:
   - Tumor core vs margin vs stroma
   - Identifies tumor heterogeneity
   ```

3. **Spatial Autocorrelation Visualization**
   ```
   Use tool: visualize_spatial_autocorrelation

   Input: Autocorrelation results from TEST_3
   Output: Bar plot showing genes with significant spatial clustering
   - Moran's I > 0.3 indicates spatial co-expression patterns
   ```

**B. Imaging Visualizations** (use mcp-openimagedata and mcp-deepcell):

4. **H&E Morphology Annotation**
   ```
   Use tool: generate_he_annotation (mcp-openimagedata)

   Input:
   - H&E image: /data/imaging/HE_PAT001.tif
   - Necrotic regions: Identified from TEST_4 (red dashed boxes)
   - High cellularity regions: Tumor regions (green solid boxes)

   Output: Annotated H&E showing:
   - Necrotic areas (tissue architecture disruption)
   - High-grade serous morphology
   - Correlation with spatial gene expression (overlay regions)
   ```

5. **Multiplex IF Composite** (TP53 + Ki67 + CD8)
   ```
   Use tool: generate_multiplex_composite (mcp-openimagedata)

   Input:
   - Channels:
     * DAPI (nuclear, blue)
     * Ki67 (proliferation, green)
     * TP53 (tumor suppressor, red)
     * CD8 (immune cells, cyan)
   - Normalize: True

   Output: RGB composite showing:
   - Yellow cells = TP53+/Ki67+ (highly proliferative mutant cells)
   - Cyan cells = CD8+ T cells (mostly at periphery)
   - Spatial relationship: immune exclusion phenotype confirmed
   ```

6. **Cell Segmentation Overlay**
   ```
   Use tool: generate_segmentation_overlay (mcp-deepcell)

   Input: Nuclear/membrane channel from IF
   Output: Segmentation masks overlaid on original image
   - Quantify cell density per region
   - Validate CD8+ cell counts (~5-15 cells/mm²)
   ```

**C. Cross-Modal Synthesis Questions:**

Answer these by analyzing all generated visualizations together:

1. **Spatial-Imaging Correlation:**
   - Do regions with high PIK3CA/AKT1 expression (spatial) correspond to TP53+/Ki67+ cells (imaging)?
   - Expected: YES - tumor core shows both high resistance gene expression AND high proliferation

2. **Immune Exclusion Validation:**
   - Compare CD8 spatial expression (TEST_3) with CD8+ IF staining (TEST_4)
   - Expected: CD8 cells EXCLUDED from tumor, located in stroma_immune region
   - Visualize: Overlay spatial CD8 heatmap on H&E annotation

3. **Pathway-Morphology Link:**
   - Do regions with activated PI3K/AKT pathway (multi-omics) show poor differentiation (H&E)?
   - Expected: YES - high-grade serous morphology correlates with pathway activation

4. **Heterogeneity Assessment:**
   - Compare spatial gene expression heterogeneity with imaging-based cell phenotype distribution
   - Identify regions of treatment-refractory tumor (high resistance signature + low immune infiltration)

**D. Generate Integrated Multi-Modal Figure:**

Create a comprehensive figure combining:
- **Panel A:** Spatial heatmap (PIK3CA/AKT1/MTOR - resistance genes)
- **Panel B:** H&E annotation (morphology + region overlays)
- **Panel C:** Multiplex IF composite (TP53/Ki67/CD8 co-localization)
- **Panel D:** Gene expression heatmap (tumor vs stroma comparison)

**Caption:** "Multi-modal analysis of Patient PAT001-OVC-2025 reveals spatially heterogeneous PI3K/AKT pathway activation (Panel A) in high-grade serous tumor regions (Panel B), with proliferative TP53-mutant cells (Panel C, yellow) and immune cell exclusion (Panel C, cyan at periphery). Gene expression patterns (Panel D) confirm resistance signature enrichment in tumor vs stroma."

### Section 6: Prognosis
Based on TCGA data and integrated findings, expected outcomes with:
- Standard platinum re-challenge: [expected response]
- Recommended targeted therapies: [expected response]
- Novel immunotherapy: [expected response]

---

## Validation Checkpoints:

✅ Synthesized: Findings from all 4 previous tests
✅ Resistance mechanisms: Identified and ranked by evidence
✅ Multi-modal consistency: Confirmed across genomics/multi-omics/spatial/imaging
✅ Targeted therapies: Prioritized based on molecular evidence
✅ Immunotherapy: Recommendation based on immune phenotype (COLD)
✅ Monitoring strategy: Molecular, imaging, and clinical biomarkers defined
✅ **Visualization Synthesis (NEW):** Generated spatial heatmaps, H&E annotations, IF composites, and segmentation overlays
✅ **Cross-Modal Correlation (NEW):** Confirmed spatial gene expression correlates with imaging phenotypes
✅ **Integrated Figure (NEW):** Created 4-panel multi-modal figure synthesizing all findings
✅ Prognosis: Realistic expectations based on TCGA cohort data
