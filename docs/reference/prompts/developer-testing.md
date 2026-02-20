# Developer Test Prompts

**Audience:** Developers, DevOps, QA Engineers
**Purpose:** Server validation, testing, debugging
**Time per prompt:** 1-15 minutes
**Expected output:** Technical verification, pass/fail results

---

## Server Validation Prompts (9+ Prompts)

### Prompt 1: List All Available Tools (Quick Smoke Test)

**Time:** <1 minute | **Servers:** Any | **Output:** Tool inventory

```
List all available tools from the spatialtools server.

Expected output:
- Tool names
- Tool descriptions
- Parameters for each tool
- Total tool count (should be 14 tools)
```

**Expected Output:**
```
Found 14 tools in mcp-spatialtools:
1. align_spatial_data - Align FASTQ to reference genome using STAR
2. perform_batch_correction - ComBat batch effect correction
3. deconvolve_cell_types - Estimate cell type proportions
4. perform_differential_expression - Compare gene expression between groups
5. perform_pathway_enrichment - Find enriched pathways
6. calculate_spatial_autocorrelation - Moran's I for spatial patterns
7. generate_spatial_heatmap - Generate spatial heatmaps
... (14 total)
```

**Pass Criteria:**
- ✅ Returns tool list without errors
- ✅ Tool count matches expected (14 for spatialtools)
- ✅ Each tool has description and parameters

---

### Prompt 2: Validate mcp-multiomics Server Configuration

**Time:** 2 minutes | **Servers:** mcp-multiomics | **Output:** Server health check

```
Test the mcp-multiomics server by:

1. List all available tools (should be 10 tools)
2. Check tool categories:
   - Data preprocessing (validate_multiomics_data, preprocess_multiomics_data, visualize_data_quality)
   - Integration (integrate_omics_data)
   - Association testing (run_halla_analysis, calculate_stouffer_meta)
   - Visualization (create_multiomics_heatmap, run_multiomics_pca)
   - Interpretation (predict_upstream_regulators)
   - Utilities (estimate_analysis_cost)

3. Verify each tool has:
   - Input parameters defined
   - Output schema documented
   - Example usage

4. Test dry-run mode:
   - Call integrate_omics_data with mock parameters
   - Should return success with simulated data
```

**Expected Output:**
```
✅ Server: mcp-multiomics (RUNNING)
✅ Tools found: 10
✅ All tools have complete documentation
✅ Dry-run mode: FUNCTIONAL

Tool Categories:
- Preprocessing: 3 tools ✓ (validate, preprocess, visualize_data_quality)
- Integration: 1 tool ✓ (integrate_omics_data)
- Association: 2 tools ✓ (run_halla_analysis, calculate_stouffer_meta)
- Visualization: 2 tools ✓ (create_multiomics_heatmap, run_multiomics_pca)
- Interpretation: 1 tool ✓ (predict_upstream_regulators)
- Utilities: 1 tool ✓ (estimate_analysis_cost)
```

**Pass Criteria:**
- ✅ 10 tools available
- ✅ All tools documented
- ✅ Dry-run call succeeds

---

### Prompt 3: Validate mcp-spatialtools with PatientOne Data

**Time:** 5 minutes | **Servers:** mcp-spatialtools | **Output:** Data processing validation

```
Validate mcp-spatialtools using PatientOne spatial data:

1. **Data Loading:**
   - Load: /data/PAT001-OVC-2025/spatial/visium_gene_expression.csv
   - Load: /data/PAT001-OVC-2025/spatial/visium_spatial_coordinates.csv
   - Verify: 900 spots, 31 genes
   - Check: No missing values in coordinates

2. **Data Processing:**
   - Run spatial autocorrelation on MKI67 gene
   - Expected: Moran's I > 0.5 (spatial clustering)
   - Runtime: <10 seconds

3. **Output Validation:**
   - Returns: {"gene": "MKI67", "morans_i": float, "p_value": float}
   - morans_i: Between 0 and 1
   - p_value: < 0.05 (significant)

4. **Error Handling:**
   - Test invalid gene name: Should return helpful error
   - Test missing coordinates: Should catch and report
```

**Expected Output:**
```
✅ Data loaded: 900 spots × 31 genes
✅ Coordinates validated: No missing values
✅ Spatial autocorrelation: MKI67, Moran's I = 0.623, p = 1.2e-08
✅ Runtime: 4.2 seconds
✅ Error handling: Graceful failure on invalid inputs

PASS: mcp-spatialtools validation complete
```

**Pass Criteria:**
- ✅ Data loads correctly
- ✅ Analysis completes without errors
- ✅ Output format matches schema
- ✅ Runtime < 10 seconds
- ✅ Error handling works

---

### Prompt 4: Validate mcp-fgbio Reference Genome Tool

**Time:** 3 minutes | **Servers:** mcp-fgbio | **Output:** Tool functionality test

```
Test mcp-fgbio reference genome functionality:

1. **List Available Resources:**
   - Query: reference://hg38
   - Expected: Metadata for GRCh38 genome
   - Check: genome_id, name, size, URL fields present

2. **Query Gene Annotations:**
   - Gene: BRCA1
   - Genome: hg38
   - Source: gencode
   - Expected: chr17:43,044,295-43,170,245

3. **Validate FASTQ Tool (dry-run):**
   - Call validate_fastq with mock file path
   - Expected: Returns quality metrics (even if simulated)

4. **Error Handling:**
   - Test: Non-existent gene name
   - Test: Invalid genome ID
   - Both should return clear error messages
```

**Expected Output:**
```
✅ Resource query: hg38 metadata retrieved
   - genome_id: hg38
   - name: Human GRCh38
   - size: 3.0 Gb
   - chromosomes: 24

✅ Gene annotation: BRCA1
   - chromosome: chr17
   - start: 43,044,295
   - end: 43,170,245
   - strand: -

✅ FASTQ validation: Tool functional (dry-run mode)

✅ Error handling: Clear error messages for invalid inputs

PASS: mcp-fgbio validation complete
```

**Pass Criteria:**
- ✅ Resource metadata correct
- ✅ Gene coordinates accurate
- ✅ Tools execute without crashes
- ✅ Error messages are helpful

---

### Prompt 5: Test mcp-epic FHIR De-identification

**Time:** 5 minutes | **Servers:** mcp-epic or mcp-mockepic | **Output:** HIPAA compliance test

```
Test FHIR de-identification (HIPAA Safe Harbor):

1. **Load Test Patient:**
   - Patient ID: PAT001-OVC-2025 (or mock patient)
   - Retrieve: Demographics, conditions, medications

2. **Verify 18 Identifiers Removed:**
   - Names: Replaced with generic "Patient [ID]"
   - Dates: Year only (no month/day)
   - Geographic: State level only (no city/zip)
   - Phone/email: Removed
   - SSN/MRN: Removed
   - Ages >89: Set to 90

3. **Output Validation:**
   - De-identified FHIR resource returned
   - Original vs de-identified comparison
   - Audit log entry created

4. **HIPAA Compliance Check:**
   - All 18 Safe Harbor identifiers removed? ✓
   - Data still usable for research? ✓
   - Re-identification risk minimal? ✓
```

**Expected Output:**
```
✅ Patient loaded: PAT001-OVC-2025

Original FHIR Resource:
- name: Sarah Elizabeth Anderson
- birthDate: 1967-04-15
- address: 123 Main St, Boston, MA 02101
- telecom: 617-555-1234, sarah.anderson@email.com

De-identified FHIR Resource:
- name: Patient PAT001
- birthDate: 1967 (year only)
- address: Massachusetts (state only)
- telecom: [removed]

✅ 18 Safe Harbor identifiers: ALL REMOVED
✅ Audit log: Entry created (timestamp, user, patient_id)
✅ HIPAA compliance: PASS

PASS: De-identification validation complete
```

**Pass Criteria:**
- ✅ All 18 identifiers removed
- ✅ Data remains usable
- ✅ Audit log created
- ✅ No PHI in output

---

### Prompt 6: Load Test - Concurrent Analysis Requests

**Time:** 5-10 minutes | **Servers:** All | **Output:** Performance metrics

```
Perform load testing with concurrent requests:

**Test Scenario:**
- Simulate 5 concurrent users
- Each requests spatial pathway enrichment analysis
- Same data (PatientOne), different gene lists

**Metrics to Measure:**
1. Response time (per request)
2. Throughput (requests/second)
3. Error rate (%)
4. Peak memory usage
5. Server auto-scaling behavior

**Expected Performance:**
- Response time: <30 seconds (p95)
- Throughput: >0.1 requests/second
- Error rate: <1%
- Memory: <2GB per request
- Auto-scaling: Triggers at >80% CPU

**Test Commands:**
Run 5 pathway enrichment queries in parallel:
- User 1: PI3K/AKT genes
- User 2: Apoptosis genes
- User 3: Cell cycle genes
- User 4: Immune response genes
- User 5: Hypoxia genes
```

**Expected Output:**
```
Load Test Results (5 concurrent users):

Response Times:
- Min: 12.3 sec
- Mean: 18.7 sec
- p95: 24.1 sec
- Max: 27.9 sec

Throughput: 0.27 requests/sec
Error Rate: 0% (0/5)
Peak Memory: 1.2 GB
Auto-scaling: Not triggered (CPU 65%)

✅ PASS: All requests completed successfully
✅ PASS: p95 response time < 30 sec
✅ PASS: Error rate < 1%
✅ PASS: Memory usage < 2GB

PASS: Load test validation complete
```

**Pass Criteria:**
- ✅ All requests complete
- ✅ p95 response < 30s
- ✅ Error rate <1%
- ✅ No crashes

---

### Prompt 7: Integration Test - Multi-Server Workflow

**Time:** 10-15 minutes | **Servers:** mcp-mockepic, mcp-multiomics, mcp-spatialtools | **Output:** End-to-end test

```
Test complete multi-server workflow for PatientOne:

**Workflow Steps:**
1. **Clinical Data** (mcp-mockepic):
   - Retrieve patient demographics
   - Get CA-125 values
   - Expected: BRCA1 mutation confirmed

2. **Multi-Omics** (mcp-multiomics):
   - Load RNA/Protein/Phospho data
   - Run Stouffer's meta-analysis
   - Expected: PIK3CA, AKT1, MTOR significant (q<0.05)

3. **Spatial Analysis** (mcp-spatialtools):
   - Load Visium data (900 spots)
   - Calculate spatial autocorrelation for PIK3CA
   - Expected: Moran's I > 0.5

4. **Cross-Server Consistency:**
   - PIK3CA appears in multi-omics AND spatial?
   - Clinical context matches molecular findings?

5. **Output Synthesis:**
   - Generate integrated summary
   - All 3 servers contributed data?
   - Findings are concordant?
```

**Expected Output:**
```
Integration Test: PatientOne Multi-Server Workflow

Step 1 - Clinical Data (mcp-mockepic):
✅ Patient demographics retrieved
✅ BRCA1 germline mutation: CONFIRMED
✅ CA-125: 1456 U/mL at diagnosis

Step 2 - Multi-Omics (mcp-multiomics):
✅ Data loaded: 14 samples (7 resistant, 7 sensitive)
✅ Stouffer's meta-analysis completed
✅ Significant genes: PIK3CA (q=0.0001), AKT1 (q<0.0001), MTOR (q=0.0003)

Step 3 - Spatial Analysis (mcp-spatialtools):
✅ Visium data loaded: 900 spots × 31 genes
✅ Spatial autocorrelation: PIK3CA, Moran's I = 0.58, p=3.2e-07

Step 4 - Cross-Server Consistency:
✅ PIK3CA: Found in multi-omics AND spatial
✅ Clinical-molecular concordance: BRCA1 (clinical) + PI3K pathway (molecular)

Step 5 - Integrated Summary:
✅ All 3 servers contributed data
✅ Findings are concordant
✅ Treatment target identified: PI3K/AKT pathway

PASS: Multi-server integration test complete
```

**Pass Criteria:**
- ✅ All 3 servers respond
- ✅ Data flows between servers
- ✅ Findings are consistent
- ✅ No errors

---

### Prompt 8: Error Handling and Recovery Test

**Time:** 5 minutes | **Servers:** All | **Output:** Robustness validation

```
Test error handling and recovery mechanisms:

**Test Cases:**

1. **Invalid Patient ID:**
   - Request data for non-existent patient "PAT999-XXX-9999"
   - Expected: Clear error message, no crash
   - Error message should include: "Patient ID not found. Valid IDs: ..."

2. **Missing Data File:**
   - Request spatial analysis with missing coordinates file
   - Expected: File not found error, suggest troubleshooting steps
   - Should not proceed with partial data

3. **Malformed Input:**
   - Pass string where integer expected (e.g., fdr_threshold="high" instead of 0.05)
   - Expected: Type error with clear correction guidance

4. **Timeout Scenario:**
   - Request very large analysis (if possible, simulate)
   - Expected: Timeout error with suggestion to reduce scope or increase timeout

5. **API Rate Limit:**
   - Make rapid successive requests (if rate-limited)
   - Expected: Rate limit error with retry-after information

**Recovery Testing:**
- After each error, verify server remains responsive
- Subsequent valid requests should succeed
```

**Expected Output:**
```
Error Handling Test Results:

Test 1 - Invalid Patient ID:
✅ Error message clear: "Patient PAT999-XXX-9999 not found"
✅ Suggested valid IDs
✅ Server remained responsive

Test 2 - Missing Data File:
✅ File not found error: "/data/missing.csv"
✅ Troubleshooting suggestion provided
✅ Did not proceed with partial data

Test 3 - Malformed Input:
✅ Type error detected: "fdr_threshold must be float, got str"
✅ Correction guidance: "Use fdr_threshold=0.05 (numeric)"
✅ Input validation prevented invalid execution

Test 4 - Timeout:
✅ Timeout error after 300 seconds
✅ Suggestion: "Reduce gene list or increase timeout parameter"
✅ Partial results saved for debugging

Test 5 - Rate Limit:
✅ Rate limit error: "Maximum 10 requests per minute exceeded"
✅ Retry-after: 42 seconds
✅ Subsequent requests succeeded after waiting

Recovery:
✅ All subsequent valid requests succeeded
✅ No lingering errors
✅ Server performance unaffected

PASS: Error handling and recovery test complete
```

**Pass Criteria:**
- ✅ All errors caught gracefully
- ✅ Error messages are helpful
- ✅ Server recovers after errors
- ✅ No crashes

---

### Prompt 9: Data Quality and Integrity Check

**Time:** 10 minutes | **Servers:** mcp-multiomics, mcp-spatialtools | **Output:** Data validation report

```
Validate data quality across all PatientOne datasets:

**Genomic Data (VCF):**
1. Check VCF format compliance
2. Verify required fields (CHROM, POS, REF, ALT, QUAL, FILTER)
3. Validate variant coordinates against hg38
4. Check for duplicate variants

**Multi-Omics Data (CSV):**
1. Check sample naming consistency (RNA, Protein, Phospho)
2. Verify metadata completeness (Batch, Condition columns)
3. Detect batch effects (PCA, PC1 vs Batch correlation)
4. Check for outlier samples (MAD > 3)
5. Report missing value percentages

**Spatial Data (Visium):**
1. Validate spot count (should be ~900 for 10x Visium)
2. Check coordinate range (x: 0-6000, y: 0-6000 typical)
3. Verify gene count (31 genes for PatientOne)
4. Check for duplicate spots (same x,y coordinates)

**Output Quality Report:**
- Data completeness: % non-missing
- Format compliance: PASS/FAIL
- Outliers detected: Count and sample IDs
- Recommendations: Any data cleaning needed
```

**Expected Output:**
```
Data Quality Report: PatientOne (PAT001-OVC-2025)

=== Genomic Data (VCF) ===
✅ VCF format: VALID (VCF 4.2)
✅ Required fields: Present
✅ Variants: 4 total, all valid hg38 coordinates
✅ Duplicates: None
PASS: Genomic data quality check

=== Multi-Omics Data ===
✅ Sample naming: Consistent across 3 modalities
✅ Metadata: Complete (Batch, Condition, Treatment)
⚠️  Batch effects detected: PC1-Batch correlation = 0.82 (proteomics)
   Recommendation: Apply ComBat batch correction
⚠️  Outlier detected: Sample_07 (MAD = 3.5)
   Recommendation: Exclude from analysis
✅ Missing values: RNA 0%, Protein 32%, Phospho 38%
   Recommendation: Apply KNN imputation
PARTIAL PASS: Multi-omics data needs preprocessing

=== Spatial Data (Visium) ===
✅ Spot count: 900 (expected for 10x Visium)
✅ Coordinates: Valid range (x: 0-5800, y: 0-5200)
✅ Genes: 31 (matches PatientOne dataset)
✅ Duplicates: None
PASS: Spatial data quality check

=== Summary ===
Overall: PASS with preprocessing recommendations
- Genomic: Ready for analysis
- Multi-omics: Apply batch correction and imputation
- Spatial: Ready for analysis

Recommended preprocessing steps:
1. Multi-omics: Run preprocess_multiomics_data tool
2. Multi-omics: Exclude Sample_07 (outlier)
3. Multi-omics: Verify batch correction (PC1-Batch < 0.3)
```

**Pass Criteria:**
- ✅ All format checks pass
- ✅ Data completeness >90%
- ✅ Issues identified and documented
- ✅ Recommendations provided

---

## Additional Developer Prompts

### Prompt 10: Cost and Token Usage Monitoring

**Time:** 2 minutes | **Output:** Usage metrics

```
Track Claude API usage for PatientOne complete workflow:

1. Run complete TEST_1 through TEST_6 workflow
2. Monitor:
   - Total tokens (input + output)
   - Cost per test
   - Total workflow cost
3. Compare to budget targets in [Cost Analysis](../shared/cost-analysis.md)

Expected:
- TEST_1: ~8K tokens, $0.05
- TEST_2: ~15K tokens, $0.10
- TEST_3: ~10K tokens, $0.06
- TEST_4: ~12K tokens, $0.08
- TEST_5: ~8K tokens, $0.05
- TEST_6: ~5K tokens, $0.03
- Total: ~58K tokens, $0.37
```

---

### Prompt 11: Audit Log Verification

**Time:** 3 minutes | **Output:** Audit compliance check

```
Verify audit logging for HIPAA compliance:

1. Trigger patient data access (PAT001-OVC-2025)
2. Check audit log entries:
   - Timestamp (ISO 8601 format)
   - User ID
   - Patient ID
   - Action (read/write/delete)
   - Resource type (FHIR Patient, VCF, etc.)
   - IP address (if available)
3. Verify 10-year retention policy configured
4. Test log immutability (attempt to modify log)

Expected: All access logged, logs immutable, retention set
```

---

## Summary: Developer Testing Checklist

### Pre-Deployment Tests
- [ ] Prompt 1: List tools (all servers)
- [ ] Prompt 2: Validate server configurations
- [ ] Prompt 3-5: Test individual servers
- [ ] Prompt 9: Data quality validation

### Integration Tests
- [ ] Prompt 7: Multi-server workflow
- [ ] Prompt 8: Error handling
- [ ] Prompt 6: Load testing

### Compliance Tests
- [ ] Prompt 5: FHIR de-identification
- [ ] Prompt 11: Audit logging

### Monitoring Tests
- [ ] Prompt 10: Cost tracking

**All tests should PASS before production deployment.**

---

**Document Version:** 1.0
**Date:** 2026-01-16
**Target Audience:** Developers, DevOps, QA Engineers
