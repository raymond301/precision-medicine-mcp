# Test Prompt: All MCP Servers

## Quick Server Validation Test

Use this prompt in Claude Desktop to verify all MCP servers are working correctly:

---

**PROMPT:**

I need to run a quick validation test for a Stage IV HGSOC patient analysis. Please perform the following checks in order and report which MCP server was called for each step:

1. **Clinical Data** - Retrieve patient demographics and treatment history for patient ID "P001"

2. **Reference Genome** - Fetch information about the hg38 reference genome

3. **Cancer Genomics** - Query TCGA for ovarian cancer (OV) cohort statistics

4. **Histology Image** - Fetch a sample histology image and extract basic features

5. **Cell Segmentation** - Get cell segmentation results for a tissue sample

6. **Spatial Transcriptomics** - Perform QC filtering on spatial transcriptomics data with default parameters

7. **ML Model** - Query available genomic ML models on Hugging Face

8. **Multi-omics Integration** - Validate multi-omics data files (check what file formats are expected)

9. **Workflow Status** - List available nf-core pipelines for genomics analysis

Please report:
- Which server was called for each step
- Whether the call succeeded
- A one-line summary of the result

---

## Expected Server Calls

If all servers are working correctly, you should see:

| Step | Expected Server | What It Does |
|------|----------------|--------------|
| 1 | **mockepic** | Retrieves clinical/EHR data |
| 2 | **fgbio** | Fetches reference genome info |
| 3 | **tcga** | Queries cancer genomics cohort |
| 4 | **openimagedata** | Fetches and processes histology images |
| 5 | **deepcell** | Performs cell segmentation |
| 6 | **spatialtools** | Processes spatial transcriptomics |
| 7 | **huggingface** | Lists available ML models |
| 8 | **multiomics** | Validates multi-omics data |
| 9 | **seqera** | Lists workflow pipelines |

## Success Criteria

✅ **PASS**: All servers are called and return synthetic/mock data
⚠️ **PARTIAL**: Some servers called but others failed
❌ **FAIL**: Servers not being called or errors occur

## Notes

- All servers are in DRY_RUN mode, so results will be synthetic/mocked
- Each result should include a warning banner indicating synthetic data
- Total execution time should be < 30 seconds
- This test does NOT require any real data files
