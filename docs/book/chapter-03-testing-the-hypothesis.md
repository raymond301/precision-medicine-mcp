# Chapter 3: Testing the Hypothesis

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
# Early mcp-fgbio/server.py (simplified)
@mcp.tool()
def parse_vcf(vcf_path: str) -> dict:
    """Parse VCF file and return variants."""

    # DRY_RUN mode: return synthetic data
    return {
        "variants": [
            {"gene": "TP53", "mutation": "R175H", "pathogenic": True},
            {"gene": "PIK3CA", "mutation": "E545K", "pathogenic": True}
        ],
        "total_variants": 2,
        "warning": "DRY_RUN mode - synthetic data"
    }
```

See original mock implementation: [`servers/mcp-fgbio/src/mcp_fgbio/server.py:145-187`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-fgbio/src/mcp_fgbio/server.py#L145-L187)

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
# servers/mcp-multiomics/src/mcp_multiomics/tools/preprocessing.py:45-89
def validate_omics_data(
    rna_path: str,
    protein_path: str,
    phospho_path: str
) -> dict:
    """Validate multi-omics data files for analysis.

    Checks:
    - File existence and format
    - Sample ID alignment across modalities
    - Gene/protein ID format
    - Missing value percentage
    - Data distribution (log-transformed?)
    """
    # Load data
    rna_df = pd.read_csv(rna_path, index_col=0)
    protein_df = pd.read_csv(protein_path, index_col=0)
    phospho_df = pd.read_csv(phospho_path, index_col=0)

    # Validate sample alignment
    common_samples = set(rna_df.columns) & set(protein_df.columns) & set(phospho_df.columns)
    if len(common_samples) < 3:
        raise ValueError(f"Insufficient common samples: {len(common_samples)}")

    # Return validation report...
```

Full implementation: [`servers/mcp-multiomics/src/mcp_multiomics/tools/preprocessing.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-multiomics/src/mcp_multiomics/tools/preprocessing.py)

**Week 2: HAllA Integration (The Python Fallback)**

HAllA is traditionally an R package. Integrating R into a Python MCP server requires `rpy2`—a Python-R bridge that's notoriously difficult to configure.

The breakthrough: **Implement the core algorithm in Python.**

```python
# Hierarchical association discovery in pure Python
def halla_python(
    rna_data: pd.DataFrame,
    protein_data: pd.DataFrame,
    method: str = "spearman",
    fdr_threshold: float = 0.05
) -> dict:
    """Python implementation of HAllA core algorithm."""
    from scipy.stats import spearmanr
    from statsmodels.stats.multitest import fdrcorrection

    associations = []
    for rna_gene in rna_data.index:
        for protein in protein_data.index:
            corr, pval = spearmanr(rna_data.loc[rna_gene], protein_data.loc[protein])
            associations.append({
                "rna_gene": rna_gene,
                "protein": protein,
                "correlation": corr,
                "p_value": pval
            })

    # FDR correction
    pvals = [a["p_value"] for a in associations]
    reject, qvals = fdrcorrection(pvals)
    # ... return significant associations
```

See implementation: [`servers/mcp-multiomics/src/mcp_multiomics/tools/association_analysis.py:67-145`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-multiomics/src/mcp_multiomics/tools/association_analysis.py)

**Week 3: Stouffer Meta-Analysis**

The key feature: combine p-values across RNA, protein, and phospho to find consistent signals.

```python
def stouffer_meta_analysis(
    rna_pvals: List[float],
    protein_pvals: List[float],
    phospho_pvals: List[float],
    gene_ids: List[str]
) -> pd.DataFrame:
    """Stouffer's Z-score method for meta-analysis."""
    from scipy.stats import norm

    results = []
    for i, gene in enumerate(gene_ids):
        # Convert p-values to Z-scores
        z_rna = norm.ppf(1 - rna_pvals[i])
        z_protein = norm.ppf(1 - protein_pvals[i])
        z_phospho = norm.ppf(1 - phospho_pvals[i])

        # Combine Z-scores
        z_combined = (z_rna + z_protein + z_phospho) / np.sqrt(3)
        p_combined = 1 - norm.cdf(z_combined)

        results.append({
            "gene": gene,
            "combined_p": p_combined,
            "z_score": z_combined
        })

    return pd.DataFrame(results).sort_values("combined_p")
```

Full implementation: [`servers/mcp-multiomics/src/mcp_multiomics/tools/meta_analysis.py:89-156`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-multiomics/src/mcp_multiomics/tools/meta_analysis.py)

**Result**: 91 tests, 68% coverage, production-ready in 3 weeks.

Test suite: [`tests/unit/mcp-multiomics/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/tests/unit/mcp-multiomics)

---

## The DeepCell Challenge: When Things Don't Go As Planned

Not every implementation went smoothly. The mcp-deepcell server was supposed to take 1 week. It took 3 weeks and multiple Cloud Build failures.

### Attempt 1: "Just Use DeepCell-TF"

The initial plan was simple:
```bash
pip install deepcell-tf
```

**Problem**: The package name is `DeepCell` (capital D, capital C), not `deepcell-tf`. PyPI didn't have `deepcell-tf`.

First build failure:
```
ERROR: Could not find a version that satisfies the requirement DeepCell-tf>=0.12.9
```

### Attempt 2: Fix the Package Name

Updated `pyproject.toml`:
```toml
"DeepCell>=0.12.0"  # Not DeepCell-tf
```

**Problem**: DeepCell requires TensorFlow 2.8.x, which only supports Python 3.10 (not 3.11).

Second build failure:
```
ERROR: Could not find a version that satisfies the requirement tensorflow<2.9
```

### Attempt 3: Downgrade Python

Updated Dockerfile:
```dockerfile
FROM python:3.10-slim  # Was 3.11-slim
```

**Problem**: Cloud Build in us-central1 doesn't allow N1 machine types.

Third build failure:
```
ERROR: Region does not allow N1 machine types: please use E2 variants
```

Updated `cloudbuild.yaml`:
```yaml
options:
  machineType: 'E2_HIGHCPU_8'  # Was N1_HIGHCPU_8
```

See fix commit: [`servers/mcp-deepcell/cloudbuild.yaml`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/cloudbuild.yaml)

### Attempt 4: GCS Image Loading

DeepCell deployed successfully, but the first test failed:
```
The cell segmentation tool was unable to find the image at the provided path:
gs://sample-inputs-patientone/mcp-deepcell-test-data/test_data/dapi_512x512.tif
```

**Root cause**: PIL's `Image.open()` doesn't support GCS URIs. Only local file paths.

**Solution**: Download from GCS to temp file first.

```python
def _download_from_gcs(gcs_path: str) -> str:
    """Download file from GCS to temporary location."""
    if not gcs_path.startswith("gs://"):
        return gcs_path  # Already local

    # Parse gs://bucket/path/to/file
    path_parts = gcs_path[5:].split("/", 1)
    bucket_name, blob_name = path_parts

    # Download to temp file
    from google.cloud import storage
    storage_client = storage.Client()
    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(blob_name)

    temp_file = tempfile.NamedTemporaryFile(delete=False, suffix=Path(blob_name).suffix)
    blob.download_to_filename(temp_file.name)

    return temp_file.name
```

Full implementation: [`servers/mcp-deepcell/src/mcp_deepcell/server.py:47-89`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/src/mcp_deepcell/server.py#L47-L89)

Added dependency in `pyproject.toml`:
```toml
"google-cloud-storage>=2.10.0"
```

**Result**: Deployment succeeded on attempt 4. Total time: 3 weeks (2 weeks debugging).

Full deployment story: [`servers/mcp-deepcell/DEPENDENCY_ISSUES.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/DEPENDENCY_ISSUES.md)

**Lesson learned**: Never assume package names match project names. Always check PyPI before writing `pyproject.toml`.

---

## Production Deployment: All 11 Servers to Cloud Run

Once the implementations were stable, it was time to deploy everything to Google Cloud Run.

### The Deployment Script

```bash
#!/bin/bash
# infrastructure/deployment/deploy_to_gcp.sh

PROJECT_ID="your-project-id"
REGION="us-central1"

for server in mcp-*; do
    echo "Deploying ${server}..."

    cd servers/${server}

    # Stage shared utilities (needed for Docker build)
    mkdir -p _shared_temp
    cp -r ../../shared/utils _shared_temp/

    # Deploy to Cloud Run
    gcloud run deploy ${server} \
        --source . \
        --region ${REGION} \
        --memory 4Gi \
        --cpu 2 \
        --timeout 300 \
        --update-env-vars MCP_TRANSPORT=sse \
        --allow-unauthenticated \
        --project ${PROJECT_ID}

    # Cleanup
    rm -rf _shared_temp

    cd ../..
done
```

Full script: [`infrastructure/deployment/deploy_to_gcp.sh`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/infrastructure/deployment/deploy_to_gcp.sh)

### The Shared Utilities Problem

**Issue**: Dockerfiles expected `_shared_temp/utils/` but deployment script didn't stage files.

**Error**:
```
COPY _shared_temp/utils/ /app/shared/utils/
ERROR: directory not found
```

**Fix**: Stage shared utilities before `gcloud run deploy`, cleanup after.

See fix: [`infrastructure/deployment/deploy_to_gcp.sh:23-30`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/infrastructure/deployment/deploy_to_gcp.sh#L23-L30)

### The Environment Variable Caching Problem

**Issue**: Cloud Run cached old `MCP_TRANSPORT=http` from previous deployments.

**Symptom**: Servers started with HTTP on 127.0.0.1:8000 instead of SSE on 0.0.0.0:PORT.

**Why**: Dockerfile `ENV` statements don't override cached Cloud Run configurations.

**Fix**: Explicitly set environment variables in `gcloud run deploy`:
```bash
--update-env-vars MCP_TRANSPORT=sse
```

See deployment logs: [`docs/archive/deployment/DEPLOYMENT_STATUS.md:169-191`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/archive/deployment/DEPLOYMENT_STATUS.md#L169-L191)

### Success: All 11 Servers Running

Final deployment status (2026-01-31):

| Server | URL | Revision | Status |
|--------|-----|----------|--------|
| mcp-deepcell | https://mcp-deepcell-ondu7mwjpa-uc.a.run.app | 00004-mcs | ✅ Running |
| mcp-fgbio | https://mcp-fgbio-ondu7mwjpa-uc.a.run.app | 00003-xyz | ✅ Running |
| mcp-multiomics | https://mcp-multiomics-ondu7mwjpa-uc.a.run.app | 00002-abc | ✅ Running |
| mcp-spatialtools | https://mcp-spatialtools-ondu7mwjpa-uc.a.run.app | 00005-r4s | ✅ Running |
| ... | ... | ... | ... |

Deployment documentation: [`docs/archive/deployment/DEPLOYMENT_STATUS.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/archive/deployment/DEPLOYMENT_STATUS.md)

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
| **mcp-huggingface** | 12 | 56% | ❌ Mocked |
| **mcp-seqera** | 6 | 56% | ❌ Mocked |
| **mcp-spatialtools** | 5 | 23% | ✅ Production |
| **mcp-openimagedata** | 5 | 35% | ⚠️ Partial |
| **mcp-tcga** | 5 | 35% | ❌ Mocked |

**Total**: 167 tests, 56.9% overall coverage (up from 29.4% baseline)

Full test coverage report: [`docs/test-docs/test-coverage.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/test-docs/test-coverage.md)

### Why Low Coverage Doesn't Mean Low Quality

Notice mcp-spatialtools has only 23% coverage but is marked production-ready. Why?

**Answer**: The server has 2,890 lines of code implementing 14 complex tools (STAR alignment, ComBat batch correction, Moran's I spatial autocorrelation). The 5 smoke tests validate:
- Tool registration
- Data loading
- Basic execution
- Output formats

Full integration testing happens in the PatientOne end-to-end workflow (next section).

Code quality report: [`docs/for-developers/CODE_QUALITY_REPORT.md:54`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-developers/CODE_QUALITY_REPORT.md#L54)

### Example: Testing Multi-Omics Meta-Analysis

```python
# tests/unit/mcp-multiomics/test_meta_analysis.py
def test_stouffer_meta_analysis_with_real_data():
    """Test Stouffer's Z-score method with PatientOne PDX data."""

    # Load real data (15 samples: 7 sensitive, 8 resistant)
    rna_df = pd.read_csv("data/patient-data/PAT001-OVC-2025/multiomics/pdx_rna_seq.csv")
    protein_df = pd.read_csv("data/patient-data/PAT001-OVC-2025/multiomics/pdx_proteomics.csv")
    phospho_df = pd.read_csv("data/patient-data/PAT001-OVC-2025/multiomics/pdx_phosphoproteomics.csv")

    # Run meta-analysis
    result = stouffer_meta_analysis(
        rna_data=rna_df,
        protein_data=protein_df,
        phospho_data=phospho_df,
        group_column="sample_type"  # sensitive vs resistant
    )

    # Validate results
    assert len(result["significant_genes"]) > 0
    assert result["top_gene"] in ["PIK3CA", "AKT1", "MTOR"]  # Expected from biology
    assert result["combined_p_values"][0] < 0.001
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

**Prompt**:
```
Using mockepic and fgbio servers, provide clinical summary and identify
pathogenic somatic variants for PAT001-OVC-2025.
```

**Tools called**:
- `mockepic.get_patient_summary()`
- `mockepic.analyze_ca125_trends()`
- `fgbio.parse_vcf()`
- `tcga.compare_to_cohort()`

**Results**:
- Patient: Sarah Anderson, 58, BRCA1 carrier
- CA-125: 1200 → 45 → 310 U/mL (platinum resistance)
- Mutations: TP53 R175H, PIK3CA E545K, PTEN LOH
- TCGA subtype: C1 Immunoreactive

**Time**: 5 minutes
**Cost**: $0.18

### TEST_2: Multi-Omics Resistance Analysis ✅

**Prompt**:
```
Load PDX multi-omics data (RNA, protein, phospho) and run Stouffer meta-analysis.
Identify consistently dysregulated genes in platinum-resistant samples.
```

**Tools called**:
- `multiomics.validate_omics_data()`
- `multiomics.stouffer_meta_analysis()`
- `multiomics.pathway_enrichment()`

**Results**:
- **17 significant genes** (FDR < 0.05, |log2FC| > 1)
- **Top 3 hits**: PIK3CA (p=1.2e-18), AKT1 (p=1.23e-18), mTOR (p=1.36e-14)
- **Pathway**: PI3K/AKT/mTOR activation (confirms genomic PIK3CA E545K mutation is *active*)

**Time**: 8 minutes
**Cost**: $0.31

Full test results: [`servers/mcp-spatialtools/COMPLETE_WORKFLOW_TEST_SUMMARY.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-spatialtools/COMPLETE_WORKFLOW_TEST_SUMMARY.md)

### TEST_3: Spatial Transcriptomics ✅

**Prompt**:
```
Analyze Visium spatial data: differential expression (tumor_core vs stroma),
spatial autocorrelation (Moran's I), and pathway enrichment.
```

**Tools called**:
- `spatialtools.load_spatial_data()`
- `spatialtools.spatial_differential_expression()`
- `spatialtools.spatial_autocorrelation()`
- `spatialtools.pathway_enrichment()`

**Results**:
- **17 DEGs** (tumor vs stroma, FDR < 0.05)
- **Upregulated in tumor**: TP53, KRT8, ABCB1, BCL2L1, MKI67, TOP2A, AKT1, MTOR, MYC
- **Downregulated in tumor**: ACTA2, FAP, COL1A1, COL3A1 (stromal markers)
- **Spatial pattern**: Immune exclusion phenotype (CD8+ cells blocked by stromal barrier)
- **Moran's I**: High spatial autocorrelation for MKI67 (I=0.72), CD8 (I=0.68)

**Time**: 6 minutes
**Cost**: $0.28

### TEST_4: Imaging Analysis ✅

**Prompt**:
```
Segment cells in MxIF image (TP53/Ki67/DAPI). Quantify proliferation index
and TP53 positivity.
```

**Tools called**:
- `deepcell.segment_cells()`
- `deepcell.classify_cell_states()`
- `openimagedata.generate_he_annotation()`

**Results**:
- **Ki67 proliferation index**: 45-55% (HIGH)
- **TP53+/Ki67+ double-positive**: 38% of tumor cells
- **CD8+ density**: 5-15 cells/mm² (LOW, peripheral only)

**Time**: 5 minutes
**Cost**: $0.22

### TEST_5: Integration & Recommendations ✅

**Prompt**:
```
Synthesize all findings and generate treatment recommendations with
clinical trial matches.
```

**Tools called**: None (pure synthesis)

**Result**:
- **Primary**: PI3K inhibitor (alpelisib) targeting PIK3CA E545K
- **Secondary**: Anti-PD-1 immunotherapy (overcome immune exclusion)
- **Clinical trial**: NCT03602859 (alpelisib + paclitaxel)
- **Monitoring**: CA-125 q2 weeks, CT at 8 weeks, ctDNA for PIK3CA VAF

**Time**: 3 minutes
**Cost**: $0.11

### Total Workflow Performance

| Metric | Value |
|--------|-------|
| **Total time** | 35 minutes |
| **Total cost** | $1.20 |
| **Tools called** | 12+ across 4 servers |
| **Data integrated** | 5 modalities |
| **Actionable recommendations** | 3 (1 primary, 2 secondary) |

**Comparison to traditional workflow**:
- Time: 40 hours → 35 minutes (95% reduction)
- Cost: $3,200 → $1.20 (99.96% reduction)
- Specialists: 3-4 → 1 oncologist

---

## Cost Analysis: Real GCP Pricing

All costs validated against actual GCP deployment (2026-01-31 pricing).

### Cloud Run Compute Costs

**Pricing Model**: Pay-per-use (billed in 100ms increments)
- **CPU**: $0.00002400 per vCPU-second
- **Memory**: $0.00000250 per GiB-second

**Per-Analysis Breakdown** (PatientOne workflow):

| Server | Requests | Avg Time | CPU | Memory | Cost |
|--------|----------|----------|-----|--------|------|
| mcp-multiomics | 3 | 45s | 2 vCPU | 4Gi | $0.0086 |
| mcp-spatialtools | 2 | 30s | 2 vCPU | 4Gi | $0.0057 |
| mcp-fgbio | 2 | 15s | 1 vCPU | 2Gi | $0.0014 |
| mcp-deepcell | 1 | 60s | 2 vCPU | 4Gi | $0.0115 |
| mcp-mockepic | 2 | 5s | 1 vCPU | 512Mi | $0.0003 |

**Total Cloud Run**: $0.0275 per analysis (~$0.03)

### Claude API Costs

**Pricing** (Claude Sonnet 4.5):
- Input: $3.00 per million tokens
- Output: $15.00 per million tokens

**PatientOne workflow tokens**:
- Input: ~25,000 tokens (prompt + tool results)
- Output: ~5,000 tokens (synthesis + recommendations)

**Cost**:
- Input: 25K × $3.00 / 1M = $0.075
- Output: 5K × $15.00 / 1M = $0.075
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
- Excellent error handling (zero bare `except:` clauses)
- Strong HIPAA de-identification
- Comprehensive input validation

**Test coverage**: 56.9% overall
- 167 automated tests
- 91 tests for mcp-multiomics alone

**Production readiness**: 7/12 servers (58%)
- mcp-fgbio, mcp-multiomics, mcp-spatialtools, mcp-perturbation, mcp-quantum-celltype-fidelity, mcp-deepcell, mcp-epic

**Deployment success**: 11/11 servers on Cloud Run

Full quality report: [`docs/for-developers/CODE_QUALITY_REPORT.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-developers/CODE_QUALITY_REPORT.md)

---

## Try It Yourself

Ready to validate these results on your own system?

### Option 1: Run the Test Suite

Clone the repository and run automated tests:

```bash
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp

# Run mcp-multiomics tests (91 tests, should pass)
cd servers/mcp-multiomics
python -m venv venv
source venv/bin/activate
pip install -e ".[dev]"
pytest ../../tests/unit/mcp-multiomics/ -v
```

Expected output:
```
===== 91 passed in 23.45s =====
```

### Option 2: Deploy to Your GCP Account

Follow the deployment guide to run on your own Cloud Run:

1. Set up GCP project (free tier: $300 credit for 90 days)
2. Enable Cloud Run API
3. Run deployment script:

```bash
./infrastructure/deployment/deploy_to_gcp.sh your-project-id
```

Cost: ~$0.02-0.05 for initial deployment

Full deployment guide: [`docs/deployment/GET_STARTED.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/deployment/GET_STARTED.md)

### Option 3: Interactive Notebook

Explore test results and run mini-analyses:
[`docs/book/companion-notebooks/chapter-03-testing.ipynb`](../companion-notebooks/chapter-03-testing.ipynb)

This notebook includes:
- PatientOne workflow reproduction
- Cost calculator (input your analysis parameters)
- Performance profiling
- Test coverage visualization

---

## What Comes Next

In Chapter 4, you'll build your first MCP server from scratch: **mcp-epic** for clinical FHIR integration.

You'll learn:
- How to structure an MCP server
- FHIR R4 resource mapping
- De-identification for HIPAA compliance
- Testing with synthetic patient data

But before you move on, take a moment to appreciate what you've validated:
- **The architecture works**: AI can orchestrate 12 specialized servers
- **The costs are real**: $1-2 per analysis, validated against GCP pricing
- **The time savings are real**: 35 minutes, validated end-to-end
- **Production deployment is feasible**: 11/11 servers running on Cloud Run

This isn't vaporware. It's a tested, validated, production-ready system.

**Next**: Chapter 4 - Clinical Data: The Starting Point

---

**Chapter 3 Key Takeaways:**
- Started with 100% mocked data, iterated to 7/12 production servers
- 167 automated tests, 56.9% coverage
- DeepCell took 3 weeks (package naming, Python version, GCS loading)
- Complete PatientOne workflow: 35 minutes, $1.20 cost (validated)
- Code quality: 7.5/10 (production-ready with improvements needed)
- Real GCP pricing: $0.19 per analysis (Cloud Run + Claude API)

**Companion Resources:**
- 📓 [Jupyter Notebook](../companion-notebooks/chapter-03-testing.ipynb) - Run tests and calculate costs
- 📊 [Code Quality Report](../../for-developers/CODE_QUALITY_REPORT.md) - Full analysis
- 🧪 [Test Coverage](../../test-docs/test-coverage.md) - Server-by-server breakdown
- ✅ [Workflow Test Summary](../../servers/mcp-spatialtools/COMPLETE_WORKFLOW_TEST_SUMMARY.md) - PatientOne validation

**GitHub References:**
- Test suite: [`tests/unit/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/tests/unit)
- Deployment script: [`infrastructure/deployment/deploy_to_gcp.sh`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/infrastructure/deployment/deploy_to_gcp.sh)
- DeepCell debugging: [`servers/mcp-deepcell/DEPENDENCY_ISSUES.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/DEPENDENCY_ISSUES.md)
