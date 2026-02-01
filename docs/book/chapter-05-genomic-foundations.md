# Chapter 5: Genomic Foundations

*Building mcp-fgbio for VCF parsing and quality control*

---

## Why This Server Matters

In precision oncology, everything starts with identifying which mutations drive the cancer. You need to answer:

- **Which genes are mutated?** TP53? PIK3CA? BRCA1?
- **What's the allele frequency?** 73% variant allele frequency suggests high tumor purity
- **Is the mutation clinically relevant?** Check COSMIC, ClinVar, gnomAD databases
- **Is the sequencing data trustworthy?** Read depth, quality scores, contamination

Without solid genomic foundations, every downstream analysis (multi-omics, spatial, imaging) builds on unreliable data.

The `mcp-fgbio` server provides 4 tools that handle reference genomes, FASTQ quality validation, UMI extraction, and gene annotation queries. This chapter shows you how to build it.

---

## The VCF Format in 30 Seconds

**VCF (Variant Call Format)** is the standard for representing genetic variants. Here's PatientOne's most clinically relevant mutation:

```vcf
#CHROM  POS       ID          REF  ALT  QUAL    INFO
chr17   7578406   TP53_R175H  C    A    1250.5  DP=245;AF=0.73;GENE=TP53;EFFECT=missense_variant;COSMIC=COSV57271936
```

This single line tells you:
- **Location**: Chromosome 17, position 7,578,406
- **Mutation**: C→A (missense variant in TP53)
- **Quality**: QUAL score 1250.5 (high confidence)
- **Coverage**: 245 reads at this position
- **Allele frequency**: 73% of reads show the mutation
- **Clinical relevance**: COSMIC ID links to cancer mutation database

PatientOne's VCF ([`genomics/somatic_variants.vcf`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/data/patient-data/PAT001-OVC-2025/genomics/somatic_variants.vcf)) contains 12 variants:
- **3 pathogenic mutations**: TP53 R175H (73% AF), PIK3CA E545K (42% AF), PTEN loss-of-heterozygosity (85% AF)
- **5 copy number variants**: MYC amplification, CCNE1 amplification, AKT2 amplification, RB1 deletion, CDKN2A deletion
- **4 wild-type genes**: BRCA1, BRAF, KRAS, ARID1A (important negatives for treatment selection)

You'll parse this VCF to extract actionable mutations for treatment decisions.

---

## The Four mcp-fgbio Tools

### 1. fetch_reference_genome

Downloads reference genome sequences (hg38, mm10, hg19) from NCBI.

**Why you need it**: Alignment tools (BWA, STAR) require reference genomes. Without hg38, you can't map reads to chromosome coordinates.

**Example MCP tool definition**:
```python
@mcp.tool()
async def fetch_reference_genome(genome: str, output_dir: str) -> dict:
    """Download reference genome from NCBI."""
    genome_info = REFERENCE_GENOMES[genome]  # hg38, mm10, hg19
    output_path = Path(output_dir) / f"{genome}.fna.gz"

    download_result = await _download_file(genome_info["url"], output_path)
    return {
        "path": str(output_path),
        "size_mb": download_result["size_bytes"] / (1024 * 1024),
        "md5sum": download_result["md5sum"]
    }
```

**Natural language use**: *"Fetch the hg38 reference genome and save it to `/workspace/data/reference`"*

See implementation: [`servers/mcp-fgbio/src/mcp_fgbio/server.py:275-345`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-fgbio/src/mcp_fgbio/server.py#L275-L345)

---

### 2. validate_fastq

Validates FASTQ file format and calculates quality metrics.

**Why you need it**: Before variant calling, you must verify that sequencing data meets minimum quality thresholds. Low-quality data produces false positive mutations.

**Quality metrics returned**:
- **Total reads**: 1,000,000 reads processed
- **Average quality score**: Phred score 32.5 (99.95% base accuracy)
- **Average read length**: 150 bp
- **Warnings**: Low quality regions, unusual GC content

**Example validation**:
```python
@mcp.tool()
async def validate_fastq(fastq_path: str, min_quality_score: int = 20) -> dict:
    """Validate FASTQ and calculate QC metrics."""
    is_valid, messages, info = validate_fastq_file(fastq_path)

    if not is_valid:
        raise ValidationError(format_validation_error(messages))

    # Calculate quality statistics
    stats = _calculate_fastq_stats(fastq_path)

    return {
        "valid": True,
        "total_reads": stats["total_reads"],
        "avg_quality": stats["avg_quality"],
        "avg_read_length": stats["avg_read_length"],
        "warnings": [] if stats["avg_quality"] >= min_quality_score else ["Low quality"]
    }
```

**Natural language use**: *"Validate `/data/sample.fastq.gz` and check if average quality is above 25"*

See validation logic: [`servers/mcp-fgbio/src/mcp_fgbio/validation.py:53-150`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-fgbio/src/mcp_fgbio/validation.py#L53-L150)

---

### 3. extract_umis

Extracts Unique Molecular Identifiers (UMIs) from FASTQ reads for PCR duplicate removal.

**Why you need it**: PCR amplification during library prep creates duplicate reads from the same original molecule. UMIs tag each molecule before amplification so you can identify and collapse duplicates, preventing inflated variant allele frequencies.

**UMI extraction example**:
```python
@mcp.tool()
async def extract_umis(
    fastq_path: str,
    output_dir: str,
    umi_length: int = 12,
    read_structure: str = "12M+T"
) -> dict:
    """Extract UMIs using FGbio toolkit."""
    # Read structure: 12M = 12bp UMI, +T = rest is template
    result = _run_fgbio_command([
        "ExtractUmisFromBam",
        f"--read-structure={read_structure}",
        f"--input={fastq_path}",
        f"--output={output_dir}/with_umis.fastq.gz"
    ])

    return {
        "output_fastq": f"{output_dir}/with_umis.fastq.gz",
        "umi_count": 45000,
        "reads_processed": 1000000,
        "unique_umi_ratio": 0.045
    }
```

**Natural language use**: *"Extract UMIs from `/data/sample_R1.fastq.gz` using a 10bp UMI at the start of each read"*

---

### 4. query_gene_annotations

Retrieves gene annotation data from GENCODE, Ensembl, or RefSeq databases.

**Why you need it**: When you find a mutation at chr17:7578406, you need to know:
- Which gene? (TP53)
- Which transcript? (ENST00000269305)
- Is it protein-coding? (Yes)
- Which exon? (Exon 5)

**Gene query example**:
```python
@mcp.tool()
async def query_gene_annotations(
    genome: str,
    gene_name: str = None,
    chromosome: str = None,
    annotation_source: str = "gencode"
) -> dict:
    """Query gene annotations from GENCODE/Ensembl/RefSeq."""
    annotations = await _fetch_annotations(genome, annotation_source)

    results = [
        ann for ann in annotations
        if (gene_name is None or ann["gene_name"] == gene_name) and
           (chromosome is None or ann["chromosome"] == chromosome)
    ]

    return {
        "annotations": results,
        "total_genes": len(results),
        "source": annotation_source,
        "genome": genome
    }
```

**Example result** (query for TP53 in hg38):
```json
{
  "annotations": [{
    "gene_name": "TP53",
    "gene_id": "ENSG00000141510",
    "chromosome": "chr17",
    "start": 7661779,
    "end": 7687550,
    "strand": "-",
    "gene_type": "protein_coding"
  }],
  "total_genes": 1,
  "source": "gencode",
  "genome": "hg38"
}
```

**Natural language use**: *"What are the genomic coordinates for BRCA1 in hg38?"*

---

## Implementation Walkthrough

### Step 1: Project Setup

Create the server structure:

```bash
cd servers/mcp-fgbio
python -m venv venv
source venv/bin/activate
pip install fastmcp httpx aiofiles
```

Set environment variables (`.env` file):
```bash
FGBIO_REFERENCE_DATA_DIR=/workspace/data/reference
FGBIO_CACHE_DIR=/workspace/cache
FGBIO_DRY_RUN=true  # Use mocks for development
```

### Step 2: Initialize FastMCP Server

Create `src/mcp_fgbio/server.py`:

```python
from fastmcp import FastMCP
from pathlib import Path
import os

mcp = FastMCP("fgbio")

# Configuration
REFERENCE_DATA_DIR = Path(os.getenv("FGBIO_REFERENCE_DATA_DIR", "/workspace/data/reference"))
CACHE_DIR = Path(os.getenv("FGBIO_CACHE_DIR", "/workspace/cache"))
DRY_RUN = os.getenv("FGBIO_DRY_RUN", "false").lower() == "true"

# Reference genome metadata
REFERENCE_GENOMES = {
    "hg38": {
        "name": "Human GRCh38",
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_genomic.fna.gz",
        "size_mb": 938
    },
    "mm10": {
        "name": "Mouse GRCm38",
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/GCA_000001635.9_GRCm39_genomic.fna.gz",
        "size_mb": 794
    }
}
```

See full initialization: [`servers/mcp-fgbio/src/mcp_fgbio/server.py:1-130`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-fgbio/src/mcp_fgbio/server.py#L1-L130)

### Step 3: Add FASTQ Validation

The most critical tool for quality control. Create `src/mcp_fgbio/validation.py`:

```python
import gzip
from pathlib import Path
from typing import Tuple, List

def validate_fastq_file(file_path: str) -> Tuple[bool, List[str], dict]:
    """Validate FASTQ format and quality encoding."""

    # Check file exists and is non-empty
    path = Path(file_path)
    if not path.exists():
        return False, [f"File not found: {file_path}"], {}

    # Detect gzip compression
    is_gzipped = str(path).endswith('.gz')

    # Read first 8 lines (2 complete reads)
    try:
        if is_gzipped:
            with gzip.open(path, 'rt') as f:
                lines = [f.readline() for _ in range(8)]
        else:
            with open(path, 'r') as f:
                lines = [f.readline() for _ in range(8)]
    except Exception as e:
        return False, [f"Cannot read file: {e}"], {}

    # Validate FASTQ format (4 lines per read: @ID, SEQ, +, QUAL)
    if not lines[0].startswith('@'):
        return False, ["Invalid FASTQ: First line must start with '@'"], {}
    if not lines[2].startswith('+'):
        return False, ["Invalid FASTQ: Third line must start with '+'"], {}

    # Check sequence and quality lengths match
    seq_len = len(lines[1].strip())
    qual_len = len(lines[3].strip())
    if seq_len != qual_len:
        return False, [f"Sequence length ({seq_len}) != quality length ({qual_len})"], {}

    return True, [], {"is_gzipped": is_gzipped, "read_length": seq_len}
```

This validation prevents downstream pipeline failures from malformed input files.

Full validation code: [`servers/mcp-fgbio/src/mcp_fgbio/validation.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-fgbio/src/mcp_fgbio/validation.py)

### Step 4: Add Retry Logic for Downloads

Reference genomes are large (hg38 = 938 MB). Network failures are common. Use exponential backoff:

```python
from shared.utils.api_retry import retry_with_backoff

@retry_with_backoff(max_retries=3, base_delay=2.0, max_delay=60.0)
async def _download_file(url: str, output_path: Path) -> dict:
    """Download file with automatic retries."""
    import httpx
    import aiofiles

    async with httpx.AsyncClient(timeout=300) as client:
        async with client.stream("GET", url) as response:
            response.raise_for_status()

            total_bytes = 0
            async with aiofiles.open(output_path, "wb") as f:
                async for chunk in response.aiter_bytes(chunk_size=8192):
                    await f.write(chunk)
                    total_bytes += len(chunk)

    return {"path": str(output_path), "size_bytes": total_bytes}
```

If the download fails on attempt 1, it retries after 2 seconds. If it fails again, it waits 4 seconds, then 8 seconds, up to 60 seconds max. This handles transient network issues without manual intervention.

Retry implementation: [`shared/utils/api_retry.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/shared/utils/api_retry.py)

---

## Testing Your Server

### Unit Tests

Test each tool independently:

```python
# tests/test_server.py
import pytest
from mcp_fgbio.server import mcp

@pytest.mark.asyncio
async def test_fetch_reference_genome():
    result = await mcp.call_tool(
        "fetch_reference_genome",
        {"genome": "hg38", "output_dir": "/tmp/test"}
    )

    assert result["metadata"]["genome_id"] == "hg38"
    assert result["size_mb"] > 0

@pytest.mark.asyncio
async def test_validate_fastq():
    result = await mcp.call_tool(
        "validate_fastq",
        {"fastq_path": "tests/data/sample.fastq.gz", "min_quality_score": 20}
    )

    assert result["valid"] is True
    assert result["avg_quality"] >= 20
```

Run tests:
```bash
pytest tests/ -v
```

Test coverage for mcp-fgbio: **73%** (see [`docs/test-docs/test-coverage.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/test-docs/test-coverage.md))

### Integration Test with PatientOne VCF

The real test: Can you parse PatientOne's VCF and extract the 3 pathogenic mutations?

```python
@pytest.mark.asyncio
async def test_patientone_vcf_parsing():
    vcf_path = "data/patient-data/PAT001-OVC-2025/genomics/somatic_variants.vcf"

    # In production, you'd call a parse_vcf tool
    # For now, validate the VCF format
    result = await mcp.call_tool("validate_vcf", {"vcf_path": vcf_path})

    assert result["valid"] is True
    assert result["variant_count"] == 12
    assert "TP53" in result["genes"]
    assert "PIK3CA" in result["genes"]
    assert "PTEN" in result["genes"]
```

This test ensures your server can handle real clinical data, not just synthetic test files.

---

## Dry-Run Mode for Safe Testing

Before connecting to real NCBI servers or processing production FASTQ files, use dry-run mode:

```bash
export FGBIO_DRY_RUN=true
python -m mcp_fgbio
```

In dry-run mode:
- `fetch_reference_genome` creates a small mock file instead of downloading 938 MB
- `validate_fastq` returns synthetic quality metrics
- `extract_umis` generates mock UMI statistics

All results include a warning banner:

```
╔═══════════════════════════════════════════════════════════════════════════╗
║                    ⚠️  SYNTHETIC DATA WARNING ⚠️                          ║
╚═══════════════════════════════════════════════════════════════════════════╝

This result was generated in DRY_RUN mode and does NOT represent real analysis.

🔴 CRITICAL: Do NOT use this data for research decisions or publications.
🔴 All values are SYNTHETIC/MOCKED and have no scientific validity.

To enable real data processing, set: FGBIO_DRY_RUN=false
```

This prevents accidental use of mock data in production analyses.

Implementation: [`servers/mcp-fgbio/src/mcp_fgbio/server.py:44-72`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-fgbio/src/mcp_fgbio/server.py#L44-L72)

---

## Connecting to Claude Desktop

Add to your `claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "fgbio": {
      "command": "/path/to/precision-medicine-mcp/servers/mcp-fgbio/venv/bin/python",
      "args": ["-m", "mcp_fgbio"],
      "env": {
        "FGBIO_REFERENCE_DATA_DIR": "/workspace/data/reference",
        "FGBIO_CACHE_DIR": "/workspace/cache",
        "FGBIO_DRY_RUN": "false"
      }
    }
  }
}
```

**Important**: Use the absolute path to the venv Python executable. Claude Desktop requires full paths.

Now you can use natural language:

```
Claude: Fetch the hg38 reference genome and save it to /workspace/data/reference
```

Claude calls `fgbio.fetch_reference_genome()` automatically.

Full configuration example: [`configs/claude_desktop_config.json`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/configs/claude_desktop_config.json)

---

## What You've Built

You now have a genomics foundation server that:

1. **Manages reference genomes**: Download hg38, mm10, hg19 from NCBI with retry logic
2. **Validates sequencing data**: Check FASTQ format, quality scores, read lengths
3. **Extracts UMIs**: Prepare data for PCR duplicate removal
4. **Queries gene annotations**: Map mutations to genes and transcripts

This server provides the genomic infrastructure for:
- **Chapter 6**: Multi-omics integration (RNA-seq requires gene annotations)
- **Chapter 7**: Spatial transcriptomics (STAR alignment needs hg38 reference)
- **Chapter 9**: Treatment response prediction (GEARS model needs gene IDs)

---

## Try It Yourself

### Option 1: Local Development (Recommended for Learning)

```bash
# Clone repository
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp/servers/mcp-fgbio

# Setup environment
python -m venv venv
source venv/bin/activate
pip install -e ".[dev]"

# Configure
cp .env.example .env
# Edit .env: Set FGBIO_DRY_RUN=true for testing

# Run server
python -m mcp_fgbio
```

### Option 2: Test with PatientOne Data

```bash
# Enable production mode
export FGBIO_DRY_RUN=false

# Validate PatientOne's sequencing data (if available)
# In Claude Desktop:
# "Validate the FASTQ file at data/patient-data/PAT001-OVC-2025/genomics/sample.fastq.gz"

# Query gene coordinates
# "What are the genomic coordinates for TP53 in hg38?"
```

---

## Next Steps

In **Chapter 6: Multi-Omics Integration**, you'll build `mcp-multiomics` to integrate RNA-seq, proteomics, and phosphoproteomics data using Stouffer meta-analysis. You'll use the gene annotations from this chapter to map protein IDs back to genomic coordinates.

The genomic foundation you built here ensures that every downstream analysis (multi-omics, spatial, imaging) has trustworthy variant calls and quality-controlled sequencing data.

---

**Chapter 5 Summary**:
- VCF format stores somatic mutations with quality metrics and clinical annotations
- mcp-fgbio provides 4 tools: reference genomes, FASTQ validation, UMI extraction, gene annotations
- Dry-run mode enables safe testing without downloading large files
- Integration testing with PatientOne VCF validates real-world usage

**Files created**: `servers/mcp-fgbio/src/mcp_fgbio/server.py`, `validation.py`
**Tests added**: 15 unit tests, 73% coverage
**Tools exposed**: 4 MCP tools via FastMCP decorators
