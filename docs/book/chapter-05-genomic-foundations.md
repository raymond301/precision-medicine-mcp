# Genomic Foundations

*Building mcp-fgbio for VCF parsing and quality control*

---

## Why This Server Matters

Precision oncology requires identifying which mutations drive cancer:
- **Which genes are mutated?** TP53? PIK3CA? BRCA1?
- **What's the allele frequency?** 73% VAF suggests high tumor purity
- **Is it clinically relevant?** Check COSMIC, ClinVar, gnomAD
- **Is sequencing trustworthy?** Read depth, quality scores, contamination

The `mcp-fgbio` server provides 4 tools for reference genomes, FASTQ validation, UMI extraction, and gene annotation queries.

---

## The VCF Format in 30 Seconds

**VCF (Variant Call Format)** is the standard for representing genetic variants. PatientOne's most clinically relevant mutation:

```vcf
#CHROM  POS       ID          REF  ALT  QUAL    INFO
chr17   7578406   TP53_R175H  C    A    1250.5  DP=245;AF=0.73;GENE=TP53;EFFECT=missense_variant
```

This tells you: chr17 position 7,578,406 C→A (TP53 missense), QUAL 1250.5, 245 reads coverage, 73% allele frequency, COSMIC ID.

PatientOne's VCF ([`genomics/somatic_variants.vcf`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/data/patient-data/PAT001-OVC-2025/genomics/somatic_variants.vcf)) contains:
- **3 pathogenic**: TP53 R175H (73% AF), PIK3CA E545K (42% AF), PTEN LOH (85% AF)
- **5 copy number variants**: MYC amplification, CCNE1 amplification, etc.
- **4 wild-type**: BRCA1, BRAF, KRAS, ARID1A (important negatives)

![VCF Parsed Output](vcf-parsed-output.png){width=100%}

**Figure 5.1: Parsed VCF with Annotations**
*Claude output showing parsed and annotated PatientOne VCF file using `fgbio_parse_vcf` and `fgbio_annotate_variants` tools. Table displays gene names, variants (protein changes), mutation types, variant allele frequencies (VAF), and ClinVar pathogenicity ratings with star ratings. Pathogenic mutations (TP53, PIK3CA, PTEN, KRAS) are highlighted for clinical decision-making.*

---

## The Four mcp-fgbio Tools

### 1. fetch_reference_genome

Downloads reference genomes (hg38, mm10, hg19) from NCBI.

```python
@mcp.tool()
async def fetch_reference_genome(genome: str, output_dir: str) -> dict:
    """Download reference genome from NCBI."""
    genome_info = REFERENCE_GENOMES[genome]
    output_path = Path(output_dir) / f"{genome}.fna.gz"
    download_result = await _download_file(genome_info["url"], output_path)
    return {"path": str(output_path), "size_mb": download_result["size_bytes"] / (1024 * 1024)}
    # Full implementation: servers/mcp-fgbio/src/mcp_fgbio/server.py:275-345
```

**Natural language use**: *"Fetch the hg38 reference genome and save to `/workspace/data/reference`"*

---

### 2. validate_fastq

Validates FASTQ format and calculates quality metrics.

```python
@mcp.tool()
async def validate_fastq(fastq_path: str, min_quality_score: int = 20) -> dict:
    """Validate FASTQ and calculate QC metrics."""
    is_valid, messages, info = validate_fastq_file(fastq_path)
    if not is_valid: raise ValidationError(messages)
    stats = _calculate_fastq_stats(fastq_path)
    return {"valid": True, "total_reads": stats["total_reads"], "avg_quality": stats["avg_quality"]}
    # Full implementation: servers/mcp-fgbio/src/mcp_fgbio/validation.py:53-150
```

**Metrics returned**: Total reads, average Phred quality score, average read length, warnings.

---

### 3. extract_umis

Extracts Unique Molecular Identifiers (UMIs) for PCR duplicate removal.

```python
@mcp.tool()
async def extract_umis(fastq_path: str, output_dir: str, umi_length: int = 12, read_structure: str = "12M+T") -> dict:
    """Extract UMIs using FGbio toolkit."""
    # Read structure: 12M = 12bp UMI, +T = rest is template
    result = _run_fgbio_command(["ExtractUmisFromBam", f"--read-structure={read_structure}", ...])
    return {"output_fastq": f"{output_dir}/with_umis.fastq.gz", "umi_count": 45000}
    # Full implementation: servers/mcp-fgbio/src/mcp_fgbio/server.py:350-425
```

**Why you need it**: PCR duplicates inflate variant allele frequencies. UMIs tag molecules before amplification for accurate duplicate removal.

---

### 4. query_gene_annotations

Retrieves gene annotation data from GENCODE, Ensembl, or RefSeq.

```python
@mcp.tool()
async def query_gene_annotations(genome: str, gene_name: str = None, chromosome: str = None, annotation_source: str = "gencode") -> dict:
    """Query gene annotations from GENCODE/Ensembl/RefSeq."""
    annotations = await _fetch_annotations(genome, annotation_source)
    results = [ann for ann in annotations if (gene_name is None or ann["gene_name"] == gene_name)]
    return {"annotations": results, "total_genes": len(results), "source": annotation_source}
    # Full implementation: servers/mcp-fgbio/src/mcp_fgbio/server.py:430-510
```

**Example result** (TP53 in hg38):
```json
{
  "annotations": [{
    "gene_name": "TP53",
    "gene_id": "ENSG00000141510",
    "chromosome": "chr17",
    "start": 7661779,
    "end": 7687550,
    "gene_type": "protein_coding"
  }]
}
```

---

## Implementation Walkthrough

### Project Setup

```bash
cd servers/mcp-fgbio
python -m venv venv && source venv/bin/activate
pip install fastmcp httpx aiofiles
```

Environment variables (`.env`):
```bash
FGBIO_REFERENCE_DATA_DIR=/workspace/data/reference
FGBIO_CACHE_DIR=/workspace/cache
FGBIO_DRY_RUN=true  # Use mocks for development
```

### Initialize FastMCP Server

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
    "hg38": {"name": "Human GRCh38", "url": "https://ftp.ncbi.nlm.nih.gov/genomes/.../GRCh38_genomic.fna.gz", "size_mb": 938},
    "mm10": {"name": "Mouse GRCm38", "url": "https://ftp.ncbi.nlm.nih.gov/genomes/.../GRCm39_genomic.fna.gz", "size_mb": 794}
}
```

Full initialization: [`servers/mcp-fgbio/src/mcp_fgbio/server.py:1-130`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-fgbio/src/mcp_fgbio/server.py#L1-L130)

### Add FASTQ Validation

```python
def validate_fastq_file(file_path: str) -> Tuple[bool, List[str], dict]:
    """Validate FASTQ format and quality encoding."""
    path = Path(file_path)
    if not path.exists(): return False, [f"File not found"], {}
    is_gzipped = str(path).endswith('.gz')
    # Read first 8 lines, validate format, check sequence/quality lengths match
    # Full implementation: servers/mcp-fgbio/src/mcp_fgbio/validation.py
```

### Add Retry Logic for Downloads

```python
from shared.utils.api_retry import retry_with_backoff

@retry_with_backoff(max_retries=3, base_delay=2.0, max_delay=60.0)
async def _download_file(url: str, output_path: Path) -> dict:
    """Download file with automatic retries."""
    # Exponential backoff: 2s, 4s, 8s, up to 60s max
    # Full implementation: shared/utils/api_retry.py
```

---

## Testing Your Server

### Unit Tests

```python
@pytest.mark.asyncio
async def test_fetch_reference_genome():
    result = await mcp.call_tool("fetch_reference_genome", {"genome": "hg38", "output_dir": "/tmp/test"})
    assert result["metadata"]["genome_id"] == "hg38"
    assert result["size_mb"] > 0
```

Test coverage: **73%** ([`docs/test-docs/test-coverage.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/test-docs/test-coverage.md))

---

## Dry-Run Mode

```bash
export FGBIO_DRY_RUN=true
python -m mcp_fgbio
```

In dry-run mode:
- `fetch_reference_genome` creates small mock file instead of downloading 938 MB
- `validate_fastq` returns synthetic quality metrics
- Results include warning banner about synthetic data

Implementation: [`servers/mcp-fgbio/src/mcp_fgbio/server.py:44-72`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-fgbio/src/mcp_fgbio/server.py#L44-L72)

---

## Connecting to Claude Desktop

```json
{
  "mcpServers": {
    "fgbio": {
      "command": "/path/to/venv/bin/python",
      "args": ["-m", "mcp_fgbio"],
      "env": {
        "FGBIO_REFERENCE_DATA_DIR": "/workspace/data/reference",
        "FGBIO_DRY_RUN": "false"
      }
    }
  }
}
```

Full configuration: [`configs/claude_desktop_config.json`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/configs/claude_desktop_config.json)

---

## What You've Built

A genomics foundation server providing:
1. **Reference genome management**: Download hg38, mm10, hg19 from NCBI with retry logic
2. **Sequencing data validation**: Check FASTQ format, quality scores, read lengths
3. **UMI extraction**: Prepare data for PCR duplicate removal
4. **Gene annotations**: Map mutations to genes and transcripts

This supports:
- **Chapter 6**: Multi-omics integration (RNA-seq requires gene annotations)
- **Chapter 7**: Spatial transcriptomics (STAR alignment needs hg38)
- **Chapter 9**: Treatment response prediction (GEARS needs gene IDs)

---

## Try It Yourself

### Local Development

```bash
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp/servers/mcp-fgbio
python -m venv venv && source venv/bin/activate
pip install -e ".[dev]"
cp .env.example .env  # Set FGBIO_DRY_RUN=true for testing
python -m mcp_fgbio
```

---

## Summary

**Chapter 5 Summary**:
- VCF format stores somatic mutations with quality metrics and clinical annotations
- mcp-fgbio provides 4 tools: reference genomes, FASTQ validation, UMI extraction, gene annotations
- Dry-run mode enables safe testing without downloading large files
- Test coverage: 73% with 15 unit tests

**Files created**: `servers/mcp-fgbio/src/mcp_fgbio/server.py`, `validation.py`
**Tools exposed**: 4 MCP tools via FastMCP decorators
