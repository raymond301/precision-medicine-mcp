# Seqera Official MCP Server Integration

This document describes how to integrate **Seqera's official MCP server** (`@seqeralabs/mcp-server-seqera`) with the Precision Medicine MCP platform for upstream pipeline orchestration.

## Overview

The Precision Medicine MCP platform uses Seqera's official MCP server for pipeline orchestration, launching and monitoring Nextflow workflows on Seqera Platform (formerly Nextflow Tower).

**Workflow:**
```
Seqera MCP (official)          mcp-genomic-results          mcp-patient-report
  Launch nf-core/sarek    -->    Parse VCF + CNS       -->    Generate report
  Monitor pipeline               Annotate mutations           Patient summary
  Retrieve outputs                HRD scoring
```

## Setup

### Prerequisites

- Node.js 18+ installed
- A Seqera Platform account ([cloud.seqera.io](https://cloud.seqera.io))
- A Seqera Platform access token
- A workspace ID

### Installation

The official Seqera MCP server runs via `npx` (no global install needed):

```bash
npx @seqeralabs/mcp-server-seqera
```

### Environment Variables

| Variable | Required | Description |
|----------|----------|-------------|
| `TOWER_ACCESS_TOKEN` | Yes | Seqera Platform API token |
| `TOWER_API_ENDPOINT` | No | API endpoint (default: `https://api.cloud.seqera.io`) |
| `TOWER_WORKSPACE_ID` | No | Default workspace ID |

## Claude Desktop Configuration

Add the Seqera official MCP server alongside the existing Precision Medicine servers in your `claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "seqera-official": {
      "command": "npx",
      "args": ["-y", "@seqeralabs/mcp-server-seqera"],
      "env": {
        "TOWER_ACCESS_TOKEN": "your-seqera-token-here",
        "TOWER_WORKSPACE_ID": "your-workspace-id"
      }
    },
    "genomic-results": {
      "command": "uv",
      "args": [
        "run", "--directory",
        "/path/to/spatial-mcp/servers/mcp-genomic-results",
        "python", "-m", "mcp_genomic_results"
      ],
      "env": {
        "GENOMIC_RESULTS_DRY_RUN": "false"
      }
    }
  }
}
```

## End-to-End Workflow

### 1. Launch Variant Calling Pipeline (Seqera Official)

Ask Claude to launch an nf-core/sarek pipeline via the Seqera MCP server:

```
Launch nf-core/sarek for Patient PAT001 with these settings:
- Input: s3://bucket/PAT001/fastq/
- Genome: GRCh38
- Tools: Mutect2, CNVkit
- Output: s3://bucket/PAT001/results/
```

### 2. Monitor Pipeline (Seqera Official)

```
Check the status of workflow wf_abc123. Is it still running?
```

### 3. Parse Results (mcp-genomic-results)

Once the pipeline completes, parse the outputs:

```
Parse the somatic variants from data/patient-data/PAT001-OVC-2025/genomics/somatic_variants.vcf
and copy number results from data/patient-data/PAT001-OVC-2025/genomics/copy_number_results.cns.
Generate a comprehensive genomic report with therapy recommendations.
```

### 4. Generate Patient Report (mcp-patient-report)

```
Create a patient-facing report for PAT001 incorporating the genomic findings.
```

## Relationship to Previous Mock Server

The mocked `mcp-seqera` server (3 tools, synthetic data) has been removed from this repository. The official `@seqeralabs/mcp-server-seqera` now fully replaces it with **real** Seqera Platform integration:

- Actual pipeline launches on cloud compute
- Real-time monitoring and log streaming
- Compute environment management
- Dataset and pipeline management

**Setup instructions:** [CONNECT_EXTERNAL_MCP.md](../for-researchers/CONNECT_EXTERNAL_MCP.md)

## References

- [Seqera MCP Server (npm)](https://www.npmjs.com/package/@seqeralabs/mcp-server-seqera)
- [Seqera Platform Documentation](https://docs.seqera.io/)
- [nf-core/sarek Pipeline](https://nf-co.re/sarek)
- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)
- [Test-9: Seqera Connector Demo](../reference/testing/patient-one/test-prompts/DRY_RUN/test-9-e2e-seqera-connector.md)
