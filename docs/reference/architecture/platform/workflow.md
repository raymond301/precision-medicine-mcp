# Workflow Orchestration Architecture

**Status:** Served by external Seqera MCP server
**Last Updated:** 2026-02-21

---

## Overview

Workflow orchestration for Nextflow bioinformatics pipelines (nf-core/sarek, rnaseq, spatial, etc.) is now provided by the **external Seqera MCP server** (`@seqeralabs/mcp-server-seqera`, 7 tools) rather than the previous 3-tool mocked `mcp-seqera` server.

**Setup instructions:** [CONNECT_EXTERNAL_MCP.md](../../../../docs/for-researchers/CONNECT_EXTERNAL_MCP.md)

**Integration details:** [SEQERA_OFFICIAL_MCP.md](../../../../docs/for-developers/SEQERA_OFFICIAL_MCP.md)

---

## Cross-Server Workflow

The external Seqera MCP server integrates with other platform servers in the same pipeline chain:

```
Seqera MCP (external)          mcp-genomic-results          mcp-patient-report
  Launch nf-core/sarek    -->    Parse VCF + CNS       -->    Generate report
  Monitor pipeline               Annotate mutations           Patient summary
  Retrieve outputs                HRD scoring
```

---

## Related

- [Architecture Overview](../README.md)
- [Server Registry](../../shared/server-registry.md)
- [Genomic Results Architecture](../dna/genomic-results.md)
