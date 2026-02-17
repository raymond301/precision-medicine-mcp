# DRY_RUN Mode (Canonical Reference)

## Purpose

DRY_RUN mode enables development and testing with realistic synthetic responses, without requiring:
- External tool dependencies (fgbio, STAR, DeepCell, etc.)
- Real data files or GCS access
- API credentials or authentication
- Large compute resources

## Activation

Set the server's DRY_RUN environment variable to `true`:

```bash
export {SERVER_NAME}_DRY_RUN=true
```

### Per-Server Variable Names

| Server | Variable |
|--------|----------|
| mcp-spatialtools | `SPATIAL_DRY_RUN` |
| mcp-fgbio | `FGBIO_DRY_RUN` |
| mcp-multiomics | `MULTIOMICS_DRY_RUN` |
| mcp-perturbation | `PERTURBATION_DRY_RUN` |
| mcp-deepcell | `DEEPCELL_DRY_RUN` |
| mcp-cell-classify | `CELL_CLASSIFY_DRY_RUN` |
| mcp-openimagedata | `IMAGE_DRY_RUN` |
| mcp-tcga | `TCGA_DRY_RUN` |
| mcp-genomic-results | `GENOMIC_RESULTS_DRY_RUN` |
| mcp-quantum-celltype-fidelity | `QUANTUM_DRY_RUN` |

## When to Use

- **Local development** — Test without installing bioinformatics tools
- **CI/CD pipelines** — Automated testing with predictable outputs
- **Demonstrations** — Show capabilities without real patient data
- **Onboarding** — New developers can explore immediately
- **Education** — Classroom use with ~$0.32 total cost

## What DRY_RUN Returns

Each server returns clinically realistic synthetic data based on the PatientOne case study. Responses include:
- Realistic data structures matching production output
- Biologically plausible values calibrated to published studies
- Appropriate metadata and annotations

## Cost

- **DRY_RUN:** ~$1 total for full PatientOne workflow (25-35 minutes, tokens only)
- **Production:** $24-104 per patient analysis (compute + API costs)

## Production

Always set `DRY_RUN=false` in production deployments. Production mode requires:
- Real data files accessible via GCS or local paths
- Bioinformatics tools installed (or Cloud Run containers)
- Appropriate API credentials

---

**Server installation:** [server-installation.md](server-installation.md)
**Cost analysis:** [cost-analysis.md](cost-analysis.md)
