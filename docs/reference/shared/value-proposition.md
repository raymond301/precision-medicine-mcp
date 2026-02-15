# Value Proposition (Canonical Reference)

## Time Savings

| Metric | Traditional | MCP Platform |
|--------|-------------|-------------|
| **Time per patient** | 40 hours manual bioinformatics | **35 minutes** AI-orchestrated |
| **Speed improvement** | — | **68x faster** |
| **Annual capacity** (1 bioinformatician) | ~50 patients | ~3,400 patients |

## Cost Savings

| Analysis Type | Traditional Cost | MCP Cost | Savings |
|---------------|-----------------|----------|---------|
| **DRY_RUN** (demo, synthetic data) | $6,000-9,000 | ~$1 | ~$6,000+ |
| **Production** (per-analysis compute) | $6,000-9,000 | $25-104 | ~$6,000+ |
| **Production** (total per patient, incl. infrastructure) | $6,000-9,000 | $324-702 | ~$3,137 avg |

**Note on cost ranges:**
- **$25-104** = per-analysis compute cost (Claude API tokens, Cloud Run, data transfer)
- **$324-702** = total per-patient cost including amortized infrastructure, storage, and overhead
- **$1** = DRY_RUN mode with synthetic data (tokens only)

## ROI Summary

- **Savings per patient:** ~$3,098-3,176 (production, total cost basis)
- **Monthly pilot** (5 users, 100 patients): ~$2,400-9,200
- **Annual production** (20 users, 500 patients): ~$12,000-51,000
- **Break-even:** ~2-3 patients per month covers infrastructure costs

## What Makes This Possible

- **15 MCP servers** with 80 specialized bioinformatics tools
- **AI orchestration** via natural language (Claude + Gemini)
- **Multi-modal integration** across 5 data types in a single workflow
- **DRY_RUN mode** for zero-risk testing and demonstrations

---

**Full cost analysis:** [cost-analysis.md](cost-analysis.md)
**Full platform overview:** [README.md](README.md)
