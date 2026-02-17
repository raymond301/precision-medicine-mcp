# Value Proposition (Canonical Reference)

## Time Savings

| Mode | Traditional | MCP Platform | Speed-up |
|------|-------------|-------------|----------|
| **DRY_RUN** (demo, synthetic data) | 40 hours | **25-35 minutes** | ~68x faster |
| **Production** (real patient data) | 40 hours | **2-5 hours** | 8-20x faster |

**Note on time modes:**
- **DRY_RUN (25-35 min):** Synthetic data, no external API calls — ideal for demos and education
- **Production (2-5 hours):** Real data analysis with full compute (sequencing alignment, spatial analysis, imaging segmentation)
- **Traditional (40 hours):** Manual bioinformatics across multiple tools and teams

**Annual capacity** (1 bioinformatician): ~50 patients (traditional) vs. ~400-1,000 patients (production MCP)

## Cost Savings

| Analysis Mode | Traditional Cost | MCP Cost | Savings |
|---------------|-----------------|----------|---------|
| **DRY_RUN** (demo, synthetic data) | N/A | ~$1 (tokens only) | — |
| **Production** (per-analysis compute) | $6,000-9,000 | $24-104 | ~$6,000+ |
| **Production** (total per patient, incl. infrastructure) | $6,000-9,000 | $324-702 | ~$3,137 avg |

**Note on cost ranges:**
- **$24-102** = per-analysis compute cost (Claude API tokens, Cloud Run, data transfer)
- **$324-702** = total per-patient cost including amortized infrastructure, storage, and overhead
- **~$1** = DRY_RUN mode with synthetic data (tokens only)

## ROI Summary

- **Savings per patient:** ~$3,098-3,176 (production, total cost basis)
- **Monthly pilot** (5 users, 100 patients): ~$2,400-9,200
- **Annual production** (20 users, 500 patients): ~$12,000-51,000
- **Break-even:** ~2-3 patients per month covers infrastructure costs

**Annual savings projections:**
- 100 patients/year: ~$313,700 saved
- 500 patients/year: ~$1,568,500 saved

## What Makes This Possible

- **Multiple MCP servers** with specialized bioinformatics tools (see [Server Registry](server-registry.md))
- **AI orchestration** via natural language (Claude + Gemini)
- **Multi-modal integration** across 5 data types in a single workflow
- **DRY_RUN mode** for zero-risk testing and demonstrations

---

**Full cost analysis:** [cost-analysis.md](cost-analysis.md)
**Full platform overview:** [README.md](README.md)
