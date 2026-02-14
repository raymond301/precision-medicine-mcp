# Cost Analysis (Canonical Reference)

## Per-Patient Cost Comparison

| Analysis Mode | Traditional | MCP Platform | Savings |
|---------------|-------------|-------------|---------|
| **DRY_RUN** (synthetic data, demo) | N/A | ~$1 | — |
| **Production compute** (API + Cloud Run) | $6,000-9,000 | $25-104 | ~$6,000+ |
| **Production total** (incl. infrastructure) | $6,000-9,000 | $324-702 | ~$3,137 avg |

## Cost Breakdown by Test (DRY_RUN Mode)

| Test | DRY_RUN Time | DRY_RUN Cost | Production Time | Production Cost |
|------|-------------|-------------|----------------|----------------|
| TEST_1: Clinical & Genomic | 4-6 min | ~$0.15 | 15-30 min | $3-8 |
| TEST_2: Multi-omics Integration | 5-8 min | ~$0.25 | 30-60 min | $5-15 |
| TEST_3: Spatial Transcriptomics | 4-6 min | ~$0.20 | 45-90 min | $8-17 |
| TEST_4: Histology & Imaging | 3-5 min | ~$0.15 | 40-90 min | $10-36 |
| TEST_5: Integrated Report | 5-8 min | ~$0.25 | 15-30 min | $5-12 |
| **TOTAL** | **25-35 min** | **~$1** | **2-5 hours** | **$31-88** |

## Monthly Infrastructure Costs

| Component | Monthly Cost | Notes |
|-----------|-------------|-------|
| Cloud Run (15 servers) | ~$50-200 | Scale-to-zero when idle |
| Cloud Storage (GCS) | ~$5-20 | Patient data + results |
| Cloud Logging | ~$10-30 | 10-year HIPAA retention |
| Networking (VPC, NAT) | ~$30-50 | VPN + private access |
| Claude API tokens | ~$500-2,000 | Usage-dependent |
| **Total infrastructure** | **~$600-2,300/mo** | — |

## Scaling Projections

| Scale | Patients/Month | Monthly Cost | Cost/Patient |
|-------|---------------|-------------|-------------|
| Pilot | 10 | ~$800 | ~$80 |
| Small | 50 | ~$1,800 | ~$36 |
| Medium | 100 | ~$3,000 | ~$30 |
| Large | 500 | ~$8,000 | ~$16 |

---

**Full value proposition:** [value-proposition.md](value-proposition.md)
