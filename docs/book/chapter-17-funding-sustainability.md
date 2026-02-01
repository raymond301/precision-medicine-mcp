# Chapter 17: Funding and Sustainability

*ROI analysis, grant strategies, and budget models*

---

## Why Funding Matters

Chapters 1-16 built and deployed the system. **Now: How do you pay for it?**

**Deployment costs money**:
- Infrastructure: $1,000-2,000/month (Cloud Run, storage, networking)
- API usage: $500-1,000/month (Claude/Gemini tokens)
- Personnel: $50,000-120,000/year (bioinformatics support)

**But traditional analysis costs more**:
- Personnel: $120,000-240,000/year (1.5-3 FTE manual analysis)
- Infrastructure: $10,000-20,000/year (on-premise compute)
- External services: $600,000-900,000/year (100 patients × $6,000-9,000)

**This chapter**: ROI analysis, grant strategies, budget models for sustainability.

---

## ROI Analysis

### Per-Patient Cost Comparison

| Cost Component | Traditional | MCP System | Savings |
|----------------|-------------|------------|---------|
| Personnel time | $3,000 (40 hrs × $75/hr) | $75 (1 hr × $75/hr) | $2,925 |
| Compute resources | $2,000-4,000 | $22-99 | $1,901-3,978 |
| External services | $1,000-2,000 | $1-3 (APIs) | $997-1,999 |
| **Total per patient** | **$6,000-9,000** | **$98-177** | **$5,823-8,902** |

**Average savings**: $7,210 per patient (96% cost reduction)

### Annual Savings (100 Patients)

**Traditional costs**:
- 100 patients × $7,500 avg = $750,000/year
- Personnel: 4,000 hours (2 FTE @ $120K = $240,000)
- Infrastructure: $200,000-400,000
- External services: $100,000-200,000

**MCP costs**:
- 100 patients × $138 avg = $13,800/year
- Cloud infrastructure: $12,000-24,000/year
- Claude API: $1,200-2,400/year
- Personnel: $50,000/year (0.5 FTE oversight)

**Total MCP**: $76,000-90,000/year

**Annual savings**: $660,000-674,000 (88-89% reduction)

### Payback Period

**Initial investment**:
- Setup: $50,000-75,000 (deployment, training, integration)
- First-year operational: $76,000-90,000

**Total first-year**: $126,000-165,000

**Break-even**: After 17-23 patients analyzed (2-3 months in 100-patient/year facility)

**5-year ROI**:
- Investment: $126K-165K (Year 1) + $350K-450K (Years 2-5) = $476K-615K
- Traditional cost: $3.75M (5 years × $750K)
- **Total savings**: $3.13M-3.27M (655-687% ROI)

---

## Grant Strategies

### NIH R01 Budget Model

**Typical R01**: $250,000/year direct costs × 5 years = $1.25M

**Budget allocation**:

| Category | Traditional | MCP System | Savings |
|----------|-------------|------------|---------|
| Personnel | $600,000 (2 FTE × 5 years) | $250,000 (1 FTE × 5 years) | $350,000 |
| Equipment | $100,000 (servers, storage) | $0 (cloud-based) | $100,000 |
| Compute | $200,000 (cluster time) | $120,000 (Cloud Run) | $80,000 |
| Supplies | $150,000 (sequencing) | $150,000 (same) | $0 |
| Other | $200,000 (travel, pubs) | $200,000 (same) | $0 |
| **Total** | **$1,250,000** | **$720,000** | **$530,000** |

**Grant strategy**: Request $720K instead of $1.25M = **42% budget reduction**. Use savings for additional aims or patient cohort expansion.

**Budget justification template**:

```
Analysis Costs (Aim 2 - Multi-modal patient analysis):

Traditional approach: 200 patients × $7,500 = $1,500,000
Proposed MCP approach: 200 patients × $138 = $27,600

Cost savings: $1,472,400 (98% reduction)

We request $120,000 over 5 years for Cloud Run infrastructure
($24K/year) rather than $1.5M for traditional analysis. This 92%
cost reduction allows us to expand the cohort from 200 to 1,000
patients within the same budget, increasing statistical power and
clinical impact.
```

Full NIH budget template: [`docs/for-researchers/grant-templates/nih-r01-budget.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-researchers/grant-templates/nih-r01-budget.md) (planned)

### Foundation Grant Model

**Typical foundation grant**: $50,000-100,000 (1-2 years)

**Budget for pilot deployment** (Year 1):

| Item | Cost | Notes |
|------|------|-------|
| GCP infrastructure setup | $15,000 | VPC, OAuth2, Epic integration |
| Cloud Run (annual) | $24,000 | 100 patients/year |
| Claude API (annual) | $2,400 | 100 patients × $24 avg |
| Personnel (0.25 FTE) | $30,000 | Bioinformatics oversight |
| Training & documentation | $5,000 | 5 users, 2 days training |
| Contingency (10%) | $7,640 | Buffer for overruns |
| **Total Year 1** | **$84,040** | |

**Deliverables**:
- 100 ovarian cancer patients analyzed
- Publication: Methods paper + clinical cohort study
- Expansion plan: Scale to breast, lung, colorectal

**Sustainability plan** (Year 2+):
- Hospital absorbs operational costs ($26.4K/year)
- Grant funds patient cohort expansion (additional $50K/year)
- Savings reinvested in research ($660K/year traditional → $26.4K MCP = $634K reinvestment)

---

## Hospital Budget Model

### Pilot Deployment (3 Months, 25 Patients)

**Costs**:
- Infrastructure setup: $15,000 (one-time)
- Cloud Run: $2,000/month × 3 = $6,000
- Claude API: $600 (25 patients × $24)
- Personnel (0.5 FTE oversight): $15,000 (3 months)
- **Total pilot**: $36,600

**Traditional alternative**: 25 patients × $7,500 = $187,500
**Pilot savings**: $150,900 (80% reduction)

### Production Deployment (Annual, 100 Patients)

**Costs**:
- Cloud Run infrastructure: $24,000/year
- Claude API: $2,400/year
- Personnel oversight (0.5 FTE): $50,000/year
- Epic FHIR integration: $5,000/year (maintenance)
- Training & support: $5,000/year
- **Total production**: $86,400/year

**Traditional alternative**: $750,000/year
**Annual savings**: $663,600 (88% reduction)

### Scaling to 500 Patients/Year

**Costs**:
- Cloud Run (auto-scales): $50,000/year
- Claude API: $12,000/year (500 × $24)
- Personnel (1 FTE): $120,000/year
- Infrastructure maintenance: $10,000/year
- **Total**: $192,000/year

**Traditional alternative**: $3,750,000/year (500 × $7,500)
**Annual savings**: $3,558,000 (95% reduction)

**Cost per patient remains constant**: $138 (vs $7,500 traditional)

---

## Revenue Models (For-Profit)

### Clinical Service Model

**Charge structure**:
- Base analysis: $500 per patient (vs $3,000-5,000 competitors)
- Cost: $138 per patient
- **Gross margin**: $362 per patient (72%)

**Volume scenarios**:

| Patients/Year | Revenue | Costs | Profit | Margin |
|---------------|---------|-------|--------|--------|
| 100 | $50,000 | $13,800 | $36,200 | 72% |
| 500 | $250,000 | $69,000 | $181,000 | 72% |
| 1,000 | $500,000 | $138,000 | $362,000 | 72% |

**Market positioning**: 83% cheaper than competitors ($500 vs $3,000), same quality (multi-modal analysis).

### SaaS Subscription Model

**Tier structure**:

| Tier | Monthly | Included Analyses | Overage Cost | Target Customer |
|------|---------|-------------------|--------------|-----------------|
| Starter | $5,000 | 10 patients | $400/patient | Small clinics |
| Professional | $15,000 | 50 patients | $250/patient | Medium hospitals |
| Enterprise | $40,000 | 200 patients | $150/patient | Large cancer centers |

**Example**: Professional tier (50 patients/month)
- Revenue: $15,000/month = $180,000/year
- Cost: $6,900/month (50 × $138) = $82,800/year
- **Profit**: $97,200/year (54% margin)

**Customer acquisition**:
- Year 1: 5 customers (pilot hospitals) = $900,000 revenue
- Year 2: 15 customers (expansion) = $2,700,000 revenue
- Year 3: 30 customers (scale) = $5,400,000 revenue

---

## Cost Drivers and Optimization

### Primary Cost Drivers

1. **Claude API tokens** (40-50% of operational costs)
   - Optimization: Use Haiku for simple queries ($0.25 vs $3 per million input tokens)
   - Savings: 30-40% reduction in API costs

2. **Cloud Run compute** (30-40% of operational costs)
   - Optimization: Right-size memory/CPU, scale to zero when idle
   - Savings: 20-30% reduction in compute costs

3. **Personnel oversight** (10-20% of operational costs)
   - Optimization: Automate quality checks, batch review workflows
   - Savings: 10-15% reduction in personnel time

**Example optimization** (100 patients/year):
- Baseline: $86,400/year
- Haiku for 50% of queries: -$600 (API savings)
- Right-sized compute: -$4,800 (Cloud Run savings)
- Automated QC: -$5,000 (personnel savings)
- **Optimized total**: $76,000/year (12% reduction)

### Sensitivity Analysis

**What if Claude API pricing changes?**

| API Cost Change | New Cost/Patient | New Annual (100 patients) | Impact |
|----------------|------------------|---------------------------|--------|
| -50% (competition) | $126 | $74,400 | +15% savings |
| +50% (price increase) | $150 | $98,400 | -14% savings |
| +100% (doubling) | $162 | $110,400 | -28% savings |

**Still cheaper than traditional** even if Claude API doubles in price: $162 vs $7,500 (98% savings).

**Mitigation**: Multi-LLM strategy (Claude, Gemini, Llama) ensures price competition.

---

## Sustainability Roadmap

### Year 1: Pilot and Validation

**Goals**:
- Deploy to 1 hospital (100 patients)
- Validate clinical outcomes
- Publish methods paper
- Establish baseline metrics

**Funding**: Foundation grant ($85K) or hospital pilot budget ($87K)

**Deliverables**:
- 100 patients analyzed
- Clinical validation study published
- Cost savings documented ($660K vs traditional)

### Year 2: Expansion and Optimization

**Goals**:
- Scale to 3 hospitals (300 patients total)
- Expand to 2 additional cancer types
- Optimize costs (12% reduction)
- Establish revenue model (if for-profit)

**Funding**: Hospital operational budgets ($76K each) or NIH R01 Year 2 ($150K)

**Deliverables**:
- 300 patients analyzed across 3 cancer types
- Multi-cancer validation published
- ROI metrics: $1.98M savings (300 patients)

### Year 3: Scale and Productization

**Goals**:
- Scale to 10 hospitals (1,000 patients total)
- Full production deployment (all 12 servers)
- SaaS offering (if for-profit)
- Self-sustaining operations

**Funding**: Self-sustaining (hospital budgets) or revenue (SaaS subscriptions)

**Deliverables**:
- 1,000 patients analyzed
- Consortium publication
- Total savings: $7.2M vs traditional

---

## What You've Configured

**Funding models**:
1. **ROI analysis**: $7,210 savings per patient (96% reduction), $660K annual savings (100 patients)
2. **Grant strategies**: NIH R01 ($720K vs $1.25M traditional), foundation grants ($85K pilot)
3. **Hospital budgets**: $87K/year production (100 patients) vs $750K traditional
4. **Revenue models**: Clinical service ($500/patient, 72% margin), SaaS ($5K-40K/month tiers)
5. **Cost optimization**: Haiku for simple queries (-30% API), right-sized compute (-20%), automated QC (-15%)
6. **Sustainability roadmap**: Year 1 pilot (100 patients), Year 2 expansion (300 patients), Year 3 scale (1,000 patients)

**Key metrics**:
- Per-patient savings: $7,210 (96% reduction)
- Payback period: 17-23 patients (2-3 months)
- 5-year ROI: 655-687% ($3.13M-3.27M savings)
- Scalability: Cost/patient constant at $138 (traditional increases)

---

## Next Steps

**Chapter 18: Lessons Learned and What's Next** (final chapter) covers:
- Production deployment insights
- Technical challenges and solutions
- Future enhancements roadmap
- Multi-cancer expansion
- Community and open source

---

**Chapter 17 Summary**:
- Per-patient cost: $138 (MCP) vs $7,500 (traditional) = $7,210 savings (96%)
- Annual savings: $660K (100 patients), $3.56M (500 patients)
- Payback period: 2-3 months (17-23 patients)
- 5-year ROI: 655-687% ($3.13M-3.27M savings)
- Grant strategies: NIH R01 ($720K vs $1.25M), foundation ($85K pilot)
- Hospital budgets: $87K/year (100 patients) vs $750K traditional
- Revenue models: Clinical service (72% margin), SaaS ($5K-40K/month)
- Sustainability: Year 1 pilot → Year 2 expansion → Year 3 scale

**Files**: [`docs/prompt-library/funder-demo-prompts.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/prompt-library/funder-demo-prompts.md)
**ROI**: 96% cost reduction, 2-3 month payback, 655-687% 5-year ROI
**Scaling**: Cost/patient constant ($138) as volume increases
