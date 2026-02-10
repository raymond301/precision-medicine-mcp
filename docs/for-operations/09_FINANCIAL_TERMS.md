# Section 09: Financial Terms

[← Back to SLA Overview](README.md)

---

## 9. Financial Terms

### 9.1 Operational Budget

- **Target:** $2,500/month
- **Infrastructure:** $1,000 (GCP Compute/Storage)
- **API Costs:** $1,500 (Anthropic Claude orchestration)

#### 9.1.1 Cost Analysis by Mode

| Mode | Est. Wait Time | Est. Cost per Run | Best For |
| :--- | :--- | :--- | :--- |
| **DRY_RUN** | 25-35 min | ~$1 | Development, demo, CI/CD |
| **Automated Report** | ~12 seconds | ~$1 | Quick analysis (pre-aligned) |
| **Real Data (Small)** | 1-3 hours | $7 - $29 | Workflow validation |
| **Production** | 4-8 hours | **$25 - $120** | Full hospital sequencing data |

*Note: Production costs vary based on whether data is pre-aligned (Space Ranger) or raw FASTQ.*

### 9.2 Service Credits

| Monthly Uptime Achieved | Service Credit |
|------------------------|----------------|
| ≥ 99.99% | 0% |
| 99.9% – < 99.99% | 10% |
| 99.0% – < 99.9% | 25% |
| 95.0% – < 99.0% | 50% |
| < 95.0% | 100% |

### 9.3 Cost Optimization Requirements

To maintain the $2,500/month target, the following strategies must be implemented:
1. **Batch Processing**: Amortize setup costs by processing ≥10 patients in sequence.
2. **Intermediate Caching**: Store aligned BAMs and segmentation masks to avoid redundant compute.
3. **Selective Analysis**: Disable imaging or omics processing for patients where data is missing.
4. **Auto-Scaling**: Ensure Cloud Run is configured to scale to zero when idle.

---

## Related Documents

- [Section 02: Service Level Objectives](02_SERVICE_LEVEL_OBJECTIVES.md)
- [Comprehensive Cost & Budget Guide](../for-hospitals/operations/cost-and-budget.md) (Full technical breakdown)
- [Section 10: Governance and Change Management](10_GOVERNANCE_CHANGE_MANAGEMENT.md)

**Next Section:** [10. Governance and Change Management →](10_GOVERNANCE_CHANGE_MANAGEMENT.md)

---

### Document History

| Date | Version | Author | Change Summary |
|:---|:---|:---|:---|
| 2026-02-08 | 2.0 | IT Operations | Initial creation |
