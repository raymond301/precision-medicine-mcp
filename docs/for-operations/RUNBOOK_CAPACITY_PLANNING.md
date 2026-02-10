# Runbook: Capacity Planning

[← Back to SLA Overview](README.md)

---

## Procedure

1. **Monthly Review:** Analyze `run.googleapis.com/request_count` trends.
2. **Forecast:** Calculate projected load for next 90 days.
3. **Triggers:** If sustained usage >80% max instances, increase Terraform max instance ceiling.
4. **Quota:** Request GCP project quota increases if forecasted load exceeds 1000 CPUS.

---

## Related Documents

- [Section 02: Service Level Objectives](02_SERVICE_LEVEL_OBJECTIVES.md)
- [Section 09: Financial Terms](09_FINANCIAL_TERMS.md)

---

### Document History

| Date | Version | Author | Change Summary |
|:---|:---|:---|:---|
| 2026-02-08 | 2.0 | IT Operations | Initial creation |
