# Section 07: Backup and Disaster Recovery

[← Back to SLA Overview](README.md)

---

## 7. Backup and Disaster Recovery

### 7.1 Recovery Objectives

| Metric | Clinical Tier (Target) | Admin Tier (Target) | Commitment |
|--------|------------------------|---------------------|---------------|
| **RTO (Recovery Time)** | **≤ 1 Hour** | ≤ 4 Hours | ✅ Achievable |
| **RPO (Data Loss)** | **≤ 15 Minutes** | ≤ 1 Hour | ✅ Achievable |
| **Data Durability** | 11 9's | 11 9's | ✅ Guaranteed |

### 7.2 Disaster Recovery Strategy

- **Primary Region:** us-central1
- **DR Region:** us-east1 (warm standby)
- **Failover:** DNS-based cutover via Cloud DNS.

---

## Related Documents

- [Runbook: Disaster Recovery Failover](RUNBOOK_DR_FAILOVER.md)
- [Section 05: Data Protection](05_DATA_PROTECTION.md)

**Next Section:** [08. Monitoring and Reporting →](08_MONITORING_REPORTING.md)

---

### Document History

| Date | Version | Author | Change Summary |
|:---|:---|:---|:---|
| 2026-02-08 | 2.0 | IT Operations | Initial creation |
