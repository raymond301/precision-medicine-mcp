# Section 05: Data Protection

[← Back to SLA Overview](README.md)

---

## 5. Storage Infrastructure & Data Protection

### 5.1 Storage Types

| Storage Type | Purpose | Encryption | Backup | HIPAA Compliant |
|--------------|---------|------------|--------|-----------------|
| Cloud Storage (uploads) | Patient input data | AES-256 (Google-managed keys) | Daily incremental | ✅ Yes (BAA signed) |
| Cloud Storage (outputs) | Analysis results | AES-256 (Google-managed keys) | Daily incremental | ✅ Yes |
| Healthcare API FHIR Store | Clinical data | AES-256 + Customer-managed keys (CMEK) | Continuous (point-in-time recovery) | ✅ Yes |
| Cloud Logging (audit logs) | HIPAA audit trail | AES-256 (immutable, tamper-proof) | Real-time export to BigQuery | ✅ Yes |

### 5.2 Data Retention Policy

- **Active Path:** Nearline after 90 days.
- **Long-term Storage:** Archive after 7 years.
- **Deletion:** Complete purge after 10 years.

---

## Related Documents

- [Section 04: Compliance and Regulatory](04_COMPLIANCE_REGULATORY.md)
- [Section 07: Backup and Disaster Recovery](07_BACKUP_DISASTER_RECOVERY.md)

**Next Section:** [06. Security Operations →](06_SECURITY_OPERATIONS.md)

---

### Document History

| Date | Version | Author | Change Summary |
|:---|:---|:---|:---|
| 2026-02-08 | 2.0 | IT Operations | Initial creation |
