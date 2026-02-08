# Section 11: Contract Termination

[← Back to SLA Overview](README.md)

---

## 11. Contract Termination

### 11.1 Data Transfer and Deletion Protocols

- **Data Export:** Upon termination, all hospital data (genomic raw data, clinical reports, and audit logs) will be provided via a secure GCP Storage bucket in industry-standard formats:
  - **Clinical Data:** FHIR R4 / NDJSON.
  - **Genomics:** VCF / BAM / FASTQ.
  - **Reports:** PDF / JSON.
- **Secure Purging:** Following confirmed receipt of data, all data will be securely purged using **NIST SP 800-88** compliant deletion methods within 30 days.

### 11.2 Access Revocation Timeline

- **Immediate Revocation:** Within **1 hour** of official contract termination, the following will be disabled:
  - Azure AD SSO integration and all user accounts.
  - Service Account keys and Secret Manager access.
  - All VPN and network-level ingress/egress.

### 11.3 Post-Contract Audit Rights

- **Compliance Window:** The hospital retains the right to audit the platform's deletion logs and infrastructure state for up to **12 months** post-termination to verify total data erasure.
- **Cost:** Such audits are hosted at the hospital's expense unless a breach or failure to delete is discovered.

### 11.4 Knowledge Transfer Procedures

- **Technical Handoff:** IT Operations will provide a final "System State" report including:
  - Final resource utilization and billing reconciliation.
  - Inventory of any hospital-owned hardware or dedicated cloud instances.
  - Documentation of current API versions and data schemas used during the contract period.

---

## Related Documents

- [Section 04: Compliance and Regulatory](04_COMPLIANCE_REGULATORY.md)
- [Section 05: Data Protection](05_DATA_PROTECTION.md)

**Next Section:** [12. Escalation and Support →](12_ESCALATION_SUPPORT.md)

---

### Document History

| Date | Version | Author | Change Summary |
|:---|:---|:---|:---|
| 2026-02-08 | 2.0 | IT Operations | Initial creation |
