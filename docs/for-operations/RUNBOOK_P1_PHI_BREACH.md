# Runbook: P1 PHI Breach Response

[← Back to SLA Overview](README.md)

---

## Procedure

1. **IMMEDIATE:** Block suspicious user or IP in Cloud Armor.
2. **Identify Scope:** Search Cloud Logging for principalEmail access to FHIR stores.
3. **Notify CISO:** Escalate within 1 hour of detection.
4. **Mandatory Report:** Determine if breach affects >500 patients (requires HHS notification).

---

## Related Documents

- [Section 04: Compliance and Regulatory](04_COMPLIANCE_REGULATORY.md)
- [Section 06: Security Operations](06_SECURITY_OPERATIONS.md)

---

### Document History

| Date | Version | Author | Change Summary |
|:---|:---|:---|:---|
| 2026-02-08 | 2.0 | IT Operations | Modular alignment and PDF spec compliance |
| 2026-01-22 | 1.0 | IT Operations | Initial monolithic document (original.md) |
