# Appendix D: HIPAA Compliance Checklist

[← Back to SLA Overview](README.md)

---

## Monthly Verification Checklist

- [ ] **Encryption at Rest:** Verify all storage buckets and FHIR stores have encryption enabled.
- [ ] **Data-in-Transit:** Verify HTTPS/TLS 1.3 enforcement on all endpoints.
- [ ] **Audit Logging:** Verify 100% of access attempts were captured in Cloud Logging.
- [ ] **MFA Compliance:** Verify all admin accounts have MFA active in Azure AD.
- [ ] **Patch Status:** Confirm all container images were updated during the maintenance window.

---

## Related Documents

- [Section 04: Compliance and Regulatory](04_COMPLIANCE_REGULATORY.md)
- [Section 06: Security Operations](06_SECURITY_OPERATIONS.md)

---

### Document History

| Date | Version | Author | Change Summary |
|:---|:---|:---|:---|
| 2026-02-08 | 2.0 | IT Operations | Initial creation |
