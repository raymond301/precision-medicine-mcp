# Section 04: Compliance and Regulatory

[← Back to SLA Overview](README.md)

---

## 4. Compliance and Regulatory

### 4.1 HIPAA Privacy & Security Safeguards

The platform adheres to the three pillars of HIPAA security:

| Category | Control Implementation | Standard |
| :--- | :--- | :--- |
| **Technical** | Encryption at rest/transit, IAM, Audit logging | §164.312 |
| **Administrative** | Annual risk assessments, HIPAA training for all staff | §164.308 |
| **Physical** | Cloud Data Center physical security (GCP/Azure) | §164.310 |

### 4.2 Business Associate Agreement (BAA)

- **Requirement:** A signed BAA must be in place before any PHI is processed.
- **Vendor Management:** The platform only utilizes "HIPAA-compliant" cloud services covered by the provider's BAA (e.g., Google Cloud BAA, Microsoft BAA).
- **Sub-contractors:** All sub-contractors with access to the platform environment must execute a BAA with the Primary Provider.

### 4.3 NIST Cybersecurity Framework (CSF) Alignment

The platform's security controls are mapped to **NIST CSF v2.0**:
- **Govern (GV):** Integrated risk management and organizational leadership oversight.
- **Identify (ID):** Automated asset inventory (Terraform/Cloud Asset Inventory).
- **Protect (PR):** Zero-Trust architecture and identity-centric access control.
- **Detect (DE):** Real-time monitoring via Cloud Guard Duty and Security Command Center.
- **Respond (RS):** Documented Incident Response runbooks for PHI breaches.
- **Recover (RC):** Multi-region backup and 1-hour RTO for clinical services.

### 4.4 Required Certifications

| Certification | Status | Target Date |
| :--- | :--- | :--- |
| **SOC 2 Type II** | ✅ Compliant (Provider Native) | Annual renewal |
| **HIPAA Attestation** | ✅ Verified | Periodic Audit |
| **ISO 27001** | 📅 In Progress | 2026 Q4 |
| **HITRUST CSF** | 🚀 Planned | 2027 Q2 |

### 4.5 State and International Compliance

- **CCPA/CPRA (California):** Privacy notices and consumer rights implemented.
- **GDPR (EU):** Data processing agreements (DPAs) and SCCs available upon request.

---

## Related Documents

- [Appendix D: HIPAA Compliance Checklist](APPENDIX_D_HIPAA_CHECKLIST.md)
- [Section 05: Data Protection](05_DATA_PROTECTION.md)

**Next Section:** [05. Data Protection →](05_DATA_PROTECTION.md)

---

### Document History

| Date | Version | Author | Change Summary |
|:---|:---|:---|:---|
| 2026-02-08 | 2.0 | IT Operations | Initial creation |
