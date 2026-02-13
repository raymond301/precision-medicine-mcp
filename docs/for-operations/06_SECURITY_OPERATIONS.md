# Section 06: Security Operations

[← Back to SLA Overview](README.md)

---

## 6. Security Operations

### 6.1 Identity and Access Management (IAM)

| Component | Policy | Enforcement | Audit |
|---------------|--------|-------------|-------|
| User Access | Azure AD SSO | MFA mandatory | Cloud Logging |
| Service Accounts | Least privilege | Terraform-managed | Key rotation: 90 days |
| Admin Access | Breakglass accounts | Just-in-time access | Alert on usage |

### 6.2 Security Monitoring

- **Vulnerability Scanning:** Continuous (Qualys/Defender).
- **Penetration Testing:** Annual external engagement.
- **Patch Management:** 
  - **Zero-day vulnerabilities:** ≤ 4 hours (immediate deployment).
  - **High-risk patches:** ≤ 24 hours.
  - **Standard updates:** Monthly maintenance window.

### 6.3 Generic Role-Based Access Control (RBAC) Matrix

To ensure portability across cloud providers (GCP, Azure, AWS), access is governed by **Logical Personas** mapped to standard cloud permission levels:

| Logical Role | Cloud IAM Equivalent | Data Permissions | Platform Capability |
| :--- | :--- | :--- | :--- |
| **Platform Viewer** | `Reader` / `Viewer` | Read-only de-identified core data | Access dashboards, view reports |
| **Platform Editor** | `Editor` / `Contributor` | Read/Write analysis & research data | Execute workflows, trigger analyses |
| **Infrastructure Admin** | `Owner` / `IAM Admin` | Metadata and service settings only | Manage secrets, deploy code, update IAM |
| **Platform Auditor** | `Security Reviewer` | Comprehensive audit & logging trails | Access immutable logs (no functional data) |

### 6.4 Authorization Enforcement Standards

Regardless of the underlying cloud infrastructure, the platform enforces:
1. **Least Privilege**: Users are assigned the lowest tier logical role required for their clinical or research task.
2. **Context-Aware Access**: Authorization is verified at the API Gateway layer using standard JWT/OIDC identity tokens.
3. **MFA Neutrality**: Multi-Factor Authentication is required at the primary Identity Provider (IdP) level before cloud resources are accessible.
4. **Agent Identity**: AI Agents are assigned unique Service Identities with scoped logical roles, preventing privilege escalation.

---

## Related Documents

- [Section 04: Compliance and Regulatory](04_COMPLIANCE_REGULATORY.md)
- [Runbook: Secret and Key Rotation](RUNBOOK_KEY_ROTATION.md)
- [Runbook: P1 PHI Breach Response](RUNBOOK_P1_PHI_BREACH.md)

**Next Section:** [07. Backup and Disaster Recovery →](07_BACKUP_DISASTER_RECOVERY.md)

---

### Document History

| Date | Version | Author | Change Summary |
|:---|:---|:---|:---|
| 2026-02-08 | 2.0 | IT Operations | Initial creation |
