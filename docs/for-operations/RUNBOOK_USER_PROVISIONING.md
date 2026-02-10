# Runbook: User Provisioning and Offboarding

[← Back to SLA Overview](README.md)

---

## Overview
This runbook defines the procedures for managing user access to the Precision Medicine MCP Platform, ensuring HIPAA-compliant "Least Privilege" access.

## 1. User Provisioning (Onboarding)

### Step 1: Identity Creation
All users must have a valid hospital identity in Azure AD (@hospital.org).

### Step 2: Group Assignment
Assign the user to the appropriate Azure AD security group:
- **`precision-medicine-viewers`**: Logical "Platform Viewer" role.
- **`precision-medicine-users`**: Logical "Platform Editor" role.
- **`precision-medicine-admins`**: Logical "Infrastructure Admin" role.

### Step 3: Secret/Certificate Enrollment
If the user requires programmatic access (Bioinformaticians):
1. Issue a temporary OIDC identity token.
2. Verify MFA enrollment in Azure AD.

## 2. User Offboarding (Termination)

### Step 1: Immediate Revocation
Upon notification of termination (or via automated HR feed sync):
1. **Azure AD**: Disable account or remove from all `precision-medicine-*` groups.
2. **GCP IAM**: Revoke any direct IAM bindings (if applicable).

### Step 2: Access Refresh
Force a sign-out of all sessions to invalidate existing JWT/OIDC tokens and session cookies.

---

## Verification checklist
- [ ] User added to Azure AD group.
- [ ] User can log in via SSO.
- [ ] Offboarded user cannot log in (verify within 1 hour).

---

### Document History

| Date | Version | Author | Change Summary |
|:---|:---|:---|:---|
| 2026-02-08 | 2.0 | IT Operations | Initial creation |
