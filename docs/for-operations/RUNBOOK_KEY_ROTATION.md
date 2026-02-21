# Runbook: Secret and Key Rotation

[← Back to SLA Overview](README.md)

---

## Overview
Procedures for rotating secrets, API keys, and service account keys used by the Precision Medicine MCP platform. Regular rotation reduces the risk of credential compromise and meets compliance requirements (HIPAA, SOC 2).

## 1. Rotation Schedule

| Secret Type | Rotation Frequency | Owner | Alert Threshold |
| :--- | :--- | :--- | :--- |
| **Service Account Keys** | Every 90 days | Infrastructure Admin | 14 days before expiry |
| **API Keys** (HuggingFace, Seqera) | Every 90 days | Infrastructure Admin | 14 days before expiry |
| **Database Credentials** | Every 90 days | Infrastructure Admin | 14 days before expiry |
| **TLS/SSL Certificates** | Every 365 days | Infrastructure Admin | 30 days before expiry |
| **MCP Server Tokens** | Every 90 days | Platform Editor | 14 days before expiry |

## 2. Rotation Workflow

### Step 1: Generate New Secret
```bash
# Example: Rotate a GCP service account key
gcloud iam service-accounts keys create new-key.json \
  --iam-account=mcp-service@PROJECT_ID.iam.gserviceaccount.com
```

### Step 2: Update Secret Store
Update the secret in your secret management system (e.g., GCP Secret Manager, Azure Key Vault):
```bash
# GCP Secret Manager example
gcloud secrets versions add MCP_SERVICE_KEY --data-file=new-key.json
```

### Step 3: Deploy Updated Configuration
1. Update `claude_desktop_config.json` or environment variables with the new secret.
2. Restart affected MCP servers to pick up the new credentials.
3. Verify connectivity for each rotated service.

### Step 4: Revoke Old Secret
```bash
# Revoke the previous key after verifying the new one works
gcloud iam service-accounts keys delete OLD_KEY_ID \
  --iam-account=mcp-service@PROJECT_ID.iam.gserviceaccount.com
```

### Step 5: Audit and Document
1. Record the rotation in the change management log.
2. Update the rotation tracking spreadsheet with new expiry dates.
3. Verify no alerts remain for the rotated credential.

## 3. MCP-Specific Secrets

| Secret | Used By | Config Location |
| :--- | :--- | :--- |
| `SEQERA_ACCESS_TOKEN` | mcp-seqera | `claude_desktop_config.json` → env |
| Service account keys | All servers (GCS access) | Environment or key file |
| TLS certificates | API Gateway | Cloud-managed |

## 4. Emergency Rotation (Suspected Compromise)

If a credential is suspected compromised:
1. **Immediately** generate a new secret (Step 1 above).
2. **Immediately** revoke the compromised secret (Step 4 above).
3. Review access logs for unauthorized usage during the exposure window.
4. File an incident report per [Runbook: P1 PHI Breach Response](RUNBOOK_P1_PHI_BREACH.md) if PHI access is suspected.

---

## Verification Checklist
- [ ] All secrets are within their rotation window (not expired).
- [ ] Old/revoked keys have been deleted from all systems.
- [ ] MCP servers restart successfully with new credentials.
- [ ] No authentication errors in server logs after rotation.

---

## Related Documents

- [Section 06: Security Operations](06_SECURITY_OPERATIONS.md)
- [Runbook: P1 PHI Breach Response](RUNBOOK_P1_PHI_BREACH.md)
- [Runbook: Monthly Maintenance](RUNBOOK_MONTHLY_MAINTENANCE.md)

---

### Document History

| Date | Version | Author | Change Summary |
|:---|:---|:---|:---|
| 2026-02-13 | 1.0 | IT Operations | Initial creation |
