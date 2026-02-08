# Appendix D: HIPAA Compliance Checklist

[← Back to SLA Overview](../README.md)

---

## Overview

This appendix provides comprehensive HIPAA compliance checklists for monthly, quarterly, and annual verification. Use these checklists to ensure continuous compliance with HIPAA Administrative, Technical, and Physical Safeguards.

---

## Monthly Compliance Checklist

**Review Frequency:** First Monday of each month  
**Owner:** Compliance Team + IT Operations  
**Deliverable:** Monthly compliance report to IT Director and CISO

### Administrative Safeguards (§164.308)

| Control | Verification Method | Status | Notes |
|---------|-------------------|--------|-------|
| ☐ **Audit log export to BigQuery** | Verify 100% success rate for past 30 days | ☐ Pass ☐ Fail | Check Cloud Logging export metrics |
| ☐ **MFA enrollment** | Verify 100% of users have MFA enabled | ☐ Pass ☐ Fail | Query Azure AD MFA status |
| ☐ **Service account key age** | All keys < 90 days old | ☐ Pass ☐ Fail | `gcloud iam service-accounts keys list` |
| ☐ **Failed login attempts** | Investigate if >50/user/day | ☐ Pass ☐ Fail | Review Cloud Logging for auth failures |
| ☐ **Unauthorized access attempts** | Review blocked requests, investigate anomalies | ☐ Pass ☐ Fail | Check Cloud Armor logs |

### Technical Safeguards (§164.312)

| Control | Verification Method | Status | Notes |
|---------|-------------------|--------|-------|
| ☐ **Encryption at rest** | Verify AES-256 on all Cloud Storage buckets | ☐ Pass ☐ Fail | `gsutil encryption ls gs://bucket` |
| ☐ **Encryption in transit** | Verify TLS 1.3 only (no TLS 1.2/1.1/1.0) | ☐ Pass ☐ Fail | SSL Labs scan results |
| ☐ **FHIR data integrity** | Run hash validation on random sample | ☐ Pass ☐ Fail | Healthcare API integrity checks |
| ☐ **Backup success rate** | Verify >95% success for daily backups | ☐ Pass ☐ Fail | Check backup job logs |
| ☐ **PHI access from unknown IP** | Review access logs, investigate new IPs | ☐ Pass ☐ Fail | Cloud Logging query |

### Data Retention (§164.316(b)(2))

| Control | Verification Method | Status | Notes |
|---------|-------------------|--------|-------|
| ☐ **Audit logs** | Verify 10-year retention in BigQuery | ☐ Pass ☐ Fail | Check BigQuery table expiration |
| ☐ **Patient data** | Verify 7-year retention, lifecycle policies | ☐ Pass ☐ Fail | Check Cloud Storage lifecycle rules |
| ☐ **Compliance documentation** | Verify 6-year retention | ☐ Pass ☐ Fail | Document management system audit |

### Automated Verification Script

```bash
#!/bin/bash
# Monthly HIPAA Compliance Check
# Run on first Monday of each month

echo "=== HIPAA Monthly Compliance Report ===" > hipaa_monthly_$(date +%Y%m).txt
echo "Generated: $(date)" >> hipaa_monthly_$(date +%Y%m).txt

# 1. Verify encryption at rest
echo "\n## Encryption at Rest" >> hipaa_monthly_$(date +%Y%m).txt
gsutil encryption ls gs://patient-data-uploads/ >> hipaa_monthly_$(date +%Y%m).txt
# Expected: Encryption: Customer-managed OR Encryption: Google-managed

# 2. Verify TLS 1.3 enforcement
echo "\n## TLS Configuration" >> hipaa_monthly_$(date +%Y%m).txt
curl -I https://mcp-fgbio-*.run.app 2>&1 | grep -i 'tls' >> hipaa_monthly_$(date +%Y%m).txt
# Expected: TLSv1.3

# 3. Verify audit log retention
echo "\n## Audit Log Retention" >> hipaa_monthly_$(date +%Y%m).txt
bq show --format=prettyjson PROJECT:audit_logs.cloudaudit_googleapis_com_activity \
  | jq '.expirationTime' >> hipaa_monthly_$(date +%Y%m).txt
# Expected: Date 10 years from now

# 4. Verify MFA enrollment (Azure AD)
echo "\n## MFA Enrollment" >> hipaa_monthly_$(date +%Y%m).txt
# Requires Azure AD admin token
az ad user list --query "[?otherMails[0]=='@hospital.org'].{UPN:userPrincipalName,MFA:strongAuthenticationMethods}" \
  >> hipaa_monthly_$(date +%Y%m).txt
# Expected: 100% enrollment

# 5. Verify service account key age
echo "\n## Service Account Key Age" >> hipaa_monthly_$(date +%Y%m).txt
for SA in mcp-fgbio mcp-multiomics mcp-spatialtools; do
  gcloud iam service-accounts keys list \
    --iam-account=${SA}@PROJECT.iam.gserviceaccount.com \
    --format='table(name,validAfterTime)' >> hipaa_monthly_$(date +%Y%m).txt
done
# Expected: All keys < 90 days old

# 6. Check backup success rate
echo "\n## Backup Success Rate (Last 30 Days)" >> hipaa_monthly_$(date +%Y%m).txt
gcloud logging read \
  'resource.type="cloud_storage_bucket" AND
   protoPayload.methodName="storage.objects.create" AND
   protoPayload.resourceName=~"backups-fhir"
  ' --limit=1000 --format=json \
  | jq '[.[] | select(.severity != "ERROR")] | length' >> hipaa_monthly_$(date +%Y%m).txt
# Expected: >95% success

# Email report
cat hipaa_monthly_$(date +%Y%m).txt | mail -s "Monthly HIPAA Compliance Report" compliance@hospital.org
```

---

## Quarterly Compliance Checklist

**Review Frequency:** End of each quarter (March 31, June 30, Sept 30, Dec 31)  
**Owner:** Compliance Team + CISO + IT Director  
**Deliverable:** Quarterly compliance report to Executive Team

### Security Testing (§164.308(a)(8))

| Control | Verification Method | Status | Notes |
|---------|-------------------|--------|-------|
| ☐ **Vulnerability scan** | Run automated scan on all infrastructure | ☐ Pass ☐ Fail | Use GCP Security Command Center |
| ☐ **Security configuration review** | Audit Cloud Run, IAM, network configs | ☐ Pass ☐ Fail | Solutions Architect review |
| ☐ **Access control review** | Remove unused roles, verify least privilege | ☐ Pass ☐ Fail | IAM role audit |
| ☐ **CMEK key rotation** | Rotate customer-managed encryption keys | ☐ Pass ☐ Fail | Cloud KMS key rotation |

### Disaster Recovery Testing (§164.308(a)(7))

| Control | Verification Method | Status | Notes |
|---------|-------------------|--------|-------|
| ☐ **DR drill** | Execute DR failover to us-east1 | ☐ Pass ☐ Fail | See [DR Failover Runbook](../runbooks/DR_FAILOVER.md) |
| ☐ **Backup restoration test** | Restore random FHIR snapshot, verify integrity | ☐ Pass ☐ Fail | Healthcare API restore test |
| ☐ **RTO/RPO verification** | Measure actual recovery time vs target | ☐ Pass ☐ Fail | RTO target: 1-4 hours |

### Business Associate Agreements (§164.308(b)(1))

| Vendor | BAA Status | Renewal Date | Contact |
|--------|-----------|--------------|---------|
| ☐ **Google Cloud Platform (GCP)** | ☐ Signed ☐ Expired | [Date] | GCP Account Manager |
| ☐ **Anthropic (Claude API)** | ☐ Signed ☐ Expired ☐ N/A | [Date] | enterprise-support@anthropic.com |
| ☐ **Microsoft (Azure AD)** | ☐ Signed ☐ Expired | [Date] | Microsoft Account Manager |

### Training & Awareness (§164.308(a)(5))

| Control | Verification Method | Status | Notes |
|---------|-------------------|--------|-------|
| ☐ **HIPAA training completion** | Verify 100% of staff completed annual training | ☐ Pass ☐ Fail | LMS completion report |
| ☐ **Security awareness training** | Phishing simulation, security best practices | ☐ Pass ☐ Fail | Training platform metrics |
| ☐ **Incident response drills** | Tabletop exercise for PHI breach | ☐ Pass ☐ Fail | Post-drill debrief |

---

## Annual Compliance Checklist

**Review Frequency:** End of fiscal year  
**Owner:** External Auditor + Compliance Team + CISO  
**Deliverable:** Annual HIPAA audit report to Board of Directors

### Security Risk Assessment (§164.308(a)(1)(ii)(A))

| Assessment Area | Method | Status | Findings |
|-----------------|--------|--------|----------|
| ☐ **Threat identification** | Review current threat landscape | ☐ Complete | Document new threats |
| ☐ **Vulnerability assessment** | Comprehensive infrastructure scan | ☐ Complete | Remediation plan for findings |
| ☐ **Risk analysis** | NIST SP 800-30 risk assessment | ☐ Complete | Risk register updated |
| ☐ **Risk mitigation** | Review controls, identify gaps | ☐ Complete | Corrective action plan |

### External Audit & Penetration Testing

| Test Type | Scope | Status | Report Date |
|-----------|-------|--------|-------------|
| ☐ **External HIPAA audit** | All HIPAA controls | ☐ Pass ☐ Fail | [Date] |
| ☐ **Penetration testing** | External-facing APIs, Cloud Run endpoints | ☐ Pass ☐ Fail | [Date] |
| ☐ **Social engineering test** | Phishing simulation (staff) | ☐ Pass ☐ Fail | [Date] |

### Certification Renewals

| Certification | Status | Expiration Date | Owner |
|---------------|--------|-----------------|-------|
| ☐ **HITRUST CSF** (GCP) | ☐ Current ☐ Expired | [Date] | GCP |
| ☐ **SOC 2 Type II** (GCP) | ☐ Current ☐ Expired | [Date] | GCP |
| ☐ **ISO 27001** (GCP) | ☐ Current ☐ Expired | [Date] | GCP |
| ☐ **Cyber Liability Insurance** | ☐ Current ☐ Expired | [Date] | Finance Team |

### Policy Review & Updates (§164.316(b)(2)(iii))

| Policy | Last Updated | Next Review | Owner |
|--------|--------------|-------------|-------|
| ☐ **HIPAA Privacy Policy** | [Date] | [Date + 1 year] | Compliance Officer |
| ☐ **HIPAA Security Policy** | [Date] | [Date + 1 year] | CISO |
| ☐ **Incident Response Policy** | [Date] | [Date + 1 year] | IT Director |
| ☐ **Breach Notification Policy** | [Date] | [Date + 1 year] | Compliance Officer |
| ☐ **Data Retention Policy** | [Date] | [Date + 1 year] | Compliance Officer |

---

## HIPAA Safeguards Summary

### Administrative Safeguards (§164.308)

| Safeguard | Implementation | Status |
|-----------|----------------|--------|
| **Security Management Process** | Annual risk assessment, documented policies | ✅ Implemented |
| **Assigned Security Responsibility** | Dedicated HIPAA Security Officer (CISO) | ✅ Implemented |
| **Workforce Security** | Authorization, supervision, termination procedures | ✅ Implemented |
| **Information Access Management** | Role-based access (Azure AD + GCP IAM) | ✅ Implemented |
| **Security Awareness Training** | Annual HIPAA training, phishing simulations | ✅ Implemented |
| **Security Incident Procedures** | Documented incident response plan, runbooks | ✅ Implemented |
| **Contingency Plan** | DR plan, backup/restore, quarterly testing | ✅ Implemented |
| **Evaluation** | Quarterly compliance reviews | ✅ Implemented |
| **Business Associate Contracts** | BAA with GCP, Anthropic, Microsoft | ✅ Implemented |

### Technical Safeguards (§164.312)

| Safeguard | Implementation | Status |
|-----------|----------------|--------|
| **Access Control** | Azure AD SSO, MFA mandatory, IAM roles | ✅ Implemented |
| **Audit Controls** | Cloud Logging (10-year retention, immutable) | ✅ Implemented |
| **Integrity** | FHIR data checksums, encryption | ✅ Implemented |
| **Person/Entity Authentication** | MFA via Azure AD (100% enrollment) | ✅ Implemented |
| **Transmission Security** | TLS 1.3 only (no TLS 1.2/1.1/1.0) | ✅ Implemented |

### Physical Safeguards (§164.310)

| Safeguard | Implementation | Status |
|-----------|----------------|--------|
| **Facility Access Controls** | GCP data center security (badge, biometrics) | ✅ GCP Responsibility |
| **Workstation Use** | Hospital-managed workstations, endpoint protection | ✅ Hospital IT |
| **Workstation Security** | Screen locks, full disk encryption | ✅ Hospital IT |
| **Device and Media Controls** | No removable media access to cloud systems | ✅ Implemented |

---

## Breach Notification Requirements

### Breach Discovery Timeline

| Timeframe | Action | Owner |
|-----------|--------|-------|
| **Day 0** | Breach discovered | IT Operations / Security Team |
| **< 1 hour** | Notify CISO (internal escalation) | On-Call Engineer |
| **< 4 hours** | Initial scope assessment, affected patient count | Security Team |
| **< 24 hours** | Forensic preservation complete | IT Security |
| **< 72 hours** | Post-incident review, root cause analysis | IT Director |
| **< 5 business days** | Corrective action plan | Change Advisory Board |

### Breach Notification Timeline (Per HIPAA)

| Breach Scope | Notification Requirement | Deadline |
|--------------|-------------------------|----------|
| **< 500 individuals** | Notify affected individuals | ≤ 60 days from discovery |
| **≥ 500 individuals** | Notify affected individuals + HHS/OCR + media | ≤ 60 days from discovery |
| **Any breach** | Document in breach log | Immediately |

### Breach Notification Content (Required by HIPAA)

Notification to affected individuals must include:

1. **What happened:** Brief description of the breach
2. **Date of breach:** When it occurred (or estimated if unknown)
3. **PHI involved:** Types of information compromised (e.g., names, SSNs, medical records)
4. **What we're doing:** Steps taken to investigate and mitigate
5. **What you can do:** Steps individuals can take to protect themselves
6. **Contact information:** Phone number and email for questions

**Example Notification Letter Template:** [Available in HIPAA policy folder]

---

## Sanctions Policy (§164.308(a)(1)(ii)(C))

### Workforce Sanctions for HIPAA Violations

| Violation | First Offense | Second Offense | Third Offense |
|-----------|---------------|----------------|---------------|
| **Unintentional PHI disclosure** (e.g., emailing unencrypted data) | Written warning + retraining | Suspension | Termination |
| **Accessing PHI without authorization** (curiosity/snooping) | Suspension | Termination | Termination + legal action |
| **Intentional PHI breach** (malicious) | Immediate termination + legal action | N/A | N/A |
| **Failure to report suspected breach** | Written warning | Suspension | Termination |
| **Sharing login credentials** | Written warning + password reset | Suspension | Termination |

**Documentation:** All sanctions must be documented in employee file and reported to Compliance Officer.

---

## Compliance Reporting

### Monthly Report Template

**To:** IT Director, CISO  
**From:** Compliance Team  
**Date:** [First Monday of month]  
**Subject:** Monthly HIPAA Compliance Report - [Month Year]

**Summary:**
- ☐ All controls passed
- ☐ [X] controls failed (see details below)

**Administrative Safeguards:**
- Audit log export: [Pass/Fail] - [Details]
- MFA enrollment: [Pass/Fail] - [100%/XX%]
- Service account key rotation: [Pass/Fail] - [All <90 days / XX keys expired]

**Technical Safeguards:**
- Encryption at rest: [Pass/Fail]
- Encryption in transit: [Pass/Fail]
- Backup success rate: [XX%]

**Findings:**
1. [Issue description]
   - Severity: [Low/Medium/High/Critical]
   - Remediation: [Action plan]
   - Owner: [Name]
   - Due Date: [Date]

**Next Steps:**
- [Action items for next month]

---

## Related Documents

- [Section 4: Compliance & Regulatory](../sections/04_COMPLIANCE_REGULATORY.md) - HIPAA safeguards detail
- [Section 5: Data Protection](../sections/05_DATA_PROTECTION.md) - Encryption standards
- [Section 6: Security Operations](../sections/06_SECURITY_OPERATIONS.md) - Security testing schedule
- [Appendix A: Incident Classification Matrix](APPENDIX_A_INCIDENT_CLASSIFICATION.md) - PHI breach response

[← Back to SLA Overview](../README.md)
