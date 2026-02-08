# Appendix E: Contact and Escalation Directory

[← Back to SLA Overview](../README.md)

---

## Overview

This appendix provides complete contact information for all IT personnel, escalation paths, vendor contacts, and emergency hotlines supporting the Precision Medicine MCP Platform.

**Last Updated:** [Date]  
**Owner:** IT Operations Team  
**Review Frequency:** Quarterly (or immediately upon personnel changes)

---

## Primary IT Contacts

### Executive Leadership

| Role | Name | Email | Phone | PagerDuty | Availability |
|------|------|-------|-------|-----------|--------------|
| **Chief Technology Officer (CTO)** | [Name] | [email@hospital.org] | +1-XXX-XXX-XXXX | N/A | Business hours |
| **Chief Information Officer (CIO)** | [Name] | [email@hospital.org] | +1-XXX-XXX-XXXX | N/A | Business hours |
| **Chief Information Security Officer (CISO)** | [Name] | [email@hospital.org] | +1-XXX-XXX-XXXX (24/7 for PHI breaches) | Security escalation | 24/7 for P0 security |
| **Chief Medical Information Officer (CMIO)** | [Name] | [email@hospital.org] | +1-XXX-XXX-XXXX | N/A | Business hours |

### IT Management

| Role | Name | Email | Phone | PagerDuty | Availability |
|------|------|-------|-------|-----------|--------------|
| **IT Director** | [Name] | [email@hospital.org] | +1-XXX-XXX-XXXX | Tertiary escalation | Business hours + P0 escalation |
| **IT Manager (Operations)** | [Name] | [email@hospital.org] | +1-XXX-XXX-XXXX | Secondary escalation | Business hours + P0/P1 escalation |
| **Solutions Architect** | [Name] | [email@hospital.org] | +1-XXX-XXX-XXXX | N/A | Business hours |
| **Security Manager** | [Name] | [email@hospital.org] | +1-XXX-XXX-XXXX | Security escalation | 24/7 for security incidents |

### IT Operations Team

| Role | Name | Email | Phone | PagerDuty | Availability |
|------|------|-------|-------|-----------|--------------|
| **On-Call Engineer (Primary)** | [Rotation] | it-oncall@hospital.org | [PagerDuty phone] | Primary on-call | 24/7 (rotating) |
| **On-Call Engineer (Backup)** | [Rotation] | it-oncall-backup@hospital.org | [PagerDuty phone] | Backup on-call | 24/7 (rotating) |
| **SRE Team Lead** | [Name] | [email@hospital.org] | +1-XXX-XXX-XXXX | N/A | Business hours |
| **DevOps Engineer** | [Name] | [email@hospital.org] | +1-XXX-XXX-XXXX | N/A | Business hours |

### Specialized Teams

| Role | Name | Email | Phone | Availability |
|------|------|-------|-------|--------------|
| **FinOps Lead** | [Name] | [email@hospital.org] | +1-XXX-XXX-XXXX | Business hours |
| **Compliance Officer** | [Name] | [email@hospital.org] | +1-XXX-XXX-XXXX | Business hours |
| **Data Privacy Officer** | [Name] | [email@hospital.org] | +1-XXX-XXX-XXXX | Business hours |
| **Legal Counsel** | [Name] | [email@hospital.org] | +1-XXX-XXX-XXXX | Business hours |

---

## On-Call Rotation Schedule

### Current On-Call (Week of [Date])

| Shift | Engineer | Contact | Backup |
|-------|----------|---------|--------|
| **Primary** | [Name] | [PagerDuty] / [Phone] | [Backup Name] |
| **Secondary** | [Name] | [PagerDuty] / [Phone] | [Backup Name] |

**PagerDuty Schedule URL:** [https://hospital.pagerduty.com/schedules/...]

### Shift Handoff Procedures

**Handoff Time:** Every Monday 9:00 AM  
**Method:** Slack #precision-medicine-oncall + Email

**Handoff Checklist:**
- [ ] Review open P0/P1 incidents from previous week
- [ ] Review upcoming maintenance windows
- [ ] Verify access to all critical systems (GCP, Anthropic, Azure AD)
- [ ] Confirm PagerDuty is receiving test alerts
- [ ] Update on-call calendar

---

## Escalation Matrix

### By Incident Severity

| Priority | L1 (On-Call) | L2 (IT Manager) | L3 (IT Director / CISO) | L4 (Executive) |
|----------|--------------|-----------------|------------------------|----------------|
| **P0 - Critical** | Immediate page | 15 min | 30 min | 1 hour |
| **P1 - High** | Immediate page | 1 hour | 4 hours (if PHI breach) | As needed |
| **P2 - Medium** | Business hours | 4 hours | As needed | Not required |
| **P3 - Low** | Ticket queue | As needed | Not required | Not required |

### By Incident Type

| Incident Type | Primary Contact | Escalation Path | Special Procedures |
|---------------|----------------|-----------------|-------------------|
| **Complete System Outage** | On-Call Engineer → IT Manager → IT Director → CTO | Standard P0 escalation | Activate incident command |
| **PHI Breach** | On-Call Engineer → CISO (within 1 hour) → Legal → Compliance | Security escalation path | HIPAA notification timeline |
| **Anthropic API Outage** | On-Call Engineer → IT Manager | No further escalation (wait for vendor) | Monitor status.anthropic.com |
| **GCP Regional Failure** | On-Call Engineer → IT Manager → IT Director | Initiate DR failover if >4 hour ETA | Execute DR runbook |
| **Audit Log Failure** | On-Call Engineer → Compliance Officer + CISO | Critical HIPAA violation | Fix within 24 hours |

---

## Emergency Hotlines (24/7)

| Purpose | Phone Number | Email | Notes |
|---------|-------------|-------|-------|
| **IT On-Call (Primary)** | [PagerDuty: +1-XXX-XXX-XXXX] | it-oncall@hospital.org | Auto-routes to current on-call |
| **HIPAA Breach Hotline** | [+1-XXX-XXX-XXXX] | hipaa-breach@hospital.org | Direct to CISO (24/7) |
| **Security Incident Hotline** | [+1-XXX-XXX-XXXX] | security-incident@hospital.org | Security Operations Center |
| **Hospital Emergency Command** | [+1-XXX-XXX-XXXX] | emergency-command@hospital.org | For P0 affecting patient care |

---

## Vendor Contacts & Escalation

### Google Cloud Platform (GCP)

**Support Plan:** Premium Support (24/7)  
**Account Manager:** [Name], [email@google.com], [Phone]

| Severity | Contact Method | Response Time | Escalation |
|----------|---------------|---------------|------------|
| **P1 (Critical)** | [GCP Support Console](https://console.cloud.google.com/support) | ≤ 1 hour | Automatic to TAM after 2 hours |
| **P2 (High)** | GCP Support Console | ≤ 4 hours | Escalate to Account Manager if needed |
| **P3 (Normal)** | GCP Support Console | ≤ 8 hours | Standard queue |
| **P4 (Low)** | GCP Support Console | ≤ 24 hours | Standard queue |

**Useful Links:**
- **GCP Status Dashboard:** https://status.cloud.google.com
- **Support Console:** https://console.cloud.google.com/support
- **TAM (Technical Account Manager):** [Name], [email@google.com]

**Escalation Path:**
1. Open support case via console
2. Call support hotline (if P1): +1-XXX-XXX-XXXX
3. Email TAM (if no response after 2 hours): [email@google.com]
4. Contact Account Manager (if critical): [email@google.com]

---

### Anthropic (Claude API)

**Support Plan:** Enterprise (email only, no phone support)  
**Contact:** enterprise-support@anthropic.com

| Severity | Contact Method | Response Time | Notes |
|----------|---------------|---------------|-------|
| **Critical Outage** | Email: enterprise-support@anthropic.com | Best-effort (no SLA) | ⚠️ No guaranteed response time |
| **Performance Issues** | Email: enterprise-support@anthropic.com | Best-effort | Monitor status.anthropic.com |
| **Billing Questions** | Email: enterprise-support@anthropic.com | 1-2 business days | |

**Useful Links:**
- **Status Page:** https://status.anthropic.com
- **Documentation:** https://docs.anthropic.com

**⚠️ Critical Note:**  
Anthropic Claude API is a **single point of failure** for this system. If api.anthropic.com is unreachable for >15 minutes, escalate to IT Manager immediately. **However, there is no phone support or escalation path beyond email.**

**Workaround During Outage:**
- Monitor https://status.anthropic.com for ETA
- Notify users of external dependency outage
- Wait for Anthropic to resolve (no alternative available)

---

### Microsoft (Azure AD)

**Support Plan:** Premier Support  
**Account Manager:** [Name], [email@microsoft.com], [Phone]

| Severity | Contact Method | Response Time | Escalation |
|----------|---------------|---------------|------------|
| **Severity A (Critical)** | Phone: 1-800-MICROSOFT | ≤ 1 hour | Automatic to TAM after 2 hours |
| **Severity B (High)** | Phone or Portal | ≤ 2 hours | Escalate to Account Manager |
| **Severity C (Normal)** | Portal | ≤ 8 hours | Standard queue |

**Useful Links:**
- **Azure Status Dashboard:** https://status.azure.com
- **Support Portal:** https://portal.azure.com/#blade/Microsoft_Azure_Support/HelpAndSupportBlade

**Escalation Path:**
1. Call 1-800-MICROSOFT (for Severity A/B)
2. Reference contract number: [CONTRACT-NUMBER]
3. Email TAM (if no response after 2 hours): [email@microsoft.com]

---

### Epic Systems (Hospital FHIR Server)

**Support:** Per hospital IT contract  
**Contact:** [Hospital Epic Team]

| Contact | Phone | Email | Availability |
|---------|-------|-------|--------------|
| **Epic IT Liaison** | [+1-XXX-XXX-XXXX] | [email@hospital.org] | Business hours |
| **Epic On-Call** | [+1-XXX-XXX-XXXX] | N/A | 24/7 (hospital IT manages) |

**Escalation:**
- Contact hospital Epic team first (not Epic Systems directly)
- Epic Systems support is managed by hospital IT contract

---

## Distribution Lists

### Operational

| List | Email | Purpose | Subscribers |
|------|-------|---------|-------------|
| **IT On-Call** | it-oncall@hospital.org | On-call engineer rotation | Current on-call engineer |
| **IT Operations Team** | it-operations@hospital.org | All IT operations staff | 12 engineers |
| **IT Leadership** | it-leadership@hospital.org | IT management | IT Director, IT Managers, CISO |
| **SRE Team** | sre-team@hospital.org | Site reliability engineers | 5 SREs |

### Incident Alerts

| List | Email | Purpose | Alert Severity |
|------|-------|---------|----------------|
| **Clinical Alerts** | clinical-alerts@hospital.org | Critical incidents affecting patient care | P0 only |
| **IT Alerts** | it-alerts@hospital.org | All incidents | P0, P1 |
| **Security Alerts** | security-alerts@hospital.org | Security incidents | P0 (security), P1 (security) |

### Compliance & Governance

| List | Email | Purpose |
|------|-------|---------|
| **Compliance Team** | compliance@hospital.org | HIPAA compliance, audits |
| **Change Advisory Board** | change-advisory-board@hospital.org | SLA changes, major updates |
| **Executive Team** | executive-team@hospital.org | C-suite (CTO, CIO, CISO, CMIO) |

---

## Slack Channels

| Channel | Purpose | Severity |
|---------|---------|----------|
| **#precision-medicine-p0-incidents** | P0 critical incidents | P0 only |
| **#precision-medicine-alerts** | All incidents and alerts | P0, P1, P2 |
| **#precision-medicine-oncall** | On-call coordination | All |
| **#precision-medicine-deployments** | Deployment notifications | Informational |
| **#precision-medicine-general** | General discussion | Informational |

---

## Conference Bridge / War Room

### Incident Command Conference Bridge

**Purpose:** P0 incidents requiring coordination across multiple teams

**Dial-In:**
- **Phone:** [+1-XXX-XXX-XXXX]
- **PIN:** [XXXX]
- **Zoom:** [https://hospital.zoom.us/j/XXXXXXXXX]

**When to Activate:**
- P0 incidents lasting >1 hour
- PHI breach incidents (any severity)
- Multi-region failures
- Incidents affecting >100 concurrent users

**Roles:**
- **Incident Commander:** IT Manager (or IT Director for P0 >2 hours)
- **Scribe:** IT Operations engineer (documents timeline)
- **Communications Lead:** IT Manager (stakeholder updates)
- **Technical Lead:** On-Call Engineer (executes fixes)

---

## External Auditor & Legal Contacts

### HIPAA Compliance Auditor

| Firm | Contact | Phone | Email |
|------|---------|-------|-------|
| **[Audit Firm Name]** | [Contact Name] | [+1-XXX-XXX-XXXX] | [email@auditfirm.com] |

**Annual Audit:** [Month] each year  
**Scope:** HIPAA compliance, SOC 2 Type II

### Legal Counsel

| Role | Name | Phone | Email |
|------|------|-------|-------|
| **General Counsel** | [Name] | [+1-XXX-XXX-XXXX] | [email@hospital.org] |
| **Privacy Counsel** | [Name] | [+1-XXX-XXX-XXXX] | [email@hospital.org] |

**When to Contact:**
- PHI breach affecting ≥500 patients (within 24 hours)
- Legal questions on HIPAA notification requirements
- Contract disputes with vendors

---

## After-Hours & Holiday Contacts

### Holiday On-Call Schedule

For major holidays (Thanksgiving, Christmas, New Year), verify on-call coverage 2 weeks in advance.

**Current Holiday Schedule:** [Link to shared calendar]

### Weekend Coverage

- **Saturdays:** Primary on-call engineer
- **Sundays:** Primary on-call engineer
- **Escalation:** Same as weekday (IT Manager → IT Director)

**Note:** Response and resolution times are the same 24/7/365. No exceptions for weekends or holidays.

---

## Contact Update Procedures

### When to Update This Directory

1. **Immediately:** Personnel changes (hire, termination, role change)
2. **Quarterly:** Vendor contact verification
3. **Annually:** Full directory review

### How to Update

1. Submit pull request to update this file
2. Get approval from IT Director
3. Merge to main branch
4. Notify stakeholders via email + Slack

### Distribution After Updates

- [ ] Email IT Operations Team
- [ ] Post to Slack #precision-medicine-general
- [ ] Update printed on-call binder (if maintained)

---

## Quick Reference Card (Print & Laminate for On-Call Binder)

```
┌─────────────────────────────────────────────────────────┐
│  PRECISION MEDICINE MCP - EMERGENCY CONTACTS            │
├─────────────────────────────────────────────────────────┤
│                                                         │
│  IT ON-CALL (24/7):         [PagerDuty: XXX-XXX-XXXX]  │
│  HIPAA BREACH HOTLINE:      [Direct CISO: XXX-XXX-XXXX]│
│  SECURITY INCIDENT:         [SOC: XXX-XXX-XXXX]        │
│                                                         │
│  ESCALATION:                                            │
│  - P0: On-Call → IT Mgr (15m) → IT Dir (30m)          │
│  - PHI Breach: On-Call → CISO (1h) → Legal            │
│                                                         │
│  VENDOR SUPPORT:                                        │
│  - GCP: console.cloud.google.com/support               │
│  - Anthropic: enterprise-support@anthropic.com         │
│  - Azure AD: 1-800-MICROSOFT                           │
│                                                         │
│  RUNBOOKS: /docs/operations/sla/runbooks/              │
└─────────────────────────────────────────────────────────┘
```

---

## Related Documents

- [Section 3: Incident Management](../sections/03_INCIDENT_MANAGEMENT.md) - Escalation procedures
- [Appendix A: Incident Classification Matrix](APPENDIX_A_INCIDENT_CLASSIFICATION.md) - Severity definitions
- [Runbook: P0 Complete System Outage](../runbooks/P0_COMPLETE_OUTAGE.md)
- [Runbook: P1 PHI Breach Response](../runbooks/P1_PHI_BREACH.md)

[← Back to SLA Overview](../README.md)
