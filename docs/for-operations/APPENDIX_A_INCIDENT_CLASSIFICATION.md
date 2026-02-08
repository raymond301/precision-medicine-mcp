# Appendix A: Incident Classification Matrix

[← Back to SLA Overview](../README.md)

---

## Incident Priority Definitions

This appendix provides detailed classification criteria for all incident severity levels (P0-P3), including response times, resolution times, escalation procedures, and examples.

---

## Priority Level Summary Table

| Priority | Response Time | Resolution Time | Escalation | Notification |
|----------|---------------|-----------------|------------|--------------|
| **P0 - Critical** | ≤ 15 minutes | ≤ 4 hours | Immediate to IT Director | PagerDuty + Email + Slack + SMS (executive) |
| **P1 - High** | ≤ 1 hour | ≤ 8 hours | IT Manager at 1 hour | PagerDuty + Email + Slack |
| **P2 - Medium** | ≤ 4 hours | ≤ 3 business days | Team lead at 4 hours | Email + Slack |
| **P3 - Low** | ≤ 24 hours | ≤ 7 business days | Standard queue | Ticket system only |

---

## P0 - Critical Incidents

### Definition
**System-wide outage or security breach affecting patient care or data integrity.**

### Criteria (Any ONE of these triggers P0)

✅ **Infrastructure Failures:**
- All 9 MCP servers unavailable
- Complete GCP regional outage (us-central1)
- Anthropic Claude API down for >15 minutes
- Healthcare API FHIR store unavailable (no read/write access)
- Complete network outage (load balancer failure)

✅ **Security Incidents:**
- **PHI breach:** Unauthorized access to patient data
- **Ransomware attack:** Encryption of critical systems
- **Data exfiltration:** Evidence of data being sent to external IPs
- **Audit log failure:** Cloud Logging export to BigQuery failing (HIPAA violation risk)

✅ **Data Loss:**
- FHIR store corruption (data integrity compromised)
- Backup failure detected during restore attempt
- Accidental deletion of production data

✅ **Service Degradation (Critical Tier 1):**
- >50% of Tier 1 servers down (2+ out of 3 clinical servers)
- Error rate >25% sustained for >15 minutes
- Response latency >10 seconds sustained for >15 minutes

### Examples

| Scenario | Impact | Response |
|----------|--------|----------|
| **All servers return 502 Bad Gateway** | Zero analysis capacity, patient care delayed | Check Anthropic API status, initiate DR failover |
| **Unauthorized user accessed 1,000+ patient FHIR records** | PHI breach, HIPAA violation | Block user, notify CISO within 1 hour, forensic preservation |
| **Cloud Logging export to BigQuery failing for 6 hours** | HIPAA audit trail gap | Critical - Fix immediately, notify Compliance Team |
| **GCP announces us-central1 outage with >4 hour ETA** | Regional failure | Execute DR failover to us-east1 |
| **Ransomware detected on mcp-fgbio server** | Potential data encryption/loss | Isolate server, block network, notify Security Team |

### Response Timeline

| Time | Action | Owner |
|------|--------|-------|
| **0-5 min** | Alert triggers, PagerDuty pages on-call engineer | Automated |
| **5-15 min** | Acknowledge incident, initial assessment | On-Call Engineer |
| **15-30 min** | IT Manager notified if not acknowledged | Automated escalation |
| **30 min** | CISO and IT Director notified | IT Manager |
| **1 hour** | Executive team notified (CTO, CMIO) | IT Director |
| **Every 30 min** | Status update to stakeholders | Incident Commander |
| **≤ 4 hours** | Service restored or workaround in place | Incident Response Team |
| **< 24 hours** | Detailed incident report | On-Call Engineer |
| **< 72 hours** | Post-mortem meeting | IT Manager |
| **< 5 business days** | Corrective action plan | IT Director |

### Escalation Path

```
On-Call Engineer (0 min)
    ↓ (15 min if not acknowledged)
IT Manager (15 min)
    ↓ (30 min if not resolved)
IT Director + CISO (30 min)
    ↓ (1 hour if not resolved)
Executive Team: CTO, CMIO (1 hour)
    ↓ (2 hours if not resolved)
Hospital Incident Command (2 hours)
```

### Communication Channels

- **PagerDuty:** Immediate page (phone + SMS + app)
- **Email:** clinical-alerts@hospital.org (immediate)
- **Slack:** #precision-medicine-p0-incidents (immediate)
- **SMS:** Executive team (within 1 hour)
- **Status Page:** Update every 30 minutes until resolved
- **Phone:** Conference bridge activated for incident command

### Special Procedures for PHI Breach (P0)

**HIPAA Timeline Requirements:**

| Timeframe | Action | Owner |
|-----------|--------|-------|
| **≤ 1 hour** | Notify CISO (internal escalation) | On-Call Engineer |
| **≤ 4 hours** | Initial scope assessment, affected patient count | Security Team |
| **≤ 24 hours** | Forensic preservation complete | IT Security |
| **≤ 72 hours** | Post-incident review, root cause analysis | IT Director |
| **≤ 5 business days** | Corrective action plan | Change Advisory Board |
| **≤ 60 days** | Notify affected individuals (if ≥500 patients) | Compliance Team + Legal |
| **≤ 60 days** | Notify HHS/OCR (if ≥500 patients) | Compliance Team |

**Financial Penalties (Internal SLA):**
- **< 500 patients:** $1,000 per patient
- **≥ 500 patients:** $50,000 base + $500 per patient
- **Audit log failure:** $5,000 per day
- **Missed notification deadline:** $10,000 per day

---

## P1 - High Priority Incidents

### Definition
**Significant service degradation affecting clinical workflows or multiple users.**

### Criteria (Any ONE of these triggers P1)

✅ **Partial Infrastructure Failures:**
- Single Tier 1 server down (mcp-fgbio, mcp-multiomics, mcp-mockepic)
- Multiple Tier 2 servers down (2+ out of 5)
- FHIR API slow (>5 second response time sustained)
- Cloud Storage unavailable (can't upload/download patient data)

✅ **Performance Degradation:**
- Error rate >5% sustained for >1 hour
- Response latency >5 seconds (p95) sustained for >1 hour
- Memory/CPU exhaustion causing throttling
- Network bandwidth saturation

✅ **Security Issues (Non-Breach):**
- Suspected unauthorized access attempt (blocked by firewall)
- MFA failures >50 per user per day
- Service account key rotation failed
- Encryption key near expiration (< 7 days)

✅ **Data Issues:**
- Backup failure (single snapshot failed)
- Data sync lag >4 hours (DR region out of sync)
- FHIR data integrity check failed (checksum mismatch)

### Examples

| Scenario | Impact | Response |
|----------|--------|----------|
| **mcp-fgbio server returning 500 errors** | Genomic analysis unavailable | Check logs, rollback recent deployment, scale up resources |
| **FHIR API response time increased from 500ms to 8s** | Slow clinical data access | Check Healthcare API quotas, optimize queries, contact GCP support |
| **3 concurrent users report "can't log in"** | Multiple user authentication failure | Check Azure AD status, verify IAM roles, review MFA issues |
| **Daily FHIR backup failed overnight** | Backup gap, potential data loss risk | Re-run backup job, verify storage permissions, alert if fails again |
| **50+ failed login attempts for admin account** | Potential brute-force attack | Lock account, notify Security Team, investigate source IPs |

### Response Timeline

| Time | Action | Owner |
|------|--------|-------|
| **0-30 min** | Alert triggers, ticket created | Automated |
| **30-60 min** | Acknowledge incident, initial triage | On-Call Engineer |
| **1 hour** | IT Manager notified if not acknowledged | Automated escalation |
| **4 hours** | IT Director notified if not resolved (50% of 8-hour window) | IT Manager |
| **Every 2 hours** | Status update in Slack | On-Call Engineer |
| **≤ 8 hours** | Service restored or workaround in place | On-Call Engineer + IT Manager |
| **< 48 hours** | Incident report documented | On-Call Engineer |

### Escalation Path

```
On-Call Engineer (0 min)
    ↓ (1 hour if not acknowledged)
IT Manager (1 hour)
    ↓ (4 hours if not resolved - 50% of resolution window)
IT Director (4 hours, if needed)
```

### Communication Channels

- **PagerDuty:** Immediate page
- **Email:** it-alerts@hospital.org
- **Slack:** #precision-medicine-alerts
- **Status Page:** Update every 2 hours

---

## P2 - Medium Priority Incidents

### Definition
**Single-user issues or non-critical service degradation.**

### Criteria (Any ONE of these triggers P2)

✅ **Individual User Issues:**
- Single user can't access system (authentication issue)
- Single user experiencing slow performance
- Single user data upload failing

✅ **Non-Critical Service Issues:**
- Tier 2 server intermittently slow (not down)
- Tier 3 server down (mcp-tcga, mcp-deepcell, etc.)
- Non-critical feature broken (e.g., UI cosmetic bug)
- Documentation error or outdated runbook

✅ **Minor Performance Issues:**
- Error rate 2-5% (above 1% target, below 5% critical)
- Response latency 3-5 seconds (above 2s target, below 5s critical)
- Cold start time 10-15 seconds

### Examples

| Scenario | Impact | Response |
|----------|--------|----------|
| **User reports "can't upload VCF file"** | Single user workflow blocked | Check user permissions, file size limits, storage quota |
| **mcp-tcga (mocked server) returning errors** | Research feature unavailable (low impact) | Check logs, restart container if needed, low priority fix |
| **Grafana dashboard showing incorrect metrics** | Monitoring data inaccurate | Verify data source, refresh dashboard, update queries |
| **User requests access to new FHIR resource type** | Feature request / permission change | Review IAM roles, grant access if approved |

### Response Timeline

| Time | Action | Owner |
|------|--------|-------|
| **0-4 hours** | Ticket created, acknowledge issue | Support Team / On-Call |
| **4 hours** | Team lead notified if not acknowledged | Automated |
| **≤ 3 business days** | Issue resolved or workaround provided | Support Team |
| **Weekly** | Status update in ticket system | Assigned Engineer |

### Escalation Path

```
Support Team / On-Call Engineer (0 hours)
    ↓ (4 hours if not acknowledged)
Team Lead (4 hours)
```

### Communication Channels

- **Email:** it-oncall@hospital.org
- **Slack:** #precision-medicine-alerts
- **Ticket System:** Jira / ServiceNow

---

## P3 - Low Priority Incidents

### Definition
**Cosmetic issues, minor bugs, or enhancement requests.**

### Criteria

✅ **Cosmetic Issues:**
- UI alignment issues
- Typos in documentation
- Non-functional buttons (no impact on workflows)

✅ **Enhancement Requests:**
- Feature requests (new capabilities)
- Performance optimization (non-urgent)
- Configuration changes (nice-to-have)

✅ **Non-Critical Bugs:**
- Error messages unclear
- Logging too verbose / not verbose enough
- Minor inconsistencies in reporting

### Examples

| Scenario | Impact | Response |
|----------|--------|----------|
| **UI button misaligned on analysis results page** | Cosmetic issue, no functional impact | Add to backlog, fix in next sprint |
| **User requests "can you add a new chart to dashboard?"** | Enhancement request | Review feasibility, prioritize in roadmap |
| **Error message says "Error 500" instead of descriptive message** | Poor UX, but system functional | Improve error handling, update in next release |

### Response Timeline

| Time | Action | Owner |
|------|--------|-------|
| **0-24 hours** | Ticket created, acknowledge receipt | Support Team |
| **≤ 7 business days** | Issue triaged, added to backlog | Product Team |
| **As scheduled** | Resolved in upcoming sprint/release | Engineering Team |

### Escalation Path

No automatic escalation. Handled through standard backlog prioritization.

### Communication Channels

- **Ticket System:** Jira / ServiceNow only
- No PagerDuty, email, or Slack alerts

---

## Incident Severity Decision Tree

Use this flowchart to determine incident severity:

```
┌─────────────────────────────────────┐
│  Is the entire system down?         │
│  OR PHI breach detected?            │
│  OR Audit log failure?              │
└────────────┬────────────────────────┘
             │
         YES │
             ↓
        ┌────────┐
        │   P0   │ ← 15 min response, 4 hour resolution
        └────────┘
             │
         NO  │
             ↓
┌─────────────────────────────────────┐
│  Is a Tier 1 server down?           │
│  OR >5% error rate?                 │
│  OR Multiple users affected?        │
└────────────┬────────────────────────┘
             │
         YES │
             ↓
        ┌────────┐
        │   P1   │ ← 1 hour response, 8 hour resolution
        └────────┘
             │
         NO  │
             ↓
┌─────────────────────────────────────┐
│  Is a single user affected?         │
│  OR Tier 2/3 server issue?          │
│  OR Minor performance degradation?  │
└────────────┬────────────────────────┘
             │
         YES │
             ↓
        ┌────────┐
        │   P2   │ ← 4 hour response, 3 days resolution
        └────────┘
             │
         NO  │
             ↓
┌─────────────────────────────────────┐
│  Is it a cosmetic issue?            │
│  OR enhancement request?            │
│  OR minor bug?                      │
└────────────┬────────────────────────┘
             │
         YES │
             ↓
        ┌────────┐
        │   P3   │ ← 24 hour response, 7 days resolution
        └────────┘
```

---

## Service Credit Eligibility by Priority

Only P0 and P1 incidents that violate SLA targets are eligible for service credits.

| Priority | SLA Violation | Service Credit |
|----------|---------------|----------------|
| **P0** | Response time >15 min OR Resolution time >4 hours | 25-100% of monthly bill |
| **P1** | Response time >1 hour OR Resolution time >8 hours | 10-50% of monthly bill |
| **P2** | Response time >4 hours OR Resolution time >3 days | No service credit (informational) |
| **P3** | Response time >24 hours OR Resolution time >7 days | No service credit (informational) |

**Service Credit Calculation Example:**

```
Scenario: P0 incident took 6 hours to resolve (exceeded 4-hour target by 2 hours)

Monthly infrastructure bill: $1,000
SLA violation: 50% (exceeded target by 50%)
Service credit: $500 (50% of $1,000)

Customer pays: $500 for that month
```

---

## After-Hours and Holiday Coverage

### On-Call Rotation

| Time Period | Coverage | Response Time | Escalation |
|-------------|----------|---------------|------------|
| **Business Hours** (Mon-Fri, 8 AM - 5 PM) | Full team available | Standard (per priority) | Standard |
| **After Hours** (Mon-Fri, 5 PM - 8 AM) | On-call engineer | Standard (per priority) | Escalate to IT Manager if needed |
| **Weekends** (Sat-Sun) | On-call engineer | Standard (per priority) | Escalate to IT Manager if needed |
| **Holidays** (US federal holidays) | On-call engineer | Standard (per priority) | Escalate to IT Director if P0 |

**Note:** Response and resolution times are the same 24/7/365. No exceptions for after-hours or holidays.

### Holiday Escalation

For P0 incidents on major holidays (Thanksgiving, Christmas, New Year):
- On-call engineer responds within 15 minutes
- IT Manager notified within 30 minutes
- IT Director notified within 1 hour
- Executive team notified within 2 hours (if not resolved)

---

## Related Documents

- [Section 3: Incident Management](../sections/03_INCIDENT_MANAGEMENT.md) - Detailed response procedures
- [Runbook: P0 Complete System Outage](../runbooks/P0_COMPLETE_OUTAGE.md) - Step-by-step response
- [Runbook: P1 PHI Breach Response](../runbooks/P1_PHI_BREACH.md) - HIPAA breach procedures
- [Appendix E: Contact and Escalation Directory](APPENDIX_E_CONTACTS.md) - On-call contacts

[← Back to SLA Overview](../README.md)
