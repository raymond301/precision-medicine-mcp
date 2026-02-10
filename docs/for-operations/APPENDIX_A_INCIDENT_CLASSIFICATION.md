# Appendix A: Incident Classification Matrix

[← Back to SLA Overview](README.md)

---

## Incident Severity Matrix

| Severity | Description | Examples | Response Target |
|----------|-------------|----------|-----------------|
| **P0 - Critical** | Complete system outage or PHI breach | api.anthropic.com unreachable, Database breach | 15 minutes |
| **P1 - High** | Significant performance degradation | Latency >10s, FHIR store errors | 1 hour |
| **P2 - Medium** | Partial outage affecting single tools | mcp-fgbio failing, Single user unable to login | 4 hours |
| **P3 - Low** | Minor UI bugs or cosmetic issues | Typo in dashboard, Icon alignment | 24 hours |

---

## Related Documents

- [Section 03: Incident Management](03_INCIDENT_MANAGEMENT.md)
- [Runbook: P0 Complete System Outage](RUNBOOK_P0_COMPLETE_OUTAGE.md)

---

### Document History

| Date | Version | Author | Change Summary |
|:---|:---|:---|:---|
| 2026-02-08 | 2.0 | IT Operations | Initial creation |
