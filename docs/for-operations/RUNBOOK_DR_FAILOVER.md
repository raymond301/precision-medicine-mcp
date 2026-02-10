# Runbook: Disaster Recovery Failover

[← Back to SLA Overview](README.md)

---

## Procedure

1. **Declare disaster:** IT Director approval required.
2. **Update DNS:** Point Cloud DNS to us-east1 load balancer IP.
3. **Restore Data:** Import latest FHIR snapshots into us-east1 FHIR store.
4. **Verify:** Run synthetic health checks against the new domain.

---

## Related Documents

- [Section 07: Backup and Disaster Recovery](07_BACKUP_DISASTER_RECOVERY.md)
- [Runbook: P0 Complete System Outage](RUNBOOK_P0_COMPLETE_OUTAGE.md)

---

### Document History

| Date | Version | Author | Change Summary |
|:---|:---|:---|:---|
| 2026-02-08 | 2.0 | IT Operations | Initial creation |
