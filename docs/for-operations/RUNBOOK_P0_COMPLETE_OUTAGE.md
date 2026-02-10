# Runbook: P0 Complete System Outage

[← Back to SLA Overview](README.md)

---

## Procedure

1. **Check Anthropic API status:**
   `curl -I https://api.anthropic.com/v1/messages`
   If 5xx or timeout → Wait for resolution.

2. **Check GCP region health:**
   `gcloud compute regions describe us-central1`
   If degraded → Initiate DR failover to us-east1.

3. **Check Cloud Run service status:**
   `gcloud run services list --region=us-central1`
   If NOT "Ready" → Check logs for crash loop.

---

## Related Documents

- [Section 03: Incident Management](03_INCIDENT_MANAGEMENT.md)
- [Appendix A: Incident Classification](APPENDIX_A_INCIDENT_CLASSIFICATION.md)

---

### Document History

| Date | Version | Author | Change Summary |
|:---|:---|:---|:---|
| 2026-02-08 | 2.0 | IT Operations | Initial creation |
