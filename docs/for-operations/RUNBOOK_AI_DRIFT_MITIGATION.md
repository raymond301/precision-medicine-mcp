# Runbook: AI Drift and Bias Mitigation

[← Back to SLA Overview](README.md)

---

## Overview
This runbook provides recovery steps when AI performance metrics (defined in `APPENDIX_C`) exceed threshold limits (Drift > 10% or Accuracy < 85%).

## 1. Detection and Triage
1. Review the **AI Health Dashboard** (Langfuse/Datadog).
2. Compare current P95 accuracy against the last 7-day baseline.
3. Identify if the drift is **Data-driven** (new variant types) or **Model-driven** (API update/version change).

## 2. Mitigation Steps

### Step 1: Version Rollback
If a recent deployment caused the drift:
```bash
# Roll back to the previous stable revision
gcloud run services update <service-name> --region=us-central1 --rollback-to-stable
```

### Step 2: Prompt Optimization
If the model behavior has changed due to vendor API updates:
1. Update the System Prompt in the `mcp-server` configuration.
2. Run the `test-6-citl-review.md` scenario to verify output alignment.

### Step 3: Human-in-the-Loop Override
Trigger a mandatory manual review for all clinical reports until the performance baseline is restored.

---

## Verification checklist
- [ ] P95 Accuracy restored to > 85%.
- [ ] Bias audit logs show zero high-confidence hallucinations.

---

### Document History

| Date | Version | Author | Change Summary |
|:---|:---|:---|:---|
| 2026-02-08 | 2.0 | IT Operations | Initial creation |
