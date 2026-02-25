# Deployment Guide

For comprehensive deployment information, see:

## [GCP Integration Guide](../../testing/gcp-integration.md)
How to test deployed servers using Claude API (HTTP/SSE remote access).

## [Security Guide](../../../for-developers/security.md)
API key and secrets management, GCP Secret Manager integration.

---

## Quick Reference

**Deployment Script:** `/infrastructure/deployment/deploy_to_gcp.sh`

```bash
# Deploy single server (development)
./infrastructure/deployment/deploy_to_gcp.sh --development --server mcp-spatialtools

# Deploy all servers
./infrastructure/deployment/deploy_to_gcp.sh --development
```

**Current Status:**
- **Deployed:** All servers on GCP Cloud Run except mcp-epic (local-only) — see [Server Registry](../../shared/server-registry.md)
- **Mode:** Development (public) + Production (epic - private)
- **Region:** us-central1

---

For deployment procedures, testing, and troubleshooting, see the [GCP Integration Guide](../../testing/gcp-integration.md).
