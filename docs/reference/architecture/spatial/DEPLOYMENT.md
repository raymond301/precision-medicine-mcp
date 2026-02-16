# Deployment Guide

For comprehensive deployment information, see:

## 🧪 [GCP Testing Guide](../../deployment/GCP_TESTING_GUIDE.md)
How to test deployed servers using Claude API (HTTP/SSE remote access).

## 🔒 [Security Guide](../../deployment/security.md)
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
- **Deployed:** 14/15 servers on GCP Cloud Run (all except mcp-epic which is local-only)
- **Mode:** Development (public) + Production (epic - private)
- **Region:** us-central1

---

For deployment procedures, testing, and troubleshooting, see the [GCP Testing Guide](../../deployment/GCP_TESTING_GUIDE.md).
