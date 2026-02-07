"""
Live MCP Server Monitor

Two live data channels for the dashboard:
  1. Health polling  - async httpx GET /health against each Cloud Run service
                       10 s timeout, 2 retries (handles cold-start)
  2. Cloud Logging   - GCP log queries for request counts, latency, error rates
                       over a user-selected time window

Both are called from streamlit_app.py when the user enables Live Mode.
"""

import asyncio
import os
import time
from datetime import datetime, timedelta, timezone
from typing import Any, Dict, List

import httpx

# ── Deployed Cloud Run base URLs ────────────────────────────────────────────
# Must stay in sync with infrastructure/deployment/deploy_to_gcp.sh
# and ui/streamlit-app/utils/mcp_config.py
MCP_SERVERS: Dict[str, str] = {
    "mcp-fgbio":                     "https://mcp-fgbio-ondu7mwjpa-uc.a.run.app",
    "mcp-multiomics":                "https://mcp-multiomics-ondu7mwjpa-uc.a.run.app",
    "mcp-spatialtools":              "https://mcp-spatialtools-ondu7mwjpa-uc.a.run.app",
    "mcp-perturbation":              "https://mcp-perturbation-ondu7mwjpa-uc.a.run.app",
    "mcp-quantum-celltype-fidelity": "https://mcp-quantum-celltype-fidelity-ondu7mwjpa-uc.a.run.app",
    "mcp-deepcell":                  "https://mcp-deepcell-ondu7mwjpa-uc.a.run.app",
    "mcp-tcga":                      "https://mcp-tcga-ondu7mwjpa-uc.a.run.app",
    "mcp-openimagedata":             "https://mcp-openimagedata-ondu7mwjpa-uc.a.run.app",
    "mcp-mockepic":                  "https://mcp-mockepic-ondu7mwjpa-uc.a.run.app",
    "mcp-seqera":                    "https://mcp-seqera-ondu7mwjpa-uc.a.run.app",
    "mcp-huggingface":               "https://mcp-huggingface-ondu7mwjpa-uc.a.run.app",
    "mcp-patient-report":            "https://mcp-patient-report-ondu7mwjpa-uc.a.run.app",
}

TIMEOUT     = 10   # seconds per attempt
MAX_RETRIES = 2    # retries after first failure
HEALTH_PATH = "/health"

GCP_PROJECT_ID = os.getenv("GCP_PROJECT_ID", "precision-medicine-poc")
GCP_REGION     = os.getenv("GCP_REGION",      "us-central1")


# ════════════════════════════════════════════════════════════════════════════
# 1. Health polling (httpx async)
# ════════════════════════════════════════════════════════════════════════════


async def _check_one(
    client: httpx.AsyncClient,
    name: str,
    base_url: str,
) -> Dict[str, Any]:
    """Health-check one server with retry.  Returns a status dict."""
    url = f"{base_url}{HEALTH_PATH}"
    for attempt in range(MAX_RETRIES + 1):
        t0 = time.perf_counter()
        try:
            resp = await client.get(url, timeout=TIMEOUT)
            ms   = round((time.perf_counter() - t0) * 1000, 1)
            return {
                "server":     name,
                "status":     "healthy" if resp.status_code == 200 else "degraded",
                "latency_ms": ms,
                "checked_at": datetime.now(timezone.utc).isoformat(),
                "error":      None if resp.status_code == 200 else f"HTTP {resp.status_code}",
            }
        except (httpx.TimeoutException, httpx.ConnectError) as exc:
            if attempt < MAX_RETRIES:
                await asyncio.sleep(1)
                continue
            return {
                "server":     name,
                "status":     "unhealthy",
                "latency_ms": round((time.perf_counter() - t0) * 1000, 1),
                "checked_at": datetime.now(timezone.utc).isoformat(),
                "error":      type(exc).__name__,
            }

    # safety fallback — should never reach here
    return {
        "server": name, "status": "unknown", "latency_ms": 0,
        "checked_at": datetime.now(timezone.utc).isoformat(), "error": "unknown",
    }


async def _poll_all() -> Dict[str, Dict[str, Any]]:
    """Fan-out health checks to all servers concurrently."""
    async with httpx.AsyncClient() as client:
        results = await asyncio.gather(*(
            _check_one(client, name, url)
            for name, url in MCP_SERVERS.items()
        ))
    return {r["server"]: r for r in results}


def get_live_health() -> Dict[str, Dict[str, Any]]:
    """Sync entry-point for Streamlit.  Spawns the async poll in a thread."""
    import concurrent.futures
    with concurrent.futures.ThreadPoolExecutor(max_workers=1) as pool:
        return pool.submit(asyncio.run, _poll_all()).result(timeout=45)


# ════════════════════════════════════════════════════════════════════════════
# 2. Cloud Logging queries (request counts, latency, errors)
# ════════════════════════════════════════════════════════════════════════════


def _parse_latency(s: str) -> float:
    """Parse Cloud Run latency string '0.123456s' -> milliseconds."""
    try:
        return float(s.rstrip("s")) * 1000
    except (ValueError, AttributeError):
        return 0.0


def _logging_client():
    from google.cloud import logging_v2
    return logging_v2.Client(project=GCP_PROJECT_ID)


def query_server_logs(
    name: str,
    hours_back: int = 1,
    client=None,
) -> Dict[str, Any]:
    """
    Query Cloud Logging for one Cloud Run service and return aggregated stats.

    Cloud Run automatically emits structured HTTP-request logs; we filter by
    service name + time window and compute request count, avg/p95 latency,
    and error rate.
    """
    start = (datetime.now(timezone.utc) - timedelta(hours=hours_back)).isoformat()
    log_filter = (
        f'resource.type="cloud_run_revision" '
        f'AND resource.labels.service_name="{name}" '
        f'AND resource.labels.location="{GCP_REGION}" '
        f'AND httpRequest.requestMethod != "" '
        f'AND timestamp >= "{start}"'
    )

    latencies: List[float] = []
    request_count = error_count = 0

    try:
        if client is None:
            client = _logging_client()
        for entry in client.list_entries(filter_=log_filter, page_size=500):
            request_count += 1
            http_req = entry.to_dict().get("httpRequest", {})
            latencies.append(_parse_latency(http_req.get("latency", "0s")))
            if int(http_req.get("status", 200)) >= 400:
                error_count += 1
    except Exception:
        pass  # graceful degradation — caller sees zeros

    avg_lat = sum(latencies) / len(latencies) if latencies else 0.0
    p95_lat = (
        sorted(latencies)[min(int(0.95 * len(latencies)), len(latencies) - 1)]
        if latencies else 0.0
    )

    return {
        "server":            name,
        "request_count":     request_count,
        "avg_latency_ms":    round(avg_lat, 1),
        "p95_latency_ms":    round(p95_lat, 1),
        "error_count":       error_count,
        "error_rate":        round(error_count / request_count, 3) if request_count else 0.0,
        "time_window_hours": hours_back,
    }


def query_all_server_logs(hours_back: int = 1) -> Dict[str, Dict[str, Any]]:
    """Query Cloud Logging for every deployed MCP server.  Single client reused."""
    try:
        client = _logging_client()
    except Exception:
        client = None  # each per-server call will gracefully return zeros
    return {name: query_server_logs(name, hours_back, client) for name in MCP_SERVERS}
