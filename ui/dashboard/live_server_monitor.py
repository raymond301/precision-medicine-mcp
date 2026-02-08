"""
Live MCP Server Monitor

Three live data channels for the dashboard:
  1. Health polling  - async httpx GET / against each Cloud Run service
                       10 s timeout, 2 retries (handles cold-start)
  2. Cloud Logging   - GCP log queries for request counts, latency, error rates
                       over a user-selected time window
  3. Token Usage     - Queries mcp-audit-log for LLM token usage and costs
                       logged by Streamlit clients via audit_logger

Monitors both MCP servers and Streamlit client apps.
Called from streamlit_app.py when the user enables Live Mode.
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

# MCP Servers (FastMCP with /sse endpoint)
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

# Streamlit Client Apps (standard HTTP with root path /)
STREAMLIT_CLIENTS: Dict[str, str] = {
    "mcp-dashboard":                 "https://mcp-dashboard-305650208648.us-central1.run.app",
    "streamlit-mcp-chat":            "https://streamlit-mcp-chat-305650208648.us-central1.run.app",
    "streamlit-mcp-chat-students":   "https://streamlit-mcp-chat-students-305650208648.us-central1.run.app",
}

# Combined for monitoring
ALL_SERVICES: Dict[str, str] = {**MCP_SERVERS, **STREAMLIT_CLIENTS}

TIMEOUT     = 10   # seconds per attempt
MAX_RETRIES = 2    # retries after first failure
# MCP servers use /sse, Streamlit apps use / - we check root for all
HEALTH_PATH = "/"

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
    """Health-check one service with retry.  Returns a status dict.

    MCP servers don't have a /health endpoint - they expose /sse for SSE transport.
    Streamlit apps respond at /. We check root path and treat any HTTP response
    (even 404/405) as "healthy" since it proves the service is running.
    Only connection failures or timeouts indicate an unhealthy service.
    """
    url = f"{base_url}{HEALTH_PATH}"
    for attempt in range(MAX_RETRIES + 1):
        t0 = time.perf_counter()
        try:
            resp = await client.get(url, timeout=TIMEOUT)
            ms   = round((time.perf_counter() - t0) * 1000, 1)
            # Any HTTP response means the service is running
            # 200 = healthy (Streamlit), 404/405 = healthy (MCP server responding)
            # Only 5xx errors indicate degraded state
            if resp.status_code < 500:
                status = "healthy"
                error = None
            else:
                status = "degraded"
                error = f"HTTP {resp.status_code}"
            return {
                "server":     name,
                "status":     status,
                "latency_ms": ms,
                "checked_at": datetime.now(timezone.utc).isoformat(),
                "error":      error,
                "http_code":  resp.status_code,
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
    """Fan-out health checks to all services (MCP servers + Streamlit clients) concurrently."""
    async with httpx.AsyncClient() as client:
        results = await asyncio.gather(*(
            _check_one(client, name, url)
            for name, url in ALL_SERVICES.items()
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
    """Query Cloud Logging for every deployed service (MCP servers + Streamlit clients).
    Single client reused."""
    try:
        client = _logging_client()
    except Exception:
        client = None  # each per-server call will gracefully return zeros
    return {name: query_server_logs(name, hours_back, client) for name in ALL_SERVICES}


# ════════════════════════════════════════════════════════════════════════════
# 3. Token Usage & Cost queries (from audit log)
# ════════════════════════════════════════════════════════════════════════════


def query_token_usage(hours_back: int = 24, client=None) -> Dict[str, Any]:
    """
    Query the mcp-audit-log for token usage and costs from Streamlit clients.

    The Streamlit apps log every LLM response with token counts and estimated costs
    via audit_logger.log_mcp_response(). This function aggregates that data.

    Args:
        hours_back: How many hours of history to query (default: 24)
        client: Optional logging client (reused for efficiency)

    Returns:
        Dictionary with aggregated token usage and costs:
        - total_input_tokens: Sum of all input tokens
        - total_output_tokens: Sum of all output tokens
        - total_tokens: Sum of all tokens
        - total_cost_usd: Sum of estimated costs
        - query_count: Number of LLM queries
        - avg_tokens_per_query: Average tokens per query
        - avg_cost_per_query: Average cost per query
        - by_model: Breakdown by model name
        - by_server: Breakdown by MCP server used
        - time_window_hours: Query time window
    """
    start = (datetime.now(timezone.utc) - timedelta(hours=hours_back)).isoformat()

    # Query audit log for mcp_response events
    log_filter = (
        f'logName="projects/{GCP_PROJECT_ID}/logs/mcp-audit-log" '
        f'AND jsonPayload.event="mcp_response" '
        f'AND timestamp >= "{start}"'
    )

    total_input = total_output = total_tokens = 0
    total_cost = 0.0
    query_count = 0
    by_model: Dict[str, Dict[str, Any]] = {}
    by_server: Dict[str, Dict[str, Any]] = {}

    try:
        if client is None:
            client = _logging_client()

        for entry in client.list_entries(filter_=log_filter, page_size=1000):
            payload = entry.to_dict().get("jsonPayload", {})

            input_tokens = payload.get("input_tokens", 0)
            output_tokens = payload.get("output_tokens", 0)
            tokens = payload.get("total_tokens", input_tokens + output_tokens)
            cost = payload.get("estimated_cost_usd", 0.0)
            servers = payload.get("servers", [])

            total_input += input_tokens
            total_output += output_tokens
            total_tokens += tokens
            total_cost += cost
            query_count += 1

            # Aggregate by MCP servers used
            for server in servers:
                if server not in by_server:
                    by_server[server] = {
                        "query_count": 0,
                        "total_tokens": 0,
                        "total_cost_usd": 0.0
                    }
                by_server[server]["query_count"] += 1
                by_server[server]["total_tokens"] += tokens
                by_server[server]["total_cost_usd"] += cost

    except Exception as e:
        # Graceful degradation - return zeros with error info
        return {
            "total_input_tokens": 0,
            "total_output_tokens": 0,
            "total_tokens": 0,
            "total_cost_usd": 0.0,
            "query_count": 0,
            "avg_tokens_per_query": 0,
            "avg_cost_per_query": 0.0,
            "by_model": {},
            "by_server": {},
            "time_window_hours": hours_back,
            "error": str(e),
        }

    return {
        "total_input_tokens": total_input,
        "total_output_tokens": total_output,
        "total_tokens": total_tokens,
        "total_cost_usd": round(total_cost, 4),
        "query_count": query_count,
        "avg_tokens_per_query": round(total_tokens / query_count, 1) if query_count else 0,
        "avg_cost_per_query": round(total_cost / query_count, 4) if query_count else 0.0,
        "by_model": by_model,
        "by_server": by_server,
        "time_window_hours": hours_back,
    }


def query_usage_by_user(hours_back: int = 24, client=None) -> List[Dict[str, Any]]:
    """
    Query token usage grouped by user (email hash).

    Args:
        hours_back: How many hours of history to query
        client: Optional logging client

    Returns:
        List of user usage summaries, sorted by total cost descending
    """
    start = (datetime.now(timezone.utc) - timedelta(hours=hours_back)).isoformat()

    log_filter = (
        f'logName="projects/{GCP_PROJECT_ID}/logs/mcp-audit-log" '
        f'AND jsonPayload.event="mcp_response" '
        f'AND timestamp >= "{start}"'
    )

    by_user: Dict[str, Dict[str, Any]] = {}

    try:
        if client is None:
            client = _logging_client()

        for entry in client.list_entries(filter_=log_filter, page_size=1000):
            payload = entry.to_dict().get("jsonPayload", {})

            user_hash = payload.get("user_email_hash", "unknown")
            tokens = payload.get("total_tokens", 0)
            cost = payload.get("estimated_cost_usd", 0.0)

            if user_hash not in by_user:
                by_user[user_hash] = {
                    "user_email_hash": user_hash,
                    "query_count": 0,
                    "total_tokens": 0,
                    "total_cost_usd": 0.0
                }

            by_user[user_hash]["query_count"] += 1
            by_user[user_hash]["total_tokens"] += tokens
            by_user[user_hash]["total_cost_usd"] += cost

    except Exception:
        return []

    # Sort by cost descending
    users = list(by_user.values())
    users.sort(key=lambda x: x["total_cost_usd"], reverse=True)

    # Round costs
    for u in users:
        u["total_cost_usd"] = round(u["total_cost_usd"], 4)

    return users


def get_combined_metrics(hours_back: int = 1) -> Dict[str, Any]:
    """
    Get combined metrics: health, request logs, and token usage.

    This is a convenience function that aggregates all monitoring data
    into a single call for the dashboard.

    Args:
        hours_back: Time window for log queries

    Returns:
        Dictionary with:
        - health: Live health status for all services
        - request_logs: Request counts, latency, errors per service
        - token_usage: Aggregated token usage and costs
        - timestamp: When this data was collected
    """
    try:
        client = _logging_client()
    except Exception:
        client = None

    return {
        "health": get_live_health(),
        "request_logs": query_all_server_logs(hours_back),
        "token_usage": query_token_usage(hours_back * 24, client),  # Token usage over longer window
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }
