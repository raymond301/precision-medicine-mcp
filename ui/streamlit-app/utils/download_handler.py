"""Download handler for file content from MCP tool results.

Extracts downloadable files (PDFs, etc.) from tool results and provides
Streamlit download buttons. Fetches files via HTTP download URL from the
MCP server rather than parsing base64 from message content.
"""

import json
import logging
import re
import streamlit as st
import httpx
from typing import List, Dict, Any

logger = logging.getLogger(__name__)


def extract_downloadable_files(trace: Any) -> List[Dict[str, Any]]:
    """Extract downloadable files from an orchestration trace.

    Note: ToolCall objects only have result_summary (truncated), not full results.
    This function is kept for interface compatibility but returns empty list.
    Real extraction happens in extract_downloadable_from_content().

    Args:
        trace: OrchestrationTrace object with tool_calls

    Returns:
        List of dicts with file_name, file_content (bytes), mime_type
    """
    return []


def extract_downloadable_from_content(content: str) -> List[Dict[str, Any]]:
    """Extract downloadable files from message content.

    Looks for JSON blocks with download_url and fetches the file content
    from the MCP server's HTTP download endpoint.

    Args:
        content: Message content string

    Returns:
        List of dicts with file info
    """
    files = []

    # Look for JSON blocks that contain download_url
    json_pattern = r'\{[^{}]*"download_url"\s*:\s*"[^"]+?"[^{}]*\}'
    code_block_pattern = r'```(?:json)?\s*(\{.*?"download_url".*?\})\s*```'

    for pattern in [json_pattern, code_block_pattern]:
        matches = re.findall(pattern, content, re.DOTALL)
        for match in matches:
            try:
                data = json.loads(match)
                download_url = data.get("download_url")
                if not download_url:
                    continue

                file_name = data.get("file_name", "report.pdf")
                output_format = data.get("output_format", "pdf")
                mime_type = "application/pdf" if output_format == "pdf" else "text/html"

                # Fetch the file from the download URL
                try:
                    resp = httpx.get(download_url, timeout=30.0)
                    resp.raise_for_status()
                    file_content = resp.content
                except httpx.HTTPError as e:
                    logger.warning(f"Failed to download {download_url}: {e}")
                    continue

                files.append({
                    "file_name": file_name,
                    "file_content": file_content,
                    "mime_type": mime_type,
                    "patient_id": data.get("patient_id", "unknown"),
                    "report_type": data.get("report_type", "report"),
                    "is_draft": data.get("is_draft", True)
                })
            except (json.JSONDecodeError, TypeError, Exception):
                continue

    return files


def render_download_buttons(files: List[Dict[str, Any]], key_prefix: str = ""):
    """Render Streamlit download buttons for files.

    Args:
        files: List of file dicts from extract_downloadable_files
        key_prefix: Prefix for Streamlit widget keys (for uniqueness)
    """
    if not files:
        return

    st.markdown("---")
    st.markdown("📥 **Generated Reports**")

    for i, file_info in enumerate(files):
        col1, col2 = st.columns([3, 1])

        with col1:
            draft_badge = "🏷️ DRAFT" if file_info.get("is_draft") else "✅ Final"
            st.markdown(
                f"**{file_info['file_name']}** {draft_badge}\n\n"
                f"Patient: `{file_info['patient_id']}` | Type: {file_info['report_type']}"
            )

        with col2:
            st.download_button(
                label="⬇️ Download",
                data=file_info["file_content"],
                file_name=file_info["file_name"],
                mime=file_info["mime_type"],
                key=f"{key_prefix}_download_{i}"
            )


def check_and_render_downloads(trace: Any, content: str, key_prefix: str = ""):
    """Check for downloadable files and render download buttons if found.

    Args:
        trace: OrchestrationTrace object (optional)
        content: Message content string
        key_prefix: Prefix for widget keys
    """
    files = []

    # Try extracting from trace first
    if trace:
        files.extend(extract_downloadable_files(trace))

    # Also try extracting from content
    files.extend(extract_downloadable_from_content(content))

    # Deduplicate by file_name
    seen = set()
    unique_files = []
    for f in files:
        if f["file_name"] not in seen:
            seen.add(f["file_name"])
            unique_files.append(f)

    if unique_files:
        render_download_buttons(unique_files, key_prefix)
