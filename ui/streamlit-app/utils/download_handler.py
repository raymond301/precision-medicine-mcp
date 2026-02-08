"""Download handler for file content from MCP tool results.

Extracts downloadable files (PDFs, etc.) from tool results and provides
Streamlit download buttons.
"""

import base64
import json
import re
import streamlit as st
from typing import Optional, List, Dict, Any


def extract_downloadable_files(trace: Any) -> List[Dict[str, Any]]:
    """Extract downloadable files from an orchestration trace.

    Looks for tool results containing file_content_base64 fields.

    Args:
        trace: OrchestrationTrace object with tool_calls

    Returns:
        List of dicts with file_name, file_content (bytes), mime_type
    """
    files = []

    if not trace or not hasattr(trace, 'tool_calls'):
        return files

    for tool_call in trace.tool_calls:
        if not tool_call.result:
            continue

        # Parse result if it's a string
        result = tool_call.result
        if isinstance(result, str):
            try:
                result = json.loads(result)
            except (json.JSONDecodeError, TypeError):
                continue

        # Check for file_content_base64
        if isinstance(result, dict) and result.get("file_content_base64"):
            try:
                file_content = base64.b64decode(result["file_content_base64"])
                file_name = result.get("file_name", "report.pdf")
                output_format = result.get("output_format", "pdf")

                # Determine MIME type
                if output_format == "pdf":
                    mime_type = "application/pdf"
                elif output_format == "html":
                    mime_type = "text/html"
                else:
                    mime_type = "application/octet-stream"

                files.append({
                    "file_name": file_name,
                    "file_content": file_content,
                    "mime_type": mime_type,
                    "patient_id": result.get("patient_id", "unknown"),
                    "report_type": result.get("report_type", "report"),
                    "is_draft": result.get("is_draft", True)
                })
            except Exception as e:
                # Log error but continue
                st.warning(f"Failed to decode file content: {e}")
                continue

    return files


def extract_downloadable_from_content(content: str) -> List[Dict[str, Any]]:
    """Extract downloadable files from message content.

    Searches for JSON blocks with file_content_base64 in the text.

    Args:
        content: Message content string

    Returns:
        List of dicts with file info
    """
    files = []

    # Look for JSON blocks that might contain file_content_base64
    # Pattern: {"status": ..., "file_content_base64": ...}
    json_pattern = r'\{[^{}]*"file_content_base64"\s*:\s*"[^"]+?"[^{}]*\}'

    # Also try to find in markdown code blocks
    code_block_pattern = r'```(?:json)?\s*(\{.*?"file_content_base64".*?\})\s*```'

    for pattern in [json_pattern, code_block_pattern]:
        matches = re.findall(pattern, content, re.DOTALL)
        for match in matches:
            try:
                data = json.loads(match)
                if data.get("file_content_base64"):
                    file_content = base64.b64decode(data["file_content_base64"])
                    file_name = data.get("file_name", "report.pdf")
                    output_format = data.get("output_format", "pdf")

                    mime_type = "application/pdf" if output_format == "pdf" else "text/html"

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
