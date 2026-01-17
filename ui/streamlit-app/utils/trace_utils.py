"""
Utilities for extracting and formatting MCP orchestration traces.

This module parses Claude API responses to extract tool calls,
providing visibility into the AI's decision-making process.
"""

from typing import List, Dict, Any, Optional
from dataclasses import dataclass
from datetime import datetime
import json


@dataclass
class ToolCall:
    """Represents a single MCP tool call."""
    step_number: int
    server_name: str
    tool_name: str
    input_params: Dict[str, Any]
    result_summary: Optional[str] = None
    duration_ms: Optional[float] = None
    tool_use_id: Optional[str] = None
    timestamp: Optional[datetime] = None


@dataclass
class OrchestrationTrace:
    """Complete trace of an orchestrated query."""
    query: str
    tool_calls: List[ToolCall]
    total_duration_ms: float
    total_tokens: int
    estimated_cost_usd: float

    @property
    def servers_used(self) -> List[str]:
        """Unique servers called in this trace."""
        return list(dict.fromkeys(tc.server_name for tc in self.tool_calls))

    @property
    def call_count(self) -> int:
        """Total number of tool calls."""
        return len(self.tool_calls)


def extract_tool_calls_from_response(response: Any) -> List[ToolCall]:
    """
    Extract tool calls from a Claude API response.

    Works with both the Anthropic Python SDK response objects
    and raw JSON responses.

    Args:
        response: Claude API response (SDK object or dict)

    Returns:
        List of ToolCall objects in order of execution
    """
    tool_calls = []
    step_number = 0

    # Handle SDK response object
    if hasattr(response, 'content'):
        content_blocks = response.content
    # Handle dict response
    elif isinstance(response, dict) and 'content' in response:
        content_blocks = response['content']
    else:
        return tool_calls

    for block in content_blocks:
        # Handle SDK ContentBlock objects
        if hasattr(block, 'type'):
            block_type = block.type
            if block_type == 'tool_use':
                step_number += 1
                server_name, tool_name = _parse_tool_name(
                    getattr(block, 'name', 'unknown')
                )
                tool_calls.append(ToolCall(
                    step_number=step_number,
                    server_name=server_name,
                    tool_name=tool_name,
                    input_params=getattr(block, 'input', {}),
                    tool_use_id=getattr(block, 'id', None)
                ))
        # Handle dict blocks
        elif isinstance(block, dict):
            if block.get('type') == 'tool_use':
                step_number += 1
                server_name, tool_name = _parse_tool_name(
                    block.get('name', 'unknown')
                )
                tool_calls.append(ToolCall(
                    step_number=step_number,
                    server_name=server_name,
                    tool_name=tool_name,
                    input_params=block.get('input', {}),
                    tool_use_id=block.get('id')
                ))

    return tool_calls


def extract_tool_results(messages: List[Dict]) -> Dict[str, str]:
    """
    Extract tool results from conversation messages.

    Args:
        messages: List of conversation messages

    Returns:
        Dict mapping tool_use_id to result summary
    """
    results = {}

    for msg in messages:
        if msg.get('role') == 'user':
            content = msg.get('content', [])
            if isinstance(content, list):
                for block in content:
                    if isinstance(block, dict) and block.get('type') == 'tool_result':
                        tool_use_id = block.get('tool_use_id')
                        result_content = block.get('content', '')
                        if tool_use_id:
                            # Truncate long results for display
                            results[tool_use_id] = _truncate_result(result_content)

    return results


def _parse_tool_name(full_name: str) -> tuple:
    """
    Parse a full tool name into server and tool components.

    Examples:
        'mcp-fgbio.analyze_variants' -> ('mcp-fgbio', 'analyze_variants')
        'get_patient' -> ('unknown', 'get_patient')
    """
    if '.' in full_name:
        parts = full_name.split('.', 1)
        return parts[0], parts[1]
    elif '-' in full_name and full_name.startswith('mcp'):
        # Handle 'mcp-fgbio_analyze_variants' format
        parts = full_name.split('_', 1)
        if len(parts) == 2:
            return parts[0], parts[1]
    return 'unknown-server', full_name


def _truncate_result(result: str, max_length: int = 200) -> str:
    """Truncate a result string for display."""
    if isinstance(result, dict):
        result = json.dumps(result)
    result = str(result)
    if len(result) > max_length:
        return result[:max_length] + "..."
    return result


def build_orchestration_trace(
    query: str,
    response: Any,
    messages: List[Dict],
    duration_ms: float = 0,
    tokens: int = 0,
    cost_usd: float = 0
) -> OrchestrationTrace:
    """
    Build a complete orchestration trace from a query and response.

    Args:
        query: The user's original query
        response: Claude API response
        messages: Full conversation messages (for extracting results)
        duration_ms: Total request duration in milliseconds
        tokens: Total tokens used
        cost_usd: Estimated cost in USD

    Returns:
        OrchestrationTrace object
    """
    tool_calls = extract_tool_calls_from_response(response)
    results = extract_tool_results(messages)

    # Match results to tool calls
    for tc in tool_calls:
        if tc.tool_use_id and tc.tool_use_id in results:
            tc.result_summary = results[tc.tool_use_id]

    return OrchestrationTrace(
        query=query,
        tool_calls=tool_calls,
        total_duration_ms=duration_ms,
        total_tokens=tokens,
        estimated_cost_usd=cost_usd
    )


# Server metadata for display
SERVER_INFO = {
    'mcp-mockepic': {
        'display_name': 'Clinical (Epic FHIR)',
        'icon': '🏥',
        'color': '#e2e3e5',
        'description': 'Patient demographics, conditions, medications'
    },
    'mcp-fgbio': {
        'display_name': 'Genomics (FGbio)',
        'icon': '🧬',
        'color': '#d4edda',
        'description': 'FASTQ/VCF QC, variant analysis'
    },
    'mcp-tcga': {
        'display_name': 'Cancer Genomics (TCGA)',
        'icon': '📊',
        'color': '#f8d7da',
        'description': 'TCGA cohort comparisons'
    },
    'mcp-multiomics': {
        'display_name': 'Multi-Omics',
        'icon': '🔬',
        'color': '#d4edda',
        'description': 'RNA, protein, phospho integration'
    },
    'mcp-spatialtools': {
        'display_name': 'Spatial Transcriptomics',
        'icon': '🗺️',
        'color': '#d4edda',
        'description': 'Visium analysis, spatial statistics'
    },
    'mcp-openimagedata': {
        'display_name': 'Imaging (Histology)',
        'icon': '🔍',
        'color': '#fff3cd',
        'description': 'H&E, multiplex IF analysis'
    },
    'mcp-deepcell': {
        'display_name': 'Cell Segmentation',
        'icon': '🔲',
        'color': '#f8d7da',
        'description': 'Cell detection and quantification'
    },
    'mcp-huggingface': {
        'display_name': 'ML Models',
        'icon': '🤗',
        'color': '#f8d7da',
        'description': 'Biomedical ML models'
    },
    'mcp-seqera': {
        'display_name': 'Workflows (Nextflow)',
        'icon': '⚙️',
        'color': '#f8d7da',
        'description': 'Pipeline orchestration'
    }
}


def get_server_info(server_name: str) -> Dict[str, str]:
    """Get display info for a server."""
    return SERVER_INFO.get(server_name, {
        'display_name': server_name,
        'icon': '🔧',
        'color': '#e9ecef',
        'description': 'MCP Server'
    })
