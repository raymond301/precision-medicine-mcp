"""Shared MCP Client utilities for Jupyter notebooks.

Contains MCPClient class, server configuration, and helper functions
used across all group-specific example notebooks.
"""

import os
import json
from typing import List, Dict, Any

import getpass

import anthropic

# ---------------------------------------------------------------------------
# Server configuration -- deployed MCP servers on GCP Cloud Run
# ---------------------------------------------------------------------------

MCP_SERVERS: Dict[str, Dict[str, Any]] = {
    # --- Imaging / Cell Analysis ---
    "deepcell": {
        "url": "https://mcp-deepcell-ondu7mwjpa-uc.a.run.app/sse",
        "description": "DeepCell-TF cell segmentation and marker quantification for MxIF",
        "group": "imaging",
        "status": "production",
        "tools_count": 3,
    },
    "cell-classify": {
        "url": "https://mcp-cell-classify-ondu7mwjpa-uc.a.run.app/sse",
        "description": "Cell phenotype classification and visualization (lightweight, no TensorFlow)",
        "group": "imaging",
        "status": "production",
        "tools_count": 3,
    },
    "openimagedata": {
        "url": "https://mcp-openimagedata-ondu7mwjpa-uc.a.run.app/sse",
        "description": "H&E/MxIF image loading, registration, feature extraction, and composite generation",
        "group": "imaging",
        "status": "production",
        "tools_count": 5,
    },
    # --- Genomics / Omics ---
    "fgbio": {
        "url": "https://mcp-fgbio-ondu7mwjpa-uc.a.run.app/sse",
        "description": "Genomic reference data and FASTQ validation",
        "group": "genomics",
        "status": "production",
        "tools_count": 4,
    },
    "multiomics": {
        "url": "https://mcp-multiomics-ondu7mwjpa-uc.a.run.app/sse",
        "description": "Multi-omics integration (RNA/Protein/Phospho)",
        "group": "genomics",
        "status": "production",
        "tools_count": 10,
    },
    "spatialtools": {
        "url": "https://mcp-spatialtools-ondu7mwjpa-uc.a.run.app/sse",
        "description": "Spatial transcriptomics analysis",
        "group": "genomics",
        "status": "production",
        "tools_count": 14,
    },
    "mocktcga": {
        "url": "https://mcp-mocktcga-ondu7mwjpa-uc.a.run.app/sse",
        "description": "Mock TCGA cancer genomics data",
        "group": "genomics",
        "status": "mock",
        "tools_count": 5,
    },
    "perturbation": {
        "url": "https://mcp-perturbation-ondu7mwjpa-uc.a.run.app/sse",
        "description": "GEARS perturbation prediction for treatment response",
        "group": "genomics",
        "status": "production",
        "tools_count": 8,
    },
    "genomic-results": {
        "url": "https://mcp-genomic-results-ondu7mwjpa-uc.a.run.app/sse",
        "description": "Somatic variant/CNV parsing with clinical annotations and HRD scoring",
        "group": "genomics",
        "status": "production",
        "tools_count": 4,
    },
    # --- Clinical ---
    "mockepic": {
        "url": "https://mcp-mockepic-ondu7mwjpa-uc.a.run.app/sse",
        "description": "Mock EHR/FHIR data",
        "group": "clinical",
        "status": "mock",
        "tools_count": 3,
    },
    "patient-report": {
        "url": "https://mcp-patient-report-ondu7mwjpa-uc.a.run.app/sse",
        "description": "Patient-facing PDF reports with plain-language summaries",
        "group": "clinical",
        "status": "production",
        "tools_count": 5,
    },
    # --- Workflow / ML ---
    "quantum-celltype-fidelity": {
        "url": "https://mcp-quantum-celltype-fidelity-ondu7mwjpa-uc.a.run.app/sse",
        "description": "Quantum computing for cell type validation and immune evasion detection",
        "group": "workflow-ml",
        "status": "production",
        "tools_count": 6,
    },
}

# Convenience groupings
SERVER_GROUPS = {
    "imaging":     ["deepcell", "cell-classify", "openimagedata"],
    "genomics":    ["fgbio", "multiomics", "spatialtools", "mocktcga", "perturbation", "genomic-results"],
    "clinical":    ["mockepic", "patient-report"],
    "workflow-ml": ["quantum-celltype-fidelity"],
}


# ---------------------------------------------------------------------------
# MCPClient
# ---------------------------------------------------------------------------

def get_api_key() -> str:
    """Prompt the user for their Anthropic API key.

    The key is requested via an interactive input prompt so it is never
    stored in the notebook, environment, or container image.
    """
    print("Enter your Anthropic API key (get one at https://console.anthropic.com/)")
    key = getpass.getpass("ANTHROPIC_API_KEY: ")
    if not key:
        raise ValueError("No API key provided.")
    return key


class MCPClient:
    """Helper class for calling MCP servers via the Anthropic Claude API."""

    def __init__(self, api_key: str = None):
        self.api_key = api_key or get_api_key()
        self.client = anthropic.Anthropic(api_key=self.api_key)
        self.conversation_history: List[Dict] = []

    def call_servers(
        self,
        prompt: str,
        servers: List[str],
        model: str = "claude-sonnet-4-5",
        max_tokens: int = 4096,
        clear_history: bool = False,
    ) -> Dict[str, Any]:
        """Call one or more MCP servers with a prompt via Claude.

        Args:
            prompt: User query / instruction.
            servers: List of server names to enable (keys in MCP_SERVERS).
            model: Claude model to use.
            max_tokens: Maximum response tokens.
            clear_history: Clear conversation history before this call.

        Returns:
            Dict with 'response', 'usage', 'model', and 'servers_used'.
        """
        if clear_history:
            self.conversation_history = []

        self.conversation_history.append({"role": "user", "content": prompt})

        mcp_servers = [
            {"type": "url", "url": MCP_SERVERS[s]["url"], "name": s}
            for s in servers if s in MCP_SERVERS
        ]
        tools = [
            {"type": "mcp_toolset", "mcp_server_name": s}
            for s in servers if s in MCP_SERVERS
        ]

        response = self.client.beta.messages.create(
            model=model,
            max_tokens=max_tokens,
            messages=self.conversation_history,
            mcp_servers=mcp_servers,
            tools=tools,
            betas=["mcp-client-2025-11-20"],
        )

        response_text = "".join(
            block.text for block in response.content if hasattr(block, "text")
        )

        self.conversation_history.append(
            {"role": "assistant", "content": response_text}
        )

        usage = {
            "input_tokens": response.usage.input_tokens,
            "output_tokens": response.usage.output_tokens,
            "total_tokens": response.usage.input_tokens + response.usage.output_tokens,
            "estimated_cost_usd": self._estimate_cost(
                response.usage.input_tokens, response.usage.output_tokens, model
            ),
        }

        return {
            "response": response_text,
            "usage": usage,
            "model": model,
            "servers_used": servers,
        }

    def clear(self):
        """Clear conversation history."""
        self.conversation_history = []

    @staticmethod
    def _estimate_cost(input_tokens: int, output_tokens: int, model: str) -> float:
        pricing = {
            "claude-sonnet-4-5": {"input": 0.003 / 1000, "output": 0.015 / 1000},
            "claude-opus-4-5":   {"input": 0.015 / 1000, "output": 0.075 / 1000},
            "claude-haiku-4-5":  {"input": 0.001 / 1000, "output": 0.005 / 1000},
        }
        rates = pricing.get(model, pricing["claude-sonnet-4-5"])
        return round(input_tokens * rates["input"] + output_tokens * rates["output"], 4)


# ---------------------------------------------------------------------------
# Display helpers
# ---------------------------------------------------------------------------

def list_servers(group: str = None):
    """Print available servers, optionally filtered by group."""
    servers = MCP_SERVERS
    if group:
        names = SERVER_GROUPS.get(group, [])
        servers = {k: v for k, v in MCP_SERVERS.items() if k in names}

    for name, cfg in servers.items():
        icon = "+" if cfg["status"] == "production" else "~"
        print(f"  [{icon}] {name:30s} {cfg['description']}")


def print_result(result: Dict[str, Any], label: str = "Response"):
    """Pretty-print an MCPClient call result."""
    print(f"--- {label} ---\n")
    print(result["response"])
    print(f"\n{'=' * 70}")
    print(f"Cost: ${result['usage']['estimated_cost_usd']}  |  "
          f"Tokens: {result['usage']['input_tokens']} in, "
          f"{result['usage']['output_tokens']} out")
