"""Mock MCP Client Manager for local development and testing.

Provides realistic tool schemas and responses without requiring Cloud Run MCP servers.
Enable with: USE_MOCK_MCP=true in .env
"""

import sys
from typing import List, Dict, Any, Optional
from contextlib import asynccontextmanager


class MockMCPClientManager:
    """Mock MCP client manager for local testing.

    Simulates MCP server connections and tool calling without requiring
    actual Cloud Run deployments. Useful for:
    - Local development
    - Fast iteration
    - Testing provider logic
    - CI/CD pipelines
    """

    def __init__(self, server_configs: List[Dict]):
        """Initialize mock MCP client manager.

        Args:
            server_configs: List of server configurations (same format as real MCPClientManager)
        """
        self.server_configs = server_configs
        self.sessions: Dict[str, str] = {}
        self.tools_cache: Dict[str, List[Dict]] = {}
        self._gemini_name_map: Dict[str, tuple] = {}

        print(f"DEBUG: Using MOCK MCP manager with {len(server_configs)} servers", file=sys.stderr)

    async def __aenter__(self):
        """Enter async context - simulate MCP connections."""
        print(f"DEBUG: [MOCK] Establishing connections to {len(self.server_configs)} MCP servers", file=sys.stderr)

        # Simulate connecting to each server
        for server in self.server_configs:
            server_name = server["name"]
            self.sessions[server_name] = f"mock_session_{server_name}"
            self.tools_cache[server_name] = self._get_mock_tools(server_name)
            print(f"DEBUG: [MOCK] Connected to {server_name} ({len(self.tools_cache[server_name])} tools)", file=sys.stderr)

        print(f"DEBUG: [MOCK] Connected to {len(self.sessions)} servers", file=sys.stderr)
        return self

    async def __aexit__(self, exc_type, exc_val, exc_tb):
        """Exit async context - simulate closing connections."""
        print(f"DEBUG: [MOCK] Closing all MCP connections", file=sys.stderr)
        self.sessions.clear()
        self.tools_cache.clear()

    def _get_mock_tools(self, server_name: str) -> List[Dict]:
        """Get mock tool definitions for a server.

        Args:
            server_name: Name of the MCP server

        Returns:
            List of mock tool definitions
        """
        # Define mock tools for each known server
        mock_tools = {
            "spatialtools": [
                {
                    "name": "spatial_autocorrelation",
                    "description": "Calculate spatial autocorrelation (Moran's I) for gene expression",
                    "input_schema": {
                        "type": "object",
                        "properties": {
                            "adata_path": {
                                "type": "string",
                                "description": "Path to AnnData h5ad file"
                            },
                            "gene": {
                                "type": "string",
                                "description": "Gene name to analyze"
                            }
                        },
                        "required": ["adata_path", "gene"]
                    }
                },
                {
                    "name": "cell_type_deconvolution",
                    "description": "Deconvolve cell types from spatial transcriptomics data",
                    "input_schema": {
                        "type": "object",
                        "properties": {
                            "adata_path": {
                                "type": "string",
                                "description": "Path to spatial AnnData file"
                            },
                            "reference_path": {
                                "type": "string",
                                "description": "Path to reference single-cell data"
                            }
                        },
                        "required": ["adata_path", "reference_path"]
                    }
                }
            ],
            "multiomics": [
                {
                    "name": "halla_analysis",
                    "description": "Run HAllA (Hierarchical All-against-All) multi-omics association analysis",
                    "input_schema": {
                        "type": "object",
                        "properties": {
                            "dataset1_path": {
                                "type": "string",
                                "description": "Path to first omics dataset (CSV/TSV)"
                            },
                            "dataset2_path": {
                                "type": "string",
                                "description": "Path to second omics dataset (CSV/TSV)"
                            },
                            "output_dir": {
                                "type": "string",
                                "description": "Output directory for results"
                            }
                        },
                        "required": ["dataset1_path", "dataset2_path"]
                    }
                },
                {
                    "name": "stouffer_meta_analysis",
                    "description": "Combine p-values from multiple omics studies using Stouffer's method",
                    "input_schema": {
                        "type": "object",
                        "properties": {
                            "pvalue_files": {
                                "type": "array",
                                "items": {"type": "string"},
                                "description": "List of files containing p-values"
                            }
                        },
                        "required": ["pvalue_files"]
                    }
                }
            ],
            "fgbio": [
                {
                    "name": "validate_fastq",
                    "description": "Validate FASTQ file format and quality metrics",
                    "input_schema": {
                        "type": "object",
                        "properties": {
                            "fastq_path": {
                                "type": "string",
                                "description": "Path to FASTQ file"
                            }
                        },
                        "required": ["fastq_path"]
                    }
                }
            ],
            "quantum-celltype-fidelity": [
                {
                    "name": "calculate_fidelity",
                    "description": "Calculate quantum fidelity between cell type states",
                    "input_schema": {
                        "type": "object",
                        "properties": {
                            "cell_type_a": {
                                "type": "string",
                                "description": "First cell type identifier"
                            },
                            "cell_type_b": {
                                "type": "string",
                                "description": "Second cell type identifier"
                            },
                            "adata_path": {
                                "type": "string",
                                "description": "Path to AnnData file with expression data"
                            }
                        },
                        "required": ["cell_type_a", "cell_type_b", "adata_path"]
                    }
                },
                {
                    "name": "quantum_state_preparation",
                    "description": "Prepare quantum state from gene expression profile",
                    "input_schema": {
                        "type": "object",
                        "properties": {
                            "expression_vector": {
                                "type": "array",
                                "items": {"type": "number"},
                                "description": "Gene expression vector"
                            },
                            "normalization": {
                                "type": "string",
                                "enum": ["l2", "minmax", "zscore"],
                                "description": "Normalization method"
                            }
                        },
                        "required": ["expression_vector"]
                    }
                }
            ]
        }

        # Return mock tools for this server, or generic tools if unknown
        if server_name in mock_tools:
            return mock_tools[server_name]
        else:
            # Generic tools for unknown servers
            return [
                {
                    "name": "generic_tool",
                    "description": f"Generic tool for {server_name} (mock)",
                    "input_schema": {
                        "type": "object",
                        "properties": {
                            "input": {
                                "type": "string",
                                "description": "Generic input parameter"
                            }
                        },
                        "required": ["input"]
                    }
                }
            ]

    async def call_tool(
        self,
        server_name: str,
        tool_name: str,
        arguments: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Mock tool execution.

        Args:
            server_name: Name of the MCP server
            tool_name: Name of the tool to call
            arguments: Tool arguments

        Returns:
            Mock tool result
        """
        print(f"DEBUG: [MOCK] Calling tool {tool_name} on {server_name}", file=sys.stderr)
        print(f"DEBUG: [MOCK] Arguments: {arguments}", file=sys.stderr)

        # Generate realistic mock response based on tool
        mock_response = self._generate_mock_response(server_name, tool_name, arguments)

        print(f"DEBUG: [MOCK] Returning mock result", file=sys.stderr)

        return {
            "content": [
                {
                    "type": "text",
                    "text": mock_response
                }
            ],
            "isError": False
        }

    def _generate_mock_response(self, server_name: str, tool_name: str, arguments: Dict) -> str:
        """Generate realistic mock response for a tool call.

        Args:
            server_name: Server name
            tool_name: Tool name
            arguments: Tool arguments

        Returns:
            Mock response text
        """
        # Tool-specific mock responses
        responses = {
            "spatial_autocorrelation": f"""Spatial Autocorrelation Analysis Results (MOCK DATA)

Gene: {arguments.get('gene', 'Unknown')}
Moran's I: 0.73
P-value: 0.0012
Z-score: 3.24

Interpretation: Strong positive spatial autocorrelation detected.
The gene shows significant spatial clustering (p < 0.05).
""",
            "cell_type_deconvolution": """Cell Type Deconvolution Results (MOCK DATA)

Estimated cell type proportions:
- T cells: 35.2%
- B cells: 18.7%
- Macrophages: 22.1%
- Endothelial cells: 12.3%
- Fibroblasts: 11.7%

RMSE: 0.089
Correlation with reference: 0.92
""",
            "halla_analysis": f"""HAllA Multi-Omics Association Analysis (MOCK DATA)

Dataset 1: {arguments.get('dataset1_path', 'Unknown')}
Dataset 2: {arguments.get('dataset2_path', 'Unknown')}

Significant associations found: 127
Top associations:
1. Gene_A ↔ Metabolite_X (p=1.2e-8, r=0.84)
2. Gene_B ↔ Protein_Y (p=3.4e-7, r=0.76)
3. Gene_C ↔ Lipid_Z (p=8.9e-6, r=0.68)

Results saved to: {arguments.get('output_dir', 'output/')}/halla_results/
""",
            "stouffer_meta_analysis": """Stouffer Meta-Analysis Results (MOCK DATA)

Combined p-values from 3 studies
Method: Stouffer's Z-score

Top significant features (FDR < 0.05):
1. Feature_1: Z=4.82, p=1.4e-6, q=0.001
2. Feature_2: Z=4.21, p=2.5e-5, q=0.008
3. Feature_3: Z=3.94, p=8.1e-5, q=0.018

Total features tested: 15,234
Significant after FDR correction: 89
""",
            "validate_fastq": f"""FASTQ Validation Results (MOCK DATA)

File: {arguments.get('fastq_path', 'Unknown')}

✓ Format: Valid FASTQ
✓ Read count: 1,234,567
✓ Average quality score: 35.2 (Excellent)
✓ GC content: 42.3%
✓ Sequence length: 150bp (uniform)

No issues detected. File is ready for downstream analysis.
""",
            "calculate_fidelity": f"""Quantum Cell Type Fidelity Analysis (MOCK DATA)

Cell Type A: {arguments.get('cell_type_a', 'Unknown')}
Cell Type B: {arguments.get('cell_type_b', 'Unknown')}

Quantum Fidelity: F = 0.8734
Trace Distance: 0.1266
Hilbert-Schmidt Distance: 0.1891

Interpretation: High fidelity (F > 0.8) indicates strong similarity
between the quantum states of these cell types. They share substantial
gene expression patterns despite being classified as distinct types.

This suggests potential:
- Common developmental origin
- Similar functional programs
- Transition state relationship
""",
            "quantum_state_preparation": """Quantum State Preparation (MOCK DATA)

Input dimension: {dim}
Normalization: {norm}

✓ State vector prepared: |ψ⟩
✓ Normalization verified: ⟨ψ|ψ⟩ = 1.0000
✓ Entanglement entropy: 2.34 bits

Quantum state successfully encoded from gene expression profile.
Ready for fidelity calculations.
""".format(
                dim=len(arguments.get('expression_vector', [])),
                norm=arguments.get('normalization', 'l2')
            )
        }

        # Return specific mock response or generic one
        if tool_name in responses:
            return responses[tool_name]
        else:
            return f"""[MOCK] Tool '{tool_name}' executed successfully on server '{server_name}'

Arguments received:
{self._format_args(arguments)}

This is a mock response for local development.
Set USE_MOCK_MCP=false to use real MCP servers.
"""

    def _format_args(self, arguments: Dict) -> str:
        """Format arguments for display."""
        lines = []
        for key, value in arguments.items():
            lines.append(f"  - {key}: {value}")
        return "\n".join(lines) if lines else "  (none)"

    def get_all_tools(self) -> List[Dict]:
        """Get all tools from all connected servers.

        Returns:
            List of all available tools with server context
        """
        all_tools = []
        for server_name, tools in self.tools_cache.items():
            for tool in tools:
                # Add server name to tool for routing
                tool_with_server = tool.copy()
                tool_with_server["server_name"] = server_name
                all_tools.append(tool_with_server)

        return all_tools

    def convert_tools_to_gemini_format(self) -> List[Dict]:
        """Convert MCP tools to Gemini function declarations.

        Returns:
            List of Gemini function declarations
        """
        gemini_tools = []
        self._gemini_name_map = {}

        for tool in self.get_all_tools():
            # Gemini function names must match [a-zA-Z_][a-zA-Z0-9_]*
            safe_name = f"{tool['server_name']}_{tool['name']}".replace("-", "_")

            # Store mapping back to original names for tool execution
            self._gemini_name_map[safe_name] = (tool['server_name'], tool['name'])

            function_declaration = {
                "name": safe_name,
                "description": tool['description'],
                "parameters": tool['input_schema']
            }

            gemini_tools.append(function_declaration)

        return gemini_tools

    def resolve_gemini_tool_name(self, gemini_name: str) -> tuple:
        """Resolve a Gemini function name back to (server_name, tool_name).

        Args:
            gemini_name: The Gemini-safe function name

        Returns:
            Tuple of (server_name, tool_name) with original names
        """
        if gemini_name in self._gemini_name_map:
            return self._gemini_name_map[gemini_name]
        # Fallback to split
        parts = gemini_name.split("_", 1)
        return (parts[0], parts[1]) if len(parts) == 2 else (gemini_name, gemini_name)

    def convert_tools_to_claude_format(self) -> List[Dict]:
        """Convert MCP tools to Claude tool format.

        Returns:
            List of Claude tool declarations
        """
        claude_tools = []

        for tool in self.get_all_tools():
            # Use same format as real MCPClientManager
            tool_declaration = {
                "name": f"{tool['server_name']}_{tool['name']}",
                "description": tool['description'],
                "input_schema": tool['input_schema']
            }

            claude_tools.append(tool_declaration)

        return claude_tools


# Alias for easier imports
MCPClientManager = MockMCPClientManager
