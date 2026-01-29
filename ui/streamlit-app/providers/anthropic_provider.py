"""Anthropic (Claude) provider implementation with MCP tool calling support.

Uses manual MCP client management similar to Gemini provider.
"""

import os
import asyncio
import sys
from typing import List, Dict, Optional, Any
import anthropic

from .base import LLMProvider, ChatMessage, ChatResponse, UsageInfo


class AnthropicProvider(LLMProvider):
    """Anthropic Claude provider with MCP support."""

    def __init__(self, api_key: Optional[str] = None):
        """Initialize Anthropic provider.

        Args:
            api_key: Anthropic API key (defaults to ANTHROPIC_API_KEY env var)
        """
        self.api_key = api_key or os.getenv("ANTHROPIC_API_KEY")
        if not self.api_key:
            raise ValueError("ANTHROPIC_API_KEY not found in environment")

        self.client = anthropic.Anthropic(api_key=self.api_key)

    def send_message(
        self,
        messages: List[ChatMessage],
        mcp_servers: List[Dict],
        model: str = "claude-sonnet-4-20250514",
        max_tokens: int = 4096,
        temperature: float = 1.0,
        uploaded_files: Optional[Dict] = None
    ) -> ChatResponse:
        """Send message to Claude with MCP tool calling.

        Args:
            messages: Chat message history
            mcp_servers: MCP server configurations
            model: Claude model to use
            max_tokens: Maximum response tokens
            temperature: Sampling temperature
            uploaded_files: Dict of uploaded files

        Returns:
            ChatResponse with standardized format
        """
        # Run async tool calling in sync context
        return asyncio.run(
            self._send_message_async(
                messages=messages,
                mcp_servers=mcp_servers,
                model=model,
                max_tokens=max_tokens,
                temperature=temperature,
                uploaded_files=uploaded_files
            )
        )

    async def _send_message_async(
        self,
        messages: List[ChatMessage],
        mcp_servers: List[Dict],
        model: str,
        max_tokens: int,
        temperature: float,
        uploaded_files: Optional[Dict]
    ) -> ChatResponse:
        """Async implementation of message sending with tool calling.

        This implements the agentic loop:
        1. Connect to MCP servers and discover tools
        2. Send message to Claude with tool declarations
        3. If Claude calls tools, execute them via MCP
        4. Feed results back to Claude
        5. Repeat until Claude responds without tool calls
        """
        try:
            # Import MCP client manager
            parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            if parent_dir not in sys.path:
                sys.path.insert(0, parent_dir)

            from utils.mcp_client import MCPClientManager

            # Connect to MCP servers using context manager
            print(f"DEBUG: Connecting to {len(mcp_servers)} MCP servers", file=sys.stderr)

            async with MCPClientManager(mcp_servers) as mcp_manager:
                # Get Claude-formatted tools
                claude_tools = mcp_manager.convert_tools_to_claude_format()
                print(f"DEBUG: Discovered {len(claude_tools)} total tools", file=sys.stderr)

                # Convert messages to Claude format
                api_messages = [
                    {"role": msg.role, "content": msg.content}
                    for msg in messages
                ]

                # Build system prompt
                system_prompt = self._build_system_prompt(mcp_servers, uploaded_files)

                # Agentic loop: keep calling until no more tool calls
                max_iterations = 30
                iteration = 0
                conversation_history = api_messages.copy()
                all_tool_calls = []  # Track all tool calls for trace

                while iteration < max_iterations:
                    iteration += 1
                    print(f"DEBUG: Iteration {iteration}", file=sys.stderr)

                    # Call Claude with current conversation and tools
                    response = self.client.messages.create(
                        model=model,
                        max_tokens=max_tokens,
                        temperature=temperature,
                        messages=conversation_history,
                        tools=claude_tools,
                        system=system_prompt
                    )

                    # Check if response contains tool calls
                    tool_calls = self._extract_tool_calls(response)

                    # Debug: log response structure
                    print(f"DEBUG: Stop reason: {response.stop_reason}", file=sys.stderr)
                    print(f"DEBUG: Content blocks: {len(response.content)}", file=sys.stderr)

                    if not tool_calls:
                        # No more tool calls - we're done
                        print(f"DEBUG: No tool calls, finishing after {iteration} iterations", file=sys.stderr)

                        # Extract final response
                        content = self._format_response(response)
                        usage = self._get_usage_info(response)

                        return ChatResponse(
                            content=content,
                            usage=usage,
                            raw_response=response,
                            tool_calls_metadata=all_tool_calls if all_tool_calls else None
                        )

                    # Execute tool calls
                    print(f"DEBUG: Executing {len(tool_calls)} tool calls", file=sys.stderr)
                    tool_results = await self._execute_tool_calls(mcp_manager, tool_calls)

                    # Store tool calls and results for trace
                    for i, tc in enumerate(tool_calls):
                        all_tool_calls.append({
                            "server_name": tc["name"].split("_", 1)[0] if "_" in tc["name"] else "unknown",
                            "tool_name": tc["name"].split("_", 1)[1] if "_" in tc["name"] else tc["name"],
                            "input": tc["input"],
                            "result": tool_results[i]["content"] if i < len(tool_results) else None
                        })

                    # Add assistant message with tool calls
                    conversation_history.append({
                        "role": "assistant",
                        "content": response.content
                    })

                    # Add tool results
                    conversation_history.append({
                        "role": "user",
                        "content": tool_results
                    })

                # Max iterations reached
                print(f"WARNING: Max iterations ({max_iterations}) reached", file=sys.stderr)
                return ChatResponse(
                    content="Error: Maximum tool calling iterations reached. Please try a simpler query.",
                    usage=None,
                    raw_response=None
                )

        except anthropic.APIError as e:
            import traceback
            traceback.print_exc()
            raise Exception(f"Claude API error: {str(e)}")
        except Exception as e:
            import traceback
            traceback.print_exc()
            raise Exception(f"Unexpected error: {str(e)}")

    def _extract_tool_calls(self, response: anthropic.types.Message) -> List[Dict]:
        """Extract tool calls from Claude response.

        Args:
            response: Claude API response

        Returns:
            List of tool calls with name and input
        """
        tool_calls = []

        if response.stop_reason == "tool_use":
            for block in response.content:
                if block.type == "tool_use":
                    tool_calls.append({
                        "id": block.id,
                        "name": block.name,
                        "input": block.input
                    })

        return tool_calls

    async def _execute_tool_calls(self, mcp_manager, tool_calls: List[Dict]) -> List[Dict]:
        """Execute MCP tool calls.

        Args:
            mcp_manager: MCP client manager instance
            tool_calls: List of tool calls from Claude

        Returns:
            List of tool results in Claude format
        """
        results = []

        for tool_call in tool_calls:
            # Parse server name and tool name (format: servername_toolname)
            full_name = tool_call["name"]
            parts = full_name.split("_", 1)

            if len(parts) != 2:
                print(f"ERROR: Invalid tool name format: {full_name}", file=sys.stderr)
                results.append({
                    "type": "tool_result",
                    "tool_use_id": tool_call["id"],
                    "content": f"Error: Invalid tool name format: {full_name}"
                })
                continue

            server_name, tool_name = parts

            # Call the tool via MCP manager
            try:
                result = await mcp_manager.call_tool(
                    server_name=server_name,
                    tool_name=tool_name,
                    arguments=tool_call["input"]
                )

                # Format result for Claude
                # Extract text from content blocks
                content_text = ""
                if isinstance(result.get("content"), list):
                    for item in result["content"]:
                        if isinstance(item, dict) and item.get("type") == "text":
                            content_text += item.get("text", "")

                results.append({
                    "type": "tool_result",
                    "tool_use_id": tool_call["id"],
                    "content": content_text or str(result.get("content", ""))
                })

            except Exception as e:
                print(f"ERROR: Tool execution failed: {e}", file=sys.stderr)
                results.append({
                    "type": "tool_result",
                    "tool_use_id": tool_call["id"],
                    "content": f"Error: {str(e)}"
                })

        return results

    def get_provider_name(self) -> str:
        """Get provider name."""
        return "Claude"

    def get_model_display_name(self, model: str) -> str:
        """Get human-readable model name."""
        model_names = {
            "claude-sonnet-4-20250514": "Claude Sonnet 4.5",
            "claude-sonnet-4-5": "Claude Sonnet 4.5",
            "claude-opus-4-5": "Claude Opus 4.5",
            "claude-haiku-4": "Claude Haiku 4"
        }
        return model_names.get(model, model)

    def format_mcp_servers(self, server_configs: List[Dict]) -> List[Dict]:
        """Format MCP servers (pass-through for manual MCP management).

        Args:
            server_configs: Standard MCP server configs

        Returns:
            Server configs (unchanged)
        """
        return server_configs

    def is_available(self) -> bool:
        """Check if Claude provider is available."""
        return bool(self.api_key)

    def _build_system_prompt(
        self,
        mcp_servers: List[Dict],
        uploaded_files: Optional[Dict] = None
    ) -> str:
        """Build system prompt with MCP server and file information."""
        server_descriptions = [
            f"- {server['name']}: Available for bioinformatics analysis"
            for server in mcp_servers
        ]

        system_prompt = f"""You are a bioinformatics assistant with access to specialized analysis tools via MCP servers.

IMPORTANT: When users ask you to perform analysis, you MUST use the appropriate MCP tool rather than providing theoretical responses.

Available MCP servers:
{chr(10).join(server_descriptions)}

Guidelines:
- When asked to perform pathway enrichment, spatial analysis, or other bioinformatics tasks, USE THE TOOLS
- Call the appropriate MCP tool to get real analysis results
- After receiving tool results, interpret and explain them to the user
- Provide actionable insights based on the actual analysis output
- If you don't have access to a required tool, let the user know

Do not simulate or describe what an analysis would show - actually perform it using the available tools."""

        # Add uploaded file information if present
        if uploaded_files and len(uploaded_files) > 0:
            file_descriptions = []
            for filename, file_info in uploaded_files.items():
                metadata = file_info['metadata']
                path = file_info['path']
                original_name = file_info.get('original_name', filename)
                is_gcs = file_info.get('source') == 'gcs'

                if is_gcs:
                    file_desc = f"""- **{original_name}** (GCS)
  - GCS URI: {path}
  - Bucket: {metadata.get('bucket', 'unknown')}
  - Path: {metadata.get('blob_path', 'unknown')}
  - Type: {metadata['extension']}"""
                else:
                    file_desc = f"""- **{original_name}**
  - File path: {path}
  - Type: {metadata['extension']}
  - Size: {metadata['size_mb']:.2f} MB"""

                if not metadata.get('is_binary', False):
                    file_desc += f"\n  - Lines: {metadata.get('line_count', 'N/A')}"

                    # Include content preview for small files
                    content = metadata.get('_gcs_content')
                    if not content and not is_gcs and metadata.get('size_bytes', 0) < 50000:
                        try:
                            with open(path, 'r') as f:
                                content = f.read()
                        except Exception:
                            pass

                    if content:
                        if len(content) > 10000:
                            content = content[:10000] + "\n\n[... truncated ...]"
                        file_desc += f"\n  - Content preview:\n```\n{content}\n```"

                file_descriptions.append(file_desc)

            system_prompt += f"""

UPLOADED FILES AVAILABLE:
The user has uploaded {len(uploaded_files)} file(s) for analysis.

{chr(10).join(file_descriptions)}

FILE ACCESS INSTRUCTIONS - READ CAREFULLY:

For GCS files (gs://bucket/path/file):
✅ MCP servers CAN and WILL access GCS URIs directly - this works and is tested
✅ Pass the full GCS URI (e.g., "gs://sample-inputs-patientone/patient-data/PAT001-OVC-2025/multiomics/pdx_rna_seq.csv") as the file parameter to MCP tools
✅ The multiomics server has fsspec and gcsfs installed and can read from GCS
✅ DO NOT apologize about file access - just call the tool with the GCS URI

For files with content preview:
- You can analyze the content directly without calling tools

CRITICAL: When the user asks you to analyze GCS files, you MUST call the MCP tools with the GCS URIs shown above. Do not provide theoretical responses or apologize about file access. The infrastructure is configured and working."""

        return system_prompt

    def _format_response(self, response: anthropic.types.Message) -> str:
        """Format Claude response for display."""
        if not response.content:
            return "No response from Claude"

        parts = []
        for block in response.content:
            # Only extract text blocks, not tool_use blocks
            if block.type == "text":
                parts.append(block.text)

        return "\n\n".join(parts) if parts else "Empty response"

    def _get_usage_info(self, response: anthropic.types.Message) -> Optional[UsageInfo]:
        """Extract usage information from response."""
        if hasattr(response, 'usage') and response.usage:
            return UsageInfo(
                input_tokens=response.usage.input_tokens,
                output_tokens=response.usage.output_tokens,
                total_tokens=response.usage.input_tokens + response.usage.output_tokens
            )
        return None
