"""MCP SSE Client Manager for Cloud Run MCP servers.

Manages SSE connections to remote MCP servers and provides tool discovery/calling.
"""

import os
import sys
from typing import List, Dict, Any, Optional
from contextlib import asynccontextmanager
import asyncio

try:
    from mcp import ClientSession
    from mcp.client.sse import sse_client
    MCP_AVAILABLE = True
except ImportError:
    MCP_AVAILABLE = False
    ClientSession = None
    sse_client = None

try:
    import google.auth
    from google.auth.transport.requests import Request
    import google.oauth2.id_token
    GOOGLE_AUTH_AVAILABLE = True
except ImportError:
    GOOGLE_AUTH_AVAILABLE = False


class MCPClientManager:
    """Manages SSE connections to Cloud Run MCP servers.

    This manager handles the lifecycle of MCP connections during a tool calling session.
    Connections are established when entering the context and closed when exiting.
    """

    def __init__(self, server_configs: List[Dict]):
        """Initialize MCP client manager.

        Args:
            server_configs: List of server configurations
                [{"name": "server1", "url": "https://..."}]
        """
        if not MCP_AVAILABLE:
            raise ImportError(
                "mcp package not installed. "
                "Install with: pip install mcp>=1.0.0"
            )

        self.server_configs = server_configs
        self.sessions: Dict[str, ClientSession] = {}
        self.tools_cache: Dict[str, List[Dict]] = {}
        self._context_managers = []

    async def __aenter__(self):
        """Enter async context - establish all MCP connections."""
        print(f"DEBUG: Establishing connections to {len(self.server_configs)} MCP servers", file=sys.stderr)

        for server in self.server_configs:
            try:
                await self._connect_to_server(server["name"], server["url"])
            except Exception as e:
                print(f"ERROR: Failed to connect to {server['name']}: {e}", file=sys.stderr)

        print(f"DEBUG: Connected to {len(self.sessions)} servers", file=sys.stderr)
        return self

    async def __aexit__(self, exc_type, exc_val, exc_tb):
        """Exit async context - close all MCP connections."""
        print(f"DEBUG: Closing all MCP connections", file=sys.stderr)

        # Close all context managers in reverse order
        for cm in reversed(self._context_managers):
            try:
                await cm.__aexit__(None, None, None)
            except Exception as e:
                print(f"WARNING: Error closing connection: {e}", file=sys.stderr)

        self._context_managers.clear()
        self.sessions.clear()
        self.tools_cache.clear()

    async def _connect_to_server(self, server_name: str, server_url: str):
        """Connect to an MCP server via SSE.

        Args:
            server_name: Name of the MCP server
            server_url: SSE endpoint URL (should end with /sse)
        """
        # Ensure URL ends with /sse
        if not server_url.endswith("/sse"):
            server_url = server_url + "/sse"

        print(f"DEBUG: Connecting to {server_name} at {server_url}", file=sys.stderr)

        # Get authentication headers for Cloud Run
        headers = self._get_auth_headers(server_url)

        # Create SSE client context manager
        sse_cm = sse_client(url=server_url, headers=headers)
        read, write = await sse_cm.__aenter__()
        self._context_managers.append(sse_cm)

        # Create client session context manager
        session_cm = ClientSession(read, write)
        session = await session_cm.__aenter__()
        self._context_managers.append(session_cm)

        # Initialize the session
        await session.initialize()

        # Cache the session
        self.sessions[server_name] = session

        # Discover and cache tools
        tools = await self._discover_tools(session)
        self.tools_cache[server_name] = tools

        print(f"DEBUG: Connected to {server_name}, found {len(tools)} tools", file=sys.stderr)

    def _get_auth_headers(self, server_url: str) -> Dict[str, str]:
        """Get authentication headers for Cloud Run service.

        Args:
            server_url: The Cloud Run service URL

        Returns:
            Headers dict with Authorization token
        """
        headers = {}

        # Only add auth for Cloud Run URLs (*.run.app)
        if ".run.app" in server_url and GOOGLE_AUTH_AVAILABLE:
            try:
                # Get ID token for Cloud Run authentication
                auth_req = Request()
                target_audience = server_url.split("/sse")[0]  # Remove /sse suffix

                id_token = google.oauth2.id_token.fetch_id_token(auth_req, target_audience)
                headers["Authorization"] = f"Bearer {id_token}"

                print(f"DEBUG: Added Cloud Run auth token for {target_audience}", file=sys.stderr)
            except Exception as e:
                print(f"WARNING: Could not get Cloud Run auth token: {e}", file=sys.stderr)

        return headers

    async def _discover_tools(self, session: ClientSession) -> List[Dict]:
        """Discover available tools from an MCP server.

        Args:
            session: Connected MCP session

        Returns:
            List of tool definitions
        """
        try:
            # List available tools
            tools_result = await session.list_tools()

            tools = []
            for tool in tools_result.tools:
                tools.append({
                    "name": tool.name,
                    "description": tool.description,
                    "input_schema": tool.inputSchema
                })

            return tools
        except Exception as e:
            print(f"ERROR: Failed to discover tools: {e}", file=sys.stderr)
            return []

    async def call_tool(
        self,
        server_name: str,
        tool_name: str,
        arguments: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Call a tool on an MCP server.

        Args:
            server_name: Name of the MCP server
            tool_name: Name of the tool to call
            arguments: Tool arguments

        Returns:
            Tool result
        """
        session = self.sessions.get(server_name)
        if not session:
            raise ValueError(f"Not connected to server: {server_name}")

        print(f"DEBUG: Calling tool {tool_name} on {server_name}", file=sys.stderr)
        print(f"DEBUG: Arguments: {arguments}", file=sys.stderr)

        try:
            # Call the tool
            result = await session.call_tool(tool_name, arguments)

            print(f"DEBUG: Tool result: {result}", file=sys.stderr)

            return {
                "content": result.content,
                "isError": result.isError if hasattr(result, 'isError') else False
            }
        except Exception as e:
            print(f"ERROR: Tool call failed: {e}", file=sys.stderr)
            import traceback
            traceback.print_exc()
            return {
                "content": [{"type": "text", "text": f"Error: {str(e)}"}],
                "isError": True
            }

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

        for tool in self.get_all_tools():
            # Clean the input schema for Gemini
            cleaned_schema = self._clean_schema_for_gemini(tool['input_schema'])

            # Convert MCP tool schema to Gemini function declaration
            function_declaration = {
                "name": f"{tool['server_name']}_{tool['name']}",  # Prefix with server name
                "description": tool['description'],
                "parameters": cleaned_schema
            }

            gemini_tools.append(function_declaration)

        return gemini_tools

    def _clean_schema_for_gemini(self, schema: Dict[str, Any]) -> Dict[str, Any]:
        """Clean JSON schema to remove properties Gemini doesn't support.

        Gemini's function calling doesn't support:
        - additional_properties
        - any_of, one_of, all_of
        - $ref, $schema, $id
        - And other advanced JSON schema features

        Args:
            schema: Original JSON schema

        Returns:
            Cleaned schema compatible with Gemini
        """
        if not isinstance(schema, dict):
            return schema

        cleaned = {}

        # List of keys to explicitly exclude
        excluded_keys = {
            "additional_properties", "additionalProperties",
            "any_of", "anyOf", "one_of", "oneOf", "all_of", "allOf",
            "$ref", "$schema", "$id", "definitions"
        }

        # List of allowed keys for Gemini
        allowed_keys = {
            "type", "properties", "required", "description",
            "items", "enum", "format", "default", "minimum",
            "maximum", "minLength", "maxLength", "pattern",
            "minItems", "maxItems", "title"
        }

        for key, value in schema.items():
            # Skip excluded keys
            if key in excluded_keys:
                continue

            # Only keep allowed keys (or pass through if we're being permissive)
            # For now, let's exclude known bad keys but allow others

            # Recursively clean nested schemas
            if key == "properties" and isinstance(value, dict):
                cleaned[key] = {
                    prop_name: self._clean_schema_for_gemini(prop_schema)
                    for prop_name, prop_schema in value.items()
                }
            elif key == "items" and isinstance(value, dict):
                cleaned[key] = self._clean_schema_for_gemini(value)
            elif isinstance(value, dict):
                # Recursively clean any nested dict
                cleaned[key] = self._clean_schema_for_gemini(value)
            elif isinstance(value, list):
                # Clean items in lists
                cleaned[key] = [
                    self._clean_schema_for_gemini(item) if isinstance(item, dict) else item
                    for item in value
                ]
            else:
                # Keep primitive values
                cleaned[key] = value

        return cleaned
