#!/usr/bin/env python3
"""
Simple functional test for all 9 GCP-deployed MCP servers.
Tests that each server responds and lists available tools.

Usage:
    export ANTHROPIC_API_KEY=your_key_here
    python tests/integration/test_all_gcp_servers.py
"""

import anthropic
import os
import sys
from pathlib import Path

# Color codes
GREEN = '\033[92m'
RED = '\033[91m'
BLUE = '\033[94m'
YELLOW = '\033[93m'
NC = '\033[0m'

SERVERS = {
    "fgbio": "https://mcp-fgbio-305650208648.us-central1.run.app/sse",
    "multiomics": "https://mcp-multiomics-305650208648.us-central1.run.app/sse",
    "spatialtools": "https://mcp-spatialtools-305650208648.us-central1.run.app/sse",
    "tcga": "https://mcp-tcga-305650208648.us-central1.run.app/sse",
    "openimagedata": "https://mcp-openimagedata-305650208648.us-central1.run.app/sse",
    "seqera": "https://mcp-seqera-305650208648.us-central1.run.app/sse",
    "huggingface": "https://mcp-huggingface-305650208648.us-central1.run.app/sse",
    "deepcell": "https://mcp-deepcell-305650208648.us-central1.run.app/sse",
    "mockepic": "https://mcp-mockepic-305650208648.us-central1.run.app/sse",
}

def test_server(name, url):
    """Test a single MCP server."""
    print(f"{BLUE}Testing {name:15}{NC} ", end="", flush=True)

    try:
        client = anthropic.Anthropic()

        response = client.beta.messages.create(
            model="claude-sonnet-4-5",
            max_tokens=512,
            messages=[{
                "role": "user",
                "content": f"List the available tools from the {name} MCP server in 1-2 sentences."
            }],
            mcp_servers=[{
                "type": "url",
                "url": url,
                "name": name,
            }],
            tools=[{"type": "mcp_toolset", "mcp_server_name": name}],
            betas=["mcp-client-2025-11-20"]
        )

        if response.content:
            print(f"{GREEN}✓ PASS{NC}")
            # Print first 100 chars of response
            text = response.content[0].text if hasattr(response.content[0], 'text') else str(response.content[0])
            print(f"        {text[:100]}...")
            return True
        else:
            print(f"{RED}✗ FAIL (no response){NC}")
            return False

    except Exception as e:
        print(f"{RED}✗ FAIL{NC}")
        print(f"        Error: {str(e)[:80]}")
        return False

def main():
    """Test all servers."""
    print(f"\n{BLUE}{'='*70}{NC}")
    print(f"{BLUE}GCP MCP Server Functional Test{NC}")
    print(f"{BLUE}{'='*70}{NC}\n")

    if not os.getenv("ANTHROPIC_API_KEY"):
        print(f"{RED}Error: ANTHROPIC_API_KEY not set{NC}")
        print(f"Set it with: export ANTHROPIC_API_KEY=your_key_here")
        sys.exit(1)

    print(f"{YELLOW}Testing {len(SERVERS)} MCP servers...{NC}\n")

    results = {}
    for name, url in SERVERS.items():
        results[name] = test_server(name, url)
        print()  # Empty line between tests

    # Summary
    passed = sum(results.values())
    failed = len(results) - passed

    print(f"{BLUE}{'='*70}{NC}")
    print(f"{BLUE}Summary{NC}")
    print(f"{BLUE}{'='*70}{NC}\n")

    for name, success in results.items():
        status = f"{GREEN}✓ PASS{NC}" if success else f"{RED}✗ FAIL{NC}"
        print(f"  {name:20} {status}")

    print(f"\n{BLUE}Results:{NC} {GREEN}{passed}{NC}/{len(results)} servers passed")

    if failed > 0:
        print(f"         {RED}{failed}{NC}/{len(results)} servers failed\n")
        return 1
    else:
        print(f"\n{GREEN}✅ All servers passed!{NC}\n")
        return 0

if __name__ == "__main__":
    sys.exit(main())
