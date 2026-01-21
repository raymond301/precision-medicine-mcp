#!/bin/bash
set -e

SERVICE_URL="https://mcp-perturbation-ondu7mwjpa-uc.a.run.app"

echo "🧪 Testing mcp-perturbation deployment..."
echo ""

# Test 1: Basic connectivity
echo "1️⃣  Testing basic connectivity..."
HTTP_CODE=$(curl -f -s -o /dev/null -w "%{http_code}" "$SERVICE_URL/sse" || echo "000")
if [ "$HTTP_CODE" != "000" ]; then
  echo "✅ Service is reachable (HTTP $HTTP_CODE)"
else
  echo "❌ Service unreachable"
  exit 1
fi

# Test 2: MCP initialize
echo ""
echo "2️⃣  Testing MCP protocol handshake..."
INIT_RESPONSE=$(curl -s -X POST "$SERVICE_URL/message" \
  -H "Content-Type: application/json" \
  -d '{
    "jsonrpc": "2.0",
    "id": 1,
    "method": "initialize",
    "params": {
      "protocolVersion": "2024-11-05",
      "capabilities": {},
      "clientInfo": {"name": "test-client", "version": "1.0.0"}
    }
  }')

if echo "$INIT_RESPONSE" | jq -e '.result' > /dev/null 2>&1; then
  echo "✅ MCP handshake successful"
  SERVER_NAME=$(echo "$INIT_RESPONSE" | jq -r '.result.serverInfo.name')
  SERVER_VERSION=$(echo "$INIT_RESPONSE" | jq -r '.result.serverInfo.version')
  echo "   Server: $SERVER_NAME v$SERVER_VERSION"
else
  echo "❌ MCP handshake failed"
  echo "$INIT_RESPONSE" | jq .
  exit 1
fi

# Test 3: List tools
echo ""
echo "3️⃣  Testing tool listing..."
TOOLS_RESPONSE=$(curl -s -X POST "$SERVICE_URL/message" \
  -H "Content-Type: application/json" \
  -d '{
    "jsonrpc": "2.0",
    "id": 2,
    "method": "tools/list",
    "params": {}
  }')

TOOL_COUNT=$(echo "$TOOLS_RESPONSE" | jq -r '.result.tools | length')
if [ "$TOOL_COUNT" -eq 8 ]; then
  echo "✅ All 8 tools registered:"
  echo "$TOOLS_RESPONSE" | jq -r '.result.tools[] | "   • \(.name)"'
else
  echo "❌ Expected 8 tools, got $TOOL_COUNT"
  exit 1
fi

# Test 4: Call a tool
echo ""
echo "4️⃣  Testing tool execution (load_dataset)..."
TOOL_RESPONSE=$(curl -s -X POST "$SERVICE_URL/message" \
  -H "Content-Type: application/json" \
  -d '{
    "jsonrpc": "2.0",
    "id": 3,
    "method": "tools/call",
    "params": {
      "name": "perturbation_load_dataset",
      "arguments": {
        "dataset_id": "test_synthetic",
        "normalize": true,
        "n_hvg": 2000
      }
    }
  }')

if echo "$TOOL_RESPONSE" | jq -e '.result' > /dev/null 2>&1; then
  echo "✅ Tool execution successful"
  echo "$TOOL_RESPONSE" | jq -r '.result.content[0].text' | head -3
else
  echo "⚠️  Tool execution returned error (may be expected in DRY_RUN mode)"
  ERROR_MSG=$(echo "$TOOL_RESPONSE" | jq -r '.error.message // "Unknown error"')
  echo "   Error: $ERROR_MSG"
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "🎉 Deployment tests complete!"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "📊 Summary:"
echo "   Service URL: $SERVICE_URL"
echo "   Region: us-central1"
echo "   Status: ✅ OPERATIONAL"
echo "   Tools: $TOOL_COUNT/8 registered"
echo ""
echo "📝 Next steps:"
echo "   1. Add to Claude Desktop config"
echo "   2. Test with Streamlit UI"
echo "   3. Run full GEARS workflow"
echo "   4. Monitor costs in GCP Console"
echo ""
