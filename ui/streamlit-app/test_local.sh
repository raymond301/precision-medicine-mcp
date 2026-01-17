#!/bin/bash
# Test the Streamlit app locally before deploying

set -e

echo "=========================================="
echo "Testing Streamlit App Locally"
echo "=========================================="
echo ""

# Check if .env exists
if [ ! -f ".env" ]; then
    echo "❌ Error: .env file not found"
    echo "Please create .env from .env.example:"
    echo "  cp .env.example .env"
    echo "  # Edit .env and set your ANTHROPIC_API_KEY"
    exit 1
fi

# Load environment variables
echo "✅ Loading environment variables from .env"
export $(grep -v '^#' .env | xargs)

# Verify API key is set
if [ -z "$ANTHROPIC_API_KEY" ]; then
    echo "❌ Error: ANTHROPIC_API_KEY not set in .env"
    exit 1
fi

# Verify ENVIRONMENT is set to development
if [ "$ENVIRONMENT" != "development" ]; then
    echo "⚠️  Warning: ENVIRONMENT is set to '$ENVIRONMENT'"
    echo "   For local testing, it should be 'development'"
    echo "   This will bypass Azure AD SSO requirement"
fi

echo "✅ API Key configured"
echo "✅ Environment mode: $ENVIRONMENT"
echo ""
echo "🚀 Starting Streamlit app on http://localhost:8501"
echo ""
echo "📝 To test the orchestration trace feature:"
echo "   1. Enable 'Show trace for responses' in the sidebar"
echo "   2. Run a query that uses MCP servers"
echo "   3. Verify the trace appears below the response"
echo ""
echo "Press Ctrl+C to stop the server"
echo ""

# Run Streamlit
streamlit run app.py
