#!/bin/bash
# Pre-deployment validation script for Streamlit MCP Chat
# Catches errors before expensive Cloud Run deployment

set -e  # Exit on any error

# Color codes for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "=========================================="
echo "  Streamlit MCP Chat - Validation"
echo "=========================================="
echo ""

# Change to script directory
cd "$(dirname "$0")"

# 1. Check Python version
echo "1. Checking Python version..."
PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}')
echo "   Found: Python ${PYTHON_VERSION}"

# Warn if Python 3.9 (EOL)
if [[ "$PYTHON_VERSION" == 3.9* ]]; then
    echo -e "   ${YELLOW}⚠️  Warning: Python 3.9 is past end-of-life${NC}"
    echo "   Cloud Run will use Python 3.11 (see Dockerfile)"
fi
echo -e "   ${GREEN}✓${NC} Python version check complete"
echo ""

# 2. Syntax check all Python files
echo "2. Running syntax checks..."
PYTHON_FILES=(
    "app.py"
    "providers/base.py"
    "providers/anthropic_provider.py"
    "providers/gemini_provider.py"
    "utils/mcp_client.py"
    "utils/mcp_config.py"
    "utils/audit_logger.py"
    "utils/trace_builder.py"
)

for file in "${PYTHON_FILES[@]}"; do
    if [ -f "$file" ]; then
        if python3 -m py_compile "$file" 2>/dev/null; then
            echo -e "   ${GREEN}✓${NC} $file"
        else
            echo -e "   ${RED}✗${NC} $file - SYNTAX ERROR"
            python3 -m py_compile "$file"  # Show error
            exit 1
        fi
    else
        echo -e "   ${YELLOW}⚠${NC}  $file - NOT FOUND"
    fi
done
echo ""

# 3. Check required dependencies
echo "3. Checking dependencies..."

# Check if mock mode is enabled
USE_MOCK_MCP=${USE_MOCK_MCP:-false}
if [ "$USE_MOCK_MCP" = "true" ]; then
    echo "   ${YELLOW}ℹ${NC}  USE_MOCK_MCP=true - MCP package not required"
fi

REQUIRED_PACKAGES=(
    "streamlit"
    "anthropic"
    "google.cloud.logging"
    "google.cloud.storage"
)

# Add MCP to required packages only if not using mock mode
if [ "$USE_MOCK_MCP" != "true" ]; then
    REQUIRED_PACKAGES+=("mcp")
fi

MISSING_PACKAGES=()
for package in "${REQUIRED_PACKAGES[@]}"; do
    if python3 -c "import ${package}" 2>/dev/null; then
        echo -e "   ${GREEN}✓${NC} ${package}"
    else
        echo -e "   ${RED}✗${NC} ${package} - NOT INSTALLED"
        MISSING_PACKAGES+=("$package")
    fi
done

# Check MCP separately if in mock mode
if [ "$USE_MOCK_MCP" = "true" ]; then
    if python3 -c "import mcp" 2>/dev/null; then
        echo -e "   ${GREEN}✓${NC} mcp (optional in mock mode)"
    else
        echo -e "   ${YELLOW}⚠${NC}  mcp (not installed - using mock mode)"
    fi
fi

# Check optional packages
OPTIONAL_PACKAGES=(
    "google.genai:google-genai"
    "pandas:pandas"
)

for entry in "${OPTIONAL_PACKAGES[@]}"; do
    IFS=':' read -r import_name pip_name <<< "$entry"
    if python3 -c "import ${import_name}" 2>/dev/null; then
        echo -e "   ${GREEN}✓${NC} ${pip_name} (optional)"
    else
        echo -e "   ${YELLOW}⚠${NC}  ${pip_name} (optional) - NOT INSTALLED"
    fi
done

if [ ${#MISSING_PACKAGES[@]} -gt 0 ]; then
    echo ""
    echo -e "${RED}Missing required packages!${NC}"
    echo "Install with: python3 -m pip install -r requirements.txt"
    exit 1
fi
echo ""

# 4. Check imports
echo "4. Checking module imports..."
if python3 -c "
import sys
sys.path.insert(0, '.')

# Try importing main modules
try:
    from providers.anthropic_provider import AnthropicProvider
    print('   ✓ AnthropicProvider')
except Exception as e:
    print(f'   ✗ AnthropicProvider: {e}')
    sys.exit(1)

try:
    from providers.gemini_provider import GeminiProvider
    print('   ✓ GeminiProvider')
except Exception as e:
    print(f'   ✗ GeminiProvider: {e}')
    sys.exit(1)

try:
    from utils.mcp_client import MCPClientManager
    print('   ✓ MCPClientManager')
except Exception as e:
    print(f'   ✗ MCPClientManager: {e}')
    sys.exit(1)

try:
    from utils.audit_logger import AuditLogger
    print('   ✓ AuditLogger')
except Exception as e:
    print(f'   ✗ AuditLogger: {e}')
    sys.exit(1)

print('All imports successful')
" 2>&1; then
    echo -e "   ${GREEN}✓${NC} All modules import successfully"
else
    echo -e "   ${RED}✗${NC} Import errors detected"
    exit 1
fi
echo ""

# 5. Check environment variables
echo "5. Checking environment configuration..."
if [ -f ".env" ]; then
    echo -e "   ${GREEN}✓${NC} .env file found"

    # Check for API keys (without showing values)
    if grep -q "ANTHROPIC_API_KEY=" .env; then
        echo -e "   ${GREEN}✓${NC} ANTHROPIC_API_KEY configured in .env"
    else
        echo -e "   ${YELLOW}⚠${NC}  ANTHROPIC_API_KEY not in .env"
    fi

    if grep -q "GEMINI_API_KEY=" .env; then
        echo -e "   ${GREEN}✓${NC} GEMINI_API_KEY configured in .env"
    else
        echo -e "   ${YELLOW}⚠${NC}  GEMINI_API_KEY not in .env (Gemini provider will be unavailable)"
    fi
else
    echo -e "   ${YELLOW}⚠${NC}  .env file not found (using environment variables)"
fi

# Check if API keys are in environment
if [ -n "$ANTHROPIC_API_KEY" ]; then
    echo -e "   ${GREEN}✓${NC} ANTHROPIC_API_KEY set in environment"
else
    echo -e "   ${YELLOW}⚠${NC}  ANTHROPIC_API_KEY not in environment"
fi

if [ -n "$GEMINI_API_KEY" ]; then
    echo -e "   ${GREEN}✓${NC} GEMINI_API_KEY set in environment"
else
    echo -e "   ${YELLOW}⚠${NC}  GEMINI_API_KEY not in environment"
fi
echo ""

# 6. Check for common issues
echo "6. Checking for common issues..."

# Check for tabs in Python files
if grep -r "$(printf '\t')" *.py providers/*.py utils/*.py 2>/dev/null; then
    echo -e "   ${YELLOW}⚠${NC}  Tabs found in Python files (should use spaces)"
else
    echo -e "   ${GREEN}✓${NC} No tabs in Python files"
fi

# Check for trailing whitespace (can cause issues)
if grep -r " $" *.py providers/*.py utils/*.py 2>/dev/null | head -5; then
    echo -e "   ${YELLOW}⚠${NC}  Trailing whitespace found (minor issue)"
else
    echo -e "   ${GREEN}✓${NC} No trailing whitespace"
fi

# Check Dockerfile exists
if [ -f "Dockerfile" ]; then
    echo -e "   ${GREEN}✓${NC} Dockerfile exists"
else
    echo -e "   ${RED}✗${NC} Dockerfile missing"
    exit 1
fi

# Check requirements.txt exists
if [ -f "requirements.txt" ]; then
    echo -e "   ${GREEN}✓${NC} requirements.txt exists"
else
    echo -e "   ${RED}✗${NC} requirements.txt missing"
    exit 1
fi
echo ""

# 7. Check provider initialization (if API keys available)
echo "7. Testing provider initialization..."
if [ -n "$ANTHROPIC_API_KEY" ] && [ -n "$GEMINI_API_KEY" ]; then
    if python3 -c "
import sys
import os
sys.path.insert(0, '.')

from providers.anthropic_provider import AnthropicProvider
from providers.gemini_provider import GeminiProvider

try:
    claude = AnthropicProvider()
    print(f'   ✓ Claude provider initialized ({claude.get_provider_name()})')
except Exception as e:
    print(f'   ✗ Claude provider error: {e}')
    sys.exit(1)

try:
    gemini = GeminiProvider()
    print(f'   ✓ Gemini provider initialized ({gemini.get_provider_name()})')
except Exception as e:
    print(f'   ✗ Gemini provider error: {e}')
    sys.exit(1)
" 2>&1; then
        echo -e "   ${GREEN}✓${NC} All providers initialized successfully"
    else
        echo -e "   ${RED}✗${NC} Provider initialization failed"
        exit 1
    fi
else
    echo -e "   ${YELLOW}⚠${NC}  Skipping provider test (API keys not set)"
fi
echo ""

# Final summary
echo "=========================================="
echo -e "${GREEN}✅ Validation Complete!${NC}"
echo "=========================================="
echo ""
echo "The app is ready to deploy."
echo "Run: ./deploy.sh"
echo ""
