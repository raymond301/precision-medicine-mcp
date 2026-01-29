# Streamlit MCP Chat - Development Guide

Quick start guide for local development and testing.

## Prerequisites

- Python 3.9+ (Python 3.11 recommended for production)
- Anthropic API key
- Optional: Gemini API key

## Quick Start

### 1. Environment Setup

```bash
cd ui/streamlit-app

# Copy example environment file
cp .env.example .env

# Edit .env and add your API keys
nano .env  # or your favorite editor
```

### 2. Install Dependencies

```bash
# Install Python packages
python3 -m pip install -r requirements.txt
```

### 3. Validate Before Running

```bash
# Run validation script
./validate.sh
```

This checks:
- ✓ Python syntax
- ✓ Module imports
- ✓ Dependencies installed
- ✓ Environment configuration
- ✓ Provider initialization

### 4. Run Locally

**Option A: Mock Mode (Recommended for Development)**

```bash
# Set USE_MOCK_MCP=true in .env, then:
streamlit run app.py
```

Benefits:
- No Cloud Run connections needed
- Instant responses
- Works offline
- No MCP server costs

**Option B: Real MCP Servers**

```bash
# Set USE_MOCK_MCP=false in .env, then:
streamlit run app.py
```

Requirements:
- MCP servers deployed to Cloud Run
- Proper authentication configured
- Python 3.10+ (for mcp package)

## Development Workflow

### Fast Iteration Cycle

```bash
# 1. Make code changes
vim app.py

# 2. Validate (10 seconds)
./validate.sh

# 3. Test locally with mock (instant)
USE_MOCK_MCP=true streamlit run app.py

# 4. When ready, deploy
./deploy.sh
```

### Before Every Deployment

```bash
# Always validate before deploying
export ANTHROPIC_API_KEY=...
export GEMINI_API_KEY=...
./validate.sh && ./deploy.sh
```

This catches errors **before** the 5-minute Cloud Run build.

## Mock vs Real MCP

### When to Use Mock Mode

✅ **Use Mock (USE_MOCK_MCP=true) for:**
- Local development
- Testing provider logic
- Testing UI changes
- Fast iteration
- CI/CD pipelines
- Python 3.9 environments

### When to Use Real MCP

✅ **Use Real (USE_MOCK_MCP=false) for:**
- Production deployments
- Integration testing
- Testing actual analysis tools
- Validating MCP server connections

## Environment Variables

### Required

```bash
ANTHROPIC_API_KEY=sk-ant-...     # Claude API key
ENVIRONMENT=development          # development or production
```

### Optional

```bash
GEMINI_API_KEY=...              # Enables Gemini provider
USE_MOCK_MCP=true               # Enable mock mode (default: false)
DEFAULT_MODEL=claude-sonnet-4-5  # Default LLM model
DEFAULT_MAX_TOKENS=4096         # Default max tokens
```

## Troubleshooting

### "Module 'mcp' not found"

**Solution 1:** Enable mock mode (recommended for local dev)
```bash
export USE_MOCK_MCP=true
```

**Solution 2:** Upgrade to Python 3.10+ and install mcp
```bash
python3.10 -m pip install mcp>=1.0.0
```

### "Validation script fails"

Run with verbose output:
```bash
bash -x ./validate.sh
```

Common issues:
- Missing dependencies: `python3 -m pip install -r requirements.txt`
- Wrong Python version: Upgrade to 3.10+
- Missing API keys: Check `.env` file

### "Provider initialization fails"

Check:
1. API keys are set: `echo $ANTHROPIC_API_KEY`
2. Keys are valid (not expired)
3. Internet connection works
4. Run validation: `./validate.sh`

## Testing

### Test Validation Script

```bash
export ANTHROPIC_API_KEY=...
export GEMINI_API_KEY=...
export USE_MOCK_MCP=true
./validate.sh
```

Expected output: `✅ Validation Complete!`

### Test Mock MCP Layer

```python
python3 -c "
import sys
import asyncio
sys.path.insert(0, '.')
from utils.mcp_mock import MockMCPClientManager

async def test():
    servers = [{'name': 'spatialtools', 'url': 'https://mock.run.app'}]
    async with MockMCPClientManager(servers) as manager:
        tools = manager.get_all_tools()
        print(f'✓ Discovered {len(tools)} tools')
        result = await manager.call_tool('spatialtools', 'spatial_autocorrelation',
            {'adata_path': '/test.h5ad', 'gene': 'CD4'})
        print(f'✓ Tool executed: {result[\"content\"][0][\"text\"][:50]}...')

asyncio.run(test())
"
```

### Test Providers

```python
python3 -c "
import sys
import os
sys.path.insert(0, '.')
os.environ['ANTHROPIC_API_KEY'] = 'your-key'
os.environ['GEMINI_API_KEY'] = 'your-key'

from providers.anthropic_provider import AnthropicProvider
from providers.gemini_provider import GeminiProvider

claude = AnthropicProvider()
print(f'✓ Claude: {claude.get_provider_name()} - Available: {claude.is_available()}')

gemini = GeminiProvider()
print(f'✓ Gemini: {gemini.get_provider_name()} - Available: {gemini.is_available()}')
"
```

## File Structure

```
ui/streamlit-app/
├── app.py                          # Main Streamlit app
├── validate.sh                     # Pre-deployment validation
├── deploy.sh                       # Cloud Run deployment
├── .env                           # Local environment (gitignored)
├── .env.example                   # Environment template
├── requirements.txt               # Python dependencies
├── Dockerfile                     # Cloud Run container
├── providers/
│   ├── base.py                    # Provider interface
│   ├── anthropic_provider.py      # Claude integration
│   └── gemini_provider.py         # Gemini integration
└── utils/
    ├── mcp_client.py              # Real MCP client
    ├── mcp_mock.py                # Mock MCP client
    ├── mcp_config.py              # Server configurations
    ├── audit_logger.py            # HIPAA logging
    └── trace_builder.py           # Orchestration traces
```

## Best Practices

1. **Always validate before deploying**
   ```bash
   ./validate.sh && ./deploy.sh
   ```

2. **Use mock mode for local development**
   ```bash
   export USE_MOCK_MCP=true
   ```

3. **Test with both providers**
   - Claude (primary)
   - Gemini (optional)

4. **Keep API keys secure**
   - Never commit `.env` to git
   - Use environment variables in CI/CD
   - Rotate keys regularly

5. **Check validation output**
   - Read warnings (yellow)
   - Fix errors (red) immediately
   - Understand what each check does

## Getting Help

- **Validation fails**: Check error messages, run `./validate.sh` with verbose mode
- **Import errors**: Reinstall dependencies: `python3 -m pip install -r requirements.txt`
- **MCP issues**: Try mock mode first: `export USE_MOCK_MCP=true`
- **Provider errors**: Verify API keys are valid and not expired

## Advanced

### Custom Mock Responses

Edit `utils/mcp_mock.py` to customize mock tool responses for your testing needs.

### CI/CD Integration

```yaml
# Example GitHub Actions workflow
- name: Validate before deploy
  run: |
    export USE_MOCK_MCP=true
    ./validate.sh
- name: Deploy to Cloud Run
  run: ./deploy.sh
```

### Local MCP Server Testing

To test with a local MCP server:

```bash
# Run MCP server locally
cd servers/mcp-spatialtools
python -m mcp_spatialtools --port 8001

# Point to local server in mcp_config.py
# Then set USE_MOCK_MCP=false
```
