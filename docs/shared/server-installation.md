# Standard MCP Server Installation (Canonical Reference)

Standard setup for Python-based MCP servers in the Precision Medicine platform.

## Prerequisites

- **Python 3.11+** (except mcp-deepcell which requires Python 3.10 for TensorFlow 2.8.x)
- **pip** package manager
- **Git** (for cloning repository)

## Virtual Environment Setup

```bash
# Navigate to the server directory
cd servers/mcp-{server-name}

# Create virtual environment
python -m venv venv

# Activate
source venv/bin/activate  # macOS/Linux
# venv\Scripts\activate   # Windows

# Install with development dependencies
pip install -e ".[dev]"
```

## Claude Desktop Configuration

Add to `~/Library/Application Support/Claude/claude_desktop_config.json` (macOS):

```json
{
  "mcpServers": {
    "{server-name}": {
      "command": "/absolute/path/to/servers/mcp-{server-name}/venv/bin/python",
      "args": ["-m", "mcp_{server_name}"],
      "cwd": "/absolute/path/to/servers/mcp-{server-name}",
      "env": {
        "PYTHONPATH": "/absolute/path/to/servers/mcp-{server-name}/src",
        "{SERVER_NAME}_DRY_RUN": "true"
      }
    }
  }
}
```

**Important:** Use absolute paths to the venv Python executable, not just `python`. Claude Desktop requires absolute paths.

## Environment Variables

Each server follows this naming convention:

| Variable | Pattern | Description |
|----------|---------|-------------|
| DRY_RUN | `{SERVER_NAME}_DRY_RUN` | Enable mock mode (`true`/`false`) |
| Data directory | `{SERVER_NAME}_DATA_DIR` | Primary data directory |
| Cache directory | `{SERVER_NAME}_CACHE_DIR` | Cache/model storage |
| Output directory | `{SERVER_NAME}_OUTPUT_DIR` | Analysis output |

See each server's README for server-specific environment variables.

## Running the Server

```bash
# Standalone (stdio)
python -m mcp_{server_name}

# Verify it starts
python -m mcp_{server_name} --help  # If supported
```

## Development Commands

```bash
# Run tests
pytest

# Run tests with coverage
pytest --cov=src --cov-report=html

# Code formatting
black src/ tests/

# Linting
ruff check src/ tests/
```

---

**DRY_RUN mode details:** [dry-run-mode.md](dry-run-mode.md)
**Server status matrix:** [platform-overview.md](platform-overview.md)
