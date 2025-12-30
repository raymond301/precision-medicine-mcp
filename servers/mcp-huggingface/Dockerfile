FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Create directory structure
RUN mkdir -p /app/servers/mcp-huggingface /app/shared/utils

# Copy shared utilities
COPY ../../shared/utils/ /app/shared/utils/

# Copy server code
COPY . /app/servers/mcp-huggingface/

# Install Python dependencies
WORKDIR /app/servers/mcp-huggingface
RUN pip install --no-cache-dir -e .

# Set environment variables for SSE transport
ENV MCP_TRANSPORT=sse
ENV MCP_PORT=3006
ENV HF_DRY_RUN=false
ENV HF_DATA_DIR=/app/data/huggingface
ENV HF_CACHE_DIR=/app/data/cache/huggingface

# Create data directories
RUN mkdir -p /app/data/huggingface /app/data/cache/huggingface

# Expose port
EXPOSE 3006

# Health check commented out - FastMCP doesn't provide /health endpoint
# HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
#     CMD curl -f http://localhost:3006/health || exit 1

# Run the server
CMD ["python", "-m", "mcp_huggingface"]
