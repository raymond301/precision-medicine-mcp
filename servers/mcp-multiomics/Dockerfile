FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Create directory structure
RUN mkdir -p /app/servers/mcp-multiomics /app/shared/utils

# Copy shared utilities
COPY ../../shared/utils/ /app/shared/utils/

# Copy server code
COPY . /app/servers/mcp-multiomics/

# Install Python dependencies
WORKDIR /app/servers/mcp-multiomics
RUN pip install --no-cache-dir -e .

# Set environment variables for SSE transport
ENV MCP_TRANSPORT=sse
ENV MCP_PORT=3001
ENV MULTIOMICS_DRY_RUN=false
ENV MULTIOMICS_DATA_DIR=/app/data/multiomics
ENV MULTIOMICS_CACHE_DIR=/app/data/cache/multiomics

# Create data directories
RUN mkdir -p /app/data/multiomics /app/data/cache/multiomics

# Expose port
EXPOSE 3001

# Health check commented out - FastMCP doesn't provide /health endpoint
# HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
#     CMD curl -f http://localhost:3001/health || exit 1

# Run the server
CMD ["python", "-m", "mcp_multiomics"]
