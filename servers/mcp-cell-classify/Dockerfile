FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install system dependencies (minimal - no TensorFlow/DeepCell needed)
RUN apt-get update && apt-get install -y \
    git \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Copy server code
COPY . /app

# Install Python dependencies
WORKDIR /app
RUN pip install --no-cache-dir -e .

# Set environment variables for SSE transport
ENV MCP_TRANSPORT=sse
ENV MCP_PORT=3009
ENV PORT=3009

# Cell classify configuration
ENV CELL_CLASSIFY_DRY_RUN=false
ENV CELL_CLASSIFY_OUTPUT_DIR=/app/data/output

# Create data directories
RUN mkdir -p /app/data/output

# Expose port
EXPOSE 3009

# Run the server
CMD ["python", "-m", "mcp_cell_classify"]
