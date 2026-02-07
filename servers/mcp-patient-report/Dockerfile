FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install system dependencies for WeasyPrint PDF generation
RUN apt-get update && apt-get install -y \
    git \
    curl \
    libpango-1.0-0 \
    libpangocairo-1.0-0 \
    libgdk-pixbuf-2.0-0 \
    libffi-dev \
    && rm -rf /var/lib/apt/lists/*

# Create directory structure
RUN mkdir -p /app/servers/mcp-patient-report /app/shared/utils /app/templates

# Copy temporarily staged shared utilities
COPY _shared_temp/utils/ /app/shared/utils/

# Copy templates from repo root
COPY _templates_temp/ /app/templates/

# Copy server code
COPY . /app/servers/mcp-patient-report/

# Install Python dependencies
WORKDIR /app/servers/mcp-patient-report
RUN pip install --no-cache-dir -e .

# Set environment variables for SSE transport
ENV MCP_TRANSPORT=sse
ENV MCP_PORT=3011
ENV PATIENT_REPORT_DRY_RUN=false
ENV PATIENT_REPORT_OUTPUT_DIR=/app/data/reports
ENV PATIENT_REPORT_TEMPLATES_DIR=/app/templates/patient_report

# Create data directories
RUN mkdir -p /app/data/reports

# Expose port
EXPOSE 3011

# Run the server
CMD ["python", "-m", "mcp_patient_report"]
