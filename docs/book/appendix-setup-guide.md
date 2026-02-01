# Appendix: MCP Server Setup Guide

*Quick reference for setting up your own copy of all MCP servers*

---

## Overview

This appendix provides a quick reference for setting up all 12 MCP servers in your environment. For detailed information, always refer to the linked documentation in the repository.

**Setup options**:
1. **Local (Claude Desktop)** - Run servers locally via stdio transport
2. **Cloud (GCP Cloud Run)** - Deploy servers to production via SSE transport
3. **Development** - Run servers in your IDE for testing and development

---

## Prerequisites

### Software Requirements

| Software | Version | Purpose |
|----------|---------|---------|
| **Python** | 3.11+ | MCP server runtime (3.10 for DeepCell) |
| **Git** | Latest | Clone repository |
| **Docker** | Latest | Cloud deployment (optional for local) |
| **Claude Desktop** | Latest | Local MCP server orchestration |
| **Google Cloud SDK** | Latest | Cloud deployment (Chapters 12-13) |

Installation instructions:
- **Python**: [python.org](https://www.python.org/downloads/)
- **Claude Desktop**: [claude.com/claude-desktop](https://claude.com/claude-desktop)
- **Google Cloud SDK**: [cloud.google.com/sdk](https://cloud.google.com/sdk/docs/install)

### Hardware Requirements

| Resource | Minimum | Recommended |
|----------|---------|-------------|
| RAM | 8GB | 16GB |
| Disk | 50GB free | 100GB free |
| CPU | 4 cores | 8 cores |
| GPU | None | Optional (DeepCell acceleration) |

---

## Quick Start: Local Setup

### Step 1: Clone Repository

```bash
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp
```

### Step 2: Install Server Dependencies

**Install all servers** (from repository root):

```bash
# Install each server in editable mode
cd servers/mcp-spatialtools && pip install -e ".[dev]" && cd ../..
cd servers/mcp-multiomics && pip install -e ".[dev]" && cd ../..
cd servers/mcp-deepcell && pip install -e ".[dev]" && cd ../..
# Repeat for all 12 servers...
```

**Or use the setup script** (if available):
```bash
./scripts/setup-all-servers.sh
```

Full server list and installation instructions: [`servers/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/README.md)

### Step 3: Configure Claude Desktop

**Location**: `~/Library/Application Support/Claude/claude_desktop_config.json` (macOS)

**Minimal configuration** (single server example):

```json
{
  "mcpServers": {
    "spatialtools": {
      "command": "/path/to/spatial-mcp/servers/mcp-spatialtools/venv/bin/python",
      "args": ["-m", "mcp_spatialtools"],
      "env": {
        "SPATIAL_DATA_DIR": "/path/to/data",
        "SPATIAL_DRY_RUN": "false"
      }
    }
  }
}
```

**Full configuration** (all 12 servers):
- Template: [`configs/claude_desktop_config.json`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/configs/claude_desktop_config.json)
- Setup guide: [`docs/test-docs/manual-testing/claude-desktop-setup.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/test-docs/manual-testing/claude-desktop-setup.md)

**Key points**:
- Use **absolute paths** to venv Python executables
- Set `DRY_RUN=false` for real analysis, `true` for demonstration mode
- Restart Claude Desktop after config changes

### Step 4: Verify Installation

**Test in Claude Desktop**:
```
List all available MCP servers and their tools.
```

**Expected response**: 12 servers listed with 124 total tools.

---

## Server-Specific Setup

Each server has unique setup requirements. Refer to individual server READMEs for details:

| Server | Setup Time | Special Requirements | README Link |
|--------|------------|---------------------|-------------|
| **mcp-epic** | 5 min | Epic FHIR credentials (or use mockepic) | [`servers/mcp-epic/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-epic/README.md) |
| **mcp-fgbio** | 10 min | Reference genome files (hg38) | [`servers/mcp-fgbio/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-fgbio/README.md) |
| **mcp-multiomics** | 5 min | None | [`servers/mcp-multiomics/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-multiomics/README.md) |
| **mcp-spatialtools** | 10 min | STAR aligner (optional) | [`servers/mcp-spatialtools/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-spatialtools/README.md) |
| **mcp-deepcell** | 15 min | Python 3.10, TensorFlow 2.8.x | [`servers/mcp-deepcell/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/README.md) |
| **mcp-perturbation** | 5 min | None (mocked) | [`servers/mcp-perturbation/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-perturbation/README.md) |
| **mcp-quantum-celltype-fidelity** | 10 min | PennyLane (for PQC) | [`servers/mcp-quantum-celltype-fidelity/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-quantum-celltype-fidelity/README.md) |
| **mcp-openimagedata** | 5 min | None | [`servers/mcp-openimagedata/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-openimagedata/README.md) |
| **mcp-tcga** | 5 min | None (mocked) | [`servers/mcp-tcga/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-tcga/README.md) |
| **mcp-huggingface** | 5 min | None (mocked) | [`servers/mcp-huggingface/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-huggingface/README.md) |
| **mcp-seqera** | 5 min | Nextflow CLI | [`servers/mcp-seqera/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-seqera/README.md) |
| **mcp-mockepic** | 2 min | None | [`servers/mcp-mockepic/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-mockepic/README.md) |

**Total setup time**: 1-2 hours for all servers (local development)

---

## Cloud Deployment (GCP Cloud Run)

For production deployment to GCP Cloud Run, see **Chapter 12: Cloud Deployment on GCP** (page 167).

### Quick Cloud Deployment

**Deploy a single server**:

```bash
cd servers/mcp-deepcell
./deploy.sh YOUR_PROJECT_ID us-central1
```

**Deploy all servers**:

```bash
./scripts/deploy-all-servers.sh YOUR_PROJECT_ID us-central1
```

Deployment scripts: [`servers/*/deploy.sh`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/servers)

### Cloud Configuration

**Environment variables** (set during deployment):

```bash
# GCP Project
export GCP_PROJECT_ID="your-project-id"
export GCP_REGION="us-central1"

# Optional: Production settings
export DEEPCELL_DRY_RUN="false"
export SPATIAL_DRY_RUN="false"
export MULTIOMICS_DRY_RUN="false"
```

**Using deployed servers** (Claude API with SSE transport):

```python
import anthropic

client = anthropic.Anthropic()

response = client.beta.messages.create(
    model="claude-sonnet-4-5",
    max_tokens=1024,
    messages=[{"role": "user", "content": "Analyze PatientOne"}],
    mcp_servers=[{
        "type": "url",
        "url": "https://mcp-spatialtools-PROJECT_ID.run.app/sse",
        "name": "spatialtools"
    }],
    tools=[{"type": "mcp_toolset", "mcp_server_name": "spatialtools"}],
    betas=["mcp-client-2025-11-20"]
)
```

Full deployment guide: **Chapter 12** (pages 167-179)

---

## PatientOne Dataset Setup

The PatientOne synthetic dataset is required for testing and exercises.

### Option 1: Use GCS Public Bucket (Recommended)

**Location**: `gs://precision-medicine-mcp-public/patient-data/PAT001-OVC-2025/`

**Access**: Public read access (no authentication required)

**Download**:
```bash
# Install gsutil (part of Google Cloud SDK)
gsutil -m cp -r gs://precision-medicine-mcp-public/patient-data/PAT001-OVC-2025/ ./data/
```

### Option 2: Use Local Copy

**Location**: [`data/patient-data/PAT001-OVC-2025/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/data/patient-data/PAT001-OVC-2025)

**Clone with data**:
```bash
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
# PatientOne data included in repository
```

### Data Structure

```
PAT001-OVC-2025/
├── clinical/           # FHIR R4 resources (2 files)
├── genomics/           # Somatic variants VCF (1 file)
├── multiomics/         # RNA/protein/phospho (4 files)
├── spatial/            # 10X Visium (3 files)
└── imaging/            # H&E, MxIF images (7 files)
```

**Total**: 17 files, 100% synthetic (safe to share, no IRB needed)

Full dataset documentation: [`data/patient-data/PAT001-OVC-2025/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/data/patient-data/PAT001-OVC-2025/README.md)

---

## Environment Variables Reference

### Common Variables (All Servers)

| Variable | Default | Description |
|----------|---------|-------------|
| `*_DRY_RUN` | `false` | Enable mock mode (true) or real analysis (false) |
| `*_DATA_DIR` | `/workspace/data` | Data directory path |
| `*_CACHE_DIR` | `/workspace/cache` | Cache directory path |
| `*_LOG_LEVEL` | `INFO` | Logging level (DEBUG, INFO, WARNING, ERROR) |

### Server-Specific Variables

**mcp-spatialtools**:
- `SPATIAL_DATA_DIR`: Spatial transcriptomics data location
- `SPATIAL_DRY_RUN`: Mock mode toggle
- `STAR_PATH`: STAR aligner executable path
- `STAR_GENOME_INDEX`: Reference genome index directory

**mcp-deepcell**:
- `DEEPCELL_DRY_RUN`: Mock mode toggle
- `DEEPCELL_CACHE_DIR`: Model cache directory
- `TF_CPP_MIN_LOG_LEVEL`: TensorFlow logging (default: 3)

**mcp-multiomics**:
- `MULTIOMICS_DATA_DIR`: Multi-omics data location
- `MULTIOMICS_DRY_RUN`: Mock mode toggle

Full environment variable reference: Each server's README.md

---

## Hospital Production Setup

For HIPAA-compliant production deployment in hospital environments, see **Chapter 13: Hospital Production Deployment** (page 180).

### Key Requirements

1. **VPC Networking**: Private deployment with Serverless VPC Connector
2. **Authentication**: Azure AD SSO via OAuth2 Proxy
3. **Compliance**: 10-year immutable audit logs, de-identification
4. **Epic Integration**: FHIR R4 with OAuth2

**Setup scripts**:
- VPC setup: [`scripts/setup-vpc.sh`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/scripts/setup-vpc.sh)
- Audit logs: [`scripts/setup-audit-logs.sh`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/scripts/setup-audit-logs.sh)
- OAuth2 Proxy: [`scripts/setup-oauth2-proxy.sh`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/scripts/setup-oauth2-proxy.sh)

Full hospital setup guide: **Chapter 13** (pages 180-194)

---

## Troubleshooting

### Issue: Claude Desktop doesn't see servers

**Symptoms**: "No MCP servers found" or tools not available

**Solutions**:
1. Verify config file location: `~/Library/Application Support/Claude/claude_desktop_config.json`
2. Check absolute paths to Python executables
3. Restart Claude Desktop after config changes
4. Check server logs: `tail -f ~/.config/Claude/logs/mcp-*.log`

### Issue: Import errors when running servers

**Symptoms**: `ModuleNotFoundError` or `ImportError`

**Solutions**:
1. Activate virtual environment: `source venv/bin/activate`
2. Reinstall dependencies: `pip install -e ".[dev]"`
3. Check Python version: `python --version` (must be 3.11+)
4. Set PYTHONPATH: `export PYTHONPATH=/path/to/server/src`

### Issue: DeepCell requires Python 3.10

**Symptom**: TensorFlow 2.8.x incompatible with Python 3.11+

**Solution**: Use Python 3.10 for mcp-deepcell only
```bash
# Install Python 3.10
brew install python@3.10  # macOS

# Create venv with Python 3.10
cd servers/mcp-deepcell
python3.10 -m venv venv
source venv/bin/activate
pip install -e ".[dev]"
```

Full story: [`servers/mcp-deepcell/DEPENDENCY_ISSUES.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/DEPENDENCY_ISSUES.md)

### Issue: Cloud Run deployment fails

**Symptoms**: Deployment timeout or resource errors

**Solutions**:
1. Check GCP project quotas: `gcloud compute project-info describe --project=PROJECT_ID`
2. Verify Docker is running: `docker ps`
3. Increase timeout: Edit `cloudbuild.yaml` (set `timeout: 1200s`)
4. Check Cloud Build logs: `gcloud builds list --limit=5`

Full deployment troubleshooting: **Chapter 12** (pages 175-178)

### Issue: "Permission denied" in MCP server

**Symptoms**: Cannot read/write files

**Solutions**:
1. Check directory permissions: `ls -la /workspace/data`
2. Create directories: `mkdir -p /workspace/data /workspace/cache`
3. Set ownership: `chown -R $USER /workspace/data`
4. Update config: Ensure `*_DATA_DIR` points to writable location

---

## Testing Your Setup

### Test 1: List Available Servers

**In Claude Desktop**:
```
List all available MCP servers and their tools.
```

**Expected**: 12 servers, 124 tools total

### Test 2: Run Simple Query

**In Claude Desktop**:
```
For PatientOne (patient-001), retrieve clinical demographics using mockepic.
```

**Expected**: Patient demographics returned (Sarah Anderson, 58yo, BRCA1+)

### Test 3: Run Complete Analysis

**In Claude Desktop**:
```
Analyze PatientOne (patient-001) spatial transcriptomics data:
1. Load data using spatialtools
2. Filter low-quality spots (min 200 UMIs, 100 genes)
3. Calculate spatial autocorrelation for MKI67, PCNA, CD8A

Use DRY_RUN=false for real analysis.
```

**Expected**: Real analysis results with Moran's I values

Full test prompts: [`docs/test-docs/manual-testing/TEST_*.txt`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/docs/test-docs/manual-testing)

---

## Development Workflow

### Running Servers in Your IDE

**For development and debugging**:

```bash
# Start server in stdio mode (for testing)
cd servers/mcp-spatialtools
python -m mcp_spatialtools

# Or use MCP Inspector (interactive testing tool)
npx @modelcontextprotocol/inspector python -m mcp_spatialtools
```

**MCP Inspector**: [modelcontextprotocol.io/docs/tools/inspector](https://modelcontextprotocol.io/docs/tools/inspector)

### Running Tests

**Run server tests**:
```bash
cd servers/mcp-spatialtools
pytest
```

**Run with coverage**:
```bash
pytest --cov=src/mcp_spatialtools --cov-report=html
```

Testing guide: [`docs/development/testing.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/development/testing.md)

---

## Additional Resources

### Book Chapters (Setup Referenced)

- **Chapter 4**: Clinical Data (mcp-epic configuration)
- **Chapter 5**: Genomic Foundations (mcp-fgbio setup)
- **Chapter 7**: Spatial Transcriptomics (mcp-spatialtools STAR setup)
- **Chapter 8**: Cell Segmentation (mcp-deepcell Python 3.10)
- **Chapter 12**: Cloud Deployment (GCP Cloud Run)
- **Chapter 13**: Hospital Production (VPC, SSO, HIPAA)
- **Chapter 16**: Teaching (Jupyter notebooks setup)

### Repository Documentation

- **Architecture overview**: [`docs/architecture/README.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/architecture/README.md)
- **Server implementation status**: [`docs/architecture/servers.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/architecture/servers.md)
- **Deployment guides**: [`docs/deployment/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/docs/deployment)
- **Getting started**: [`docs/getting-started/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/docs/getting-started)

### External Resources

- **MCP Specification**: [modelcontextprotocol.io](https://modelcontextprotocol.io)
- **FastMCP Framework**: [github.com/jlowin/fastmcp](https://github.com/jlowin/fastmcp)
- **Claude Desktop**: [claude.com/claude-desktop](https://claude.com/claude-desktop)
- **Claude API Docs**: [docs.anthropic.com](https://docs.anthropic.com)

---

## Summary

**Local setup steps**:
1. Install prerequisites (Python 3.11, Git, Claude Desktop)
2. Clone repository
3. Install server dependencies (`pip install -e ".[dev]"` for each server)
4. Configure Claude Desktop (`claude_desktop_config.json`)
5. Download PatientOne data (GCS or local)
6. Test setup (list servers, run simple query)

**Cloud setup steps** (Chapter 12):
1. Install Google Cloud SDK
2. Authenticate (`gcloud auth login`)
3. Deploy servers (`./deploy.sh PROJECT_ID REGION`)
4. Test with Claude API (SSE transport)

**Production setup steps** (Chapter 13):
1. Run setup scripts (VPC, audit logs, OAuth2 Proxy)
2. Configure Epic FHIR integration
3. Deploy with HIPAA compliance
4. Validate security controls

**Total time**:
- Local: 1-2 hours (all servers)
- Cloud: 2-3 hours (deployment + testing)
- Production: 1 week (security review + compliance)

---

**This appendix consolidates setup information from all chapters. For detailed instructions, always refer to the linked repository documentation.**

---

**Appendix complete**: Quick reference for setting up all 12 MCP servers (local, cloud, production)
**Key sections**: Prerequisites, local setup, cloud deployment, PatientOne data, troubleshooting, testing
**Links**: 25+ links to detailed documentation in repository (chapters, READMEs, scripts, guides)
