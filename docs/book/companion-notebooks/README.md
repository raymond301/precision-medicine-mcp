# Companion Jupyter Notebooks

**Hands-on exercises for _Building AI-Orchestrated Precision Oncology Systems_**

---

## Overview

Each chapter in the book has a corresponding Jupyter notebook with:
- Executable code examples
- Interactive exercises
- "Try changing this parameter..." experiments
- Links to deployed MCP servers (once you deploy them)

**IMPORTANT**: These notebooks require you to deploy your own MCP servers to GCP Cloud Run. They will NOT work without your own infrastructure.

---

## Prerequisites

### Required

1. **Python 3.11+** installed
2. **Git** for cloning the repository
3. **GCP Project** with Cloud Run enabled
4. **MCP Servers deployed** to your GCP project (see Setup below)
5. **API Key** - Either Anthropic API key (for Claude) or Google AI API key (for Gemini)

### Optional (Recommended)

- PatientOne dataset (100% synthetic, included in repository)
- JupyterLab or VS Code with Jupyter extension
- Google Cloud SDK (for deployment)

---

## Setup

### Step 1: Clone Repository

```bash
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp/docs/book/companion-notebooks
```

### Step 2: Install Dependencies

```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install required packages
pip install -r requirements.txt
```

### Step 3: Deploy MCP Servers

**CRITICAL**: You must deploy the MCP servers to your own GCP project.

**Quick deployment**:
```bash
# Deploy all servers
./infrastructure/deployment/deploy_to_gcp.sh YOUR_PROJECT_ID us-central1
```

**Detailed setup**: See **Appendix: MCP Server Setup Guide** in the book

### Step 4: Configure API Keys

Create a `.env` file in this directory:

```bash
# Anthropic (Claude)
ANTHROPIC_API_KEY=your_anthropic_api_key_here

# OR Google AI (Gemini)
GOOGLE_API_KEY=your_google_api_key_here

# Your MCP server URLs (replace YOUR_PROJECT_ID)
MCP_FGBIO_URL=https://mcp-fgbio-YOUR_PROJECT_ID.run.app/sse
MCP_MULTIOMICS_URL=https://mcp-multiomics-YOUR_PROJECT_ID.run.app/sse
MCP_SPATIALTOOLS_URL=https://mcp-spatialtools-YOUR_PROJECT_ID.run.app/sse
```

### Step 5: Launch Jupyter

```bash
jupyter lab
```

---

## Available Notebooks

### Part 1: Why This Matters
- `chapter-01-patientone-story.ipynb` - End-to-end workflow demo
- `chapter-02-architecture.ipynb` - MCP server discovery and orchestration
- `chapter-03-testing.ipynb` - Test coverage and cost analysis

### Part 2: Building the Foundation
- `chapter-04-clinical-data.ipynb` - FHIR integration and de-identification
- `chapter-05-genomics.ipynb` - VCF parsing and variant annotation
- `chapter-06-multiomics.ipynb` - HAllA and Stouffer meta-analysis
- `chapter-07-spatial.ipynb` - Spatial transcriptomics analysis

### Part 3: Advanced Capabilities
- `chapter-08-deepcell.ipynb` - Cell segmentation with DeepCell
- `chapter-09-treatment-response.ipynb` - GEARS GNN prediction
- `chapter-10-quantum.ipynb` - Quantum fidelity and Bayesian UQ
- `chapter-11-imaging.ipynb` - H&E and MxIF analysis

### Part 4: Deployment and Operations
- `chapter-12-cloud-deployment.ipynb` - Docker and Cloud Run
- `chapter-13-hospital-deployment.ipynb` - HIPAA compliance and VPC
- `chapter-14-operations.ipynb` - Logging, alerts, runbooks

### Part 5: Research and Education
- `chapter-15-research.ipynb` - Exploratory analysis and prompts
- `chapter-16-education.ipynb` - Educational workflows
- `chapter-17-funding.ipynb` - ROI calculator and grant budgets
- `chapter-18-lessons-learned.ipynb` - Production insights

**Total**: 18 Jupyter notebooks

---

## Cost Expectations

### Development/Learning
- **Local execution**: Free
- **Claude API**: ~$5-10 for all 18 notebooks
- **Gemini API**: ~$2-5 for all 18 notebooks

### Your MCP Server Deployment
- **Cloud Run**: ~$0.02-0.21 per analysis
- **100 analyses**: ~$2-20 total
- **GCP Free Tier**: $300 credit for 90 days

**Total cost to complete all exercises**: ~$10-20

---

## Support

**Issues**: https://github.com/lynnlangit/precision-medicine-mcp/issues
**Documentation**: https://github.com/lynnlangit/precision-medicine-mcp
**Book**: [`../README.md`](../README.md)
**Appendix**: [`../appendix-b-installation-setup.md`](../appendix-b-installation-setup.md)

---

**Note**: These notebooks are companions to the book. Deploy your own MCP servers (see Appendix) before running notebooks.

**Last Updated**: 2026-02-01
**Status**: Template structure complete, individual notebooks in development
