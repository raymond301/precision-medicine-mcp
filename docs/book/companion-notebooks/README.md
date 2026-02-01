# Companion Jupyter Notebooks

**For the book**: *Building AI-Orchestrated Precision Oncology Systems*

These Jupyter notebooks provide hands-on, executable examples for each chapter of the book. You can run these notebooks locally or in Google Colab to interact with the deployed MCP servers on GCP Cloud Run.

---

## Setup

### Local Installation

```bash
# Clone the repository
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp/docs/book/companion-notebooks

# Create virtual environment
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Launch Jupyter
jupyter notebook
```

### Google Colab

Click the "Open in Colab" badge at the top of each notebook to run in the cloud without local setup.

---

## Authentication

Most notebooks require API keys for Claude or Gemini:

**Claude (Anthropic)**:
```python
import os
os.environ["ANTHROPIC_API_KEY"] = "your_key_here"
```

**Gemini (Google AI)**:
```python
import os
os.environ["GOOGLE_API_KEY"] = "your_key_here"
```

Get your keys:
- Claude: https://console.anthropic.com/
- Gemini: https://aistudio.google.com/app/apikey

---

## Notebook Structure

Each chapter has a corresponding notebook:

| Chapter | Notebook | Description |
|---------|----------|-------------|
| 1 | `chapter-01-patientone-story.ipynb` | PatientOne scenario walkthrough |
| 2 | `chapter-02-architecture.ipynb` | MCP architecture exploration |
| 3 | `chapter-03-testing.ipynb` | Testing servers and tools |
| 4 | `chapter-04-clinical-data.ipynb` | FHIR data with mcp-epic |
| 5 | `chapter-05-genomics.ipynb` | Genomic QC with mcp-fgbio |
| 6 | `chapter-06-multiomics.ipynb` | Multi-omics integration |
| 7 | `chapter-07-spatial.ipynb` | Spatial transcriptomics analysis |
| 8 | `chapter-08-deepcell.ipynb` | Cell segmentation |
| 9 | `chapter-09-perturbation.ipynb` | Treatment response prediction |
| 10 | `chapter-10-quantum.ipynb` | Quantum fidelity analysis |
| 11 | `chapter-11-imaging.ipynb` | Histopathology imaging |
| 12 | `chapter-12-cloud-deployment.ipynb` | Deploy to Cloud Run |
| 13 | `chapter-13-hospital-deployment.ipynb` | HIPAA-compliant setup |
| 14 | `chapter-14-operations.ipynb` | Monitoring and operations |
| 15 | `chapter-15-research.ipynb` | Research workflows |
| 16 | `chapter-16-education.ipynb` | Educational use cases |
| 17 | `chapter-17-funding.ipynb` | ROI analysis |
| 18 | `chapter-18-lessons.ipynb` | Lessons learned exploration |

---

## Sample Data

Notebooks use the PatientOne dataset stored in GCS:

```
gs://sample-inputs-patientone/
├── fhir/           # Clinical data
├── genomics/       # VCF files
├── multiomics/     # RNA, protein, phospho
├── spatial/        # 10X Visium data
└── imaging/        # H&E and MxIF images
```

Public read access is enabled for educational use.

---

## MCP Server Endpoints

Deployed Cloud Run services (as of 2026-01-31):

- **mcp-epic**: Local only (HIPAA compliance)
- **mcp-fgbio**: https://mcp-fgbio-ondu7mwjpa-uc.a.run.app
- **mcp-multiomics**: https://mcp-multiomics-ondu7mwjpa-uc.a.run.app
- **mcp-spatialtools**: https://mcp-spatialtools-ondu7mwjpa-uc.a.run.app
- **mcp-perturbation**: https://mcp-perturbation-ondu7mwjpa-uc.a.run.app
- **mcp-quantum-celltype-fidelity**: https://mcp-quantum-celltype-fidelity-ondu7mwjpa-uc.a.run.app
- **mcp-deepcell**: https://mcp-deepcell-ondu7mwjpa-uc.a.run.app
- **mcp-openimagedata**: https://mcp-openimagedata-ondu7mwjpa-uc.a.run.app
- **mcp-tcga**: https://mcp-tcga-ondu7mwjpa-uc.a.run.app
- **mcp-huggingface**: https://mcp-huggingface-ondu7mwjpa-uc.a.run.app
- **mcp-seqera**: https://mcp-seqera-ondu7mwjpa-uc.a.run.app

All servers use SSE transport: `{base_url}/sse`

---

## Troubleshooting

**"API key not found"**: Set environment variables before running cells

**"Server connection failed"**: Check Cloud Run service status:
```bash
gcloud run services list --region us-central1
```

**"Out of memory"**: Some notebooks require 4GB+ RAM. Use Colab Pro or run locally.

**"Rate limit exceeded"**: Add delays between API calls or reduce iteration count.

---

## Contributing

Found a bug or have an improvement? Open an issue or PR:
- Issues: https://github.com/lynnlangit/precision-medicine-mcp/issues
- Pull Requests: https://github.com/lynnlangit/precision-medicine-mcp/pulls

---

**License**: Apache 2.0
**Maintained by**: Precision Medicine MCP Team
**Last Updated**: 2026-01-31
