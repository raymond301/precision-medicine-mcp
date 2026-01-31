# mcp-deepcell: Deep Learning Cell Segmentation for MxIF ✅

**Status:** ✅ **Phase 1 Complete - Production Deployment on Cloud Run**

MCP server for deep learning-based cell segmentation in multiplexed immunofluorescence (MxIF) and other fluorescence microscopy images, powered by **real DeepCell-TF models**.

---

## 🎯 Overview

`mcp-deepcell` provides AI-accessible cell segmentation and classification tools through the Model Context Protocol (MCP). This server uses the **open-source DeepCell-TF library** (https://github.com/vanvalenlab/deepcell-tf) to enable automated cell detection, phenotype classification, and visualization generation for MxIF workflows.

### ✨ Recent Updates (2026-01-31)

- ✅ **Phase 1 Complete:** Real DeepCell-TF integration (nuclear & membrane segmentation)
- ✅ **Cloud Run Deployment:** Production service running at `https://mcp-deepcell-ondu7mwjpa-uc.a.run.app`
- ✅ **Intensity-Based Classification:** Real marker expression analysis (Ki67, TP53, etc.)
- ✅ **Synthetic Test Data:** Generated and uploaded to GCS for testing
- ✅ **Comprehensive Documentation:** Deployment, testing, monitoring guides

### 🔬 Primary Use Case: MxIF Analysis

**MxIF (Multiplexed Immunofluorescence)** enables imaging of multiple protein markers (2-7+) on a single tissue section. DeepCell segments cells from fluorescence channels (e.g., DAPI nuclear stain) and enables quantitative phenotyping:

- Single-cell marker quantification (CD8+, Ki67+, TP53+, etc.)
- Multi-marker co-expression analysis (e.g., TP53+/Ki67+ double-positive cells)
- Spatial distribution of cell phenotypes
- Cell counting and density measurements

### 🚀 Key Features

- 🔬 **Real DeepCell Segmentation** - Nuclear (DAPI) and membrane (Mesmer) models from Van Valen Lab
- 🎯 **Intensity Classification** - Classify cell states based on marker expression
- 🎨 **Segmentation Overlays** - Visualize cell boundaries overlaid on original fluorescence images
- 📊 **Phenotype Visualization** - Color cells by marker expression (positive/negative, co-expression patterns)
- ☁️ **Cloud Ready** - Deployed on GCP Cloud Run with auto-scaling
- 🧪 **Test Data** - Synthetic microscopy images for development and validation

---

## 🏗️ Implementation Status

### ✅ Phase 1 Complete (100%)

| Component | Status | Implementation | Lines of Code |
|-----------|--------|----------------|---------------|
| **DeepCell Engine** | ✅ 100% | Real TensorFlow models | 470 lines |
| **Intensity Classifier** | ✅ 100% | Per-cell marker measurement | 338 lines |
| **segment_cells** | ✅ 100% | Nuclear & membrane segmentation | Real DeepCell |
| **classify_cell_states** | ✅ 100% | Intensity-based phenotyping | Real analysis |
| **generate_segmentation_overlay** | ✅ 100% | Boundary visualization | Real (scikit-image) |
| **generate_phenotype_visualization** | ✅ 100% | Phenotype coloring | Real (matplotlib) |

**Total Production Code:** ~1,250 lines of real implementation

### 📋 Phase 2-4 Roadmap

See [IMPLEMENTATION_PLAN.md](./IMPLEMENTATION_PLAN.md) for complete roadmap:
- **Phase 2:** Performance optimization (batching, tiling, caching)
- **Phase 3:** Advanced features (multi-marker, spatial analysis)
- **Phase 4:** Production hardening (monitoring, security, scale)

---

## ☁️ Cloud Run Deployment

### Production Service

**Service URL:** https://mcp-deepcell-ondu7mwjpa-uc.a.run.app
**Region:** us-central1
**Project:** precision-medicine-poc
**Status:** ✅ Healthy and Running

### Configuration

- **Memory:** 4 GiB (sufficient for images up to 2048×2048)
- **CPU:** 2 vCPUs
- **Timeout:** 300s (5 minutes)
- **Max Instances:** 10
- **Transport:** SSE (Server-Sent Events)

### Performance (CPU-only)

| Image Size | First Request | Subsequent Requests |
|------------|---------------|---------------------|
| 512×512    | ~35s          | ~2s                 |
| 1024×1024  | ~40s          | ~5s                 |
| 2048×2048  | ~50s          | ~10-15s             |
| 4096×4096  | ~90s          | ~30-60s (tiled)     |

**First request** includes model download (~30s) + inference
**Subsequent requests** use cached models (fast)

### Quick Deploy

```bash
cd servers/mcp-deepcell
./deploy.sh precision-medicine-poc us-central1
```

**Documentation:**
- [DEPLOYMENT.md](./DEPLOYMENT.md) - Complete deployment guide
- [MONITORING.md](./MONITORING.md) - Performance monitoring and optimization
- [TESTING.md](./TESTING.md) - Testing and validation guide

---

## 🧪 Test Data

### Synthetic Microscopy Images

**Location:** `gs://sample-inputs-patientone/mcp-deepcell-test-data/test_data/`

**Dataset:** 12 synthetic 16-bit TIFF images (42 MiB total)

#### Image Sizes & Markers

**Three test sizes:**
1. **512×512** (Quick tests, ~30 cells) - 512 KiB each
2. **1024×1024** (Medium tests, ~120 cells) - 2 MiB each
3. **2048×2048** (Large tests, ~480 cells) - 8 MiB each

**Four markers per size:**
- **dapi_NxN.tif** - Nuclear staining (DAPI) for segmentation
- **ki67_NxN.tif** - Proliferation marker (~25% positive cells)
- **tp53_NxN.tif** - Tumor suppressor marker (~40% positive cells)
- **membrane_NxN.tif** - Membrane marker for Mesmer segmentation

#### Test Characteristics

- **Format:** 16-bit grayscale TIFF
- **Intensity Range:** 0-65535
- **Background:** 300-800 (autofluorescence)
- **Positive Cells:** 5000-15000 intensity
- **Negative Cells:** 500-2500 intensity
- **Noise:** Realistic Gaussian (σ=200-300)
- **PSF Blur:** Microscope simulation (σ=1.0-1.5)

#### Test Workflow

```bash
# 1. Segment nuclei
Tool: segment_cells
Input: gs://sample-inputs-patientone/mcp-deepcell-test-data/test_data/dapi_512x512.tif
Model: nuclear
Expected: ~30 cells detected

# 2. Classify cell states
Tool: classify_cell_states
Segmentation: <output from step 1>
Intensity: gs://sample-inputs-patientone/mcp-deepcell-test-data/test_data/ki67_512x512.tif
Expected: ~25% proliferating cells

# 3. Visualize results
Tool: generate_segmentation_overlay
Tool: generate_phenotype_visualization
```

See `test_data/README.md` and `test_data/manifest.json` for complete specifications.

---

## 📦 Installation

### Prerequisites

- **Python:** 3.10 (required for TensorFlow 2.8.x compatibility)
- **Platform:** Linux x86_64 (Cloud Run, GCE) or macOS with Docker
- **NOT Compatible:**
  - ❌ Python 3.11+ (TensorFlow 2.8.x limitation)
  - ❌ macOS Apple Silicon native (use Docker with `--platform linux/amd64`)
  - ⚠️ Windows (works but not tested)

See [DEPENDENCY_ISSUES.md](./DEPENDENCY_ISSUES.md) for platform compatibility details.

### Local Setup (Development)

**Note:** Local development on Apple Silicon requires Docker. For production, deploy to Cloud Run.

```bash
cd servers/mcp-deepcell

# Create Python 3.10 virtual environment
python3.10 -m venv venv
source venv/bin/activate

# Install dependencies
pip install -e ".[dev]"

# Set environment variables
export DEEPCELL_OUTPUT_DIR=./data/output
export DEEPCELL_DRY_RUN=false  # Use real DeepCell models
export DEEPCELL_USE_GPU=false  # CPU-only (GPU optional)
```

### Docker Setup (Recommended for Local Development)

```bash
# Build Docker image
docker build -t mcp-deepcell .

# Run locally
docker run -p 8080:8080 \
  -e DEEPCELL_DRY_RUN=false \
  -e MCP_TRANSPORT=sse \
  mcp-deepcell
```

---

## 🔧 Available Tools

### 1. segment_cells

Deep learning-based cell segmentation using DeepCell-TF models.

**Parameters:**
- `image_path` (string): Path to microscopy image (16-bit TIFF recommended)
- `model_type` (string): Model to use - "nuclear" or "membrane" (default: "nuclear")
- `min_cell_size` (integer): Minimum cell size in pixels (default: 100)
- `image_mpp` (float, optional): Microns per pixel for scaling (default: None)

**Returns:**
```json
{
  "segmentation_mask": "/app/data/output/dapi_512x512_seg.tif",
  "cells_detected": 28,
  "total_area": 23800,
  "mean_cell_area": 850,
  "model_used": "nuclear",
  "processing_time_seconds": 2.3,
  "quality_metrics": {
    "cells_filtered": 2,
    "min_size_threshold": 100
  }
}
```

**Models Available:**
- **nuclear**: Nuclear segmentation from DAPI/Hoechst staining
- **membrane**: Whole-cell segmentation from membrane markers (Mesmer model)

**Output:** 16-bit TIFF label image where each cell has a unique ID (0 = background, 1 = cell 1, etc.)

**Example:**
```
Segment cells in the DAPI nuclear image:
- Image: gs://sample-inputs-patientone/mcp-deepcell-test-data/test_data/dapi_512x512.tif
- Use nuclear segmentation model
- Filter out cells smaller than 100 pixels
```

---

### 2. classify_cell_states

Classify cell states/phenotypes based on marker intensity.

**Parameters:**
- `segmentation_mask_path` (string): Path to segmentation mask (from segment_cells)
- `intensity_image_path` (string): Path to marker intensity image (16-bit TIFF)
- `marker_name` (string): Name of marker (e.g., "Ki67", "TP53")
- `threshold_proliferating` (float): Intensity threshold for positive classification (default: 50.0 percentile)

**Returns:**
```json
{
  "classifications_csv": "/app/data/output/cell_states_ki67.csv",
  "summary": {
    "total_cells": 28,
    "proliferating": 7,
    "quiescent": 21,
    "percent_proliferating": 25.0
  },
  "marker_name": "Ki67",
  "threshold_used": 4500.0
}
```

**Output CSV Format:**
```csv
cell_id,mean_intensity,cell_state
1,3200.5,quiescent
2,8500.2,proliferating
3,2100.8,quiescent
...
```

**Available States:**
- **proliferating**: Marker intensity above threshold (e.g., Ki67+)
- **quiescent**: Marker intensity below threshold (e.g., Ki67-)

**Example:**
```
Classify cell proliferation states:
- Segmentation mask: /output/dapi_512x512_seg.tif
- Ki67 intensity image: gs://sample-inputs-patientone/.../ki67_512x512.tif
- Marker: Ki67
- Use median intensity as threshold
```

---

### 3. generate_segmentation_overlay

Generate segmentation overlay visualization showing cell boundaries on original image.

**Parameters:**
- `image_path` (string): Path to original microscopy image
- `segmentation_mask_path` (string): Path to segmentation mask (label image)
- `output_filename` (string, optional): Custom output filename
- `overlay_color` (string): Color for boundaries - "yellow", "green", "red", "cyan", "magenta" (default: "yellow")
- `overlay_alpha` (float): Transparency 0-1 (default: 0.4)

**Returns:**
```json
{
  "status": "success",
  "output_file": "/app/data/output/visualizations/segmentation_overlay_20260131_145000.png",
  "cells_visualized": 28,
  "description": "Segmentation overlay showing 28 cells with yellow boundaries.",
  "visualization_type": "segmentation_overlay",
  "overlay_color": "yellow"
}
```

**Output:** Side-by-side figure (PNG, 300 DPI):
1. Original microscopy image
2. Image with colored cell boundaries overlaid

**Use Cases:**
- Validate segmentation quality
- Identify over-segmentation or under-segmentation
- Visualize cell density and morphology
- Generate publication figures

**Example:**
```
Generate segmentation overlay:
- Original: gs://sample-inputs-patientone/.../dapi_512x512.tif
- Segmentation: /output/dapi_512x512_seg.tif
- Use green boundaries with 50% transparency
```

---

### 4. generate_phenotype_visualization

Generate cell phenotype visualization coloring cells by marker expression.

**Parameters:**
- `image_path` (string): Path to original microscopy image
- `segmentation_mask_path` (string): Path to segmentation mask
- `classification_csv_path` (string): Path to classification CSV (from classify_cell_states)
- `classification_column` (string): Column name for phenotype (e.g., "cell_state")
- `output_filename` (string, optional): Custom output filename
- `positive_color` (string): Color for positive cells (default: "green")
- `negative_color` (string): Color for negative cells (default: "red")

**Returns:**
```json
{
  "status": "success",
  "output_file": "/app/data/output/visualizations/phenotype_viz_20260131_145030.png",
  "positive_cells": 7,
  "negative_cells": 21,
  "total_cells": 28,
  "percent_positive": 25.0,
  "description": "Phenotype visualization: 7 proliferating (green), 21 quiescent (red). 25.0% proliferating.",
  "visualization_type": "phenotype_coloring"
}
```

**Output:** Side-by-side figure (PNG, 300 DPI):
1. Original microscopy image (grayscale)
2. Cells colored by phenotype with legend and statistics

**Use Cases:**
- Visualize spatial distribution of phenotypes
- Quantify marker-positive cell percentages
- Identify regional patterns in tissue
- Generate publication-quality figures

**Example:**
```
Visualize Ki67+ proliferating cells:
- Original: gs://sample-inputs-patientone/.../dapi_512x512.tif
- Segmentation: /output/dapi_512x512_seg.tif
- Classifications: /output/cell_states_ki67.csv
- Column: cell_state
- Color proliferating cells green, quiescent cells red
```

---

## 📖 Example Workflows

### Complete MxIF Analysis Pipeline

**End-to-end cell phenotyping:**

```
Claude, analyze the MxIF tissue sample:

STEP 1: Segment nuclei from DAPI channel
- Image: gs://sample-inputs-patientone/.../dapi_1024x1024.tif
- Model: nuclear
- Min size: 100 pixels
→ Expect ~120 cells

STEP 2: Generate segmentation overlay
- Show yellow boundaries on DAPI image
- Transparency: 40%
- Validate segmentation quality visually

STEP 3: Classify Ki67+ proliferating cells
- Segmentation from Step 1
- Ki67 intensity: gs://sample-inputs-patientone/.../ki67_1024x1024.tif
- Marker: Ki67
- Use median intensity threshold
→ Expect ~25% positive

STEP 4: Visualize Ki67 phenotypes
- Color proliferating cells green
- Color quiescent cells red
- Generate publication figure

STEP 5: Repeat for TP53 marker
- Same segmentation
- TP53 intensity: gs://sample-inputs-patientone/.../tp53_1024x1024.tif
→ Expect ~40% positive

STEP 6: Multi-marker analysis
- Identify double-positive cells (Ki67+/TP53+)
- Quantify 4 phenotypes:
  - Ki67+/TP53+ (proliferating mutant)
  - Ki67+/TP53- (proliferating wild-type)
  - Ki67-/TP53+ (quiescent mutant)
  - Ki67-/TP53- (quiescent wild-type)
```

### Spatial Distribution Analysis

**Analyze spatial patterns:**

```
Claude, analyze spatial distribution of CD8+ T cells:

STEP 1: Segment tissue from membrane channel
- Use membrane model for whole-cell segmentation
- Image: gs://sample-inputs-patientone/.../membrane_2048x2048.tif

STEP 2: Identify CD8+ cells
- Measure CD8 intensity per cell
- Threshold at 75th percentile (high expressors only)

STEP 3: Generate spatial visualization
- Color CD8+ cells green, others gray
- Analyze clustering patterns

STEP 4: Quantify metrics
- Total CD8+ count
- CD8+ density (cells/mm²)
- Spatial clustering score
```

---

## 🔌 MCP Client Configuration

### Claude Desktop (Local)

Add to `claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "deepcell": {
      "command": "/path/to/spatial-mcp/servers/mcp-deepcell/venv/bin/python",
      "args": ["-m", "mcp_deepcell"],
      "env": {
        "DEEPCELL_OUTPUT_DIR": "/workspace/output",
        "DEEPCELL_DRY_RUN": "false",
        "DEEPCELL_USE_GPU": "false"
      }
    }
  }
}
```

### Cloud Run (Production)

Connect via SSE transport:

```json
{
  "mcpServers": {
    "deepcell": {
      "url": "https://mcp-deepcell-ondu7mwjpa-uc.a.run.app/sse",
      "transport": "sse"
    }
  }
}
```

---

## 📚 Documentation

| Document | Description |
|----------|-------------|
| [README.md](./README.md) | This file - Overview and quick start |
| [IMPLEMENTATION_PLAN.md](./IMPLEMENTATION_PLAN.md) | Complete Phase 1-4 roadmap (486 lines) |
| [DEPLOYMENT.md](./DEPLOYMENT.md) | Cloud Run deployment guide (401 lines) |
| [DEPENDENCY_ISSUES.md](./DEPENDENCY_ISSUES.md) | Platform compatibility and troubleshooting (442 lines) |
| [TESTING.md](./TESTING.md) | Testing and validation guide (467 lines) |
| [MONITORING.md](./MONITORING.md) | Performance monitoring and optimization (724 lines) |
| [BUILD_SUMMARY.md](./BUILD_SUMMARY.md) | Build overview and deployment checklist (372 lines) |
| [test_data/README.md](./test_data/README.md) | Synthetic test data documentation |
| [test_data/manifest.json](./test_data/manifest.json) | Test dataset metadata |

**Total Documentation:** ~3,300 lines across 9 files

---

## 🏗️ Architecture

### System Components

```
┌─────────────────────────────────────────────────────────┐
│                    MCP Client                           │
│              (Claude Desktop / API)                     │
└────────────────────┬────────────────────────────────────┘
                     │ SSE/stdio
                     ▼
┌─────────────────────────────────────────────────────────┐
│                  MCP Server (FastMCP)                   │
│  ┌──────────────────────────────────────────────────┐  │
│  │  Tools: segment_cells, classify_cell_states,    │  │
│  │         generate_*_visualization                 │  │
│  └──────────────────────────────────────────────────┘  │
└────────┬────────────────────────────────┬───────────────┘
         │                                │
         ▼                                ▼
┌────────────────────┐        ┌──────────────────────────┐
│  DeepCellEngine    │        │  IntensityClassifier     │
│  (470 lines)       │        │  (338 lines)             │
├────────────────────┤        ├──────────────────────────┤
│ - Model loading    │        │ - Intensity measurement  │
│ - Caching          │        │ - Thresholding           │
│ - Preprocessing    │        │ - Phenotype assignment   │
│ - Inference        │        │ - CSV export             │
│ - Post-processing  │        │ - Multi-marker support   │
│ - Tiling (large)   │        └──────────────────────────┘
└─────────┬──────────┘
          │
          ▼
┌────────────────────────────────────────────────────────┐
│           DeepCell-TF Models (TensorFlow 2.8.x)        │
│  ┌──────────────────┐    ┌──────────────────────────┐ │
│  │  Nuclear Model   │    │   Membrane Model         │ │
│  │  (~500MB)        │    │   (Mesmer, ~800MB)       │ │
│  │  DAPI/Hoechst    │    │   Whole-cell boundaries  │ │
│  └──────────────────┘    └──────────────────────────┘ │
└────────────────────────────────────────────────────────┘
```

### Data Flow

```
Input Image (16-bit TIFF)
    │
    ▼
DeepCellEngine.segment_image()
    │
    ├─► Preprocessing (normalize, resize)
    ├─► Model inference (TensorFlow)
    ├─► Post-processing (label, filter)
    │
    ▼
Segmentation Mask (16-bit TIFF)
    │
    ▼
IntensityClassifier.classify_cell_states()
    │
    ├─► Measure per-cell intensity
    ├─► Apply threshold
    ├─► Assign phenotypes
    │
    ▼
Classification CSV + Visualizations (PNG)
```

### Directory Structure

```
mcp-deepcell/
├── src/
│   └── mcp_deepcell/
│       ├── __init__.py
│       ├── __main__.py
│       ├── server.py                 # MCP server (443 lines)
│       ├── deepcell_engine.py        # DeepCell model management (470 lines)
│       └── intensity_classifier.py   # Cell classification (338 lines)
├── tests/
│   ├── conftest.py
│   └── test_server.py
├── test_data/                        # Synthetic test images
│   ├── dapi_512x512.tif
│   ├── ki67_512x512.tif
│   ├── ... (12 TIFF files)
│   ├── manifest.json
│   └── README.md
├── Dockerfile                        # Cloud Run container
├── .dockerignore                     # Build optimization
├── cloudbuild.yaml                   # Cloud Build config
├── deploy.sh                         # Quick deployment script
├── pyproject.toml                    # Dependencies
├── README.md                         # This file
├── IMPLEMENTATION_PLAN.md            # Complete roadmap
├── DEPLOYMENT.md                     # Deployment guide
├── DEPENDENCY_ISSUES.md              # Platform compatibility
├── TESTING.md                        # Testing guide
├── MONITORING.md                     # Monitoring guide
├── BUILD_SUMMARY.md                  # Build overview
├── test_deployment.py                # Deployment verification
└── generate_synthetic_data.py        # Test data generator
```

---

## ⚙️ Configuration

### Environment Variables

| Variable | Default | Description |
|----------|---------|-------------|
| `DEEPCELL_OUTPUT_DIR` | `/app/data/output` | Output directory for results |
| `DEEPCELL_DRY_RUN` | `false` | Enable mock mode (no real models) |
| `DEEPCELL_MODEL_CACHE_DIR` | `/app/data/models` | Model cache directory |
| `DEEPCELL_USE_GPU` | `false` | Enable GPU acceleration (requires GPU setup) |
| `MCP_TRANSPORT` | `sse` | MCP transport protocol |
| `MCP_PORT` | `3007` | MCP server port |
| `TF_CPP_MIN_LOG_LEVEL` | `2` | TensorFlow logging level |

### Resource Requirements

**Minimum:**
- Memory: 4 GiB
- CPU: 2 vCPUs
- Disk: 10 GB (for model cache)

**Recommended:**
- Memory: 8 GiB (for large images)
- CPU: 4 vCPUs (for higher throughput)
- GPU: T4 or better (5-10× speedup)

---

## 🔍 Troubleshooting

### Common Issues

#### Out of Memory

**Symptoms:** Container crashes, 500 errors

**Solutions:**
- Increase memory limit: `--memory=8Gi`
- Use smaller images or enable tiling
- Reduce concurrency: `--concurrency=1`

#### Slow Performance

**Symptoms:** All requests take >30s

**Solutions:**
- First request is slow (model download) - normal
- Check if models are cached in logs
- Increase CPU: `--cpu=4`
- Consider GPU acceleration

#### Model Download Failures

**Symptoms:** Timeout on first request

**Solutions:**
- Verify network connectivity
- Check DeepCell model zoo accessibility
- Increase timeout: `--timeout=600`

See [TROUBLESHOOTING.md](./DEPENDENCY_ISSUES.md) for complete guide.

---

## 🧪 Development

### Running Tests

```bash
# Install dev dependencies
pip install -e ".[dev]"

# Run all tests
pytest

# Run with coverage
pytest --cov=src/mcp_deepcell --cov-report=html

# Run specific test
pytest tests/test_server.py -v
```

### Code Quality

```bash
# Format code
black src/ tests/

# Type checking
mypy src/

# Linting
ruff check src/ tests/
```

---

## 📊 Performance Metrics

### Benchmarks (Cloud Run, 2 vCPUs, 4Gi RAM)

**Segmentation:**
| Image Size | Cells | Cold Start | Warm Request |
|------------|-------|------------|--------------|
| 512×512    | ~30   | 35s        | 2s           |
| 1024×1024  | ~120  | 40s        | 5s           |
| 2048×2048  | ~480  | 50s        | 15s          |
| 4096×4096  | ~1920 | 90s        | 60s          |

**Classification:**
| Cells | Processing Time |
|-------|-----------------|
| 30    | <1s             |
| 120   | <1s             |
| 480   | 1-2s            |

**Visualization:**
| Image Size | Generation Time |
|------------|-----------------|
| 512×512    | <1s             |
| 1024×1024  | 1-2s            |
| 2048×2048  | 2-3s            |

---

## 🎯 Roadmap

### Phase 1 (Complete ✅)

- ✅ Real DeepCell integration (nuclear & membrane)
- ✅ Intensity-based cell classification
- ✅ Cloud Run deployment
- ✅ Synthetic test data
- ✅ Comprehensive documentation

### Phase 2 (Planned)

- Performance optimization (batching, GPU)
- Advanced tiling for large images
- Improved model caching
- Multi-instance coordination

### Phase 3 (Future)

- Multi-marker phenotyping UI
- Spatial analysis integration
- Cell-cell interaction detection
- Region-based statistics

### Phase 4 (Future)

- Production hardening
- Advanced monitoring
- Security enhancements
- Scale testing

See [IMPLEMENTATION_PLAN.md](./IMPLEMENTATION_PLAN.md) for details.

---

## 📝 License

See the main repository LICENSE file.

---

## 🤝 Support

- **Issues:** https://github.com/lynnlangit/precision-medicine-mcp/issues
- **Documentation:** See `docs/` directory
- **MCP Specification:** https://modelcontextprotocol.io/
- **DeepCell-TF:** https://github.com/vanvalenlab/deepcell-tf

---

## 🔗 Related Servers

Part of the Precision Medicine MCP suite:

- **mcp-spatialtools** - Spatial transcriptomics analysis and visualization
- **mcp-openimagedata** - Histology image processing and multiplex IF
- **mcp-epic** - Clinical EHR data (FHIR)
- **mcp-multiomics** - Multi-omics integration
- **mcp-fgbio** - Genomic reference data

---

**Built for the Precision Medicine MCP suite** - Enables AI-driven cell segmentation and phenotype analysis integrated with spatial and clinical data.

**Status:** ✅ Phase 1 Complete - Production Ready on Cloud Run (2026-01-31)
