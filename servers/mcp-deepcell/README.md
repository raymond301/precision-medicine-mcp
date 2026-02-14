# mcp-deepcell: Deep Learning Cell Segmentation + Quantification for MxIF

**Status:** Production Deployment on Cloud Run

MCP server for deep learning-based cell segmentation and per-cell marker quantification in multiplexed immunofluorescence (MxIF) and other fluorescence microscopy images, powered by **real DeepCell-TF models**.

> **Note:** Classification tools (classify_cell_states, generate_phenotype_visualization) have moved to the lightweight **mcp-cell-classify** server. See [../mcp-cell-classify/README.md](../mcp-cell-classify/README.md).

---

## Overview

`mcp-deepcell` provides AI-accessible cell segmentation and marker quantification tools through the Model Context Protocol (MCP). This server uses the **open-source DeepCell-TF library** (https://github.com/vanvalenlab/deepcell-tf) to enable automated cell detection, per-cell intensity measurement, and visualization generation for MxIF workflows.

### Primary Use Case: MxIF Analysis

**MxIF (Multiplexed Immunofluorescence)** enables imaging of multiple protein markers (2-7+) on a single tissue section. DeepCell segments cells from fluorescence channels (e.g., DAPI nuclear stain) and `quantify_markers` measures per-cell intensities:

- Single-cell marker quantification (CD8, Ki67, TP53, etc.)
- Per-cell intensity CSV for downstream classification
- Cell counting and density measurements
- Segmentation quality visualization

### Key Features

- **Real DeepCell Segmentation** - Nuclear (DAPI) and membrane (Mesmer) models from Van Valen Lab
- **Marker Quantification** - Per-cell mean/max/min intensity for any number of markers
- **Segmentation Overlays** - Visualize cell boundaries overlaid on original fluorescence images
- **Cloud Ready** - Deployed on GCP Cloud Run with auto-scaling

---

## Available Tools

### 1. segment_cells

Deep learning-based cell segmentation using DeepCell-TF models.

**Parameters:**
- `image_path` (string): Path to microscopy image (16-bit TIFF recommended)
- `model_type` (string): Model to use - "nuclear" or "membrane" (default: "nuclear")
- `min_cell_size` (integer): Minimum cell size in pixels (default: 100)
- `image_mpp` (float, optional): Microns per pixel for scaling

**Returns:**
```json
{
  "segmentation_mask": "/output/segmentations/dapi_512x512_segmentation_nuclear.tif",
  "cells_detected": 28,
  "mean_cell_area": 850,
  "model_used": "nuclear",
  "processing_time_seconds": 2.3
}
```

**Models Available:**
- **nuclear**: Nuclear segmentation from DAPI/Hoechst staining
- **membrane**: Whole-cell segmentation from membrane markers (Mesmer model)

**Output:** 16-bit TIFF label image where each cell has a unique ID (0 = background, 1 = cell 1, etc.)

---

### 2. quantify_markers

Measure per-cell marker intensities from a segmentation mask and one or more marker images.

**Parameters:**
- `segmentation_mask_path` (string): Path to segmentation mask (16-bit TIFF from segment_cells)
- `marker_images` (list): List of `{"path": "...", "name": "..."}` dicts for each marker
- `output_filename` (string, optional): Custom output CSV filename

**Returns:**
```json
{
  "quantification_csv": "/output/quantifications/quantification_Ki67_TP53_seg.csv",
  "total_cells": 28,
  "markers_quantified": ["Ki67", "TP53"],
  "marker_summaries": {
    "Ki67": {"mean": 3200.5, "std": 2100.3, "min": 120.0, "max": 12500.0},
    "TP53": {"mean": 4100.2, "std": 1800.1, "min": 200.0, "max": 15000.0}
  }
}
```

**Output CSV Format:**
```csv
cell_id,Ki67_mean_intensity,Ki67_max_intensity,Ki67_min_intensity,TP53_mean_intensity,TP53_max_intensity,TP53_min_intensity
1,3200.5,8500.0,120.0,4100.2,12000.0,200.0
2,8500.2,15000.0,500.0,2100.8,6000.0,150.0
```

**Use Cases:**
- Feed into mcp-cell-classify for phenotype classification
- Export for analysis in FlowSOM, Leiden, scikit-learn, or other tools
- Per-cell quantification for downstream spatial analysis

---

### 3. generate_segmentation_overlay

Generate segmentation overlay visualization showing cell boundaries on original image.

**Parameters:**
- `original_image_path` (string): Path to original microscopy image
- `segmentation_mask_path` (string): Path to segmentation mask
- `output_filename` (string, optional): Custom output filename
- `overlay_color` (string): Color for boundaries (default: "yellow")
- `overlay_alpha` (float): Transparency 0-1 (default: 0.4)

**Returns:**
```json
{
  "status": "success",
  "output_file": "/output/visualizations/segmentation_overlay.png",
  "cells_visualized": 28,
  "visualization_type": "segmentation_overlay"
}
```

---

## Workflow

### Segmentation + Quantification (mcp-deepcell only)

```
Step 1: Segment nuclei
  Tool: segment_cells
  Input: dapi_1024x1024.tif, model=nuclear
  Output: segmentation_mask.tif

Step 2: Quantify markers
  Tool: quantify_markers
  Input: segmentation_mask.tif + [ki67.tif, tp53.tif]
  Output: per_cell_intensities.csv

Step 3: Validate segmentation
  Tool: generate_segmentation_overlay
  Input: original_image + segmentation_mask
  Output: overlay_visualization.png
```

### Full Pipeline (with mcp-cell-classify)

```
mcp-deepcell                          mcp-cell-classify
┌──────────────────────┐              ┌──────────────────────────┐
│  segment_cells       │──mask.tif──►│  classify_cell_states    │
│  quantify_markers    │──csv──────►  │  classify_multi_marker   │
│  generate_seg_overlay│              │  generate_phenotype_viz  │
└──────────────────────┘              └──────────────────────────┘
```

---

## Cloud Run Deployment

### Production Service

**Service URL:** https://mcp-deepcell-ondu7mwjpa-uc.a.run.app
**Region:** us-central1
**Project:** precision-medicine-poc

### Configuration

- **Memory:** 4 GiB (sufficient for images up to 2048x2048)
- **CPU:** 2 vCPUs
- **Timeout:** 300s (5 minutes)

### Quick Deploy

```bash
cd servers/mcp-deepcell
./deploy.sh precision-medicine-poc us-central1
```

---

## Test Data

### Synthetic Microscopy Images

**Location:** `gs://sample-inputs-patientone/mcp-deepcell-test-data/test_data/`

**Dataset:** 12 synthetic 16-bit TIFF images (42 MiB total)

**Three test sizes:** 512x512, 1024x1024, 2048x2048
**Four markers per size:** dapi, ki67, tp53, membrane

### Test Workflow

```bash
# 1. Segment nuclei
Tool: segment_cells
Input: gs://sample-inputs-patientone/mcp-deepcell-test-data/test_data/dapi_512x512.tif
Model: nuclear

# 2. Quantify markers
Tool: quantify_markers
Segmentation: <output from step 1>
Markers: [ki67_512x512.tif, tp53_512x512.tif]

# 3. Visualize segmentation
Tool: generate_segmentation_overlay

# 4. Classify (use mcp-cell-classify)
Tool: classify_multi_marker (on mcp-cell-classify server)
```

---

## Installation

> **Standard setup:** See [Server Installation Guide](../../docs/shared/server-installation.md) for common setup steps. Note: this server requires Python 3.10 (not 3.11+) due to TensorFlow compatibility.

### Prerequisites

- **Python:** 3.10 (required for TensorFlow 2.8.x compatibility)
- **Platform:** Linux x86_64 (Cloud Run, GCE) or macOS with Docker

### Local Setup (Development)

```bash
cd servers/mcp-deepcell
python3.10 -m venv venv
source venv/bin/activate
pip install -e ".[dev]"

export DEEPCELL_OUTPUT_DIR=./data/output
export DEEPCELL_DRY_RUN=false
export DEEPCELL_USE_GPU=false
```

### Docker Setup

```bash
docker build -t mcp-deepcell .
docker run -p 8080:8080 \
  -e DEEPCELL_DRY_RUN=false \
  -e MCP_TRANSPORT=sse \
  mcp-deepcell
```

---

## MCP Client Configuration

### Claude Desktop (Local)

```json
{
  "mcpServers": {
    "deepcell": {
      "command": "uv",
      "args": ["run", "--directory", "/path/to/servers/mcp-deepcell", "python", "-m", "mcp_deepcell"],
      "env": {
        "DEEPCELL_OUTPUT_DIR": "/workspace/output",
        "DEEPCELL_DRY_RUN": "true",
        "DEEPCELL_USE_GPU": "false"
      }
    }
  }
}
```

---

## Configuration

| Variable | Default | Description |
|----------|---------|-------------|
| `DEEPCELL_OUTPUT_DIR` | `/workspace/output` | Output directory for results |
| `DEEPCELL_DRY_RUN` | `true` | Enable mock mode (no real models) |
| `DEEPCELL_MODEL_CACHE_DIR` | `~/.deepcell/models` | Model cache directory |
| `DEEPCELL_USE_GPU` | `true` | Enable GPU acceleration |
| `MCP_TRANSPORT` | `stdio` | MCP transport protocol |
| `MCP_PORT` | `8000` | MCP server port |

---

## Architecture

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
│  │  Tools: segment_cells, quantify_markers,         │  │
│  │         generate_segmentation_overlay             │  │
│  └──────────────────────────────────────────────────┘  │
└────────┬────────────────────────────────┬───────────────┘
         │                                │
         ▼                                ▼
┌────────────────────┐        ┌──────────────────────────┐
│  DeepCellEngine    │        │  IntensityClassifier     │
│  (segmentation)    │        │  (quantification)        │
├────────────────────┤        ├──────────────────────────┤
│ - Model loading    │        │ - Intensity measurement  │
│ - Preprocessing    │        │ - Per-cell stats         │
│ - Inference        │        │ - CSV export             │
│ - Post-processing  │        └──────────────────────────┘
└─────────┬──────────┘
          │
          ▼
┌────────────────────────────────────────────────────────┐
│           DeepCell-TF Models (TensorFlow 2.8.x)        │
│  ┌──────────────────┐    ┌──────────────────────────┐ │
│  │  Nuclear Model   │    │   Membrane Model         │ │
│  │  (~500MB)        │    │   (Mesmer, ~800MB)       │ │
│  └──────────────────┘    └──────────────────────────┘ │
└────────────────────────────────────────────────────────┘
```

---

## Related Servers

- **mcp-cell-classify** - Lightweight cell phenotype classification (downstream of this server)
- **mcp-spatialtools** - Spatial transcriptomics analysis and visualization
- **mcp-openimagedata** - Histology image processing and multiplex IF

---

**Built for the Precision Medicine MCP suite** - Enables AI-driven cell segmentation and per-cell marker quantification integrated with spatial and clinical data.
