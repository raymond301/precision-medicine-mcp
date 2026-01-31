# Synthetic Microscopy Test Data

## Overview

This directory contains synthetic microscopy images for testing mcp-deepcell.

## Files

### Image Sizes
- **512x512**: Small images for quick testing (~30 cells)
- **1024x1024**: Medium images (~120 cells)
- **2048x2048**: Large images (~480 cells)

### Markers
- **dapi_NxN.tif**: Nuclear staining (DAPI) - for segmentation
- **ki67_NxN.tif**: Proliferation marker (~25% positive cells)
- **tp53_NxN.tif**: Tumor suppressor marker (~40% positive cells)
- **membrane_NxN.tif**: Membrane marker - for Mesmer segmentation

### Metadata
- **manifest.json**: Dataset metadata and expected results

## Testing Workflow

### 1. Nuclear Segmentation
```
Tool: segment_cells
Input: dapi_512x512.tif
Model: nuclear
Expected: ~30 cells
```

### 2. Cell State Classification
```
Tool: classify_cell_states
Segmentation: From step 1
Intensity: ki67_512x512.tif
Expected: ~25% proliferating
```

### 3. Visualization
```
Tool: generate_segmentation_overlay
Tool: generate_phenotype_visualization
```

## Image Characteristics

- **Format**: 16-bit grayscale TIFF
- **Bit depth**: 0-65535 intensity range
- **Background**: 300-800 (autofluorescence)
- **Signal**: 5000-15000 (positive cells)
- **Noise**: Gaussian, σ=200-300
- **PSF Blur**: σ=1.0-1.5

## Generated

Date: 2026-01-31
Tool: generate_synthetic_data.py
Purpose: mcp-deepcell Cloud Run deployment testing
