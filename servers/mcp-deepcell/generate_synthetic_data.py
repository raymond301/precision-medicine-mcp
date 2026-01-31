#!/usr/bin/env python3
"""
Generate synthetic microscopy test data for mcp-deepcell testing
Creates realistic DAPI nuclear images and marker images (Ki67, TP53, etc.)
"""

import numpy as np
from PIL import Image
from pathlib import Path
import json
from scipy import ndimage
from skimage import filters, morphology

def add_gaussian_noise(image, sigma=200):
    """Add Gaussian noise to simulate microscopy background"""
    noise = np.random.normal(0, sigma, image.shape)
    return np.clip(image + noise, 0, 65535).astype(np.uint16)

def create_synthetic_nuclei(size=512, num_cells=30, seed=42):
    """
    Create realistic synthetic nuclear image (DAPI-like)

    Args:
        size: Image dimensions (square)
        num_cells: Number of nuclei to generate
        seed: Random seed for reproducibility

    Returns:
        16-bit TIFF image with synthetic nuclei
    """
    np.random.seed(seed)
    image = np.zeros((size, size), dtype=np.float32)

    # Store cell positions for later
    cell_positions = []

    for i in range(num_cells):
        # Random position (avoid edges)
        x = np.random.randint(40, size - 40)
        y = np.random.randint(40, size - 40)

        # Check for overlap with existing cells
        too_close = False
        for cx, cy in cell_positions:
            dist = np.sqrt((x - cx)**2 + (y - cy)**2)
            if dist < 35:  # Minimum distance between nuclei
                too_close = True
                break

        if too_close:
            continue

        cell_positions.append((x, y))

        # Nucleus parameters
        radius_x = np.random.randint(12, 20)
        radius_y = np.random.randint(12, 20)
        intensity = np.random.randint(8000, 15000)

        # Create elliptical nucleus with smooth gradient
        yy, xx = np.ogrid[:size, :size]
        ellipse = ((xx - x) / radius_x)**2 + ((yy - y) / radius_y)**2
        mask = ellipse <= 1.0

        # Smooth gradient from center
        gradient = np.exp(-ellipse * 0.5) * mask
        image += gradient * intensity

        # Add chromatin texture (brighter spots within nucleus)
        num_chromatin = np.random.randint(3, 8)
        for _ in range(num_chromatin):
            cx = x + np.random.randint(-radius_x//2, radius_x//2)
            cy = y + np.random.randint(-radius_y//2, radius_y//2)
            c_radius = np.random.randint(2, 5)
            c_intensity = intensity * np.random.uniform(0.3, 0.6)

            chromatin = ((xx - cx)**2 + (yy - cy)**2) <= c_radius**2
            image[chromatin] += c_intensity

    # Add background autofluorescence
    background = np.random.uniform(300, 800, (size, size))
    image += background

    # Add Gaussian blur to simulate microscope PSF
    image = filters.gaussian(image, sigma=1.2)

    # Add realistic noise
    image = add_gaussian_noise(image, sigma=250)

    # Clip to 16-bit range
    image = np.clip(image, 0, 65535).astype(np.uint16)

    return image, cell_positions

def create_synthetic_marker(nuclear_image, cell_positions, marker_name="Ki67",
                           positive_fraction=0.25, seed=43):
    """
    Create synthetic marker image (Ki67, TP53, etc.) based on nuclear positions

    Args:
        nuclear_image: Reference nuclear image for cell locations
        cell_positions: List of (x, y) cell center positions
        marker_name: Name of marker (affects intensity distribution)
        positive_fraction: Fraction of cells that are marker-positive
        seed: Random seed

    Returns:
        16-bit TIFF marker image
    """
    np.random.seed(seed)
    size = nuclear_image.shape[0]
    image = np.zeros((size, size), dtype=np.float32)

    # Randomly select positive cells
    num_positive = int(len(cell_positions) * positive_fraction)
    positive_indices = np.random.choice(len(cell_positions), num_positive, replace=False)

    for i, (x, y) in enumerate(cell_positions):
        is_positive = i in positive_indices

        # Determine intensity based on marker status
        if is_positive:
            # High intensity for positive cells
            intensity = np.random.randint(7000, 14000)
            radius = np.random.randint(10, 16)
        else:
            # Low intensity for negative cells
            intensity = np.random.randint(500, 2500)
            radius = np.random.randint(8, 14)

        # Create circular marker signal
        yy, xx = np.ogrid[:size, :size]
        circle = (xx - x)**2 + (yy - y)**2
        mask = circle <= radius**2

        # Smooth gradient
        gradient = np.exp(-circle / (2 * (radius * 0.5)**2)) * mask
        image += gradient * intensity

    # Add background
    background = np.random.uniform(200, 600, (size, size))
    image += background

    # Blur
    image = filters.gaussian(image, sigma=1.5)

    # Add noise
    image = add_gaussian_noise(image, sigma=200)

    # Clip to 16-bit
    image = np.clip(image, 0, 65535).astype(np.uint16)

    return image

def create_synthetic_membrane(size=512, num_cells=25, seed=44):
    """
    Create synthetic membrane marker image for Mesmer model testing

    Args:
        size: Image dimensions
        num_cells: Number of cells
        seed: Random seed

    Returns:
        16-bit TIFF membrane image
    """
    np.random.seed(seed)
    image = np.zeros((size, size), dtype=np.float32)

    for i in range(num_cells):
        x = np.random.randint(50, size - 50)
        y = np.random.randint(50, size - 50)
        radius = np.random.randint(18, 28)
        thickness = np.random.randint(2, 4)
        intensity = np.random.randint(6000, 12000)

        # Create ring (membrane)
        yy, xx = np.ogrid[:size, :size]
        dist = np.sqrt((xx - x)**2 + (yy - y)**2)
        ring = (dist >= radius - thickness) & (dist <= radius + thickness)

        image[ring] += intensity

    # Background
    background = np.random.uniform(300, 700, (size, size))
    image += background

    # Blur
    image = filters.gaussian(image, sigma=1.0)

    # Noise
    image = add_gaussian_noise(image, sigma=300)

    # Clip
    image = np.clip(image, 0, 65535).astype(np.uint16)

    return image

def save_tiff(image, filepath):
    """Save 16-bit TIFF image"""
    Image.fromarray(image.astype(np.uint16)).save(filepath)
    print(f"✅ Saved: {filepath}")

def main():
    """Generate synthetic test dataset"""
    output_dir = Path("./test_data")
    output_dir.mkdir(exist_ok=True)

    print("=" * 60)
    print("Generating Synthetic Microscopy Test Data")
    print("=" * 60)

    # Create different image sizes
    sizes = [512, 1024, 2048]
    datasets = []

    for size in sizes:
        print(f"\n📊 Generating {size}x{size} images...")

        # Adjust cell count based on size
        num_cells = int(30 * (size / 512)**2)

        # Nuclear image (DAPI)
        print(f"  Creating nuclear image ({num_cells} cells)...")
        nuclear_img, positions = create_synthetic_nuclei(size, num_cells, seed=42+size)
        nuclear_path = output_dir / f"dapi_{size}x{size}.tif"
        save_tiff(nuclear_img, nuclear_path)

        # Ki67 marker (proliferation)
        print(f"  Creating Ki67 marker (25% positive)...")
        ki67_img = create_synthetic_marker(nuclear_img, positions, "Ki67",
                                          positive_fraction=0.25, seed=43+size)
        ki67_path = output_dir / f"ki67_{size}x{size}.tif"
        save_tiff(ki67_img, ki67_path)

        # TP53 marker (tumor suppressor)
        print(f"  Creating TP53 marker (40% positive)...")
        tp53_img = create_synthetic_marker(nuclear_img, positions, "TP53",
                                          positive_fraction=0.40, seed=44+size)
        tp53_path = output_dir / f"tp53_{size}x{size}.tif"
        save_tiff(tp53_img, tp53_path)

        # Membrane marker
        print(f"  Creating membrane marker...")
        membrane_img = create_synthetic_membrane(size, num_cells, seed=45+size)
        membrane_path = output_dir / f"membrane_{size}x{size}.tif"
        save_tiff(membrane_img, membrane_path)

        # Create metadata
        dataset = {
            "size": f"{size}x{size}",
            "num_cells": num_cells,
            "files": {
                "nuclear": nuclear_path.name,
                "ki67": ki67_path.name,
                "tp53": tp53_path.name,
                "membrane": membrane_path.name
            },
            "expected_results": {
                "segmentation": {
                    "approximate_cell_count": num_cells,
                    "min_cell_size": 100
                },
                "classification": {
                    "ki67_positive_fraction": 0.25,
                    "tp53_positive_fraction": 0.40
                }
            }
        }
        datasets.append(dataset)

    # Save manifest
    manifest = {
        "description": "Synthetic microscopy test dataset for mcp-deepcell",
        "generated": "2026-01-31",
        "image_format": "16-bit TIFF",
        "pixel_type": "grayscale",
        "datasets": datasets
    }

    manifest_path = output_dir / "manifest.json"
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)
    print(f"\n✅ Saved manifest: {manifest_path}")

    # Create README
    readme_path = output_dir / "README.md"
    with open(readme_path, 'w') as f:
        f.write("""# Synthetic Microscopy Test Data

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
""")
    print(f"✅ Saved README: {readme_path}")

    print("\n" + "=" * 60)
    print("✅ Synthetic Data Generation Complete")
    print("=" * 60)
    print(f"\nOutput directory: {output_dir.absolute()}")
    print(f"Total files: {len(list(output_dir.glob('*.tif')))} TIFFs + manifest + README")
    print(f"\nImage sizes: {', '.join([f'{s}x{s}' for s in sizes])}")
    print(f"Total dataset size: ~{sum([s*s*2*4 for s in sizes])/1024/1024:.1f} MB")
    print("\n📦 Ready to upload to GCS!")

if __name__ == "__main__":
    main()
