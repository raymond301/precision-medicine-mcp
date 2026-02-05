#!/usr/bin/env python3
"""
Generate synthetic microscopy imaging data for PAT002-BC-2026 breast cancer patient.
Creates realistic H&E, immunofluorescence (ER, PR, HER2, Ki67, CD8), and DAPI images.
"""

import numpy as np
from PIL import Image
from pathlib import Path
from scipy import ndimage, ndimage as ndi
from skimage import filters

def add_gaussian_noise(image, sigma=200):
    """Add Gaussian noise to simulate microscopy background"""
    noise = np.random.normal(0, sigma, image.shape)
    return np.clip(image + noise, 0, 65535).astype(np.uint16)


def create_synthetic_nuclei(size=1024, num_cells=150, seed=43):
    """Create realistic synthetic nuclear image (DAPI-like) for breast cancer.

    Args:
        size: Image dimensions (square)
        num_cells: Number of nuclei to generate
        seed: Random seed for reproducibility

    Returns:
        16-bit TIFF image with synthetic nuclei, cell positions
    """
    np.random.seed(seed)
    image = np.zeros((size, size), dtype=np.float32)

    cell_positions = []

    for i in range(num_cells):
        # Random position (avoid edges)
        x = np.random.randint(40, size - 40)
        y = np.random.randint(40, size - 40)

        # Check for overlap
        too_close = False
        for cx, cy, _ in cell_positions:
            dist = np.sqrt((x - cx)**2 + (y - cy)**2)
            if dist < 30:  # Minimum distance
                too_close = True
                break

        if too_close:
            continue

        # Nucleus parameters
        radius_x = np.random.randint(10, 18)
        radius_y = np.random.randint(10, 18)
        intensity = np.random.randint(8000, 14000)

        cell_positions.append((x, y, radius_x))

        # Create elliptical nucleus with smooth gradient
        yy, xx = np.ogrid[:size, :size]
        ellipse = ((xx - x) / radius_x)**2 + ((yy - y) / radius_y)**2
        mask = ellipse <= 1.0

        # Smooth gradient from center
        gradient = np.exp(-ellipse * 0.5) * mask
        image += gradient * intensity

        # Add chromatin texture
        num_chromatin = np.random.randint(3, 7)
        for _ in range(num_chromatin):
            cx = x + np.random.randint(-radius_x//2, radius_x//2)
            cy = y + np.random.randint(-radius_y//2, radius_y//2)
            c_radius = np.random.randint(2, 4)
            c_intensity = intensity * np.random.uniform(0.3, 0.5)

            chromatin = ((xx - cx)**2 + (yy - cy)**2) <= c_radius**2
            image[chromatin] += c_intensity

    # Add background autofluorescence
    background = np.random.uniform(300, 700, (size, size))
    image += background

    # Gaussian blur to simulate microscope PSF
    image = filters.gaussian(image, sigma=1.3)

    # Add noise
    image = add_gaussian_noise(image, sigma=250)

    # Clip to 16-bit range
    image = np.clip(image, 0, 65535).astype(np.uint16)

    return image, cell_positions


def create_nuclear_marker(cell_positions, size=1024, marker_name="ER",
                          positive_fraction=0.85, seed=44):
    """Create nuclear marker (ER, PR, Ki67) image.

    Args:
        cell_positions: List of (x, y, radius) tuples
        size: Image dimensions
        marker_name: Marker name (ER, PR, Ki67)
        positive_fraction: Fraction of cells that are marker-positive
        seed: Random seed

    Returns:
        16-bit TIFF marker image
    """
    np.random.seed(seed)
    image = np.zeros((size, size), dtype=np.float32)

    # Randomly select positive cells
    num_positive = int(len(cell_positions) * positive_fraction)
    positive_indices = np.random.choice(len(cell_positions), num_positive, replace=False)

    for i, (x, y, radius) in enumerate(cell_positions):
        is_positive = i in positive_indices

        if is_positive:
            # High intensity for positive cells
            intensity = np.random.randint(9000, 15000)
            radius_mult = 0.9  # Nuclear staining
        else:
            # Low intensity for negative cells
            intensity = np.random.randint(400, 1500)
            radius_mult = 0.7

        # Create circular nuclear marker signal
        yy, xx = np.ogrid[:size, :size]
        circle = (xx - x)**2 + (yy - y)**2
        marker_radius = radius * radius_mult
        mask = circle <= marker_radius**2

        # Smooth gradient
        gradient = np.exp(-circle / (2 * (marker_radius * 0.6)**2)) * mask
        image += gradient * intensity

    # Background
    background = np.random.uniform(250, 550, (size, size))
    image += background

    # Blur
    image = filters.gaussian(image, sigma=1.4)

    # Noise
    image = add_gaussian_noise(image, sigma=220)

    # Clip to 16-bit
    image = np.clip(image, 0, 65535).astype(np.uint16)

    return image


def create_membrane_marker(cell_positions, size=1024, marker_name="HER2",
                           positive_fraction=0.05, seed=45):
    """Create membrane marker (HER2, CD8) image.

    Args:
        cell_positions: List of (x, y, radius) tuples
        size: Image dimensions
        marker_name: Marker name
        positive_fraction: Fraction of cells that are marker-positive
        seed: Random seed

    Returns:
        16-bit TIFF membrane marker image
    """
    np.random.seed(seed)
    image = np.zeros((size, size), dtype=np.float32)

    # For HER2-negative breast cancer, all cells have low signal
    # For CD8, only scattered T cells are positive

    if marker_name == "HER2":
        # HER2 negative (IHC 0): weak membrane staining all cells
        positive_intensity_range = (800, 2000)  # Very low
    elif marker_name == "CD8":
        # Scattered CD8+ T cells
        positive_intensity_range = (8000, 14000)  # High for positive cells
    else:
        positive_intensity_range = (6000, 12000)

    # Select positive cells
    num_positive = int(len(cell_positions) * positive_fraction)
    positive_indices = np.random.choice(len(cell_positions), num_positive, replace=False)

    for i, (x, y, radius) in enumerate(cell_positions):
        is_positive = i in positive_indices

        if is_positive or marker_name == "HER2":
            # Membrane ring
            membrane_radius = radius + np.random.randint(2, 5)
            thickness = np.random.randint(2, 4)

            if is_positive and marker_name != "HER2":
                intensity = np.random.randint(*positive_intensity_range)
            else:
                # HER2 negative or CD8 negative
                intensity = np.random.randint(500, 1800)

            # Create ring (membrane)
            yy, xx = np.ogrid[:size, :size]
            dist = np.sqrt((xx - x)**2 + (yy - y)**2)
            ring = (dist >= membrane_radius - thickness) & (dist <= membrane_radius + thickness)

            image[ring] += intensity

    # Background
    background = np.random.uniform(300, 600, (size, size))
    image += background

    # Blur
    image = filters.gaussian(image, sigma=1.2)

    # Noise
    image = add_gaussian_noise(image, sigma=250)

    # Clip
    image = np.clip(image, 0, 65535).astype(np.uint16)

    return image


def create_synthetic_he(size=1024, num_cells=150, seed=46):
    """Create synthetic H&E histology image.

    H&E staining:
    - Hematoxylin (purple/blue): Nuclei
    - Eosin (pink): Cytoplasm, stroma

    Args:
        size: Image dimensions
        num_cells: Number of cells
        seed: Random seed

    Returns:
        RGB 8-bit image
    """
    np.random.seed(seed)

    # Create RGB channels
    r_channel = np.ones((size, size), dtype=np.float32) * 240  # Base pink/red
    g_channel = np.ones((size, size), dtype=np.float32) * 200  # Base pink/green
    b_channel = np.ones((size, size), dtype=np.float32) * 230  # Base pink/blue

    # Add tissue texture (collagen, cytoplasm)
    tissue_texture = np.random.uniform(0.8, 1.0, (size, size))
    r_channel *= tissue_texture
    g_channel *= tissue_texture * 0.85
    b_channel *= tissue_texture * 0.95

    # Generate nuclei (hematoxylin staining - purple/blue)
    for i in range(num_cells):
        x = np.random.randint(30, size - 30)
        y = np.random.randint(30, size - 30)
        radius = np.random.randint(8, 14)

        # Create nucleus
        yy, xx = np.ogrid[:size, :size]
        ellipse = ((xx - x) / radius)**2 + ((yy - y) / (radius * 0.9))**2
        mask = ellipse <= 1.0

        # Hematoxylin: reduce R/G, increase B (purple/blue)
        nucleus_intensity = np.random.uniform(0.15, 0.35)
        gradient = np.exp(-ellipse * 0.5) * mask

        r_channel -= gradient * (240 * nucleus_intensity)
        g_channel -= gradient * (200 * nucleus_intensity)
        b_channel += gradient * (150 * nucleus_intensity * 0.3)  # Slight blue increase

        # Cytoplasm around nucleus (eosin - pink)
        cyto_radius = radius + np.random.randint(3, 8)
        cyto_mask = ((xx - x) / cyto_radius)**2 + ((yy - y) / cyto_radius)**2 <= 1.0
        cyto_mask = cyto_mask & ~mask  # Ring around nucleus

        cyto_intensity = np.random.uniform(0.9, 1.0)
        r_channel[cyto_mask] *= cyto_intensity
        g_channel[cyto_mask] *= cyto_intensity * 0.8
        b_channel[cyto_mask] *= cyto_intensity * 0.85

    # Add stromal features (fibrous texture)
    for _ in range(20):
        x1, y1 = np.random.randint(0, size, 2)
        x2, y2 = x1 + np.random.randint(-50, 50), y1 + np.random.randint(-50, 50)

        # Draw line (collagen fiber)
        rr, cc = np.meshgrid(range(size), range(size), indexing='ij')
        line_dist = np.abs((y2-y1)*rr - (x2-x1)*cc + x2*y1 - y2*x1) / np.sqrt((y2-y1)**2 + (x2-x1)**2 + 1e-10)
        fiber = line_dist < 3

        fiber_color = np.random.uniform(1.0, 1.1)
        r_channel[fiber] *= fiber_color * 0.95
        g_channel[fiber] *= fiber_color * 0.9
        b_channel[fiber] *= fiber_color

    # Clip and convert to 8-bit RGB
    r_channel = np.clip(r_channel, 0, 255).astype(np.uint8)
    g_channel = np.clip(g_channel, 0, 255).astype(np.uint8)
    b_channel = np.clip(b_channel, 0, 255).astype(np.uint8)

    # Stack channels
    rgb_image = np.stack([r_channel, g_channel, b_channel], axis=-1)

    return rgb_image


def save_tiff_16bit(image, filepath):
    """Save 16-bit TIFF image"""
    Image.fromarray(image.astype(np.uint16)).save(filepath)
    print(f"   ✓ Saved: {filepath.name}")


def save_tiff_rgb(image, filepath):
    """Save 8-bit RGB TIFF image"""
    Image.fromarray(image.astype(np.uint8), mode='RGB').save(filepath)
    print(f"   ✓ Saved: {filepath.name}")


def main():
    """Generate all imaging files for PAT002-BC-2026."""
    print("="*70)
    print("Generating Synthetic Microscopy Imaging Data")
    print("Patient: PAT002-BC-2026 (ER+/PR+/HER2- Breast Cancer)")
    print("="*70)

    output_dir = Path(__file__).parent
    he_dir = output_dir / "he"
    if_dir = output_dir / "if"
    he_dir.mkdir(exist_ok=True)
    if_dir.mkdir(exist_ok=True)

    # Image sizes and cell counts
    configs = [
        {"size": 512, "num_cells": 80, "suffix": ""},
        {"size": 1024, "num_cells": 200, "suffix": ""},
        {"size": 2048, "num_cells": 450, "suffix": "_high"},
    ]

    for config in configs:
        size = config["size"]
        num_cells = config["num_cells"]
        suffix = config["suffix"]

        print(f"\n{'='*70}")
        print(f"Generating {size}x{size} images ({num_cells} cells)")
        print(f"{'='*70}")

        # Generate nuclei (DAPI)
        print(f"\n1. Nuclear staining (DAPI)")
        dapi_img, cell_positions = create_synthetic_nuclei(size, num_cells, seed=43+size)

        if size == 1024:
            # Main directory files (1024x1024)
            dapi_path = output_dir / "PAT002_tumor_IF_DAPI.tiff"
            save_tiff_16bit(dapi_img, dapi_path)

            # ER marker (85% positive, nuclear)
            print(f"\n2. ER (Estrogen Receptor) - 85% positive")
            er_img = create_nuclear_marker(cell_positions, size, "ER",
                                          positive_fraction=0.85, seed=44)
            er_path = output_dir / "PAT002_tumor_IF_ER.tiff"
            save_tiff_16bit(er_img, er_path)

            # PR marker (70% positive, nuclear)
            print(f"\n3. PR (Progesterone Receptor) - 70% positive")
            pr_img = create_nuclear_marker(cell_positions, size, "PR",
                                          positive_fraction=0.70, seed=45)
            pr_path = output_dir / "PAT002_tumor_IF_PR.tiff"
            save_tiff_16bit(pr_img, pr_path)

            # HER2 marker (negative, weak membrane all cells)
            print(f"\n4. HER2 - Negative (IHC 0)")
            her2_img = create_membrane_marker(cell_positions, size, "HER2",
                                             positive_fraction=1.0, seed=46)
            her2_path = output_dir / "PAT002_tumor_IF_HER2.tiff"
            save_tiff_16bit(her2_img, her2_path)

            # Ki67 marker (20% positive, nuclear, proliferation)
            print(f"\n5. Ki67 (Proliferation) - 20% positive")
            ki67_img = create_nuclear_marker(cell_positions, size, "Ki67",
                                            positive_fraction=0.20, seed=47)
            ki67_path = output_dir / "PAT002_tumor_IF_KI67.tiff"
            save_tiff_16bit(ki67_img, ki67_path)

            # CD8 marker (5-10% positive, membrane, T cells)
            print(f"\n6. CD8 (T cells) - 8% positive")
            cd8_img = create_membrane_marker(cell_positions, size, "CD8",
                                            positive_fraction=0.08, seed=48)
            cd8_path = output_dir / "PAT002_tumor_IF_CD8.tiff"
            save_tiff_16bit(cd8_img, cd8_path)

            # H&E staining
            print(f"\n7. H&E Histology")
            he_img = create_synthetic_he(size, num_cells, seed=49)
            he_path = output_dir / "PAT002_tumor_HE_20x.tiff"
            save_tiff_rgb(he_img, he_path)

        elif size == 2048:
            # High-resolution files for he/ and if/ subdirectories
            print(f"\nHigh-resolution images:")

            # H&E high-resolution
            print(f"  - H&E high-resolution")
            he_high = create_synthetic_he(size, num_cells, seed=50)
            he_high_path = he_dir / "PAT002_tumor_HE_20x_high.tif"
            save_tiff_rgb(he_high, he_high_path)

            he_high_path2 = he_dir / "sample_002_high.tif"
            save_tiff_rgb(he_high, he_high_path2)

            # H&E low (resampled from high)
            print(f"  - H&E low-resolution (downsampled)")
            he_low = he_high[::4, ::4]  # Downsample to 512x512
            he_low_path = he_dir / "sample_002_low.tif"
            save_tiff_rgb(he_low, he_low_path)

            # ER high-resolution
            print(f"  - ER high-resolution")
            dapi_high, positions_high = create_synthetic_nuclei(size, num_cells, seed=51)
            er_high = create_nuclear_marker(positions_high, size, "ER",
                                           positive_fraction=0.85, seed=52)
            er_high_path = if_dir / "PAT002_tumor_IF_ER_high.tif"
            save_tiff_16bit(er_high, er_high_path)

            # Ki67 high-resolution
            print(f"  - Ki67 high-resolution")
            ki67_high = create_nuclear_marker(positions_high, size, "Ki67",
                                             positive_fraction=0.20, seed=53)
            ki67_high_path = if_dir / "PAT002_tumor_IF_KI67_high.tif"
            save_tiff_16bit(ki67_high, ki67_high_path)

    # Summary
    print("\n" + "="*70)
    print("Imaging Data Generation Complete")
    print("="*70)

    # Count files
    main_files = list(output_dir.glob("*.tiff")) + list(output_dir.glob("*.tif"))
    he_files = list(he_dir.glob("*.tif"))
    if_files = list(if_dir.glob("*.tif"))

    print(f"\nFiles generated:")
    print(f"  Main directory: {len([f for f in main_files if f.is_file()])} files")
    print(f"    - PAT002_tumor_HE_20x.tiff (H&E, RGB 8-bit)")
    print(f"    - PAT002_tumor_IF_*.tiff (6 markers, 16-bit)")
    print(f"  he/ subdirectory: {len(he_files)} files (high/low resolution)")
    print(f"  if/ subdirectory: {len(if_files)} files (high resolution)")
    print(f"\nTotal files: {len(main_files) + len(he_files) + len(if_files)}")

    print(f"\nBreast cancer imaging features:")
    print(f"  - ER: 85% positive (strong nuclear staining)")
    print(f"  - PR: 70% positive (moderate nuclear staining)")
    print(f"  - HER2: Negative (IHC 0, weak membrane all cells)")
    print(f"  - Ki67: 20% positive (post-treatment proliferation)")
    print(f"  - CD8: 8% positive (tumor-infiltrating lymphocytes)")
    print(f"  - H&E: Synthetic histology (pink eosin, purple hematoxylin)")
    print(f"  - DAPI: Nuclear staining (all cells)")

    print("\n✓ Generation successful!")


if __name__ == "__main__":
    main()
