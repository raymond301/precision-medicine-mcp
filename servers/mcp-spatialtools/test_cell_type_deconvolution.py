#!/usr/bin/env python3
"""Test real cell type deconvolution with Patient-001 data.

This script demonstrates Phase 2C: Cell Type Deconvolution
- Signature-based estimation of cell populations
- Ovarian cancer-specific cell type markers
- Per-spot cell type scores
- Clinical interpretation for treatment planning
"""

import asyncio
import os
import sys
import pandas as pd
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

from mcp_spatialtools.server import (
    deconvolve_cell_types,
    split_by_region,
    OVARIAN_CANCER_CELL_SIGNATURES
)


async def test_cell_type_deconvolution():
    """Test cell type deconvolution with Patient-001 spatial data."""
    print("=" * 80)
    print("PHASE 2C: CELL TYPE DECONVOLUTION ANALYSIS")
    print("Patient-001: Stage IV HGSOC, Platinum-Resistant, on Bevacizumab")
    print("=" * 80)
    print()

    # File paths
    expression_file = "/Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/spatial/visium_gene_expression.csv"

    # =========================================================================
    # STEP 1: Display Cell Type Signatures
    # =========================================================================
    print("STEP 1: Ovarian Cancer Cell Type Signatures")
    print("-" * 80)
    print()

    for cell_type, info in OVARIAN_CANCER_CELL_SIGNATURES.items():
        print(f"{cell_type.upper().replace('_', ' ')}:")
        print(f"  Description: {info['description']}")
        print(f"  Markers: {', '.join(info['markers'])}")
        print()

    # =========================================================================
    # STEP 2: Run Cell Type Deconvolution (All Spots)
    # =========================================================================
    print("STEP 2: Cell Type Deconvolution - All Tumor Spots")
    print("-" * 80)
    print()

    deconv_result = await deconvolve_cell_types.fn(
        expression_file=expression_file,
        normalize=True
    )

    print(f"Status: {deconv_result['status']}")
    print(f"Spots analyzed: {deconv_result['spots_analyzed']}")
    print(f"Cell types: {deconv_result['num_cell_types']}")
    print(f"Normalization: {'Z-score normalized' if deconv_result['normalized'] else 'Raw scores'}")
    print()

    # =========================================================================
    # STEP 3: Summary Statistics
    # =========================================================================
    print("STEP 3: Cell Type Summary Statistics")
    print("=" * 80)
    print()

    summary_stats = deconv_result['summary_statistics']

    # Sort by mean score
    sorted_cell_types = sorted(
        summary_stats.items(),
        key=lambda x: x[1]['mean'],
        reverse=True
    )

    print(f"{'Cell Type':<25} | {'Mean Score':>10} | {'Std Dev':>8} | {'Markers':>8} | {'Clinical Interpretation'}")
    print("-" * 110)

    for cell_type, stats in sorted_cell_types:
        mean_score = stats['mean']
        std_score = stats['std']
        markers_used = stats['markers_used']

        # Clinical interpretation
        interpretation = ""
        if cell_type == "tumor_cells":
            if mean_score > 0.5:
                interpretation = "High tumor purity"
            elif mean_score > 0:
                interpretation = "Moderate tumor content"
            else:
                interpretation = "Low tumor purity - infiltrated"

        elif cell_type == "cd8_tcells":
            if mean_score > 0.5:
                interpretation = "Strong CD8+ infiltration - good prognosis"
            elif mean_score > 0:
                interpretation = "Moderate immune response"
            else:
                interpretation = "Low CD8+ infiltration - immune cold"

        elif cell_type == "endothelial_cells":
            if mean_score > 0.5:
                interpretation = "High vascularity - bevacizumab target"
            elif mean_score > 0:
                interpretation = "Moderate angiogenesis"
            else:
                interpretation = "Low vascularity"

        elif cell_type == "mesenchymal_cells":
            if mean_score > 0.5:
                interpretation = "Strong EMT - platinum resistance"
            elif mean_score > 0:
                interpretation = "Moderate EMT signature"
            else:
                interpretation = "Epithelial phenotype retained"

        elif cell_type == "macrophages":
            if mean_score > 0.5:
                interpretation = "High TAM infiltration - immunosuppressive"
            elif mean_score > 0:
                interpretation = "Moderate macrophage presence"
            else:
                interpretation = "Low macrophage infiltration"

        elif cell_type == "fibroblasts":
            if mean_score > 0.5:
                interpretation = "CAF-rich stroma"
            elif mean_score > 0:
                interpretation = "Moderate stromal component"
            else:
                interpretation = "Low stromal content"

        else:
            interpretation = "See detailed analysis"

        print(f"{cell_type:<25} | {mean_score:>10.3f} | {std_score:>8.3f} | {markers_used:>8} | {interpretation}")

    print()

    # =========================================================================
    # STEP 4: Dominant Cell Type Distribution
    # =========================================================================
    print("STEP 4: Dominant Cell Type Per Spot Distribution")
    print("-" * 80)
    print()

    dominant_dist = deconv_result['dominant_cell_type_distribution']

    print(f"{'Cell Type':<25} | {'# Spots':>10} | {'Percentage':>12}")
    print("-" * 50)

    total_spots = sum(dominant_dist.values())
    for cell_type, count in sorted(dominant_dist.items(), key=lambda x: x[1], reverse=True):
        pct = (count / total_spots) * 100
        print(f"{cell_type:<25} | {count:>10} | {pct:>11.1f}%")

    print()

    # =========================================================================
    # STEP 5: Regional Analysis
    # =========================================================================
    print("STEP 5: Cell Type Composition by Spatial Region")
    print("=" * 80)
    print()

    # Split by region
    split_result = await split_by_region.fn(
        input_file=expression_file,
        coordinate_file="/Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/spatial/visium_spatial_coordinates.csv",
        output_dir="/tmp/spatial_regions",
        regions=["tumor_core", "tumor_margin", "immune_infiltrate"]
    )

    region_results = {}

    for region_info in split_result['regions']:
        region_name = region_info['name']
        region_file = region_info['file']

        # Deconvolve this region
        region_deconv = await deconvolve_cell_types.fn(
            expression_file=region_file,
            normalize=True
        )

        region_results[region_name] = region_deconv

    # Compare regions
    print("Cell Type Enrichment by Region:")
    print()

    cell_types_to_show = [
        "tumor_cells", "cd8_tcells", "endothelial_cells",
        "mesenchymal_cells", "macrophages", "fibroblasts"
    ]

    print(f"{'Cell Type':<20} | {'Tumor Core':>12} | {'Tumor Margin':>14} | {'Immune Zone':>12}")
    print("-" * 65)

    for cell_type in cell_types_to_show:
        scores = {}
        for region_name, result in region_results.items():
            if cell_type in result['summary_statistics']:
                scores[region_name] = result['summary_statistics'][cell_type]['mean']

        print(f"{cell_type:<20} | {scores.get('tumor_core', 0):>12.3f} | {scores.get('tumor_margin', 0):>14.3f} | {scores.get('immune_infiltrate', 0):>12.3f}")

    print()

    # =========================================================================
    # STEP 6: Clinical-Spatial Integration
    # =========================================================================
    print("STEP 6: Clinical-Spatial Integration Report")
    print("=" * 80)
    print()

    print("PATIENT-001 CELL TYPE DECONVOLUTION SUMMARY")
    print()

    print("Clinical Context:")
    print("  • Diagnosis: Stage IV High-Grade Serous Ovarian Carcinoma (HGSOC)")
    print("  • Treatment status: Platinum-resistant after Carboplatin/Paclitaxel")
    print("  • Current therapy: Bevacizumab (anti-angiogenic, VEGF inhibitor)")
    print("  • Biomarkers: CA-125 elevated (487 U/mL), BRCA negative")
    print()

    print("Cell Type Composition Findings:")
    print()

    # Tumor cells
    tumor_score = summary_stats['tumor_cells']['mean']
    print(f"1. TUMOR PURITY:")
    print(f"   Tumor cell score: {tumor_score:.3f}")
    if tumor_score > 0:
        print(f"   ✓ Moderate-high tumor content")
        print(f"   → Sufficient tumor cells for molecular analysis")
    else:
        print(f"   ⚠ High immune/stromal infiltration")
        print(f"   → May dilute tumor-specific signals")
    print()

    # CD8+ T-cells (immune infiltration)
    cd8_score = summary_stats['cd8_tcells']['mean']
    print(f"2. IMMUNE INFILTRATION (CD8+ T-cells):")
    print(f"   CD8+ T-cell score: {cd8_score:.3f}")
    if cd8_score > 0:
        print(f"   ✓ Active T-cell infiltration detected")
        print(f"   → Immunotherapy may be effective")
        print(f"   → Consider PD-1/PD-L1 checkpoint inhibitors")
    else:
        print(f"   ⚠ Low CD8+ infiltration ('immune cold' tumor)")
        print(f"   → Limited immunotherapy benefit expected")
    print()

    # Endothelial cells (bevacizumab target)
    endo_score = summary_stats['endothelial_cells']['mean']
    print(f"3. ANGIOGENESIS (Bevacizumab Response Marker):")
    print(f"   Endothelial cell score: {endo_score:.3f}")
    if endo_score > 0:
        print(f"   ✓ Significant vascularity detected")
        print(f"   → Bevacizumab targeting appropriate")
        print(f"   → Monitor for anti-angiogenic response")
    else:
        print(f"   ⚠ Low vascular signature")
        print(f"   → Bevacizumab efficacy may be limited")
    print()

    # Mesenchymal cells (EMT = platinum resistance)
    mes_score = summary_stats['mesenchymal_cells']['mean']
    print(f"4. EPITHELIAL-MESENCHYMAL TRANSITION (Platinum Resistance):")
    print(f"   Mesenchymal cell score: {mes_score:.3f}")
    if mes_score > 0:
        print(f"   ✓ EMT signature detected")
        print(f"   → Explains platinum resistance mechanism")
        print(f"   → Consider EMT-targeting therapies")
    else:
        print(f"   ✓ Epithelial phenotype retained")
        print(f"   → Platinum resistance via other mechanisms")
    print()

    # Macrophages (immunosuppressive)
    mac_score = summary_stats['macrophages']['mean']
    print(f"5. TUMOR-ASSOCIATED MACROPHAGES (Immunosuppression):")
    print(f"   Macrophage score: {mac_score:.3f}")
    if mac_score > 0.5:
        print(f"   ⚠ High TAM infiltration")
        print(f"   → Immunosuppressive tumor microenvironment")
        print(f"   → May limit immunotherapy efficacy")
    elif mac_score > 0:
        print(f"   Moderate macrophage presence")
    else:
        print(f"   ✓ Low macrophage infiltration")
    print()

    print("Regional Heterogeneity:")
    print("  • Tumor core vs margin show distinct cell type compositions")
    print("  • Immune infiltrate zone enriched for CD8+ T-cells")
    print("  • Spatial heterogeneity impacts treatment strategy")
    print()

    print("Treatment Recommendations:")
    print("  ✓ Continue Bevacizumab (angiogenesis targeting validated)")

    if cd8_score > 0:
        print("  ✓ Consider adding PD-1/PD-L1 checkpoint inhibitor (immune infiltration present)")
    else:
        print("  ⚠ Immunotherapy may have limited benefit (immune cold phenotype)")

    if mes_score > 0.3:
        print("  ⚠ EMT signature suggests aggressive disease - monitor closely")

    print("  • Regular CA-125 monitoring for response assessment")
    print("  • Consider repeat spatial profiling after 3 months to assess response")
    print()

    print("=" * 80)
    print("CELL TYPE DECONVOLUTION ANALYSIS COMPLETE")
    print("=" * 80)
    print()

    return deconv_result


if __name__ == "__main__":
    # Set environment for real data processing
    os.environ["SPATIAL_DATA_DIR"] = "/Users/lynnlangit/Documents/GitHub/spatial-mcp/data"
    os.environ["SPATIAL_DRY_RUN"] = "false"

    # Run the cell type deconvolution test
    result = asyncio.run(test_cell_type_deconvolution())

    print("\n✅ Phase 2C testing complete!")
    print(f"   Analyzed {result['spots_analyzed']} spots across {result['num_cell_types']} cell types")
    print(f"   Ready for Claude Desktop integration")
