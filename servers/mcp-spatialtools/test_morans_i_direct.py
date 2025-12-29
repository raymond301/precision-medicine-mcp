#!/usr/bin/env python3
"""Direct test of Moran's I calculation without FastMCP.

Tests the _calculate_morans_i function with Patient-001 data.
"""

import sys
from pathlib import Path
import pandas as pd
import numpy as np

# Add src to path
src_path = Path(__file__).parent / "src"
sys.path.insert(0, str(src_path))

# Import the Moran's I calculation function
from mcp_spatialtools.server import _calculate_morans_i


def test_morans_i():
    """Test Moran's I with Patient-001 spatial data."""

    print("=" * 80)
    print("MORAN'S I TEST - Direct Calculation")
    print("=" * 80)
    print()

    # Load Patient-001 data
    data_dir = Path("/Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/spatial")
    expression_file = data_dir / "visium_gene_expression.csv"
    coordinates_file = data_dir / "visium_spatial_coordinates.csv"

    print(f"Loading data from: {data_dir}")
    print()

    # Load files
    expr_data = pd.read_csv(expression_file, index_col=0)
    coord_data = pd.read_csv(coordinates_file, index_col=0)

    print(f"✅ Expression: {expr_data.shape[0]} spots × {expr_data.shape[1]} genes")
    print(f"✅ Coordinates: {coord_data.shape[0]} spots")
    print()

    # Extract coordinates
    coordinates = coord_data[["array_row", "array_col"]].values

    # Test 1: Proliferation markers
    print("Test 1: Proliferation Markers")
    print("-" * 80)

    proliferation_genes = ["MKI67", "PCNA", "TOP2A", "CCND1"]
    print(f"Genes: {', '.join(proliferation_genes)}")
    print()

    print(f"{'Gene':15s} | {'Moran I':>8s} | {'Z-score':>8s} | {'P-value':>8s} | {'Interpretation':30s}")
    print("-" * 80)

    for gene in proliferation_genes:
        if gene not in expr_data.columns:
            print(f"{gene:15s} | NOT FOUND")
            continue

        expression_values = expr_data[gene].values

        # Calculate Moran's I
        morans_i, z_score, p_value = _calculate_morans_i(
            expression_values,
            coordinates,
            distance_threshold=1.5
        )

        # Interpret
        if p_value < 0.05:
            if morans_i > 0.3:
                interpretation = "significantly clustered"
            elif morans_i < -0.3:
                interpretation = "significantly dispersed"
            else:
                interpretation = "weakly patterned"
        else:
            interpretation = "random (not significant)"

        sig_marker = "***" if p_value < 0.05 else ""

        print(f"{gene:15s} | {morans_i:>8.3f} | {z_score:>8.2f} | {p_value:>8.4f} | "
              f"{interpretation:30s} {sig_marker}")

    print()
    print("=" * 80)
    print()

    # Test 2: Immune markers
    print("Test 2: Immune Markers")
    print("-" * 80)

    immune_genes = ["CD3D", "CD3E", "CD8A", "CD4"]
    print(f"Genes: {', '.join(immune_genes)}")
    print()

    print(f"{'Gene':15s} | {'Moran I':>8s} | {'Z-score':>8s} | {'P-value':>8s} | {'Interpretation':30s}")
    print("-" * 80)

    for gene in immune_genes:
        if gene not in expr_data.columns:
            print(f"{gene:15s} | NOT FOUND")
            continue

        expression_values = expr_data[gene].values
        morans_i, z_score, p_value = _calculate_morans_i(
            expression_values,
            coordinates,
            distance_threshold=1.5
        )

        if p_value < 0.05:
            if morans_i > 0.3:
                interpretation = "significantly clustered"
            elif morans_i < -0.3:
                interpretation = "significantly dispersed"
            else:
                interpretation = "weakly patterned"
        else:
            interpretation = "random (not significant)"

        sig_marker = "***" if p_value < 0.05 else ""

        print(f"{gene:15s} | {morans_i:>8.3f} | {z_score:>8.2f} | {p_value:>8.4f} | "
              f"{interpretation:30s} {sig_marker}")

    print()
    print("=" * 80)
    print()

    # Test 3: Stromal markers
    print("Test 3: Stromal/CAF Markers")
    print("-" * 80)

    stromal_genes = ["FAP", "ACTA2", "COL1A1", "COL3A1", "VIM"]
    print(f"Genes: {', '.join(stromal_genes)}")
    print()

    print(f"{'Gene':15s} | {'Moran I':>8s} | {'Z-score':>8s} | {'P-value':>8s} | {'Interpretation':30s}")
    print("-" * 80)

    for gene in stromal_genes:
        if gene not in expr_data.columns:
            print(f"{gene:15s} | NOT FOUND")
            continue

        expression_values = expr_data[gene].values
        morans_i, z_score, p_value = _calculate_morans_i(
            expression_values,
            coordinates,
            distance_threshold=1.5
        )

        if p_value < 0.05:
            if morans_i > 0.3:
                interpretation = "significantly clustered"
            elif morans_i < -0.3:
                interpretation = "significantly dispersed"
            else:
                interpretation = "weakly patterned"
        else:
            interpretation = "random (not significant)"

        sig_marker = "***" if p_value < 0.05 else ""

        print(f"{gene:15s} | {morans_i:>8.3f} | {z_score:>8.2f} | {p_value:>8.4f} | "
              f"{interpretation:30s} {sig_marker}")

    print()
    print("=" * 80)
    print()

    # Test 4: Test all genes and identify top SVGs
    print("Test 4: Identify Top Spatially Variable Genes (SVGs)")
    print("-" * 80)

    all_genes = expr_data.columns.tolist()
    print(f"Analyzing {len(all_genes)} genes...")
    print()

    svg_results = []
    for gene in all_genes:
        expression_values = expr_data[gene].values
        morans_i, z_score, p_value = _calculate_morans_i(
            expression_values,
            coordinates,
            distance_threshold=1.5
        )

        if p_value < 0.05:  # Significant
            svg_results.append({
                'gene': gene,
                'morans_i': morans_i,
                'z_score': z_score,
                'p_value': p_value
            })

    # Sort by absolute Moran's I
    svg_results.sort(key=lambda x: abs(x['morans_i']), reverse=True)

    print(f"✅ Analyzed: {len(all_genes)} genes")
    print(f"✅ Significantly variable: {len(svg_results)} genes")
    print()

    print("Top 15 Spatially Variable Genes:")
    print(f"{'Rank':>4s} | {'Gene':15s} | {'Moran I':>8s} | {'Z-score':>8s} | {'P-value':>8s} | {'Pattern':15s}")
    print("-" * 85)

    for i, res in enumerate(svg_results[:15], 1):
        pattern = "Clustered" if res['morans_i'] > 0 else "Dispersed"
        print(f"{i:>4d} | {res['gene']:15s} | {res['morans_i']:>8.3f} | {res['z_score']:>8.2f} | "
              f"{res['p_value']:>8.4f} | {pattern:15s}")

    print()
    print("=" * 80)
    print()

    # Test 5: Distance threshold sensitivity
    print("Test 5: Distance Threshold Sensitivity (TP53)")
    print("-" * 80)

    gene_test = "TP53"
    thresholds = [1.2, 1.5, 2.0, 3.0, 5.0]

    print(f"Gene: {gene_test}")
    print(f"Testing thresholds: {thresholds}")
    print()

    expression_values = expr_data[gene_test].values

    print(f"{'Threshold':>10s} | {'Moran I':>8s} | {'Z-score':>8s} | {'P-value':>8s} | {'Interpretation':30s}")
    print("-" * 80)

    for threshold in thresholds:
        morans_i, z_score, p_value = _calculate_morans_i(
            expression_values,
            coordinates,
            distance_threshold=threshold
        )

        if p_value < 0.05:
            if morans_i > 0.3:
                interpretation = "significantly clustered"
            elif morans_i < -0.3:
                interpretation = "significantly dispersed"
            else:
                interpretation = "weakly patterned"
        else:
            interpretation = "random (not significant)"

        sig_marker = "***" if p_value < 0.05 else ""

        print(f"{threshold:>10.1f} | {morans_i:>8.3f} | {z_score:>8.2f} | {p_value:>8.4f} | "
              f"{interpretation:30s} {sig_marker}")

    print()
    print("=" * 80)
    print("✅ Moran's I tests completed successfully!")
    print("=" * 80)
    print()


if __name__ == "__main__":
    test_morans_i()
