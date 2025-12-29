#!/usr/bin/env python3
"""Test script for spatial autocorrelation (Moran's I) implementation.

This script tests the calculate_spatial_autocorrelation function with Patient-001
data to identify spatially variable genes.
"""

import os
import sys
from pathlib import Path
import asyncio

# Add src to path
src_path = Path(__file__).parent / "src"
sys.path.insert(0, str(src_path))


async def test_spatial_autocorrelation():
    """Test spatial autocorrelation with Patient-001 data."""

    # Import the server module
    from mcp_spatialtools import server

    # Ensure NOT in DRY_RUN mode
    os.environ["SPATIAL_DRY_RUN"] = "false"

    print("=" * 80)
    print("SPATIAL AUTOCORRELATION TEST - Moran's I")
    print("=" * 80)
    print()

    # Patient-001 data paths
    data_dir = Path("/Users/lynnlangit/Documents/GitHub/spatial-mcp/data/patient-data/PAT001-OVC-2025/spatial")
    expression_file = str(data_dir / "visium_gene_expression.csv")
    coordinates_file = str(data_dir / "visium_spatial_coordinates.csv")

    print(f"Data directory: {data_dir}")
    print(f"Expression file: {expression_file}")
    print(f"Coordinates file: {coordinates_file}")
    print()

    # Test 1: Proliferation markers (expected to be clustered in tumor regions)
    print("Test 1: Proliferation Markers (Expected: Clustered)")
    print("-" * 80)

    proliferation_genes = ["MKI67", "PCNA", "TOP2A", "CCND1"]
    print(f"Genes: {', '.join(proliferation_genes)}")
    print()

    result1 = await server.calculate_spatial_autocorrelation(
        expression_file=expression_file,
        coordinates_file=coordinates_file,
        genes=proliferation_genes,
        method="morans_i",
        distance_threshold=150.0  # Neighbors within 150 units
    )

    print(f"✅ Status: {result1['status']}")
    print(f"✅ Method: {result1['method']}")
    print(f"✅ Spots analyzed: {result1['num_spots']}")
    print(f"✅ Distance threshold: {result1['distance_threshold']}")
    print(f"✅ Genes analyzed: {result1['genes_analyzed']}")
    print()

    print("Summary:")
    print(f"  Significantly clustered: {result1['summary']['significantly_clustered']}")
    print(f"  Significantly dispersed: {result1['summary']['significantly_dispersed']}")
    print(f"  Random pattern: {result1['summary']['random_pattern']}")
    print()

    print("Results per gene:")
    for res in result1['results']:
        if 'morans_i' in res:
            sig_marker = "***" if res['significant'] else ""
            print(f"  {res['gene']:15s} | Moran's I: {res['morans_i']:6.3f} | "
                  f"Z-score: {res['z_score']:6.2f} | P-value: {res['p_value']:.4f} | "
                  f"{res['interpretation']:30s} {sig_marker}")
        else:
            print(f"  {res['gene']:15s} | {res['message']}")
    print()

    # Test 2: Immune markers (expected spatial patterns)
    print()
    print("=" * 80)
    print()
    print("Test 2: Immune Markers (Expected: Variable)")
    print("-" * 80)

    immune_genes = ["CD3D", "CD3E", "CD8A", "CD4"]
    print(f"Genes: {', '.join(immune_genes)}")
    print()

    result2 = await server.calculate_spatial_autocorrelation(
        expression_file=expression_file,
        coordinates_file=coordinates_file,
        genes=immune_genes,
        distance_threshold=150.0
    )

    print(f"✅ Genes analyzed: {result2['genes_analyzed']}")
    print(f"✅ Significantly clustered: {result2['summary']['significantly_clustered']}")
    print()

    print("Results:")
    for res in result2['results']:
        if 'morans_i' in res:
            sig_marker = "***" if res['significant'] else ""
            print(f"  {res['gene']:15s} | Moran's I: {res['morans_i']:6.3f} | "
                  f"Z-score: {res['z_score']:6.2f} | P-value: {res['p_value']:.4f} | "
                  f"{res['interpretation']:30s} {sig_marker}")
        else:
            print(f"  {res['gene']:15s} | {res['message']}")
    print()

    # Test 3: Stromal markers (expected to cluster in stroma)
    print()
    print("=" * 80)
    print()
    print("Test 3: Stromal/CAF Markers (Expected: Clustered in Stroma)")
    print("-" * 80)

    stromal_genes = ["FAP", "ACTA2", "COL1A1", "COL3A1", "VIM"]
    print(f"Genes: {', '.join(stromal_genes)}")
    print()

    result3 = await server.calculate_spatial_autocorrelation(
        expression_file=expression_file,
        coordinates_file=coordinates_file,
        genes=stromal_genes,
        distance_threshold=150.0
    )

    print(f"✅ Genes analyzed: {result3['genes_analyzed']}")
    print(f"✅ Significantly clustered: {result3['summary']['significantly_clustered']}")
    print()

    print("Results:")
    for res in result3['results']:
        if 'morans_i' in res:
            sig_marker = "***" if res['significant'] else ""
            print(f"  {res['gene']:15s} | Moran's I: {res['morans_i']:6.3f} | "
                  f"Z-score: {res['z_score']:6.2f} | P-value: {res['p_value']:.4f} | "
                  f"{res['interpretation']:30s} {sig_marker}")
        else:
            print(f"  {res['gene']:15s} | {res['message']}")
    print()

    print()
    print("=" * 80)
    print("✅ Spatial autocorrelation tests completed!")
    print("=" * 80)
    print()


if __name__ == "__main__":
    asyncio.run(test_spatial_autocorrelation())
