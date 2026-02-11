#!/usr/bin/env python3
"""Test script for pathway enrichment implementation.

This script tests the perform_pathway_enrichment function with example gene lists
to verify Fisher's exact test and FDR correction are working correctly.
"""

import asyncio
import os
import sys
from pathlib import Path
import pytest

# Add src to path
src_path = Path(__file__).parent / "src"
sys.path.insert(0, str(src_path))


@pytest.mark.asyncio
async def test_pathway_enrichment():
    """Test pathway enrichment with different gene lists."""

    # Import the server module
    from mcp_spatialtools import server

    # Ensure NOT in DRY_RUN mode
    os.environ["SPATIAL_DRY_RUN"] = "false"

    print("=" * 80)
    print("PATHWAY ENRICHMENT TEST")
    print("=" * 80)
    print()

    # Test 1: DNA repair genes (should enrich DNA repair pathways)
    print("Test 1: DNA Repair Gene List")
    print("-" * 80)
    dna_repair_genes = [
        "BRCA1", "BRCA2", "TP53", "ATM", "ATR", "RAD51", "XRCC1", "ERCC1", "MSH2", "MLH1"
    ]
    print(f"Genes: {', '.join(dna_repair_genes)}")
    print()

    result1 = await server.perform_pathway_enrichment.fn(
        gene_list=dna_repair_genes,
        database="GO_BP",
        p_value_cutoff=0.05
    )

    print(f"Database: {result1['database']}")
    print(f"Genes analyzed: {result1['genes_analyzed']}")
    print(f"Background size: {result1['background_size']}")
    print(f"Pathways tested: {result1['pathways_tested']}")
    print(f"Pathways enriched: {result1['pathways_enriched']}")
    print(f"Mode: {result1['mode']}")
    print()

    if result1['pathways_enriched'] > 0:
        print(f"Top enriched pathway: {result1['top_pathway']}")
        print()
        print("Top 5 enriched pathways:")
        for i, pathway in enumerate(result1['pathways'][:5], 1):
            print(f"{i}. {pathway['pathway_name']}")
            print(f"   ID: {pathway['pathway_id']}")
            print(f"   Overlap: {pathway['genes_overlapping']}/{pathway['genes_in_pathway']} genes")
            print(f"   P-value: {pathway['p_value']:.2e}")
            print(f"   P-adj (FDR): {pathway['p_adj']:.4f}")
            print(f"   Fold enrichment: {pathway['fold_enrichment']}")
            print(f"   Genes: {', '.join(pathway['overlapping_genes'][:5])}")
            if len(pathway['overlapping_genes']) > 5:
                print(f"          ... and {len(pathway['overlapping_genes']) - 5} more")
            print()
    else:
        print("⚠️  No significant pathways found")

    print()
    print("=" * 80)
    print()

    # Test 2: Cell cycle/proliferation genes
    print("Test 2: Cell Cycle & Proliferation Gene List")
    print("-" * 80)
    cell_cycle_genes = [
        "CCND1", "CCNE1", "CDK4", "CDK6", "MYC", "E2F1", "PCNA", "Ki67", "TOP2A", "AURKA"
    ]
    print(f"Genes: {', '.join(cell_cycle_genes)}")
    print()

    result2 = await server.perform_pathway_enrichment.fn(
        gene_list=cell_cycle_genes,
        database="Hallmark",
        p_value_cutoff=0.05
    )

    print(f"Database: {result2['database']}")
    print(f"Pathways enriched: {result2['pathways_enriched']}")
    print()

    if result2['pathways_enriched'] > 0:
        print(f"Top enriched pathway: {result2['top_pathway']}")
        print()
        print("Top 5 enriched pathways:")
        for i, pathway in enumerate(result2['pathways'][:5], 1):
            print(f"{i}. {pathway['pathway_name']}")
            print(f"   Overlap: {pathway['genes_overlapping']}/{pathway['genes_in_pathway']} genes")
            print(f"   P-adj (FDR): {pathway['p_adj']:.4f}")
            print(f"   Fold enrichment: {pathway['fold_enrichment']}")
            print()
    else:
        print("⚠️  No significant pathways found")

    print()
    print("=" * 80)
    print()

    # Test 3: Platinum resistance genes (Patient-001 relevant)
    print("Test 3: Platinum Resistance Gene List (Patient-001)")
    print("-" * 80)
    resistance_genes = [
        "PIK3CA", "AKT1", "PTEN", "ABCB1", "BCL2L1", "ERCC1", "GSTP1", "BRCA1"
    ]
    print(f"Genes: {', '.join(resistance_genes)}")
    print()

    result3 = await server.perform_pathway_enrichment.fn(
        gene_list=resistance_genes,
        database="Drug_Resistance",
        p_value_cutoff=0.1  # More lenient for small gene list
    )

    print(f"Database: {result3['database']}")
    print(f"Pathways enriched: {result3['pathways_enriched']}")
    print()

    if result3['pathways_enriched'] > 0:
        print(f"Top enriched pathway: {result3['top_pathway']}")
        print()
        print("Enriched pathways:")
        for i, pathway in enumerate(result3['pathways'], 1):
            print(f"{i}. {pathway['pathway_name']}")
            print(f"   Overlap: {pathway['genes_overlapping']}/{pathway['genes_in_pathway']} genes")
            print(f"   P-adj (FDR): {pathway['p_adj']:.4f}")
            print(f"   Fold enrichment: {pathway['fold_enrichment']}")
            print(f"   Genes: {', '.join(pathway['overlapping_genes'])}")
            print()
    else:
        print("⚠️  No significant pathways found")

    print()
    print("=" * 80)
    print()

    # Test 4: KEGG pathways
    print("Test 4: KEGG Pathway Database Test")
    print("-" * 80)
    pi3k_genes = ["PIK3CA", "AKT1", "AKT2", "MTOR", "PTEN", "RPS6KB1", "TSC1", "TSC2"]
    print(f"Genes: {', '.join(pi3k_genes)}")
    print()

    result4 = await server.perform_pathway_enrichment.fn(
        gene_list=pi3k_genes,
        database="KEGG",
        p_value_cutoff=0.05
    )

    print(f"Database: {result4['database']}")
    print(f"Pathways enriched: {result4['pathways_enriched']}")
    print()

    if result4['pathways_enriched'] > 0:
        print(f"Top enriched pathway: {result4['top_pathway']}")
        print()
        print("Top 3 enriched pathways:")
        for i, pathway in enumerate(result4['pathways'][:3], 1):
            print(f"{i}. {pathway['pathway_name']}")
            print(f"   ID: {pathway['pathway_id']}")
            print(f"   Overlap: {pathway['genes_overlapping']}/{pathway['genes_in_pathway']} genes")
            print(f"   P-adj (FDR): {pathway['p_adj']:.6f}")
            print(f"   Fold enrichment: {pathway['fold_enrichment']}")
            print()
    else:
        print("⚠️  No significant pathways found")

    print()
    print("=" * 80)
    print()
    print("✅ Pathway enrichment tests completed!")
    print()


if __name__ == "__main__":
    asyncio.run(test_pathway_enrichment())
