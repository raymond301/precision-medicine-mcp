"""Generate synthetic Visium spatial transcriptomics data for PAT002-BC-2026 breast cancer."""

import numpy as np
import pandas as pd
from pathlib import Path

# Set random seed for reproducibility
np.random.seed(43)

# Grid dimensions (10x Visium standard)
GRID_SIZE = 30  # 30x30 = 900 spots
N_SPOTS = GRID_SIZE * GRID_SIZE

# Breast cancer-specific gene panel (35 genes)
BREAST_CANCER_GENES = [
    # Hormone receptors and luminal markers
    "ESR1", "PGR", "GATA3", "FOXA1",

    # Proliferation markers
    "MKI67", "PCNA", "TOP2A", "CCND1", "MYC",

    # Tumor suppressors
    "BRCA2", "GATA3", "MAP3K1",

    # Oncogenes
    "PIK3CA", "AKT1", "MTOR",

    # Luminal cytokeratins
    "KRT8", "KRT18", "KRT19", "EPCAM",

    # Basal markers (low in Luminal)
    "KRT5", "KRT14", "EGFR",

    # Stromal markers
    "COL1A1", "COL3A1", "ACTA2", "FAP", "VIM",

    # Immune markers
    "CD3D", "CD8A", "CD4", "FOXP3", "CD68", "CD163",

    # Hypoxia
    "HIF1A", "VEGFA", "CA9",
]

# Breast tissue regions (7 types)
BREAST_REGIONS = [
    "tumor_core_luminal",      # ER+ tumor cells, low proliferation
    "tumor_proliferative",      # Ki67-high tumor zone
    "tumor_invasive_front",     # Tumor-stroma boundary
    "stroma_fibrous",           # Dense collagen, CAFs
    "stroma_immune",            # Lymphocyte infiltration (TILs)
    "adipose",                  # Normal breast adipose tissue
    "ductal_normal",            # Adjacent normal ductal epithelium
]


def generate_spatial_coordinates(grid_size: int) -> pd.DataFrame:
    """Generate Visium spatial coordinates.

    Args:
        grid_size: Grid dimension (e.g., 30 for 30x30 grid)

    Returns:
        DataFrame with spatial coordinates
    """
    spots = []
    for row in range(grid_size):
        for col in range(grid_size):
            barcode = f"SPOT_{row:02d}_{col:02d}"
            # Most spots are in tissue (95% coverage), edges may be out
            is_edge = (row == 0 or row == grid_size - 1 or
                      col == 0 or col == grid_size - 1)
            in_tissue = 0 if (is_edge and np.random.rand() < 0.3) else 1

            spots.append({
                "barcode": barcode,
                "in_tissue": in_tissue,
                "array_row": row,
                "array_col": col,
                "pxl_row_in_fullres": row * 1000,  # 1000 pixels per spot
                "pxl_col_in_fullres": col * 1000,
            })

    return pd.DataFrame(spots)


def assign_breast_tissue_regions(coords: pd.DataFrame) -> pd.DataFrame:
    """Assign breast tissue regions to spatial coordinates.

    Breast cancer spatial organization:
    - Center: tumor core (luminal, proliferative)
    - Periphery: tumor invasive front
    - Surrounding: stroma (fibrous, immune)
    - Corners: adipose tissue
    - Scattered: normal ductal epithelium

    Args:
        coords: DataFrame with spatial coordinates

    Returns:
        DataFrame with region annotations
    """
    annotations = []
    center_row, center_col = GRID_SIZE // 2, GRID_SIZE // 2

    for _, row in coords.iterrows():
        barcode = row["barcode"]
        r, c = row["array_row"], row["array_col"]

        # Calculate distance from center
        dist_from_center = np.sqrt((r - center_row)**2 + (c - center_col)**2)

        # Region assignment based on spatial location
        if dist_from_center < 6:
            # Central tumor core: luminal subtype (ER+/PR+)
            if np.random.rand() < 0.6:
                region = "tumor_core_luminal"
            else:
                # Proliferative hotspot (Ki67-high)
                region = "tumor_proliferative"

        elif dist_from_center < 10:
            # Tumor invasive front
            region = "tumor_invasive_front"

        elif dist_from_center < 14:
            # Stromal regions
            if np.random.rand() < 0.6:
                region = "stroma_fibrous"  # CAFs, collagen
            else:
                region = "stroma_immune"   # TILs

        else:
            # Periphery: adipose or normal ductal
            if (r < 5 or r > GRID_SIZE - 5) or (c < 5 or c > GRID_SIZE - 5):
                region = "adipose"  # Breast adipose tissue
            else:
                region = "ductal_normal"  # Adjacent normal

        annotations.append({
            "barcode": barcode,
            "region": region,
        })

    return pd.DataFrame(annotations)


def generate_gene_expression(
    coords: pd.DataFrame,
    annotations: pd.DataFrame,
    genes: list,
) -> pd.DataFrame:
    """Generate Visium gene expression counts by region.

    Args:
        coords: Spatial coordinates
        annotations: Region annotations
        genes: List of gene names

    Returns:
        DataFrame with gene expression (spots x genes)
    """
    n_spots = len(coords)
    n_genes = len(genes)

    # Initialize expression matrix
    expression = np.zeros((n_spots, n_genes))

    # Define baseline expression per gene (Visium UMI counts)
    baseline_expr = {
        # Luminal markers: High expression
        "ESR1": (150, 80), "PGR": (120, 70), "GATA3": (180, 90), "FOXA1": (100, 60),

        # Proliferation: Moderate (post-treatment, Ki67 ~20%)
        "MKI67": (80, 50), "PCNA": (60, 40), "TOP2A": (50, 30), "CCND1": (90, 50),
        "MYC": (70, 40),

        # Tumor suppressors
        "BRCA2": (60, 30),  # Reduced due to germline mutation
        "MAP3K1": (50, 25),

        # PI3K pathway (active)
        "PIK3CA": (100, 60), "AKT1": (80, 50), "MTOR": (70, 40),

        # Luminal cytokeratins: High
        "KRT8": (200, 100), "KRT18": (180, 90), "KRT19": (150, 80), "EPCAM": (140, 70),

        # Basal markers: Very low (Luminal type)
        "KRT5": (10, 5), "KRT14": (8, 4), "EGFR": (15, 8),

        # Stromal markers
        "COL1A1": (120, 70), "COL3A1": (100, 60), "ACTA2": (80, 50),
        "FAP": (70, 40), "VIM": (90, 50),

        # Immune markers
        "CD3D": (50, 30), "CD8A": (40, 25), "CD4": (45, 28), "FOXP3": (30, 20),
        "CD68": (60, 35), "CD163": (55, 32),

        # Hypoxia: Moderate
        "HIF1A": (80, 50), "VEGFA": (100, 60), "CA9": (60, 40),
    }

    # Generate expression for each spot
    for i, (_, coord_row) in enumerate(coords.iterrows()):
        barcode = coord_row["barcode"]
        region = annotations[annotations["barcode"] == barcode]["region"].values[0]

        for j, gene in enumerate(genes):
            base_mean, base_std = baseline_expr.get(gene, (50, 30))

            # Modulate expression by region
            if region == "tumor_core_luminal":
                # High ER/PR, low proliferation
                if gene in ["ESR1", "PGR", "GATA3", "FOXA1", "KRT8", "KRT18", "KRT19", "EPCAM"]:
                    multiplier = np.random.uniform(2.5, 4.0)
                elif gene in ["MKI67", "PCNA", "TOP2A"]:
                    multiplier = np.random.uniform(0.6, 1.0)  # Low proliferation
                elif gene == "BRCA2":
                    multiplier = 0.5  # Germline mutation
                else:
                    multiplier = np.random.uniform(0.8, 1.2)

            elif region == "tumor_proliferative":
                # Ki67-high zone
                if gene in ["MKI67", "PCNA", "TOP2A", "CCND1", "MYC"]:
                    multiplier = np.random.uniform(3.0, 5.0)
                elif gene in ["ESR1", "PGR"]:
                    multiplier = np.random.uniform(1.5, 2.5)
                else:
                    multiplier = np.random.uniform(0.8, 1.2)

            elif region == "tumor_invasive_front":
                # Mesenchymal transition, EMT markers
                if gene in ["VIM", "ACTA2", "FAP"]:
                    multiplier = np.random.uniform(1.5, 2.5)
                elif gene in ["CD8A", "CD3D", "CD4"]:
                    multiplier = np.random.uniform(1.3, 2.0)  # TILs at boundary
                else:
                    multiplier = np.random.uniform(0.7, 1.0)

            elif region == "stroma_fibrous":
                # High collagen, CAFs
                if gene in ["COL1A1", "COL3A1", "ACTA2", "FAP", "VIM"]:
                    multiplier = np.random.uniform(3.0, 5.0)
                elif gene in ["ESR1", "PGR", "KRT8", "KRT18"]:
                    multiplier = 0.1  # No epithelial markers
                else:
                    multiplier = np.random.uniform(0.2, 0.5)

            elif region == "stroma_immune":
                # TILs (tumor-infiltrating lymphocytes)
                if gene in ["CD3D", "CD8A", "CD4", "FOXP3", "CD68", "CD163"]:
                    multiplier = np.random.uniform(3.0, 5.0)
                elif gene in ["ESR1", "PGR", "KRT8"]:
                    multiplier = 0.1
                else:
                    multiplier = np.random.uniform(0.3, 0.6)

            elif region == "adipose":
                # Breast adipose tissue (low cellularity)
                multiplier = np.random.uniform(0.05, 0.15)

            elif region == "ductal_normal":
                # Normal ductal epithelium
                if gene in ["KRT8", "KRT18", "KRT19", "EPCAM"]:
                    multiplier = np.random.uniform(1.0, 1.5)
                elif gene in ["ESR1", "PGR"]:
                    multiplier = np.random.uniform(0.5, 1.0)  # Lower than tumor
                elif gene in ["MKI67", "PCNA"]:
                    multiplier = 0.2  # Minimal proliferation
                else:
                    multiplier = np.random.uniform(0.3, 0.7)

            else:
                multiplier = 1.0

            # Generate count with log-normal distribution
            count = np.random.lognormal(np.log(base_mean * multiplier), base_std / base_mean)
            count = max(0, int(count))  # Ensure non-negative integer

            expression[i, j] = count

    # Create DataFrame
    expr_df = pd.DataFrame(
        expression,
        columns=genes,
    )
    # Add barcode as first column (unnamed index)
    expr_df.insert(0, "", coords["barcode"].values)

    return expr_df


def filter_low_count_genes(expr_df: pd.DataFrame, min_counts: int = 10) -> pd.DataFrame:
    """Filter genes with low total counts across all spots.

    Args:
        expr_df: Gene expression DataFrame
        min_counts: Minimum total counts threshold

    Returns:
        Filtered DataFrame
    """
    # Skip first column (barcodes)
    gene_totals = expr_df.iloc[:, 1:].sum(axis=0)
    keep_genes = gene_totals[gene_totals >= min_counts].index.tolist()

    # Keep barcode column + filtered genes
    filtered_df = expr_df[[""] + keep_genes].copy()

    return filtered_df


def main():
    """Generate all spatial transcriptomics files for PAT002-BC-2026."""
    print("="*60)
    print("Generating Visium Spatial Transcriptomics Data")
    print("Patient: PAT002-BC-2026 (Breast Cancer)")
    print("="*60)

    output_dir = Path(__file__).parent

    # 1. Generate spatial coordinates
    print(f"\n1. Generating spatial coordinates: {N_SPOTS} spots ({GRID_SIZE}x{GRID_SIZE} grid)")
    coords = generate_spatial_coordinates(GRID_SIZE)
    coords_path = output_dir / "visium_spatial_coordinates.csv"
    coords.to_csv(coords_path, index=False)
    print(f"   Saved to {coords_path}")
    print(f"   In-tissue spots: {coords['in_tissue'].sum()} / {N_SPOTS}")

    # 2. Assign tissue regions
    print(f"\n2. Generating region annotations: {len(BREAST_REGIONS)} tissue types")
    annotations = assign_breast_tissue_regions(coords)
    annot_path = output_dir / "visium_region_annotations.csv"
    annotations.to_csv(annot_path, index=False)
    print(f"   Saved to {annot_path}")
    print("   Region distribution:")
    for region in BREAST_REGIONS:
        count = (annotations["region"] == region).sum()
        pct = 100 * count / N_SPOTS
        print(f"     - {region}: {count} spots ({pct:.1f}%)")

    # 3. Generate gene expression
    print(f"\n3. Generating gene expression: {len(BREAST_CANCER_GENES)} genes x {N_SPOTS} spots")
    expression = generate_gene_expression(coords, annotations, BREAST_CANCER_GENES)
    expr_path = output_dir / "visium_gene_expression.csv"
    expression.to_csv(expr_path, index=False)
    print(f"   Saved to {expr_path}")

    # Calculate statistics (skip first column which is barcodes)
    expr_values = expression.iloc[:, 1:].values.flatten()
    print(f"   Expression range: {expr_values.min():.0f} - {expr_values.max():.0f} UMI counts")
    print(f"   Mean counts per spot: {expr_values.mean():.0f}")
    print(f"   Median counts per spot: {np.median(expr_values):.0f}")

    # 4. Generate filtered version
    print(f"\n4. Generating QC-filtered expression (min 10 counts)")
    filtered_dir = output_dir / "filtered"
    filtered_dir.mkdir(exist_ok=True)

    filtered_expression = filter_low_count_genes(expression, min_counts=10)
    filtered_path = filtered_dir / "visium_gene_expression_filtered.csv"
    filtered_expression.to_csv(filtered_path, index=False)
    print(f"   Saved to {filtered_path}")
    print(f"   Genes retained: {filtered_expression.shape[1] - 1} / {len(BREAST_CANCER_GENES)}")

    # Summary
    print("\n" + "="*60)
    print("Spatial Transcriptomics Data Generation Complete")
    print("="*60)
    print(f"\nFiles generated:")
    print(f"  - visium_spatial_coordinates.csv ({N_SPOTS} spots)")
    print(f"  - visium_region_annotations.csv ({len(BREAST_REGIONS)} regions)")
    print(f"  - visium_gene_expression.csv ({len(BREAST_CANCER_GENES)} genes)")
    print(f"  - filtered/visium_gene_expression_filtered.csv (QC-filtered)")
    print(f"\nBreast cancer spatial features:")
    print(f"  - ER+/PR+ luminal tumor core (high hormone receptor expression)")
    print(f"  - Ki67-high proliferative zones (post-treatment, moderate)")
    print(f"  - Tumor-stroma interface with TILs")
    print(f"  - Fibrous stroma with CAFs (collagen-rich)")
    print(f"  - Immune infiltration zones (CD8+ T cells)")
    print(f"  - Adjacent normal breast tissue (adipose, ductal)")
    print(f"  - BRCA2 haploinsufficiency (50% expression)")
    print("\n✓ Generation successful!")


if __name__ == "__main__":
    main()
