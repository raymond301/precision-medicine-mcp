"""Generate synthetic multi-omics data for PAT002-BC-2026 breast cancer patient."""

import numpy as np
import pandas as pd
from pathlib import Path

# Set random seed for reproducibility
np.random.seed(43)  # Different from PAT001's 42

# Sample information for pre/post treatment comparison
N_SAMPLES = 12
N_PRE_TREATMENT = 6
N_POST_TREATMENT = 6

# Feature counts
N_GENES = 1000
N_PROTEINS = 400
N_PHOSPHO = 250

# Sample names
sample_names = (
    [f"Pre_Sample_{i:02d}" for i in range(1, N_PRE_TREATMENT + 1)] +
    [f"Post_Sample_{i:02d}" for i in range(1, N_POST_TREATMENT + 1)]
)

# Treatment timepoints
timepoints = (
    ["pre_treatment"] * N_PRE_TREATMENT +
    ["post_treatment"] * N_POST_TREATMENT
)

# Breast cancer-specific gene panel
breast_cancer_genes = [
    # Hormone receptors and luminal markers
    "ESR1", "PGR", "GATA3", "FOXA1", "XBP1", "TFF1", "PDZK1",

    # HER2 pathway
    "ERBB2", "EGFR", "GRB7",

    # Proliferation markers
    "MKI67", "PCNA", "TOP2A", "CCND1", "MYC", "CCNE1", "CDK4", "CDK6",

    # Tumor suppressors
    "BRCA2", "BRCA1", "TP53", "RB1", "PTEN", "MAP3K1", "CDH1",

    # PI3K/AKT/mTOR pathway
    "PIK3CA", "AKT1", "AKT2", "MTOR", "PTEN", "PIK3R1",

    # Luminal cytokeratins
    "KRT8", "KRT18", "KRT19", "EPCAM", "MUC1",

    # Basal markers (low in Luminal type)
    "KRT5", "KRT14", "KRT17", "EGFR",

    # Stromal markers
    "COL1A1", "COL3A1", "ACTA2", "FAP", "VIM", "FN1", "PDGFRA",

    # Immune markers
    "CD3D", "CD3E", "CD8A", "CD8B", "CD4", "FOXP3", "CD68", "CD163",
    "PTPRC", "PDCD1", "CTLA4", "LAG3", "HAVCR2",

    # Hypoxia and angiogenesis
    "HIF1A", "VEGFA", "CA9", "LDHA",

    # Metastasis and invasion
    "SNAI1", "TWIST1", "CDH2", "MMP2", "MMP9",

    # Tamoxifen response
    "CYP2D6", "SULT1A1", "UGT2B15",

    # DNA repair (HRD signature genes)
    "RAD51", "BRCA2", "PALB2", "ATM", "ATR", "CHEK1", "CHEK2",

    # Cell cycle
    "CCNA2", "CCNB1", "AURKA", "AURKB", "PLK1", "BUB1",
]

# Additional generic genes
generic_genes = [f"Gene_{i:04d}" for i in range(1, N_GENES - len(breast_cancer_genes) + 1)]
gene_names = breast_cancer_genes + generic_genes

# Protein names (subset of genes)
protein_names = gene_names[:N_PROTEINS]

# Phosphosite names (key signaling pathways)
phospho_names = [
    # ER pathway
    "ESR1_S118", "ESR1_S167", "PGR_S294",

    # PI3K/AKT/mTOR
    "AKT1_S473", "AKT1_T308", "MTOR_S2448", "MTOR_S2481",
    "PIK3CA_S361", "GSK3B_S9", "TSC2_T1462",

    # MAPK pathway
    "ERK1_T202", "ERK2_Y204", "MEK1_S217", "MEK1_S221",

    # Cell cycle
    "RB1_S807", "RB1_S811", "CCND1_T286", "CDK4_T172",

    # DNA damage response
    "ATM_S1981", "CHEK2_T68", "TP53_S15", "TP53_S20",
    "BRCA1_S1524", "RAD51_Y315",
]

# Add generic phosphosites
remaining_phospho = N_PHOSPHO - len(phospho_names)
for i in range(remaining_phospho):
    gene = np.random.choice(gene_names[:200])
    site = f"S{np.random.randint(1, 500)}"
    phospho_names.append(f"{gene}_{site}")


def generate_expression_data(
    n_features: int,
    feature_names: list,
    n_pre: int,
    n_post: int,
    n_differential: int = 50,
    baseline_mean: float = 5.0,
    baseline_std: float = 1.5,
    fold_change_range: tuple = (1.5, 3.0),
    post_treatment_direction: str = "down",  # "down" = decreased in post-treatment
) -> pd.DataFrame:
    """Generate log-normal expression data with differential features for pre/post treatment.

    Args:
        n_features: Number of features (genes/proteins)
        feature_names: List of feature names
        n_pre: Number of pre-treatment samples
        n_post: Number of post-treatment samples
        n_differential: Number of differentially expressed features
        baseline_mean: Mean of log-normal distribution
        baseline_std: Std of log-normal distribution
        fold_change_range: Range of fold changes for differential features
        post_treatment_direction: "down" for decreased post-treatment (proliferation),
                                  "up" for increased post-treatment (immune response)

    Returns:
        DataFrame with features as rows, samples as columns
    """
    n_samples = n_pre + n_post

    # Initialize data matrix
    data = np.zeros((n_features, n_samples))

    # Generate baseline expression (log-normal)
    for i in range(n_features):
        baseline = np.random.lognormal(baseline_mean, baseline_std, n_samples)
        data[i, :] = baseline

    # Add differential expression for subset of features
    differential_indices = np.random.choice(n_features, n_differential, replace=False)

    for idx in differential_indices:
        # Random fold change
        fold_change = np.random.uniform(*fold_change_range)

        if post_treatment_direction == "down":
            # Downregulated in post-treatment (e.g., proliferation markers)
            data[idx, n_pre:] /= fold_change
        else:
            # Upregulated in post-treatment (e.g., immune markers)
            data[idx, n_pre:] *= fold_change

    # Modulate specific breast cancer markers
    for i, gene in enumerate(feature_names[:n_features]):
        if gene in ["MKI67", "PCNA", "TOP2A", "CCNA2", "CCNB1", "AURKA", "PLK1"]:
            # Proliferation markers: high pre-treatment, low post-treatment
            data[i, :n_pre] *= np.random.uniform(2.0, 3.5)
            data[i, n_pre:] /= np.random.uniform(1.5, 2.5)

        elif gene in ["CD8A", "CD8B", "CD3D", "GZMB", "PRF1", "IFNG"]:
            # Immune markers: increased post-treatment
            data[i, n_pre:] *= np.random.uniform(1.3, 2.0)

        elif gene in ["ESR1", "PGR", "GATA3", "FOXA1"]:
            # ER+ markers: high and stable (minor decrease due to tamoxifen)
            data[i, :] *= np.random.uniform(3.0, 5.0)
            data[i, n_pre:] *= np.random.uniform(0.8, 1.0)  # Slight decrease

        elif gene == "BRCA2":
            # BRCA2: reduced expression due to germline mutation
            data[i, :] *= 0.5  # Heterozygous loss

        elif gene in ["PIK3CA", "AKT1", "MTOR"]:
            # PI3K pathway: active in ER+ breast cancer
            data[i, :] *= np.random.uniform(1.5, 2.5)

        elif gene in ["KRT5", "KRT14", "KRT17"]:
            # Basal markers: very low in Luminal type
            data[i, :] *= 0.2

    # Add small amount of noise
    noise = np.random.normal(0, 0.1, data.shape)
    data = data * (1 + noise)

    # Ensure positive values
    data = np.maximum(data, 0.1)

    # Create DataFrame
    df = pd.DataFrame(
        data,
        index=feature_names[:n_features],
        columns=sample_names,
    )

    return df


def generate_metadata(sample_names: list, timepoints: list) -> pd.DataFrame:
    """Generate sample metadata.

    Args:
        sample_names: List of sample names
        timepoints: List of timepoint labels (pre_treatment / post_treatment)

    Returns:
        Metadata DataFrame
    """
    metadata = pd.DataFrame({
        "sample_id": sample_names,
        "timepoint": timepoints,
        "treatment_response": ["responsive"] * len(sample_names),  # All responsive (disease-free)
        "batch": [1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 3, 3],
        "patient_id": ["PAT002"] * len(sample_names),
        "tissue_type": ["primary_tumor"] * len(sample_names),
    })

    return metadata


def main():
    """Generate all PAT002 multi-omics fixtures."""
    print("Generating synthetic multi-omics data for PAT002-BC-2026...")

    # Get output directory
    output_dir = Path(__file__).parent

    # Generate RNA-seq expression data (FPKM/TPM values)
    print(f"\n1. Generating RNA-seq data: {N_GENES} genes x {N_SAMPLES} samples")
    rna_data = generate_expression_data(
        n_features=N_GENES,
        feature_names=gene_names,
        n_pre=N_PRE_TREATMENT,
        n_post=N_POST_TREATMENT,
        n_differential=120,  # 12% differential
        baseline_mean=6.5,
        baseline_std=2.0,
        post_treatment_direction="down",  # Overall decreased proliferation
    )
    rna_path = output_dir / "tumor_rna_seq.csv"
    rna_data.to_csv(rna_path)
    print(f"   Saved to {rna_path}")
    print(f"   Range: {rna_data.min().min():.2f} - {rna_data.max().max():.2f}")

    # Generate protein abundance data (TMT/iTRAQ values)
    print(f"\n2. Generating Proteomics data: {N_PROTEINS} proteins x {N_SAMPLES} samples")
    protein_data = generate_expression_data(
        n_features=N_PROTEINS,
        feature_names=protein_names,
        n_pre=N_PRE_TREATMENT,
        n_post=N_POST_TREATMENT,
        n_differential=50,  # 12% differential
        baseline_mean=5.0,
        baseline_std=1.5,
        post_treatment_direction="down",
    )
    protein_path = output_dir / "tumor_proteomics.csv"
    protein_data.to_csv(protein_path)
    print(f"   Saved to {protein_path}")
    print(f"   Range: {protein_data.min().min():.2f} - {protein_data.max().max():.2f}")

    # Generate phosphoproteomics data
    print(f"\n3. Generating Phosphoproteomics data: {N_PHOSPHO} sites x {N_SAMPLES} samples")
    phospho_data = generate_expression_data(
        n_features=N_PHOSPHO,
        feature_names=phospho_names,
        n_pre=N_PRE_TREATMENT,
        n_post=N_POST_TREATMENT,
        n_differential=35,  # 14% differential
        baseline_mean=4.0,
        baseline_std=1.8,
        fold_change_range=(2.0, 4.0),  # More dramatic changes in phosphorylation
        post_treatment_direction="down",  # Decreased PI3K/AKT signaling
    )
    phospho_path = output_dir / "tumor_phosphoproteomics.csv"
    phospho_data.to_csv(phospho_path)
    print(f"   Saved to {phospho_path}")
    print(f"   Range: {phospho_data.min().min():.2f} - {phospho_data.max().max():.2f}")

    # Generate sample metadata
    print(f"\n4. Generating sample metadata: {N_SAMPLES} samples")
    metadata = generate_metadata(sample_names, timepoints)
    metadata_path = output_dir / "sample_metadata.csv"
    metadata.to_csv(metadata_path, index=False)
    print(f"   Saved to {metadata_path}")

    # Print summary
    print("\n" + "="*60)
    print("PAT002-BC-2026 Multi-Omics Data Generation Complete")
    print("="*60)
    print(f"\nPatient: Michelle Thompson, 42F")
    print(f"Diagnosis: Stage IIA ER+/PR+/HER2- Breast Cancer, BRCA2+")
    print(f"Treatment: Post-adjuvant (AC-T chemo + radiation + tamoxifen)")
    print(f"\nData files generated:")
    print(f"  - tumor_rna_seq.csv: {N_GENES} genes x {N_SAMPLES} samples")
    print(f"  - tumor_proteomics.csv: {N_PROTEINS} proteins x {N_SAMPLES} samples")
    print(f"  - tumor_phosphoproteomics.csv: {N_PHOSPHO} phosphosites x {N_SAMPLES} samples")
    print(f"  - sample_metadata.csv: {N_SAMPLES} samples metadata")
    print(f"\nSample design:")
    print(f"  - Pre-treatment: {N_PRE_TREATMENT} samples (tumor biopsies)")
    print(f"  - Post-treatment: {N_POST_TREATMENT} samples (surveillance biopsies)")
    print(f"  - All samples: Responsive (disease-free status)")
    print("\nBiological features:")
    print(f"  - High ER/PR expression (Luminal subtype)")
    print(f"  - Low proliferation post-treatment (Ki67 ↓)")
    print(f"  - BRCA2 haploinsufficiency (50% expression)")
    print(f"  - PI3K/AKT pathway activation (PIK3CA mutation)")
    print(f"  - Decreased phospho-signaling post-treatment")
    print(f"  - Increased immune markers post-treatment (TILs)")
    print("\n✓ Generation successful!")


if __name__ == "__main__":
    main()
