#!/usr/bin/env python3
"""
Generate synthetic perturbation data for PAT002-BC-2026.
T cell exhaustion CRISPR screen: PD-1 (PDCD1) knockout vs control.
"""

import numpy as np
import pandas as pd
import anndata as ad
from pathlib import Path

# Set random seed
np.random.seed(44)

# Configuration
N_CELLS = 500
N_GENES = 100
OUTPUT_DIR = Path(__file__).parent
OUTPUT_FILE = OUTPUT_DIR / "pat002_tcells.h5ad"

print("="*70)
print("Generating T Cell Perturbation Data for PAT002-BC-2026")
print("CRISPR Screen: PD-1 (PDCD1) Knockout vs Control")
print("="*70)

# ============================================================================
# 1. Define experimental design
# ============================================================================
print("\n1. Setting up experimental design...")

# Two conditions: Control (non-targeting guide) vs PDCD1 knockout
n_control = 250
n_pdcd1_ko = 250

conditions = ['control'] * n_control + ['PDCD1_KO'] * n_pdcd1_ko
perturbations = ['non-targeting'] * n_control + ['PDCD1_knockout'] * n_pdcd1_ko

print(f"   ✓ {n_control} control cells, {n_pdcd1_ko} PDCD1 KO cells")

# ============================================================================
# 2. Define T cell marker genes
# ============================================================================
print("\n2. Defining T cell gene panel...")

# T cell core markers
tcell_markers = ['CD8A', 'CD8B', 'CD3D', 'CD3E', 'CD3G', 'TCR', 'TRAC', 'TRBC1']

# Effector functions
effector_genes = [
    'GZMB', 'GZMA', 'GZMK', 'PRF1', 'GNLY', 'NKG7',  # Cytotoxicity
    'IFNG', 'TNF', 'IL2', 'FASLG', 'LTA'  # Cytokines
]

# Exhaustion markers (higher in control, lower in PDCD1 KO)
exhaustion_genes = [
    'PDCD1',  # PD-1 (knocked out in treatment group)
    'HAVCR2',  # TIM-3
    'LAG3',  # LAG-3
    'TIGIT',  # TIGIT
    'CTLA4',  # CTLA-4
    'TOX',  # TOX (exhaustion TF)
    'TOX2',
    'ENTPD1',  # CD39
    'NT5E'  # CD73
]

# Activation/memory markers
activation_genes = [
    'CD69', 'CD25', 'IL2RA', 'ICOS', 'CD28', 'CD27', 'CCR7', 'CD62L', 'SELL'
]

# Proliferation markers (increased in PDCD1 KO)
proliferation_genes = [
    'MKI67', 'PCNA', 'TOP2A', 'CCND1', 'CCNB1', 'AURKA'
]

# Transcription factors
tf_genes = [
    'EOMES', 'TBX21', 'PRDM1', 'TCF7', 'LEF1', 'ID2', 'ID3', 'BATF', 'IRF4'
]

# Metabolism genes
metabolism_genes = [
    'LDHA', 'PKM', 'HK2', 'PFKP', 'GAPDH', 'ENO1', 'SLC2A1'
]

# Migration/homing
migration_genes = [
    'ITGAE', 'ITGAL', 'SELP', 'CXCR3', 'CCR5', 'S1PR1'
]

# Survival/apoptosis
survival_genes = [
    'BCL2', 'BCL2L1', 'MCL1', 'BIM', 'BMF', 'BAX', 'CASP3', 'FAS'
]

# Combine all genes
all_genes = (
    tcell_markers + effector_genes + exhaustion_genes +
    activation_genes + proliferation_genes + tf_genes +
    metabolism_genes + migration_genes + survival_genes
)

# Fill to N_GENES with generic names if needed
if len(all_genes) < N_GENES:
    remaining = N_GENES - len(all_genes)
    all_genes += [f'TCELL_GENE_{i:03d}' for i in range(remaining)]

gene_names = all_genes[:N_GENES]

print(f"   ✓ {N_GENES} T cell genes defined")
print(f"   - Core markers: {len(tcell_markers)}")
print(f"   - Effector: {len(effector_genes)}")
print(f"   - Exhaustion: {len(exhaustion_genes)}")
print(f"   - Proliferation: {len(proliferation_genes)}")

# ============================================================================
# 3. Generate gene expression matrix
# ============================================================================
print("\n3. Generating gene expression matrix...")

# Initialize expression matrix (negative binomial for count data)
expression = np.zeros((N_CELLS, N_GENES))

for i in range(N_CELLS):
    condition = conditions[i]

    for j, gene in enumerate(gene_names):
        # Baseline expression for all T cells
        baseline = np.random.negative_binomial(8, 0.4)

        # Modulate expression based on gene and condition
        if gene in tcell_markers:
            # High expression in all T cells
            expression[i, j] = np.random.negative_binomial(50, 0.2)

        elif gene == 'PDCD1':
            # PDCD1 (PD-1): High in control, knocked out in PDCD1_KO
            if condition == 'control':
                expression[i, j] = np.random.negative_binomial(60, 0.15)
            else:
                # Near-zero in knockout (residual low expression)
                expression[i, j] = np.random.negative_binomial(2, 0.7)

        elif gene in exhaustion_genes and gene != 'PDCD1':
            # Other exhaustion markers: High in control, reduced in PDCD1 KO
            if condition == 'control':
                expression[i, j] = np.random.negative_binomial(40, 0.2)
            else:
                # Reduced exhaustion in PDCD1 KO
                expression[i, j] = np.random.negative_binomial(15, 0.35)

        elif gene in effector_genes:
            # Effector functions: Enhanced in PDCD1 KO (reinvigorated T cells)
            if condition == 'control':
                expression[i, j] = np.random.negative_binomial(20, 0.3)
            else:
                # Increased effector function in PDCD1 KO
                expression[i, j] = np.random.negative_binomial(45, 0.18)

        elif gene in proliferation_genes:
            # Proliferation: Increased in PDCD1 KO
            if condition == 'control':
                expression[i, j] = np.random.negative_binomial(10, 0.4)
            else:
                expression[i, j] = np.random.negative_binomial(30, 0.25)

        elif gene in activation_genes:
            # Activation markers: Increased in PDCD1 KO
            if condition == 'control':
                expression[i, j] = np.random.negative_binomial(15, 0.35)
            else:
                expression[i, j] = np.random.negative_binomial(35, 0.22)

        elif gene in metabolism_genes:
            # Glycolysis: Increased in activated PDCD1 KO cells
            if condition == 'control':
                expression[i, j] = np.random.negative_binomial(25, 0.28)
            else:
                expression[i, j] = np.random.negative_binomial(45, 0.2)

        elif gene == 'TOX' or gene == 'TOX2':
            # TOX: Exhaustion TF, reduced in PDCD1 KO
            if condition == 'control':
                expression[i, j] = np.random.negative_binomial(50, 0.18)
            else:
                expression[i, j] = np.random.negative_binomial(12, 0.4)

        else:
            # Generic genes: slight variation by condition
            if condition == 'control':
                expression[i, j] = baseline
            else:
                expression[i, j] = baseline * np.random.uniform(0.8, 1.3)

print(f"   ✓ Generated expression matrix: {N_CELLS} cells × {N_GENES} genes")

# ============================================================================
# 4. Create AnnData object
# ============================================================================
print("\n4. Creating AnnData object...")

# Create observation metadata
obs = pd.DataFrame({
    'cell_type': ['CD8_T_cells'] * N_CELLS,
    'condition': conditions,
    'perturbation': perturbations,
    'patient_id': ['PAT002'] * N_CELLS,
    'batch': np.random.choice(['Batch1', 'Batch2', 'Batch3'], N_CELLS),
    'n_genes': (expression > 0).sum(axis=1),
    'total_counts': expression.sum(axis=1)
})

# Create gene metadata
var = pd.DataFrame({
    'gene_name': gene_names,
    'gene_category': [
        'tcell_marker' if g in tcell_markers else
        'effector' if g in effector_genes else
        'exhaustion' if g in exhaustion_genes else
        'proliferation' if g in proliferation_genes else
        'activation' if g in activation_genes else
        'transcription_factor' if g in tf_genes else
        'metabolism' if g in metabolism_genes else
        'migration' if g in migration_genes else
        'survival' if g in survival_genes else
        'other'
        for g in gene_names
    ]
})
var.index = gene_names

# Create AnnData
adata = ad.AnnData(
    X=expression,
    obs=obs,
    var=var
)

# Add metadata
adata.uns['patient_id'] = 'PAT002'
adata.uns['experiment_type'] = 'CRISPR_screen'
adata.uns['perturbation_target'] = 'PDCD1 (PD-1)'
adata.uns['cell_type'] = 'CD8+ T cells (tumor-infiltrating)'
adata.uns['hypothesis'] = 'PD-1 blockade reverses T cell exhaustion and restores effector function'
adata.uns['disease_context'] = 'ER+/PR+/HER2- breast cancer, BRCA2+, post-adjuvant'
adata.uns['conditions'] = ['control', 'PDCD1_KO']

print(f"   ✓ Created AnnData object with metadata")

# ============================================================================
# 5. Save AnnData object
# ============================================================================
print(f"\n5. Saving to {OUTPUT_FILE}...")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(OUTPUT_FILE)

# ============================================================================
# 6. Print summary
# ============================================================================
print("\n" + "="*70)
print("✓ T Cell Perturbation Data Created Successfully!")
print("="*70)
print(f"\nFile: {OUTPUT_FILE}")
print(f"Size: {OUTPUT_FILE.stat().st_size / 1024:.2f} KB")
print(f"\nData dimensions:")
print(f"  - Cells: {adata.n_obs}")
print(f"  - Genes: {adata.n_vars}")
print(f"\nExperimental design:")
print(f"  - Control: {(adata.obs['condition'] == 'control').sum()} cells")
print(f"  - PDCD1 KO: {(adata.obs['condition'] == 'PDCD1_KO').sum()} cells")
print(f"\nGene categories:")
gene_category_counts = adata.var['gene_category'].value_counts()
for cat, count in gene_category_counts.items():
    print(f"  - {cat}: {count} genes")

# Calculate mean expression for key genes
print(f"\nKey gene expression (mean counts):")
key_genes_check = ['PDCD1', 'GZMB', 'IFNG', 'TOX', 'MKI67']
for gene in key_genes_check:
    if gene in gene_names:
        gene_idx = gene_names.index(gene)
        control_mean = expression[:n_control, gene_idx].mean()
        pdcd1_ko_mean = expression[n_control:, gene_idx].mean()
        fold_change = pdcd1_ko_mean / (control_mean + 1)
        print(f"  - {gene}: Control={control_mean:.1f}, PDCD1_KO={pdcd1_ko_mean:.1f} (FC={fold_change:.2f})")

print(f"\nExpected phenotype:")
print(f"  - Control: High PD-1, high exhaustion markers (TOX, LAG3)")
print(f"  - PDCD1 KO: No PD-1, restored effector function (GZMB, IFNG↑)")
print(f"  - PDCD1 KO: Increased proliferation (MKI67↑), reduced exhaustion")

print("\n" + "="*70)
print("Ready for perturbation analysis and machine learning!")
print("="*70)
