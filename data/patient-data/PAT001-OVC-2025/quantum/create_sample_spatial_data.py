#!/usr/bin/env python3
"""
Generate sample spatial transcriptomics data for quantum-celltype-fidelity testing.
Creates realistic ovarian cancer tumor microenvironment with spatial patterns.
"""

import numpy as np
import pandas as pd
import anndata as ad
from pathlib import Path

# Set random seed for reproducibility
np.random.seed(42)

# Configuration
N_CELLS = 500
N_GENES = 2000
OUTPUT_DIR = Path("/home/claude")
OUTPUT_FILE = OUTPUT_DIR / "PAT001_tumor_spatial.h5ad"

print("Generating sample spatial transcriptomics data...")

# ============================================================================
# 1. Generate spatial coordinates with realistic tumor architecture
# ============================================================================
print("  - Creating spatial coordinates...")

# Create distinct spatial regions
cell_types = []
spatial_x = []
spatial_y = []

# Tumor core (center, high density)
n_tumor_core = 150
x_tumor = np.random.normal(50, 10, n_tumor_core)
y_tumor = np.random.normal(50, 10, n_tumor_core)
cell_types.extend(['tumor_cells'] * n_tumor_core)
spatial_x.extend(x_tumor)
spatial_y.extend(y_tumor)

# CAFs (cancer-associated fibroblasts) - surrounding tumor
n_cafs = 100
x_cafs = np.random.normal(50, 20, n_cafs)
y_cafs = np.random.normal(50, 20, n_cafs)
cell_types.extend(['CAFs'] * n_cafs)
spatial_x.extend(x_cafs)
spatial_y.extend(y_cafs)

# CD8+ T cells (tumor margin, clustered)
n_cd8 = 80
x_cd8 = np.random.normal(70, 8, n_cd8)
y_cd8 = np.random.normal(30, 8, n_cd8)
cell_types.extend(['CD8_T_cells'] * n_cd8)
spatial_x.extend(x_cd8)
spatial_y.extend(y_cd8)

# Macrophages (scattered throughout)
n_macro = 70
x_macro = np.random.uniform(20, 80, n_macro)
y_macro = np.random.uniform(20, 80, n_macro)
cell_types.extend(['macrophages'] * n_macro)
spatial_x.extend(x_macro)
spatial_y.extend(y_macro)

# B cells (tertiary lymphoid structure candidate)
n_bcells = 50
x_bcells = np.random.normal(75, 5, n_bcells)
y_bcells = np.random.normal(75, 5, n_bcells)
cell_types.extend(['B_cells'] * n_bcells)
spatial_x.extend(x_bcells)
spatial_y.extend(y_bcells)

# Endothelial cells (vasculature)
n_endo = 30
x_endo = np.random.uniform(20, 80, n_endo)
y_endo = np.random.uniform(20, 80, n_endo)
cell_types.extend(['endothelial'] * n_endo)
spatial_x.extend(x_endo)
spatial_y.extend(y_endo)

# Exhausted T cells (near tumor, dysfunctional)
n_exhausted = 20
x_exhausted = np.random.normal(55, 12, n_exhausted)
y_exhausted = np.random.normal(55, 12, n_exhausted)
cell_types.extend(['exhausted_T_cells'] * n_exhausted)
spatial_x.extend(x_exhausted)
spatial_y.extend(y_exhausted)

# Ensure we have exactly N_CELLS
current_n = len(cell_types)
if current_n < N_CELLS:
    # Add more tumor cells to reach N_CELLS
    remaining = N_CELLS - current_n
    x_extra = np.random.normal(50, 15, remaining)
    y_extra = np.random.normal(50, 15, remaining)
    cell_types.extend(['tumor_cells'] * remaining)
    spatial_x.extend(x_extra)
    spatial_y.extend(y_extra)

spatial_x = np.array(spatial_x[:N_CELLS])
spatial_y = np.array(spatial_y[:N_CELLS])
cell_types = cell_types[:N_CELLS]

# ============================================================================
# 2. Generate gene expression matrix with realistic patterns
# ============================================================================
print("  - Generating gene expression matrix...")

# Define marker genes for each cell type
marker_genes = {
    'tumor_cells': ['TP53', 'MUC16', 'PAX8', 'WT1', 'CCND1', 'MYC', 'ERBB2'],
    'CAFs': ['COL1A1', 'VIM', 'ACTA2', 'FAP', 'PDGFRA', 'TGFB1'],
    'CD8_T_cells': ['CD8A', 'CD8B', 'GZMB', 'PRF1', 'CD3E', 'CD3D'],
    'macrophages': ['CD68', 'CD163', 'CSF1R', 'MARCO', 'CD14'],
    'B_cells': ['CD19', 'CD20', 'MS4A1', 'CD79A', 'IGHM', 'IGHG1'],
    'endothelial': ['PECAM1', 'VWF', 'CDH5', 'CD34', 'ENG'],
    'exhausted_T_cells': ['PDCD1', 'HAVCR2', 'LAG3', 'TIGIT', 'CTLA4', 'TOX']
}

# Proliferation markers
proliferation_genes = ['Ki67', 'MKI67', 'PCNA', 'TOP2A', 'CCNB1']

# Immune checkpoint genes
checkpoint_genes = ['PD1', 'PDL1', 'CD274', 'PDCD1LG2']

# Angiogenesis markers
angio_genes = ['VEGFA', 'VEGFR2', 'KDR', 'FLT1']

# All marker genes combined
all_markers = (
    list(set([g for genes in marker_genes.values() for g in genes])) +
    proliferation_genes + checkpoint_genes + angio_genes
)

# Add random background genes
background_genes = [f'GENE_{i:04d}' for i in range(N_GENES - len(all_markers))]
gene_names = all_markers + background_genes

# Initialize expression matrix
expression = np.random.negative_binomial(5, 0.3, size=(N_CELLS, N_GENES))

# Add cell-type-specific expression patterns
for i, cell_type in enumerate(cell_types):
    if cell_type in marker_genes:
        for marker in marker_genes[cell_type]:
            if marker in gene_names:
                gene_idx = gene_names.index(marker)
                # High expression of markers in corresponding cell type
                expression[i, gene_idx] = np.random.negative_binomial(50, 0.2)
    
    # Add spatial gradients (e.g., hypoxia in tumor core)
    if cell_type == 'tumor_cells':
        # Distance from center
        dist_from_center = np.sqrt((spatial_x[i] - 50)**2 + (spatial_y[i] - 50)**2)
        if dist_from_center < 10:  # Core hypoxia
            # Reduce proliferation markers in hypoxic core
            for marker in proliferation_genes:
                if marker in gene_names:
                    gene_idx = gene_names.index(marker)
                    expression[i, gene_idx] = expression[i, gene_idx] * 0.3

# Add proliferation to tumor cells
for i, cell_type in enumerate(cell_types):
    if cell_type == 'tumor_cells':
        for marker in proliferation_genes:
            if marker in gene_names:
                gene_idx = gene_names.index(marker)
                expression[i, gene_idx] = np.random.negative_binomial(30, 0.25)

# Add checkpoint expression to exhausted T cells
for i, cell_type in enumerate(cell_types):
    if cell_type == 'exhausted_T_cells':
        for marker in checkpoint_genes:
            if marker in gene_names:
                gene_idx = gene_names.index(marker)
                expression[i, gene_idx] = np.random.negative_binomial(40, 0.2)

# ============================================================================
# 3. Create AnnData object
# ============================================================================
print("  - Creating AnnData object...")

# Create observation metadata
obs = pd.DataFrame({
    'cell_type': cell_types,
    'spatial_x': spatial_x,
    'spatial_y': spatial_y,
    'sample_id': 'PAT001',
    'tissue_region': ['tumor_core' if ct in ['tumor_cells', 'CAFs'] else 'margin' 
                      for ct in cell_types],
    'n_genes': (expression > 0).sum(axis=1),
    'total_counts': expression.sum(axis=1)
})

# Create gene metadata
var = pd.DataFrame({
    'gene_name': gene_names,
    'highly_variable': [True if g in all_markers else False for g in gene_names]
})
var.index = gene_names

# Create AnnData
adata = ad.AnnData(
    X=expression,
    obs=obs,
    var=var
)

# Add spatial coordinates as obsm (standard location)
adata.obsm['spatial'] = np.column_stack([spatial_x, spatial_y])

# Add metadata
adata.uns['patient_id'] = 'PAT001'
adata.uns['diagnosis'] = 'High-grade serous ovarian carcinoma (HGSOC)'
adata.uns['sample_type'] = 'Primary tumor biopsy'
adata.uns['spatial_technology'] = '10x Visium'
adata.uns['ca125_level'] = 487
adata.uns['treatment_status'] = 'treatment_naive'

# ============================================================================
# 4. Save AnnData object
# ============================================================================
print(f"  - Saving to {OUTPUT_FILE}...")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(OUTPUT_FILE)

# ============================================================================
# 5. Print summary
# ============================================================================
print("\n" + "="*70)
print("✓ Sample spatial transcriptomics data created successfully!")
print("="*70)
print(f"\nFile: {OUTPUT_FILE}")
print(f"Size: {OUTPUT_FILE.stat().st_size / 1024 / 1024:.2f} MB")
print(f"\nData dimensions:")
print(f"  - Cells: {adata.n_obs}")
print(f"  - Genes: {adata.n_vars}")
print(f"\nCell type distribution:")
cell_type_counts = adata.obs['cell_type'].value_counts()
for ct, count in cell_type_counts.items():
    print(f"  - {ct}: {count} cells ({count/adata.n_obs*100:.1f}%)")
print(f"\nSpatial coordinates:")
print(f"  - X range: [{spatial_x.min():.1f}, {spatial_x.max():.1f}]")
print(f"  - Y range: [{spatial_y.min():.1f}, {spatial_y.max():.1f}]")
print(f"\nMarker genes included: {len(all_markers)}")
print(f"Example markers: {', '.join(all_markers[:10])}")
print(f"\nMetadata keys: {list(adata.uns.keys())}")
print("\n" + "="*70)
print("Ready for quantum embedding training!")
print("="*70)
