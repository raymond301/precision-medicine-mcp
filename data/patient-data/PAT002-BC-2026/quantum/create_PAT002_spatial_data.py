#!/usr/bin/env python3
"""
Generate sample spatial transcriptomics data for PAT002-BC-2026.
Creates realistic breast cancer tumor microenvironment with spatial patterns.
"""

import numpy as np
import pandas as pd
import anndata as ad
from pathlib import Path

# Set random seed for reproducibility (different from PAT001's 42)
np.random.seed(43)

# Configuration
N_CELLS = 500
N_GENES = 2000
OUTPUT_DIR = Path(__file__).parent
OUTPUT_FILE = OUTPUT_DIR / "PAT002_tumor_spatial.h5ad"

print("="*70)
print("Generating Spatial Transcriptomics Data for PAT002-BC-2026")
print("Patient: Michelle Thompson, ER+/PR+/HER2- Breast Cancer")
print("="*70)

# ============================================================================
# 1. Generate spatial coordinates with realistic breast tumor architecture
# ============================================================================
print("\n1. Creating spatial coordinates...")

# Create distinct spatial regions for breast cancer
cell_types = []
spatial_x = []
spatial_y = []

# Tumor core - luminal subtype (ER+/PR+ epithelial, center)
n_tumor = 180
x_tumor = np.random.normal(50, 12, n_tumor)
y_tumor = np.random.normal(50, 12, n_tumor)
cell_types.extend(['tumor_cells_luminal'] * n_tumor)
spatial_x.extend(x_tumor)
spatial_y.extend(y_tumor)

# CAFs (cancer-associated fibroblasts) - surrounding tumor
n_cafs = 100
x_cafs = np.random.normal(50, 22, n_cafs)
y_cafs = np.random.normal(50, 22, n_cafs)
cell_types.extend(['CAFs'] * n_cafs)
spatial_x.extend(x_cafs)
spatial_y.extend(y_cafs)

# CD8+ T cells (tumor-infiltrating lymphocytes, clustered at margins)
n_cd8 = 70
x_cd8 = np.random.normal(30, 10, n_cd8)
y_cd8 = np.random.normal(70, 10, n_cd8)
cell_types.extend(['CD8_T_cells'] * n_cd8)
spatial_x.extend(x_cd8)
spatial_y.extend(y_cd8)

# CD4+ T cells (helper T cells, scattered)
n_cd4 = 50
x_cd4 = np.random.normal(35, 12, n_cd4)
y_cd4 = np.random.normal(65, 12, n_cd4)
cell_types.extend(['CD4_T_cells'] * n_cd4)
spatial_x.extend(x_cd4)
spatial_y.extend(y_cd4)

# Macrophages (scattered throughout tumor and stroma)
n_macro = 40
x_macro = np.random.uniform(20, 80, n_macro)
y_macro = np.random.uniform(20, 80, n_macro)
cell_types.extend(['macrophages'] * n_macro)
spatial_x.extend(x_macro)
spatial_y.extend(y_macro)

# B cells (tertiary lymphoid structure)
n_bcells = 30
x_bcells = np.random.normal(75, 6, n_bcells)
y_bcells = np.random.normal(75, 6, n_bcells)
cell_types.extend(['B_cells'] * n_bcells)
spatial_x.extend(x_bcells)
spatial_y.extend(y_bcells)

# Endothelial cells (vasculature)
n_endo = 20
x_endo = np.random.uniform(20, 80, n_endo)
y_endo = np.random.uniform(20, 80, n_endo)
cell_types.extend(['endothelial'] * n_endo)
spatial_x.extend(x_endo)
spatial_y.extend(y_endo)

# Adipocytes (breast-specific, periphery)
n_adipo = 10
x_adipo_1 = np.random.normal(10, 5, n_adipo // 2)
y_adipo_1 = np.random.normal(10, 5, n_adipo // 2)
x_adipo_2 = np.random.normal(90, 5, n_adipo // 2)
y_adipo_2 = np.random.normal(90, 5, n_adipo // 2)
cell_types.extend(['adipocytes'] * n_adipo)
spatial_x.extend(list(x_adipo_1) + list(x_adipo_2))
spatial_y.extend(list(y_adipo_1) + list(y_adipo_2))

# Ensure we have exactly N_CELLS
current_n = len(cell_types)
if current_n < N_CELLS:
    # Add more tumor cells to reach N_CELLS
    remaining = N_CELLS - current_n
    x_extra = np.random.normal(50, 15, remaining)
    y_extra = np.random.normal(50, 15, remaining)
    cell_types.extend(['tumor_cells_luminal'] * remaining)
    spatial_x.extend(x_extra)
    spatial_y.extend(y_extra)

spatial_x = np.array(spatial_x[:N_CELLS])
spatial_y = np.array(spatial_y[:N_CELLS])
cell_types = cell_types[:N_CELLS]

print(f"   ✓ Generated {N_CELLS} cells across 8 cell types")

# ============================================================================
# 2. Generate gene expression matrix with realistic patterns
# ============================================================================
print("\n2. Generating gene expression matrix...")

# Define marker genes for each cell type (breast cancer-specific)
marker_genes = {
    'tumor_cells_luminal': [
        'ESR1', 'PGR', 'GATA3', 'FOXA1', 'XBP1', 'TFF1',  # Hormone receptors & luminal
        'KRT8', 'KRT18', 'KRT19', 'EPCAM', 'MUC1',  # Epithelial/luminal cytokeratins
        'CCND1', 'MYC', 'PIK3CA', 'AKT1', 'MTOR',  # Oncogenes
        'BRCA2'  # Tumor suppressor (reduced expression)
    ],
    'CAFs': ['COL1A1', 'COL3A1', 'VIM', 'ACTA2', 'FAP', 'PDGFRA', 'TGFB1'],
    'CD8_T_cells': ['CD8A', 'CD8B', 'GZMB', 'PRF1', 'CD3E', 'CD3D', 'IFNG'],
    'CD4_T_cells': ['CD4', 'CD3E', 'CD3D', 'IL2', 'IL4', 'FOXP3'],
    'macrophages': ['CD68', 'CD163', 'CSF1R', 'MARCO', 'CD14', 'CD206'],
    'B_cells': ['CD19', 'CD20', 'MS4A1', 'CD79A', 'IGHM', 'IGHG1'],
    'endothelial': ['PECAM1', 'VWF', 'CDH5', 'CD34', 'ENG', 'VEGFR2'],
    'adipocytes': ['ADIPOQ', 'FABP4', 'PLIN1', 'LEP', 'PPARG', 'CEBPA']  # Breast-specific
}

# Proliferation markers (reduced post-treatment)
proliferation_genes = ['MKI67', 'PCNA', 'TOP2A', 'CCNA2', 'CCNB1', 'AURKA']

# Immune checkpoint genes
checkpoint_genes = ['PDCD1', 'HAVCR2', 'LAG3', 'TIGIT', 'CTLA4', 'CD274']

# Tamoxifen response genes
tamoxifen_genes = ['CYP2D6', 'SULT1A1', 'UGT2B15']

# DNA repair (HRD signature)
dna_repair_genes = ['RAD51', 'PALB2', 'ATM', 'ATR', 'CHEK2']

# All marker genes combined
all_markers = (
    list(set([g for genes in marker_genes.values() for g in genes])) +
    proliferation_genes + checkpoint_genes + tamoxifen_genes + dna_repair_genes
)

# Add random background genes
background_genes = [f'GENE_{i:04d}' for i in range(N_GENES - len(all_markers))]
gene_names = all_markers + background_genes

# Initialize expression matrix (negative binomial for count data)
expression = np.random.negative_binomial(5, 0.3, size=(N_CELLS, N_GENES))

# Add cell-type-specific expression patterns
for i, cell_type in enumerate(cell_types):
    if cell_type in marker_genes:
        for marker in marker_genes[cell_type]:
            if marker in gene_names:
                gene_idx = gene_names.index(marker)
                # High expression of markers in corresponding cell type
                expression[i, gene_idx] = np.random.negative_binomial(50, 0.2)

    # Special handling for BRCA2 (germline mutation, reduced expression)
    if 'BRCA2' in gene_names:
        brca2_idx = gene_names.index('BRCA2')
        if cell_type == 'tumor_cells_luminal':
            # 50% reduced expression due to heterozygous germline mutation
            expression[i, brca2_idx] = np.random.negative_binomial(25, 0.3)

# Add proliferation to tumor cells (moderate, post-treatment)
for i, cell_type in enumerate(cell_types):
    if cell_type == 'tumor_cells_luminal':
        for marker in proliferation_genes:
            if marker in gene_names:
                gene_idx = gene_names.index(marker)
                # Moderate proliferation (post-adjuvant, Ki67 ~20%)
                expression[i, gene_idx] = np.random.negative_binomial(15, 0.3)

# Add high ER/PR expression to luminal tumor cells
for i, cell_type in enumerate(cell_types):
    if cell_type == 'tumor_cells_luminal':
        for marker in ['ESR1', 'PGR', 'GATA3', 'FOXA1']:
            if marker in gene_names:
                gene_idx = gene_names.index(marker)
                # Very high expression (ER+ 85%, PR+ 70%)
                expression[i, gene_idx] = np.random.negative_binomial(70, 0.15)

# Add tamoxifen response genes to tumor cells
for i, cell_type in enumerate(cell_types):
    if cell_type == 'tumor_cells_luminal':
        for marker in tamoxifen_genes:
            if marker in gene_names:
                gene_idx = gene_names.index(marker)
                expression[i, gene_idx] = np.random.negative_binomial(30, 0.25)

print(f"   ✓ Generated expression matrix: {N_CELLS} cells × {N_GENES} genes")

# ============================================================================
# 3. Create AnnData object
# ============================================================================
print("\n3. Creating AnnData object...")

# Create observation metadata
obs = pd.DataFrame({
    'cell_type': cell_types,
    'spatial_x': spatial_x,
    'spatial_y': spatial_y,
    'sample_id': 'PAT002',
    'tissue_region': ['tumor' if ct in ['tumor_cells_luminal', 'CAFs'] else
                      'immune' if ct in ['CD8_T_cells', 'CD4_T_cells', 'B_cells', 'macrophages'] else
                      'stroma' for ct in cell_types],
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
adata.uns['patient_id'] = 'PAT002'
adata.uns['diagnosis'] = 'Stage IIA (T2N0M0) ER+/PR+/HER2- Invasive Ductal Carcinoma'
adata.uns['sample_type'] = 'Post-adjuvant tumor bed biopsy'
adata.uns['spatial_technology'] = '10x Visium'
adata.uns['receptor_status'] = 'ER+ (85%), PR+ (70%), HER2-'
adata.uns['brca2_status'] = 'Germline pathogenic variant (c.5946delT)'
adata.uns['treatment_status'] = 'post_adjuvant'
adata.uns['current_therapy'] = 'Tamoxifen 20mg daily'
adata.uns['disease_status'] = 'Disease-free (NED)'

print(f"   ✓ Created AnnData object with metadata")

# ============================================================================
# 4. Save AnnData object
# ============================================================================
print(f"\n4. Saving to {OUTPUT_FILE}...")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(OUTPUT_FILE)

# ============================================================================
# 5. Print summary
# ============================================================================
print("\n" + "="*70)
print("✓ Spatial transcriptomics data created successfully!")
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
print(f"Example markers: {', '.join(all_markers[:15])}")
print(f"\nMetadata keys: {list(adata.uns.keys())}")
print(f"\nBreast cancer features:")
print(f"  - ER+/PR+ luminal tumor cells (high ESR1, PGR)")
print(f"  - BRCA2 haploinsufficiency (50% expression)")
print(f"  - Moderate proliferation (post-treatment, Ki67 ~20%)")
print(f"  - Tumor-infiltrating lymphocytes (CD8+, CD4+)")
print(f"  - CAFs and fibrous stroma (COL1A1, ACTA2)")
print(f"  - Breast adipocytes (ADIPOQ, FABP4)")
print("\n" + "="*70)
print("Ready for quantum embedding training!")
print("="*70)
