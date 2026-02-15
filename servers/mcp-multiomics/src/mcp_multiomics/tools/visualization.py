"""Multi-omics visualization: heatmap and PCA."""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from .utils import load_integrated_data

logger = logging.getLogger(__name__)


def create_multiomics_heatmap_impl(
    data_path: str,
    features: Optional[List[str]] = None,
    cluster_rows: bool = True,
    cluster_cols: bool = True,
    output_path: Optional[str] = None,
) -> Dict[str, Any]:
    """Create clustered heatmap from integrated multi-omics data.

    Args:
        data_path: Path to integrated data pickle (from integrate_omics_data)
        features: Specific features to include (default: top 50 by variance)
        cluster_rows: Hierarchical clustering on rows (features)
        cluster_cols: Hierarchical clustering on columns (samples)
        output_path: Path to save PNG

    Returns:
        Dictionary with plot path, cluster info, and figure metadata.
    """
    logger.info(f"Loading integrated data from {data_path}")
    omics_data, metadata = load_integrated_data(data_path)

    # Concatenate all modalities into a single matrix (features x samples)
    dfs = []
    modality_labels = []
    for modality, df in omics_data.items():
        # Prefix features with modality to avoid name collisions
        prefixed = df.copy()
        prefixed.index = [f"{modality}:{f}" for f in prefixed.index]
        dfs.append(prefixed)
        modality_labels.extend([modality] * len(prefixed))
    combined = pd.concat(dfs, axis=0)

    # Fill NaN with 0 for clustering
    combined = combined.fillna(0)

    # Select features
    if features:
        # Match with or without modality prefix
        mask = combined.index.isin(features)
        if mask.sum() == 0:
            # Try matching without prefix
            base_names = [idx.split(":", 1)[-1] for idx in combined.index]
            mask = pd.Series(base_names, index=combined.index).isin(features)
        combined = combined[mask]
        if combined.empty:
            raise ValueError(f"None of the requested features found in data")
    else:
        # Top 50 features by variance across samples
        variances = combined.var(axis=1)
        top_n = min(50, len(variances))
        combined = combined.loc[variances.nlargest(top_n).index]

    n_features, n_samples = combined.shape
    logger.info(f"Heatmap: {n_features} features x {n_samples} samples")

    # Build cluster info
    cluster_info = {
        "row_clusters": 0,
        "col_clusters": 0,
        "row_linkage": "complete" if cluster_rows else "none",
        "col_linkage": "complete" if cluster_cols else "none",
    }

    # Create the heatmap
    g = sns.clustermap(
        combined,
        cmap="RdBu_r",
        center=0,
        row_cluster=cluster_rows,
        col_cluster=cluster_cols,
        method="complete",
        figsize=(max(8, n_samples * 0.6), max(6, n_features * 0.3)),
        xticklabels=True,
        yticklabels=True,
        linewidths=0.5,
        cbar_kws={"label": "Z-score"},
    )
    g.ax_heatmap.set_xlabel("Samples")
    g.ax_heatmap.set_ylabel("Features")
    g.figure.suptitle("Multi-omics Integrated Heatmap", y=1.02, fontsize=14)

    if cluster_rows:
        cluster_info["row_clusters"] = len(set(leaves_list(g.dendrogram_row.linkage)))
    if cluster_cols:
        cluster_info["col_clusters"] = len(set(leaves_list(g.dendrogram_col.linkage)))

    # Save
    if not output_path:
        output_path = str(Path(data_path).parent / "multiomics_heatmap.png")
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    g.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close("all")
    logger.info(f"Heatmap saved to {output_path}")

    return {
        "plot_path": output_path,
        "cluster_info": cluster_info,
        "figure_data": {
            "type": "heatmap",
            "rows": n_features,
            "cols": n_samples,
            "format": "png",
        },
        "status": "success",
    }


def run_multiomics_pca_impl(
    data_path: str,
    modalities: Optional[List[str]] = None,
    n_components: int = 3,
    scale_features: bool = True,
    output_path: Optional[str] = None,
) -> Dict[str, Any]:
    """Run PCA on integrated multi-omics data.

    Args:
        data_path: Path to integrated data pickle (from integrate_omics_data)
        modalities: Which modalities to include (default: all)
        n_components: Number of PCs (default: 3)
        scale_features: StandardScaler before PCA
        output_path: Path to save scatter plot PNG

    Returns:
        Dictionary with variance explained, loadings, sample coordinates, and plot path.
    """
    logger.info(f"Loading integrated data from {data_path}")
    omics_data, metadata = load_integrated_data(data_path)

    # Filter modalities
    if modalities:
        omics_data = {k: v for k, v in omics_data.items() if k in modalities}
        if not omics_data:
            raise ValueError(f"None of the requested modalities found: {modalities}")

    # Concatenate features across modalities (features x samples -> transpose for PCA)
    dfs = []
    for modality, df in omics_data.items():
        prefixed = df.copy()
        prefixed.index = [f"{modality}:{f}" for f in prefixed.index]
        dfs.append(prefixed)
    combined = pd.concat(dfs, axis=0).fillna(0)

    # PCA expects samples x features
    X = combined.T.values
    sample_names = list(combined.columns)
    feature_names = list(combined.index)

    n_components = min(n_components, X.shape[0], X.shape[1])

    if scale_features:
        X = StandardScaler().fit_transform(X)

    pca = PCA(n_components=n_components)
    coords = pca.fit_transform(X)

    # Build results
    variance_explained = [round(float(v), 4) for v in pca.explained_variance_ratio_]

    # Top loadings per component
    loadings = {}
    for i in range(n_components):
        pc_label = f"PC{i + 1}"
        abs_loadings = np.abs(pca.components_[i])
        top_idx = abs_loadings.argsort()[-5:][::-1]
        loadings[pc_label] = {
            "top_features": [feature_names[j] for j in top_idx],
            "n_features": len(feature_names),
        }

    # Sample coordinates
    sample_coords = {
        name: {f"PC{j + 1}": round(float(coords[i, j]), 4) for j in range(n_components)}
        for i, name in enumerate(sample_names)
    }

    # Create scatter plot (PC1 vs PC2)
    fig, ax = plt.subplots(figsize=(8, 6))

    # Color by metadata group if available
    colors = None
    if metadata is not None:
        for col in ["Response", "Treatment", "Group"]:
            if col in metadata.columns:
                aligned_meta = metadata.loc[metadata.index.intersection(sample_names)]
                groups = aligned_meta[col]
                unique_groups = groups.unique()
                cmap = plt.cm.Set1(np.linspace(0, 1, len(unique_groups)))
                group_colors = {g: cmap[i] for i, g in enumerate(unique_groups)}
                colors = [group_colors.get(groups.get(s, "unknown"), (0.5, 0.5, 0.5, 1)) for s in sample_names]
                # Add legend
                for g, c in group_colors.items():
                    ax.scatter([], [], color=c, label=g, s=60)
                ax.legend(title=col)
                break

    ax.scatter(coords[:, 0], coords[:, 1], c=colors, s=60, edgecolors="black", linewidths=0.5)
    for i, name in enumerate(sample_names):
        ax.annotate(name, (coords[i, 0], coords[i, 1]), fontsize=7, alpha=0.7,
                    xytext=(5, 5), textcoords="offset points")

    ax.set_xlabel(f"PC1 ({variance_explained[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({variance_explained[1] * 100:.1f}%)" if n_components > 1 else "PC2")
    ax.set_title("Multi-omics PCA")
    ax.grid(alpha=0.3)

    if not output_path:
        output_path = str(Path(data_path).parent / "multiomics_pca.png")
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"PCA plot saved to {output_path}")

    return {
        "variance_explained": variance_explained,
        "loadings": loadings,
        "sample_coordinates": sample_coords,
        "plot_path": output_path,
        "statistics": {
            "total_variance": round(float(sum(variance_explained)), 4),
            "n_features": len(feature_names),
            "modalities_used": list(omics_data.keys()),
        },
        "status": "success",
    }
