"""MCP DeepCell server - Deep learning cell segmentation."""

import json
import logging
import os
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional

import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for server use
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from fastmcp import FastMCP
from google.cloud import storage
from PIL import Image
from skimage import color, measure

from .deepcell_engine import DeepCellEngine
from .intensity_classifier import IntensityClassifier

# Configure logging
logger = logging.getLogger(__name__)

mcp = FastMCP("deepcell")

def _is_dry_run() -> bool:
    """Check if DRY_RUN mode is enabled."""
    return os.getenv("DEEPCELL_DRY_RUN", "false").lower() == "true"

DRY_RUN = _is_dry_run()

# DRY_RUN warning wrapper
def add_dry_run_warning(result: Any) -> Any:
    """Add warning banner to results when in DRY_RUN mode."""
    if not DRY_RUN:
        return result

    warning = """
╔═══════════════════════════════════════════════════════════════════════════╗
║                    ⚠️  SYNTHETIC DATA WARNING ⚠️                          ║
╚═══════════════════════════════════════════════════════════════════════════╝

This result was generated in DRY_RUN mode and does NOT represent real analysis.

🔴 CRITICAL: Do NOT use this data for research decisions or publications.
🔴 All values are SYNTHETIC/MOCKED and have no scientific validity.

To enable real data processing, set: DEEPCELL_DRY_RUN=false

═══════════════════════════════════════════════════════════════════════════

"""

    if isinstance(result, dict):
        result["_DRY_RUN_WARNING"] = "SYNTHETIC DATA - NOT FOR RESEARCH USE"
        result["_message"] = warning.strip()
    elif isinstance(result, str):
        result = warning + result

    return result


DRY_RUN = os.getenv("DEEPCELL_DRY_RUN", "true").lower() == "true"
OUTPUT_DIR = Path(os.getenv("DEEPCELL_OUTPUT_DIR", "/workspace/output"))
MODEL_CACHE_DIR = Path(os.getenv("DEEPCELL_MODEL_CACHE_DIR", Path.home() / ".deepcell" / "models"))
USE_GPU = os.getenv("DEEPCELL_USE_GPU", "true").lower() == "true"

def _ensure_output_dir() -> None:
    """Ensure output directory exists."""
    try:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        (OUTPUT_DIR / "visualizations").mkdir(exist_ok=True)
        (OUTPUT_DIR / "segmentations").mkdir(exist_ok=True)
        (OUTPUT_DIR / "classifications").mkdir(exist_ok=True)
    except (OSError, PermissionError):
        pass

_ensure_output_dir()

def _download_from_gcs(gcs_path: str) -> str:
    """Download file from GCS to temporary location.

    Args:
        gcs_path: GCS path in format gs://bucket/path/to/file

    Returns:
        Local temporary file path

    Raises:
        ValueError: If path is not a valid GCS path
        Exception: If download fails
    """
    if not gcs_path.startswith("gs://"):
        # Return as-is if not a GCS path
        return gcs_path

    try:
        # Parse GCS path: gs://bucket/path/to/file
        path_parts = gcs_path[5:].split("/", 1)
        if len(path_parts) != 2:
            raise ValueError(f"Invalid GCS path format: {gcs_path}")

        bucket_name, blob_name = path_parts

        logger.info(f"Downloading from GCS: {gcs_path}")

        # Initialize GCS client (uses default credentials)
        storage_client = storage.Client()
        bucket = storage_client.bucket(bucket_name)
        blob = bucket.blob(blob_name)

        # Create temp file with same extension
        suffix = Path(blob_name).suffix
        temp_file = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
        temp_path = temp_file.name
        temp_file.close()

        # Download file
        blob.download_to_filename(temp_path)
        logger.info(f"Downloaded to: {temp_path}")

        return temp_path

    except Exception as e:
        logger.error(f"Failed to download from GCS: {e}")
        raise Exception(
            f"The cell segmentation tool was unable to find the image at the provided path: {gcs_path}.\n\n"
            f"Please verify that the path is correct and accessible. "
            f"If you are running this in an environment that requires local file paths, "
            f"you may need to download the file first or provide a local path to the image."
        )

# Initialize DeepCell engine and classifier (lazy loading)
_deepcell_engine: Optional[DeepCellEngine] = None
_intensity_classifier: Optional[IntensityClassifier] = None

def _get_deepcell_engine() -> DeepCellEngine:
    """Get or create DeepCell engine instance (singleton pattern)."""
    global _deepcell_engine
    if _deepcell_engine is None:
        logger.info("Initializing DeepCell engine...")
        _deepcell_engine = DeepCellEngine(model_cache_dir=MODEL_CACHE_DIR, use_gpu=USE_GPU)
    return _deepcell_engine

def _get_intensity_classifier() -> IntensityClassifier:
    """Get or create IntensityClassifier instance (singleton pattern)."""
    global _intensity_classifier
    if _intensity_classifier is None:
        logger.info("Initializing IntensityClassifier...")
        _intensity_classifier = IntensityClassifier()
    return _intensity_classifier

@mcp.tool()
async def segment_cells(
    image_path: str,
    model_type: str = "nuclear",
    min_cell_size: int = 100,
    image_mpp: Optional[float] = None
) -> Dict[str, Any]:
    """Deep learning-based cell segmentation using DeepCell models.

    Args:
        image_path: Path to microscopy image (TIFF, PNG, or other formats)
        model_type: Model to use - "nuclear" or "membrane" (default: "nuclear")
        min_cell_size: Minimum cell size in pixels (default: 100)
        image_mpp: Microns per pixel for physical size calculation (optional)

    Returns:
        Dictionary with:
        - segmentation_mask: Path to saved segmentation mask (TIFF)
        - cells_detected: Number of cells found
        - cells_filtered: Number of small cells filtered out
        - mean_cell_area: Mean cell area in pixels
        - median_cell_area: Median cell area in pixels
        - processing_time_seconds: Time taken for segmentation
        - model_used: Which model was used
        - gpu_used: Whether GPU acceleration was used

    Example:
        >>> result = await segment_cells(
        ...     image_path="/data/DAPI_nuclear.tiff",
        ...     model_type="nuclear",
        ...     min_cell_size=100
        ... )
        >>> print(f"Found {result['cells_detected']} cells")
        >>> # Use result['segmentation_mask'] for downstream analysis
    """
    if DRY_RUN:
        return add_dry_run_warning({
            "segmentation_mask": f"{image_path}.seg.tif",
            "cells_detected": 1247,
            "mean_cell_area": 850,
            "model_used": model_type,
            "processing_time_seconds": 15.3,
            "quality_metrics": {
                "mean_confidence": 0.89,
                "cells_filtered": 23
            },
            "mode": "dry_run"
        })

    try:
        # Download from GCS if needed
        local_image_path = _download_from_gcs(image_path)

        # Load image
        logger.info(f"Loading image: {local_image_path}")
        image = np.array(Image.open(local_image_path))
        logger.info(f"Image loaded: shape={image.shape}, dtype={image.dtype}")

        # Get DeepCell engine
        engine = _get_deepcell_engine()

        # Segment image
        logger.info(f"Segmenting with {model_type} model...")
        segmentation_mask, metadata = engine.segment_image(
            image=image,
            model_type=model_type,
            min_cell_size=min_cell_size,
            image_mpp=image_mpp
        )

        # Save segmentation mask as TIFF
        mask_filename = f"{Path(image_path).stem}_segmentation_{model_type}.tif"
        mask_path = OUTPUT_DIR / "segmentations" / mask_filename
        mask_path.parent.mkdir(parents=True, exist_ok=True)

        # Save as 16-bit TIFF (supports up to 65535 cells)
        mask_img = Image.fromarray(segmentation_mask.astype(np.uint16))
        mask_img.save(mask_path)
        logger.info(f"Segmentation mask saved to: {mask_path}")

        result = {
            "status": "success",
            "segmentation_mask": str(mask_path),
            "cells_detected": metadata["cells_detected"],
            "cells_filtered": metadata["cells_filtered"],
            "mean_cell_area": metadata["mean_cell_area"],
            "median_cell_area": metadata["median_cell_area"],
            "min_cell_area": metadata["min_cell_area"],
            "max_cell_area": metadata["max_cell_area"],
            "processing_time_seconds": metadata["processing_time_seconds"],
            "model_used": metadata["model_used"],
            "gpu_used": metadata["gpu_used"],
            "image_mpp": metadata.get("image_mpp"),
            "description": f"Segmented {metadata['cells_detected']} cells using {model_type} model in {metadata['processing_time_seconds']:.1f}s"
        }

        logger.info(f"✅ Segmentation complete: {result['cells_detected']} cells detected")
        return result

    except FileNotFoundError:
        error_msg = f"Image file not found: {image_path}"
        logger.error(error_msg)
        return {
            "status": "error",
            "error": error_msg,
            "message": "Please check the file path and try again."
        }
    except ValueError as e:
        error_msg = f"Invalid parameter: {str(e)}"
        logger.error(error_msg)
        return {
            "status": "error",
            "error": error_msg,
            "message": "Please check your parameters (model_type should be 'nuclear' or 'membrane')"
        }
    except Exception as e:
        error_msg = f"Segmentation failed: {str(e)}"
        logger.error(error_msg, exc_info=True)
        return {
            "status": "error",
            "error": error_msg,
            "message": "An unexpected error occurred during segmentation. Check logs for details."
        }

@mcp.tool()
async def classify_cell_states(
    segmentation_mask_path: str,
    intensity_image_path: str,
    marker_name: str = "Ki67",
    threshold_proliferating: float = 50.0,
    threshold_quiescent: float = 20.0
) -> Dict[str, Any]:
    """Classify cells into functional states based on marker intensity.

    Classifies cells as proliferating, quiescent, or intermediate based on
    marker expression (e.g., Ki67 for proliferation).

    Args:
        segmentation_mask_path: Path to segmentation mask (TIFF from segment_cells)
        intensity_image_path: Path to marker intensity image (TIFF, PNG, etc.)
        marker_name: Name of the marker (default: "Ki67")
        threshold_proliferating: Intensity threshold for proliferating cells (default: 50)
        threshold_quiescent: Intensity threshold for quiescent cells (default: 20)

    Returns:
        Dictionary with:
        - classifications: List of cell states with cell_id, state, confidence
        - total_cells: Total number of cells classified
        - state_counts: Number of cells in each state
        - classifications_csv: Path to exported CSV file
        - description: Summary description

    Example:
        >>> result = await classify_cell_states(
        ...     segmentation_mask_path="/output/segmentations/DAPI_nuclear_seg.tif",
        ...     intensity_image_path="/data/Ki67_intensity.tiff",
        ...     marker_name="Ki67",
        ...     threshold_proliferating=50,
        ...     threshold_quiescent=20
        ... )
        >>> print(f"{result['state_counts']['proliferating']} proliferating cells")
    """
    if DRY_RUN:
        return add_dry_run_warning({
            "classifications": [
                {"cell_id": 1, "state": "proliferating", "confidence": 0.92, "Ki67_intensity": 85.3},
                {"cell_id": 2, "state": "quiescent", "confidence": 0.87, "Ki67_intensity": 12.1},
                {"cell_id": 3, "state": "intermediate", "confidence": 0.75, "Ki67_intensity": 35.2}
            ],
            "total_cells": 1247,
            "state_counts": {
                "proliferating": 523,
                "quiescent": 612,
                "intermediate": 112
            },
            "mode": "dry_run"
        })

    try:
        # Download from GCS if needed
        local_seg_mask_path = _download_from_gcs(segmentation_mask_path)
        local_intensity_path = _download_from_gcs(intensity_image_path)

        # Load segmentation mask
        logger.info(f"Loading segmentation mask: {local_seg_mask_path}")
        seg_mask = np.array(Image.open(local_seg_mask_path))

        # Load intensity image
        logger.info(f"Loading intensity image: {local_intensity_path}")
        intensity_image = np.array(Image.open(local_intensity_path))

        logger.info(f"Mask shape: {seg_mask.shape}, Intensity shape: {intensity_image.shape}")

        # Get classifier
        classifier = _get_intensity_classifier()

        # Measure cell intensities
        logger.info("Measuring cell intensities...")
        intensities_df = classifier.measure_cell_intensities(
            image=intensity_image,
            segmentation_mask=seg_mask
        )

        # Classify cell states
        logger.info(f"Classifying cells using {marker_name} marker...")
        classifications = classifier.classify_cell_states(
            intensities=intensities_df,
            marker_name=marker_name,
            threshold_proliferating=threshold_proliferating,
            threshold_quiescent=threshold_quiescent
        )

        # Save classifications to CSV
        csv_filename = f"classifications_{marker_name}_{Path(segmentation_mask_path).stem}.csv"
        csv_path = OUTPUT_DIR / "classifications" / csv_filename
        csv_path.parent.mkdir(parents=True, exist_ok=True)

        classifications_df = pd.DataFrame(classifications)
        classifications_df.to_csv(csv_path, index=False)
        logger.info(f"Classifications saved to: {csv_path}")

        # Calculate state counts
        state_counts = {}
        for c in classifications:
            state = c["state"]
            state_counts[state] = state_counts.get(state, 0) + 1

        # Create summary
        total_cells = len(classifications)
        states_list = list(state_counts.keys())

        result = {
            "status": "success",
            "classifications": classifications[:100],  # Return first 100 for brevity
            "total_cells": total_cells,
            "state_counts": state_counts,
            "states_identified": states_list,
            "classifications_csv": str(csv_path),
            "marker_used": marker_name,
            "threshold_proliferating": threshold_proliferating,
            "threshold_quiescent": threshold_quiescent,
            "description": f"Classified {total_cells} cells: {state_counts.get('proliferating', 0)} proliferating, {state_counts.get('quiescent', 0)} quiescent, {state_counts.get('intermediate', 0)} intermediate"
        }

        logger.info(f"✅ Classification complete: {total_cells} cells classified")
        return result

    except FileNotFoundError as e:
        error_msg = f"File not found: {str(e)}"
        logger.error(error_msg)
        return {
            "status": "error",
            "error": error_msg,
            "message": "Please check that both the segmentation mask and intensity image paths are correct."
        }
    except ValueError as e:
        error_msg = f"Invalid data: {str(e)}"
        logger.error(error_msg)
        return {
            "status": "error",
            "error": error_msg,
            "message": "Image and segmentation mask shapes don't match or invalid parameters."
        }
    except Exception as e:
        error_msg = f"Classification failed: {str(e)}"
        logger.error(error_msg, exc_info=True)
        return {
            "status": "error",
            "error": error_msg,
            "message": "An unexpected error occurred during classification. Check logs for details."
        }


# ============================================================================
# VISUALIZATION TOOLS
# ============================================================================


@mcp.tool()
async def generate_segmentation_overlay(
    original_image_path: str,
    segmentation_mask_path: str,
    output_filename: Optional[str] = None,
    overlay_color: str = "yellow",
    overlay_alpha: float = 0.4
) -> Dict[str, Any]:
    """Generate segmentation overlay visualization.

    Overlays segmentation boundaries on the original microscopy image.
    Useful for validating segmentation quality and visualizing cell boundaries.

    Args:
        original_image_path: Path to original microscopy image
        segmentation_mask_path: Path to segmentation mask (label image)
        output_filename: Custom output filename (default: segmentation_overlay_TIMESTAMP.png)
        overlay_color: Color for cell boundaries (default: "yellow")
        overlay_alpha: Transparency of overlay (0=transparent, 1=opaque, default: 0.4)

    Returns:
        Dictionary with:
        - output_file: Path to saved overlay image
        - cells_visualized: Number of cells shown
        - description: Text description

    Example:
        >>> result = await generate_segmentation_overlay(
        ...     original_image_path="/data/IF_CD8.tiff",
        ...     segmentation_mask_path="/data/IF_CD8.tiff.seg.tif",
        ...     overlay_color="green"
        ... )
    """
    if DRY_RUN:
        return add_dry_run_warning({
            "output_file": str(OUTPUT_DIR / "visualizations" / "segmentation_overlay_dryrun.png"),
            "cells_visualized": 1247,
            "description": "DRY_RUN: Would generate segmentation overlay",
            "message": "Set DEEPCELL_DRY_RUN=false to generate real visualizations"
        })

    try:
        # Download from GCS if needed
        local_original_path = _download_from_gcs(original_image_path)
        local_seg_mask_path = _download_from_gcs(segmentation_mask_path)

        # Load original image
        original_img = Image.open(local_original_path)
        original_arr = np.array(original_img)

        # Convert grayscale to RGB if needed
        if len(original_arr.shape) == 2:
            original_rgb = color.gray2rgb(original_arr)
        elif original_arr.shape[2] == 1:
            original_rgb = color.gray2rgb(original_arr[:, :, 0])
        else:
            original_rgb = original_arr[:, :, :3]

        # Load segmentation mask
        seg_mask = np.array(Image.open(local_seg_mask_path))

        # Find cell boundaries
        boundaries = measure.find_boundaries(seg_mask, mode='outer')

        # Create color map for overlay
        color_map = {
            'yellow': [255, 255, 0],
            'green': [0, 255, 0],
            'red': [255, 0, 0],
            'cyan': [0, 255, 255],
            'magenta': [255, 0, 255]
        }
        overlay_rgb = color_map.get(overlay_color.lower(), [255, 255, 0])

        # Create overlay
        overlay = original_rgb.copy().astype(float)
        overlay[boundaries] = overlay_rgb

        # Blend original and overlay
        blended = (1 - overlay_alpha) * original_rgb + overlay_alpha * overlay
        blended = np.clip(blended, 0, 255).astype(np.uint8)

        # Plot
        fig, axes = plt.subplots(1, 2, figsize=(12, 6))

        axes[0].imshow(original_rgb)
        axes[0].set_title('Original Image', fontsize=12, fontweight='bold')
        axes[0].axis('off')

        axes[1].imshow(blended)
        axes[1].set_title(f'Segmentation Overlay ({len(np.unique(seg_mask)) - 1} cells)',
                         fontsize=12, fontweight='bold')
        axes[1].axis('off')

        plt.tight_layout()

        # Save
        timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
        if output_filename is None:
            output_filename = f"segmentation_overlay_{timestamp}.png"

        output_path = OUTPUT_DIR / "visualizations" / output_filename
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close(fig)

        num_cells = len(np.unique(seg_mask)) - 1  # Subtract background
        description = f"Segmentation overlay showing {num_cells} segmented cells with {overlay_color} boundaries overlaid on original image."

        return {
            "status": "success",
            "output_file": str(output_path),
            "cells_visualized": int(num_cells),
            "description": description,
            "visualization_type": "segmentation_overlay",
            "overlay_color": overlay_color
        }

    except Exception as e:
        logger.error(f"Error generating segmentation overlay: {e}")
        return {
            "status": "error",
            "error": str(e),
            "message": "Failed to generate segmentation overlay. Check file paths and formats."
        }


@mcp.tool()
async def generate_phenotype_visualization(
    original_image_path: str,
    segmentation_mask_path: str,
    marker_positive_cells: List[int],
    output_filename: Optional[str] = None,
    positive_color: str = "green",
    negative_color: str = "red"
) -> Dict[str, Any]:
    """Generate cell phenotype visualization.

    Colors cells based on marker expression (positive vs negative).
    Useful for visualizing marker-positive cell distribution (e.g., CD8+, Ki67+).

    Args:
        original_image_path: Path to original microscopy image
        segmentation_mask_path: Path to segmentation mask
        marker_positive_cells: List of cell IDs that are marker-positive
        output_filename: Custom output filename (default: phenotype_viz_TIMESTAMP.png)
        positive_color: Color for marker-positive cells (default: "green")
        negative_color: Color for marker-negative cells (default: "red")

    Returns:
        Dictionary with:
        - output_file: Path to saved visualization
        - positive_cells: Number of marker-positive cells
        - negative_cells: Number of marker-negative cells
        - percent_positive: Percentage of positive cells
        - description: Text description

    Example:
        >>> result = await generate_phenotype_visualization(
        ...     original_image_path="/data/IF_Ki67.tiff",
        ...     segmentation_mask_path="/data/IF_Ki67.tiff.seg.tif",
        ...     marker_positive_cells=[1, 5, 12, 18, 25, ...],  # Ki67+ cell IDs
        ...     positive_color="yellow",
        ...     negative_color="blue"
        ... )
    """
    if DRY_RUN:
        return add_dry_run_warning({
            "output_file": str(OUTPUT_DIR / "visualizations" / "phenotype_viz_dryrun.png"),
            "positive_cells": 623,
            "negative_cells": 624,
            "percent_positive": 50.0,
            "description": "DRY_RUN: Would generate phenotype visualization",
            "message": "Set DEEPCELL_DRY_RUN=false to generate real visualizations"
        })

    try:
        # Download from GCS if needed
        local_original_path = _download_from_gcs(original_image_path)
        local_seg_mask_path = _download_from_gcs(segmentation_mask_path)

        # Load images
        original_img = np.array(Image.open(local_original_path))
        seg_mask = np.array(Image.open(local_seg_mask_path))

        # Get unique cell IDs
        cell_ids = np.unique(seg_mask)[1:]  # Exclude background (0)
        positive_set = set(marker_positive_cells)

        # Create colored phenotype mask
        phenotype_mask = np.zeros(seg_mask.shape + (3,), dtype=np.uint8)

        color_map = {
            'green': [0, 255, 0],
            'red': [255, 0, 0],
            'yellow': [255, 255, 0],
            'blue': [0, 0, 255],
            'cyan': [0, 255, 255],
            'magenta': [255, 0, 255]
        }
        pos_color = color_map.get(positive_color.lower(), [0, 255, 0])
        neg_color = color_map.get(negative_color.lower(), [255, 0, 0])

        for cell_id in cell_ids:
            mask = seg_mask == cell_id
            if cell_id in positive_set:
                phenotype_mask[mask] = pos_color
            else:
                phenotype_mask[mask] = neg_color

        # Plot
        fig, axes = plt.subplots(1, 2, figsize=(12, 6))

        # Original image
        if len(original_img.shape) == 2:
            axes[0].imshow(original_img, cmap='gray')
        else:
            axes[0].imshow(original_img)
        axes[0].set_title('Original Image', fontsize=12, fontweight='bold')
        axes[0].axis('off')

        # Phenotype visualization
        axes[1].imshow(phenotype_mask)
        num_positive = len(positive_set.intersection(set(cell_ids)))
        num_negative = len(cell_ids) - num_positive
        pct_positive = 100 * num_positive / len(cell_ids) if len(cell_ids) > 0 else 0

        axes[1].set_title(f'Cell Phenotypes ({pct_positive:.1f}% Positive)',
                         fontsize=12, fontweight='bold')
        axes[1].axis('off')

        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=np.array(pos_color)/255, label=f'Marker+ ({num_positive} cells)'),
            Patch(facecolor=np.array(neg_color)/255, label=f'Marker- ({num_negative} cells)')
        ]
        axes[1].legend(handles=legend_elements, loc='upper right')

        plt.tight_layout()

        # Save
        import pandas as pd
        timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
        if output_filename is None:
            output_filename = f"phenotype_viz_{timestamp}.png"

        output_path = OUTPUT_DIR / "visualizations" / output_filename
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close(fig)

        description = f"Phenotype visualization showing {num_positive} marker-positive cells ({positive_color}) and {num_negative} marker-negative cells ({negative_color}). {pct_positive:.1f}% of cells are positive."

        return {
            "status": "success",
            "output_file": str(output_path),
            "positive_cells": int(num_positive),
            "negative_cells": int(num_negative),
            "total_cells": int(len(cell_ids)),
            "percent_positive": round(float(pct_positive), 2),
            "description": description,
            "visualization_type": "phenotype_coloring",
            "positive_color": positive_color,
            "negative_color": negative_color
        }

    except Exception as e:
        logger.error(f"Error generating phenotype visualization: {e}")
        return {
            "status": "error",
            "error": str(e),
            "message": "Failed to generate phenotype visualization. Check file paths and formats."
        }


@mcp.resource("model://deepcell/membrane")
def get_membrane_model_info() -> str:
    """Membrane segmentation model information."""
    return json.dumps({
        "model": "DeepCell Membrane Segmentation",
        "description": "Deep learning model for cell membrane detection",
        "architecture": "ResNet50 + Feature Pyramid Network",
        "training_data": "10,000+ annotated cell images",
        "accuracy": "95% IoU on test set",
        "use_cases": ["Cell counting", "Morphology analysis", "Cell tracking"]
    }, indent=2)

def main() -> None:
    """Run the MCP mcp-deepcell server."""
    logger.info("Starting mcp-deepcell server...")

    if DRY_RUN:
        logger.warning("=" * 80)
        logger.warning("⚠️  DRY_RUN MODE ENABLED - RETURNING SYNTHETIC DATA")
        logger.warning("⚠️  Results are MOCKED and do NOT represent real analysis")
        logger.warning("⚠️  Set DEEPCELL_DRY_RUN=false for production use")
        logger.warning("=" * 80)
    else:
        logger.info("✅ Real data processing mode enabled (DEEPCELL_DRY_RUN=false)")

    # Get transport and port from environment
    transport = os.getenv("MCP_TRANSPORT", "stdio")
    port = int(os.getenv("PORT", os.getenv("MCP_PORT", "8000")))

    # Run the server with appropriate transport
    if transport in ("sse", "streamable-http"):
        mcp.run(transport=transport, port=port, host="0.0.0.0")
    else:
        mcp.run(transport=transport)

if __name__ == "__main__":
    main()
