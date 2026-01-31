"""DeepCell model management and inference engine.

Handles loading, caching, and running DeepCell models for cell segmentation.
Supports nuclear and membrane segmentation with GPU acceleration.
"""

import logging
import os
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import numpy as np
from PIL import Image

logger = logging.getLogger(__name__)


class DeepCellEngine:
    """Manages DeepCell models and performs cell segmentation inference.

    Features:
    - Model caching (load once, reuse)
    - GPU acceleration (automatic if available)
    - Nuclear and membrane segmentation
    - Image preprocessing and post-processing
    - Quality metrics and filtering
    """

    def __init__(self, model_cache_dir: Optional[Path] = None, use_gpu: bool = True):
        """Initialize DeepCell engine.

        Args:
            model_cache_dir: Directory to cache downloaded models (default: ~/.deepcell/models)
            use_gpu: Enable GPU acceleration if available (default: True)
        """
        self.model_cache_dir = model_cache_dir or Path.home() / ".deepcell" / "models"
        self.model_cache_dir.mkdir(parents=True, exist_ok=True)

        self.use_gpu = use_gpu
        self._models: Dict[str, Any] = {}  # Cached models

        # Configure TensorFlow/GPU
        self._configure_tensorflow()

        logger.info(f"DeepCell engine initialized (GPU: {self.use_gpu})")
        logger.info(f"Model cache: {self.model_cache_dir}")

    def _configure_tensorflow(self) -> None:
        """Configure TensorFlow for GPU/CPU usage."""
        try:
            import tensorflow as tf

            if not self.use_gpu:
                # Force CPU-only mode
                tf.config.set_visible_devices([], 'GPU')
                logger.info("TensorFlow configured for CPU-only mode")
            else:
                # Check GPU availability
                gpus = tf.config.list_physical_devices('GPU')
                if gpus:
                    try:
                        # Enable memory growth to avoid OOM
                        for gpu in gpus:
                            tf.config.experimental.set_memory_growth(gpu, True)
                        logger.info(f"TensorFlow configured with {len(gpus)} GPU(s)")
                    except RuntimeError as e:
                        logger.warning(f"GPU configuration failed: {e}")
                else:
                    logger.info("No GPUs detected, using CPU")

        except ImportError:
            logger.warning("TensorFlow not available, using CPU fallback")

    def load_model(self, model_type: str) -> Any:
        """Load DeepCell model (with caching).

        Args:
            model_type: Model type - "nuclear" or "membrane"

        Returns:
            Loaded DeepCell application model

        Raises:
            ValueError: If model_type is invalid
        """
        if model_type not in ["nuclear", "membrane"]:
            raise ValueError(f"Invalid model_type: {model_type}. Must be 'nuclear' or 'membrane'")

        # Check cache
        if model_type in self._models:
            logger.debug(f"Using cached {model_type} model")
            return self._models[model_type]

        logger.info(f"Loading {model_type} segmentation model...")

        try:
            from deepcell.applications import NuclearSegmentation, Mesmer

            if model_type == "nuclear":
                model = NuclearSegmentation()
            else:  # membrane
                # Mesmer handles both nuclear and whole-cell segmentation
                model = Mesmer()

            self._models[model_type] = model
            logger.info(f"✅ {model_type} model loaded successfully")

            return model

        except Exception as e:
            logger.error(f"Failed to load {model_type} model: {e}")
            raise

    def preprocess_image(self, image: np.ndarray) -> np.ndarray:
        """Preprocess image for DeepCell inference.

        Args:
            image: Input image (H, W) or (H, W, C)

        Returns:
            Preprocessed image (1, H, W, 1) for nuclear or (1, H, W, 2) for membrane
        """
        # Convert to grayscale if needed
        if len(image.shape) == 3:
            # Take first channel or convert RGB to grayscale
            if image.shape[2] == 3:
                from skimage.color import rgb2gray
                image = rgb2gray(image)
            else:
                image = image[:, :, 0]

        # Ensure float32
        if image.dtype != np.float32:
            image = image.astype(np.float32)

        # Normalize to [0, 1]
        if image.max() > 1.0:
            image = image / image.max()

        # Add batch and channel dimensions: (H, W) -> (1, H, W, 1)
        image = np.expand_dims(image, axis=0)  # Add batch
        image = np.expand_dims(image, axis=-1)  # Add channel

        return image

    def segment_image(
        self,
        image: np.ndarray,
        model_type: str = "nuclear",
        min_cell_size: int = 100,
        image_mpp: Optional[float] = None
    ) -> Tuple[np.ndarray, Dict[str, Any]]:
        """Segment cells using DeepCell model.

        Args:
            image: Input image (H, W) or (H, W, C)
            model_type: "nuclear" or "membrane"
            min_cell_size: Minimum cell size in pixels (cells smaller are filtered)
            image_mpp: Microns per pixel (optional, for physical size calculation)

        Returns:
            Tuple of:
            - segmentation_mask: Label image (H, W) where each cell has unique ID
            - metadata: Dict with quality metrics and statistics
        """
        import time
        start_time = time.time()

        logger.info(f"Segmenting image with {model_type} model...")
        logger.debug(f"Image shape: {image.shape}, dtype: {image.dtype}")

        # Load model
        model = self.load_model(model_type)

        # Preprocess
        preprocessed = self.preprocess_image(image)
        logger.debug(f"Preprocessed shape: {preprocessed.shape}")

        # Run inference
        try:
            if model_type == "nuclear":
                # Nuclear segmentation
                segmentation_mask = model.predict(
                    preprocessed,
                    image_mpp=image_mpp,
                    batch_size=1
                )
            else:  # membrane (Mesmer)
                # Mesmer expects (batch, H, W, 2) with nuclear + membrane channels
                # For single channel, duplicate it
                if preprocessed.shape[-1] == 1:
                    preprocessed = np.concatenate([preprocessed, preprocessed], axis=-1)

                segmentation_mask = model.predict(
                    preprocessed,
                    image_mpp=image_mpp,
                    batch_size=1,
                    compartment='whole-cell'  # Get whole cell, not just nucleus
                )

            # Remove batch dimension
            segmentation_mask = np.squeeze(segmentation_mask, axis=0)

            logger.debug(f"Segmentation mask shape: {segmentation_mask.shape}")

        except Exception as e:
            logger.error(f"Segmentation failed: {e}")
            raise

        # Post-process: filter small cells
        filtered_mask, filter_stats = self._filter_small_cells(
            segmentation_mask,
            min_size=min_cell_size
        )

        # Calculate metrics
        cells_detected = len(np.unique(filtered_mask)) - 1  # Exclude background

        # Calculate cell areas
        from skimage.measure import regionprops
        props = regionprops(filtered_mask.astype(int))
        cell_areas = [prop.area for prop in props]

        mean_cell_area = np.mean(cell_areas) if cell_areas else 0
        median_cell_area = np.median(cell_areas) if cell_areas else 0

        processing_time = time.time() - start_time

        metadata = {
            "cells_detected": int(cells_detected),
            "cells_filtered": int(filter_stats["cells_removed"]),
            "mean_cell_area": float(mean_cell_area),
            "median_cell_area": float(median_cell_area),
            "min_cell_area": float(min(cell_areas)) if cell_areas else 0,
            "max_cell_area": float(max(cell_areas)) if cell_areas else 0,
            "processing_time_seconds": float(processing_time),
            "model_used": model_type,
            "image_mpp": image_mpp,
            "gpu_used": self.use_gpu
        }

        logger.info(f"✅ Segmentation complete: {cells_detected} cells in {processing_time:.1f}s")

        return filtered_mask, metadata

    def _filter_small_cells(
        self,
        segmentation_mask: np.ndarray,
        min_size: int
    ) -> Tuple[np.ndarray, Dict[str, int]]:
        """Filter out small cells from segmentation mask.

        Args:
            segmentation_mask: Label image
            min_size: Minimum cell size in pixels

        Returns:
            Tuple of (filtered_mask, stats_dict)
        """
        from skimage.morphology import remove_small_objects

        original_cells = len(np.unique(segmentation_mask)) - 1

        # Remove small objects (cells)
        filtered_mask = remove_small_objects(
            segmentation_mask.astype(int),
            min_size=min_size
        )

        # Relabel to ensure consecutive IDs
        from skimage.measure import label
        filtered_mask = label(filtered_mask)

        filtered_cells = len(np.unique(filtered_mask)) - 1
        cells_removed = original_cells - filtered_cells

        if cells_removed > 0:
            logger.debug(f"Filtered {cells_removed} small cells (< {min_size} pixels)")

        return filtered_mask, {
            "original_cells": original_cells,
            "filtered_cells": filtered_cells,
            "cells_removed": cells_removed
        }

    def get_model_info(self, model_type: str) -> Dict[str, Any]:
        """Get metadata about a DeepCell model.

        Args:
            model_type: "nuclear" or "membrane"

        Returns:
            Dictionary with model information
        """
        model_info = {
            "nuclear": {
                "name": "DeepCell Nuclear Segmentation",
                "description": "Deep learning model for nuclear segmentation in fluorescence microscopy",
                "architecture": "ResNet50 + Feature Pyramid Network (FPN)",
                "training_data": "Thousands of annotated nuclear images across multiple tissue types",
                "input_channels": 1,
                "output": "Nuclear instance segmentation",
                "use_cases": [
                    "Cell counting",
                    "Nuclear morphology analysis",
                    "Cell density quantification",
                    "DAPI/Hoechst stained images"
                ],
                "reference": "https://www.nature.com/articles/s41592-021-01249-6"
            },
            "membrane": {
                "name": "DeepCell Mesmer (Membrane + Nuclear)",
                "description": "Whole-cell segmentation combining nuclear and membrane signals",
                "architecture": "ResNet50 + FPN with multi-scale features",
                "training_data": "Multiplexed imaging data with nuclear + membrane markers",
                "input_channels": 2,
                "output": "Whole-cell instance segmentation",
                "use_cases": [
                    "Whole-cell segmentation",
                    "Cell shape analysis",
                    "Membrane marker quantification",
                    "Multiplexed IF (MxIF) analysis"
                ],
                "reference": "https://www.nature.com/articles/s41587-021-01094-0"
            }
        }

        return model_info.get(model_type, {})

    def segment_image_tiled(
        self,
        image: np.ndarray,
        model_type: str = "nuclear",
        tile_size: int = 2048,
        overlap: int = 128,
        **kwargs
    ) -> Tuple[np.ndarray, Dict[str, Any]]:
        """Segment large images using tiling strategy.

        For images larger than tile_size, splits into overlapping tiles,
        segments each tile, and stitches results together.

        Args:
            image: Input image (H, W) or (H, W, C)
            model_type: "nuclear" or "membrane"
            tile_size: Size of each tile (default: 2048)
            overlap: Overlap between tiles in pixels (default: 128)
            **kwargs: Additional arguments passed to segment_image

        Returns:
            Tuple of (segmentation_mask, metadata)
        """
        h, w = image.shape[:2]

        # If image is small enough, process directly
        if h <= tile_size and w <= tile_size:
            return self.segment_image(image, model_type=model_type, **kwargs)

        logger.info(f"Image too large ({h}×{w}), using tiling with {tile_size}×{tile_size} tiles")

        # Calculate tile positions
        tiles = self._calculate_tile_positions(h, w, tile_size, overlap)
        logger.info(f"Processing {len(tiles)} tiles...")

        # Segment each tile
        tile_masks = []
        for i, (y_start, y_end, x_start, x_end) in enumerate(tiles):
            logger.debug(f"Tile {i+1}/{len(tiles)}: [{y_start}:{y_end}, {x_start}:{x_end}]")

            tile_image = image[y_start:y_end, x_start:x_end]
            tile_mask, _ = self.segment_image(tile_image, model_type=model_type, **kwargs)

            tile_masks.append({
                "mask": tile_mask,
                "position": (y_start, y_end, x_start, x_end)
            })

        # Stitch tiles together
        logger.info("Stitching tiles...")
        final_mask = self._stitch_tiles(tile_masks, (h, w), overlap)

        # Calculate final metadata
        cells_detected = len(np.unique(final_mask)) - 1
        metadata = {
            "cells_detected": int(cells_detected),
            "tiles_processed": len(tiles),
            "tile_size": tile_size,
            "overlap": overlap,
            "model_used": model_type,
            "tiled": True
        }

        logger.info(f"✅ Tiled segmentation complete: {cells_detected} cells from {len(tiles)} tiles")

        return final_mask, metadata

    def _calculate_tile_positions(
        self,
        height: int,
        width: int,
        tile_size: int,
        overlap: int
    ) -> list:
        """Calculate tile positions with overlap."""
        tiles = []
        stride = tile_size - overlap

        for y in range(0, height, stride):
            y_end = min(y + tile_size, height)
            if y_end == height and y > 0:
                y = height - tile_size  # Adjust last tile

            for x in range(0, width, stride):
                x_end = min(x + tile_size, width)
                if x_end == width and x > 0:
                    x = width - tile_size  # Adjust last tile

                tiles.append((y, y_end, x, x_end))

        return tiles

    def _stitch_tiles(
        self,
        tile_masks: list,
        output_shape: Tuple[int, int],
        overlap: int
    ) -> np.ndarray:
        """Stitch tile segmentation masks together.

        Simple strategy: Take center region of each tile to avoid boundary artifacts.
        """
        from skimage.measure import label

        h, w = output_shape
        stitched = np.zeros((h, w), dtype=np.int32)
        current_label = 1

        for tile_data in tile_masks:
            mask = tile_data["mask"]
            y_start, y_end, x_start, x_end = tile_data["position"]

            # For overlapping regions, use center portion
            margin = overlap // 2

            # Relabel tile to avoid ID conflicts
            unique_labels = np.unique(mask)[1:]  # Exclude background
            for old_label in unique_labels:
                mask[mask == old_label] = current_label
                current_label += 1

            # Place tile in stitched image
            # In overlap regions, prefer non-zero values
            tile_region = stitched[y_start:y_end, x_start:x_end]
            stitched[y_start:y_end, x_start:x_end] = np.where(
                tile_region == 0,
                mask,
                tile_region
            )

        # Final relabeling to ensure consecutive IDs
        stitched = label(stitched > 0)

        return stitched

    def __del__(self):
        """Cleanup on engine destruction."""
        if hasattr(self, '_models'):
            self._models.clear()
        logger.debug("DeepCell engine cleaned up")
