"""Cell state classification based on marker intensity.

Classifies cells into phenotypes (positive/negative, proliferating/quiescent, etc.)
based on fluorescence marker intensity measurements.
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from PIL import Image
from skimage.measure import regionprops_table

logger = logging.getLogger(__name__)


class IntensityClassifier:
    """Classify cell states/phenotypes from segmentation masks and marker intensity.

    Features:
    - Measure per-cell marker intensities
    - Threshold-based classification (positive/negative)
    - Multi-marker phenotyping (e.g., Ki67+/TP53+)
    - Statistical classification methods
    """

    def __init__(self):
        """Initialize intensity classifier."""
        logger.info("IntensityClassifier initialized")

    def measure_cell_intensities(
        self,
        image: np.ndarray,
        segmentation_mask: np.ndarray,
        intensity_properties: List[str] = None
    ) -> pd.DataFrame:
        """Measure intensity statistics for each segmented cell.

        Args:
            image: Intensity image (H, W) or (H, W, C) for multiple channels
            segmentation_mask: Label image where each cell has unique ID
            intensity_properties: Properties to measure (default: ["mean_intensity", "max_intensity"])

        Returns:
            DataFrame with columns: cell_id, mean_intensity, max_intensity, min_intensity, etc.
        """
        if intensity_properties is None:
            intensity_properties = ["mean_intensity", "max_intensity", "min_intensity"]

        logger.debug(f"Measuring {len(intensity_properties)} intensity properties for cells")

        # Handle multi-channel images
        if len(image.shape) == 3:
            # Take first channel or average
            if image.shape[2] > 1:
                logger.debug(f"Multi-channel image detected ({image.shape[2]} channels), using channel 0")
                intensity_image = image[:, :, 0]
            else:
                intensity_image = image[:, :, 0]
        else:
            intensity_image = image

        # Ensure same shape
        if intensity_image.shape != segmentation_mask.shape:
            raise ValueError(
                f"Image shape {intensity_image.shape} doesn't match "
                f"segmentation mask shape {segmentation_mask.shape}"
            )

        # Measure properties using regionprops
        props = regionprops_table(
            segmentation_mask.astype(int),
            intensity_image=intensity_image,
            properties=["label"] + intensity_properties
        )

        df = pd.DataFrame(props)
        df = df.rename(columns={"label": "cell_id"})

        logger.info(f"✅ Measured intensities for {len(df)} cells")

        return df

    def classify_by_threshold(
        self,
        intensities: pd.DataFrame,
        marker_name: str,
        threshold: float,
        intensity_column: str = "mean_intensity"
    ) -> pd.DataFrame:
        """Classify cells as marker-positive or marker-negative based on intensity threshold.

        Args:
            intensities: DataFrame from measure_cell_intensities()
            marker_name: Name of marker (e.g., "Ki67", "CD8", "TP53")
            threshold: Intensity threshold for positive classification
            intensity_column: Column to use for thresholding (default: "mean_intensity")

        Returns:
            DataFrame with added columns: {marker_name}_positive (bool), {marker_name}_intensity (float)
        """
        if intensity_column not in intensities.columns:
            raise ValueError(f"Column '{intensity_column}' not found in intensities DataFrame")

        df = intensities.copy()

        # Add marker intensity column
        df[f"{marker_name}_intensity"] = df[intensity_column]

        # Classify as positive/negative
        df[f"{marker_name}_positive"] = df[intensity_column] > threshold

        n_positive = df[f"{marker_name}_positive"].sum()
        n_total = len(df)
        pct_positive = 100 * n_positive / n_total if n_total > 0 else 0

        logger.info(
            f"✅ {marker_name} classification: {n_positive}/{n_total} positive ({pct_positive:.1f}%)"
        )

        return df

    def classify_multi_marker(
        self,
        intensities: pd.DataFrame,
        thresholds: Dict[str, float],
        intensity_column: str = "mean_intensity"
    ) -> pd.DataFrame:
        """Classify cells by multiple markers simultaneously.

        Args:
            intensities: DataFrame from measure_cell_intensities()
            thresholds: Dict mapping marker_name -> threshold value
            intensity_column: Column to use for thresholding

        Returns:
            DataFrame with boolean columns for each marker and phenotype column

        Example:
            >>> thresholds = {"Ki67": 50, "TP53": 100}
            >>> classified = classifier.classify_multi_marker(intensities, thresholds)
            >>> # Result has columns: Ki67_positive, TP53_positive, phenotype
            >>> # phenotype values: "Ki67+/TP53+", "Ki67+/TP53-", "Ki67-/TP53+", "Ki67-/TP53-"
        """
        df = intensities.copy()

        # Classify each marker
        for marker_name, threshold in thresholds.items():
            df = self.classify_by_threshold(
                df,
                marker_name=marker_name,
                threshold=threshold,
                intensity_column=intensity_column
            )

        # Create phenotype label (e.g., "Ki67+/TP53+")
        phenotype_parts = []
        for marker_name in thresholds.keys():
            pos_col = f"{marker_name}_positive"
            phenotype_parts.append(
                df[pos_col].apply(lambda x: f"{marker_name}+" if x else f"{marker_name}-")
            )

        df["phenotype"] = phenotype_parts[0]
        for part in phenotype_parts[1:]:
            df["phenotype"] = df["phenotype"] + "/" + part

        # Log phenotype distribution
        phenotype_counts = df["phenotype"].value_counts()
        logger.info(f"✅ Multi-marker phenotypes identified:")
        for phenotype, count in phenotype_counts.items():
            logger.info(f"   {phenotype}: {count} cells ({100*count/len(df):.1f}%)")

        return df

    def classify_cell_states(
        self,
        intensities: pd.DataFrame,
        marker_name: str = "Ki67",
        threshold_proliferating: float = 50,
        threshold_quiescent: float = 20
    ) -> List[Dict[str, Any]]:
        """Classify cells into functional states (proliferating, quiescent, etc.).

        Args:
            intensities: DataFrame from measure_cell_intensities()
            marker_name: Proliferation marker (default: "Ki67")
            threshold_proliferating: Intensity above which cells are proliferating
            threshold_quiescent: Intensity below which cells are quiescent

        Returns:
            List of dicts with cell_id, state, confidence
        """
        results = []

        for _, row in intensities.iterrows():
            cell_id = int(row["cell_id"])
            intensity = row["mean_intensity"]

            # Classify state
            if intensity > threshold_proliferating:
                state = "proliferating"
                confidence = min(1.0, (intensity - threshold_proliferating) / 100)
            elif intensity < threshold_quiescent:
                state = "quiescent"
                confidence = min(1.0, (threshold_quiescent - intensity) / 100)
            else:
                # Intermediate - unclear state
                state = "intermediate"
                mid_point = (threshold_proliferating + threshold_quiescent) / 2
                distance = abs(intensity - mid_point)
                confidence = 1.0 - min(1.0, distance / mid_point)

            results.append({
                "cell_id": cell_id,
                "state": state,
                "confidence": float(max(0.5, min(1.0, confidence))),  # Clamp to [0.5, 1.0]
                f"{marker_name}_intensity": float(intensity)
            })

        # Count states
        state_counts = {}
        for r in results:
            state_counts[r["state"]] = state_counts.get(r["state"], 0) + 1

        logger.info(f"✅ Cell states identified:")
        for state, count in state_counts.items():
            logger.info(f"   {state}: {count} cells")

        return results

    def auto_threshold_otsu(self, intensities: pd.DataFrame) -> float:
        """Calculate automatic threshold using Otsu's method.

        Args:
            intensities: DataFrame with intensity measurements

        Returns:
            Optimal threshold value
        """
        from skimage.filters import threshold_otsu

        values = intensities["mean_intensity"].values
        threshold = threshold_otsu(values)

        logger.info(f"Otsu threshold: {threshold:.2f}")

        return float(threshold)

    def auto_threshold_percentile(
        self,
        intensities: pd.DataFrame,
        percentile: float = 75
    ) -> float:
        """Calculate threshold as a percentile of intensities.

        Args:
            intensities: DataFrame with intensity measurements
            percentile: Percentile to use (default: 75 = top 25% are positive)

        Returns:
            Threshold value at specified percentile
        """
        values = intensities["mean_intensity"].values
        threshold = np.percentile(values, percentile)

        logger.info(f"{percentile}th percentile threshold: {threshold:.2f}")

        return float(threshold)

    def get_marker_positive_cells(
        self,
        classified_df: pd.DataFrame,
        marker_name: str
    ) -> List[int]:
        """Get list of cell IDs that are positive for a marker.

        Args:
            classified_df: DataFrame from classify_by_threshold() or classify_multi_marker()
            marker_name: Marker name

        Returns:
            List of cell IDs (integers)
        """
        pos_col = f"{marker_name}_positive"

        if pos_col not in classified_df.columns:
            raise ValueError(f"Marker '{marker_name}' not found in classified DataFrame")

        positive_cells = classified_df[classified_df[pos_col]]["cell_id"].tolist()

        logger.debug(f"Found {len(positive_cells)} {marker_name}+ cells")

        return [int(c) for c in positive_cells]

    def export_classifications(
        self,
        classified_df: pd.DataFrame,
        output_path: Path
    ) -> None:
        """Export classifications to CSV file.

        Args:
            classified_df: DataFrame with classifications
            output_path: Path to save CSV file
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        classified_df.to_csv(output_path, index=False)

        logger.info(f"✅ Exported classifications to {output_path}")

    def load_expression_matrix(self, expression_path: str) -> pd.DataFrame:
        """Load expression matrix from CSV file.

        Expected format:
        - CSV with columns: cell_id, marker1, marker2, ...
        - Or: cell_id, mean_intensity (for single marker)

        Args:
            expression_path: Path to CSV file

        Returns:
            DataFrame with cell expression data
        """
        df = pd.read_csv(expression_path)

        if "cell_id" not in df.columns:
            # Assume first column is cell_id
            df = df.rename(columns={df.columns[0]: "cell_id"})

        logger.info(f"✅ Loaded expression data: {len(df)} cells, {len(df.columns)-1} markers")

        return df
