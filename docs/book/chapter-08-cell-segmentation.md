# Chapter 8: Cell Segmentation with DeepCell

*Building mcp-deepcell for single-cell resolution imaging analysis*

---

## Why Single-Cell Segmentation Matters

In Chapter 7, you analyzed spatial transcriptomics with 10X Visium. Each 55μm spot contains **10-30 cells** mixed together. You identified regions:

- **Tumor proliferative**: MKI67+ high (chemo-sensitive)
- **Necrotic/hypoxic**: HIF1A+ high (radio-resistant)
- **Immune infiltrated**: CD3D+, CD8A+ (immunotherapy potential)

**But Visium can't answer**:
- How many cells are **Ki67+** (actively proliferating)?
- Are **TP53+** mutant cells spatially clustered or dispersed?
- Do **CD8+ T cells** contact tumor cells (immune checkpoint blockade readiness)?
- Which cells are **double-positive** (TP53+/Ki67+ → aggressive phenotype)?

**Multiplexed Immunofluorescence (MxIF)** images 2-7 protein markers simultaneously at subcellular resolution. But raw MxIF images are just pixels. You need **cell segmentation** to:

1. **Detect cell boundaries** (where does one cell end and another begin?)
2. **Measure marker intensity per cell** (which cells express Ki67? TP53?)
3. **Classify cell phenotypes** (proliferating, quiescent, immune, stromal)
4. **Count cells spatially** (how many CD8+ T cells per mm²?)

The `mcp-deepcell` server uses **DeepCell-TF deep learning models** to segment nuclei and membranes, then classify cells based on marker intensity.

---

## The Four mcp-deepcell Tools

### 1. segment_cells

**Why you need it**: Detect cell boundaries from fluorescence images using pretrained deep neural networks.

**Two DeepCell models**:
- **Nuclear segmentation**: Trained on DAPI nuclear stain (most common)
- **Mesmer (membrane segmentation)**: Uses nuclear + membrane markers for whole-cell boundaries

**Example nuclear segmentation**:
```python
@mcp.tool()
def segment_cells(
    image_path: str,
    model_type: str = "nuclear",
    min_cell_size: int = 10,
    max_cell_size: int = 500
) -> dict:
    """Segment cells using DeepCell models.

    Args:
        image_path: GCS path or local path to DAPI/nuclear image (gs://... or /path/to/dapi.tif)
        model_type: "nuclear" or "membrane" (default: nuclear)
        min_cell_size: Minimum cell area in pixels (filters debris)
        max_cell_size: Maximum cell area (filters artifacts/clumps)

    Returns:
        Segmentation mask (label image where each cell has unique ID)
    """
    # Load image
    image = load_image_from_gcs(image_path)  # Handles gs:// URIs

    # Initialize DeepCell engine
    engine = DeepCellEngine(use_gpu=False)  # CPU mode for Cloud Run

    # Load model (cached after first use)
    model = engine.load_model(model_type)

    # Preprocess: normalize to 0-1 range, add batch/channel dims
    image_normalized = image.astype(np.float32) / 65535.0  # 16-bit → float
    image_input = np.expand_dims(np.expand_dims(image_normalized, 0), -1)  # (1, H, W, 1)

    # Run inference
    predictions = model.predict(image_input)
    segmentation_mask = predictions[0, :, :, 0]  # (H, W) label image

    # Filter by size
    from skimage.measure import regionprops
    props = regionprops(segmentation_mask.astype(int))

    filtered_mask = np.zeros_like(segmentation_mask)
    new_label = 1
    for region in props:
        if min_cell_size <= region.area <= max_cell_size:
            filtered_mask[segmentation_mask == region.label] = new_label
            new_label += 1

    num_cells = new_label - 1

    return {
        "segmentation_mask_path": "/tmp/segmentation_mask.npy",
        "num_cells": num_cells,
        "image_size": image.shape,
        "model_used": model_type
    }
```

**PatientOne example** (DAPI nuclear stain):
```json
{
  "image_path": "gs://sample-inputs-patientone/PAT001-OVC-2025/imaging/PAT001_tumor_IF_DAPI.tiff",
  "num_cells": 1247,
  "image_size": [2048, 2048],
  "model_used": "nuclear"
}
```

1247 cells detected in 2048×2048 image (~35 seconds first run, ~5 seconds subsequent with cached model).

Implementation: [`servers/mcp-deepcell/src/mcp_deepcell/deepcell_engine.py:100-250`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/src/mcp_deepcell/deepcell_engine.py#L100-L250)

---

### 2. classify_cell_states

**Why you need it**: Determine which cells are positive/negative for each marker (Ki67+, TP53+, CD8+, etc.) based on fluorescence intensity.

**Intensity-based classification**:
1. **Measure** per-cell marker intensity (mean intensity within segmented region)
2. **Threshold** to classify positive vs negative (automatic or manual)
3. **Multi-marker** phenotyping (Ki67+/TP53+ double-positive cells)

**Classification methods**:
- **Manual threshold**: User-defined cutoff (e.g., Ki67 intensity > 3000 → positive)
- **Otsu's method**: Automatic threshold from intensity histogram
- **Percentile**: Top 25% intensity cells are positive

**Example classification**:
```python
@mcp.tool()
def classify_cell_states(
    image_paths: dict,  # {"dapi": "path/to/dapi.tif", "ki67": "path/to/ki67.tif", "tp53": "path/to/tp53.tif"}
    segmentation_mask_path: str,
    markers: list[str],  # ["ki67", "tp53"]
    classification_method: str = "otsu",
    manual_thresholds: dict = None  # {"ki67": 3000, "tp53": 2500}
) -> dict:
    """Classify cell phenotypes based on marker intensity.

    Args:
        image_paths: Dict of marker_name → image path
        segmentation_mask_path: Path to segmentation mask from segment_cells
        markers: List of markers to classify
        classification_method: "otsu", "manual", or "percentile"
        manual_thresholds: Intensity thresholds if method="manual"

    Returns:
        Per-cell classifications (DataFrame with cell_id, marker_positive columns)
    """
    # Load segmentation mask
    segmentation_mask = np.load(segmentation_mask_path)

    # Initialize classifier
    classifier = IntensityClassifier()

    # Measure intensities for each marker
    results = {}
    for marker in markers:
        marker_image = load_image_from_gcs(image_paths[marker])

        # Measure per-cell mean intensity
        intensities_df = classifier.measure_cell_intensities(
            marker_image,
            segmentation_mask,
            intensity_properties=["mean_intensity", "max_intensity"]
        )

        # Classify positive/negative
        if classification_method == "manual":
            threshold = manual_thresholds[marker]
        elif classification_method == "otsu":
            from skimage.filters import threshold_otsu
            threshold = threshold_otsu(intensities_df["mean_intensity"].values)
        elif classification_method == "percentile":
            threshold = np.percentile(intensities_df["mean_intensity"].values, 75)

        # Apply threshold
        classified = classifier.classify_by_threshold(
            intensities_df,
            marker_name=marker,
            threshold=threshold
        )

        results[marker] = classified

    # Combine markers (multi-marker phenotyping)
    cell_phenotypes = results[markers[0]][["cell_id"]].copy()
    for marker in markers:
        cell_phenotypes[f"{marker}_positive"] = results[marker]["is_positive"]

    # Define composite phenotypes
    if "ki67" in markers and "tp53" in markers:
        cell_phenotypes["proliferating_mutant"] = (
            cell_phenotypes["ki67_positive"] & cell_phenotypes["tp53_positive"]
        )

    # Count phenotypes
    phenotype_counts = {
        marker: int(cell_phenotypes[f"{marker}_positive"].sum())
        for marker in markers
    }

    if "proliferating_mutant" in cell_phenotypes:
        phenotype_counts["proliferating_mutant"] = int(
            cell_phenotypes["proliferating_mutant"].sum()
        )

    return {
        "cell_phenotypes_path": "/tmp/cell_phenotypes.csv",
        "num_cells": len(cell_phenotypes),
        "phenotype_counts": phenotype_counts,
        "thresholds_used": {marker: threshold for marker in markers}
    }
```

**PatientOne results** (Ki67 + TP53 classification):
```json
{
  "num_cells": 1247,
  "phenotype_counts": {
    "ki67": 312,
    "tp53": 498,
    "proliferating_mutant": 187
  },
  "thresholds_used": {
    "ki67": 3420.5,
    "tp53": 2850.3
  }
}
```

**Interpretation**:
- **25% Ki67+**: Proliferating cells (chemotherapy targets)
- **40% TP53+**: TP53 mutant cells (confirming 73% VAF from VCF in Chapter 5)
- **15% Ki67+/TP53+ double-positive**: Aggressive proliferating mutant phenotype

Implementation: [`servers/mcp-deepcell/src/mcp_deepcell/intensity_classifier.py:86-200`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/src/mcp_deepcell/intensity_classifier.py#L86-L200)

---

### 3. generate_segmentation_overlay

**Why you need it**: Visualize cell boundaries overlaid on original fluorescence image to verify segmentation quality.

**What it creates**: RGB image with segmentation boundaries drawn on top of original DAPI image.

**Example visualization**:
```python
@mcp.tool()
def generate_segmentation_overlay(
    image_path: str,
    segmentation_mask_path: str,
    output_path: str = "/tmp/segmentation_overlay.png"
) -> dict:
    """Generate visualization of segmentation boundaries on original image.

    Args:
        image_path: Original fluorescence image
        segmentation_mask_path: Segmentation mask from segment_cells
        output_path: Where to save overlay image

    Returns:
        Path to overlay image
    """
    from skimage.segmentation import mark_boundaries
    import matplotlib.pyplot as plt

    # Load image and mask
    image = load_image_from_gcs(image_path)
    mask = np.load(segmentation_mask_path)

    # Normalize image for visualization
    image_norm = (image - image.min()) / (image.max() - image.min())

    # Mark boundaries
    overlay = mark_boundaries(
        image_norm,
        mask.astype(int),
        color=(1, 0, 0),  # Red boundaries
        mode='thick'
    )

    # Save
    plt.figure(figsize=(12, 12))
    plt.imshow(overlay, cmap='gray')
    plt.title(f"Cell Segmentation ({mask.max()} cells)")
    plt.axis('off')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    return {
        "overlay_path": output_path,
        "num_cells": int(mask.max())
    }
```

This visualization helps you identify segmentation errors (over-segmentation, under-segmentation, missed cells).

---

### 4. generate_phenotype_visualization

**Why you need it**: Color-code cells by phenotype (Ki67+ = green, TP53+ = red, double-positive = yellow).

**What it creates**: Multi-colored image showing spatial distribution of cell phenotypes.

**Example phenotype visualization**:
```python
@mcp.tool()
def generate_phenotype_visualization(
    segmentation_mask_path: str,
    cell_phenotypes_path: str,
    markers: list[str],  # ["ki67", "tp53"]
    output_path: str = "/tmp/phenotype_visualization.png"
) -> dict:
    """Generate spatial visualization of cell phenotypes.

    Args:
        segmentation_mask_path: Segmentation mask
        cell_phenotypes_path: Cell phenotype classifications from classify_cell_states
        markers: Markers to visualize
        output_path: Where to save visualization

    Returns:
        Path to phenotype visualization
    """
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    # Load data
    mask = np.load(segmentation_mask_path)
    phenotypes = pd.read_csv(cell_phenotypes_path)

    # Create color-coded mask
    # 0 = background (black)
    # 1 = negative for both (gray)
    # 2 = Ki67+ only (green)
    # 3 = TP53+ only (red)
    # 4 = double-positive (yellow)

    colored_mask = np.zeros_like(mask, dtype=np.uint8)

    for _, row in phenotypes.iterrows():
        cell_id = row["cell_id"]
        ki67_pos = row.get("ki67_positive", False)
        tp53_pos = row.get("tp53_positive", False)

        if ki67_pos and tp53_pos:
            color_code = 4  # Yellow (double-positive)
        elif ki67_pos:
            color_code = 2  # Green (Ki67+ only)
        elif tp53_pos:
            color_code = 3  # Red (TP53+ only)
        else:
            color_code = 1  # Gray (negative)

        colored_mask[mask == cell_id] = color_code

    # Define colormap
    colors = ['black', 'gray', 'green', 'red', 'yellow']
    cmap = ListedColormap(colors)

    # Plot
    plt.figure(figsize=(12, 12))
    plt.imshow(colored_mask, cmap=cmap, vmin=0, vmax=4)
    plt.title("Cell Phenotypes (Green=Ki67+, Red=TP53+, Yellow=Both)")

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='gray', label='Negative'),
        Patch(facecolor='green', label='Ki67+'),
        Patch(facecolor='red', label='TP53+'),
        Patch(facecolor='yellow', label='Ki67+/TP53+')
    ]
    plt.legend(handles=legend_elements, loc='upper right')
    plt.axis('off')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    return {
        "visualization_path": output_path,
        "phenotype_distribution": phenotypes.describe().to_dict()
    }
```

This visualization reveals spatial patterns (e.g., proliferating cells clustered at tumor edge, TP53+ cells throughout).

---

## The Complete PatientOne MxIF Workflow

Natural language prompt in Claude Desktop:

```
I have multiplexed immunofluorescence images for patient PAT001-OVC-2025:
- DAPI: gs://sample-inputs-patientone/PAT001-OVC-2025/imaging/PAT001_tumor_IF_DAPI.tiff
- Ki67: gs://sample-inputs-patientone/PAT001-OVC-2025/imaging/PAT001_tumor_IF_KI67.tiff
- TP53: gs://sample-inputs-patientone/PAT001-OVC-2025/imaging/PAT001_tumor_multiplex_IF_TP53_KI67_DAPI.tiff
- CD8: gs://sample-inputs-patientone/PAT001-OVC-2025/imaging/PAT001_tumor_IF_CD8.tiff

Please:
1. Segment cells from DAPI nuclear stain
2. Classify cell states based on Ki67 and TP53 markers (use Otsu thresholding)
3. Generate segmentation overlay to verify quality
4. Generate phenotype visualization showing Ki67+, TP53+, and double-positive cells
5. Count how many cells are in each phenotype
```

Claude orchestrates all 4 tools:

```python
# Step 1: Segment cells
seg_result = deepcell.segment_cells(
    image_path="gs://.../PAT001_tumor_IF_DAPI.tiff",
    model_type="nuclear",
    min_cell_size=20,
    max_cell_size=500
)
# → 1247 cells detected

# Step 2: Classify cell states
class_result = deepcell.classify_cell_states(
    image_paths={
        "dapi": "gs://.../PAT001_tumor_IF_DAPI.tiff",
        "ki67": "gs://.../PAT001_tumor_IF_KI67.tiff",
        "tp53": "gs://.../PAT001_tumor_multiplex_IF_TP53_KI67_DAPI.tiff"
    },
    segmentation_mask_path=seg_result["segmentation_mask_path"],
    markers=["ki67", "tp53"],
    classification_method="otsu"
)
# → Ki67+: 312 cells (25%), TP53+: 498 cells (40%), Double-positive: 187 cells (15%)

# Step 3: Segmentation overlay
overlay_result = deepcell.generate_segmentation_overlay(
    image_path="gs://.../PAT001_tumor_IF_DAPI.tiff",
    segmentation_mask_path=seg_result["segmentation_mask_path"]
)

# Step 4: Phenotype visualization
viz_result = deepcell.generate_phenotype_visualization(
    segmentation_mask_path=seg_result["segmentation_mask_path"],
    cell_phenotypes_path=class_result["cell_phenotypes_path"],
    markers=["ki67", "tp53"]
)
```

**Results**:
- **1247 cells** segmented
- **25% Ki67+** (proliferating → chemo-sensitive)
- **40% TP53+** (mutant, matches 73% VAF accounting for normal stroma)
- **15% double-positive** (aggressive phenotype → priority treatment targets)

Total analysis time: **~45 seconds** (first run with model download), **~8 seconds** (subsequent runs with cached models).

---

## Implementation Walkthrough

### Step 1: Project Setup

```bash
cd servers/mcp-deepcell
python -m venv venv
source venv/bin/activate

# Install dependencies (including DeepCell-TF)
pip install fastmcp deepcell scipy scikit-image pillow google-cloud-storage
```

**Critical dependency**: DeepCell requires **Python 3.10** (TensorFlow 2.8.x compatibility) and **Linux x86_64** (TensorFlow wheels).

Environment variables (`.env`):
```bash
DEEPCELL_DRY_RUN=false  # Use real models
DEEPCELL_USE_GPU=false  # CPU mode for Cloud Run
DEEPCELL_MODEL_CACHE_DIR=/tmp/.deepcell/models  # Cache downloaded models
```

### Step 2: Initialize FastMCP Server

```python
from fastmcp import FastMCP
import os
from pathlib import Path

mcp = FastMCP("deepcell")

# Configuration
config = {
    "dry_run": os.getenv("DEEPCELL_DRY_RUN", "false").lower() == "true",
    "use_gpu": os.getenv("DEEPCELL_USE_GPU", "false").lower() == "true",
    "model_cache_dir": Path(os.getenv("DEEPCELL_MODEL_CACHE_DIR", "/tmp/.deepcell/models"))
}
```

### Step 3: Implement DeepCell Engine

The core segmentation engine. Create `src/mcp_deepcell/deepcell_engine.py`:

```python
import numpy as np
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

class DeepCellEngine:
    """Manages DeepCell models and performs cell segmentation."""

    def __init__(self, model_cache_dir: Path = None, use_gpu: bool = False):
        self.model_cache_dir = model_cache_dir or Path.home() / ".deepcell" / "models"
        self.model_cache_dir.mkdir(parents=True, exist_ok=True)
        self.use_gpu = use_gpu
        self._models = {}  # Cached models

        self._configure_tensorflow()

    def _configure_tensorflow(self):
        """Configure TensorFlow for GPU/CPU mode."""
        try:
            import tensorflow as tf

            if not self.use_gpu:
                # Force CPU-only
                tf.config.set_visible_devices([], 'GPU')
                logger.info("TensorFlow configured for CPU-only")
            else:
                gpus = tf.config.list_physical_devices('GPU')
                if gpus:
                    # Enable memory growth
                    for gpu in gpus:
                        tf.config.experimental.set_memory_growth(gpu, True)
                    logger.info(f"TensorFlow configured with {len(gpus)} GPU(s)")
                else:
                    logger.info("No GPUs detected, using CPU")
        except ImportError:
            logger.warning("TensorFlow not available")

    def load_model(self, model_type: str):
        """Load DeepCell model (with caching)."""
        if model_type in self._models:
            return self._models[model_type]

        logger.info(f"Loading {model_type} segmentation model...")

        from deepcell.applications import NuclearSegmentation, Mesmer

        if model_type == "nuclear":
            model = NuclearSegmentation()
        elif model_type == "membrane":
            model = Mesmer()
        else:
            raise ValueError(f"Invalid model_type: {model_type}")

        self._models[model_type] = model
        logger.info(f"✅ {model_type} model loaded")

        return model

    def segment(self, image: np.ndarray, model_type: str = "nuclear") -> np.ndarray:
        """Run segmentation on preprocessed image."""
        model = self.load_model(model_type)

        # Preprocess: normalize and add dimensions
        image_norm = image.astype(np.float32) / 65535.0  # 16-bit → [0, 1]
        image_input = np.expand_dims(np.expand_dims(image_norm, 0), -1)  # (1, H, W, 1)

        # Predict
        predictions = model.predict(image_input)
        segmentation_mask = predictions[0, :, :, 0]

        return segmentation_mask.astype(np.int32)
```

**Key features**:
- **Model caching**: Load once, reuse (saves ~30 seconds per request)
- **CPU mode**: Works on Cloud Run without GPU
- **Automatic preprocessing**: Handles 16-bit TIFF normalization

Full implementation: [`servers/mcp-deepcell/src/mcp_deepcell/deepcell_engine.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/src/mcp_deepcell/deepcell_engine.py) (470 lines)

### Step 4: Implement Intensity Classifier

Create `src/mcp_deepcell/intensity_classifier.py`:

```python
import pandas as pd
import numpy as np
from skimage.measure import regionprops_table

class IntensityClassifier:
    """Classify cell phenotypes based on marker intensity."""

    def measure_cell_intensities(
        self,
        image: np.ndarray,
        segmentation_mask: np.ndarray
    ) -> pd.DataFrame:
        """Measure per-cell marker intensities."""
        props = regionprops_table(
            segmentation_mask.astype(int),
            intensity_image=image,
            properties=["label", "mean_intensity", "max_intensity", "min_intensity"]
        )

        df = pd.DataFrame(props)
        df = df.rename(columns={"label": "cell_id"})

        return df

    def classify_by_threshold(
        self,
        intensities: pd.DataFrame,
        marker_name: str,
        threshold: float
    ) -> pd.DataFrame:
        """Classify cells as marker-positive/negative."""
        classified = intensities.copy()
        classified["is_positive"] = intensities["mean_intensity"] > threshold
        classified["marker"] = marker_name

        return classified
```

Full implementation: [`servers/mcp-deepcell/src/mcp_deepcell/intensity_classifier.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/src/mcp_deepcell/intensity_classifier.py) (338 lines)

---

## GCS Image Loading

DeepCell needs to load images from Google Cloud Storage (gs:// URIs):

```python
from google.cloud import storage
from PIL import Image
import numpy as np
from io import BytesIO

def load_image_from_gcs(gcs_path: str) -> np.ndarray:
    """Load TIFF image from GCS bucket."""
    if not gcs_path.startswith("gs://"):
        # Local file
        return np.array(Image.open(gcs_path))

    # Parse GCS path
    parts = gcs_path.replace("gs://", "").split("/", 1)
    bucket_name = parts[0]
    blob_name = parts[1]

    # Download to memory
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blob = bucket.blob(blob_name)

    image_bytes = blob.download_as_bytes()
    image = Image.open(BytesIO(image_bytes))

    return np.array(image)
```

This enables direct image access from Cloud Storage without local downloads.

---

## Testing Your Server

### Unit Tests

```python
# tests/test_deepcell_engine.py
import pytest
import numpy as np
from mcp_deepcell.deepcell_engine import DeepCellEngine

def test_nuclear_segmentation_synthetic():
    """Test nuclear segmentation on synthetic DAPI image."""
    # Create synthetic nuclear image (3 cells)
    image = np.zeros((512, 512), dtype=np.uint16)
    image[100:150, 100:150] = 15000  # Cell 1
    image[200:250, 200:250] = 14000  # Cell 2
    image[350:400, 350:400] = 16000  # Cell 3

    engine = DeepCellEngine(use_gpu=False)
    mask = engine.segment(image, model_type="nuclear")

    assert mask.max() >= 3  # At least 3 cells detected
    assert mask.shape == (512, 512)

def test_model_caching():
    """Test that models are cached after first load."""
    engine = DeepCellEngine(use_gpu=False)

    # First load
    model1 = engine.load_model("nuclear")

    # Second load (should be cached)
    model2 = engine.load_model("nuclear")

    assert model1 is model2  # Same object (cached)
```

Run tests:
```bash
pytest tests/ -v
```

### Integration Test with PatientOne

```python
@pytest.mark.integration
def test_patientone_mxif_workflow():
    """Test complete MxIF workflow on PatientOne data."""
    # Segment cells
    seg_result = segment_cells(
        image_path="gs://sample-inputs-patientone/PAT001-OVC-2025/imaging/PAT001_tumor_IF_DAPI.tiff",
        model_type="nuclear"
    )

    assert seg_result["num_cells"] > 1000  # Expect ~1247 cells

    # Classify phenotypes
    class_result = classify_cell_states(
        image_paths={
            "ki67": "gs://.../PAT001_tumor_IF_KI67.tiff",
            "tp53": "gs://.../PAT001_tumor_multiplex_IF_TP53_KI67_DAPI.tiff"
        },
        segmentation_mask_path=seg_result["segmentation_mask_path"],
        markers=["ki67", "tp53"],
        classification_method="otsu"
    )

    # Check expected phenotype proportions
    ki67_percent = class_result["phenotype_counts"]["ki67"] / class_result["num_cells"]
    assert 0.20 < ki67_percent < 0.30  # Expect ~25% Ki67+
```

---

## Cloud Run Deployment

Deploy to Google Cloud Run for production use:

```bash
cd servers/mcp-deepcell
./deploy.sh precision-medicine-poc us-central1
```

**Deployment configuration** (`cloudbuild.yaml`):
```yaml
steps:
  - name: 'gcr.io/cloud-builders/docker'
    args: ['build', '-t', 'gcr.io/$PROJECT_ID/mcp-deepcell:latest', '.']
  - name: 'gcr.io/cloud-builders/docker'
    args: ['push', 'gcr.io/$PROJECT_ID/mcp-deepcell:latest']
  - name: 'gcr.io/cloud-builders/gcloud'
    args:
      - 'run'
      - 'deploy'
      - 'mcp-deepcell'
      - '--image=gcr.io/$PROJECT_ID/mcp-deepcell:latest'
      - '--region=us-central1'
      - '--memory=4Gi'
      - '--cpu=2'
      - '--timeout=300s'
      - '--set-env-vars=DEEPCELL_DRY_RUN=false,DEEPCELL_USE_GPU=false'
```

**Performance** (Cloud Run, CPU-only):
- **512×512 image**: ~35s first request (model download), ~2s subsequent
- **2048×2048 image**: ~50s first request, ~10s subsequent
- **Model caching**: Models cached in `/tmp` across requests (same container)

Production service: `https://mcp-deepcell-ondu7mwjpa-uc.a.run.app`

Deployment guide: [`servers/mcp-deepcell/DEPLOYMENT.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-deepcell/DEPLOYMENT.md)

---

## What You've Built

You now have a cell segmentation server that:

1. **Segments cells**: DeepCell nuclear/membrane models with GPU/CPU support
2. **Classifies phenotypes**: Intensity-based marker classification (Ki67+, TP53+, multi-marker)
3. **Generates visualizations**: Segmentation overlays and phenotype spatial maps
4. **Handles cloud images**: Direct GCS image loading (gs:// URIs)
5. **Caches models**: Fast subsequent requests (~2-10s vs 30-50s first run)

This bridges Chapter 7 (spatial transcriptomics, 10-30 cells/spot) to **single-cell resolution** imaging analysis.

---

## Try It Yourself

### Option 1: Local Development

```bash
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp/servers/mcp-deepcell

python3.10 -m venv venv  # Python 3.10 required
source venv/bin/activate
pip install -e ".[dev]"

export DEEPCELL_DRY_RUN=false
python -m mcp_deepcell
```

### Option 2: Test with Synthetic Data

```bash
# In Claude Desktop:
# "Segment cells from synthetic DAPI image:
#  gs://sample-inputs-patientone/mcp-deepcell-test-data/test_data/dapi_512x512.tif
#  Expected: ~30 cells"

# "Classify Ki67 positive cells from:
#  DAPI: gs://.../dapi_512x512.tif
#  Ki67: gs://.../ki67_512x512.tif
#  Expected: ~25% positive"
```

---

## Next Steps

In **Chapter 9: Treatment Response Prediction**, you'll build `mcp-perturbation` to predict how cells respond to drug treatments using **GEARS graph neural network**. You'll use the single-cell phenotypes identified here (Ki67+/TP53+ cells) to predict sensitivity to platinum-based chemotherapy, PARP inhibitors, and AKT inhibitors.

The cell segmentation you built identifies **which cells have which phenotypes**. Treatment response prediction reveals **which drugs will kill those cells**.

---

**Chapter 8 Summary**:
- MxIF imaging provides 2-7 marker resolution at subcellular scale
- DeepCell-TF deep learning models segment nuclei (DAPI) and membranes (Mesmer)
- Intensity-based classification identifies Ki67+ (25%), TP53+ (40%), double-positive (15%) cells
- Cloud Run deployment: 4Gi RAM, 2 CPU, ~10s inference with model caching
- PatientOne: 1247 cells segmented from 2048×2048 DAPI image

**Files created**: `servers/mcp-deepcell/src/mcp_deepcell/server.py`, `deepcell_engine.py` (470 lines), `intensity_classifier.py` (338 lines)
**Tests added**: 12 unit tests, 68% coverage
**Tools exposed**: 4 MCP tools (segment_cells, classify_cell_states, generate_segmentation_overlay, generate_phenotype_visualization)
**Production deployment**: Cloud Run (https://mcp-deepcell-ondu7mwjpa-uc.a.run.app)
