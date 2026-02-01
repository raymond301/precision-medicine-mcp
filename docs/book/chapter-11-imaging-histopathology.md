# Chapter 11: Imaging and Histopathology

*Building mcp-openimagedata for H&E and MxIF image analysis*

---

## Why Imaging Matters

Chapters 7-10 analyzed cells from different perspectives:
- **Chapter 7**: Spatial transcriptomics (gene expression in tissue regions)
- **Chapter 8**: Cell segmentation (single-cell phenotypes from fluorescence)
- **Chapter 9**: Treatment prediction (which drugs will work)
- **Chapter 10**: Quantum fidelity (classification confidence)

**Missing piece**: Tissue morphology. Pathologists diagnose cancer by visual assessment:
- **Necrotic regions**: Dead tissue (pink/pale areas on H&E)
- **Nuclear atypia**: Abnormal nuclei size/shape (cancer hallmark)
- **Cellularity**: Dense vs sparse cell packing
- **Immune infiltration**: Lymphocytes surrounding tumor

The `mcp-openimagedata` server bridges computational analysis with traditional pathology imaging.

---

## Two Imaging Modalities

### H&E (Hematoxylin and Eosin) Staining

**Standard brightfield microscopy**:
- **Hematoxylin**: Blue/purple nuclear stain
- **Eosin**: Pink cytoplasm and extracellular matrix stain
- **Use**: Morphological assessment by pathologists

**PatientOne H&E** shows:
- High-grade serous carcinoma features
- Papillary architecture with necrotic cores
- High nuclear-to-cytoplasmic ratio
- Mitotic figures (proliferation)

### MxIF (Multiplexed Immunofluorescence)

**Fluorescence microscopy with 2-7 antibody markers**:
- **DAPI**: Nuclear stain (blue)
- **Ki67**: Proliferation marker (green)
- **CD8**: Cytotoxic T cells (red)
- **PanCK**: Epithelial/tumor cells (yellow)
- **TP53**: Mutant protein accumulation (cyan)

**Use**: Quantitative single-cell phenotyping (Chapter 8 used MxIF for segmentation).

---

## The Five mcp-openimagedata Tools

### 1. fetch_histology_image

Load H&E or MxIF images from storage.

```python
@mcp.tool()
def fetch_histology_image(image_id: str, stain_type: str = "he") -> dict:
    """Fetch histology image."""
    # Implementation: servers/mcp-openimagedata/src/mcp_openimagedata/image_loader.py
```

**Returns**: Image path, dimensions (4096×4096), metadata (magnification, format).

Full implementation: [`servers/mcp-openimagedata/src/mcp_openimagedata/image_loader.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-openimagedata/src/mcp_openimagedata/image_loader.py)

---

### 2. register_image_to_spatial

Align histology image with spatial transcriptomics spot coordinates.

**Why you need it**: Overlay H&E tissue morphology on 10X Visium gene expression data (Chapter 7).

```python
@mcp.tool()
def register_image_to_spatial(
    image_path: str,
    spatial_coords_path: str,
    registration_method: str = "affine"
) -> dict:
    """Align histology to spatial transcriptomics coordinates.

    Returns: Transform matrix, alignment quality metrics.
    """
```

**Registration methods**:
- **Affine**: Rotation, scaling, translation (handles tissue orientation differences)
- **Deformable**: Non-linear warping (for tissue deformation)

**Output**: Transformation matrix to map spot coordinates → image pixels.

Full implementation: [`servers/mcp-openimagedata/src/mcp_openimagedata/registration.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-openimagedata/src/mcp_openimagedata/registration.py)

---

### 3. extract_image_features

Compute texture and morphology features from H&E regions.

**Features extracted**:
- **Texture**: Haralick features (contrast, correlation, energy, homogeneity)
- **Color**: RGB/HSV histograms (distinguish necrosis from viable tissue)
- **Morphology**: Cellularity estimation (nuclei density per mm²)

```python
@mcp.tool()
def extract_image_features(image_path: str, region_mask: str = None) -> dict:
    """Extract texture/morphology features from H&E image.

    Returns: Feature vectors for ML classification.
    """
```

**Use case**: Train ML classifier to predict necrosis from texture features, validate against spatial HIF1A expression.

Full implementation: [`servers/mcp-openimagedata/src/mcp_openimagedata/feature_extraction.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-openimagedata/src/mcp_openimagedata/feature_extraction.py)

---

### 4. generate_multiplex_composite

Combine MxIF channels into RGB composite for visualization.

**Why you need it**: Raw MxIF has 5-7 separate grayscale channels. Pathologists need false-color composites.

```python
@mcp.tool()
def generate_multiplex_composite(
    channel_paths: dict,  # {"dapi": "path", "ki67": "path", "cd8": "path"}
    output_path: str,
    color_scheme: dict = None  # {"dapi": "blue", "ki67": "green", "cd8": "red"}
) -> dict:
    """Composite MxIF channels into RGB visualization.

    Returns: Composite image path, channel assignments.
    """
```

**Example output**: Blue nuclei (DAPI) + green proliferating cells (Ki67) + red T cells (CD8).

Full implementation: [`servers/mcp-openimagedata/src/mcp_openimagedata/mxif_composite.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-openimagedata/src/mcp_openimagedata/mxif_composite.py)

---

### 5. generate_he_annotation

Annotate H&E images with necrotic regions and cellularity zones.

**Why you need it**: Pathologists mark regions of interest. This tool digitizes those annotations for computational analysis.

```python
@mcp.tool()
def generate_he_annotation(
    image_path: str,
    annotation_type: str = "necrosis"  # "necrosis", "cellularity", "tumor_boundary"
) -> dict:
    """Generate H&E region annotations.

    Returns: Annotation masks, region statistics.
    """
```

**Annotation methods**:
- **Intensity-based**: Necrosis = pale pink regions (low eosin intensity)
- **Texture-based**: High cellularity = high texture contrast (many nuclei)
- **Manual overlay**: Pathologist drawings imported as masks

Full implementation: [`servers/mcp-openimagedata/src/mcp_openimagedata/he_annotation.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-openimagedata/src/mcp_openimagedata/he_annotation.py)

---

## PatientOne Imaging Workflow

**Natural language prompt**:

```
Analyze PatientOne's tumor histology:
- H&E image: data/patient-data/PAT001-OVC-2025/imaging/PAT001_tumor_HE_20x.tiff
- MxIF DAPI: data/patient-data/PAT001-OVC-2025/imaging/PAT001_tumor_IF_DAPI.tiff
- MxIF Ki67: data/patient-data/PAT001-OVC-2025/imaging/PAT001_tumor_IF_KI67.tiff
- MxIF CD8: data/patient-data/PAT001-OVC-2025/imaging/PAT001_tumor_IF_CD8.tiff
- Spatial coordinates: data/patient-data/PAT001-OVC-2025/spatial/visium_spatial_coordinates.csv

Please:
1. Fetch H&E image and annotate necrotic regions
2. Create MxIF composite (DAPI=blue, Ki67=green, CD8=red)
3. Register H&E to spatial coordinates
4. Extract texture features from necrotic vs viable regions
5. Correlate H&E necrosis with spatial HIF1A expression
```

**Claude orchestrates**:

```python
# Step 1: H&E annotation
he_result = openimagedata.generate_he_annotation(
    image_path="data/.../PAT001_tumor_HE_20x.tiff",
    annotation_type="necrosis"
)
# → 23% of tissue area is necrotic

# Step 2: MxIF composite
composite = openimagedata.generate_multiplex_composite(
    channel_paths={
        "dapi": "data/.../PAT001_tumor_IF_DAPI.tiff",
        "ki67": "data/.../PAT001_tumor_IF_KI67.tiff",
        "cd8": "data/.../PAT001_tumor_IF_CD8.tiff"
    },
    color_scheme={"dapi": "blue", "ki67": "green", "cd8": "red"}
)

# Step 3: Register to spatial
registration = openimagedata.register_image_to_spatial(
    image_path="data/.../PAT001_tumor_HE_20x.tiff",
    spatial_coords_path="data/.../visium_spatial_coordinates.csv"
)

# Step 4: Extract features
features = openimagedata.extract_image_features(
    image_path="data/.../PAT001_tumor_HE_20x.tiff",
    region_mask=he_result["necrosis_mask"]
)
# → Necrotic regions: low contrast, low cellularity

# Step 5: Correlate with spatial data
# Compare H&E necrosis locations with HIF1A expression from Chapter 7
# Result: 94% overlap (H&E necrosis matches spatial hypoxia signature)
```

**Analysis time**: ~5 minutes (image loading, registration, feature extraction).

---

## Integration with Other Servers

**mcp-openimagedata connects to**:

### Chapter 7 (Spatial Transcriptomics)
- **Register H&E to Visium spots**: Overlay morphology on gene expression
- **Validate regions**: Necrotic areas (H&E) match HIF1A+ spots (spatial RNA)

### Chapter 8 (Cell Segmentation)
- **MxIF → DeepCell**: Use DAPI channel for nuclear segmentation
- **Phenotype visualization**: Color cells by marker expression on tissue image

### Chapter 10 (Quantum Fidelity)
- **Spatial context**: Tissue architecture informs cell type embeddings
- **TLS identification**: H&E lymphoid aggregates confirm quantum TLS signatures

---

## Example: Necrosis Validation

**Question**: Does H&E-detected necrosis match spatial HIF1A expression (hypoxia signature)?

**Analysis**:
1. **H&E annotation**: 23% tissue area necrotic (pale pink, low cellularity)
2. **Spatial HIF1A expression**: 150 spots (17% of 900 total) show HIF1A > 500
3. **Registration**: Overlay H&E necrosis mask on Visium spots
4. **Correlation**: 141/150 HIF1A+ spots (94%) fall within necrotic regions

**Conclusion**: H&E morphology and molecular markers are concordant → multi-modal validation.

---

## Try It Yourself

### Option 1: Local Testing

```bash
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp/servers/mcp-openimagedata

python -m venv venv
source venv/bin/activate
pip install -e ".[dev]"

# Test with dry-run mode
export IMAGE_DRY_RUN=true
python -m mcp_openimagedata
```

### Option 2: PatientOne Analysis

```bash
# In Claude Desktop:
# "Analyze PAT001 H&E image, annotate necrotic regions, register to spatial coordinates"
```

---

## What You've Built

**mcp-openimagedata provides**:
1. **Image loading**: H&E and MxIF histology
2. **Spatial registration**: Align images to transcriptomics coordinates
3. **Feature extraction**: Texture, color, morphology quantification
4. **MxIF compositing**: Multi-channel visualization
5. **H&E annotation**: Automated necrosis/cellularity detection

**Integration**: Bridges traditional pathology (visual assessment) with computational biology (gene expression, cell phenotypes, spatial patterns).

---

## Next Steps

**Part 3 Complete!** You've built all advanced capability servers:
- Chapter 8: Cell segmentation (DeepCell)
- Chapter 9: Treatment prediction (GEARS)
- Chapter 10: Quantum fidelity (Bayesian UQ)
- Chapter 11: Imaging (H&E + MxIF)

**Part 4: Deployment and Operations** (Chapters 12-14) covers:
- Chapter 12: Cloud deployment on GCP (Cloud Run, Docker, SSE transport)
- Chapter 13: Hospital production deployment (HIPAA compliance, VPC, de-identification)
- Chapter 14: Operations and monitoring (logging, alerts, cost tracking)

---

**Chapter 11 Summary**:
- H&E staining shows tissue morphology (necrosis, cellularity, nuclear atypia)
- MxIF provides quantitative multi-marker imaging (2-7 antibodies)
- Spatial registration overlays images on Visium gene expression coordinates
- PatientOne: 94% concordance between H&E necrosis and HIF1A expression
- Integration validates findings across modalities (imaging + spatial + molecular)

**Files**: [`servers/mcp-openimagedata/src/mcp_openimagedata/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/servers/mcp-openimagedata/src/mcp_openimagedata)
**Tools exposed**: 5 MCP tools (fetch_image, register_to_spatial, extract_features, generate_composite, generate_annotation)
**PatientOne images**: 7 files (H&E + 6 MxIF channels)
