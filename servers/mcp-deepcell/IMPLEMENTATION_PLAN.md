# mcp-deepcell Full Implementation Plan

**Date:** 2026-01-30
**Current Status:** 30% real (visualizations working, segmentation mocked)
**Goal:** 100% real implementation with DeepCell-TF integration

---

## Current State Analysis

### ✅ What's Working (100% Real)
1. **Visualization Tools** (2/4 tools)
   - `generate_segmentation_overlay`: Creates boundary overlays ✅
   - `generate_phenotype_visualization`: Colors cells by marker expression ✅
   - Using: scikit-image, matplotlib, PIL
   - Output: Publication-quality 300 DPI PNG files

2. **Infrastructure**
   - FastMCP server framework ✅
   - DRY_RUN mode with warnings ✅
   - Output directory management ✅
   - Logging and error handling ✅

### ❌ What's Mocked (0% Real)
1. **segment_cells** - Returns synthetic segmentation
   - Currently: Hardcoded mock data
   - Needs: Real DeepCell-TF model inference

2. **classify_cell_states** - Returns synthetic classifications
   - Currently: Hardcoded mock classifications
   - Needs: Real intensity-based phenotyping

3. **Model Management**
   - Currently: No models loaded
   - Needs: Model downloading, caching, loading

---

## Implementation Phases

### Phase 1: DeepCell-TF Integration (Core Functionality) ⭐ START HERE

**Goal:** Real cell segmentation using DeepCell models

**Dependencies to Add:**
```toml
[project]
dependencies = [
    # Existing
    "fastmcp>=0.2.0",
    "pydantic>=2.0.0",
    "numpy>=1.26.0",
    "pandas>=2.2.0",
    "matplotlib>=3.8.0",
    "pillow>=10.0.0",
    "scikit-image>=0.22.0",

    # NEW - DeepCell Integration
    "deepcell>=0.12.0",           # DeepCell library
    "tensorflow>=2.13.0,<2.16.0", # TensorFlow (required by DeepCell)
    "scipy>=1.11.0",               # Required for DeepCell processing
]
```

**Files to Create:**

1. **`src/mcp_deepcell/deepcell_engine.py`** (NEW - 300 lines)
   ```python
   """DeepCell model management and inference."""

   class DeepCellEngine:
       """Manages DeepCell models and performs inference."""

       def __init__(self, model_cache_dir: Path):
           """Initialize with model cache directory."""

       def load_model(self, model_type: str) -> Any:
           """Load DeepCell model (nuclear, membrane, cytoplasm)."""

       def segment_image(
           self,
           image: np.ndarray,
           model_type: str,
           min_cell_size: int = 100
       ) -> np.ndarray:
           """Segment cells using DeepCell model."""

       def get_model_info(self, model_type: str) -> Dict:
           """Get model metadata."""
   ```

2. **`src/mcp_deepcell/intensity_classifier.py`** (NEW - 200 lines)
   ```python
   """Cell state classification based on marker intensity."""

   class IntensityClassifier:
       """Classify cell states from segmentation + intensity data."""

       def measure_cell_intensities(
           self,
           image: np.ndarray,
           segmentation_mask: np.ndarray
       ) -> pd.DataFrame:
           """Measure mean/median intensity per cell."""

       def classify_by_threshold(
           self,
           intensities: pd.DataFrame,
           marker_name: str,
           threshold: float
       ) -> List[Dict]:
           """Classify cells as positive/negative."""

       def classify_multi_marker(
           self,
           intensities: pd.DataFrame,
           thresholds: Dict[str, float]
       ) -> List[Dict]:
           """Classify cells by multiple markers (e.g., Ki67+/TP53+)."""
   ```

**Implementation Steps:**

1. **Add Dependencies** (30 min)
   - Update `pyproject.toml`
   - Test installation: `pip install -e .`
   - Verify TensorFlow/DeepCell imports

2. **Create DeepCellEngine** (4 hours)
   - Model downloading from DeepCell model zoo
   - Model caching (avoid re-downloading)
   - Preprocessing (normalization, resizing)
   - Inference with batch support
   - Post-processing (watershed, size filtering)

3. **Implement segment_cells** (2 hours)
   - Replace mock with real DeepCellEngine calls
   - Handle different model types (nuclear, membrane, cytoplasm)
   - Add quality metrics (mean confidence, filtered cells)
   - Save segmentation masks as TIFF

4. **Create IntensityClassifier** (3 hours)
   - Measure per-cell intensities using regionprops
   - Threshold-based classification
   - Multi-marker phenotyping
   - Export classifications as CSV

5. **Implement classify_cell_states** (1 hour)
   - Replace mock with IntensityClassifier
   - Load expression data (CSV or from image channels)
   - Return real classifications

**Testing:**
- Unit tests for DeepCellEngine
- Integration tests with real images
- Validate against known good segmentations

**Time Estimate:** 10-12 hours

---

### Phase 2: Performance Optimization (Production Ready)

**Goal:** Fast, memory-efficient segmentation for large images

**Optimizations:**

1. **Model Caching** (1 hour)
   - Load models once at server startup
   - Keep in memory for fast inference
   - Singleton pattern for shared engine

2. **Image Tiling** (2 hours)
   - Tile large images (>4096×4096) into smaller patches
   - Process tiles in parallel
   - Stitch results with overlap handling

3. **Batch Processing** (1 hour)
   - Process multiple images in one API call
   - GPU batch inference if available

4. **Memory Management** (1 hour)
   - Clear unused models
   - Process images in chunks
   - Add memory limits

**Files to Modify:**
- `server.py`: Add model preloading at startup
- `deepcell_engine.py`: Add tiling support

**Time Estimate:** 5 hours

---

### Phase 3: Enhanced Features

**Goal:** Advanced capabilities beyond basic segmentation

**New Tools to Add:**

1. **`measure_cell_morphology`** (2 hours)
   ```python
   @mcp.tool()
   async def measure_cell_morphology(
       segmentation_mask: str
   ) -> Dict[str, Any]:
       """Extract morphological features per cell.

       Features: area, perimeter, eccentricity, solidity, etc.
       Useful for identifying abnormal cell shapes.
       """
   ```

2. **`track_cells_over_time`** (4 hours)
   ```python
   @mcp.tool()
   async def track_cells_over_time(
       image_sequence: List[str],
       model_type: str = "nuclear"
   ) -> Dict[str, Any]:
       """Track cells across time series.

       Useful for live-cell imaging, proliferation tracking.
       """
   ```

3. **`export_cell_data`** (1 hour)
   ```python
   @mcp.tool()
   async def export_cell_data(
       segmentation_mask: str,
       output_format: str = "csv"
   ) -> Dict[str, Any]:
       """Export per-cell data (coordinates, area, intensity).

       Formats: CSV, AnnData (.h5ad), CellProfiler
       """
   ```

**Time Estimate:** 7 hours

---

### Phase 4: Integration & Testing

**Goal:** Comprehensive testing and PatientOne integration

**Testing Strategy:**

1. **Unit Tests** (3 hours)
   - Test each tool independently
   - Mock DeepCell models for speed
   - Validate outputs

2. **Integration Tests** (2 hours)
   - Test complete workflows
   - Use real test images
   - Validate visualization outputs

3. **PatientOne Integration** (4 hours)
   - Create MxIF test data for PatientOne
   - Write prompt templates
   - Document complete workflow

**Test Files to Create:**
- `tests/test_deepcell_engine.py`
- `tests/test_intensity_classifier.py`
- `tests/test_segment_cells.py`
- `tests/test_classify_cells.py`
- `tests/test_integration.py`

**Time Estimate:** 9 hours

---

## Total Implementation Timeline

| Phase | Focus | Time | Priority |
|-------|-------|------|----------|
| **Phase 1** | DeepCell Integration | 10-12 hours | ⭐⭐⭐ Critical |
| **Phase 2** | Performance Optimization | 5 hours | ⭐⭐ High |
| **Phase 3** | Enhanced Features | 7 hours | ⭐ Medium |
| **Phase 4** | Testing & Integration | 9 hours | ⭐⭐ High |
| **TOTAL** | | **31-33 hours** | |

**Recommended Approach:** Implement phases sequentially, deploying after each phase.

---

## Dependencies & Requirements

### Python Packages (Updated)
```toml
dependencies = [
    "fastmcp>=0.2.0",
    "pydantic>=2.0.0",
    "numpy>=1.26.0",
    "pandas>=2.2.0",
    "matplotlib>=3.8.0",
    "pillow>=10.0.0",
    "scikit-image>=0.22.0",
    "deepcell>=0.12.0",           # Core DeepCell library
    "tensorflow>=2.13.0,<2.16.0", # TensorFlow backend
    "scipy>=1.11.0",               # Scientific computing
]
```

### System Requirements
- **CPU:** 4+ cores recommended
- **RAM:** 8GB minimum, 16GB recommended
- **GPU:** Optional but recommended for speed (NVIDIA CUDA 11.8+)
- **Storage:** 2GB for model weights

### Environment Variables
```bash
# Existing
DEEPCELL_OUTPUT_DIR=/workspace/output
DEEPCELL_DRY_RUN=false

# New
DEEPCELL_MODEL_CACHE_DIR=/workspace/models  # Model weights cache
DEEPCELL_USE_GPU=true                       # Enable GPU if available
DEEPCELL_MAX_IMAGE_SIZE=4096                # Max image dimension before tiling
```

---

## Technical Challenges & Solutions

### Challenge 1: Model Download Size
**Problem:** DeepCell models are 100-500MB each
**Solution:**
- Download on first use, cache locally
- Include model checksum validation
- Provide pre-downloaded models in Docker image

### Challenge 2: TensorFlow Version Conflicts
**Problem:** TensorFlow can conflict with other packages
**Solution:**
- Pin TensorFlow version: `>=2.13.0,<2.16.0`
- Use separate virtual environment
- Document known conflicts

### Challenge 3: GPU Support
**Problem:** CUDA setup is complex
**Solution:**
- Make GPU optional (CPU fallback)
- Use TensorFlow's automatic device selection
- Document GPU setup separately

### Challenge 4: Large Image Processing
**Problem:** Images >4K×4K can OOM
**Solution:**
- Implement image tiling
- Process tiles with overlap
- Stitch results intelligently

### Challenge 5: Multi-channel Images
**Problem:** MxIF has 4-7 channels, DeepCell expects 1-2
**Solution:**
- Extract nuclear channel (DAPI) for segmentation
- Use segmentation mask to measure other channels
- Support channel selection parameter

---

## Success Criteria

### Phase 1 Complete When:
✅ `segment_cells` returns real segmentation masks
✅ `classify_cell_states` returns real intensity-based classifications
✅ Works with nuclear, membrane, cytoplasm models
✅ Produces valid TIFF label images
✅ Server runs without mocking (DEEPCELL_DRY_RUN=false)

### Production Ready When:
✅ Processes 2K×2K images in <30 seconds
✅ Handles 4K×4K images with tiling
✅ Model caching reduces subsequent calls to <5 seconds
✅ 95%+ test coverage
✅ Comprehensive error handling
✅ Deployed to Cloud Run successfully

### PatientOne Integration When:
✅ MxIF test data available
✅ Complete workflow documented
✅ Prompt templates created
✅ Visualization examples generated

---

## Deployment Considerations

### Docker Image Updates
```dockerfile
# Add to Dockerfile
RUN pip install deepcell tensorflow-cpu  # or tensorflow-gpu

# Pre-download models (optional, speeds up first run)
RUN python -c "from deepcell.applications import NuclearSegmentation; NuclearSegmentation()"
```

### Cloud Run Configuration
```yaml
# Increase memory for TensorFlow
memory: 4Gi
cpu: 2

# Add GPU (optional)
gpu:
  type: nvidia-tesla-t4
  count: 1
```

### Performance Metrics to Monitor
- Model load time (should be <30s)
- Segmentation time per image (should be <1min for 2K×2K)
- Memory usage (should be <4GB without GPU, <8GB with GPU)
- Model cache hit rate (should be >80% after warmup)

---

## Alternative Approaches Considered

### Option A: Use Cellpose Instead of DeepCell
**Pros:** Newer, potentially better performance
**Cons:** Different API, less documentation
**Decision:** Stick with DeepCell (better documented, proven for nuclear/membrane)

### Option B: Use StarDist for Nuclear Segmentation
**Pros:** Fast, accurate for nuclei
**Cons:** Only nuclei, no membrane/cytoplasm
**Decision:** Keep DeepCell for flexibility (supports all 3 compartments)

### Option C: Implement Custom U-Net Model
**Pros:** Full control, optimized for our data
**Cons:** Requires training data, months of work
**Decision:** Use DeepCell (pre-trained, production-ready)

---

## References

- **DeepCell Documentation:** https://deepcell.readthedocs.io/
- **DeepCell GitHub:** https://github.com/vanvalenlab/deepcell-tf
- **DeepCell Paper:** Van Valen et al., Nature Methods 2016
- **Model Zoo:** https://github.com/vanvalenlab/deepcell-tf/releases
- **Example Notebooks:** https://github.com/vanvalenlab/intro-to-deepcell

---

## Next Steps for Implementation

### Immediate (Phase 1 - Week 1)
1. ✅ Review this plan with team
2. ⏭️ Set up development environment with TensorFlow
3. ⏭️ Create `deepcell_engine.py` skeleton
4. ⏭️ Test basic DeepCell model loading
5. ⏭️ Implement simple segmentation (nuclear model only)

### Short-term (Phase 1-2 - Week 2)
6. ⏭️ Complete all 3 model types (nuclear, membrane, cytoplasm)
7. ⏭️ Implement intensity classifier
8. ⏭️ Replace mocks in `segment_cells` and `classify_cell_states`
9. ⏭️ Add model caching and optimization
10. ⏭️ Write unit tests

### Medium-term (Phase 3-4 - Week 3-4)
11. ⏭️ Add advanced features (morphology, tracking, export)
12. ⏭️ Create comprehensive test suite
13. ⏭️ PatientOne MxIF integration
14. ⏭️ Performance benchmarking
15. ⏭️ Deploy to Cloud Run

---

**Questions for Review:**

1. **Priority:** Should we implement Phase 1 only first, or all phases together?
2. **GPU:** Do we want GPU support, or CPU-only is acceptable?
3. **Models:** Nuclear + membrane sufficient, or need cytoplasm too?
4. **Timeline:** Is 31-33 hours reasonable, or should we scope down?
5. **Testing:** What level of test coverage do we need before deploying?

**Decision Needed:** Approve plan and prioritize phases for implementation.
