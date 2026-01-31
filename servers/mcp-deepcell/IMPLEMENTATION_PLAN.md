# mcp-deepcell Implementation Plan

**Date:** 2026-01-31
**Status:** ✅ Phase 1 Complete (Production Ready)

---

## ✅ Phase 1: DeepCell-TF Integration (COMPLETE)

**Real cell segmentation and intensity-based phenotyping**

**Implemented:**
- `deepcell_engine.py` (470 lines): Model management, inference, preprocessing, tiling
- `intensity_classifier.py` (338 lines): Marker-based cell state classification
- Real segmentation: Nuclear, membrane models with GPU/CPU support
- Real classification: Intensity thresholds for multi-marker phenotyping
- GCS support: Direct image loading from `gs://` URIs
- Quality metrics: Cell counts, filtering, confidence scores

**Dependencies Added:**
```toml
"deepcell>=0.12.0"
"scipy>=1.11.0"
"google-cloud-storage>=2.10.0"
```

**Platform:** Python 3.10, TensorFlow 2.8.x, Linux x86_64 (Cloud Run ready)

---

## Phase 2: Advanced Features (Optional)

**Enhanced capabilities beyond basic segmentation**

**Potential Tools:**

1. **`measure_cell_morphology`** - Extract shape features (area, eccentricity, solidity)
2. **`track_cells_over_time`** - Cell tracking across time series
3. **`export_cell_data`** - Export to CSV, AnnData (.h5ad), CellProfiler formats

**Time Estimate:** 6-8 hours

---

## Phase 3: Comprehensive Testing (Recommended)

**Production-grade testing and validation**

**Test Coverage:**
- Unit tests for DeepCellEngine and IntensityClassifier
- Integration tests with real GCS test images
- PatientOne MxIF workflow validation
- Performance benchmarks (512x512: ~5s, 2048x2048: ~30s)

**Test Data:** `gs://sample-inputs-patientone/mcp-deepcell-test-data/`

**Time Estimate:** 4-6 hours

---

## Success Criteria

**Production Ready (Phase 1):** ✅
- Real segmentation masks from DeepCell models
- Real intensity-based classifications
- GCS image loading support
- Deployed to Cloud Run
- Works with nuclear and membrane models
- Handles large images (tiling for >2048x2048)

**Fully Tested (Phase 3):**
- 90%+ test coverage
- Validated against known segmentations
- PatientOne workflow documented with prompts

---

## Technical Decisions

**DeepCell vs Alternatives:**
- Chose DeepCell-TF over Cellpose/StarDist for proven nuclear+membrane support
- Pre-trained models (no custom training required)

**Performance:**
- Model caching in `/tmp` (fast subsequent requests)
- CPU-only mode (GPU optional)
- Automatic tiling for large images

**Deployment:**
- Cloud Run: 4Gi RAM, 2 CPU, 300s timeout
- Environment: `DEEPCELL_DRY_RUN=false` for production

---

## Next Steps

1. **Optional:** Implement Phase 2 advanced features (morphology, tracking, export)
2. **Recommended:** Add comprehensive test suite (Phase 3)
3. **Documentation:** Update imaging.md with real usage examples

---

## References

- **DeepCell-TF:** https://github.com/vanvalenlab/deepcell-tf
- **Models:** Nuclear (Mesmer), Membrane (Mesmer whole-cell)
- **Test Data:** `gs://sample-inputs-patientone/mcp-deepcell-test-data/test_data/`
