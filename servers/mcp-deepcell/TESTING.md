# mcp-deepcell Testing Guide

## Deployment Verification ✅

**Completed:** 2026-01-31
**Service URL:** https://mcp-deepcell-ondu7mwjpa-uc.a.run.app
**Service Status:** Healthy and Running

### Verified Components

1. **Docker Build** ✅
   - Python 3.10 base image
   - All system dependencies installed (gcc, libgomp1, libhdf5-dev)
   - Python packages installed successfully (DeepCell, FastMCP, etc.)

2. **Cloud Run Deployment** ✅
   - Service deployed to us-central1
   - Resource allocation: 4Gi RAM, 2 vCPUs, 300s timeout
   - SSE transport configured on port 8080
   - Environment variables set correctly

3. **Service Health** ✅
   - MCP server started successfully
   - Uvicorn running on http://0.0.0.0:8080
   - SSE endpoint available at /sse
   - Startup probe succeeded

### From Cloud Run Logs
```
INFO: Started server process [1]
INFO: Waiting for application startup.
INFO: Application startup complete.
INFO: Uvicorn running on http://0.0.0.0:8080 (Press CTRL+C to quit)
INFO: Starting MCP server 'deepcell' with transport 'sse' on http://0.0.0.0:8080/sse
Default STARTUP TCP probe succeeded after 1 attempt for container "mcp-deepcell-1" on port 8080.
```

---

## Testing Requirements

### Quick Test (Completed ✅)

**Script:** `test_deployment.py`

Tests basic service availability and creates synthetic test data.

```bash
python3 test_deployment.py
```

**Results:**
- Service is responding
- Synthetic 512x512 nuclear image created
- SSE endpoint available (timeouts are expected without client)

---

## Full Integration Testing (Pending)

### Prerequisites

1. **MCP Client**
   - Claude Desktop with MCP server configuration
   - OR MCP client library (Python/TypeScript)

2. **Test Images**
   - Real microscopy images required:
     - DAPI nuclear staining (16-bit TIFF, 512x512 or larger)
     - Membrane markers for Mesmer model (optional)
     - Cell state markers (Ki67, TP53, etc.) for classification

### MCP Client Configuration

Add to Claude Desktop config (`claude_desktop_config.json`):

```json
{
  "mcpServers": {
    "deepcell": {
      "command": "curl",
      "args": [
        "-N",
        "-H",
        "Accept: text/event-stream",
        "https://mcp-deepcell-ondu7mwjpa-uc.a.run.app/sse"
      ],
      "transport": "sse"
    }
  }
}
```

### Test Cases

#### Test 1: Nuclear Segmentation
**Tool:** `segment_cells`

**Input:**
- `image_path`: Path to DAPI nuclear image (16-bit TIFF)
- `model_type`: "nuclear"
- `min_cell_size`: 100 (pixels)
- `image_mpp`: 0.325 (microns per pixel, optional)

**Expected Output:**
- Segmentation mask saved as 16-bit TIFF
- Metadata with cell count, total area, average cell size
- Processing time (first run ~35s with model download, subsequent ~2-5s)

**Success Criteria:**
- Segmentation mask created successfully
- Cell count reasonable for image (10-500 cells for typical field of view)
- No errors or crashes

#### Test 2: Membrane Segmentation
**Tool:** `segment_cells`

**Input:**
- `image_path`: Path to membrane marker image (16-bit TIFF)
- `model_type`: "membrane"
- `min_cell_size`: 100

**Expected Output:**
- Cell boundary segmentation mask
- Higher cell counts than nuclear segmentation (membrane more accurate)

**Success Criteria:**
- Different segmentation than nuclear model
- Cell boundaries clearly defined

#### Test 3: Cell State Classification
**Tool:** `classify_cell_states`

**Input:**
- `segmentation_mask_path`: Path to segmentation mask from Test 1
- `intensity_image_path`: Path to Ki67 marker image
- `marker_name`: "Ki67"
- `threshold_proliferating`: 50.0 (median intensity)

**Expected Output:**
- CSV file with cell classifications
- State counts (proliferating vs quiescent)
- Per-cell intensity measurements

**Success Criteria:**
- CSV contains all segmented cells
- Classifications reasonable (~10-30% proliferating for normal tissue)
- Intensity values in expected range

#### Test 4: Segmentation Overlay
**Tool:** `generate_segmentation_overlay`

**Input:**
- `image_path`: Original DAPI image
- `segmentation_mask_path`: Segmentation mask from Test 1
- `alpha`: 0.5 (transparency)

**Expected Output:**
- RGB overlay image (PNG or TIFF)
- Cell boundaries visible on original image

**Success Criteria:**
- Visual verification of segmentation quality
- Boundaries align with nuclei in original image

#### Test 5: Phenotype Visualization
**Tool:** `generate_phenotype_visualization`

**Input:**
- `image_path`: Original image
- `segmentation_mask_path`: Segmentation mask
- `classification_csv_path`: Classification CSV from Test 3
- `classification_column`: "cell_state"

**Expected Output:**
- Color-coded visualization (green=proliferating, red=quiescent)

**Success Criteria:**
- Cells colored according to classification
- Visual inspection shows reasonable distribution

---

## Performance Testing

### Benchmarks (Expected on Cloud Run 2 vCPU, 4Gi RAM)

| Image Size | First Request | Subsequent Requests |
|------------|---------------|---------------------|
| 512×512    | ~35s          | ~2s                 |
| 1024×1024  | ~40s          | ~5s                 |
| 2048×2048  | ~50s          | ~10-15s             |
| 4096×4096  | ~90s          | ~30-60s (tiled)     |

**First request** includes model download (~30s) + inference
**Subsequent requests** use cached models

### Test Performance

1. **Cold Start Test**
   - Restart service (scale to zero, wait 5 minutes)
   - Send first request
   - Measure time to first response
   - **Expected:** 30-40s (model download + inference)

2. **Warm Request Test**
   - Send request immediately after first request
   - Measure time to response
   - **Expected:** 2-5s for 512x512 image

3. **Large Image Test**
   - Test with 4096×4096 image
   - Verify tiling strategy works
   - **Expected:** <120s, no out-of-memory errors

4. **Concurrent Requests Test**
   - Send 5 concurrent requests
   - Verify all complete successfully
   - **Expected:** All complete, Cloud Run auto-scales if needed

---

## Model Download Verification

### On First Request

**Check Cloud Run Logs:**
```bash
gcloud logging read "resource.type=cloud_run_revision AND resource.labels.service_name=mcp-deepcell" \
  --limit=100 --format=json --project=precision-medicine-poc
```

**Expected Log Messages:**
- "Downloading DeepCell model: nuclear (or membrane)"
- "Model cached to: /app/data/models/..."
- "Model loaded successfully"
- HTTP requests to DeepCell model zoo

**Timing:**
- Download time: 15-25s (depends on model size and network)
- Load time: 5-10s (model initialization)
- Total overhead: 20-35s on first request only

### On Subsequent Requests

**Expected:**
- "Loading cached model from: /app/data/models/..."
- No download messages
- Fast model loading (<2s)

---

## Troubleshooting

### Service Not Responding

**Check:**
```bash
gcloud run services describe mcp-deepcell \
  --region=us-central1 \
  --project=precision-medicine-poc
```

**Look for:**
- Service status: "Ready"
- URL: https://mcp-deepcell-...
- Latest revision traffic: 100%

### Out of Memory

**Symptoms:**
- 500 errors
- Container crashes in logs
- "Out of memory" messages

**Solutions:**
1. Increase memory limit in cloudbuild.yaml (currently 4Gi)
2. Reduce image size (use tiling for large images)
3. Use smaller batch size

### Slow Performance

**Check:**
1. First request is slow: Normal (model download)
2. All requests slow: Check if model caching is working
3. Large images slow: Expected, verify tiling is enabled

**Optimize:**
- Enable GPU (requires GPU-enabled Dockerfile)
- Increase CPU allocation (currently 2 vCPUs)
- Use minimum instance count (keeps models warm)

### Model Download Fails

**Symptoms:**
- Timeouts on first request
- "Failed to download model" errors
- HTTP 403/404 errors

**Solutions:**
1. Check internet connectivity from Cloud Run
2. Verify DeepCell model zoo is accessible
3. Check firewall rules allow outbound HTTPS

---

## Next Steps

### Immediate (Pending)

1. **Obtain Real Microscopy Images**
   - Source: MxIF imaging pipeline
   - Format: 16-bit TIFF
   - Markers: DAPI, Ki67, TP53, CD45, etc.

2. **Configure MCP Client**
   - Add mcp-deepcell to Claude Desktop config
   - Test connection to SSE endpoint

3. **Run Integration Tests**
   - Execute all 5 test cases above
   - Document results and performance
   - Verify model download and caching

### Short-term

4. **Create Test Dataset**
   - Add sample images to repository
   - Document image sources and preprocessing
   - Create expected outputs for validation

5. **Automated Testing**
   - Write pytest integration tests
   - Add to CI/CD pipeline
   - Test on every deployment

6. **Performance Monitoring**
   - Set up Cloud Monitoring alerts
   - Track request latency, error rate
   - Monitor memory usage patterns

---

## Test Results Log

### 2026-01-31: Initial Deployment Verification ✅

**Tested:**
- Docker build and deployment
- Service startup and health
- Basic service availability

**Results:**
- All systems operational
- Service responding correctly
- Ready for integration testing

**Pending:**
- Full MCP tool testing with real images
- Performance benchmarking
- Model download verification on first use

---

## Contact

**Issues:** Report at [GitHub repository]
**Questions:** See DEPLOYMENT.md for operations guide
