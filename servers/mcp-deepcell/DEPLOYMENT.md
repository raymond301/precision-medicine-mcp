# Cloud Run Deployment Guide

Deploy mcp-deepcell to Google Cloud Run with real DeepCell-TF segmentation.

## Prerequisites

- Google Cloud Project with billing enabled
- gcloud CLI installed and authenticated
- Docker installed (for local builds)
- Cloud Run API enabled
- Container Registry API enabled

## Quick Deploy

### Option 1: Using deploy.sh Script (Recommended)

```bash
cd servers/mcp-deepcell

# Deploy to your project
./deploy.sh YOUR_PROJECT_ID us-central1
```

The script will:
1. Build the Docker image
2. Push to Google Container Registry
3. Deploy to Cloud Run with optimized settings

### Option 2: Using Cloud Build

```bash
# From repository root
gcloud builds submit servers/mcp-deepcell \
  --config=servers/mcp-deepcell/cloudbuild.yaml \
  --project=YOUR_PROJECT_ID
```

This uses Cloud Build to build and deploy automatically.

### Option 3: Manual Deployment

```bash
cd servers/mcp-deepcell

# 1. Build image
docker build -t gcr.io/YOUR_PROJECT_ID/mcp-deepcell:latest .

# 2. Push to registry
docker push gcr.io/YOUR_PROJECT_ID/mcp-deepcell:latest

# 3. Deploy to Cloud Run
gcloud run deploy mcp-deepcell \
  --image=gcr.io/YOUR_PROJECT_ID/mcp-deepcell:latest \
  --platform=managed \
  --region=us-central1 \
  --allow-unauthenticated \
  --memory=4Gi \
  --cpu=2 \
  --timeout=300 \
  --set-env-vars=DEEPCELL_DRY_RUN=false
```

## Configuration

### Resource Requirements

**Minimum (Light usage):**
- Memory: 2Gi
- CPU: 1
- Timeout: 120s

**Recommended (Production):**
- Memory: 4Gi (TensorFlow + DeepCell models)
- CPU: 2 (faster inference)
- Timeout: 300s (large images)

**High-performance:**
- Memory: 8Gi
- CPU: 4
- Consider GPU (T4) for production workloads

### Environment Variables

```bash
# Production settings (recommended)
DEEPCELL_DRY_RUN=false           # Enable real segmentation
MCP_TRANSPORT=sse                 # Use SSE transport
DEEPCELL_USE_GPU=false            # CPU-only on Cloud Run
DEEPCELL_OUTPUT_DIR=/tmp/output   # Temporary storage
DEEPCELL_MODEL_CACHE_DIR=/tmp/models  # Model cache

# Development/Testing
DEEPCELL_DRY_RUN=true            # Mock data for testing
```

### Enable GPU (Optional)

For production workloads with high throughput:

```bash
gcloud run deploy mcp-deepcell \
  --image=gcr.io/YOUR_PROJECT_ID/mcp-deepcell:latest \
  --platform=managed \
  --region=us-central1 \
  --memory=8Gi \
  --cpu=4 \
  --gpu=1 \
  --gpu-type=nvidia-tesla-t4 \
  --set-env-vars=DEEPCELL_USE_GPU=true
```

**Note:** GPU support requires:
- Cloud Run GPU-enabled regions (us-central1, us-west1, etc.)
- Higher costs (~$0.35/hour for T4)
- Modified Dockerfile with tensorflow-gpu

## Verify Deployment

```bash
# Get service URL
SERVICE_URL=$(gcloud run services describe mcp-deepcell \
  --region=us-central1 \
  --format='value(status.url)')

echo "Service URL: $SERVICE_URL"

# Test health check
curl "${SERVICE_URL}/health"

# View logs
gcloud run logs read mcp-deepcell --region=us-central1 --limit=50
```

## Usage After Deployment

### Connect from Claude Desktop

Add to your Claude Desktop MCP config:

```json
{
  "mcpServers": {
    "deepcell": {
      "url": "https://YOUR_SERVICE_URL",
      "transport": "sse"
    }
  }
}
```

### API Endpoints

- `POST /segment_cells` - Segment cells in microscopy images
- `POST /classify_cell_states` - Classify cell phenotypes
- `POST /generate_segmentation_overlay` - Visualize segmentation
- `POST /generate_phenotype_visualization` - Visualize phenotypes

### Example Request

```bash
curl -X POST "${SERVICE_URL}/segment_cells" \
  -H "Content-Type: application/json" \
  -d '{
    "image_path": "/path/to/image.tiff",
    "model_type": "nuclear",
    "min_cell_size": 100
  }'
```

## Monitoring

### View Metrics

```bash
# Cloud Console
https://console.cloud.google.com/run/detail/us-central1/mcp-deepcell/metrics

# CLI
gcloud run services describe mcp-deepcell \
  --region=us-central1
```

### Key Metrics to Monitor

- **Request latency:** Should be <30s for 2K×2K images
- **Memory usage:** Should stay <4Gi for normal workloads
- **Error rate:** Should be <1%
- **Cold start time:** First model load takes ~30s

### Logs

```bash
# Tail logs in real-time
gcloud run logs tail mcp-deepcell --region=us-central1

# Filter errors only
gcloud run logs read mcp-deepcell \
  --region=us-central1 \
  --filter="severity>=ERROR"
```

## Troubleshooting

### Issue: Out of Memory

**Symptoms:** Container crashes, "OOM killed" errors

**Solution:**
```bash
# Increase memory
gcloud run services update mcp-deepcell \
  --memory=8Gi \
  --region=us-central1
```

### Issue: Timeout on Large Images

**Symptoms:** "Deadline exceeded" errors

**Solution:**
```bash
# Increase timeout
gcloud run services update mcp-deepcell \
  --timeout=600 \
  --region=us-central1
```

### Issue: Slow Model Loading

**Symptoms:** First requests take 30-60s

**Solution:**
- This is normal (model download + load)
- Models are cached in /tmp/models after first load
- Use minimum instances to keep warm:

```bash
gcloud run services update mcp-deepcell \
  --min-instances=1 \
  --region=us-central1
```

### Issue: Cannot Install Dependencies

**Symptoms:** Build fails with TensorFlow errors

**Solution:**
- Ensure using Python 3.10 (not 3.11+)
- Check Dockerfile uses `python:3.10-slim`
- Verify system dependencies (libgomp1, libhdf5-dev)

## Cost Optimization

### Estimated Costs

**CPU-only (no minimum instances):**
- $0.00002400/vCPU-second + $0.00000250/GiB-second
- ~$0.05-0.20 per 1000 requests (2-3s each)
- Free tier: 2 million requests/month

**With 1 minimum instance (always warm):**
- ~$50-80/month (4Gi, 2 CPU)

**With GPU (T4):**
- ~$250-350/month (continuous) + request costs

### Reduce Costs

1. **Use DRY_RUN for development:**
   ```bash
   gcloud run services update mcp-deepcell \
     --set-env-vars=DEEPCELL_DRY_RUN=true
   ```

2. **Set max instances:**
   ```bash
   gcloud run services update mcp-deepcell \
     --max-instances=3
   ```

3. **Use smaller resources for testing:**
   ```bash
   gcloud run services update mcp-deepcell \
     --memory=2Gi --cpu=1
   ```

## Security

### Restrict Access

For production, enable authentication:

```bash
gcloud run services update mcp-deepcell \
  --no-allow-unauthenticated \
  --region=us-central1
```

Then use service account or API keys for access.

### VPC Integration

For private deployment:

```bash
gcloud run services update mcp-deepcell \
  --vpc-connector=YOUR_CONNECTOR \
  --vpc-egress=private-ranges-only \
  --region=us-central1
```

## Next Steps

- [ ] Set up CI/CD with Cloud Build triggers
- [ ] Configure custom domain
- [ ] Set up monitoring alerts
- [ ] Enable Cloud Armor for DDoS protection
- [ ] Implement request authentication
- [ ] Add Cloud Storage integration for image uploads

## Support

For issues, see:
- `DEPENDENCY_ISSUES.md` - Platform compatibility
- `IMPLEMENTATION_PLAN.md` - Feature roadmap
- GitHub Issues: https://github.com/YOUR_REPO/issues
