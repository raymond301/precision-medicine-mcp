#!/bin/bash
# Deploy mcp-deepcell to Cloud Run
# Usage: ./deploy.sh [PROJECT_ID] [REGION]

set -e

# Configuration
PROJECT_ID="${1:-your-project-id}"
REGION="${2:-us-central1}"
SERVICE_NAME="mcp-deepcell"
IMAGE_NAME="gcr.io/${PROJECT_ID}/${SERVICE_NAME}"

echo "🚀 Deploying ${SERVICE_NAME} to Cloud Run..."
echo "   Project: ${PROJECT_ID}"
echo "   Region: ${REGION}"
echo ""

# Build and push image
echo "📦 Building Docker image..."
docker build -t "${IMAGE_NAME}:latest" .

echo "📤 Pushing to Container Registry..."
docker push "${IMAGE_NAME}:latest"

# Deploy to Cloud Run
echo "☁️  Deploying to Cloud Run..."
gcloud run deploy "${SERVICE_NAME}" \
  --image="${IMAGE_NAME}:latest" \
  --platform=managed \
  --region="${REGION}" \
  --project="${PROJECT_ID}" \
  --allow-unauthenticated \
  --memory=4Gi \
  --cpu=2 \
  --timeout=300 \
  --max-instances=10 \
  --set-env-vars="DEEPCELL_DRY_RUN=false,MCP_TRANSPORT=sse,DEEPCELL_USE_GPU=false,DEEPCELL_OUTPUT_DIR=/tmp/output,DEEPCELL_MODEL_CACHE_DIR=/tmp/models"

echo ""
echo "✅ Deployment complete!"
echo ""
echo "Get service URL:"
echo "  gcloud run services describe ${SERVICE_NAME} --region=${REGION} --format='value(status.url)'"
echo ""
echo "View logs:"
echo "  gcloud run logs read ${SERVICE_NAME} --region=${REGION}"
