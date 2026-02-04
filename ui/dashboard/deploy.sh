#!/bin/bash
# Deploy MCP Cost & Performance Observatory Dashboard to GCP Cloud Run

set -e

# Configuration
PROJECT_ID="${1:-precision-medicine-poc}"
REGION="${2:-us-central1}"
SERVICE_NAME="mcp-dashboard"

SA_NAME="mcp-dashboard-sa"
SA_EMAIL="${SA_NAME}@${PROJECT_ID}.iam.gserviceaccount.com"

echo "=========================================="
echo "Deploying MCP Dashboard to Cloud Run"
echo "=========================================="
echo "Project:  $PROJECT_ID"
echo "Region:   $REGION"
echo "Service:  $SERVICE_NAME"
echo "SA:       $SA_EMAIL"
echo ""

# ── Provision service account (idempotent) ──────────────────────────────
if ! gcloud iam service-accounts describe "$SA_EMAIL" --project "$PROJECT_ID" &>/dev/null; then
    echo "Creating service account: $SA_NAME"
    gcloud iam service-accounts create "$SA_NAME" \
        --project "$PROJECT_ID" \
        --display-name "MCP Dashboard – Cloud Logging reader" \
        --quiet
fi

# Grant logging.logViewer so the dashboard can query Cloud Run request logs
echo "Granting roles/logging.logViewer to $SA_NAME"
gcloud projects add-iam-policy-binding "$PROJECT_ID" \
    --member "serviceAccount:${SA_EMAIL}" \
    --role "roles/logging.viewer" \
    --quiet

# ── Deploy to Cloud Run ──────────────────────────────────────────────────
echo "Building and deploying..."
gcloud run deploy "$SERVICE_NAME" \
    --source . \
    --platform managed \
    --region "$REGION" \
    --project "$PROJECT_ID" \
    --allow-unauthenticated \
    --service-account "$SA_EMAIL" \
    --set-env-vars "GCP_PROJECT_ID=${PROJECT_ID},GCP_REGION=${REGION}" \
    --memory 1Gi \
    --cpu 1 \
    --min-instances 0 \
    --max-instances 5 \
    --timeout 300 \
    --port 8501 \
    --quiet

echo ""
echo "=========================================="
echo "✅ Deployment Complete!"
echo "=========================================="
echo ""

# Get the service URL
SERVICE_URL=$(gcloud run services describe "$SERVICE_NAME" \
    --region "$REGION" \
    --project "$PROJECT_ID" \
    --format 'value(status.url)')

echo "🌐 Your MCP Dashboard is live at:"
echo "   $SERVICE_URL"
echo ""
echo "📊 View logs:"
echo "   gcloud logging read \"resource.type=cloud_run_revision AND resource.labels.service_name=$SERVICE_NAME\" --limit=50 --project=$PROJECT_ID"
echo ""
echo "🔧 Manage service:"
echo "   https://console.cloud.google.com/run/detail/$REGION/$SERVICE_NAME?project=$PROJECT_ID"
echo ""
echo "💡 Usage:"
echo "   Navigate to the URL above to monitor MCP server costs and performance"
echo "   Upload your own metrics JSON or use the included sample data"
echo ""
