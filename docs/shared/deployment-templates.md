# Deployment Templates (Canonical Reference)

Common deployment patterns for the Precision Medicine MCP platform UI applications.

---

## Local Development

### Prerequisites
- Python 3.11+
- pip or uv package manager
- API key(s): `ANTHROPIC_API_KEY` and/or `GEMINI_API_KEY`

### Setup

```bash
# Navigate to app directory
cd ui/{app-name}

# Create virtual environment
python -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Set API keys
export ANTHROPIC_API_KEY=your_key
export GEMINI_API_KEY=your_key  # Optional

# Run
streamlit run app.py        # Streamlit apps
jupyter notebook            # Jupyter notebook
```

### Tips
- **Hot reload:** Streamlit auto-reloads on code changes
- **Port conflicts:** Use `--server.port 8502` if 8501 is taken
- **Debug:** Add `st.write(variable)` to inspect values

---

## Streamlit Cloud (Free Hosting)

Best for public demos and education.

1. Fork the repository on GitHub
2. Go to [share.streamlit.io](https://share.streamlit.io)
3. Connect your GitHub account
4. Select the repo and `ui/{app-name}/app.py`
5. Add API keys in **Secrets** (Settings > Secrets):
   ```toml
   ANTHROPIC_API_KEY = "your_key"
   GEMINI_API_KEY = "your_key"
   ```
6. Deploy

---

## GCP Cloud Run (Production)

Best for hospital and production use.

### Build and Deploy

```bash
# Build container
gcloud builds submit --tag gcr.io/{PROJECT_ID}/{app-name}

# Deploy to Cloud Run
gcloud run deploy {app-name} \
  --image gcr.io/{PROJECT_ID}/{app-name} \
  --platform managed \
  --region us-central1 \
  --allow-unauthenticated \
  --set-env-vars "ANTHROPIC_API_KEY=your_key"
```

### Security Best Practices

- Use **Secret Manager** for API keys (not environment variables in production)
- Enable **Cloud Armor** WAF for public-facing apps
- Use **VPC connector** for hospital deployments
- Set `--no-allow-unauthenticated` and add OAuth2 Proxy for hospital use
- Set concurrency and scaling limits to control costs

---

**Cost analysis:** [cost-analysis.md](cost-analysis.md)
**HIPAA deployment:** [`docs/for-hospitals/DEPLOYMENT_CHECKLIST.md`](../for-hospitals/DEPLOYMENT_CHECKLIST.md)
