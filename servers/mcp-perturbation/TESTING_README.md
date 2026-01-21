# Testing mcp-perturbation Deployment

**Quick Start**: Once deployment completes, run the test script to verify everything works.

---

## 🚀 Quick Test (30 seconds)

```bash
cd servers/mcp-perturbation
./test_deployment.sh
```

**What it tests:**
- ✅ Service connectivity
- ✅ MCP protocol handshake
- ✅ Tool registration (8 tools)
- ✅ Tool execution

---

## 📚 Full Testing Guide

See [DEPLOYMENT_TESTING.md](DEPLOYMENT_TESTING.md) for comprehensive tests including:
- Health checks
- SSE endpoint testing
- MCP protocol validation
- Error handling
- Performance benchmarks
- Concurrent request handling

---

## 🖥️ Claude Desktop Configuration

### Option 1: Local Server (stdio)

Add to your `~/Library/Application Support/Claude/claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "mcp-perturbation": {
      "command": "python",
      "args": ["-m", "mcp_perturbation"],
      "env": {
        "MCP_TRANSPORT": "stdio",
        "PERTURBATION_DRY_RUN": "false",
        "PERTURBATION_LOG_LEVEL": "INFO"
      }
    }
  }
}
```

**Requirements:**
- mcp-perturbation installed: `cd servers/mcp-perturbation && pip install -e .`
- Python 3.10+ with all dependencies

**Pros:**
- Faster (no network latency)
- No GCP costs
- Full access to local files

**Cons:**
- Requires local installation
- Limited to your machine's resources

---

### Option 2: Remote Server (SSE via Cloud Run)

Add to your `~/Library/Application Support/Claude/claude_desktop_config.json`:

```json
{
  "mcpServers": {
    "mcp-perturbation": {
      "command": "npx",
      "args": [
        "-y",
        "@modelcontextprotocol/client-sse",
        "https://mcp-perturbation-ondu7mwjpa-uc.a.run.app/sse"
      ],
      "env": {
        "NODE_OPTIONS": "--max-old-space-size=4096"
      }
    }
  }
}
```

**Requirements:**
- Node.js installed
- Internet connection
- GCP deployment complete

**Pros:**
- No local installation needed
- Access from anywhere
- Scales automatically

**Cons:**
- Network latency (~100-500ms)
- GCP costs ($0.09/hour when active)
- Requires internet connection

---

## 🧪 Testing Scenarios

### Scenario 1: Test Here (Claude Code CLI)

```bash
# Start local server
cd servers/mcp-perturbation
python -m mcp_perturbation

# In another terminal, test with MCP client or curl
```

### Scenario 2: Test with Claude Desktop

1. Add configuration (see above)
2. Restart Claude Desktop
3. Verify server appears in available tools
4. Test: "Load dataset GSE184880 using perturbation server"

### Scenario 3: Test with Streamlit UI (Local)

```bash
# Configure Streamlit to use remote server
export MCP_SERVER_URL="https://mcp-perturbation-ondu7mwjpa-uc.a.run.app/sse"

# Start Streamlit
cd ui/jupyter-notebook
streamlit run streamlit_app.py
```

### Scenario 4: Test End-to-End (Remote)

Both Streamlit UI and mcp-perturbation running on GCP - full cloud deployment.

---

## 🔍 Troubleshooting

### Issue: "Connection refused"

**Solution:**
```bash
# Check if deployment is complete
gcloud run services describe mcp-perturbation \
  --region us-central1 \
  --project precision-medicine-poc
```

### Issue: "Tool execution failed"

**Possible causes:**
1. DRY_RUN mode enabled (expected for testing)
2. Missing data files
3. Resource limits exceeded

**Check logs:**
```bash
gcloud logging read \
  "resource.type=cloud_run_revision AND resource.labels.service_name=mcp-perturbation" \
  --limit 50 \
  --project precision-medicine-poc
```

### Issue: "Timeout after 60s"

**Solution:** Increase timeout (first request has cold start):
```bash
curl --max-time 120 ...
```

### Issue: Claude Desktop doesn't see server

**Solutions:**
1. Restart Claude Desktop after config change
2. Check config file syntax (valid JSON)
3. Verify command/args are correct
4. Check logs: `~/Library/Logs/Claude/mcp*.log`

---

## 📊 Expected Performance

| Operation | Local (stdio) | Remote (SSE) |
|-----------|---------------|--------------|
| Initialize | <1s | 1-3s (cold start: 5-15s) |
| Load dataset | 2-5s | 3-10s |
| Train model (20 epochs) | 5-10 min | 5-10 min |
| Predict response | 1-2s | 2-5s |

**Note:** First request to Cloud Run has cold start delay. Subsequent requests are fast.

---

## ✅ Deployment Verification Checklist

- [ ] Run `./test_deployment.sh` successfully
- [ ] All 8 tools listed:
  - [ ] perturbation_load_dataset
  - [ ] perturbation_setup_model
  - [ ] perturbation_train_model
  - [ ] perturbation_compute_delta
  - [ ] perturbation_predict_response
  - [ ] perturbation_differential_expression
  - [ ] perturbation_get_latent
  - [ ] perturbation_visualize
- [ ] No ERROR logs in Cloud Run
- [ ] Environment variables set correctly
- [ ] Resources allocated (4Gi memory, 2 CPU)
- [ ] SSE endpoint responds
- [ ] MCP handshake completes

---

## 🎯 Next Steps

1. **Immediate** (once deployment completes):
   - [ ] Run `./test_deployment.sh`
   - [ ] Check Cloud Run logs for errors

2. **Short-term** (within 1 hour):
   - [ ] Test with Claude Desktop
   - [ ] Run a full GEARS workflow
   - [ ] Verify prediction outputs

3. **Medium-term** (within 1 day):
   - [ ] Monitor GCP costs
   - [ ] Test with real patient data
   - [ ] Performance benchmarking

4. **Long-term** (within 1 week):
   - [ ] Consider GPU upgrade (see GPU_UPGRADE_PLAN.md)
   - [ ] Set up monitoring/alerting
   - [ ] Production deployment planning

---

**Created**: 2026-01-20
**Deployment**: GCP Cloud Run (us-central1)
**Status**: Ready for testing
