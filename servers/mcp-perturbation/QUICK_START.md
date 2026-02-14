# mcp-perturbation Quick Start

**Status**: ✅ **DEPLOYED AND READY**
**Service URL**: https://mcp-perturbation-ondu7mwjpa-uc.a.run.app/sse

---

## 🚀 Start Using in 2 Minutes

### For Claude Desktop Users

1. **Open your Claude Desktop config file**:
   ```bash
   open ~/Library/Application\ Support/Claude/claude_desktop_config.json
   ```

2. **Add the mcp-perturbation server**:
   ```json
   {
     "mcpServers": {
       "perturbation": {
         "command": "npx",
         "args": [
           "-y",
           "@modelcontextprotocol/client-sse",
           "https://mcp-perturbation-ondu7mwjpa-uc.a.run.app/sse"
         ]
       }
     }
   }
   ```

3. **Restart Claude Desktop**

4. **Test it**:
   - Open a new conversation
   - Type: *"List the available perturbation tools"*
   - You should see 8 GEARS tools listed!

---

## 🧬 Example: Predict Treatment Response

Once connected, try this workflow:

```
User: "Load the GSE184880 ovarian cancer dataset using the perturbation server"

[Server loads dataset with ~12,000 cells from ovarian cancer patients]

User: "Setup a GEARS model for predicting checkpoint inhibitor response"

[Server initializes graph neural network]

User: "Train the model for 20 epochs"

[Server trains GEARS GNN, takes ~5-10 minutes]

User: "Predict how T cells from patient PAT001-OVC-2025 would respond to anti-PD1 therapy"

[Server returns predicted gene expression changes and top affected genes like GZMB↑, PRF1↑, IFNG↑]
```

---

## 🛠️ Available Tools

1. **perturbation_load_dataset** - Load scRNA-seq from GEO or .h5ad file
2. **perturbation_setup_model** - Initialize GEARS graph neural network
3. **perturbation_train_model** - Train model on reference data
4. **perturbation_compute_delta** - Calculate perturbation effect vector
5. **perturbation_predict_response** - Predict treatment response for patient cells
6. **perturbation_differential_expression** - Find top changed genes
7. **perturbation_get_latent** - Extract GNN embeddings
8. **perturbation_visualize** - Generate PCA/UMAP plots

---

## 📊 What It Does

The mcp-perturbation server uses **GEARS** (Graph-Enhanced Gene Activation and Repression Simulator), a state-of-the-art graph neural network from Nature Biotechnology 2024, to:

- Predict how patient cells will respond to treatments **without experiments**
- Model multi-gene perturbations (e.g., checkpoint inhibitor combinations)
- Identify which genes change in response to therapy
- Help prioritize which treatments to test in the clinic

**Key advantages**:
- 40% more accurate than older VAE-based methods
- Handles combinatorial drug effects
- Integrates biological knowledge (gene networks)
- Fast predictions (~2-5 seconds per patient)

---

## 💡 Use Cases

### 1. **Immunotherapy Response Prediction**
Predict if a patient's T cells will respond to checkpoint inhibitors (anti-PD1, anti-CTLA4)

### 2. **Drug Screening**
Test multiple treatments *in silico* before expensive experiments

### 3. **Personalized Medicine**
Identify optimal therapies based on patient-specific cell profiles

### 4. **Clinical Trial Design**
Pre-screen patients likely to respond to experimental therapies

---

## 🧪 Test Data Available

- **GSE184880**: Ovarian cancer (HGSOC) scRNA-seq (5 controls, 7 patients)
- **PatientOne**: Synthetic ovarian cancer patient (500 T cells)
- **GEARS datasets**: Norman, Adamson, Dixit benchmark datasets

You can also upload your own .h5ad files to Google Cloud Storage and reference them by GCS path.

---

## 💰 Cost

- **Per prediction**: ~$0.02 (10 min training + prediction)
- **Idle time**: $0 (scales to zero)
- **Always-on**: ~$65/month (not recommended)

The server automatically scales down when idle, so you only pay for actual usage.

---

## ❓ Troubleshooting

### Issue: Claude Desktop doesn't see the server

**Solution**:
1. Check the config file syntax (must be valid JSON)
2. Make sure you restarted Claude Desktop completely
3. Check logs: `~/Library/Logs/Claude/mcp*.log`

### Issue: "Server connection failed"

**Solution**:
1. Verify you have internet connection
2. Check the URL is correct: `https://mcp-perturbation-ondu7mwjpa-uc.a.run.app/sse`
3. First request may take 5-15 seconds (cold start)

### Issue: Tools are listed but execution fails

**Possible causes**:
1. Dataset not found (check GEO ID or GCS path)
2. Model not trained yet (must train before predicting)
3. Invalid parameters (check tool documentation)

---

## 📚 More Information

- [README.md](README.md) - Comprehensive documentation
- [ARCHITECTURE.md](ARCHITECTURE.md) - System architecture with diagrams
- [GCP Deployment Guide](../../docs/reference/deployment/GCP_TESTING_GUIDE.md) - Deployment and testing

---

## 🎯 Next Steps

1. ✅ Add to Claude Desktop (see above)
2. ✅ Test with "List perturbation tools"
3. ✅ Try the example workflow
4. ✅ Predict treatment response for your own patient data
5. ✅ Monitor costs in GCP Console

---

**Ready to predict treatment responses? Add it to Claude Desktop now!** 🚀

**Service URL**: https://mcp-perturbation-ondu7mwjpa-uc.a.run.app/sse
