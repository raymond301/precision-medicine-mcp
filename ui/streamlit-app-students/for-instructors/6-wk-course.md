# Study Group Design: "AI-Powered Bioinformatics Workflows"

## Format Overview

- Duration: 3 months (6 sessions).
- Meetings: 1 hour, 2x/month.
- Homework: 2-4 hours between sessions.
- Focus: Practical application to real analysis workflows.

---

## Learning Objectives

By the end, participants will:
1. Use AI to accelerate common bioinformatics tasks.
2. Understand when MCP helps vs traditional scripting.
3. Create custom MCP tools for their specific needs.
4. Integrate AI into existing analysis pipelines.

---

## 6-Session Curriculum

### Month 1: Understanding the Value

#### Session 1: "AI as Your Bioinformatics Assistant" (Week 1)

**Goal:** Show immediate practical value. No setup required — students only use the deployed student app.

**Pre-session (Instructor):** Warm up servers 5 minutes before class:
```bash
curl -s https://mcp-spatialtools-ondu7mwjpa-uc.a.run.app/sse &
curl -s https://mcp-multiomics-ondu7mwjpa-uc.a.run.app/sse &
curl -s https://mcp-fgbio-ondu7mwjpa-uc.a.run.app/sse &
wait
```

**Meeting (60 min):**

1. Demo (15 min): Instructor demos the student app solving real problems
   - Open https://streamlit-mcp-chat-students-ondu7mwjpa-uc.a.run.app
   - Run "Warm Up Servers" — show tool discovery across 3 MCP servers
   - Run "Genomic QC" — FASTQ validation for PAT001 exome data from GCS
   - Run "Spatial Analysis" — Moran's I on Visium data for CD3D, CD8A, EPCAM, MKI67
   - Explain: each prompt calls real Python analysis code on GCS data via MCP

2. Hands-on (30 min): Everyone runs the 6 built-in prompts
   - Students open the student app URL in their browser (no install needed)
   - Run prompts in order: Warm Up Servers → Genomic QC → Spatial Analysis → Multi-omics Integration → Pathway Enrichment → PatientOne Mini Workflow
   - All 6 prompts use PatientOne (PAT001-OVC-2025) sample data from GCS
   - Note: if a prompt times out, click "Warm Up Servers" again and retry
   - Observe: which MCP server and tool was called for each prompt

3. Discussion (15 min): When would this help YOUR work?
   - Each person shares one repetitive analysis task they do
   - Group discusses which tasks could benefit from AI + MCP

**Homework (2-3 hours):**

Assignment 1: "Explore the Platform"

1. Open the student app (URL provided in session)
2. Run all 6 built-in prompts again — this time read the full responses carefully
3. Write 3 custom prompts related to YOUR work using the same GCS data paths:
   - `gs://sample-inputs-patientone/patient-data/PAT001-OVC-2025/spatial/` (Visium CSV files)
   - `gs://sample-inputs-patientone/patient-data/PAT001-OVC-2025/multiomics/` (RNA/Protein/Phospho CSVs)
   - `gs://sample-inputs-patientone/patient-data/PAT001-OVC-2025/genomics/fastq/` (FASTQ files)
4. Document:
   - Which built-in prompt was most relevant to your work?
   - What analysis do you do weekly that could be automated?
   - What data formats do you work with most?

**Deliverable:** 1-page document with your top 3 use cases

---

#### Session 2: "Understanding Your Data Through AI" (Week 3)

**Goal:** Look under the hood — understand how MCP servers work, clone the repo, explore the code.

**Meeting (60 min):**

1. Show & Tell (15 min): 3 participants share homework findings
   - Which custom prompts worked? Which didn't?
   - Surprises or errors encountered?

2. Deep Dive (25 min): How MCP tools actually work
   - Clone the repo: `git clone https://github.com/lynnlangit/precision-medicine-mcp`
   - Walk through `servers/mcp-spatialtools/` code together
   - Show: "This is just Python calling scanpy/squidpy — MCP wraps it with an AI interface"
   - Show: `server.py` tool definitions → how the student app calls them
   - Explain: the 3 active servers (spatialtools, multiomics, fgbio) each run on Cloud Run

3. Hands-on (20 min): Run deeper analyses
   - Use the student app with new custom prompts
   - Try different gene lists for spatial autocorrelation
   - Try pathway enrichment with your own gene sets
   - Discuss: what worked, what hit limits

**Homework (3-4 hours):**

Assignment 2: "Analyze Test Data"

1. Browse the PatientOne sample data in the repo:
   `tests/manual_testing/PatientOne-OvarianCancer/`

2. Using the student app, run 5 analyses:
   - Spatial autocorrelation (pick 3 different genes from the expression CSV)
   - Pathway enrichment (pick genes from your own research)
   - Multi-omics integration (observe sample counts and QC metrics)
   - Genomic QC (check FASTQ quality metrics)
   - PatientOne Mini Workflow (observe 2-step orchestration)

3. Compare AI results to what you'd do manually:
   - Did AI catch something you'd miss?
   - Any hallucinations or errors in the output?
   - Estimated time saved?

**Deliverable:** Jupyter notebook or document with your findings

---

### Month 2: Hands-On Integration

#### Session 3: "Integrating AI into Your Pipeline" (Week 5)

**Goal:** Understand the full architecture and design your own integration.

**Meeting (60 min):**

1. Case Study (20 min): Real pipeline integration
   - Show: Single-cell RNA-seq pipeline
   - Before: 10 manual steps in Jupyter
   - After: Natural language → automated workflow via MCP

2. Architecture Review (15 min): How the system is built
   - Streamlit app → LLM (Gemini) → MCP servers → bioinformatics tools
   - Show `servers/README.md` — how servers are deployed to Cloud Run
   - Show `ui/streamlit-app-students/` — how the student app is configured
   - Explain: you can add your own tools to existing servers

3. Group Activity (25 min): Design YOUR integration
   - Teams of 2-3
   - Pick one person's workflow from Session 1 homework
   - Sketch how MCP could automate it (which tools, what inputs/outputs)
   - Present to group (5 min each team)

**Homework (4 hours):**

Assignment 3: "Build a Mini-Pipeline"

Choose one:

Option A — For Computational Biologists:
1. Take your most common analysis (e.g., DE analysis)
2. Write a Python function that does it
3. Add it to mcp-spatialtools or create a new server
4. Test it locally with `python server.py`

Option B — For Wet-Lab Focused:
1. Document your standard QC workflow
2. Identify which steps could be automated
3. Write pseudocode for an MCP tool
4. Test existing student app tools on your QC data

**Deliverable:** Working code OR detailed specification + 5-minute demo video

---

#### Session 4: "Tool Show & Tell + Debugging" (Week 7)

**Goal:** Learn from each other, troubleshoot together.

**Meeting (60 min):**

1. Presentations (30 min): 3-4 people demo homework
   - What you built
   - What worked / what didn't
   - Questions for group

2. Debug Session (20 min): Fix issues together
   - Common problems:
     - Cloud Run cold starts (use warm-up prompt first)
     - File format mismatches (GCS paths must be exact)
     - Token limits (simplify prompts if responses truncate)
     - Memory issues with large datasets

3. Best Practices (10 min): Lessons learned
   - When to use AI vs traditional scripts
   - Data size limits for Cloud Run (300s timeout)
   - Cost considerations (Gemini flash pricing)

**Homework (3 hours):**

Assignment 4: "Comparative Analysis"

1. Pick ONE analysis task you do regularly
2. Do it THREE ways:
   a) Your traditional method (R/Python script)
   b) Using the student app AI interface
   c) Hybrid (script + AI for interpretation)

3. Measure and compare:
   - Time spent
   - Lines of code written
   - Accuracy/quality of results
   - Reproducibility

**Deliverable:** Comparison table + recommendation

---

### Month 3: Advanced Topics & Independence

#### Session 5: "Building Custom MCP Servers" (Week 9)

**Goal:** Create tools for your specific needs.

**Meeting (60 min):**

1. Advanced Demo (15 min): Custom server walkthrough
   - Show: `mcp-quantum-celltype-fidelity` server code
   - Explain: when to build custom vs use existing tools

2. Workshop (35 min): Build a simple MCP tool
   - Template provided (based on existing server structure)
   - Everyone writes a basic tool with `@mcp.tool` decorator
   - Test it locally with `python server.py`

3. Deployment Options (10 min): Where to run your server
   - Local only (for development)
   - Shared Cloud Run instance (for team use)
   - Your lab's infrastructure

**Homework (4-5 hours):**

Assignment 5: "Create Your Custom Tool"

Build an MCP server for your specific use case:

Ideas:
- Lab-specific QC tool
- Your favorite R package wrapped for MCP
- Integration with your LIMS system
- Custom visualization tool
- Your proprietary analysis method

Requirements:
1. Follows MCP protocol (uses `@mcp.tool` decorators)
2. Has 3+ tools/functions
3. Handles errors gracefully
4. Documented with example prompts
5. Works when tested locally

**Deliverable:** GitHub repo with your MCP server

---

#### Session 6: "Final Projects & Future Plans" (Week 11)

**Goal:** Showcase work, plan next steps.

**Meeting (60 min):**

1. Project Demos (40 min): Everyone presents
   - 5-7 minutes each
   - Show custom tool in action
   - Explain real-world impact

2. Group Discussion (15 min): What's next?
   - What worked well in the study group?
   - Should we continue? What topics?
   - Opportunities to contribute to the open-source repo

3. Resources & Next Steps (5 min)
   - Developer Streamlit app (all 14 prompts, Claude + Gemini)
   - Jupyter notebooks for deeper analysis
   - How to contribute to the repo

**Final Assignment (Optional):**

Capstone Project: "Production Deployment"

Deploy your custom tool to production:
1. Set up Cloud Run deployment with `deploy.sh`
2. Add your server to the Streamlit app config
3. Train 2 colleagues to use it
4. Gather feedback and iterate

**Deliverable:** Production tool + usage documentation
