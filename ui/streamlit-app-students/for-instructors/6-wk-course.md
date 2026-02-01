 # 📚 Study Group Design: "AI-Powered Bioinformatics Workflows"

  ## Format Overview

  - Duration: 3 months (6 sessions). 
  - Meetings: 1 hour, 2x/month. 
  - Homework: 2-4 hours between sessions. 
  - Focus: Practical application to real analysis workflows. 

  ---
  
  ## 🎯 Learning Objectives

  By the end, participants will:  
  1. Use AI to accelerate common bioinformatics tasks. 
  2. Understand when MCP helps vs traditional scripting. 
  3. Create custom MCP tools for their specific needs. 
  4. Integrate AI into existing analysis pipelines. 

  ---

  ## 📅 6-Session Curriculum

  ### Month 1: Understanding the Value. 

  #### Session 1: "AI as Your Bioinformatics Assistant" (Week 1)  

  #### Goal: Show immediate practical value. 

  #### Meeting (60 min):
  1. Demo (15 min): Live demo of spatial-mcp solving real problems
    - "Calculate spatial autocorrelation for marker genes"
    - "Run multi-omics integration on my data"
    - "Validate my FASTQ files"
    - Show how it's faster than writing scripts
  2. Hands-on (30 min): Everyone runs the Streamlit app
    - Use mock mode (no infrastructure setup needed!)
    - Try 3-5 pre-written queries
    - See instant "analysis" results
  3. Discussion (15 min): When would this help YOUR work?
    - Each person shares one repetitive task they do
    - Group discusses if MCP could help

  #### Homework (2-3 hours):
  Assignment 1: "Find Your Pain Points"

  1. Install the repo locally (DEVELOPMENT.md guide)
  2. Run app in mock mode (USE_MOCK_MCP=true)
  3. Try 10 different queries related to YOUR work
  4. Document:
     - Which 3 queries gave useful results?
     - What analysis do you do weekly that could be automated?
     - What data formats do you work with most?

  #### Deliverable: 1-page document with your top 3 use cases

  ---
  #### Session 2: "Understanding Your Data Through AI" (Week 3)

  #### Goal: Use AI for exploratory data analysis

  #### Meeting (60 min):
  1. Show & Tell (20 min): 3 participants share homework findings
  2. Deep Dive (25 min): How MCP tools actually work
    - Walk through servers/mcp-spatialtools/ code
    - Show: "This is just Python calling scanpy/squidpy"
    - Explain: MCP wraps existing tools, adds AI interface
  3. Hands-on (15 min): Upload YOUR data
    - Use GCS file upload feature
    - Try analysis on real data (if small dataset)
    - Or use provided PatientOne test data

  Homework (3-4 hours):
  Assignment 2: "Analyze Test Data"

  1. Download PatientOne ovarian cancer dataset:
     tests/manual_testing/PatientOne-OvarianCancer/

  2. Using the Streamlit app, run 5 analyses:
     - Spatial autocorrelation (pick 3 genes)
     - Cell type deconvolution
     - Pathway enrichment
     - Multi-omics correlation
     - Your choice

  3. Compare AI results to what you'd do manually:
     - Did AI catch something you'd miss?
     - Any hallucinations/errors?
     - Time saved estimate

  Deliverable: Jupyter notebook with your findings

  ---
  ### Month 2: Hands-On Integration

  Session 3: "Integrating AI into Your Pipeline" (Week 5)

  Goal: Add MCP to existing workflows

  Meeting (60 min):
  1. Case Study (20 min): Real pipeline integration
    - Show: Single-cell RNA-seq pipeline
    - Before: 10 manual steps in Jupyter
    - After: Natural language → automated workflow
  2. Architecture Review (15 min): How it works
    - Streamlit app → MCP servers → Your tools
    - Show server deployment (servers/README.md)
    - Explain: You can add your own tools
  3. Group Activity (25 min): Design YOUR integration
    - Teams of 2-3
    - Pick one person's workflow
    - Sketch how MCP could help
    - Present to group (5 min each team)

  Homework (4 hours):
  Assignment 3: "Build a Mini-Pipeline"

  Choose one:

  Option A - For Computational Biologists:
  1. Take your most common analysis (e.g., DE analysis)
  2. Write a Python function that does it
  3. Add it to mcp-spatialtools or create new server
  4. Test it via Streamlit app

  Option B - For Wet-Lab Focused:
  1. Document your standard QC workflow
  2. Identify which steps could be automated
  3. Write pseudocode for an MCP tool
  4. Test existing tools on your QC data

  Deliverable:
  - Working code OR detailed specification
  - 5-minute demo video for next session

  ---
  Session 4: "Tool Show & Tell + Debugging" (Week 7)

  Goal: Learn from each other, troubleshoot together

  Meeting (60 min):
  1. Presentations (30 min): 3-4 people demo homework
    - What you built
    - What worked / what didn't
    - Questions for group
  2. Debug Session (20 min): Fix issues together
    - Common problems:
        - File format mismatches
      - Memory issues with large datasets
      - API rate limits
  3. Best Practices (10 min): Lessons learned
    - When to use AI vs traditional scripts
    - Data size limits
    - Cost considerations

  Homework (3 hours):
  Assignment 4: "Comparative Analysis"

  1. Pick ONE analysis task you do regularly
  2. Do it THREE ways:
     a) Your traditional method (R/Python script)
     b) Using spatial-mcp AI interface
     c) Hybrid (script + AI for interpretation)

  3. Measure and compare:
     - Time spent
     - Lines of code written
     - Accuracy/quality
     - Reproducibility

  Deliverable: Comparison table + recommendation

  ---
  ### Month 3: Advanced Topics & Independence

  Session 5: "Building Custom MCP Servers" (Week 9)

  Goal: Create tools for your specific needs

  Meeting (60 min):
  1. Advanced Demo (15 min): Custom server walkthrough
    - Show: mcp-quantum-celltype-fidelity (your newest server)
    - Explain: When to build custom vs use existing
  2. Workshop (35 min): Build a simple MCP tool
    - Template provided
    - Everyone writes a basic tool
    - Test it locally (mock mode)
  3. Deployment Options (10 min): Where to run
    - Local only
    - Shared Cloud Run instance
    - Your lab's infrastructure

  Homework (4-5 hours):
  Assignment 5: "Create Your Custom Tool"

  Build an MCP server for your specific use case:

  Ideas:
  - Lab-specific QC tool
  - Your favorite R package wrapped for MCP
  - Integration with your LIMS system
  - Custom visualization tool
  - Your proprietary analysis method

  Requirements:
  1. Follows MCP protocol
  2. Works with Streamlit UI
  3. Handles errors gracefully
  4. Has 3+ tools/functions
  5. Documented with examples

  Deliverable: GitHub repo with your MCP server

  ---
  Session 6: "Final Projects & Future Plans" (Week 11)

  Goal: Showcase work, plan next steps

  Meeting (60 min):
  1. Project Demos (40 min): Everyone presents
    - 5-7 minutes each
    - Show custom tool in action
    - Explain real-world impact
  2. Group Discussion (15 min): What's next?
    - What worked well in study group?
    - Should we continue?
    - What topics for next session?
  3. Resources & Next Steps (5 min)
    - Advanced topics to explore
    - How to contribute to repo
    - Community resources

  Final Assignment (Optional):
  Capstone Project: "Production Deployment"

  Deploy your custom tool to production:
  1. Set up Cloud Run deployment
  2. Add to your lab's Streamlit instance
  3. Train 2 colleagues to use it
  4. Gather feedback and iterate

  Deliverable: Production tool + usage documentation