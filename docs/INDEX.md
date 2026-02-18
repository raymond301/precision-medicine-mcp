# Documentation Index

**Welcome to the Precision Medicine MCP Documentation**

This is your central navigation hub for all documentation. Use this index to quickly find guides, references, and resources.

---

## 📋 Shared Reference Files

Canonical single-source-of-truth documents referenced across the repository:

- **[PatientOne Profile](reference/shared/patientone-profile.md)** - Clinical profile, genomic findings, data locations
- **[Platform Overview](reference/shared/README.md)** - Servers, tools, architecture summary
- **[Server Registry](reference/shared/server-registry.md)** - Server/tool counts, production readiness
- **[Value Proposition](reference/shared/value-proposition.md)** - 40 hours → 2-5 hours (production), cost savings, ROI
- **[Cost Analysis](reference/shared/cost-analysis.md)** - Per-patient costs, infrastructure, scaling projections
- **[HIPAA Summary](reference/shared/hipaa-summary.md)** - Compliance checklist, de-identification, retention
- **[Server Installation](reference/shared/server-installation.md)** - Standard Python MCP server setup
- **[DRY_RUN Mode](reference/shared/dry-run-mode.md)** - Mock mode explanation and per-server variables
- **[Deployment Templates](reference/shared/deployment-templates.md)** - Local, Streamlit Cloud, GCP Cloud Run

---

## 🚀 Getting Started

**New to the Precision Medicine MCP system? Start here:**

- **[Installation Guide](./getting-started/installation.md)** - Complete setup instructions (Quick Start: 5 minutes)
- **[Executive Summary](for-funders/EXECUTIVE_SUMMARY.md)** - High-level overview for decision-makers
- **[Why MCP for Healthcare?](reference/architecture/WHY_MCP_FOR_HEALTHCARE.md)** - Understand MCP architecture advantages

---

## 📚 User Guides

**Role-specific guides for different users:**

- **[For Funders & Grant Reviewers](./for-funders/README.md)** - ROI analysis, competitive landscape, grant talking points
- **[For Hospitals & IT Teams](./for-hospitals/README.md)** - Security, deployment, HIPAA compliance
- **[For Developers](./for-developers/README.md)** - Architecture, contributing, building new servers
- **[For Researchers](./for-researchers/README.md)** - Research workflows, analysis examples, production servers
- **[For Educators](./for-educators/README.md)** - Teaching materials, course integration, classroom activities
- **[For Patients & Families](./for-patients/README.md)** - Understanding precision medicine results (plain language)

**Specialized Guides:**
- **[Automated Patient Reports](./for-developers/automation-guides/AUTOMATED_PATIENT_REPORTS.md)** - Generate reports automatically
- **[Patient Summaries](./for-developers/automation-guides/GENERATE_PATIENT_SUMMARIES.md)** - Create patient-friendly summaries
- **[Add New Modality Server](./for-developers/ADD_NEW_MODALITY_SERVER.md)** - Extend the system with new data types

---

## 🏗️ Architecture

**Technical system design and implementation:**

- **[Architecture Overview](reference/architecture/README.md)** - System design principles
- **[Server Status](reference/shared/server-registry.md)** - Implementation status and capabilities matrix
- **[Clinical-Spatial Bridge](reference/architecture/clinical-spatial-bridge.md)** - Integration between clinical and spatial data
- **[Error Handling](reference/architecture/error-handling.md)** - Error handling and retry logic
- **[References](reference/architecture/references.md)** - Technical references and citations

---

## 📋 Operations

**Day-to-day operations and management:**

- **[Cost and Budget Management](./for-hospitals/operations/cost-and-budget.md)** - Cost estimation, tracking, and optimization
- **[Data Governance](./for-hospitals/compliance/data-governance.md)** - Data handling policies and procedures
- **[Live Monitoring Dashboard](../ui/dashboard/README.md)** - Real-time health monitoring for MCP servers + Streamlit clients, token usage tracking, cost optimization

---

## 🔒 Compliance

**Regulatory compliance and data security:**

- **[Compliance Overview](./for-hospitals/compliance/README.md)** - Compliance framework summary
- **[HIPAA Compliance](./for-hospitals/compliance/hipaa.md)** - De-identification, audit logging, encryption
- **[Data Governance](./for-hospitals/compliance/data-governance.md)** - GDPR, Common Rule, IRB requirements
- **[Risk Assessment](./for-hospitals/compliance/risk-assessment.md)** - Risk mitigation strategies
- **[Disclaimers](./for-hospitals/compliance/disclaimers.md)** - Legal disclaimers and limitations

---

## 🚀 Deployment

**Deployment guides and infrastructure:**

- **[Security Guide](reference/deployment/security.md)** - API keys, secrets management, GCP Secret Manager
- **[Hospital Deployment](./for-hospitals/)** - Enterprise deployment
  - [Operations Manual](./for-hospitals/OPERATIONS_MANUAL.md)
  - [Admin Guide](./for-hospitals/ADMIN_GUIDE.md)
  - [User Guide](./for-hospitals/USER_GUIDE.md)
  - [Audit Log Guide](./for-hospitals/AUDIT_LOG_GUIDE.md)
  - [Runbooks](./for-hospitals/RUNBOOKS/) - Incident response procedures

---

## 🧪 Testing

**Test documentation, strategies, and test data:**

- **[Test Documentation Index](reference/testing/README.md)** - Overview of all test documentation
- **[Test Coverage & Guidelines](reference/testing/test-coverage.md)** - Test structure and best practices
- **[Manual Testing](reference/testing)** - Quick test prompts and verification
- **[PatientOne Scenario](reference/testing/patient-one)** - Complete end-to-end testing scenario
  - [Quick Reference](reference/testing/patient-one/quick-reference.md)
  - [CITL Quick Test](reference/testing/patient-one/citl-quick-test.md)
  - [Test Prompts](reference/testing/patient-one/test-prompts) - Ready-to-use test prompts (6 tests)
- **[Integration Testing](reference/testing)** - GCP and API testing

---

## 🏥 Clinical Workflows

**Clinical decision support and review processes:**

- **[Clinical Overview](./for-hospitals/citl-workflows/)** - Clinical workflows overview
- **[CITL Workflow](./for-hospitals/citl-workflows/CITL_WORKFLOW_GUIDE.md)** - Clinician-in-the-Loop workflow
- **[CITL Review Template](./for-hospitals/citl-workflows/CITL_REVIEW_TEMPLATE.md)** - Review form template
- **[CITL Examples](./for-hospitals/citl-workflows/CITL_EXAMPLES.md)** - Example reviews

---

## ⚖️ Ethics

**Ethical considerations and bias auditing:**

- **[Ethics Overview](./for-hospitals/ethics/README.md)** - Ethics framework
- **[Bias Framework](./for-hospitals/ethics/ETHICS_AND_BIAS.md)** - Bias detection and mitigation
- **[Audit Checklist](./for-hospitals/ethics/BIAS_AUDIT_CHECKLIST.md)** - Step-by-step bias audit
- **[Implementation](./for-hospitals/ethics/IMPLEMENTATION_PLAN.md)** - Implementation guidelines

---

## 📦 Examples & Templates

**Example data, prompts, and templates:**

- **[Example Patients](./for-developers/automation-guides/examples/)** - Sample patient data and analyses
- **[Prompt Templates](./for-developers/automation-guides/prompts/)** - Reusable analysis prompts

---

## 📖 Additional Resources

### Quick Links

- **Installation:** [5-Minute Quick Start](./getting-started/installation.md#quick-start-5-minutes)
- **Test It:** [PatientOne Quick Test](reference/testing/patient-one/quick-reference.md)
- **Costs:** [Cost Overview](./for-hospitals/operations/cost-and-budget.md#cost-overview)
- **HIPAA:** [HIPAA Quick Reference](./for-hospitals/compliance/hipaa.md#executive-summary)

### By Task

**I want to:**
- **Install the system** → [Installation Guide](./getting-started/installation.md)
- **Run my first analysis** → [Quick Test Prompts](reference/testing/quick-test-prompts.md)
- **Understand costs** → [Cost and Budget Guide](./for-hospitals/operations/cost-and-budget.md)
- **Deploy to production** → [GCP Deployment Guide](reference/deployment/GCP_TESTING_GUIDE.md)
- **Ensure HIPAA compliance** → [HIPAA Compliance](./for-hospitals/compliance/hipaa.md)
- **Add a new server** → [Add New Modality Server](./for-developers/ADD_NEW_MODALITY_SERVER.md)
- **Review patient results** → [CITL Workflow](./for-hospitals/citl-workflows/CITL_WORKFLOW_GUIDE.md)
- **Test the system** → [PatientOne Scenario](reference/testing/patient-one/README.md)

### By Role

- **Funder/Grant Reviewer** → [For Funders](./for-funders/README.md) - ROI, competitive landscape, grant materials
- **Hospital Administrator** → [For Hospitals](./for-hospitals/README.md) - Security, deployment checklist, HIPAA
- **Developer** → [For Developers](./for-developers/README.md) - Architecture, contributing, building servers
- **Researcher/Bioinformatician** → [For Researchers](./for-researchers/README.md) - Workflows, analysis, production servers
- **Educator/Professor** → [For Educators](./for-educators/README.md) - Teaching materials, course integration
- **Patient/Family** → [For Patients](./for-patients/README.md) - Understanding results (plain language)
- **Clinician** → [CITL Workflow](./for-hospitals/citl-workflows/CITL_WORKFLOW_GUIDE.md) - Clinical workflow integration

---

## 🗂️ Repository Structure

Complete directory structure showing all major components:

```
precision-medicine-mcp/
├── ACKNOWLEDGMENTS.md      # Credits & scientific references
├── CLAUDE.md               # Claude Code project context
├── LICENSE                 # Apache 2.0 License
├── README.md               # Main repository README
├── llms.txt                # LLM-readable project summary
├── data/                   # Synthetic patient data (100% safe for demos)
├── docs/                   # Documentation organized by audience
│   ├── for-funders/        # ROI, competitive landscape, grant talking points, demos
│   ├── for-hospitals/      # Deployment checklist, security overview, operations
│   ├── for-developers/     # Architecture, contributing guide, quick reference
│   ├── for-researchers/    # Analysis workflows, bioinformatics methods
│   ├── for-educators/      # Classroom guides, learning objectives
│   ├── for-operations/     # SLA, implementation guide, operations
│   ├── for-patients/       # Patient-friendly resources
│   ├── getting-started/    # Installation, quick start, desktop-configs
│   ├── book/               # Quarto book: AI-Orchestrated Precision Oncology
│   └── reference/          # Technical reference (architecture, deployment, prompts, tests)
│       ├── architecture/   # System design & modality workflows
│       ├── deployment/     # GCP deployment status & guides
│       ├── prompts/ # 20+ ready-to-use clinical prompts
│       ├── shared/         # Canonical single-source-of-truth files
│       └── testing/      # Testing guides & PatientOne scenarios
├── infrastructure/         # Deployment, audit, environment setup
│   ├── deployment/         # GCP deployment scripts
│   ├── docker/             # Base Docker images for Cloud Run
│   ├── audit/              # Bias detection and audit tools
│   └── hospital-deployment/  # Hospital-specific infrastructure
├── servers/                # MCP servers + boilerplate template (Python)
│   ├── mcp-cell-classify/  # Cell phenotype classification
│   ├── mcp-deepcell/       # Cell segmentation
│   ├── mcp-epic/           # Epic FHIR integration
│   ├── mcp-fgbio/          # Reference genomes, FASTQ QC
│   ├── mcp-genomic-results/ # Somatic variant/CNV/HRD
│   ├── mcp-huggingface/    # AI/ML inference
│   ├── mcp-mockepic/       # Mock Epic for testing
│   ├── mcp-multiomics/     # Multi-omics integration
│   ├── mcp-openimagedata/  # Imaging data (H&E, MxIF)
│   ├── mcp-patient-report/ # Patient-facing reports
│   ├── mcp-perturbation/   # GEARS treatment prediction
│   ├── mcp-quantum-celltype-fidelity/ # Quantum fidelity
│   ├── mcp-seqera/         # Workflow orchestration
│   ├── mcp-server-boilerplate/ # Template for new servers
│   ├── mcp-spatialtools/   # Spatial transcriptomics
│   └── mcp-tcga/           # TCGA cohort data
├── shared/                 # Shared Python packages
│   ├── common/             # Common utilities
│   ├── models/             # Data models
│   ├── schemas/            # JSON schemas (CitL review, etc.)
│   └── utils/              # Helper functions
├── results/                # Analysis output files
├── tests/                  # Automated tests (unit, integration, verification)
└── ui/                     # User interfaces
    ├── streamlit-app/      # Main Streamlit chat interface
    ├── streamlit-app-students/ # Student/classroom version
    ├── dashboard/          # Live monitoring dashboard
    └── jupyter-notebook/   # Jupyter notebook interface
```

---

