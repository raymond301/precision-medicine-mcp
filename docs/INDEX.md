# Documentation Index

**Welcome to the Precision Medicine MCP Documentation**

This is your central navigation hub for all documentation. Use this index to quickly find guides, references, and resources.

---

## 🚀 Getting Started

**New to the Precision Medicine MCP system? Start here:**

- **[Installation Guide](./getting-started/installation.md)** - Complete setup instructions (Quick Start: 5 minutes)
- **[Executive Summary](./EXECUTIVE_SUMMARY.md)** - High-level overview for decision-makers
- **[Why MCP for Healthcare?](./WHY_MCP_FOR_HEALTHCARE.md)** - Understand MCP architecture advantages

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
- **[Server Status](./SERVER_REGISTRY.md)** - Implementation status and capabilities matrix
- **[Clinical-Spatial Bridge](reference/architecture/clinical-spatial-bridge.md)** - Integration between clinical and spatial data
- **[Error Handling](reference/architecture/error-handling.md)** - Error handling and retry logic
- **[References](reference/architecture/references.md)** - Technical references and citations

---

## 📋 Operations

**Day-to-day operations and management:**

- **[Cost and Budget Management](./for-hospitals/operations/cost-and-budget.md)** - Cost estimation, tracking, and optimization
- **[Data Governance](./for-hospitals/compliance/data-governance.md)** - Data handling policies and procedures

---

## 🔒 Compliance

**Regulatory compliance and data security:**

- **[Compliance Overview](./for-hospitals/compliance/README.md)** - Compliance framework summary
- **[HIPAA Compliance](./for-hospitals/compliance/hipaa.md)** - De-identification, audit logging, encryption
- **[Data Governance](./for-hospitals/compliance/data-governance.md)** - GDPR, Common Rule, IRB requirements
- **[Risk Assessment](./for-hospitals/compliance/risk-assessment.md)** - Risk mitigation strategies (if exists)
- **[Disclaimers](./for-hospitals/compliance/disclaimers.md)** - Legal disclaimers and limitations

---

## 🚀 Deployment

**Deployment guides and infrastructure:**

- **[Deployment Roadmap](./reference/archive/deployment/roadmap.md)** - Production deployment planning
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

- **[Test Documentation Index](reference/test-docs/README.md)** - Overview of all test documentation
- **[Test Coverage & Guidelines](reference/test-docs/test-coverage.md)** - Test structure and best practices
- **[Manual Testing](reference/test-docs/manual-testing)** - Quick test prompts and verification
- **[PatientOne Scenario](reference/test-docs/patient-one-scenario)** - Complete end-to-end testing scenario
  - [Quick Reference](reference/test-docs/patient-one-scenario/quick-reference.md)
  - [CITL Quick Test](reference/test-docs/patient-one-scenario/citl-quick-test.md)
  - [Test Prompts](reference/test-docs/patient-one-scenario/test-prompts) - Ready-to-use test prompts (6 tests)
- **[Integration Testing](reference/test-docs/integration-testing)** - GCP and API testing

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

## 🗄️ Archive

**Outdated documentation preserved for historical reference:**

- **[Archive Index](reference/archive/README.md)** - Why docs were archived
- **[2025 Q3/Q4 Archive](reference/archive/2025-q3-q4)** - Pre-October 2025 docs

⚠️ **Warning:** Archived documentation is outdated and should not be used for current implementations.

---

## 📖 Additional Resources

### Quick Links

- **Installation:** [5-Minute Quick Start](./getting-started/installation.md#quick-start-5-minutes)
- **Test It:** [PatientOne Quick Test](reference/test-docs/patient-one-scenario/quick-reference.md)
- **Costs:** [Cost Overview](./for-hospitals/operations/cost-and-budget.md#cost-overview)
- **HIPAA:** [HIPAA Quick Reference](./for-hospitals/compliance/hipaa.md#executive-summary)

### By Task

**I want to:**
- **Install the system** → [Installation Guide](./getting-started/installation.md)
- **Run my first analysis** → [Quick Test Prompts](reference/test-docs/manual-testing/quick-test-prompts.md)
- **Understand costs** → [Cost and Budget Guide](./for-hospitals/operations/cost-and-budget.md)
- **Deploy to production** → [Deployment Roadmap](./reference/archive/deployment/roadmap.md)
- **Ensure HIPAA compliance** → [HIPAA Compliance](./for-hospitals/compliance/hipaa.md)
- **Add a new server** → [Add New Modality Server](./for-developers/ADD_NEW_MODALITY_SERVER.md)
- **Review patient results** → [CITL Workflow](./for-hospitals/citl-workflows/CITL_WORKFLOW_GUIDE.md)
- **Test the system** → [PatientOne Scenario](reference/test-docs/patient-one-scenario/README.md)

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
├── LICENSE                 # Apache 2.0 License
├── README.md               # Main repository README
├── data/                   # Synthetic patient data (100% safe for demos)
├── docs/                   # Documentation organized by audience
│   ├── for-funders/        # ROI, competitive landscape, grant talking points, demos
│   ├── for-hospitals/      # Deployment checklist, security overview, operations
│   ├── for-developers/     # Architecture, contributing guide, quick reference
│   ├── for-researchers/    # Analysis workflows, bioinformatics methods
│   ├── for-educators/      # Classroom guides, learning objectives
│   ├── for-patients/       # Patient-friendly resources
│   ├── getting-started/    # Installation, quick start, desktop-configs
│   ├── book/               # Quarto book: AI-Orchestrated Precision Oncology
│   └── reference/          # Technical reference (architecture, deployment, prompts, tests)
│       ├── architecture/   # System design & modality workflows
│       ├── deployment/     # GCP deployment status & guides
│       ├── prompt-library/ # 20+ ready-to-use clinical prompts
│       ├── test-docs/      # Testing guides & PatientOne scenarios
│       └── archive/        # Historical documentation
├── infrastructure/         # Deployment, audit, environment setup
│   ├── deployment/         # GCP deployment scripts
│   ├── audit/              # Bias detection and audit tools
│   └── hospital-deployment/  # Hospital-specific infrastructure
├── servers/                # 10 MCP servers (Python)
│   ├── mcp-deepcell/       # Cell segmentation
│   ├── mcp-epic/           # Epic FHIR integration
│   ├── mcp-fgbio/          # Reference genomes, FASTQ QC
│   ├── mcp-huggingface/    # AI/ML inference
│   ├── mcp-mockepic/       # Mock Epic for testing
│   ├── mcp-multiomics/     # Multi-omics integration
│   ├── mcp-openimagedata/  # Imaging data (H&E, MxIF)
│   ├── mcp-seqera/         # Workflow orchestration
│   ├── mcp-spatialtools/   # Spatial transcriptomics
│   └── mcp-tcga/           # TCGA cohort data
├── shared/                 # Shared Python packages
│   ├── common/             # Common utilities
│   ├── models/             # Data models
│   ├── schemas/            # JSON schemas (CitL review, etc.)
│   └── utils/              # Helper functions
├── tests/                  # 167 automated tests
├── tools/                  # Automation & reporting tools
│   └── reports/            # Patient report generation, CitL submission
└── ui/                     # Streamlit chat, Jupyter notebook
```

---

