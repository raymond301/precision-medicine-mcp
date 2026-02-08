# Inventory: Software Bill of Materials (SBOM)

[← Back to SLA Overview](README.md)

---

## 1. Core Platform Components

| Component | Version | Source | HIPAA Status |
| :--- | :--- | :--- | :--- |
| **Precision Medicine Framework** | v2.0.1 | Internal Git | ✅ Verified |
| **Python** | 3.11-slim | Docker Hub | ✅ Verified |
| **fgbio** | 2.0.1 | Bioconda | ✅ Verified |
| **samtools** | 1.17 | Bioconda | ✅ Verified |
| **pysam** | 0.21.0 | PyPI | ✅ Verified |

## 2. Infrastructure Layer

| Tool | Version | Purpose |
| :--- | :--- | :--- |
| **Terraform** | 1.5+ | Infrastructure as Code |
| **gcloud SDK** | 440+ | GCP CLI Management |
| **Docker** | 24.0+ | Container Engine |

## 3. Identity & Orchestration

| Component | Version | Purpose |
| :--- | :--- | :--- |
| **Azure AD SSO** | Enterprise | Identity Provider |
| **Anthropic API** | Claude 3.5 Sonnet | LLM Orchestration |
| **OAuth2 Proxy** | 7.4.0 | Authentication Proxy |

---

### Document History

| Date | Version | Author | Change Summary |
|:---|:---|:---|:---|
| 2026-02-08 | 2.0 | IT Operations | Initial creation |
