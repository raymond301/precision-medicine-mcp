# Project-Specific AI Skills

This directory contains **AI Skills** designed for the Precision Medicine MCP project. Each skill is a structured Markdown file that provides domain-specific guidance to AI coding assistants.

## What are AI Skills?

AI Skills are instruction sets (Markdown files) that help AI coding assistants (Claude Code, Cursor, Gemini CLI) follow project-specific patterns. They codify "tribal knowledge" about:
- Testing conventions (`.fn` pattern, DRY_RUN modes)
- Data handling (patient data paths, scRNA-seq preprocessing)
- Infrastructure deployment (Cloud Run, hospital environments)
- Clinical compliance (HIPAA, HUGO nomenclature)
- Documentation standards (persona-based writing)
- UI development (Streamlit, LLM providers)

## How to Use These Skills

Skills are **not** automatically detected by AI tools. You must explicitly reference them:

### In Claude Code
```
Read docs/reference/skills/precision-bio-tests/SKILL.md and use it to guide writing tests for my new server tool.
```

### In a prompt or chat
```
Follow the patterns in docs/reference/skills/patient-data-ops/SKILL.md when preprocessing this dataset.
```

### As project context
Reference this directory in your AI tool's project configuration (e.g., CLAUDE.md) so the assistant knows these skills exist.

## Available Skills

| Skill | Purpose | Key Topics |
| :--- | :--- | :--- |
| **precision-bio-tests** | Testing guidance | `.fn` pattern, DRY_RUN modes, PatientOne scenarios |
| **precision-bio-infra** | GCP deployment | Cloud Run, IAM, VPC, hospital deployment |
| **precision-bio-docs** | Documentation standards | Persona-based writing, Mermaid, server registry |
| **precision-bio-ui** | Frontend development | Streamlit, LLM providers, trace visualization |
| **clinical-compliance-doc** | Compliance & security | HIPAA, audit trails, HUGO nomenclature |
| **patient-data-ops** | Data management | Patient data paths, scRNA-seq, synthetic data |

## Project Context

This project has multiple custom MCP servers (plus a boilerplate template) in `servers/` — see the [Server Registry](../shared/server-registry.md) for current counts. Each server exposes bioinformatics tools via the FastMCP framework. All servers use `uv` for dependency management and support DRY_RUN mode for safe testing with synthetic data.

## Directory Structure

```
docs/reference/skills/
├── README.md                          # This file
├── precision-bio-tests/SKILL.md       # Testing guidance
├── precision-bio-infra/SKILL.md       # Infrastructure & deployment
├── precision-bio-docs/SKILL.md        # Documentation standards
├── precision-bio-ui/SKILL.md          # UI development
├── clinical-compliance-doc/SKILL.md   # Compliance & security
└── patient-data-ops/SKILL.md          # Data operations
```
