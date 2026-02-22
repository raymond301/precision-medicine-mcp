"""MCP Server Configuration

Contains URLs and metadata for MCP servers configured in the Streamlit app
(excludes mcp-epic which is local-only and mcp-genomic-results).

Source of truth for server/tool counts: docs/reference/shared/server-registry.md
"""

from typing import Dict, List, TypedDict


class MCPServerConfig(TypedDict):
    """MCP server configuration."""
    name: str
    url: str
    description: str
    status: str
    tools_count: int


# MCP Server URLs - Set to local for testing, Cloud Run for production
# LOCAL: http://localhost:PORT/sse (requires ngrok/tunnel for Claude API to reach)
# CLOUD RUN: https://mcp-SERVER-ondu7mwjpa-uc.a.run.app/sse
USE_LOCAL_SERVERS = False  # Must be False when using Claude API beta MCP client

MCP_SERVERS: Dict[str, MCPServerConfig] = {
    "fgbio": {
        "name": "fgbio",
        "url": "http://localhost:8001/sse" if USE_LOCAL_SERVERS else "https://mcp-fgbio-ondu7mwjpa-uc.a.run.app/sse",
        "description": "Genomic reference data and FASTQ validation",
        "status": "production",
        "tools_count": 4
    },
    "multiomics": {
        "name": "multiomics",
        "url": "http://localhost:8002/sse" if USE_LOCAL_SERVERS else "https://mcp-multiomics-ondu7mwjpa-uc.a.run.app/sse",
        "description": "Multi-omics integration (RNA/Protein/Phospho)",
        "status": "production",
        "tools_count": 10
    },
    "spatialtools": {
        "name": "spatialtools",
        "url": "http://localhost:8000/sse" if USE_LOCAL_SERVERS else "https://mcp-spatialtools-ondu7mwjpa-uc.a.run.app/sse",
        "description": "Spatial transcriptomics analysis",
        "status": "production",
        "tools_count": 14
    },
    "tcga": {
        "name": "tcga",
        "url": "https://mcp-tcga-ondu7mwjpa-uc.a.run.app/sse",
        "description": "TCGA cancer genomics data",
        "status": "mock",
        "tools_count": 5
    },
    "openimagedata": {
        "name": "openimagedata",
        "url": "https://mcp-openimagedata-ondu7mwjpa-uc.a.run.app/sse",
        "description": "H&E/MxIF image loading, registration, feature extraction, and composite generation",
        "status": "production",
        "tools_count": 5
    },
    "deepcell": {
        "name": "deepcell",
        "url": "https://mcp-deepcell-ondu7mwjpa-uc.a.run.app/sse",
        "description": "DeepCell-TF cell segmentation and marker quantification for MxIF",
        "status": "production",
        "tools_count": 3
    },
    "cell-classify": {
        "name": "cell-classify",
        "url": "https://mcp-cell-classify-ondu7mwjpa-uc.a.run.app/sse",
        "description": "Cell phenotype classification and visualization (lightweight, no TensorFlow)",
        "status": "production",
        "tools_count": 3
    },
    "mockepic": {
        "name": "mockepic",
        "url": "https://mcp-mockepic-ondu7mwjpa-uc.a.run.app/sse",
        "description": "Mock EHR/FHIR data",
        "status": "mock",
        "tools_count": 3
    },
    "perturbation": {
        "name": "perturbation",
        "url": "http://localhost:8080/sse" if USE_LOCAL_SERVERS else "https://mcp-perturbation-ondu7mwjpa-uc.a.run.app/sse",
        "description": "GEARS perturbation prediction for treatment response",
        "status": "production",
        "tools_count": 8
    },
    "quantum-celltype-fidelity": {
        "name": "quantum-celltype-fidelity",
        "url": "http://localhost:3010/sse" if USE_LOCAL_SERVERS else "https://mcp-quantum-celltype-fidelity-ondu7mwjpa-uc.a.run.app/sse",
        "description": "Quantum computing for cell type validation and immune evasion detection",
        "status": "production",
        "tools_count": 6
    },
    "patient-report": {
        "name": "patient-report",
        "url": "http://localhost:3011/sse" if USE_LOCAL_SERVERS else "https://mcp-patient-report-ondu7mwjpa-uc.a.run.app/sse",
        "description": "Patient-facing PDF reports with plain-language summaries",
        "status": "production",
        "tools_count": 5
    }
}


def get_server_config(server_names: List[str]) -> List[Dict]:
    """Get MCP server configs for Claude API.

    Args:
        server_names: List of server names to include

    Returns:
        List of MCP server configurations for Claude API
    """
    configs = []
    for name in server_names:
        if name in MCP_SERVERS:
            configs.append({
                "type": "url",
                "url": MCP_SERVERS[name]["url"],
                "name": name
            })
    return configs


def get_tools_config(server_names: List[str]) -> List[Dict]:
    """Get tools configuration for Claude API.

    Args:
        server_names: List of server names to include

    Returns:
        List of tool configurations for Claude API
    """
    return [
        {"type": "mcp_toolset", "mcp_server_name": name}
        for name in server_names
        if name in MCP_SERVERS
    ]


def get_server_categories() -> Dict[str, List[str]]:
    """Get servers organized by category.

    Returns:
        Dict of category -> list of server names
    """
    return {
        "Production Servers (Real Analysis)": [
            "fgbio",
            "multiomics",
            "spatialtools",
            "perturbation",
            "quantum-celltype-fidelity",
            "deepcell",
            "cell-classify",
            "openimagedata",
            "patient-report"
        ],
        "Mock Servers (Workflow Demo)": [
            "tcga",
            "mockepic"
        ]
    }


# Example prompts for different use cases
EXAMPLE_PROMPTS = {
    "Warm Up Servers": "List all available tools from the connected MCP servers. For each tool, show its name and a one-line description.",

    "Spatial Analysis": "Use the spatial_autocorrelation tool with expression_file=gs://sample-inputs-patientone/patient-data/PAT001-OVC-2025/spatial/visium_gene_expression.csv and coordinates_file=gs://sample-inputs-patientone/patient-data/PAT001-OVC-2025/spatial/visium_spatial_coordinates.csv to calculate spatial autocorrelation for genes CD3D, CD8A, EPCAM, MKI67.",

    "Multi-omics Integration": "Use the integrate_omics_data tool with rna_path=gs://sample-inputs-patientone/patient-data/PAT001-OVC-2025/multiomics/pdx_rna_seq.csv, protein_path=gs://sample-inputs-patientone/patient-data/PAT001-OVC-2025/multiomics/pdx_proteomics.csv, phospho_path=gs://sample-inputs-patientone/patient-data/PAT001-OVC-2025/multiomics/pdx_phosphoproteomics.csv. Report the number of common samples, features per modality, and QC metrics.",

    "Genomic QC": "Validate the FASTQ file for Patient-001 (PAT001-OVC-2025) at gs://sample-inputs-patientone/patient-data/PAT001-OVC-2025/genomics/fastq/PAT001_OVC_exome_R1.fastq.gz and check quality metrics. What is the average quality score and read length?",

    "Pathway Enrichment": "For the upregulated genes [TP53, BRCA1, MYC, KRAS], perform pathway enrichment analysis using GO_BP database.",

    "Complete PatientOne Workflow": """For Patient-001 (PAT001-OVC-2025, ovarian cancer). Sample data is in GCS bucket gs://sample-inputs-patientone/patient-data/PAT001-OVC-2025/:
1. Get clinical data from FHIR using query_patient_records with patient_id=PAT001-OVC-2025
2. Use spatial_autocorrelation with expression_file=gs://sample-inputs-patientone/patient-data/PAT001-OVC-2025/spatial/visium_gene_expression.csv and coordinates_file=gs://sample-inputs-patientone/patient-data/PAT001-OVC-2025/spatial/visium_spatial_coordinates.csv for genes CD3D, CD8A, EPCAM, MKI67
3. Use integrate_omics_data with rna_path=gs://sample-inputs-patientone/patient-data/PAT001-OVC-2025/multiomics/pdx_rna_seq.csv, protein_path=gs://sample-inputs-patientone/patient-data/PAT001-OVC-2025/multiomics/pdx_proteomics.csv, phospho_path=gs://sample-inputs-patientone/patient-data/PAT001-OVC-2025/multiomics/pdx_phosphoproteomics.csv
4. Summarize findings across clinical, spatial, and multi-omics data""",

    "Batch Correction": "I have 3 batches of proteomics data with batch effects. Apply ComBat batch correction and verify PC1 no longer correlates with batch.",

    "Predict Treatment Response": "Load the GSE184880 ovarian cancer dataset, setup a GEARS model, train it, and predict how Patient-001's T cells will respond to checkpoint inhibitor therapy.",

    "Immunotherapy Prediction": "Predict the response of T cells to anti-PD1 and anti-CTLA4 combination therapy. Which genes are most upregulated?",

    "Drug Screening": "For Patient-001 with ovarian cancer, test responses to: 1) Checkpoint inhibitors (PD1/CTLA4), 2) PARP inhibitors, 3) Platinum therapy. Which shows the best predicted response?",

    "Quantum Cell Type Fidelity": "Train quantum embeddings on the T-cell spatial transcriptomics data. Compute fidelity scores and identify cells with immune evasion states near the tumor boundary.",

    "Immune Evasion Detection": "Using quantum fidelity analysis, identify T-cells that are evading immune surveillance. What is the evasion score for cells near the tumor margin?",

    "TLS Analysis": "Analyze the quantum signatures of tertiary lymphoid structures in the spatial data. Which TLS candidates show the highest quantum coherence?",

    "Quantum + GEARS Validation": "First predict T-cell response to checkpoint inhibitors using GEARS. Then encode the predicted gene expression changes into quantum states and compute fidelity changes. Do the quantum and GEARS predictions agree?",
}
