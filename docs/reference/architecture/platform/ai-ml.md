# AI/ML Model Inference Architecture

**Status:** Served by external Hugging Face MCP server
**Last Updated:** 2026-02-21

---

## Overview

AI/ML model inference for genomic analysis (cell type prediction, sequence embedding, foundation model inference) is now provided by the **external Hugging Face MCP server** rather than a local mocked server.

The external Hugging Face MCP server provides 7 real tools with direct Hub access, replacing the previous 3-tool mocked `mcp-huggingface` server.

**Setup instructions:** [CONNECT_EXTERNAL_MCP.md](../../../../docs/for-researchers/CONNECT_EXTERNAL_MCP.md)

---

## Supported Foundation Models

The following genomic foundation models are accessible via the external Hugging Face MCP server:

| Model | Architecture | Parameters | Use Case |
|-------|-------------|------------|----------|
| **DNABERT-2** | BERT-based transformer | 117M | DNA sequence embeddings, promoter/enhancer prediction |
| **Geneformer** | GPT-based transformer | 10M+ cells pre-trained | Cell type annotation from scRNA-seq |
| **Nucleotide-Transformer** | Transformer | 2.5B | Long-range regulatory interactions, splice site prediction |
| **scGPT** | Generative transformer | 33M | Single-cell data analysis |

---

## Integration with Other Modalities

### With Spatial Transcriptomics (mcp-spatialtools)

**Integration Point:** Cell type deconvolution validation
- Run cell type deconvolution on spatial data (signature-based)
- Validate predicted cell types using foundation model predictions via the external HF server
- Map cell types to spatial coordinates

### With Multi-omics (mcp-multiomics)

**Integration Point:** Feature extraction for multi-modal integration
- Generate sequence embeddings for genomic regions
- Combine with RNA/protein expression features from mcp-multiomics

### With Genomic Analysis (mcp-fgbio)

**Integration Point:** Variant annotation enhancement
- Extract variant context sequences
- Generate embeddings for variant effect prediction

---

## External Resources

- **Hugging Face Models:** [https://huggingface.co/models](https://huggingface.co/models)
- **DNABERT-2:** [Paper](https://arxiv.org/abs/2306.15006)
- **Geneformer:** [Paper](https://www.nature.com/articles/s41586-023-06139-9)

---

## Related Workflows

- [Spatial Transcriptomics](../spatial/README.md) - Cell type deconvolution validation
- [Multi-omics Integration](../rna/multiomics.md) - Multi-modal feature fusion
- [Genomic Analysis](../dna/genomic-cohorts.md) - Variant effect prediction
- [External MCP Servers](../../../../docs/for-researchers/CONNECT_EXTERNAL_MCP.md) - Setup guide

---

**See also:** [Main Architecture](../README.md) | [Hugging Face Hub](https://huggingface.co/)
