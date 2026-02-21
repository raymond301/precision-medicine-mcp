# Appendix B: Infrastructure Components

[← Back to SLA Overview](README.md)

---

## Platform Infrastructure Inventory

### Compute (Cloud Run)
- **MCP Servers** ([Server Registry](../reference/shared/server-registry.md)): fgbio, multiomics, spatialtools, tcga, openimagedata, deepcell, seqera, mockepic, epic, cell-classify, genomic-results, patient-report, perturbation, quantum-celltype-fidelity.
- **Auto-scaling:** Min 1, Max 100 instances.
- **Resource Limits:** 2 vCPU, 4GB Memory per instance.

### Storage (GCP)
- **Cloud Storage:** patient-data-uploads, analysis-outputs.
- **Healthcare API:** FHIR Store with dataset precision-medicine.
- **Encryption:** AES-256 for all data-at-rest.

### Networking
- **Load Balancer:** Global HTTPS LB with Cloud Armor WAF.
- **Firewall:** Egress allowlist for 13 specific domains.

---

## Related Documents

- [Section 02: Service Level Objectives](02_SERVICE_LEVEL_OBJECTIVES.md)
- [Section 05: Data Protection](05_DATA_PROTECTION.md)

---

### Document History

| Date | Version | Author | Change Summary |
|:---|:---|:---|:---|
| 2026-02-08 | 2.0 | IT Operations | Initial creation |
