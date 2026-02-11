"""Unit tests for mcp-genomic-results server."""

import os
import pytest


# ---------------------------------------------------------------------------
# TestServerImport - verify the module loads without errors
# ---------------------------------------------------------------------------


class TestServerImport:
    """Verify the server module can be imported."""

    def test_import_server(self):
        from mcp_genomic_results.server import mcp
        assert mcp.name == "genomic-results"

    def test_import_annotations(self):
        from mcp_genomic_results.annotations import (
            VARIANT_ANNOTATIONS,
            CN_ANNOTATIONS,
            HRD_GENES,
            OVC_GENE_PANEL,
            get_variant_annotation,
            get_cn_annotation,
        )
        assert "TP53_R175H" in VARIANT_ANNOTATIONS
        assert "MYC" in CN_ANNOTATIONS
        assert "BRCA1" in HRD_GENES
        assert "TP53" in OVC_GENE_PANEL

    def test_annotation_lookup(self):
        from mcp_genomic_results.annotations import get_variant_annotation, get_cn_annotation

        tp53 = get_variant_annotation("TP53_R175H")
        assert tp53 is not None
        assert tp53["gene"] == "TP53"
        assert tp53["cosmic_id"] == "COSV57271936"

        myc = get_cn_annotation("MYC")
        assert myc is not None
        assert myc["alteration"] == "Amplification"

        assert get_variant_annotation("NONEXISTENT") is None
        assert get_cn_annotation("NONEXISTENT") is None


# ---------------------------------------------------------------------------
# TestVcfParsing - verify pure-Python VCF parser
# ---------------------------------------------------------------------------


class TestVcfParsing:
    """Test VCF file parsing."""

    def test_parse_vcf_all_variants(self, sample_vcf_file):
        from mcp_genomic_results.server import _parse_vcf_file

        variants = _parse_vcf_file(sample_vcf_file, min_af=0.0)
        assert len(variants) == 7

    def test_parse_vcf_filter_by_af(self, sample_vcf_file):
        from mcp_genomic_results.server import _parse_vcf_file

        variants = _parse_vcf_file(sample_vcf_file, min_af=0.40)
        assert len(variants) == 3  # TP53 (0.73), PIK3CA (0.42), PTEN (0.85)
        genes = {v["gene"] for v in variants}
        assert "TP53" in genes
        assert "PIK3CA" in genes
        assert "PTEN" in genes

    def test_parse_vcf_variant_fields(self, sample_vcf_file):
        from mcp_genomic_results.server import _parse_vcf_file

        variants = _parse_vcf_file(sample_vcf_file, min_af=0.0)
        tp53 = [v for v in variants if v["gene"] == "TP53"][0]

        assert tp53["chrom"] == "chr17"
        assert tp53["pos"] == 7578406
        assert tp53["id"] == "TP53_R175H"
        assert tp53["allele_frequency"] == 0.73
        assert tp53["depth"] == 245
        assert tp53["effect"] == "missense_variant"
        assert tp53["cosmic_id"] == "COSV57271936"

    def test_parse_vcf_annotations_attached(self, sample_vcf_file):
        from mcp_genomic_results.server import _parse_vcf_file

        variants = _parse_vcf_file(sample_vcf_file, min_af=0.0)
        tp53 = [v for v in variants if v["gene"] == "TP53"][0]
        assert "annotation" in tp53
        assert tp53["annotation"]["classification"] == "Pathogenic"

    def test_parse_vcf_missing_file(self):
        from mcp_genomic_results.server import _parse_vcf_file

        with pytest.raises(FileNotFoundError):
            _parse_vcf_file("/nonexistent/path.vcf")


# ---------------------------------------------------------------------------
# TestCnsParsing - verify pure-Python CNS parser
# ---------------------------------------------------------------------------


class TestCnsParsing:
    """Test CNS file parsing."""

    def test_parse_cns_classification(self, sample_cns_file):
        from mcp_genomic_results.server import _parse_cns_file

        result = _parse_cns_file(sample_cns_file)
        assert result["total_segments"] == 6
        assert len(result["amplifications"]) == 2  # MYC, CCNE1
        assert len(result["deletions"]) == 2  # RB1, CDKN2A
        assert len(result["neutral"]) == 2  # TP53, BRAF

    def test_parse_cns_amplification_genes(self, sample_cns_file):
        from mcp_genomic_results.server import _parse_cns_file

        result = _parse_cns_file(sample_cns_file)
        amp_genes = {s["gene"] for s in result["amplifications"]}
        assert amp_genes == {"MYC", "CCNE1"}

    def test_parse_cns_deletion_genes(self, sample_cns_file):
        from mcp_genomic_results.server import _parse_cns_file

        result = _parse_cns_file(sample_cns_file)
        del_genes = {s["gene"] for s in result["deletions"]}
        assert del_genes == {"RB1", "CDKN2A"}

    def test_parse_cns_annotations_attached(self, sample_cns_file):
        from mcp_genomic_results.server import _parse_cns_file

        result = _parse_cns_file(sample_cns_file)
        myc = [s for s in result["amplifications"] if s["gene"] == "MYC"][0]
        assert "annotation" in myc
        assert myc["annotation"]["alteration"] == "Amplification"

    def test_parse_cns_custom_thresholds(self, sample_cns_file):
        from mcp_genomic_results.server import _parse_cns_file

        # Very high threshold - nothing qualifies as amplification
        result = _parse_cns_file(sample_cns_file, amp_threshold=2.0, del_threshold=-2.5)
        assert len(result["amplifications"]) == 0
        assert len(result["deletions"]) == 0

    def test_parse_cns_missing_file(self):
        from mcp_genomic_results.server import _parse_cns_file

        with pytest.raises(FileNotFoundError):
            _parse_cns_file("/nonexistent/path.cns")


# ---------------------------------------------------------------------------
# TestToolsDryRun - verify DRY_RUN returns synthetic data
# ---------------------------------------------------------------------------


class TestToolsDryRun:
    """Test tools in DRY_RUN mode return synthetic data."""

    @pytest.fixture(autouse=True)
    def enable_dry_run(self, monkeypatch):
        monkeypatch.setenv("GENOMIC_RESULTS_DRY_RUN", "true")
        # Reload the module to pick up env change
        import importlib
        import mcp_genomic_results.server as srv
        monkeypatch.setattr(srv, "DRY_RUN", True)

    @pytest.mark.asyncio
    async def test_parse_somatic_variants_dry_run(self):
        from mcp_genomic_results.server import parse_somatic_variants

        result = await parse_somatic_variants.fn(vcf_path="/any/path.vcf")
        assert "_DRY_RUN_WARNING" in result
        assert result["total_variants"] == 12
        assert len(result["somatic_mutations"]) == 3

    @pytest.mark.asyncio
    async def test_parse_cnv_calls_dry_run(self):
        from mcp_genomic_results.server import parse_cnv_calls

        result = await parse_cnv_calls.fn(cns_path="/any/path.cns")
        assert "_DRY_RUN_WARNING" in result
        assert len(result["amplifications"]) == 3

    @pytest.mark.asyncio
    async def test_calculate_hrd_dry_run(self):
        from mcp_genomic_results.server import calculate_hr_deficiency_score

        result = await calculate_hr_deficiency_score.fn(
            vcf_path="/any/path.vcf", cns_path="/any/path.cns"
        )
        assert "_DRY_RUN_WARNING" in result
        assert result["hrd_score"] == 44
        assert result["parp_eligible"] is True

    @pytest.mark.asyncio
    async def test_generate_report_dry_run(self):
        from mcp_genomic_results.server import generate_genomic_report

        result = await generate_genomic_report.fn(
            vcf_path="/any/path.vcf", cns_path="/any/path.cns"
        )
        assert "_DRY_RUN_WARNING" in result
        assert result["report_type"] == "Comprehensive Genomic Report"
        assert len(result["therapy_recommendations"]) > 0


# ---------------------------------------------------------------------------
# TestToolsRealData - verify tools with real fixture data
# ---------------------------------------------------------------------------


class TestToolsRealData:
    """Test tools with real fixture files (DRY_RUN=false)."""

    @pytest.fixture(autouse=True)
    def disable_dry_run(self, monkeypatch):
        monkeypatch.setenv("GENOMIC_RESULTS_DRY_RUN", "false")
        import mcp_genomic_results.server as srv
        monkeypatch.setattr(srv, "DRY_RUN", False)

    @pytest.mark.asyncio
    async def test_parse_somatic_variants_real(self, sample_vcf_file):
        from mcp_genomic_results.server import parse_somatic_variants

        # Use min_allele_frequency=0.0 to include all variants (incl. wild-type/CN)
        result = await parse_somatic_variants.fn(
            vcf_path=sample_vcf_file, min_allele_frequency=0.0
        )
        assert "_DRY_RUN_WARNING" not in result
        assert result["total_variants"] == 7
        assert len(result["somatic_mutations"]) == 3
        mut_genes = {m["gene"] for m in result["somatic_mutations"]}
        assert mut_genes == {"TP53", "PIK3CA", "PTEN"}

    @pytest.mark.asyncio
    async def test_parse_cnv_calls_real(self, sample_cns_file):
        from mcp_genomic_results.server import parse_cnv_calls

        result = await parse_cnv_calls.fn(cns_path=sample_cns_file)
        assert "_DRY_RUN_WARNING" not in result
        assert len(result["amplifications"]) == 2
        assert len(result["deletions"]) == 2
        assert result["total_segments"] == 6
