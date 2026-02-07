"""Tests for PatientReportData Pydantic models."""

import json
from pathlib import Path
import pytest

from mcp_patient_report.models import (
    PatientReportData,
    PatientInfo,
    DiagnosisSummary,
    GenomicFinding,
    TreatmentOption,
    MonitoringPlan,
    EvidenceLevel,
)


FIXTURES_DIR = Path(__file__).parent / "fixtures"


class TestPatientInfo:
    """Tests for PatientInfo model."""

    def test_valid_patient_info(self):
        """Test creating valid patient info."""
        info = PatientInfo(
            name="Sarah Anderson",
            age=58,
            sex="Female",
            patient_id="PAT001-OVC-2025",
            diagnosis="Stage IV Ovarian Cancer"
        )
        assert info.name == "Sarah Anderson"
        assert info.age == 58

    def test_age_validation(self):
        """Test age must be within valid range."""
        with pytest.raises(ValueError):
            PatientInfo(
                name="Test",
                age=-1,
                sex="Female",
                patient_id="TEST",
                diagnosis="Test"
            )

        with pytest.raises(ValueError):
            PatientInfo(
                name="Test",
                age=200,
                sex="Female",
                patient_id="TEST",
                diagnosis="Test"
            )


class TestGenomicFinding:
    """Tests for GenomicFinding model."""

    def test_valid_finding(self):
        """Test creating valid genomic finding."""
        finding = GenomicFinding(
            gene="BRCA1",
            variant="Pathogenic germline mutation",
            significance="Pathogenic",
            plain_language="You carry a BRCA1 mutation.",
            actionability="PARP inhibitors available",
            icon="positive"
        )
        assert finding.gene == "BRCA1"
        assert finding.significance == "Pathogenic"

    def test_vaf_validation(self):
        """Test VAF must be between 0 and 1."""
        finding = GenomicFinding(
            gene="TP53",
            variant="R175H",
            significance="Pathogenic",
            variant_allele_frequency=0.85,
            plain_language="Test"
        )
        assert finding.variant_allele_frequency == 0.85

        with pytest.raises(ValueError):
            GenomicFinding(
                gene="TP53",
                variant="R175H",
                significance="Pathogenic",
                variant_allele_frequency=1.5,
                plain_language="Test"
            )


class TestTreatmentOption:
    """Tests for TreatmentOption model."""

    def test_valid_treatment(self):
        """Test creating valid treatment option."""
        treatment = TreatmentOption(
            name="Olaparib",
            type="Targeted Therapy",
            evidence_level=EvidenceLevel.FDA_APPROVED,
            plain_language_description="A pill that blocks PARP enzyme.",
            why_recommended="Your BRCA1 mutation makes you a good candidate."
        )
        assert treatment.name == "Olaparib"
        assert treatment.evidence_level == EvidenceLevel.FDA_APPROVED

    def test_evidence_levels(self):
        """Test all evidence levels are valid."""
        for level in EvidenceLevel:
            treatment = TreatmentOption(
                name="Test",
                type="Test",
                evidence_level=level,
                plain_language_description="Test",
                why_recommended="Test"
            )
            assert treatment.evidence_level == level


class TestPatientReportData:
    """Tests for complete PatientReportData model."""

    def test_load_pat001_fixture(self):
        """Test loading and validating PAT001 fixture data."""
        fixture_path = FIXTURES_DIR / "pat001_report_data.json"

        with open(fixture_path) as f:
            data = json.load(f)

        report = PatientReportData(**data)

        assert report.patient_info.patient_id == "PAT001-OVC-2025"
        assert report.patient_info.name == "Sarah Anderson"
        assert report.diagnosis_summary.cancer_type == "Ovarian Cancer"
        assert len(report.genomic_findings) == 4
        assert len(report.treatment_options) == 3
        assert report.spatial_findings is not None
        assert report.histology_findings is not None

    def test_minimal_valid_report(self):
        """Test creating report with only required fields."""
        report = PatientReportData(
            patient_info=PatientInfo(
                name="Test Patient",
                age=50,
                sex="Female",
                patient_id="TEST-001",
                diagnosis="Test Diagnosis"
            ),
            diagnosis_summary=DiagnosisSummary(
                cancer_type="Test Cancer",
                stage="Stage I",
                plain_language_description="This is a test."
            ),
            monitoring_plan=MonitoringPlan(
                warning_signs=["Test warning"]
            )
        )

        assert report.patient_info.patient_id == "TEST-001"
        assert report.spatial_findings is None
        assert report.histology_findings is None

    def test_json_schema_generation(self):
        """Test that JSON schema can be generated."""
        schema = PatientReportData.model_json_schema()

        assert "properties" in schema
        assert "patient_info" in schema["properties"]
        assert "diagnosis_summary" in schema["properties"]

    def test_json_serialization(self):
        """Test report can be serialized to JSON."""
        fixture_path = FIXTURES_DIR / "pat001_report_data.json"

        with open(fixture_path) as f:
            data = json.load(f)

        report = PatientReportData(**data)

        # Serialize and deserialize
        json_str = report.model_dump_json()
        parsed = json.loads(json_str)

        assert parsed["patient_info"]["patient_id"] == "PAT001-OVC-2025"
