"""Fixtures for mcp-genomic-results tests."""

import os
import tempfile
from pathlib import Path

import pytest

FIXTURES_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture
def sample_vcf_file() -> str:
    """Path to the test somatic_variants.vcf fixture."""
    return str(FIXTURES_DIR / "test_somatic_variants.vcf")


@pytest.fixture
def sample_cns_file() -> str:
    """Path to the test copy_number_results.cns fixture."""
    return str(FIXTURES_DIR / "test_copy_number_results.cns")


@pytest.fixture
def temp_dir():
    """Provide a temporary directory for test outputs."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


@pytest.fixture(autouse=True)
def set_dry_run_env(monkeypatch):
    """Ensure DRY_RUN is disabled for unit tests so real parsing is exercised."""
    monkeypatch.setenv("GENOMIC_RESULTS_DRY_RUN", "false")
