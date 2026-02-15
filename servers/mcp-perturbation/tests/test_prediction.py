"""Tests for prediction functionality."""

import pytest
from mcp_perturbation.prediction import (
    DifferentialExpressionAnalyzer,
    get_top_changed_genes
)
from mcp_perturbation.data_loader import DatasetLoader


@pytest.fixture
def test_adata():
    """Create test dataset."""
    loader = DatasetLoader()
    return loader.load_geo_dataset("GSE184880", normalize=True, n_hvg=1000)


class TestDifferentialExpressionAnalyzer:
    """Test differential expression analysis."""

    def test_compute_de(self, test_adata):
        """Test computing differential expression."""
        # Split into baseline and predicted (mock)
        n_half = test_adata.n_obs // 2
        baseline = test_adata[:n_half]
        predicted = test_adata[n_half:]

        analyzer = DifferentialExpressionAnalyzer()
        results = analyzer.compute_differential_expression(
            baseline,
            predicted,
            n_top_genes=20,
            method="wilcoxon"
        )

        assert "upregulated_genes" in results
        assert "downregulated_genes" in results
        assert len(results["upregulated_genes"]) <= 20
        assert results["method"] == "wilcoxon"


def test_get_top_changed_genes(test_adata):
    """Test getting top changed genes."""
    # Split data
    n_half = test_adata.n_obs // 2
    baseline = test_adata[:n_half]
    predicted = test_adata[n_half:]

    upregulated, downregulated = get_top_changed_genes(
        baseline,
        predicted,
        n_genes=10
    )

    assert len(upregulated) == 10
    assert len(downregulated) == 10
    assert all(isinstance(g, str) for g in upregulated)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
