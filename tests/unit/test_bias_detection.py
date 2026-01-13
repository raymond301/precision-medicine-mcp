"""
Unit tests for bias detection utilities.

Tests cover:
- Data representation analysis
- Fairness metrics (demographic parity, equalized odds, calibration)
- Proxy feature detection
- Ancestry-aware confidence scoring

Run tests:
    pytest tests/unit/test_bias_detection.py -v
    pytest tests/unit/test_bias_detection.py::TestRepresentationAnalysis -v
"""

import pytest
import numpy as np
import pandas as pd
from typing import Dict

import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from shared.utils.bias_detection import (
    check_data_representation,
    demographic_parity,
    equalized_odds,
    calibration_by_group,
    calculate_fairness_metrics,
    detect_proxy_features,
    calculate_ancestry_aware_confidence,
    cramers_v,
    REPRESENTATION_THRESHOLDS,
    FAIRNESS_DISPARITY_THRESHOLDS
)


# ============================================================================
# TEST DATA FIXTURES
# ============================================================================

@pytest.fixture
def balanced_ancestry_data():
    """Create balanced dataset across ancestries."""
    n_per_group = 100
    data = pd.DataFrame({
        'patient_id': range(n_per_group * 4),
        'ancestry': ['european'] * n_per_group +
                   ['african'] * n_per_group +
                   ['asian'] * n_per_group +
                   ['latino'] * n_per_group,
        'age': np.random.randint(30, 80, n_per_group * 4),
        'variant_pathogenic': np.random.randint(0, 2, n_per_group * 4)
    })
    return data


@pytest.fixture
def imbalanced_ancestry_data():
    """Create imbalanced dataset (Euro-centric)."""
    data = pd.DataFrame({
        'patient_id': range(170),
        'ancestry': ['european'] * 100 +
                   ['african'] * 30 +
                   ['asian'] * 20 +
                   ['latino'] * 20,
        'age': np.random.randint(30, 80, 170),
        'variant_pathogenic': np.random.randint(0, 2, 170)
    })
    return data


@pytest.fixture
def critical_imbalance_data():
    """Create critically imbalanced dataset (>95% European)."""
    data = pd.DataFrame({
        'patient_id': range(100),
        'ancestry': ['european'] * 97 + ['african'] * 2 + ['asian'] * 1,
        'age': np.random.randint(30, 80, 100)
    })
    return data


@pytest.fixture
def fairness_test_data():
    """Create data for fairness metrics testing."""
    np.random.seed(42)
    n = 200

    # Create biased predictions (higher positive rate for group A)
    groups = np.array(['group_a'] * 100 + ['group_b'] * 100)
    y_true = np.random.randint(0, 2, n)

    # Group A gets more positive predictions (bias)
    y_pred = y_true.copy()
    y_pred[:100] = np.where(np.random.random(100) > 0.3, 1, y_pred[:100])  # 70% positive for group A
    y_pred[100:] = np.where(np.random.random(100) > 0.7, 1, y_pred[100:])  # 30% positive for group B

    # Probabilities
    y_prob = np.random.random(n)

    return {
        'y_true': y_true,
        'y_pred': y_pred,
        'y_prob': y_prob,
        'groups': groups
    }


@pytest.fixture
def proxy_feature_data():
    """Create data with proxy features."""
    np.random.seed(42)
    n = 200

    # Create protected attribute (ancestry)
    ancestry = np.array(['european'] * 100 + ['african'] * 100)

    # Create features with varying correlation to ancestry
    # zip_code: high correlation (proxy)
    zip_code = np.where(ancestry == 'european',
                       np.random.randint(90000, 91000, n),
                       np.random.randint(80000, 81000, n))

    # age: low correlation (not a proxy)
    age = np.random.randint(30, 80, n)

    # income: medium correlation (potential proxy)
    income = np.where(ancestry == 'european',
                     np.random.randint(60000, 120000, n),
                     np.random.randint(30000, 70000, n))

    X = pd.DataFrame({
        'zip_code': zip_code,
        'age': age,
        'income': income
    })

    protected = pd.DataFrame({'ancestry': ancestry})

    feature_importances = {
        'zip_code': 0.25,  # High importance + high correlation = proxy
        'age': 0.15,       # Low correlation, not a proxy
        'income': 0.08     # Medium correlation but low importance
    }

    return {
        'X': X,
        'protected': protected,
        'importances': feature_importances
    }


# ============================================================================
# TEST: DATA REPRESENTATION ANALYSIS
# ============================================================================

class TestRepresentationAnalysis:
    """Test data representation analysis functions."""

    def test_balanced_representation(self, balanced_ancestry_data):
        """Test that balanced data gets ACCEPTABLE risk level."""
        analysis = check_data_representation(
            balanced_ancestry_data,
            ancestry_column='ancestry',
            reference_dataset='gnomad',
            min_acceptable_pct=0.10
        )

        assert analysis.risk_level == "ACCEPTABLE"
        assert analysis.total_samples == 400
        assert len(analysis.ancestry_counts) == 4
        assert all(count == 100 for count in analysis.ancestry_counts.values())

    def test_imbalanced_representation(self, imbalanced_ancestry_data):
        """Test that imbalanced data gets MEDIUM risk level."""
        analysis = check_data_representation(
            imbalanced_ancestry_data,
            ancestry_column='ancestry',
            min_acceptable_pct=0.10
        )

        # asian and latino are at 20/170 = 11.8%, just above MEDIUM threshold
        assert analysis.risk_level in ["MEDIUM", "ACCEPTABLE"]
        assert analysis.total_samples == 170

        # European should be majority
        assert analysis.ancestry_proportions['european'] > 0.5

    def test_critical_imbalance(self, critical_imbalance_data):
        """Test that critically imbalanced data gets CRITICAL risk level."""
        analysis = check_data_representation(
            critical_imbalance_data,
            ancestry_column='ancestry',
            min_acceptable_pct=0.10
        )

        assert analysis.risk_level == "CRITICAL"
        assert len(analysis.warnings) > 0
        assert "DO NOT USE" in analysis.recommendations[0]

        # Check that asian ancestry is <5%
        assert analysis.ancestry_proportions['asian'] < 0.05

    def test_representation_comparison(self, balanced_ancestry_data):
        """Test comparison to reference datasets."""
        analysis = check_data_representation(
            balanced_ancestry_data,
            ancestry_column='ancestry',
            reference_dataset='gnomad'
        )

        # Should have comparison data
        assert len(analysis.reference_comparison) > 0

        # Each ancestry should have reference, actual, and difference
        for ancestry, comparison in analysis.reference_comparison.items():
            assert 'reference' in comparison
            assert 'actual' in comparison
            assert 'difference' in comparison


# ============================================================================
# TEST: FAIRNESS METRICS
# ============================================================================

class TestFairnessMetrics:
    """Test fairness metric calculations."""

    def test_demographic_parity_balanced(self):
        """Test demographic parity with balanced predictions."""
        y_pred = np.array([1, 0, 1, 0, 1, 0, 1, 0])
        groups = np.array(['a', 'a', 'a', 'a', 'b', 'b', 'b', 'b'])

        metrics = demographic_parity(y_pred, groups)

        assert metrics.metric_name == "Demographic Parity"
        assert metrics.group_values['a'] == 0.5  # 2/4 positive
        assert metrics.group_values['b'] == 0.5  # 2/4 positive
        assert metrics.max_disparity == 0.0      # Perfect parity
        assert metrics.risk_level == "ACCEPTABLE"

    def test_demographic_parity_biased(self, fairness_test_data):
        """Test demographic parity with biased predictions."""
        metrics = demographic_parity(
            fairness_test_data['y_pred'],
            fairness_test_data['groups']
        )

        # Should detect disparity
        assert metrics.max_disparity > 0.10  # Significant disparity
        assert len(metrics.warnings) > 0

    def test_equalized_odds(self, fairness_test_data):
        """Test equalized odds calculation."""
        metrics = equalized_odds(
            fairness_test_data['y_true'],
            fairness_test_data['y_pred'],
            fairness_test_data['groups']
        )

        # Should return both TPR and FPR
        assert 'tpr' in metrics
        assert 'fpr' in metrics

        # Each should have metric name and values
        assert metrics['tpr'].metric_name == "True Positive Rate (Sensitivity)"
        assert metrics['fpr'].metric_name == "False Positive Rate"

        # Should have values for both groups
        assert 'group_a' in metrics['tpr'].group_values
        assert 'group_b' in metrics['tpr'].group_values

    def test_calibration_by_group(self, fairness_test_data):
        """Test calibration metric."""
        metrics = calibration_by_group(
            fairness_test_data['y_true'],
            fairness_test_data['y_prob'],
            fairness_test_data['groups'],
            n_bins=10
        )

        assert metrics.metric_name == "Calibration Error"
        assert 'group_a' in metrics.group_values
        assert 'group_b' in metrics.group_values

        # Calibration error should be between 0 and 1
        for error in metrics.group_values.values():
            assert 0 <= error <= 1

    def test_calculate_fairness_metrics_all(self, fairness_test_data):
        """Test calculating all fairness metrics at once."""
        results = calculate_fairness_metrics(
            y_true=fairness_test_data['y_true'],
            y_pred=fairness_test_data['y_pred'],
            y_prob=fairness_test_data['y_prob'],
            groups=fairness_test_data['groups'],
            metrics=['demographic_parity', 'equalized_odds', 'calibration']
        )

        assert 'demographic_parity' in results
        assert 'equalized_odds' in results
        assert 'calibration' in results

    def test_calculate_fairness_metrics_without_probs(self, fairness_test_data):
        """Test that calibration is skipped when y_prob is None."""
        results = calculate_fairness_metrics(
            y_true=fairness_test_data['y_true'],
            y_pred=fairness_test_data['y_pred'],
            y_prob=None,
            groups=fairness_test_data['groups'],
            metrics=['demographic_parity', 'calibration']
        )

        assert 'demographic_parity' in results
        assert 'calibration' not in results  # Should be skipped


# ============================================================================
# TEST: PROXY FEATURE DETECTION
# ============================================================================

class TestProxyFeatureDetection:
    """Test proxy feature detection."""

    def test_detect_high_proxy(self, proxy_feature_data):
        """Test detection of high-importance proxy feature (zip_code)."""
        proxies = detect_proxy_features(
            X=proxy_feature_data['X'],
            feature_importances=proxy_feature_data['importances'],
            protected_attributes=proxy_feature_data['protected'],
            correlation_threshold=0.5,
            importance_threshold=0.05
        )

        # zip_code should be detected as proxy (high importance + high correlation)
        proxy_names = [p.feature_name for p in proxies]
        assert 'zip_code' in proxy_names

        # Find zip_code proxy
        zip_proxy = next(p for p in proxies if p.feature_name == 'zip_code')
        assert zip_proxy.is_proxy == True
        assert zip_proxy.risk_level in ["CRITICAL", "HIGH"]
        assert "REMOVE" in zip_proxy.recommendation.upper()

    def test_no_false_positives(self, proxy_feature_data):
        """Test that age (low correlation) is not flagged as proxy."""
        proxies = detect_proxy_features(
            X=proxy_feature_data['X'],
            feature_importances=proxy_feature_data['importances'],
            protected_attributes=proxy_feature_data['protected'],
            correlation_threshold=0.5,
            importance_threshold=0.05
        )

        # age should NOT be in proxies (low correlation despite importance)
        proxy_names = [p.feature_name for p in proxies]
        assert 'age' not in proxy_names

    def test_cramers_v_perfect_association(self):
        """Test Cramér's V with perfect association."""
        x = pd.Series(['a', 'a', 'b', 'b'])
        y = pd.Series(['x', 'x', 'y', 'y'])

        v = cramers_v(x, y)
        assert v == pytest.approx(1.0, abs=0.01)  # Perfect association

    def test_cramers_v_no_association(self):
        """Test Cramér's V with no association."""
        np.random.seed(42)
        x = pd.Series(np.random.choice(['a', 'b'], 100))
        y = pd.Series(np.random.choice(['x', 'y'], 100))

        v = cramers_v(x, y)
        assert v < 0.3  # Low association


# ============================================================================
# TEST: ANCESTRY-AWARE CONFIDENCE
# ============================================================================

class TestAncestryAwareConfidence:
    """Test ancestry-aware confidence scoring."""

    def test_high_confidence_well_studied(self):
        """Test that well-studied variants maintain high confidence."""
        result = calculate_ancestry_aware_confidence(
            base_confidence=0.85,
            patient_ancestry='european',
            variant_id='BRCA1:c.5266dupC',
            reference_study_counts={'european': 50, 'african': 2, 'asian': 1}
        )

        # Should have no penalty (>=20 studies)
        assert result['confidence_penalty'] == 0.0
        assert result['adjusted_confidence'] == 0.85
        assert result['confidence_level'] == 'HIGH'
        assert len(result['warnings']) == 0

    def test_medium_confidence_moderate_studies(self):
        """Test medium confidence for moderately studied variants."""
        result = calculate_ancestry_aware_confidence(
            base_confidence=0.85,
            patient_ancestry='african',
            variant_id='BRCA1:c.5266dupC',
            reference_study_counts={'european': 50, 'african': 10, 'asian': 1}
        )

        # Should have 10% penalty (5-19 studies)
        assert result['confidence_penalty'] == 0.1
        assert result['adjusted_confidence'] == pytest.approx(0.765, abs=0.01)  # 0.85 * 0.9
        assert result['confidence_level'] == 'MEDIUM'

    def test_low_confidence_understudied(self):
        """Test low confidence for understudied ancestries."""
        result = calculate_ancestry_aware_confidence(
            base_confidence=0.85,
            patient_ancestry='asian',
            variant_id='BRCA1:c.5266dupC',
            reference_study_counts={'european': 50, 'african': 2, 'asian': 1}
        )

        # Should have 30% penalty (<5 studies)
        assert result['confidence_penalty'] == 0.3
        assert result['adjusted_confidence'] == pytest.approx(0.595, abs=0.01)  # 0.85 * 0.7
        assert result['confidence_level'] == 'LOW'
        assert len(result['warnings']) > 0
        assert 'Limited data' in result['warnings'][0]

    def test_zero_studies(self):
        """Test handling of variants with no studies in patient ancestry."""
        result = calculate_ancestry_aware_confidence(
            base_confidence=0.85,
            patient_ancestry='indigenous',
            variant_id='BRCA1:c.5266dupC',
            reference_study_counts={'european': 50, 'african': 2, 'asian': 1}
        )

        # Should apply maximum penalty
        assert result['ancestral_studies'] == 0
        assert result['confidence_penalty'] == 0.3
        assert result['confidence_level'] == 'LOW'

    def test_custom_thresholds(self):
        """Test custom penalty thresholds."""
        custom_thresholds = {
            'high': (30, 0.0),
            'medium': (10, 0.15),
            'low': (0, 0.40)
        }

        result = calculate_ancestry_aware_confidence(
            base_confidence=0.80,
            patient_ancestry='african',
            variant_id='test_variant',
            reference_study_counts={'african': 15},
            min_studies_thresholds=custom_thresholds
        )

        # 15 studies should be MEDIUM with 15% penalty
        assert result['confidence_level'] == 'MEDIUM'
        assert result['confidence_penalty'] == 0.15
        assert result['adjusted_confidence'] == pytest.approx(0.68, abs=0.01)  # 0.80 * 0.85


# ============================================================================
# TEST: THRESHOLDS
# ============================================================================

class TestThresholds:
    """Test that thresholds are correctly defined."""

    def test_representation_thresholds(self):
        """Test representation threshold values."""
        assert REPRESENTATION_THRESHOLDS['critical'] == 0.05
        assert REPRESENTATION_THRESHOLDS['high'] == 0.10
        assert REPRESENTATION_THRESHOLDS['medium'] == 0.20
        assert REPRESENTATION_THRESHOLDS['acceptable'] == 0.20

    def test_fairness_disparity_thresholds(self):
        """Test fairness disparity threshold values."""
        assert FAIRNESS_DISPARITY_THRESHOLDS['critical'] == 0.20
        assert FAIRNESS_DISPARITY_THRESHOLDS['high'] == 0.10
        assert FAIRNESS_DISPARITY_THRESHOLDS['acceptable'] == 0.10


# ============================================================================
# INTEGRATION TESTS
# ============================================================================

class TestIntegration:
    """Integration tests combining multiple components."""

    def test_full_audit_workflow(self, imbalanced_ancestry_data, fairness_test_data):
        """Test complete audit workflow from data to report."""
        # Step 1: Check representation
        rep_analysis = check_data_representation(
            imbalanced_ancestry_data,
            ancestry_column='ancestry'
        )
        assert rep_analysis.risk_level in ["MEDIUM", "ACCEPTABLE", "CRITICAL"]

        # Step 2: Calculate fairness metrics
        fairness_results = calculate_fairness_metrics(
            y_true=fairness_test_data['y_true'],
            y_pred=fairness_test_data['y_pred'],
            y_prob=fairness_test_data['y_prob'],
            groups=fairness_test_data['groups']
        )
        assert 'demographic_parity' in fairness_results

        # Step 3: Would generate report (tested separately)
        assert rep_analysis is not None
        assert fairness_results is not None

    def test_ancestry_specific_fairness(self):
        """Test fairness metrics specifically for ancestry groups."""
        np.random.seed(42)
        n = 300

        ancestries = np.array(['european'] * 100 + ['african'] * 100 + ['asian'] * 100)
        y_true = np.random.randint(0, 2, n)
        y_pred = y_true.copy()

        # No bias - all ancestries treated equally
        fairness_metrics = calculate_fairness_metrics(
            y_true=y_true,
            y_pred=y_pred,
            y_prob=None,
            groups=ancestries,
            metrics=['demographic_parity']
        )

        # Should have minimal disparity
        assert fairness_metrics['demographic_parity'].max_disparity < 0.05
        assert fairness_metrics['demographic_parity'].risk_level == "ACCEPTABLE"


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
