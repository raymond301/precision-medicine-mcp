"""Bayesian uncertainty quantification for quantum fidelity predictions.

Provides confidence intervals and uncertainty estimates for fidelity scores
using Bayesian parameter distributions and Monte Carlo sampling.
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from scipy import stats


@dataclass
class UncertaintyEstimate:
    """Container for uncertainty quantification results.

    Attributes:
        mean: Mean prediction value
        std: Standard deviation of predictions
        confidence_interval_95: 95% credible interval (lower, upper)
        confidence_interval_90: 90% credible interval (lower, upper)
        samples: Raw samples from posterior (optional)
        epistemic_uncertainty: Model uncertainty (from parameters)
        aleatoric_uncertainty: Data uncertainty (inherent noise)
    """
    mean: float
    std: float
    confidence_interval_95: Tuple[float, float]
    confidence_interval_90: Tuple[float, float]
    samples: Optional[np.ndarray] = None
    epistemic_uncertainty: float = 0.0
    aleatoric_uncertainty: float = 0.0

    def to_dict(self) -> Dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "mean": float(self.mean),
            "std": float(self.std),
            "confidence_interval_95": {
                "lower": float(self.confidence_interval_95[0]),
                "upper": float(self.confidence_interval_95[1])
            },
            "confidence_interval_90": {
                "lower": float(self.confidence_interval_90[0]),
                "upper": float(self.confidence_interval_90[1])
            },
            "epistemic_uncertainty": float(self.epistemic_uncertainty),
            "aleatoric_uncertainty": float(self.aleatoric_uncertainty),
            "total_uncertainty": float(self.std)
        }


class BayesianParameterDistribution:
    """Bayesian distribution over quantum circuit parameters.

    Maintains mean and covariance of parameter distribution for
    variational quantum circuits. Enables sampling for uncertainty
    quantification.

    Attributes:
        param_mean: Mean parameter values (point estimate)
        param_cov: Covariance matrix (full or diagonal)
        use_diagonal: Whether to use diagonal covariance approximation
    """

    def __init__(
        self,
        param_mean: np.ndarray,
        param_cov: Optional[np.ndarray] = None,
        param_std: Optional[np.ndarray] = None,
        use_diagonal: bool = True
    ):
        """Initialize Bayesian parameter distribution.

        Args:
            param_mean: Mean parameter values (n_params,)
            param_cov: Covariance matrix (n_params, n_params) - optional
            param_std: Standard deviations (n_params,) - alternative to cov
            use_diagonal: Use diagonal covariance (faster, less expressive)
        """
        self.param_mean = np.array(param_mean)
        self.n_params = len(self.param_mean)
        self.use_diagonal = use_diagonal

        # Initialize covariance
        if param_cov is not None:
            self.param_cov = np.array(param_cov)
        elif param_std is not None:
            # Diagonal covariance from standard deviations
            self.param_cov = np.diag(np.array(param_std) ** 2)
        else:
            # Default: Small diagonal covariance (10% of mean)
            default_std = np.abs(self.param_mean) * 0.1 + 0.01
            self.param_cov = np.diag(default_std ** 2)

        # Validate shape
        if self.param_cov.shape != (self.n_params, self.n_params):
            raise ValueError(
                f"Covariance shape {self.param_cov.shape} doesn't match "
                f"parameter shape ({self.n_params},)"
            )

    def sample(self, n_samples: int = 100, random_state: Optional[int] = None) -> np.ndarray:
        """Sample parameter values from posterior distribution.

        Args:
            n_samples: Number of samples to draw
            random_state: Random seed for reproducibility

        Returns:
            Parameter samples (n_samples, n_params)
        """
        rng = np.random.RandomState(random_state)

        if self.use_diagonal:
            # Fast diagonal sampling
            std = np.sqrt(np.diag(self.param_cov))
            samples = rng.normal(
                loc=self.param_mean,
                scale=std,
                size=(n_samples, self.n_params)
            )
        else:
            # Full covariance sampling
            samples = rng.multivariate_normal(
                mean=self.param_mean,
                cov=self.param_cov,
                size=n_samples
            )

        return samples

    def update_from_gradient_history(
        self,
        gradient_history: List[np.ndarray],
        learning_rate: float
    ):
        """Update covariance from gradient history (empirical Bayes).

        Uses gradient variability to estimate parameter uncertainty.

        Args:
            gradient_history: List of gradient vectors from training
            learning_rate: Learning rate used in training
        """
        if len(gradient_history) < 2:
            return

        gradients = np.array(gradient_history)  # (n_steps, n_params)

        # Estimate covariance from gradient variance
        # Higher gradient variance → higher parameter uncertainty
        gradient_cov = np.cov(gradients.T)

        # Scale by learning rate (approximate Hessian^-1)
        self.param_cov = learning_rate ** 2 * gradient_cov

        # Add small diagonal for numerical stability
        self.param_cov += np.eye(self.n_params) * 1e-6


class BayesianFidelityEstimator:
    """Estimate fidelity with uncertainty quantification.

    Uses Monte Carlo sampling from parameter posterior to compute
    confidence intervals on fidelity predictions.
    """

    def __init__(
        self,
        param_distribution: BayesianParameterDistribution,
        n_samples: int = 100
    ):
        """Initialize Bayesian fidelity estimator.

        Args:
            param_distribution: Posterior distribution over parameters
            n_samples: Number of Monte Carlo samples for UQ
        """
        self.param_distribution = param_distribution
        self.n_samples = n_samples

    def estimate_fidelity_with_uncertainty(
        self,
        fidelity_fn,
        features_a: np.ndarray,
        features_b: np.ndarray,
        random_state: Optional[int] = None
    ) -> UncertaintyEstimate:
        """Estimate fidelity with confidence intervals.

        Args:
            fidelity_fn: Function that computes fidelity given parameters
                         Signature: fidelity_fn(params, features_a, features_b) -> float
            features_a: Features for first cell/state
            features_b: Features for second cell/state
            random_state: Random seed

        Returns:
            UncertaintyEstimate with mean, std, and confidence intervals
        """
        # Sample parameters from posterior
        param_samples = self.param_distribution.sample(
            n_samples=self.n_samples,
            random_state=random_state
        )

        # Compute fidelity for each parameter sample
        fidelity_samples = []
        for params in param_samples:
            fidelity = fidelity_fn(params, features_a, features_b)
            fidelity_samples.append(fidelity)

        fidelity_samples = np.array(fidelity_samples)

        # Compute statistics
        mean_fidelity = np.mean(fidelity_samples)
        std_fidelity = np.std(fidelity_samples)

        # Confidence intervals
        ci_95_lower = np.percentile(fidelity_samples, 2.5)
        ci_95_upper = np.percentile(fidelity_samples, 97.5)
        ci_90_lower = np.percentile(fidelity_samples, 5.0)
        ci_90_upper = np.percentile(fidelity_samples, 95.0)

        # Epistemic uncertainty (from parameter uncertainty)
        epistemic = std_fidelity

        # Aleatoric uncertainty (could be estimated from data noise, set to 0 for now)
        aleatoric = 0.0

        return UncertaintyEstimate(
            mean=mean_fidelity,
            std=std_fidelity,
            confidence_interval_95=(ci_95_lower, ci_95_upper),
            confidence_interval_90=(ci_90_lower, ci_90_upper),
            samples=fidelity_samples,
            epistemic_uncertainty=epistemic,
            aleatoric_uncertainty=aleatoric
        )

    def estimate_classification_confidence(
        self,
        fidelity_estimate: UncertaintyEstimate,
        threshold: float = 0.5
    ) -> float:
        """Estimate confidence in binary classification decision.

        Args:
            fidelity_estimate: Fidelity estimate with uncertainty
            threshold: Classification threshold

        Returns:
            Confidence in [0, 1] - probability that decision is correct
        """
        if fidelity_estimate.samples is None:
            # Use Gaussian approximation
            z_score = abs(fidelity_estimate.mean - threshold) / fidelity_estimate.std
            confidence = 2 * stats.norm.cdf(z_score) - 1  # Convert to [0, 1]
        else:
            # Empirical: fraction of samples on correct side of threshold
            if fidelity_estimate.mean >= threshold:
                confidence = np.mean(fidelity_estimate.samples >= threshold)
            else:
                confidence = np.mean(fidelity_estimate.samples < threshold)

        return float(confidence)


class EnsembleFidelityEstimator:
    """Ensemble-based uncertainty quantification.

    Uses multiple trained models (ensemble) to estimate epistemic uncertainty.
    Simpler than full Bayesian inference but effective in practice.
    """

    def __init__(self, param_distributions: List[BayesianParameterDistribution]):
        """Initialize ensemble estimator.

        Args:
            param_distributions: List of parameter distributions from ensemble members
        """
        self.param_distributions = param_distributions
        self.n_ensemble = len(param_distributions)

    def estimate_fidelity_with_uncertainty(
        self,
        fidelity_fn,
        features_a: np.ndarray,
        features_b: np.ndarray
    ) -> UncertaintyEstimate:
        """Estimate fidelity using ensemble disagreement.

        Args:
            fidelity_fn: Fidelity computation function
            features_a: Features for first cell
            features_b: Features for second cell

        Returns:
            UncertaintyEstimate from ensemble predictions
        """
        fidelity_predictions = []

        for param_dist in self.param_distributions:
            # Use mean parameters from each ensemble member
            fidelity = fidelity_fn(
                param_dist.param_mean,
                features_a,
                features_b
            )
            fidelity_predictions.append(fidelity)

        fidelity_predictions = np.array(fidelity_predictions)

        # Statistics from ensemble
        mean_fidelity = np.mean(fidelity_predictions)
        std_fidelity = np.std(fidelity_predictions)

        # Confidence intervals
        ci_95_lower = np.percentile(fidelity_predictions, 2.5)
        ci_95_upper = np.percentile(fidelity_predictions, 97.5)
        ci_90_lower = np.percentile(fidelity_predictions, 5.0)
        ci_90_upper = np.percentile(fidelity_predictions, 95.0)

        return UncertaintyEstimate(
            mean=mean_fidelity,
            std=std_fidelity,
            confidence_interval_95=(ci_95_lower, ci_95_upper),
            confidence_interval_90=(ci_90_lower, ci_90_upper),
            samples=fidelity_predictions,
            epistemic_uncertainty=std_fidelity,
            aleatoric_uncertainty=0.0
        )


def calibrate_uncertainty(
    predictions: List[float],
    uncertainties: List[float],
    ground_truth: List[float],
    n_bins: int = 10
) -> Dict[str, float]:
    """Calibrate uncertainty estimates against ground truth.

    Checks if predicted confidence intervals contain true values
    at the expected frequency.

    Args:
        predictions: Predicted mean values
        uncertainties: Predicted standard deviations
        ground_truth: True values
        n_bins: Number of bins for calibration plot

    Returns:
        Calibration metrics (expected vs observed coverage)
    """
    predictions = np.array(predictions)
    uncertainties = np.array(uncertainties)
    ground_truth = np.array(ground_truth)

    # Check coverage at different confidence levels
    coverage_levels = [0.68, 0.90, 0.95, 0.99]
    coverage_results = {}

    for level in coverage_levels:
        # Z-score for confidence level
        z = stats.norm.ppf((1 + level) / 2)

        # Compute intervals
        lower = predictions - z * uncertainties
        upper = predictions + z * uncertainties

        # Check coverage
        within_interval = (ground_truth >= lower) & (ground_truth <= upper)
        observed_coverage = np.mean(within_interval)

        coverage_results[f"coverage_{int(level*100)}"] = float(observed_coverage)

    return coverage_results
