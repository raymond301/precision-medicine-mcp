"""Token usage benchmarking for comparing Claude vs Gemini across prompts.

Runs example prompts on one or more providers and collects per-run token
usage, cache metrics, iteration counts, latency, and estimated cost.
"""

import time
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Optional

from providers.base import ChatMessage, UsageInfo


@dataclass
class BenchmarkResult:
    """Result of a single benchmark run."""
    prompt_name: str
    provider: str
    model: str
    run_type: str  # "cold", "warm", or "repeat"
    input_tokens: int = 0
    output_tokens: int = 0
    total_tokens: int = 0
    cache_read_tokens: int = 0
    cache_creation_tokens: int = 0
    iterations: int = 1
    duration_ms: float = 0.0
    estimated_cost: float = 0.0
    error: Optional[str] = None

    def to_dict(self) -> Dict:
        return asdict(self)


class TokenBenchmark:
    """Batch-runs prompts on providers and collects token usage data.

    Usage:
        benchmark = TokenBenchmark(
            prompts={"Spatial Analysis": "Use the spatial_autocorrelation ..."},
            providers=[("claude", provider_instance, "claude-sonnet-4-6")],
            mcp_servers=server_configs
        )
        results = benchmark.run_cold()
        results += benchmark.run_warm()
    """

    def __init__(
        self,
        prompts: Dict[str, str],
        providers: List[tuple],  # [(name, instance, model), ...]
        mcp_servers: List[Dict],
        max_tokens: int = 4096,
        uploaded_files: Optional[Dict] = None,
    ):
        """Initialize benchmark.

        Args:
            prompts: Dict of prompt_name -> prompt_text
            providers: List of (provider_key, provider_instance, model_id)
            mcp_servers: MCP server configurations
            max_tokens: Max tokens per response
            uploaded_files: Optional uploaded files dict
        """
        self.prompts = prompts
        self.providers = providers
        self.mcp_servers = mcp_servers
        self.max_tokens = max_tokens
        self.uploaded_files = uploaded_files

    def _estimate_cost(self, usage: UsageInfo, provider_key: str) -> float:
        """Estimate cost from usage info with cache-aware pricing."""
        cache_read = usage.cache_read_tokens
        cache_creation = usage.cache_creation_tokens
        raw_input = usage.input_tokens

        if provider_key == "gemini":
            uncached = max(0, raw_input - cache_read)
            input_cost = uncached * 0.00015 / 1000 + cache_read * 0.000015 / 1000
            output_cost = usage.output_tokens * 0.0006 / 1000
        else:
            uncached = max(0, raw_input - cache_read - cache_creation)
            input_cost = (uncached * 0.003 / 1000
                          + cache_read * 0.0003 / 1000
                          + cache_creation * 0.00375 / 1000)
            output_cost = usage.output_tokens * 0.015 / 1000

        return input_cost + output_cost

    def _run_single(
        self,
        prompt_name: str,
        prompt_text: str,
        provider_key: str,
        provider_instance,
        model: str,
        run_type: str,
        progress_callback=None,
    ) -> BenchmarkResult:
        """Run a single prompt on a single provider.

        Args:
            prompt_name: Display name of the prompt
            prompt_text: Full prompt text
            provider_key: Provider key ("claude" or "gemini")
            provider_instance: Provider instance
            model: Model ID
            run_type: "cold", "warm", or "repeat"
            progress_callback: Optional callable(status_str)

        Returns:
            BenchmarkResult
        """
        if progress_callback:
            progress_callback(f"{run_type.title()} | {provider_key} | {prompt_name}")

        try:
            messages = [ChatMessage(role="user", content=prompt_text)]
            start = time.time()

            chat_response = provider_instance.send_message(
                messages=messages,
                mcp_servers=self.mcp_servers,
                model=model,
                max_tokens=self.max_tokens,
                uploaded_files=self.uploaded_files,
            )

            duration_ms = (time.time() - start) * 1000
            usage = chat_response.usage

            if usage:
                return BenchmarkResult(
                    prompt_name=prompt_name,
                    provider=provider_key,
                    model=model,
                    run_type=run_type,
                    input_tokens=usage.input_tokens,
                    output_tokens=usage.output_tokens,
                    total_tokens=usage.total_tokens,
                    cache_read_tokens=usage.cache_read_tokens,
                    cache_creation_tokens=usage.cache_creation_tokens,
                    iterations=usage.iterations,
                    duration_ms=duration_ms,
                    estimated_cost=self._estimate_cost(usage, provider_key),
                )
            else:
                return BenchmarkResult(
                    prompt_name=prompt_name,
                    provider=provider_key,
                    model=model,
                    run_type=run_type,
                    duration_ms=duration_ms,
                    error="No usage info returned",
                )

        except Exception as e:
            return BenchmarkResult(
                prompt_name=prompt_name,
                provider=provider_key,
                model=model,
                run_type=run_type,
                error=str(e),
            )

    def run_cold(self, progress_callback=None) -> List[BenchmarkResult]:
        """Run each prompt once on each provider (fresh, no cache).

        Args:
            progress_callback: Optional callable(status_str)

        Returns:
            List of BenchmarkResult
        """
        results = []
        for prompt_name, prompt_text in self.prompts.items():
            for provider_key, provider_instance, model in self.providers:
                result = self._run_single(
                    prompt_name, prompt_text,
                    provider_key, provider_instance, model,
                    run_type="cold",
                    progress_callback=progress_callback,
                )
                results.append(result)
        return results

    def run_warm(self, progress_callback=None) -> List[BenchmarkResult]:
        """Run each prompt again (same session, should hit cache).

        Args:
            progress_callback: Optional callable(status_str)

        Returns:
            List of BenchmarkResult
        """
        results = []
        for prompt_name, prompt_text in self.prompts.items():
            for provider_key, provider_instance, model in self.providers:
                result = self._run_single(
                    prompt_name, prompt_text,
                    provider_key, provider_instance, model,
                    run_type="warm",
                    progress_callback=progress_callback,
                )
                results.append(result)
        return results

    def run_repeat(
        self, prompt_names: List[str], progress_callback=None
    ) -> List[BenchmarkResult]:
        """Run selected prompts a third time (cache stability check).

        Args:
            prompt_names: Subset of prompt names to re-run
            progress_callback: Optional callable(status_str)

        Returns:
            List of BenchmarkResult
        """
        results = []
        for prompt_name in prompt_names:
            prompt_text = self.prompts.get(prompt_name)
            if not prompt_text:
                continue
            for provider_key, provider_instance, model in self.providers:
                result = self._run_single(
                    prompt_name, prompt_text,
                    provider_key, provider_instance, model,
                    run_type="repeat",
                    progress_callback=progress_callback,
                )
                results.append(result)
        return results

    @staticmethod
    def results_to_csv_string(results: List[BenchmarkResult]) -> str:
        """Convert results to CSV string for download.

        Args:
            results: List of BenchmarkResult

        Returns:
            CSV string
        """
        if not results:
            return ""

        headers = [
            "prompt_name", "provider", "model", "run_type",
            "input_tokens", "output_tokens", "total_tokens",
            "cache_read_tokens", "cache_creation_tokens",
            "iterations", "duration_ms", "estimated_cost", "error"
        ]

        lines = [",".join(headers)]
        for r in results:
            d = r.to_dict()
            row = []
            for h in headers:
                val = d.get(h, "")
                # Quote strings that might contain commas
                if isinstance(val, str) and ("," in val or '"' in val):
                    val = '"' + val.replace('"', '""') + '"'
                row.append(str(val) if val is not None else "")
            lines.append(",".join(row))

        return "\n".join(lines)
