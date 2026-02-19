# Token Usage Benchmarking & Cache Optimization

## Context

We want to compare Claude vs Gemini token usage across common bioinformatics prompts, understand where tokens are consumed (especially in multi-step tool-calling loops), and leverage built-in caching to reduce costs. The Streamlit app already has 14 example prompts, per-message usage display, session totals, and an audit logger — but it had three gaps:

1. **No cache metrics** — Claude returns `cache_read_input_tokens` / `cache_creation_input_tokens`; Gemini returns `cached_tokens`. Neither was captured.
2. **Only final-iteration usage** — Both providers loop up to 30 tool-calling iterations but only reported the *last* response's tokens, not cumulative across all iterations.
3. **No automated benchmark mode** — Couldn't batch-run all 14 prompts on both providers and export a comparison report.

Additionally, **Claude prompt caching was not enabled**. Gemini has implicit caching (automatic since May 2025) but we weren't tracking whether it was hitting.

## How Caching Works

### Claude — Automatic Caching (February 2026)

As of February 2026, Claude supports **automatic caching** via a single top-level `cache_control` parameter on the API request. This replaces the previous approach of placing explicit `cache_control` breakpoints on individual content blocks.

- Add `cache_control: {"type": "ephemeral"}` at the **top level** of `messages.create()`
- The system automatically applies the cache breakpoint to the last cacheable block
- In multi-turn conversations and agentic loops, the cache point **moves forward automatically** as the conversation grows
- No need to manually mark system prompts, tool definitions, or message blocks
- 5-minute TTL (refreshed on use), or 1-hour TTL at 2x write cost
- **Reads cost 0.1x** (90% discount); writes cost 1.25x (25% premium)
- Minimum cacheable: 1,024 tokens (Sonnet), 4,096 tokens (Opus/Haiku)
- Up to 4 cache breakpoints supported (automatic uses one slot)
- Can be combined with explicit block-level breakpoints for fine-grained control

**Why automatic caching is ideal for our agentic loop:**

| Request | Content | Cache behavior |
|---------|---------|----------------|
| Iteration 1 | tools + system + user message | Everything written to cache |
| Iteration 2 | tools + system + user + assistant + tool_result | Iteration 1 content read from cache; new content written |
| Iteration 3 | tools + system + user + ... + tool_result_2 | Iterations 1-2 read from cache; new content written |

Each agentic loop iteration benefits from the cache of all prior iterations without any manual bookkeeping.

**Reference:** [Prompt Caching docs](https://platform.claude.com/docs/en/build-with-claude/prompt-caching#automatic-caching)

### Gemini — Implicit (automatic)

- Enabled by default on Gemini 2.5+ models
- Automatic detection of repeated context prefixes
- **Reads cost 0.1x** (90% discount on 2.5 models)
- Minimum: 1,024 tokens (Flash), 2,048 tokens (Pro)
- No developer action needed — just need to *track* hits via `cached_tokens` field

## What Changed

### Step 1: Extended UsageInfo with cache fields

**File:** `ui/streamlit-app/providers/base.py`

Added optional cache fields to `UsageInfo`:
```python
@dataclass
class UsageInfo:
    input_tokens: int
    output_tokens: int
    total_tokens: int
    cache_read_tokens: int = 0       # tokens served from cache
    cache_creation_tokens: int = 0   # tokens written to cache (Claude only)
    iterations: int = 1              # number of agentic loop iterations
```

### Step 2: Cumulative usage tracking across tool-calling iterations

**Files:** `ui/streamlit-app/providers/anthropic_provider.py`, `ui/streamlit-app/providers/gemini_provider.py`

In both providers' agentic loops, we now accumulate usage from every iteration (not just the final one):
- Initialize `cumulative_input = 0, cumulative_output = 0, cumulative_cache_read = 0, cumulative_cache_creation = 0` before the while loop
- After each `messages.create()` / `generate_content()` call, add that response's usage to the running totals
- Return the cumulative totals in the final `ChatResponse`

### Step 3: Enabled Claude automatic prompt caching

**File:** `ui/streamlit-app/providers/anthropic_provider.py`

Added a single top-level `cache_control` parameter to the `messages.create()` call. This uses Claude's automatic caching (released February 2026) which automatically caches everything up to the last cacheable block and moves the cache breakpoint forward as the conversation grows during the agentic loop.

```python
# Automatic caching — one parameter handles tools, system, and conversation
response = self.client.messages.create(
    model=model,
    max_tokens=max_tokens,
    temperature=temperature,
    messages=conversation_history,
    tools=claude_tools,
    system=system_prompt,        # plain string, no content blocks needed
    cache_control={"type": "ephemeral"}  # auto-caches the entire prefix
)
```

This replaced an earlier explicit approach that required wrapping the system prompt in content blocks and mutating tool definitions. The automatic approach is simpler and also caches the growing conversation history between agentic loop iterations (which the explicit approach did not).

**Requires:** `anthropic>=0.75.0` (updated in `requirements.txt`)

### Step 4: Extracted cache fields from both provider responses

**File:** `ui/streamlit-app/providers/anthropic_provider.py`
```python
# Claude usage object has: input_tokens, output_tokens,
# cache_read_input_tokens, cache_creation_input_tokens
cache_read = getattr(response.usage, 'cache_read_input_tokens', 0)
cache_creation = getattr(response.usage, 'cache_creation_input_tokens', 0)
```

**File:** `ui/streamlit-app/providers/gemini_provider.py`
```python
# Gemini usage_metadata has: prompt_token_count, candidates_token_count,
# total_token_count, cached_content_token_count (implicit cache)
cache_read = getattr(usage, 'cached_content_token_count', 0)
```

### Step 5: UI cache metrics display

**File:** `ui/streamlit-app/app.py`

Updated the per-message "Token Usage" expander to:
- Show a 4th column "Cached" when cache reads are present
- Display iterations count, cache write tokens, and cache hit rate below metrics
- Use cache-aware cost calculation:
  - Claude: cached reads at 0.1x rate, cache writes at 1.25x rate, uncached at 1.0x
  - Gemini: cached reads at 0.1x rate, uncached at 1.0x

### Step 6: Automated benchmark mode

**New file:** `ui/streamlit-app/utils/benchmark.py`

A `TokenBenchmark` class that:
1. Takes a dict of prompt names/texts and a list of providers/models to test
2. Runs each prompt on each provider sequentially (cold run)
3. Runs each prompt again (warm run) to measure cache hits
4. Collects per-run data: prompt_name, provider, model, input_tokens, output_tokens, cache_read_tokens, cache_creation_tokens, iterations, duration_ms, estimated_cost
5. Returns results as `BenchmarkResult` dataclass instances

**File:** `ui/streamlit-app/app.py`

Added a "Token Benchmark" section in the sidebar:
- Multi-select for prompts (default: all 14)
- Checkboxes for providers (Claude, Gemini — only if API keys set)
- Model selectors for each provider
- "Run Benchmark" button with progress bar
- Results displayed as a Streamlit dataframe with cache_hit_rate column
- Bar chart of token usage by prompt (grouped by provider + run type)
- "Export CSV" download button

### Step 7: Benchmark results in audit log

**File:** `ui/streamlit-app/utils/audit_logger.py`

Added `log_benchmark_run()` method that logs the full benchmark results to GCP Cloud Logging with a `benchmark` event type, including summary stats (total tokens, cost, cache reads, error count).

## Files Modified (Summary)

| File | Change |
|------|--------|
| `providers/base.py` | Add cache fields + iterations to UsageInfo |
| `providers/anthropic_provider.py` | Enable automatic caching, cumulative usage tracking, extract cache fields |
| `providers/gemini_provider.py` | Cumulative usage tracking, extract cached_content_token_count |
| `requirements.txt` | Bump `anthropic>=0.75.0` for automatic caching support |
| `app.py` | Display cache metrics, incremental benchmark UI, cache-aware cost calc |
| `utils/benchmark.py` | **New** — benchmark runner class |
| `utils/audit_logger.py` | Add log_benchmark_run() method |
| `deploy.sh`, `deploy_now.sh` | Increase Cloud Run memory from 1Gi to 2Gi |

## Test Protocol

### Prompts to use

All 14 existing prompts from `EXAMPLE_PROMPTS` in `utils/mcp_config.py`:

1. Warm Up Servers (lightweight, tool discovery)
2. Spatial Analysis (single tool call)
3. Multi-omics Integration (single tool call, larger data)
4. Genomic QC (single tool call)
5. Pathway Enrichment (single tool call)
6. Complete PatientOne Workflow (multi-step, 3+ tool calls)
7. Batch Correction (single tool call)
8. Predict Treatment Response (multi-step)
9. Immunotherapy Prediction (single tool call)
10. Drug Screening (multi-step, 3 tool calls)
11. Quantum Cell Type Fidelity (single tool call)
12. Immune Evasion Detection (single tool call)
13. TLS Analysis (single tool call)
14. Quantum + GEARS Validation (multi-step, 2+ tool calls)

### Run structure

For each provider (Claude Sonnet 4.6, Gemini 3 Flash):
- **Run A (cold):** Fresh session, run each prompt once — measures baseline tokens + cache creation
- **Run B (warm):** Same session, run each prompt again — measures cache hits
- **Run C (repeat):** Run 3 multi-step prompts (#6, #8, #14) a third time — confirms cache stability

### What to measure

| Metric | Why |
|--------|-----|
| `input_tokens` | Baseline cost driver |
| `output_tokens` | Response verbosity |
| `cache_read_tokens` | Cache effectiveness |
| `cache_creation_tokens` | Cache overhead (Claude) |
| `iterations` | Tool-calling depth |
| `duration_ms` | Latency impact |
| `estimated_cost` | Bottom-line comparison |

### Expected insights

- **Automatic caching across agentic loop iterations**: With automatic caching, each iteration of the tool-calling loop caches the entire prefix (tools + system + conversation so far). Iterations 2-N read the growing prefix from cache at 0.1x cost. This is especially impactful for multi-step prompts with 3+ iterations.
- **System prompt + tools caching**: With 5-6 servers, tool definitions are likely 2,000-5,000 tokens. These are cached automatically and reused across all iterations and across prompts within the 5-minute TTL window.
- **Multi-step prompts**: Prompts like #6 (Complete Workflow) make 3+ API calls in the loop. Cumulative tracking will reveal the true token cost vs. the previously reported final-iteration-only number.
- **Gemini implicit caching**: May already be providing savings we weren't seeing because `cached_content_token_count` wasn't tracked.
- **Prompt pattern guidance**: If caching is effective, guidance is: "Keep the same servers selected across prompts in a session" and "Run related analyses together rather than switching contexts."

## Verification

1. Run a single prompt on Claude — verify `cache_creation_tokens > 0` on first call, `cache_read_tokens > 0` on second call
2. Run same prompt on Gemini — verify `cached_tokens` field appears (may be 0 on first, > 0 on repeat)
3. Run the full benchmark — verify CSV export contains all expected columns
4. Compare `iterations` count against tool_calls_metadata length — should match
5. Verify cumulative tokens > final-iteration tokens for multi-step prompts
6. Check GCP Cloud Logging for benchmark event entries
