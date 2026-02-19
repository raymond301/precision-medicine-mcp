# Token Benchmark Findings — 2026-02-19

## What We Tested

Ran all 14 example prompts across Claude Sonnet 4.6 and Gemini 3 Flash in 3 phases:
- **Cold** — fresh cache, baseline token usage
- **Warm** — same session, measures cache hits
- **Repeat** — 3 multi-step prompts a 3rd time, confirms cache stability

62 total runs (28 cold + 28 warm + 6 repeat). 47 succeeded, 15 timed out (90s limit).

## Key Findings

### Claude automatic caching works

- System prompt + tool definitions (~14,600 tokens) cached after the first call and reused by every subsequent call at 0.1x cost
- Cache grows across agentic loop iterations: multi-step prompts cached 80K–96K tokens of conversation prefix
- Cache is **stable** — warm and repeat phases showed identical cache reads for the same prompts
- First prompt (Warm Up Servers) cost $0.074 cold vs $0.023 warm — **69% cost reduction** from caching

### Gemini implicit caching is unreliable

- Only triggered on 6 of 32 successful runs (19%)
- When it hits, savings are significant (TLS Analysis: 162K tokens cached, 53% of input)
- Cannot be relied on for cost planning — no developer control

### Cost and speed

| | Claude Sonnet 4.6 | Gemini 3 Flash | Ratio |
|---|---|---|---|
| Per-prompt avg cost | ~$0.04 | ~$0.005 | Claude is **8x** more |
| Simple prompt speed | ~17s | ~9s | Gemini is **2x** faster |
| Complex prompt speed | ~60s | ~45s | Gemini is **1.3x** faster |

Gemini Flash is cheaper due to base pricing ($0.15/M vs $3/M input). Claude's caching reduces cost but doesn't close the 20x pricing gap.

### Iteration variability

The same prompt can take 1–7 agentic loop iterations across runs. This causes unpredictable latency and occasional timeouts. Not a caching issue — it's the model deciding how many tool calls to make.

### Claude's `input_tokens` field is misleading

Claude reports `input_tokens` as only the uncached portion (3–9 tokens). Real input = `input_tokens + cache_read_tokens + cache_creation_tokens`. Gemini reports the full count. The benchmark UI now computes `real_input_tokens` to normalize this.

## Implementation Notes

- **Benchmark split into 3 phases** to avoid Streamlit websocket timeouts. Each phase runs ~14 prompts with one prompt per `st.rerun()` cycle.
- **Per-prompt timeout** set to 180s (increased from 90s after first run showed too many timeouts).
- **Per-prompt audit logging** — each result is written to GCP Cloud Logging as it completes, so partial results survive interruptions.
- **Claude automatic caching** enabled via a single `cache_control={"type": "ephemeral"}` parameter on `messages.create()`. No manual breakpoints needed.

## Recommendations

1. **Use Claude for complex multi-step workflows** — caching compounds across iterations, results are more consistent
2. **Use Gemini for high-volume simple queries** — 8x cheaper, 2x faster on single-tool prompts
3. **Keep servers stable within a session** — switching servers invalidates the ~14,600-token tool definition cache
4. **Run related analyses together** — cache TTL is 5 minutes, so sequential prompts benefit from the warm cache
