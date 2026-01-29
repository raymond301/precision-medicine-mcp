# Token Usage Analysis - Why Costs Are Higher Than Expected

## The Problem

A simple query like "which servers are here?" used **25K tokens** instead of the expected 1-2K tokens.

## Root Cause: Tool Schema Overhead

Every Claude API call includes **ALL tool schemas** from selected servers, regardless of whether they'll be used.

### Default Configuration (3 servers)

| Server | Tools | Estimated Schema Size |
|--------|-------|----------------------|
| spatialtools | 10 tools | ~5K tokens |
| multiomics | 9 tools | ~4.5K tokens |
| fgbio | 4 tools | ~2K tokens |
| **TOTAL** | **23 tools** | **~12K tokens** |

### Token Breakdown Per Query

| Component | Tokens | % of Total |
|-----------|--------|------------|
| **Tool schemas** (23 tools) | 10-12K | 48% |
| **System prompt** (server descriptions) | 2-3K | 12% |
| **Conversation history** (previous turns) | 8-10K | 40% |
| **User query** | 100-500 | <2% |
| **Response** | 500-2K | <8% |
| **TOTAL per query** | **20-27K** | **100%** |

## Why Token Usage Accumulates

### First Query in Session
```
Tool schemas:     12K
System prompt:     2K
User query:      0.2K
Response:        0.5K
-------------------
TOTAL:          ~15K tokens
```

### Second Query (conversation history included)
```
Tool schemas:     12K
System prompt:     2K
Previous msg 1:    1K  (user)
Previous msg 2:   0.5K (assistant)
User query:      0.2K
Response:        0.5K
-------------------
TOTAL:          ~16K tokens
```

### Fifth Query (lots of history)
```
Tool schemas:     12K
System prompt:     2K
Previous msgs:   8-10K  (4 exchanges)
User query:      0.2K
Response:        0.5K
-------------------
TOTAL:          ~23K tokens
```

## Actual Token Usage by Query Type

### Simple Informational Query
"Which servers are here?"
- **First query**: 15K tokens
- **After 2-3 exchanges**: 20-25K tokens

### Medium Analysis Query
"Calculate spatial autocorrelation for gene CD8A"
- **First query**: 15-18K tokens (tool execution adds conversation)
- **After tool call**: 25-30K tokens (tool response in history)
- **Total for 1 analysis**: ~40-50K tokens

### Complex Multi-Step Query
"Load data, run QC, identify spatially variable genes, visualize"
- **Initial query**: 15K tokens
- **Tool call 1**: 20K tokens
- **Tool call 2**: 25K tokens
- **Tool call 3**: 30K tokens
- **Final response**: 35K tokens
- **Total for 1 workflow**: ~125K+ tokens (exceeds 50K student limit!)

## Impact on Cost Forecast

### Original Estimates (WRONG ❌)
- Simple query: 1-2K tokens → $0.009-0.018
- Medium query: 5-8K tokens → $0.045-0.072
- Complex query: 15-25K tokens → $0.135-0.225

### Actual Usage (CORRECT ✅)
- Simple query: 15-25K tokens → $0.135-0.225
- Medium query: 40-50K tokens → $0.36-0.45
- Complex query: 50K+ tokens → $0.45+ (hits student limit!)

## Revised Per-Student Session Estimates

| Session Activity | Queries | Tokens | Cost |
|-----------------|---------|--------|------|
| **Week 1-2** (Mock mode) | Unlimited | 0 | $0.00 |
| **Week 3** (First real analysis) | 2-3 queries | 40-50K | $0.36-0.45 |
| **Week 4** (Tool integration) | 2-3 queries | 45-50K | $0.41-0.45 |
| **Week 5** (Advanced workflows) | 1-2 complex | 50K (hit limit) | $0.45 |
| **Week 6** (Final project) | Multiple resets | 50K × 2-3 | $0.90-1.35 |

## Student Limit Protection

The 50K token limit per session means:
- **Maximum cost per session**: ~$0.45
- **Students will hit limit** after 2-3 queries (not 30-40 as estimated)
- **Requires frequent resets** during active sessions
- **This is actually GOOD** - prevents runaway costs!

## Implications for 6-Week Program

### Per Student (4 students)
- Week 1-2: $0 (mock mode)
- Week 3: $0.45 (1 session + homework)
- Week 4: $0.90 (likely hit limit twice)
- Week 5: $0.90 (multiple resets for complex work)
- Week 6: $1.35 (final project with iterations)
- **Total per student**: ~$3.60

### Instructor (1)
- Testing and demos: Higher usage, ~$2-3
- Likely hits limits frequently: $5-8
- **Total instructor**: ~$5-8

### Total Program Cost
- **4 students**: 4 × $3.60 = $14.40
- **1 instructor**: $5-8
- **Buffer (20%)**: $4
- **TOTAL**: **$23-26**

## Recommendations

### Budget
- **Minimum**: $25 (covers expected usage)
- **Comfortable**: $30-35 (allows extra exploration)
- **Safe**: $40 (handles unexpected resets)

### For Students
1. **Expect to hit 50K limit** after 2-3 queries
2. **Clear conversation frequently** - this is normal!
3. **Don't panic at reset** - limit is working as designed
4. **Use mock mode for practice** - saves tokens

### For Instructor
1. **Set expectations** - explain why limits hit quickly
2. **Demo the reset process** - make it feel normal
3. **Encourage focused queries** - less back-and-forth
4. **Consider reducing default servers** - fewer tool schemas

## Optimization Options (Future)

### Reduce Tool Overhead
Currently: 23 tools (3 servers) = ~12K tokens per query

**Option 1**: Start with 1 server (spatialtools only)
- 10 tools = ~5K tokens
- Saves ~7K per query
- ~30% cost reduction

**Option 2**: Dynamic tool loading (only relevant tools)
- Requires code changes
- Could reduce to 2-3 tools per query
- ~80% cost reduction

**Option 3**: Shorter conversation history
- Keep only last 2 exchanges
- Saves ~5-10K tokens after several queries
- ~20-40% cost reduction

## Conclusion

**The 50K token limit is appropriate** - it caps costs at $0.45/session even with high tool overhead.

**Students will need to reset 2-3 times** per active session, but this is acceptable and prevents accidental overspending.

**Budget $25-30 for the program** instead of the original $5-10 estimate.
