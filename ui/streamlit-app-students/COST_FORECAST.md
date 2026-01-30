# Student App - Token Usage & Cost Forecast (REVISED)

**Study Group**: 4 students + 1 instructor  
**Duration**: 6 weeks (2 sessions/month + homework)  
**Pricing**: Claude Sonnet 4.5 (~$0.009/1K tokens average)  

⚠️ **IMPORTANT**: Token usage is ~10x higher than initially estimated due to tool schema overhead. See TOKEN_ANALYSIS.md for details.

## Weekly Breakdown (Revised)

| Week | Session Focus | Mode | Example Activity | Tokens/Student | Cost/Student | Weekly Total (4 students + instructor) |
|------|---------------|------|-----------------|----------------|--------------|----------------------------------------|
| **1** | Understanding MCP | Mock | "List tools", "What does spatialtools do?" | 0 (mock) | $0.00 | **$0.00** |
| **2** | Writing Prompts | Mock | "Explain HAllA", "Show quantum example" | 0 (mock) | $0.00 | **$0.00** |
| **3** | First Real Analysis | Real | 2-3 simple queries (hit 50K limit once) | 50K | $0.45 | **$2.70** |
| **4** | Tool Integration | Real | 2-3 queries + homework (2 resets) | 50K × 2 | $0.90 | **$5.40** |
| **5** | Advanced Workflows | Real | Complex multi-step (2 resets) | 50K × 2 | $0.90 | **$5.40** |
| **6** | Final Project | Real | Iterations + debugging (3 resets) | 50K × 3 | $1.35 | **$8.10** |

### Cost Breakdown by Week

**Weeks 1-2 (Mock Mode)**:
- Free unlimited practice
- Subtotal: **$0.00**

**Weeks 3-6 (Real Mode)**:
- Week 3: $2.70 (students learning, 1 reset each)
- Week 4: $5.40 (more complex, 2 resets each)
- Week 5: $5.40 (multi-step workflows, 2 resets each)
- Week 6: $8.10 (final projects, 3 resets each)
- Subtotal: **$21.60**

## Total Program Cost (Revised)

| Category | Tokens | Cost |
|----------|--------|------|
| **Students (4)** | ~750K (50K × 3.75 avg resets × 4 students) | $14.40 |
| **Instructor (1)** | ~500K (higher usage, more testing) | $7.00 |
| **Buffer (20%)** | - | $4.30 |
| **TOTAL** | 1.25M | **$25-26** |

## Per-Student Safety Limits

Each student has:
- **50K tokens per session** (enforced by app)
- **50 requests per session** (enforced by app)
- **Easy reset**: Clear conversation to start fresh

### Token Usage Reality Check (REVISED)

⚠️ **Token usage is ~10x higher than expected** due to tool schema overhead (~12K tokens) sent with EVERY query.

| Activity | Actual Tokens | Queries Before Limit |
|----------|---------------|---------------------|
| Simple queries (no tool calls) | 15-25K each | 2-3 queries |
| Medium analysis (1-2 tool calls) | 40-50K each | 1-2 queries |
| Complex workflow (multi-step) | 50K+ | Hits limit immediately |
| Accidentally hit limit | 50K | 🛑 Auto-stops at $0.45 |

## Example Prompts by Token Usage (REVISED)

### Simple Queries (~15-25K tokens each)
Includes ~12K tool schema overhead + conversation history
```
"What tools are available in spatialtools?"        → 15-20K tokens
"Explain what quantum fidelity measures"           → 15-20K tokens
"List the steps in a HAllA analysis"               → 16-22K tokens
```
**Impact**: 2-3 queries before hitting 50K limit

### Medium Analysis (~40-50K tokens total)
Single analysis with tool execution
```
"Calculate spatial autocorrelation for gene CD8A"  → 40-45K tokens
  (includes tool call and result)

"Run differential expression between regions"       → 42-48K tokens
  (includes data loading + analysis)
```
**Impact**: 1-2 queries before hitting 50K limit, requires reset

### Complex Workflows (50K+ tokens - EXCEEDS LIMIT)
Multi-step processes hit limit before completing
```
"Load Visium data, perform QC, identify spatially variable genes,
 then cluster them by expression pattern"
  → Hits 50K limit partway through, needs reset to continue

"Complete spatial analysis: load → normalize →
 spatial autocorrelation → cell type deconvolution → visualize"
  → Requires 2-3 resets (total ~100-150K tokens)
```
**Impact**: Requires multiple resets per workflow

## Budget Recommendations (REVISED)

### Minimum Budget
- **$25-30** - Covers expected usage with small buffer
- Assumes ~3-4 resets per student across 4 weeks
- Good for: Structured curriculum, students follow assignments closely

### Comfortable Budget
- **$35-40** - Allows extra exploration and mistakes
- Students can reset freely and experiment
- Good for: Encouraging hands-on learning and debugging

### Safe Budget
- **$50** - Generous buffer for unexpected usage
- Covers instructor demos, student retries, and edge cases
- Good for: First-time teaching, uncertain group dynamics, complex projects

## Cost Control Features

The student app has these protections:

1. **Per-Session Limits** (CRITICAL)
   - 50K tokens max (~$0.45 per session)
   - Auto-stops before runaway costs
   - **Students will hit this after 2-3 queries** - this is expected!

2. **Real-Time Tracking**
   - Visible usage counter in sidebar
   - Warnings at 80% (40K tokens) - usually means 1-2 queries left
   - Shows cost estimate

3. **Easy Reset** (Use frequently!)
   - Clear conversation = fresh 50K tokens
   - No penalty, just restart
   - **Expected 3-4 resets per student** across the program

4. **Mock Mode First** (Essential!)
   - Weeks 1-2 completely free
   - Students learn interface without ANY token usage
   - Practice writing effective prompts before using real mode

## Worst-Case Scenario (REVISED)

**If students hit limits frequently (5 resets each):**
- 4 students × 5 resets × 50K × $0.009 = **$9.00**
- Instructor (heavy usage): 10 resets × 50K × $0.009 = **$4.50**
- Week 6 intensive final projects: +$8-10
- **Total worst case: ~$22-24**

Even in worst case, the **50K limit prevents runaway costs** - maximum $0.45 per reset.

## Key Insights

### Why Token Usage Is High
See **TOKEN_ANALYSIS.md** for full details:
- **Tool schema overhead**: ~12K tokens (23 tools) sent with EVERY query
- **Conversation history**: Grows with each exchange (~10K by 5th query)
- **System prompt**: ~2-3K tokens describing servers
- **Actual query + response**: Only ~1-3K tokens

### What This Means
- Students will hit 50K limit after **2-3 queries** (not 30-40)
- **Resets are normal and expected** during active learning
- Complex workflows need **multiple resets** to complete
- **50K limit is appropriate** - caps cost at $0.45/session

## Recommendations (REVISED)  

✅ **Budget $30-35** - realistic for 4 students + instructor  
✅ **Use mock mode weeks 1-2** - build confidence for free (crucial!)  
✅ **Set expectations** - explain students will hit limits after 2-3 queries  
✅ **Normalize resets** - "Clear conversation" is a normal part of the workflow  
✅ **Start with simple queries** - build up to complex workflows gradually  
✅ **Consider fewer servers** - reduce from 3 default to 1-2 (saves ~7K tokens/query)  

With 4 students following the curriculum, expect **$25-30 total cost** for the entire 6-week program.
