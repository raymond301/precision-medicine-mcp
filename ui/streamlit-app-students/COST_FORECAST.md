# Student App - Token Usage & Cost Forecast

**Study Group**: 4 students + 1 instructor
**Duration**: 6 weeks (2 sessions/month + homework)
**Pricing**: Claude Sonnet 4.5 (~$0.009/1K tokens average)

## Weekly Breakdown

| Week | Session Focus | Mode | Example Prompts | Tokens/Student | Tokens/Instructor | Session Cost | Homework Cost | Weekly Total |
|------|---------------|------|-----------------|----------------|-------------------|--------------|---------------|--------------|
| **1** | Understanding MCP | Mock | "List all available tools", "What does spatialtools do?" | 0 (mock) | 0 (mock) | $0.00 | $0.00 | **$0.00** |
| **2** | Writing Prompts | Mock | "Explain HAllA", "Show me quantum fidelity example" | 0 (mock) | 0 (mock) | $0.00 | $0.00 | **$0.00** |
| **3** | First Real Analysis | Real | "Calculate spatial autocorrelation for CD8A", "Run HAllA on test data" | 10K | 20K | $0.54 | $0.54 | **$1.08** |
| **4** | Tool Integration | Real | "Identify spatially variable genes", "Integrate proteomics + metabolomics" | 15K | 25K | $0.77 | $0.68 | **$1.45** |
| **5** | Advanced Workflows | Real | "Multi-step: load → QC → analyze → visualize", "Batch process 3 samples" | 20K | 30K | $0.99 | $0.90 | **$1.89** |
| **6** | Final Project | Real | "Complete analysis pipeline for ovarian cancer dataset" | 25K | 35K | $1.22 | $1.13 | **$2.35** |

### Cost Breakdown by Week

**Weeks 1-2 (Mock Mode)**:
- Session cost: $0
- Homework cost: $0
- Subtotal: **$0.00**

**Weeks 3-6 (Real Mode)**:
- Session costs: $0.54 + $0.77 + $0.99 + $1.22 = **$3.52**
- Homework costs: $0.54 + $0.68 + $0.90 + $1.13 = **$3.25**
- Subtotal: **$6.77**

## Total Program Cost

| Category | Tokens | Cost |
|----------|--------|------|
| **Students (4)** | 280K | $2.52 |
| **Instructor (1)** | 110K | $0.99 |
| **Buffer (20%)** | - | $0.70 |
| **TOTAL** | 390K | **$4.21** |

## Per-Student Safety Limits

Each student has:
- **50K tokens per session** (enforced by app)
- **50 requests per session** (enforced by app)
- **Easy reset**: Clear conversation to start fresh

### Token Usage Reality Check

| Activity | Estimated Tokens | Within Limit? |
|----------|-----------------|---------------|
| 1 hour session (10-15 queries) | 10-15K | ✅ Yes (30% of limit) |
| Homework (5-10 queries) | 5-10K | ✅ Yes (20% of limit) |
| Complex final project | 25-30K | ✅ Yes (60% of limit) |
| Accidentally hit limit | 50K | 🛑 Auto-stops, cost capped |

## Example Prompts by Complexity

### Simple Queries (~1-2K tokens)
```
"What tools are available in spatialtools?"
"Explain what quantum fidelity measures"
"List the steps in a HAllA analysis"
```

### Medium Queries (~5-8K tokens)
```
"Calculate spatial autocorrelation for gene CD8A in my dataset"
"Run differential expression between tumor and normal regions"
"Find associations between proteomics and metabolomics using HAllA"
```

### Complex Queries (~15-25K tokens)
```
"Load my Visium data, perform QC, identify spatially variable genes,
 then cluster them by expression pattern"

"Integrate RNA-seq, proteomics, and metabolomics data for
 ovarian cancer samples, run HAllA, and interpret results"

"Complete spatial analysis workflow: load → normalize →
 spatial autocorrelation → cell type deconvolution → visualize"
```

## Budget Recommendations

### Minimum Budget
- **$5** - Covers expected usage with small buffer
- Assumes students stay mostly within session limits
- Good for: Structured curriculum following the plan above

### Comfortable Budget
- **$10** - Allows extra exploration and mistakes
- Students can reset and retry failed queries
- Good for: Encouraging experimentation

### Safe Budget
- **$20** - Generous buffer for unexpected usage
- Covers if students hit limits multiple times
- Good for: First-time teaching, uncertain group dynamics

## Cost Control Features

The student app has these protections:

1. **Per-Session Limits**
   - 50K tokens max (~$0.45 per session)
   - Auto-stops before runaway costs

2. **Real-Time Tracking**
   - Visible usage counter in sidebar
   - Warnings at 80% (40K tokens)

3. **Easy Reset**
   - Clear conversation = fresh 50K tokens
   - No penalty, just restart

4. **Mock Mode First**
   - Weeks 1-2 completely free
   - Students learn interface without cost

## Worst-Case Scenario

**If every student hits limit every week:**
- 4 students × 4 weeks × 50K tokens × $0.009 = **$7.20**
- Plus instructor: 4 weeks × 50K × $0.009 = **$1.80**
- **Total worst case: ~$9.00**

Even in the worst case, cost is under $10 due to safety limits.

## Recommendations

✅ **Start with $10 budget** - comfortable margin
✅ **Use mock mode weeks 1-2** - build confidence for free
✅ **Monitor week 3 usage** - adjust if needed
✅ **Remind students about limits** - they're protective, not punitive

With 4 students following the curriculum, expect **$4-7 total cost** for the entire 6-week program.
