# Student MCP Chat App - Quick Start Guide

Welcome to the **Student MCP Chat App**! This app is designed for the bioinformatics study group with **built-in safety guardrails** to prevent accidental costs.

## 🎓 What Makes This App Special?

This is **not** the full production app. It has safety features:

✅ **Token Limits**: Max 50,000 tokens per session (~$1.50 max cost)  
✅ **Request Limits**: Max 50 requests per session  
✅ **Usage Tracking**: See your usage in real-time  
✅ **Auto-Warnings**: Alerts at 80% usage  
✅ **Easy Reset**: Clear conversation to reset limits  

**Your instructor's app has NO limits** - this lets them test freely while you learn safely.

---

## 🚀 Quick Start (5 Minutes)

### Option 1: Use the Deployed App (Recommended)

Your instructor will provide:
```
Student App URL: https://streamlit-mcp-chat-students-XXXXX.run.app
```

Just click and start using! No setup needed.

### Option 2: Run Locally

```bash
# 1. Clone the repo
git clone https://github.com/yourusername/spatial-mcp
cd spatial-mcp/ui/streamlit-app-students

# 2. Copy environment file
cp .env.example .env

# 3. Add API keys (instructor will provide)
# Edit .env and paste the keys

# 4. Install dependencies
python3 -m pip install -r ../streamlit-app/requirements.txt

# 5. Run the app
streamlit run app.py
```

---

## 📊 Understanding the Safety Limits

### Session Limits

| Limit | Value | Why? |
|-------|-------|------|
| **Tokens per request** | 4,096 | Prevents single huge queries |
| **Tokens per session** | 50,000 | ~$1.50 max cost per session |
| **Requests per session** | 50 | Prevents runaway loops |

### What Happens When You Hit a Limit?

You'll see a friendly message:
```
🛑 Session Token Limit Reached
Great job exploring! You've used 50,000 tokens this session.

To continue:
- Click "Clear conversation" in sidebar
- Or refresh the page

Cost estimate: $1.50
```

**No money is wasted** - limits stop you before costs grow.

---

## 🎯 Study Group Workflow

### **Weeks 1-2: Mock Mode (FREE!)**

```bash
# In .env
USE_MOCK_MCP=true
```

**What you get:**
- ✅ Practice with instant fake responses
- ✅ Learn the interface
- ✅ Test prompts without any cost
- ✅ Unlimited usage

**Best for:**
- Getting comfortable with the UI
- Writing effective prompts
- Understanding MCP servers
- Homework assignments 1-2

### **Weeks 3+: Real Mode (With Safety Limits)**

```bash
# In .env
USE_MOCK_MCP=false
```

**What you get:**
- ✅ Real bioinformatics analysis
- ✅ Actual tool execution
- ✅ Meaningful results
- ✅ Safety limits protect you

**Best for:**
- Real data analysis
- Building custom tools
- Final projects
- Production workflows

---

## 📈 Monitoring Your Usage

### Sidebar Dashboard

```
🎓 Student Mode
┌─────────────────────────┐
│ Tokens Used             │
│ 12,450 / 50,000        │
│ [████████░░░░░] 25%     │
│                         │
│ Requests                │
│ 8 / 50                  │
│ [████░░░░░░░░] 16%     │
│                         │
│ 💰 Est. cost: $0.11    │
│ ⏱️ Session: 15m 32s     │
└─────────────────────────┘
```

### After Each Query

```
📊 Session: 8/50 requests, 12,450/50,000 tokens
```

### Warnings

At 80% usage:
```
⚠️ Token Usage Warning: 40,000 / 50,000 tokens used (10,000 remaining)
```

---

## 💡 Best Practices

### 1. **Start with Mock Mode**

Week 1-2, use mock mode to:
- Practice writing good prompts
- Learn what each server does
- Experiment freely (it's free!)

### 2. **Switch to Real Mode Gradually**

Week 3+, use real mode for:
- Homework requiring actual analysis
- Custom tool testing
- Real data exploration

### 3. **Monitor Your Usage**

Keep an eye on the sidebar:
- If tokens spike → your query might be too broad
- If requests stack up → you might have a loop

### 4. **Reset When Needed**

Hit a limit? No problem:
- Click "Clear conversation" (sidebar)
- Or refresh the page
- Start fresh!

### 5. **Ask Focused Questions**

Good prompts are:
- ✅ Specific: "Calculate spatial autocorrelation for CD4"
- ✅ Clear: "Run HAllA on proteomics vs metabolomics"

Bad prompts are:
- ❌ Vague: "Analyze everything"
- ❌ Open-ended: "What can you tell me?"

---

## 🆘 Troubleshooting

### "Session Request Limit Reached"

**Cause**: You made 50 queries in this session

**Solution**:
1. Click "Clear conversation" in sidebar
2. Or refresh the page
3. Your usage resets

### "Session Token Limit Reached"

**Cause**: You used 50,000 tokens (~$1.50)

**Solution**:
1. Great job exploring!
2. Clear conversation or refresh
3. Continue learning

### "API Key Missing"

**Cause**: Environment variable not set

**Solution**:
1. Check .env file has the key
2. Restart the app
3. Contact instructor if still broken

### "Mock mode but I want real results"

**Cause**: USE_MOCK_MCP=true in .env

**Solution**:
1. Edit .env: `USE_MOCK_MCP=false`
2. Restart the app
3. Now you get real analysis!

### "Too slow / not responding"

**Possible causes**:
- Real MCP servers starting up (first query ~30s)
- Large dataset processing
- Network issues

**Solutions**:
- Wait 30-60 seconds for first query
- Try smaller test datasets
- Switch to mock mode for practice

---

## 📚 Example Queries

### Spatial Analysis

```
Calculate spatial autocorrelation for gene CD8A in my ovarian cancer dataset
```

```
Identify spatially variable genes in my Visium data
```

```
Deconvolve cell types from my spatial transcriptomics data
```

### Multi-Omics

```
Run HAllA to find associations between my proteomics and metabolomics data
```

```
Combine p-values from my three omics studies using Stouffer's method
```

```
Integrate RNA-seq and protein abundance data
```

### Quantum Analysis

```
Calculate quantum fidelity between T cells and B cells
```

```
Prepare quantum states from my gene expression profiles
```

### File Validation

```
Validate my FASTQ files for quality issues
```

```
Check if my BAM file is properly formatted
```

---

## 🎓 Homework Tips

### Assignment 1-2 (Mock Mode)

Use mock mode for practice:
- Free unlimited usage
- Instant responses
- Perfect for learning prompts

### Assignment 3+ (Real Mode)

Switch to real mode:
- You have 50,000 tokens (~30-40 queries)
- Monitor usage as you go
- Reset if you hit limits

### Managing Costs

Each student session costs ~$0.50-1.50 max:
- Your 50K token limit caps cost
- Most assignments use <20K tokens
- Reset lets you start fresh

**Study group total cost: ~$30-60 for entire 3-month program**

---

## 🔒 What's Different from Instructor's App?

| Feature | Student App | Instructor App |
|---------|------------|----------------|
| **Token Limits** | ✅ 50K per session | ❌ Unlimited |
| **Request Limits** | ✅ 50 per session | ❌ Unlimited |
| **Authentication** | ✅ None (public) | ℹ️ Optional SSO |
| **Usage Tracking** | ✅ Visible sidebar | ℹ️ Optional |
| **Warnings** | ✅ At 80% usage | ❌ None |
| **Cost Protection** | ✅ Auto-stops | ❌ None |
| **MCP Servers** | ✅ Same servers | ✅ Same servers |
| **Providers** | ✅ Gemini only (Flash/Pro) | ✅ Claude + Gemini |

**You share the same MCP infrastructure** - only the safety limits differ!

---

## 💬 Getting Help

### During Study Group

1. **Ask in Slack/Discord** - Other students might have same question
2. **Pair with a buddy** - Two heads better than one
3. **Check DEVELOPMENT.md** - Technical details there

### Common Questions

**Q: Can I increase my limits?**
A: No - limits are for your protection. Reset conversation instead.

**Q: What if I need more tokens for homework?**
A: Most assignments fit easily. If not, split into multiple sessions.

**Q: Why mock mode first?**
A: Free practice! Learn without spending tokens.

**Q: Can I use instructor's app?**
A: Ask your instructor - they might grant access for specific projects.

---

## 🎯 Success Metrics

Track your progress:

**Week 1-2:**
- ✅ Completed 10+ mock queries
- ✅ Understand what each server does
- ✅ Written effective prompts

**Week 3-4:**
- ✅ Switched to real mode
- ✅ Ran actual analyses on test data
- ✅ Stayed within token limits

**Week 5-6:**
- ✅ Built custom MCP tool
- ✅ Integrated into workflow
- ✅ Helped another student

---

## 🚀 Ready to Start?

1. Open the student app URL (or run locally)
2. Check sidebar shows "🎓 Student Mode"
3. Try your first query: `"What tools are available in spatialtools?"`
4. Watch your usage counters
5. Have fun learning! 🧬

**Remember**: Limits are your friend - they let you explore freely without worry!

---

