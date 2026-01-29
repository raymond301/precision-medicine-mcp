# Student MCP Chat App - For Study Group

This is a **safety-limited version** of the MCP Chat app designed for the bioinformatics study group.

## 🎓 Key Differences from Main App

| Feature | This App (Students) | Main App (Instructor) |
|---------|-------------------|---------------------|
| **Purpose** | Learning with safety limits | Unrestricted testing |
| **Token Limit** | 50,000 per session | Unlimited |
| **Request Limit** | 50 per session | Unlimited |
| **Cost Cap** | ~$1.50 per session | No limit |
| **Usage Display** | Always visible | Optional |
| **API Keys** | Shared student keys | Instructor's keys |
| **Service Name** | `streamlit-mcp-chat-students` | `streamlit-mcp-chat` |

## 🔗 Shared Infrastructure

Both apps use:
- ✅ **Same MCP servers** (no duplication)
- ✅ **Same provider code** (copied from main app)
- ✅ **Same utilities** (copied from main app)

Only `app.py` differs (~100 lines of safety code).

**Note**: `providers/` and `utils/` are copied from `../streamlit-app/` to ensure Docker builds work correctly (symlinks don't work in Docker build context).

## 📁 Directory Structure

```
streamlit-app-students/
├── app.py                    # Student app with guardrails
├── providers/                # Copied from ../streamlit-app/providers/
├── utils/                    # Copied from ../streamlit-app/utils/
├── .env.example              # Student-specific config
├── deploy.sh                 # Deploys to different service
├── STUDENT_GUIDE.md          # Student documentation
└── README.md                 # This file
```

**Maintenance Note**: When updating providers or utils in the main app, copy changes to student app:
```bash
cp -r ../streamlit-app/providers .
cp -r ../streamlit-app/utils .
```

## 🚀 For Students

See **[STUDENT_GUIDE.md](STUDENT_GUIDE.md)** for complete instructions.

## 🚀 For Instructors

### Deploy Student App

```bash
export ANTHROPIC_API_KEY="sk-ant-student-..."
export GEMINI_API_KEY="..."
cd ui/streamlit-app-students
./deploy.sh
```

### Adjust Safety Limits

Edit `app.py`:
```python
MAX_TOKENS_PER_SESSION = 50000  # Change as needed
MAX_REQUESTS_PER_SESSION = 50
```

Then redeploy: `./deploy.sh`

## 📊 Expected Costs

- **Per student per session**: ~$0.50-1.50 (capped)
- **Study group (10 students, 6 weeks)**: ~$50-100 total

## ⚙️ Safety Features

- ✅ Token limits per session
- ✅ Request limits per session
- ✅ Real-time usage tracking
- ✅ Automatic warnings at 80%
- ✅ Easy reset (clear conversation)

See README for full details.
