# Test Documentation

This directory contains all documentation related to testing the Precision Medicine MCP system.

**📍 Note:** Test **code** (`.py` files, test fixtures, scripts) remains in `/tests`. This directory contains test **documentation** only.

## Contents

### [Test Coverage](./test-coverage.md)
Overview of test coverage, test structure, and testing guidelines for the project.

### [GCP Integration Testing](./gcp-integration.md)
Testing deployed servers via Claude API.

### Manual Testing
- [Quick Test Prompts](./quick-test-prompts.md) - 10 copy-paste prompts for Claude Desktop
- [Claude Desktop Setup](./claude-desktop-setup.md) - File access configuration

### [PatientOne Scenario](./patient-one/)
Complete testing scenario using synthetic ovarian cancer patient data.

- [Overview](./patient-one/README.md) - PatientOne testing scenario
- [Architecture Overview](./patient-one/overview.md) - Architecture overview
- [Quick Reference](./patient-one/quick-reference.md) - Quick test reference
- [CITL Quick Test](./patient-one/citl-quick-test.md) - Clinician-in-the-Loop workflow test
- [Data Modes Guide](./patient-one/data-modes-guide.md) - Data mode configuration
- [Test All Servers](./patient-one/test-all-servers.md) - Complete server testing
- [End-to-End Test](./patient-one/end-to-end-test.md) - Full workflow test

#### Test Prompts
Ready-to-use test prompts for the complete PatientOne workflow:

1. [Test 1: Clinical Genomic](./patient-one/test-prompts/test-1-clinical-genomic.md)
2. [Test 2: Multiomics Enhanced](./patient-one/test-prompts/test-2-multiomics-enhanced.md)
3. [Test 3: Spatial](./patient-one/test-prompts/test-3-spatial.md)
4. [Test 4: Imaging](./patient-one/test-prompts/test-4-imaging.md)
5. [Test 5: Integration](./patient-one/test-prompts/test-5-integration.md)
6. [Test 6: CitL Review](./patient-one/test-prompts/test-6-citl-review.md)


---

## Quick Start

**For automated testing:** See `/tests` directory for pytest unit and integration tests.

**For manual testing:**
1. Start with [Quick Test Prompts](./quick-test-prompts.md) for rapid verification
2. Run [PatientOne scenario](./patient-one/README.md) for comprehensive end-to-end testing
3. Follow [CitL Quick Test](./patient-one/citl-quick-test.md) to validate clinical review workflow

---

## Manual Testing Setup

This section covers scripts and documentation for manually testing the Precision Medicine MCP servers.

### Scripts (Executable)

| File | Purpose | Usage |
|------|---------|-------|
| `install_dependencies.sh` | Install all dependencies for all MCP servers | `./install_dependencies.sh` |
| `verify_servers.sh` | Verify all servers can be imported | `./verify_servers.sh` |
| `setup_and_test_servers.sh` | Combined setup and verification | `./setup_and_test_servers.sh --install` |
| `test_all_servers.py` | Python-based server verification | `python3 test_all_servers.py` |

### Install All Server Dependencies

```bash
./install_dependencies.sh
```

This will:
- Create Python 3.11 virtual environments for each server
- Install FastMCP and all dependencies
- Set up all servers in development mode

**Time:** ~5-10 minutes

### Verify All Servers

```bash
./verify_servers.sh
```

Expected output:
```
Servers working: 15/15
Total tools: 80
🎉 All MCP servers are operational!
```

**Time:** ~10-30 seconds

### Claude Code vs Claude Desktop

**Scripts run in Claude Code** (VSCode extension):
- ✅ Can install dependencies
- ✅ Can verify server code
- ❌ Cannot orchestrate MCP protocol

**To test MCP workflows:**
- ✅ Use Claude Desktop (standalone app)
- ✅ Configure with [`docs/getting-started/desktop-configs/`](../../getting-started/desktop-configs/)

### Python Version Requirement

All servers require **Python 3.11+**. The install script automatically uses `python3.11`.

### Server Test Status

| Server | Tools | Status |
|--------|-------|--------|
| mcp-fgbio | 4 | ✅ |
| mcp-multiomics | 10 | ✅ |
| mcp-spatialtools | 14 | ✅ |
| mcp-epic | 4 | ✅ |
| mcp-mockepic | 3 | ✅ |
| mcp-perturbation | 8 | ✅ |
| mcp-quantum-celltype-fidelity | 6 | ✅ |
| mcp-openimagedata | 5 | ✅ |
| mcp-deepcell | 3 | ✅ |
| mcp-cell-classify | 3 | ✅ |
| mcp-patient-report | 5 | ✅ |
| mcp-genomic-results | 4 | ✅ |
| mcp-tcga | 5 | ✅ |
| mcp-huggingface | 3 | ✅ |
| mcp-seqera | 3 | ✅ |
| **TOTAL** | **80** | **✅** |

---

**Last Updated:** 2026-01-13
