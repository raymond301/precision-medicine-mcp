# Tests

This directory contains all **test code** for the Precision Medicine MCP system.

## Contents

- **Unit Tests** (`/tests/unit/`) - 42 pytest unit test files
- **Integration Tests** (`/tests/integration/`) - End-to-end integration tests
- **Verification** (`/tests/verification/`) - Server verification scripts
- **Test Fixtures** (`/tests/unit/mcp-multiomics/fixtures/`) - Test data (CSV files)
- **Test Scripts** (`/tests/manual_testing/Solution-Testing/`) - Setup and verification shell scripts

## Test Documentation

📚 **Test documentation has moved to `/docs/test-docs/`**

For testing guides, test prompts, and testing strategies, see:
- [Test Documentation Index](../docs/reference/test-docs/README.md)
- [Test Coverage & Guidelines](../docs/reference/test-docs/test-coverage.md)
- [Manual Testing Guides](../docs/reference/test-docs/manual-testing/)
- [PatientOne Testing Scenario](../docs/reference/test-docs/patient-one-scenario/)

## Running Tests

### Unit Tests
```bash
# Run all unit tests
pytest tests/unit/ -v

# Run specific server tests
pytest tests/unit/mcp-spatialtools/ -v
pytest tests/unit/mcp-multiomics/ -v
```

### Integration Tests
```bash
# Run integration tests
pytest tests/integration/ -v

# Run CitL end-to-end tests
pytest tests/integration/test_citl_end_to_end.py -v
```

### Manual Testing
See [Manual Testing Documentation](../docs/reference/test-docs/manual-testing/) for copy-paste prompts and testing procedures.

---

**Test Code:** This directory
**Test Documentation:** `/docs/test-docs/`
