# Chapter 4: Clinical Data—The Starting Point

> *"Every analysis begins with a patient. How do we integrate EHR data?"*

---

## Why Clinical Data Comes First

Before you analyze genomics, transcriptomics, or imaging, you need to answer one fundamental question:

**Who is this patient?**

Not their name (that's protected health information). But essential clinical context:
- Age and demographics (risk factors vary by age and ancestry)
- Diagnosis and staging (IIIA vs. IV makes treatment completely different)
- Treatment history (platinum-sensitive vs. platinum-resistant)
- Lab markers (CA-125 trends show response or progression)
- Medications (what's working? what failed?)

This information lives in the **Electronic Health Record (EHR)**—typically Epic, Cerner, or similar hospital systems. Without it, your genomic analysis is guesswork. You might identify a PIK3CA mutation, but without knowing the patient is platinum-resistant with rising CA-125, you can't prioritize treatment recommendations effectively.

Sarah's case from Chapter 1 illustrates this perfectly. Her TP53 and PIK3CA mutations are important, but the *clinical context*—8-month platinum-free interval, CA-125 trajectory from 45 to 310 U/mL, BRCA1 germline mutation status—is what transforms those mutations into actionable intelligence.

This chapter teaches you how to build **mcp-epic**: an MCP server that connects to Epic FHIR APIs, retrieves clinical data, and automatically de-identifies it for HIPAA compliance.

---

## The FHIR Standard: Healthcare's Common Language

Healthcare data is notoriously fragmented. Different EHR systems use different formats, terminologies, and APIs. HL7 FHIR (Fast Healthcare Interoperability Resources) is the modern standard that solves this.

### What Is FHIR?

FHIR is a REST API standard for healthcare data exchange. It defines:
- **Resources**: Standardized data structures (Patient, Condition, Observation, Medication)
- **Terminologies**: Standard codes (ICD-10 for diagnoses, LOINC for labs, RxNorm for drugs)
- **Operations**: HTTP verbs for CRUD (Create, Read, Update, Delete)

**FHIR R4** is the current production version (released 2019), supported by Epic, Cerner, Allscripts, and athenahealth.

### Example: Patient Resource

Here's what a FHIR Patient resource looks like:

```json
{
  "resourceType": "Patient",
  "id": "RESEARCH-PAT001",
  "identifier": [{
    "system": "urn:oid:hospital.mrn",
    "value": "MRN12345"
  }],
  "name": [{
    "family": "Anderson",
    "given": ["Sarah"]
  }],
  "gender": "female",
  "birthDate": "1966-03-15",
  "address": [{
    "line": ["123 Main St"],
    "city": "Boston",
    "state": "MA",
    "postalCode": "02115"
  }]
}
```

FHIR specification: https://www.hl7.org/fhir/patient.html

### The Problem: HIPAA Identifiers

Notice what's in that Patient resource:
- **Name**: "Anderson, Sarah"
- **Birth date**: Full date (1966-03-15)
- **Address**: Street, city, zip code
- **Medical record number**: MRN12345

All of these are **HIPAA identifiers** that must be removed before use in research or AI systems.

---

## HIPAA Safe Harbor: The De-identification Standard

The Health Insurance Portability and Accountability Act (HIPAA) defines two methods for de-identifying patient data:
1. **Safe Harbor**: Remove 18 specific identifiers (mechanical, clear rules)
2. **Expert Determination**: Statistical disclosure risk analysis (requires expert, expensive)

We use **Safe Harbor** because it's deterministic and automatable.

### The 18 HIPAA Identifiers

You must remove:

1. Names
2. Geographic subdivisions smaller than state (addresses, zip codes)
3. Dates (except year) — special rule: ages >89 aggregated to ">89"
4. Telephone numbers
5. Fax numbers
6. Email addresses
7. Social Security numbers
8. Medical record numbers
9. Health plan beneficiary numbers
10. Account numbers
11. Certificate/license numbers
12. Vehicle identifiers
13. Device identifiers and serial numbers
14. Web URLs
15. IP addresses
16. Biometric identifiers
17. Full-face photos
18. Any other unique identifying characteristics

HIPAA de-identification guidance: https://www.hhs.gov/hipaa/for-professionals/privacy/special-topics/de-identification/index.html

### Example: De-identified Patient

After Safe Harbor de-identification:

```json
{
  "resourceType": "Patient",
  "id": "deidentified-a4f9c82b1e3d5f",
  "identifier": [{
    "system": "urn:oid:research.hashed-mrn",
    "value": "HASH-d5a7f8c3b1e9"
  }],
  "gender": "female",
  "birthDate": "1966",  // Year only
  "_deidentified": true,
  "_method": "HIPAA Safe Harbor"
}
```

What's removed:
- ❌ Name
- ❌ Full birth date (kept year: 1966)
- ❌ Address
- ✅ Gender (not an identifier under HIPAA)
- ✅ ID (hashed using SHA-256)

---

## Building mcp-epic: Architecture

The mcp-epic server has four layers:

```
┌─────────────────────────────────────┐
│  MCP Tools (FastMCP decorator)     │  ← Natural language interface
├─────────────────────────────────────┤
│  FHIR Client (Epic API wrapper)    │  ← HTTP requests to Epic
├─────────────────────────────────────┤
│  De-identification Layer            │  ← HIPAA Safe Harbor
├─────────────────────────────────────┤
│  Epic FHIR API (OAuth 2.0)          │  ← Hospital EHR system
└─────────────────────────────────────┘
```

Server structure:
```
servers/mcp-epic/
├── src/mcp_epic/
│   ├── server.py              # MCP tools (4 tools)
│   ├── epic_fhir_client.py    # Epic API client
│   ├── deidentify.py          # HIPAA Safe Harbor implementation
│   └── __main__.py            # Entry point
├── pyproject.toml
└── Dockerfile
```

Repository: [`servers/mcp-epic/`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/servers/mcp-epic)

---

## Implementation: The De-identification Layer

Let's build the core de-identification function.

### Step 1: Hash Identifiers

Patient IDs and medical record numbers must be hashed (one-way transformation):

```python
import hashlib

def hash_identifier(value: str) -> str:
    """Hash identifier using SHA-256."""
    hashed = hashlib.sha256(value.encode()).hexdigest()[:16]
    return f"HASH-{hashed}"

# Example
original_id = "MRN12345"
hashed_id = hash_identifier(original_id)
# Result: "HASH-d5a7f8c3b1e9"
```

See implementation: [`servers/mcp-epic/src/mcp_epic/deidentify.py:69-83`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-epic/src/mcp_epic/deidentify.py#L69-L83)

### Step 2: Remove Direct Identifiers

Names, addresses, contact info must be deleted:

```python
def deidentify_patient(patient: dict) -> dict:
    """Apply HIPAA Safe Harbor de-identification."""
    deidentified = patient.copy()

    # Remove direct identifiers
    identifiers_to_remove = [
        "name",          # Names
        "telecom",       # Phone, fax, email
        "address",       # Geographic subdivisions
        "photo",         # Full-face photos
        "contact",       # Contact persons
    ]

    for identifier in identifiers_to_remove:
        if identifier in deidentified:
            del deidentified[identifier]

    return deidentified
```

Full implementation: [`servers/mcp-epic/src/mcp_epic/deidentify.py:19-110`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-epic/src/mcp_epic/deidentify.py#L19-L110)

### Step 3: Date Reduction

Dates must be reduced to year only, with special handling for ages >89:

```python
from datetime import datetime

def reduce_date_to_year(date_string: str) -> str:
    """Reduce date to year, aggregate ages >89."""
    birth_date = datetime.fromisoformat(date_string)
    birth_year = birth_date.year
    current_year = datetime.utcnow().year
    age = current_year - birth_year

    if age > 89:
        # HIPAA: aggregate ages >89
        return f"{current_year - 90}"  # Shows as ">89"
    else:
        return str(birth_year)

# Examples
reduce_date_to_year("1966-03-15")  # "1966" (age 60)
reduce_date_to_year("1920-01-01")  # "1936" (>89, aggregated)
```

See implementation: [`servers/mcp-epic/src/mcp_epic/deidentify.py:84-100`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-epic/src/mcp_epic/deidentify.py#L84-L100)

---

## The Four MCP Tools

The mcp-epic server exposes four tools for clinical data retrieval.

### Tool 1: get_patient_demographics

Retrieves patient age, gender, and hashed identifiers:

```python
from fastmcp import FastMCP

mcp = FastMCP("epic")

@mcp.tool()
async def get_patient_demographics(patient_id: str) -> dict:
    """Retrieve patient demographics from Epic FHIR API.

    All data is automatically de-identified using HIPAA Safe Harbor.

    Args:
        patient_id: Patient identifier (e.g., "RESEARCH-PAT001")

    Returns:
        De-identified demographics including gender, birth year,
        hashed ID, and de-identification metadata.
    """
    client = get_epic_client()
    patient = await client.get_patient(patient_id)

    return {
        "status": "success",
        "data": patient,
        "deidentified": True
    }
```

Full implementation: [`servers/mcp-epic/src/mcp_epic/server.py:40-77`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-epic/src/mcp_epic/server.py#L40-L77)

**Example output**:
```json
{
  "status": "success",
  "data": {
    "id": "deidentified-a4f9c82b1e3d5f",
    "gender": "female",
    "birthDate": "1966",
    "_deidentified": true
  },
  "deidentified": true
}
```

### Tool 2: get_patient_conditions

Retrieves diagnoses with ICD-10 codes and staging:

```python
@mcp.tool()
async def get_patient_conditions(
    patient_id: str,
    category: str = None
) -> dict:
    """Retrieve patient conditions/diagnoses.

    Args:
        patient_id: Patient identifier
        category: Optional filter (e.g., "encounter-diagnosis")

    Returns:
        De-identified list of Condition resources including:
        - ICD-10 diagnosis codes
        - Clinical status (active, resolved)
        - Cancer staging (if applicable)
        - Year of diagnosis
    """
    # Implementation...
```

Full implementation: [`servers/mcp-epic/src/mcp_epic/server.py:79-120`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-epic/src/mcp_epic/server.py#L79-L120)

**Example output**:
```json
{
  "status": "success",
  "data": [{
    "code": "C56.9",  // ICD-10: Ovarian cancer
    "display": "Malignant neoplasm of ovary",
    "clinicalStatus": "active",
    "stage": {
      "summary": "Stage IV",
      "type": "TNM",
      "assessment": "M1 (distant metastasis)"
    },
    "recordedDate": "2023"  // Year only
  }],
  "count": 1
}
```

### Tool 3: get_patient_observations

Retrieves lab results and vital signs:

```python
@mcp.tool()
async def get_patient_observations(
    patient_id: str,
    category: str = None,
    code: str = None,
    limit: int = 100
) -> dict:
    """Retrieve patient observations (labs, vitals).

    Args:
        patient_id: Patient identifier
        category: Optional filter ("laboratory", "vital-signs")
        code: Optional specific observation code (LOINC)
        limit: Maximum number of results

    Returns:
        De-identified observations with values, units,
        reference ranges, and year of observation.
    """
    # Implementation...
```

**Example output (CA-125 tumor marker)**:
```json
{
  "status": "success",
  "data": [{
    "code": "10334-1",  // LOINC: CA-125
    "display": "CA 125 [Units/volume] in Serum",
    "valueQuantity": {
      "value": 310,
      "unit": "U/mL"
    },
    "referenceRange": {
      "low": 0,
      "high": 35
    },
    "interpretation": "High",
    "effectiveDateTime": "2025"  // Year only
  }]
}
```

### Tool 4: get_patient_medications

Retrieves current and historical medications:

```python
@mcp.tool()
async def get_patient_medications(
    patient_id: str,
    status: str = None
) -> dict:
    """Retrieve patient medications.

    Args:
        patient_id: Patient identifier
        status: Optional filter ("active", "completed")

    Returns:
        De-identified medication list with drug names,
        dosages, and treatment year.
    """
    # Implementation...
```

**Example output**:
```json
{
  "status": "success",
  "data": [{
    "medicationCodeableConcept": {
      "coding": [{
        "system": "http://www.nlm.nih.gov/research/umls/rxnorm",
        "code": "198042",
        "display": "Carboplatin"
      }]
    },
    "dosage": "AUC 5, IV every 3 weeks",
    "status": "completed",
    "effectivePeriod": {
      "start": "2023",
      "end": "2023"
    }
  }]
}
```

Full implementation: [`servers/mcp-epic/src/mcp_epic/server.py:177-210`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-epic/src/mcp_epic/server.py#L177-L210)

---

## Epic FHIR Client: Authentication and API Calls

The FHIR client handles OAuth 2.0 authentication and HTTP requests to Epic.

### OAuth 2.0 Client Credentials Flow

Epic uses OAuth 2.0 for API access:

```python
import httpx
import os

class EpicFHIRClient:
    def __init__(self):
        self.base_url = os.getenv("EPIC_FHIR_ENDPOINT")
        self.client_id = os.getenv("EPIC_CLIENT_ID")
        self.client_secret = os.getenv("EPIC_CLIENT_SECRET")
        self.token = None

    async def get_access_token(self) -> str:
        """Obtain OAuth 2.0 access token."""
        async with httpx.AsyncClient() as client:
            response = await client.post(
                f"{self.base_url}/oauth2/token",
                data={
                    "grant_type": "client_credentials",
                    "client_id": self.client_id,
                    "client_secret": self.client_secret
                }
            )
            self.token = response.json()["access_token"]
            return self.token
```

Full implementation: [`servers/mcp-epic/src/mcp_epic/epic_fhir_client.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-epic/src/mcp_epic/epic_fhir_client.py)

### Example: Fetching a Patient

```python
async def get_patient(self, patient_id: str) -> dict:
    """Fetch Patient resource from Epic FHIR API."""
    if not self.token:
        await self.get_access_token()

    url = f"{self.base_url}/Patient/{patient_id}"
    headers = {"Authorization": f"Bearer {self.token}"}

    async with httpx.AsyncClient() as client:
        response = await client.get(url, headers=headers)
        patient = response.json()

        # De-identify before returning
        from .deidentify import deidentify_patient
        return deidentify_patient(patient)
```

Note: De-identification happens *before* data leaves the FHIR client. This ensures no PHI is ever stored in logs or returned to callers.

---

## Configuration: Environment Variables

The server requires four environment variables:

```bash
# Epic FHIR endpoint
EPIC_FHIR_ENDPOINT="https://hospital.epic.com/api/FHIR/R4/"

# OAuth 2.0 credentials (from Epic App Orchard)
EPIC_CLIENT_ID="abc123-your-client-id"
EPIC_CLIENT_SECRET="your-secret-here"

# De-identification toggle (default: enabled)
DEIDENTIFY_ENABLED="true"
```

Get Epic credentials: https://fhir.epic.com/Documentation?docId=epiconfhirrequestprocess

---

## Testing with Synthetic Data: mcp-mockepic

Before connecting to a real Epic system, test with **mcp-mockepic**—a server that returns synthetic FHIR data.

### Why a Mock Server?

1. **No credentials needed**: Works without Epic access
2. **Instant testing**: No OAuth setup, no VPN
3. **Reproducible**: Same synthetic patients every time
4. **Public demos**: Safe for presentations and training

### Using mcp-mockepic

```python
# Same tool interface as mcp-epic
@mcp.tool()
async def query_patient_records(patient_id: str) -> dict:
    """Retrieve synthetic patient from mock database.

    Returns same structure as mcp-epic but with
    pre-generated synthetic data (Synthea-based).
    """
    # Returns synthetic patient PAT001-OVC-2025
```

Deployed endpoint: https://mcp-mockepic-ondu7mwjpa-uc.a.run.app

**Example: Testing PatientOne**

```bash
# Using Claude Desktop with mcp-mockepic
prompt: "Retrieve clinical summary for PAT001-OVC-2025"

# Returns:
{
  "patient_id": "PAT001-OVC-2025",
  "demographics": {
    "age": 58,
    "gender": "female",
    "brca_status": "BRCA1 pathogenic variant carrier"
  },
  "diagnosis": {
    "code": "C56.9",
    "description": "Stage IV High-Grade Serous Ovarian Carcinoma",
    "year": "2023"
  },
  "ca125_trend": [
    {"year": 2023, "value": 1200, "interpretation": "High (at diagnosis)"},
    {"year": 2023, "value": 45, "interpretation": "Near-normal (post-chemo)"},
    {"year": 2025, "value": 310, "interpretation": "Rising (platinum-resistant recurrence)"}
  ]
}
```

Mock server implementation: [`servers/mcp-mockepic/src/mcp_mockepic/server.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-mockepic/src/mcp_mockepic/server.py)

---

## Integration with Other Servers

Clinical data flows to downstream analysis servers:

**mcp-fgbio** (genomics):
- Uses diagnosis (ovarian cancer) to select relevant gene panels
- Filters variants by cancer type (TP53, PIK3CA, PTEN for HGSOC)

**mcp-multiomics**:
- Uses treatment history (platinum-sensitive vs. resistant) to group samples
- Correlates drug response with molecular signatures

**mcp-spatialtools** (spatial transcriptomics):
- Links clinical outcomes (progression-free survival) to spatial patterns
- Maps treatment response to tissue regions

**mcp-tcga** (cohort comparison):
- Uses diagnosis and staging to select matching TCGA cohort
- Compares patient's molecular profile to reference population

Example integration prompt:
```
Retrieve clinical data for PAT001-OVC-2025 using mcp-epic, then:
1. Use diagnosis to select ovarian cancer gene panel (mcp-fgbio)
2. Group PDX samples by platinum sensitivity from medication history (mcp-multiomics)
3. Correlate CA-125 trajectory with spatial immune infiltration (mcp-spatialtools)
```

---

## Deployment: Local Only for HIPAA

**Critical**: mcp-epic uses **STDIO transport**, not SSE over HTTP. This means:
- ✅ Runs on local machine (hospital workstation or VPN)
- ✅ No network exposure of PHI
- ✅ Complies with HIPAA requirement for local processing
- ❌ Cannot deploy to Cloud Run (would expose Epic credentials and patient data)

### Claude Desktop Configuration

```json
{
  "mcpServers": {
    "epic": {
      "command": "python",
      "args": ["-m", "mcp_epic"],
      "env": {
        "EPIC_FHIR_ENDPOINT": "https://hospital.epic.com/api/FHIR/R4/",
        "EPIC_CLIENT_ID": "your-client-id",
        "EPIC_CLIENT_SECRET": "your-secret",
        "DEIDENTIFY_ENABLED": "true"
      }
    }
  }
}
```

Setup guide: [`servers/mcp-epic/CLAUDE_DESKTOP_TESTING.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-epic/CLAUDE_DESKTOP_TESTING.md)

---

## Validation: HIPAA Compliance Checklist

Before production use, validate de-identification:

**✅ Checklist**:
- [ ] All 18 HIPAA identifiers removed (run unit tests)
- [ ] Dates reduced to year only (check birthDate fields)
- [ ] Ages >89 aggregated (verify _ageAggregated field)
- [ ] IDs hashed with SHA-256 (check id and identifier fields)
- [ ] No PHI in logs (check application logs for names/MRNs)
- [ ] BAA in place with Google Cloud (hospital legal)
- [ ] Audit logging enabled (10-year retention)

**Unit test example**:
```python
def test_deidentify_patient_removes_name():
    """Verify names are removed during de-identification."""
    patient = {
        "id": "PAT001",
        "name": [{"family": "Anderson", "given": ["Sarah"]}],
        "birthDate": "1966-03-15",
        "gender": "female"
    }

    deidentified = deidentify_patient(patient)

    assert "name" not in deidentified
    assert deidentified["birthDate"] == "1966"  # Year only
    assert deidentified["id"].startswith("deidentified-")
```

Test suite: [`tests/unit/mcp-epic/test_deidentification.py`](https://github.com/lynnlangit/precision-medicine-mcp/tree/main/tests/unit/mcp-epic)

HIPAA compliance documentation: [`docs/for-hospitals/compliance/hipaa.md`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/docs/for-hospitals/compliance/hipaa.md)

---

## Try It Yourself

Ready to retrieve clinical data from Epic?

### Option 1: Test with mcp-mockepic (No Credentials)

Deploy the mock server locally:

```bash
git clone https://github.com/lynnlangit/precision-medicine-mcp.git
cd precision-medicine-mcp/servers/mcp-mockepic

# Install dependencies
python -m venv venv
source venv/bin/activate
pip install -e .

# Run server
python -m mcp_mockepic
```

Configure Claude Desktop to use it (STDIO transport).

### Option 2: Connect to Epic Sandbox

Epic provides a public sandbox for testing:
- Endpoint: https://fhir.epic.com/interconnect-fhir-oauth/api/FHIR/R4/
- Register for credentials: https://fhir.epic.com/Developer/Apps

Follow Epic's sandbox tutorial, then configure mcp-epic with your credentials.

### Option 3: Interactive Notebook

Explore FHIR resources and de-identification:
[`docs/book/companion-notebooks/chapter-04-clinical-data.ipynb`](../companion-notebooks/chapter-04-clinical-data.ipynb)

This notebook includes:
- FHIR resource parsing examples
- De-identification validation tests
- CA-125 trend visualization
- Integration with genomics prompts

---

## What Comes Next

In Chapter 5, you'll build **mcp-fgbio** for genomic QC and variant calling. You'll learn to:
- Parse VCF files (Variant Call Format)
- Validate genomic data quality (depth, allele frequency)
- Annotate variants with ClinVar and gnomAD
- Filter for clinically relevant mutations

But first, appreciate what you've built. You now have:
- A production FHIR client for Epic EHR systems
- HIPAA-compliant de-identification (Safe Harbor method)
- Four clinical data tools (demographics, conditions, observations, medications)
- Integration with downstream analysis servers
- A mock server for testing without credentials

Clinical data is the foundation. Everything else builds on this.

**Next**: Chapter 5 - Genomic Foundations

---

**Chapter 4 Key Takeaways:**
- Clinical data provides essential patient context for precision medicine
- FHIR R4 is the standard for healthcare data interoperability
- HIPAA Safe Harbor removes 18 identifiers (names, dates, addresses, etc.)
- mcp-epic: 4 tools for demographics, conditions, observations, medications
- De-identification happens automatically before data leaves FHIR client
- Local STDIO deployment (not Cloud Run) for HIPAA compliance
- mcp-mockepic provides synthetic data for testing

**Companion Resources:**
- 📓 [Jupyter Notebook](../companion-notebooks/chapter-04-clinical-data.ipynb) - Hands-on FHIR examples
- 🏥 [HIPAA Compliance Guide](../../for-hospitals/compliance/hipaa.md) - Full compliance checklist
- 📋 [Clinical Architecture](../../architecture/clinical/README.md) - Technical deep dive
- 🔧 [Epic Setup Guide](../../servers/mcp-epic/CLAUDE_DESKTOP_TESTING.md) - Local configuration

**GitHub References:**
- mcp-epic server: [`servers/mcp-epic/src/mcp_epic/server.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-epic/src/mcp_epic/server.py)
- De-identification: [`servers/mcp-epic/src/mcp_epic/deidentify.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-epic/src/mcp_epic/deidentify.py)
- FHIR client: [`servers/mcp-epic/src/mcp_epic/epic_fhir_client.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-epic/src/mcp_epic/epic_fhir_client.py)
- mcp-mockepic: [`servers/mcp-mockepic/src/mcp_mockepic/server.py`](https://github.com/lynnlangit/precision-medicine-mcp/blob/main/servers/mcp-mockepic/src/mcp_mockepic/server.py)
