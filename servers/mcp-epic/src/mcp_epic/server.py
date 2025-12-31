"""MCP server for Epic FHIR integration with HIPAA de-identification.

This server connects to Epic FHIR API (research endpoint) and provides tools for:
- Retrieving patient demographics
- Querying patient conditions
- Accessing patient observations
- Retrieving patient medications

All data is automatically de-identified using HIPAA Safe Harbor method.

Configuration:
    EPIC_FHIR_ENDPOINT: Epic FHIR base URL (e.g., https://hospital.epic.com/api/FHIR/R4/)
    EPIC_CLIENT_ID: OAuth 2.0 client ID
    EPIC_CLIENT_SECRET: OAuth 2.0 client secret
    DEIDENTIFY_ENABLED: Enable/disable de-identification (default: true)
"""

import logging
import os
from typing import Any, Dict, Optional

from fastmcp import FastMCP

from .epic_fhir_client import get_epic_client
from .deidentify import deidentify_patient, deidentify_bundle

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Initialize FastMCP server
mcp = FastMCP("epic")

# Configuration
DEIDENTIFY_ENABLED = os.getenv("DEIDENTIFY_ENABLED", "true").lower() == "true"
EPIC_ENDPOINT = os.getenv("EPIC_FHIR_ENDPOINT", "")



@mcp.tool()
async def get_patient_demographics(patient_id: str) -> Dict[str, Any]:
    """Retrieve patient demographics from Epic FHIR API.

    Fetches patient demographic information including age, gender, and identifiers.
    All data is automatically de-identified using HIPAA Safe Harbor method.

    Args:
        patient_id: Patient identifier (e.g., "RESEARCH-PAT001")

    Returns:
        De-identified patient demographics including:
        - Hashed patient ID
        - Gender
        - Birth year (not full date, ages >89 aggregated)
        - Active status
        - Resource type
        - De-identification metadata
    """
    try:
        client = get_epic_client()
        patient = await client.get_patient(patient_id)

        logger.info(f"Retrieved and de-identified patient from Epic: {patient_id}")
        return {
            "status": "success",
            "data": patient,
            "source": "Epic FHIR API",
            "deidentified": patient.get("_deidentified", False),
        }
    except Exception as e:
        logger.error(f"Error retrieving patient {patient_id} from Epic: {e}")
        return {
            "status": "error",
            "error": "Failed to retrieve patient from Epic FHIR",
            "message": str(e),
        }


@mcp.tool()
async def get_patient_conditions(patient_id: str, category: Optional[str] = None) -> Dict[str, Any]:
    """Retrieve patient conditions/diagnoses from Epic FHIR API.

    Fetches all active and historical conditions for the patient.
    All data is automatically de-identified.

    Args:
        patient_id: Patient identifier (e.g., "RESEARCH-PAT001")
        category: Optional category filter (e.g., "encounter-diagnosis")

    Returns:
        De-identified list of Condition resources including:
        - Diagnosis codes (ICD-10) and descriptions
        - Clinical status (active, resolved, etc.)
        - Verification status
        - Severity
        - Stage information (for cancer diagnoses)
        - Year of diagnosis (dates reduced to year)
    """
    try:
        client = get_epic_client()
        conditions = await client.get_conditions(patient_id, category=category)

        count = len(conditions)
        logger.info(f"Retrieved {count} condition(s) from Epic for patient: {patient_id}")

        return {
            "status": "success",
            "data": conditions,
            "count": count,
            "source": "Epic FHIR API",
            "deidentified": all(c.get("_deidentified", False) for c in conditions),
        }
    except Exception as e:
        logger.error(f"Error retrieving conditions from Epic for {patient_id}: {e}")
        return {
            "status": "error",
            "error": "Failed to retrieve conditions from Epic FHIR",
            "message": str(e),
        }


@mcp.tool()
async def get_patient_observations(
    patient_id: str,
    category: Optional[str] = None,
    code: Optional[str] = None,
    limit: int = 100,
) -> Dict[str, Any]:
    """Retrieve patient observations (lab results, vitals, etc.) from Epic FHIR API.

    Fetches laboratory results, vital signs, and other clinical observations.
    All data is automatically de-identified.

    Args:
        patient_id: Patient identifier (e.g., "RESEARCH-PAT001")
        category: Optional category filter (e.g., "laboratory", "vital-signs")
        code: Optional specific observation code (LOINC, etc.)
        limit: Maximum number of results (default: 100)

    Returns:
        De-identified list of Observation resources including:
        - Test/observation codes (LOINC, SNOMED)
        - Values and units (e.g., "120 mmHg", "5.5 mg/dL")
        - Reference ranges
        - Interpretation (e.g., high, low, normal)
        - Year of observation (dates reduced to year)
    """
    try:
        client = get_epic_client()
        observations = await client.get_observations(
            patient_id, category=category, code=code, limit=limit
        )

        count = len(observations)
        logger.info(f"Retrieved {count} observation(s) from Epic for patient: {patient_id}")

        return {
            "status": "success",
            "data": observations,
            "count": count,
            "source": "Epic FHIR API",
            "deidentified": all(o.get("_deidentified", False) for o in observations),
        }
    except Exception as e:
        logger.error(f"Error retrieving observations from Epic for {patient_id}: {e}")
        return {
            "status": "error",
            "error": "Failed to retrieve observations from Epic FHIR",
            "message": str(e),
        }


@mcp.tool()
async def get_patient_medications(
    patient_id: str,
    status: Optional[str] = None,
    limit: int = 100,
) -> Dict[str, Any]:
    """Retrieve patient medications from Epic FHIR API.

    Fetches current and historical medication information.
    All data is automatically de-identified.

    Args:
        patient_id: Patient identifier (e.g., "RESEARCH-PAT001")
        status: Optional status filter (e.g., "active", "completed", "stopped")
        limit: Maximum number of results (default: 100)

    Returns:
        De-identified list of MedicationStatement resources including:
        - Medication names (generic and brand)
        - Status (active, completed, stopped)
        - Dosage and route (e.g., "50mg oral daily")
        - Treatment dates (year only)
        - Reason for medication
        - Prescriber information (de-identified)
    """
    try:
        client = get_epic_client()
        medications = await client.get_medications(
            patient_id, status=status, limit=limit
        )

        count = len(medications)
        logger.info(f"Retrieved {count} medication(s) from Epic for patient: {patient_id}")

        return {
            "status": "success",
            "data": medications,
            "count": count,
            "source": "Epic FHIR API",
            "deidentified": all(m.get("_deidentified", False) for m in medications),
        }
    except Exception as e:
        logger.error(f"Error retrieving medications from Epic for {patient_id}: {e}")
        return {
            "status": "error",
            "error": "Failed to retrieve medications from Epic FHIR",
            "message": str(e),
        }


def main():
    """Run the mcp-epic server."""
    logger.info("Starting mcp-epic server with Epic FHIR integration")
    logger.info(f"Epic FHIR Endpoint: {EPIC_ENDPOINT or 'NOT CONFIGURED'}")
    logger.info(f"De-identification: {'ENABLED' if DEIDENTIFY_ENABLED else 'DISABLED'}")

    # Verify Epic client is properly configured
    client = get_epic_client()
    if client.is_configured():
        logger.info("✅ Epic FHIR client configured successfully")
    else:
        logger.warning("⚠️  Epic FHIR client NOT properly configured")
        logger.warning("Set EPIC_FHIR_ENDPOINT, EPIC_CLIENT_ID, EPIC_CLIENT_SECRET")

    mcp.run()


if __name__ == "__main__":
    main()
