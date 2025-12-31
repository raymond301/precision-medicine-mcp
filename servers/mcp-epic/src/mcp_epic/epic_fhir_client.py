"""
Epic FHIR Client for Research Hospital Deployment

Handles OAuth 2.0 authentication and API calls to Epic FHIR endpoints.
Implements automatic de-identification using HIPAA Safe Harbor method.
"""

import os
import httpx
import asyncio
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Any
import logging

from .deidentify import (
    deidentify_patient,
    deidentify_observation,
    deidentify_condition,
    deidentify_medication_statement
)

logger = logging.getLogger(__name__)


class EpicFHIRClient:
    """
    Client for Epic FHIR API with OAuth 2.0 authentication.

    Supports:
    - Client credentials OAuth 2.0 flow
    - Automatic token refresh
    - Rate limiting and retry logic
    - Automatic de-identification of all patient data
    - Comprehensive error handling and logging
    """

    def __init__(
        self,
        endpoint: Optional[str] = None,
        client_id: Optional[str] = None,
        client_secret: Optional[str] = None,
        token_url: Optional[str] = None,
        deidentify_enabled: bool = True
    ):
        """
        Initialize Epic FHIR client.

        Args:
            endpoint: Epic FHIR base URL (e.g., https://hospital.epic.com/api/FHIR/R4/)
            client_id: OAuth 2.0 client ID
            client_secret: OAuth 2.0 client secret
            token_url: OAuth 2.0 token endpoint (defaults to endpoint + /oauth2/token)
            deidentify_enabled: Whether to de-identify patient data (default: True)
        """
        # Get configuration from environment if not provided
        self.epic_endpoint = endpoint or os.getenv("EPIC_FHIR_ENDPOINT", "")
        self.client_id = client_id or os.getenv("EPIC_CLIENT_ID", "")
        self.client_secret = client_secret or os.getenv("EPIC_CLIENT_SECRET", "")

        # Remove trailing slash from endpoint
        if self.epic_endpoint.endswith("/"):
            self.epic_endpoint = self.epic_endpoint[:-1]

        # Default token URL if not provided
        if token_url:
            self.token_url = token_url
        elif self.epic_endpoint:
            # Epic typically uses /oauth2/token
            self.token_url = f"{self.epic_endpoint}/oauth2/token"
        else:
            self.token_url = ""

        # De-identification setting
        self.deidentify_enabled = deidentify_enabled

        # Token management
        self.access_token = None
        self.token_expiry = None

        # Rate limiting
        self.request_count = 0
        self.rate_limit_window_start = datetime.utcnow()
        self.max_requests_per_hour = 1000  # Adjust based on Epic's limits

        # Validate configuration
        self._validate_config()

        logger.info(
            f"Epic FHIR client initialized: {self.epic_endpoint}",
            extra={
                "endpoint": self.epic_endpoint,
                "deidentify_enabled": self.deidentify_enabled
            }
        )

    def _validate_config(self):
        """Validate that required configuration is present."""
        missing = []

        if not self.epic_endpoint:
            missing.append("EPIC_FHIR_ENDPOINT")
        if not self.client_id:
            missing.append("EPIC_CLIENT_ID")
        if not self.client_secret:
            missing.append("EPIC_CLIENT_SECRET")

        if missing:
            raise ValueError(
                f"Missing required Epic FHIR configuration: {', '.join(missing)}"
            )

    async def get_access_token(self) -> str:
        """
        Get Epic FHIR OAuth 2.0 access token.

        Uses client credentials flow. Caches token and refreshes when expired.

        Returns:
            str: Valid access token

        Raises:
            httpx.HTTPError: If token request fails
        """
        # Return cached token if still valid (with 60 second buffer)
        if self.access_token and self.token_expiry:
            if self.token_expiry > datetime.utcnow() + timedelta(seconds=60):
                return self.access_token

        logger.info("Requesting new Epic FHIR access token")

        # Request new access token using client credentials flow
        async with httpx.AsyncClient() as client:
            try:
                response = await client.post(
                    self.token_url,
                    data={
                        "grant_type": "client_credentials",
                        "client_id": self.client_id,
                        "client_secret": self.client_secret,
                        "scope": "patient/*.read"  # Read-only access to patient data
                    },
                    timeout=30.0
                )

                response.raise_for_status()
                token_data = response.json()

                # Extract token and expiry
                self.access_token = token_data.get("access_token")
                expires_in = token_data.get("expires_in", 3600)  # Default 1 hour

                if not self.access_token:
                    raise ValueError("No access_token in Epic OAuth response")

                # Set expiry time
                self.token_expiry = datetime.utcnow() + timedelta(seconds=expires_in)

                logger.info(
                    f"Epic FHIR access token obtained, expires in {expires_in}s",
                    extra={"expires_in": expires_in}
                )

                return self.access_token

            except httpx.HTTPStatusError as e:
                logger.error(
                    f"Epic OAuth token request failed: {e.response.status_code} {e.response.text}",
                    extra={
                        "status_code": e.response.status_code,
                        "response": e.response.text[:500]
                    }
                )
                raise

            except Exception as e:
                logger.error(f"Epic OAuth token request error: {e}")
                raise

    async def _make_request(
        self,
        method: str,
        path: str,
        params: Optional[Dict] = None
    ) -> Dict:
        """
        Make authenticated request to Epic FHIR API.

        Args:
            method: HTTP method (GET, POST, etc.)
            path: FHIR resource path (e.g., /Patient/123)
            params: Query parameters

        Returns:
            dict: FHIR response data

        Raises:
            httpx.HTTPError: If request fails
        """
        # Get valid access token
        token = await self.get_access_token()

        # Build full URL
        if path.startswith("/"):
            url = f"{self.epic_endpoint}{path}"
        else:
            url = f"{self.epic_endpoint}/{path}"

        # Check rate limiting
        await self._check_rate_limit()

        # Make request
        async with httpx.AsyncClient() as client:
            try:
                response = await client.request(
                    method,
                    url,
                    headers={
                        "Authorization": f"Bearer {token}",
                        "Accept": "application/fhir+json",
                        "Content-Type": "application/fhir+json"
                    },
                    params=params,
                    timeout=30.0
                )

                response.raise_for_status()

                # Increment request counter
                self.request_count += 1

                # Parse FHIR response
                data = response.json()

                logger.info(
                    f"Epic FHIR request successful: {method} {path}",
                    extra={
                        "method": method,
                        "path": path,
                        "status_code": response.status_code
                    }
                )

                return data

            except httpx.HTTPStatusError as e:
                logger.error(
                    f"Epic FHIR request failed: {e.response.status_code} {e.response.text}",
                    extra={
                        "method": method,
                        "path": path,
                        "status_code": e.response.status_code,
                        "response": e.response.text[:500]
                    }
                )
                raise

            except Exception as e:
                logger.error(
                    f"Epic FHIR request error: {e}",
                    extra={"method": method, "path": path}
                )
                raise

    async def _check_rate_limit(self):
        """Check and enforce rate limiting."""
        # Reset counter every hour
        if datetime.utcnow() - self.rate_limit_window_start > timedelta(hours=1):
            self.request_count = 0
            self.rate_limit_window_start = datetime.utcnow()

        # Check if rate limit exceeded
        if self.request_count >= self.max_requests_per_hour:
            wait_time = 3600 - (datetime.utcnow() - self.rate_limit_window_start).total_seconds()
            logger.warning(
                f"Rate limit reached, waiting {wait_time:.0f}s",
                extra={"request_count": self.request_count}
            )
            await asyncio.sleep(wait_time)
            # Reset after waiting
            self.request_count = 0
            self.rate_limit_window_start = datetime.utcnow()

    # ========================================================================
    # FHIR Resource Methods
    # ========================================================================

    async def get_patient(self, patient_id: str) -> Dict:
        """
        Get patient demographics.

        Args:
            patient_id: Patient ID (e.g., "RESEARCH-PAT001")

        Returns:
            dict: De-identified patient resource

        Raises:
            httpx.HTTPError: If request fails
        """
        # Request patient from Epic
        patient = await self._make_request("GET", f"Patient/{patient_id}")

        # De-identify before returning
        if self.deidentify_enabled:
            patient = deidentify_patient(patient)

        return patient

    async def get_observations(
        self,
        patient_id: str,
        category: Optional[str] = None,
        code: Optional[str] = None,
        limit: int = 100
    ) -> List[Dict]:
        """
        Get observations (labs, vitals, etc.) for a patient.

        Args:
            patient_id: Patient ID
            category: Observation category (e.g., "laboratory", "vital-signs")
            code: Specific observation code (LOINC, etc.)
            limit: Maximum number of results

        Returns:
            list: De-identified observation resources
        """
        params = {
            "patient": patient_id,
            "_count": limit
        }

        if category:
            params["category"] = category
        if code:
            params["code"] = code

        # Request observations
        bundle = await self._make_request("GET", "Observation", params=params)

        # Extract observations from bundle
        observations = []
        if "entry" in bundle:
            for entry in bundle["entry"]:
                if "resource" in entry:
                    obs = entry["resource"]

                    # De-identify
                    if self.deidentify_enabled:
                        obs = deidentify_observation(obs)

                    observations.append(obs)

        return observations

    async def get_conditions(
        self,
        patient_id: str,
        category: Optional[str] = None,
        limit: int = 100
    ) -> List[Dict]:
        """
        Get conditions (diagnoses) for a patient.

        Args:
            patient_id: Patient ID
            category: Condition category (e.g., "encounter-diagnosis")
            limit: Maximum number of results

        Returns:
            list: De-identified condition resources
        """
        params = {
            "patient": patient_id,
            "_count": limit
        }

        if category:
            params["category"] = category

        # Request conditions
        bundle = await self._make_request("GET", "Condition", params=params)

        # Extract conditions from bundle
        conditions = []
        if "entry" in bundle:
            for entry in bundle["entry"]:
                if "resource" in entry:
                    cond = entry["resource"]

                    # De-identify
                    if self.deidentify_enabled:
                        cond = deidentify_condition(cond)

                    conditions.append(cond)

        return conditions

    async def get_medications(
        self,
        patient_id: str,
        status: Optional[str] = None,
        limit: int = 100
    ) -> List[Dict]:
        """
        Get medication statements for a patient.

        Args:
            patient_id: Patient ID
            status: Medication status (e.g., "active", "completed")
            limit: Maximum number of results

        Returns:
            list: De-identified medication statement resources
        """
        params = {
            "patient": patient_id,
            "_count": limit
        }

        if status:
            params["status"] = status

        # Request medication statements
        bundle = await self._make_request("GET", "MedicationStatement", params=params)

        # Extract medications from bundle
        medications = []
        if "entry" in bundle:
            for entry in bundle["entry"]:
                if "resource" in entry:
                    med = entry["resource"]

                    # De-identify
                    if self.deidentify_enabled:
                        med = deidentify_medication_statement(med)

                    medications.append(med)

        return medications

    async def search_patients(
        self,
        given: Optional[str] = None,
        family: Optional[str] = None,
        identifier: Optional[str] = None,
        limit: int = 20
    ) -> List[Dict]:
        """
        Search for patients.

        IMPORTANT: Name-based search should NOT be used in production.
        Use identifier-based search with research patient IDs instead.

        Args:
            given: Given name
            family: Family name
            identifier: Patient identifier (RECOMMENDED for research)
            limit: Maximum number of results

        Returns:
            list: De-identified patient resources
        """
        params = {"_count": limit}

        if identifier:
            params["identifier"] = identifier
        if given:
            params["given"] = given
        if family:
            params["family"] = family

        # Request patients
        bundle = await self._make_request("GET", "Patient", params=params)

        # Extract patients from bundle
        patients = []
        if "entry" in bundle:
            for entry in bundle["entry"]:
                if "resource" in entry:
                    patient = entry["resource"]

                    # De-identify
                    if self.deidentify_enabled:
                        patient = deidentify_patient(patient)

                    patients.append(patient)

        return patients

    async def get_capability_statement(self) -> Dict:
        """
        Get Epic FHIR capability statement.

        Useful for verifying connection and discovering supported resources.

        Returns:
            dict: CapabilityStatement resource
        """
        return await self._make_request("GET", "metadata")

    # ========================================================================
    # Helper Methods
    # ========================================================================

    def is_configured(self) -> bool:
        """Check if client is properly configured."""
        return bool(self.epic_endpoint and self.client_id and self.client_secret)

    async def test_connection(self) -> bool:
        """
        Test Epic FHIR connection.

        Returns:
            bool: True if connection successful
        """
        try:
            await self.get_capability_statement()
            logger.info("Epic FHIR connection test successful")
            return True
        except Exception as e:
            logger.error(f"Epic FHIR connection test failed: {e}")
            return False


# Global client instance
_epic_client = None


def get_epic_client() -> EpicFHIRClient:
    """
    Get or create global Epic FHIR client instance.

    Returns:
        EpicFHIRClient: Configured Epic FHIR client
    """
    global _epic_client

    if _epic_client is None:
        _epic_client = EpicFHIRClient()

    return _epic_client
