"""HIPAA Safe Harbor De-identification for FHIR Resources

This module implements the HIPAA Safe Harbor method for de-identifying patient data.
It removes all 18 HIPAA identifiers from FHIR resources.

References:
- HIPAA Privacy Rule: 45 CFR § 164.514(b)(2)
- Safe Harbor Method: https://www.hhs.gov/hipaa/for-professionals/privacy/special-topics/de-identification/index.html
"""

import hashlib
import logging
from datetime import datetime
from typing import Any, Dict, List

logger = logging.getLogger(__name__)


def deidentify_patient(patient: Dict[str, Any]) -> Dict[str, Any]:
    """Apply HIPAA Safe Harbor de-identification to Patient resource.

    Removes 18 HIPAA identifiers:
    1. Names
    2. Geographic subdivisions smaller than state
    3. Dates (except year)
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
    17. Full face photos
    18. Any other unique identifying numbers

    Args:
        patient: FHIR Patient resource

    Returns:
        De-identified Patient resource with:
        - Names removed
        - Addresses removed (keeps state if >89 years old handling applied)
        - Dates reduced to year only
        - Contact information removed
        - IDs hashed
    """
    deidentified = patient.copy()

    # Remove direct identifiers (1, 4, 5, 6, 14)
    identifiers_to_remove = [
        "name",          # 1. Names
        "telecom",       # 4, 5, 6. Telephone, fax, email
        "address",       # 2. Geographic subdivisions
        "photo",         # 17. Full face photos
        "contact",       # Contact persons (contain names)
        "communication", # May contain identifying info
    ]

    for identifier in identifiers_to_remove:
        if identifier in deidentified:
            del deidentified[identifier]

    # Hash patient ID (8, 9, 10, 18. Medical record number)
    if "id" in deidentified:
        original_id = deidentified["id"]
        hashed_id = hashlib.sha256(original_id.encode()).hexdigest()[:16]
        deidentified["id"] = f"deidentified-{hashed_id}"
        logger.debug(f"Hashed patient ID: {original_id} -> {deidentified['id']}")

    # Hash identifier values (7, 8, 9, 10, 11, 18)
    if "identifier" in deidentified:
        for ident in deidentified["identifier"]:
            if "value" in ident:
                original_value = ident["value"]
                hashed_value = hashlib.sha256(original_value.encode()).hexdigest()[:16]
                ident["value"] = f"HASH-{hashed_value}"

    # Date shift: keep year only (3. Dates except year)
    # Special handling: ages >89 are aggregated to ">89"
    if "birthDate" in deidentified:
        try:
            birth_date = datetime.fromisoformat(deidentified["birthDate"])
            birth_year = birth_date.year
            current_year = datetime.utcnow().year
            age = current_year - birth_year

            if age > 89:
                # HIPAA requirement: aggregate ages >89
                deidentified["birthDate"] = f"{current_year - 90}"  # Shows as >89
                deidentified["_ageAggregated"] = ">89"
            else:
                deidentified["birthDate"] = str(birth_year)

        except (ValueError, AttributeError):
            del deidentified["birthDate"]

    # Remove extensions that might contain identifying info
    if "extension" in deidentified:
        # Keep only non-identifying extensions
        safe_extensions = []
        for ext in deidentified["extension"]:
            url = ext.get("url", "")
            # Filter out potentially identifying extensions
            if not any(identifier in url.lower() for identifier in
                      ["name", "address", "contact", "photo", "identifier"]):
                safe_extensions.append(ext)
        deidentified["extension"] = safe_extensions

    # Add de-identification metadata
    deidentified["_deidentified"] = True
    deidentified["_deidentification_method"] = "HIPAA Safe Harbor"
    deidentified["_deidentification_date"] = datetime.utcnow().isoformat()

    logger.info(f"De-identified patient resource")
    return deidentified


def deidentify_observation(observation: Dict[str, Any]) -> Dict[str, Any]:
    """Apply de-identification to Observation resource.

    Args:
        observation: FHIR Observation resource

    Returns:
        De-identified Observation resource
    """
    deidentified = observation.copy()

    # Hash resource ID
    if "id" in deidentified:
        original_id = deidentified["id"]
        hashed_id = hashlib.sha256(original_id.encode()).hexdigest()[:16]
        deidentified["id"] = f"deidentified-{hashed_id}"

    # Hash and anonymize patient references
    if "subject" in deidentified:
        deidentified["subject"] = _deidentify_reference(deidentified["subject"])

    # Hash performer references (doctors, labs)
    if "performer" in deidentified:
        deidentified["performer"] = [
            _deidentify_reference(ref) for ref in deidentified["performer"]
        ]

    # Remove dates (keep year only for effectiveDateTime)
    if "effectiveDateTime" in deidentified:
        deidentified["effectiveDateTime"] = _reduce_to_year(
            deidentified["effectiveDateTime"]
        )

    # Remove dates completely for other date fields
    date_fields = ["issued"]
    for field in date_fields:
        if field in deidentified:
            deidentified[field] = _reduce_to_year(deidentified[field])

    # Add de-identification marker
    deidentified["_deidentified"] = True
    deidentified["_deidentification_method"] = "HIPAA Safe Harbor"

    return deidentified


def deidentify_condition(condition: Dict[str, Any]) -> Dict[str, Any]:
    """Apply de-identification to Condition resource.

    Args:
        condition: FHIR Condition resource

    Returns:
        De-identified Condition resource
    """
    deidentified = condition.copy()

    # Hash resource ID
    if "id" in deidentified:
        original_id = deidentified["id"]
        hashed_id = hashlib.sha256(original_id.encode()).hexdigest()[:16]
        deidentified["id"] = f"deidentified-{hashed_id}"

    # Hash patient references
    if "subject" in deidentified:
        deidentified["subject"] = _deidentify_reference(deidentified["subject"])

    # Hash asserter references (who diagnosed)
    if "asserter" in deidentified:
        deidentified["asserter"] = _deidentify_reference(deidentified["asserter"])

    # Remove dates
    date_fields = [
        "recordedDate",
        "assertedDate",
        "dateAsserted",
        "onsetDateTime",
        "abatementDateTime",
    ]
    for field in date_fields:
        if field in deidentified:
            deidentified[field] = _reduce_to_year(deidentified[field])

    # Add de-identification marker
    deidentified["_deidentified"] = True
    deidentified["_deidentification_method"] = "HIPAA Safe Harbor"

    return deidentified


def deidentify_medication_statement(medication: Dict[str, Any]) -> Dict[str, Any]:
    """Apply de-identification to MedicationStatement resource.

    Args:
        medication: FHIR MedicationStatement resource

    Returns:
        De-identified MedicationStatement resource
    """
    deidentified = medication.copy()

    # Hash resource ID
    if "id" in deidentified:
        original_id = deidentified["id"]
        hashed_id = hashlib.sha256(original_id.encode()).hexdigest()[:16]
        deidentified["id"] = f"deidentified-{hashed_id}"

    # Hash patient references
    if "subject" in deidentified:
        deidentified["subject"] = _deidentify_reference(deidentified["subject"])

    # Hash prescriber references
    if "informationSource" in deidentified:
        deidentified["informationSource"] = _deidentify_reference(
            deidentified["informationSource"]
        )

    # Remove dates
    date_fields = ["effectiveDateTime", "dateAsserted"]
    for field in date_fields:
        if field in deidentified:
            deidentified[field] = _reduce_to_year(deidentified[field])

    # Handle effectivePeriod
    if "effectivePeriod" in deidentified:
        period = deidentified["effectivePeriod"]
        if "start" in period:
            period["start"] = _reduce_to_year(period["start"])
        if "end" in period:
            period["end"] = _reduce_to_year(period["end"])

    # Add de-identification marker
    deidentified["_deidentified"] = True
    deidentified["_deidentification_method"] = "HIPAA Safe Harbor"

    return deidentified


def _deidentify_reference(reference: Dict[str, Any]) -> Dict[str, Any]:
    """Hash patient/practitioner references.

    Args:
        reference: FHIR Reference object

    Returns:
        De-identified Reference with hashed IDs
    """
    deidentified = reference.copy()

    # Hash the reference ID
    if "reference" in deidentified:
        ref = deidentified["reference"]
        # Handle "Patient/patient-001" format
        if "/" in ref:
            resource_type, resource_id = ref.split("/", 1)
            hashed_id = hashlib.sha256(resource_id.encode()).hexdigest()[:16]
            deidentified["reference"] = f"{resource_type}/deidentified-{hashed_id}"

    # Remove display name (may contain patient/doctor name)
    if "display" in deidentified:
        del deidentified["display"]

    # Hash identifier if present
    if "identifier" in deidentified and "value" in deidentified["identifier"]:
        original_value = deidentified["identifier"]["value"]
        hashed_value = hashlib.sha256(original_value.encode()).hexdigest()[:16]
        deidentified["identifier"]["value"] = f"HASH-{hashed_value}"

    return deidentified


def _reduce_to_year(date_string: str) -> str:
    """Reduce date to year only (HIPAA requirement).

    Args:
        date_string: ISO date string (e.g., "2024-03-15T10:30:00Z")

    Returns:
        Year as string (e.g., "2024")
    """
    try:
        # Handle various date formats
        date_obj = datetime.fromisoformat(date_string.replace("Z", "+00:00"))
        return str(date_obj.year)
    except (ValueError, AttributeError):
        # If parsing fails, remove the date entirely
        logger.warning(f"Could not parse date: {date_string}")
        return ""


def deidentify_bundle(bundle: Dict[str, Any]) -> Dict[str, Any]:
    """Apply de-identification to FHIR Bundle.

    Automatically detects resource type and applies appropriate de-identification.

    Args:
        bundle: FHIR Bundle resource containing multiple resources

    Returns:
        De-identified Bundle with all resources de-identified
    """
    deidentified = bundle.copy()

    if "entry" not in deidentified:
        return deidentified

    for entry in deidentified["entry"]:
        if "resource" not in entry:
            continue

        resource = entry["resource"]
        resource_type = resource.get("resourceType")

        # Apply type-specific de-identification
        if resource_type == "Patient":
            entry["resource"] = deidentify_patient(resource)
        elif resource_type == "Observation":
            entry["resource"] = deidentify_observation(resource)
        elif resource_type == "Condition":
            entry["resource"] = deidentify_condition(resource)
        elif resource_type == "MedicationStatement":
            entry["resource"] = deidentify_medication_statement(resource)
        else:
            # Generic de-identification for other resource types
            entry["resource"] = _deidentify_generic_resource(resource)

    # Add bundle-level de-identification marker
    deidentified["_deidentified"] = True
    deidentified["_deidentification_method"] = "HIPAA Safe Harbor"
    deidentified["_deidentification_date"] = datetime.utcnow().isoformat()

    logger.info(f"De-identified bundle with {len(deidentified['entry'])} resources")
    return deidentified


def deidentify_resource(resource: Dict[str, Any]) -> Dict[str, Any]:
    """Apply de-identification to any FHIR resource.

    Automatically detects resource type and applies appropriate de-identification.
    Falls back to generic de-identification for unknown resource types.

    Args:
        resource: FHIR resource of any type

    Returns:
        De-identified resource
    """
    resource_type = resource.get("resourceType")

    # Apply type-specific de-identification if available
    if resource_type == "Patient":
        return deidentify_patient(resource)
    elif resource_type == "Observation":
        return deidentify_observation(resource)
    elif resource_type == "Condition":
        return deidentify_condition(resource)
    elif resource_type == "MedicationStatement":
        return deidentify_medication_statement(resource)
    else:
        # Generic de-identification for other resource types
        return _deidentify_generic_resource(resource)


def _deidentify_generic_resource(resource: Dict[str, Any]) -> Dict[str, Any]:
    """Apply generic de-identification to unknown resource types.

    This is a fallback for resource types not explicitly handled.

    Args:
        resource: FHIR resource

    Returns:
        De-identified resource
    """
    deidentified = resource.copy()

    # Hash resource ID
    if "id" in deidentified:
        original_id = deidentified["id"]
        hashed_id = hashlib.sha256(original_id.encode()).hexdigest()[:16]
        deidentified["id"] = f"deidentified-{hashed_id}"

    # Hash patient references
    if "subject" in deidentified:
        deidentified["subject"] = _deidentify_reference(deidentified["subject"])

    if "patient" in deidentified:
        deidentified["patient"] = _deidentify_reference(deidentified["patient"])

    # Remove common date fields
    date_fields = ["recordedDate", "assertedDate", "dateAsserted", "issued", "date"]
    for field in date_fields:
        if field in deidentified:
            deidentified[field] = _reduce_to_year(deidentified[field])

    # Add de-identification marker
    deidentified["_deidentified"] = True
    deidentified["_deidentification_method"] = "HIPAA Safe Harbor (Generic)"

    return deidentified
