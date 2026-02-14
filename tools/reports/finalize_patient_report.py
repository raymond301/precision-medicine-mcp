#!/usr/bin/env python3
"""
Finalize patient report after Clinician-in-the-Loop (CitL) approval.

This script takes a draft clinical report and a signed clinician review,
validates the approval status, and generates a final approved report that
is ready for clinical decision-making. The final report includes the
clinician's attestation and guideline compliance assessment.

Usage:
    python finalize_patient_report.py --patient-id PAT001-OVC-2025 --output-dir ./results

Author: Claude Sonnet 4.5
Date: 2026-01-13
Part of: Precision Medicine MCP - CitL Validation Workflow
"""

import argparse
import json
import sys
from pathlib import Path
from datetime import datetime
from typing import Dict, Any


def load_draft_report(patient_output_dir: Path) -> Dict[str, Any]:
    """
    Load draft report JSON.

    Args:
        patient_output_dir: Path to patient output directory

    Returns:
        Draft report dictionary

    Raises:
        FileNotFoundError: If draft report not found
    """
    draft_file = patient_output_dir / "draft_report.json"
    if not draft_file.exists():
        raise FileNotFoundError(f"Draft report not found: {draft_file}")

    with open(draft_file) as f:
        return json.load(f)


def load_signed_review(patient_output_dir: Path) -> Dict[str, Any]:
    """
    Load signed review JSON (most recent if multiple exist).

    Args:
        patient_output_dir: Path to patient output directory

    Returns:
        Signed review dictionary

    Raises:
        FileNotFoundError: If no signed review found
    """
    review_files = list(patient_output_dir.glob("*_signed.json"))
    if not review_files:
        raise FileNotFoundError(
            f"No signed review found in {patient_output_dir}\n"
            "   Expected file pattern: *_signed.json\n"
            "   Run citl_submit_review.py first to create signed review"
        )

    # Get most recent review (by filename sort, or could use file mtime)
    review_file = sorted(review_files)[-1]

    print(f"📂 Loading signed review: {review_file.name}")

    with open(review_file) as f:
        return json.load(f)


def generate_final_report(draft_report: Dict[str, Any],
                          signed_review: Dict[str, Any]) -> Dict[str, Any]:
    """
    Generate final approved report from draft and signed review.

    Combines the draft report with clinician validation data to create
    a final approved report ready for clinical use.

    Args:
        draft_report: Draft report dictionary
        signed_review: Signed review dictionary

    Returns:
        Final approved report dictionary
    """
    # Start with draft report
    final_report = draft_report.copy()

    # Update metadata with approval information
    final_report['report_metadata']['status'] = 'clinically_approved'
    final_report['report_metadata']['approval_date'] = datetime.now().isoformat()
    final_report['report_metadata']['reviewer'] = signed_review['reviewer']['name']
    final_report['report_metadata']['reviewer_credentials'] = signed_review['reviewer']['credentials']
    final_report['report_metadata']['review_decision'] = signed_review['decision']['status']
    final_report['report_metadata']['review_rationale'] = signed_review['decision']['rationale']

    # Add version tracking
    final_report['report_metadata']['draft_version'] = draft_report['report_metadata'].get('version', '1.0')
    final_report['report_metadata']['final_version'] = '1.0-approved'

    # Add guideline compliance from clinician review
    final_report['guideline_compliance'] = signed_review.get('guideline_compliance', {})

    # Add per-finding validation results
    final_report['findings_validation'] = signed_review.get('per_finding_validation', [])

    # Add treatment recommendations review
    final_report['treatment_recommendations_review'] = signed_review.get('treatment_recommendations_review', [])

    # Add quality flags assessment
    final_report['quality_flags_clinician_assessment'] = signed_review.get('quality_flags_assessment', [])

    # Add clinical attestation (digital signature)
    final_report['clinical_attestation'] = signed_review.get('attestation', {})

    # Add audit trail metadata
    final_report['audit_trail'] = {
        'review_date': signed_review.get('review_date'),
        'review_id': signed_review.get('attestation', {}).get('signature_hash'),
        'reviewer_email_hash': signed_review.get('reviewer', {}).get('email', 'unknown'),
        'approval_timestamp': datetime.now().isoformat()
    }

    return final_report


def save_final_report(final_report: Dict[str, Any],
                      patient_output_dir: Path) -> Path:
    """
    Save final approved report to file.

    Args:
        final_report: Final approved report dictionary
        patient_output_dir: Path to patient output directory

    Returns:
        Path to saved final report file
    """
    final_file = patient_output_dir / "final_report_approved.json"
    with open(final_file, 'w') as f:
        json.dump(final_report, f, indent=2)

    return final_file


def print_summary(final_report: Dict[str, Any], final_file: Path):
    """Print a summary of the finalized report."""
    metadata = final_report['report_metadata']

    print(f"\n{'='*70}")
    print(f"Final Report Summary")
    print(f"{'='*70}")
    print(f"Patient ID:           {metadata['patient_id']}")
    print(f"Status:               {metadata['status']}")
    print(f"Reviewer:             {metadata['reviewer']} ({metadata.get('reviewer_credentials', 'N/A')})")
    print(f"Approval Date:        {metadata['approval_date']}")
    print(f"Review Decision:      {metadata['review_decision']}")
    print(f"Draft Version:        {metadata.get('draft_version', 'N/A')}")
    print(f"Final Version:        {metadata.get('final_version', 'N/A')}")

    # Guideline compliance
    gc = final_report.get('guideline_compliance', {})
    if gc:
        print(f"\nGuideline Compliance:")
        print(f"  NCCN Aligned:       {gc.get('nccn_aligned', 'Not assessed')}")
        print(f"  Institutional:      {gc.get('institutional_aligned', 'Not assessed')}")

    # Findings validation summary
    findings_val = final_report.get('findings_validation', [])
    if findings_val:
        confirmed = sum(1 for f in findings_val if f.get('validation_status') == 'CONFIRMED')
        uncertain = sum(1 for f in findings_val if f.get('validation_status') == 'UNCERTAIN')
        incorrect = sum(1 for f in findings_val if f.get('validation_status') == 'INCORRECT')
        print(f"\nFindings Validation:")
        print(f"  Total Validated:    {len(findings_val)}")
        print(f"  - Confirmed:        {confirmed}")
        print(f"  - Uncertain:        {uncertain}")
        print(f"  - Incorrect:        {incorrect}")

    # Quality checks
    qc = final_report.get('quality_checks', {})
    if qc:
        print(f"\nQuality Checks:")
        print(f"  All Passed:         {qc.get('all_checks_passed', 'N/A')}")
        flags_count = len(qc.get('flags', []))
        if flags_count > 0:
            print(f"  Flags Raised:       {flags_count}")

    # Audit trail
    audit = final_report.get('audit_trail', {})
    if audit.get('review_id'):
        print(f"\nAudit Trail:")
        print(f"  Review ID:          {audit['review_id'][:32]}...")
        print(f"  Approval Timestamp: {audit.get('approval_timestamp', 'N/A')}")

    print(f"\nFinal Report Saved: {final_file}")
    print(f"{'='*70}\n")


def main():
    """Main entry point for report finalization."""
    parser = argparse.ArgumentParser(
        description='Finalize patient report after CitL approval',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Finalize report after approval
    python finalize_patient_report.py --patient-id PAT001-OVC-2025 --output-dir ./results

    # Use custom output directory
    python finalize_patient_report.py --patient-id PAT001-OVC-2025 \\
        --output-dir /path/to/custom/results
        """
    )
    parser.add_argument('--patient-id', required=True,
                       help='Patient identifier (e.g., PAT001-OVC-2025)')
    parser.add_argument('--output-dir', default='./results',
                       help='Base output directory containing patient subdirectories (default: ./results)')
    parser.add_argument('--force', action='store_true',
                       help='Force finalization even if review status is not APPROVE (use with caution)')

    args = parser.parse_args()

    patient_output_dir = Path(args.output_dir) / args.patient_id

    print(f"📋 Report Finalization")
    print(f"{'='*70}")
    print(f"Patient ID:         {args.patient_id}")
    print(f"Output Directory:   {patient_output_dir}")
    print(f"{'='*70}\n")

    # Verify output directory exists
    if not patient_output_dir.exists():
        print(f"❌ Error: Patient output directory not found: {patient_output_dir}", file=sys.stderr)
        print(f"   Expected: {patient_output_dir}", file=sys.stderr)
        print(f"   Run generate_patient_report.py first to create draft report", file=sys.stderr)
        sys.exit(1)

    # Load draft report
    print(f"📂 Loading draft report...")
    try:
        draft_report = load_draft_report(patient_output_dir)
        print(f"   Status: {draft_report['report_metadata']['status']}")
    except FileNotFoundError as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        sys.exit(1)
    except (json.JSONDecodeError, KeyError) as e:
        print(f"❌ Error loading draft report: {e}", file=sys.stderr)
        sys.exit(1)

    # Load signed review
    print(f"📂 Loading signed review...")
    try:
        signed_review = load_signed_review(patient_output_dir)
    except FileNotFoundError as e:
        print(f"❌ Error: {e}", file=sys.stderr)
        sys.exit(1)
    except (json.JSONDecodeError, KeyError) as e:
        print(f"❌ Error loading signed review: {e}", file=sys.stderr)
        sys.exit(1)

    # Verify approval status
    review_status = signed_review.get('decision', {}).get('status')
    reviewer_name = signed_review.get('reviewer', {}).get('name', 'Unknown')

    print(f"   Reviewer: {reviewer_name}")
    print(f"   Decision: {review_status}")

    if review_status != 'APPROVE':
        print(f"\n⚠️  Warning: Review status is '{review_status}', not 'APPROVE'", file=sys.stderr)

        if review_status == 'REVISE':
            print(f"\n📋 Revision Instructions:", file=sys.stderr)
            revision_inst = signed_review.get('revision_instructions', {})
            issues = revision_inst.get('issues_to_address', [])
            if issues:
                for i, issue in enumerate(issues, 1):
                    print(f"   {i}. {issue}", file=sys.stderr)
            else:
                print(f"   No specific issues documented", file=sys.stderr)

            print(f"\n❌ Cannot finalize: Address revision instructions and resubmit for review", file=sys.stderr)
            sys.exit(1)

        elif review_status == 'REJECT':
            print(f"\n❌ Cannot finalize: Report was rejected by clinician", file=sys.stderr)
            print(f"   Review rationale: {signed_review.get('decision', {}).get('rationale', 'Not provided')}", file=sys.stderr)
            sys.exit(1)

        else:
            if not args.force:
                print(f"\n❌ Cannot finalize: Unknown review status '{review_status}'", file=sys.stderr)
                print(f"   Use --force to finalize anyway (not recommended)", file=sys.stderr)
                sys.exit(1)
            else:
                print(f"⚠️  Proceeding with finalization (--force flag used)", file=sys.stderr)

    print(f"\n✅ Review status: APPROVED")
    print(f"   Reviewer: {reviewer_name}")

    # Generate final report
    print(f"\n📝 Generating final approved report...")
    final_report = generate_final_report(draft_report, signed_review)

    # Save final report
    print(f"💾 Saving final report...")
    final_file = save_final_report(final_report, patient_output_dir)

    # Print summary
    print_summary(final_report, final_file)

    print(f"✅ Report finalization complete!")
    print(f"\n📍 Next Steps:")
    print(f"   1. Review final report: {final_file}")
    print(f"   2. Present findings at tumor board")
    print(f"   3. Document clinical decision in patient chart")
    print(f"   4. Archive audit trail (10-year retention)")


if __name__ == "__main__":
    main()
