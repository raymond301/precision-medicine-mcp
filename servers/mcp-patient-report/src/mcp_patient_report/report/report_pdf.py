"""
PDF Generator - WeasyPrint HTML to PDF conversion.

Converts rendered HTML reports to print-ready PDF documents.
"""

import logging
from pathlib import Path
from typing import Optional
from datetime import datetime

logger = logging.getLogger(__name__)

# Check if WeasyPrint is available
try:
    from weasyprint import HTML, CSS
    from weasyprint.text.fonts import FontConfiguration
    WEASYPRINT_AVAILABLE = True
except ImportError:
    WEASYPRINT_AVAILABLE = False
    logger.warning(
        "WeasyPrint not available. PDF generation will be disabled. "
        "Install with: pip install weasyprint"
    )


class PDFGenerator:
    """Generate PDF documents from HTML using WeasyPrint."""

    def __init__(self, output_dir: Optional[Path] = None):
        """
        Initialize the PDF generator.

        Args:
            output_dir: Directory for output files. Defaults to current directory.
        """
        if not WEASYPRINT_AVAILABLE:
            raise ImportError(
                "WeasyPrint is required for PDF generation. "
                "Install with: pip install weasyprint\n"
                "Note: WeasyPrint requires system libraries (libpango, libcairo, libgdk-pixbuf). "
                "See: https://doc.courtbouillon.org/weasyprint/stable/first_steps.html#installation"
            )

        self.output_dir = output_dir or Path.cwd()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.font_config = FontConfiguration()

        logger.info(f"PDFGenerator initialized. Output directory: {self.output_dir}")

    def generate_pdf(
        self,
        html_content: str,
        output_filename: str,
        base_url: Optional[str] = None,
    ) -> Path:
        """
        Generate a PDF from HTML content.

        Args:
            html_content: Rendered HTML string
            output_filename: Output PDF filename (without path)
            base_url: Base URL for resolving relative paths (images, CSS)

        Returns:
            Path to the generated PDF file
        """
        # Ensure filename ends with .pdf
        if not output_filename.endswith('.pdf'):
            output_filename = f"{output_filename}.pdf"

        output_path = self.output_dir / output_filename

        logger.info(f"Generating PDF: {output_path}")

        # Create HTML document
        html_doc = HTML(
            string=html_content,
            base_url=base_url,
        )

        # Write PDF
        html_doc.write_pdf(
            output_path,
            font_config=self.font_config,
        )

        logger.info(f"PDF generated successfully: {output_path}")
        return output_path

    def generate_report_pdf(
        self,
        html_content: str,
        patient_id: str,
        report_type: str = "full",
        is_draft: bool = True,
    ) -> Path:
        """
        Generate a patient report PDF with standardized naming.

        Args:
            html_content: Rendered HTML string
            patient_id: Patient identifier
            report_type: Type of report (full, onepage)
            is_draft: Whether this is a draft (adds DRAFT to filename)

        Returns:
            Path to the generated PDF file
        """
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        draft_suffix = "_DRAFT" if is_draft else ""

        filename = f"{patient_id}_{report_type}_report{draft_suffix}_{timestamp}.pdf"

        return self.generate_pdf(html_content, filename)


class PDFGeneratorFallback:
    """
    Fallback PDF generator when WeasyPrint is not available.

    Saves HTML output instead of PDF.
    """

    def __init__(self, output_dir: Optional[Path] = None):
        """Initialize fallback generator."""
        self.output_dir = output_dir or Path.cwd()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        logger.warning(
            "Using PDFGeneratorFallback - PDF generation not available. "
            "HTML files will be saved instead."
        )

    def generate_pdf(
        self,
        html_content: str,
        output_filename: str,
        base_url: Optional[str] = None,
    ) -> Path:
        """Save HTML content to file (fallback when PDF not available)."""
        # Change extension to .html
        if output_filename.endswith('.pdf'):
            output_filename = output_filename[:-4] + '.html'
        elif not output_filename.endswith('.html'):
            output_filename = f"{output_filename}.html"

        output_path = self.output_dir / output_filename

        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_content)

        logger.info(f"HTML saved (PDF fallback): {output_path}")
        return output_path

    def generate_report_pdf(
        self,
        html_content: str,
        patient_id: str,
        report_type: str = "full",
        is_draft: bool = True,
    ) -> Path:
        """Generate report file with standardized naming (HTML fallback)."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        draft_suffix = "_DRAFT" if is_draft else ""

        filename = f"{patient_id}_{report_type}_report{draft_suffix}_{timestamp}.html"

        return self.generate_pdf(html_content, filename)


def get_pdf_generator(output_dir: Optional[Path] = None):
    """
    Get the appropriate PDF generator based on available dependencies.

    Args:
        output_dir: Output directory for generated files

    Returns:
        PDFGenerator if WeasyPrint is available, otherwise PDFGeneratorFallback
    """
    if WEASYPRINT_AVAILABLE:
        return PDFGenerator(output_dir)
    else:
        return PDFGeneratorFallback(output_dir)
