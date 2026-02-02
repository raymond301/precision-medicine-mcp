# Building the Book PDF

This document explains how to build a unified PDF from all book chapters using GitHub Actions.

---

## Overview

Two GitHub Action workflows automatically build PDFs:

1. **`build-book-pdf.yml`**: Builds PDF on every push to `main` (for review)
2. **`release-book.yml`**: Creates GitHub Release with PDF when you tag a version

---

## Automatic Build (On Push)

**Trigger**: Every push to `main` branch that modifies book content

**What it does**:
1. Concatenates all 18 chapters + 4 appendices in correct order
2. Converts markdown to PDF using Pandoc
3. Uploads PDF as artifact (downloadable for 90 days)

**Files monitored**:
- `docs/book/chapter-*.md`
- `docs/book/appendix-*.md`
- `docs/book/chapter-00-introduction.md`

**To download the PDF**:
1. Go to: https://github.com/lynnlangit/precision-medicine-mcp/actions
2. Click on the latest "Build Book PDF" workflow run
3. Scroll down to "Artifacts"
4. Download `book-pdf`

---

## Release Build (On Git Tag)

**Trigger**: Push a version tag like `v1.0.0`

**What it does**:
1. Builds PDF with version number in filename
2. Creates GitHub Release with release notes
3. Attaches PDF to the release
4. Uploads PDF as backup artifact (365 days retention)

### How to Create a Release

**Step 1: Create and push a tag**
```bash
# Tag the current commit
git tag v1.0.0 -m "Book v1.0.0 - First complete draft"

# Push the tag to GitHub
git push origin v1.0.0
```

**Step 2: Wait for workflow**
- GitHub Actions will automatically start the "Release Book PDF" workflow
- Build takes ~5-10 minutes
- Check progress: https://github.com/lynnlangit/precision-medicine-mcp/actions

**Step 3: Release is published**
- Go to: https://github.com/lynnlangit/precision-medicine-mcp/releases
- The PDF is attached: `precision-medicine-mcp-book-v1.0.0.pdf`
- Release notes are auto-generated with book statistics

### Version Numbering

Use semantic versioning:
- **Major** (v2.0.0): Major content changes, new chapters
- **Minor** (v1.1.0): Significant updates, new sections
- **Patch** (v1.0.1): Typo fixes, minor edits

Examples:
```bash
git tag v1.0.0 -m "First complete draft with all 18 chapters"
git tag v1.1.0 -m "Added advanced exercises to Chapter 16"
git tag v1.0.1 -m "Fixed typos in Chapters 3, 7, 10"
git push origin v1.0.0 v1.1.0 v1.0.1
```

---

## Manual Build (Local)

If you want to build the PDF locally without GitHub Actions:

### Prerequisites

Install Pandoc and LaTeX:

**macOS**:
```bash
brew install pandoc
brew install --cask mactex-no-gui
```

**Linux**:
```bash
sudo apt-get install pandoc texlive-xetex texlive-fonts-recommended
```

**Windows**:
- Download Pandoc: https://pandoc.org/installing.html
- Download MiKTeX: https://miktex.org/download

### Build Script

```bash
cd docs/book

# Concatenate all chapters
cat chapter-00-introduction.md \
    chapter-01-the-patientone-story.md \
    chapter-02-the-architecture-problem.md \
    chapter-03-testing-the-hypothesis.md \
    chapter-04-clinical-data.md \
    chapter-05-genomic-foundations.md \
    chapter-06-multiomics-integration.md \
    chapter-07-spatial-transcriptomics.md \
    chapter-08-cell-segmentation.md \
    chapter-09-treatment-response-prediction.md \
    chapter-10-quantum-celltype-fidelity.md \
    chapter-11-imaging-histopathology.md \
    chapter-12-cloud-deployment-gcp.md \
    chapter-13-hospital-production-deployment.md \
    chapter-14-operations-monitoring.md \
    chapter-15-for-researchers.md \
    chapter-16-teaching-precision-medicine.md \
    chapter-17-funding-sustainability.md \
    chapter-18-lessons-learned.md \
    appendix-a-quick-reference.md \
    appendix-b-installation-setup.md \
    appendix-c-patientone-dataset.md \
    appendix-d-bias-and-ethics.md > book-full.md

# Convert to PDF
pandoc book-full.md \
  -o precision-medicine-mcp-book.pdf \
  --pdf-engine=xelatex \
  --toc \
  --toc-depth=2 \
  --number-sections \
  --highlight-style=tango \
  --variable linkcolor=blue \
  --variable urlcolor=blue \
  --variable geometry:margin=1in \
  --variable fontsize=11pt \
  --metadata title="Building AI-Orchestrated Precision Oncology Systems" \
  --metadata author="Lynn Langit" \
  --metadata date="$(date +'%Y-%m-%d')"

echo "✅ PDF created: precision-medicine-mcp-book.pdf"
open precision-medicine-mcp-book.pdf  # macOS
```

---

## Troubleshooting

### Workflow Fails: "No such file"

**Cause**: Chapter filename mismatch in workflow

**Fix**: Check that all chapter files exist with exact names:
```bash
cd docs/book
ls -1 chapter-*.md appendix-*.md
```

### PDF Rendering Issues

**Mermaid diagrams not showing**:
- Mermaid diagrams are rendered as code blocks in PDF
- For production, consider pre-rendering mermaid to PNG images

**Images not found**:
- Ensure image paths are relative: `../../data/images/file.png`
- GitHub Actions runs from repo root, not `docs/book/`

**LaTeX errors**:
- Check for special characters that need escaping: `_`, `$`, `#`
- Use backticks for code: `` `example` ``

### Workflow Doesn't Trigger

**On push but no trigger**:
- Check that changed files match the `paths:` filter
- Verify you pushed to `main` branch: `git branch`

**On tag but no trigger**:
- Check tag format: must match `v*.*.*` (e.g., `v1.0.0`)
- Verify tag was pushed: `git push origin v1.0.0`

---

## PDF Customization

### Change PDF Styling

Edit the YAML frontmatter in the workflow files:

```yaml
---
title: "Your Custom Title"
subtitle: "Your Subtitle"
author: "Your Name"
fontsize: 12pt              # Change font size
geometry: margin=1.5in      # Change margins
mainfont: "Times New Roman" # Change font
---
```

### Add Cover Page

Create `docs/book/cover.md`:
```markdown
\begin{titlepage}
\centering
\vspace*{2cm}
{\Huge\bfseries Building AI-Orchestrated\\Precision Oncology Systems\par}
\vspace{1cm}
{\Large Lynn Langit\par}
\vspace{2cm}
{\large Version 1.0.0\par}
{\large February 2026\par}
\end{titlepage}
```

Update workflow to include:
```bash
cat cover.md chapter-00-introduction.md ... > book-full.md
```

### Generate EPUB

Change the conversion command in workflow:
```bash
pandoc book-full.md \
  -o precision-medicine-mcp-book.epub \
  --toc \
  --toc-depth=2 \
  --epub-cover-image=cover.jpg
```

---

## Book Statistics

After build completes, the workflow outputs:
- Total lines
- Total words
- PDF file size
- Number of chapters
- Number of appendices

Example output:
```
📊 Book Statistics:
Total lines: 15234
Total words: 98567
PDF size: 4.2M
Chapters: 18
Appendices: 4
```

---

## Next Steps

1. **Test the workflow**: Push a small change to a chapter and verify PDF builds
2. **Create first release**: Tag v1.0.0 when ready for publication
3. **Share the PDF**: Download from GitHub Releases and distribute

---

## Additional Resources

- **Pandoc Documentation**: https://pandoc.org/MANUAL.html
- **GitHub Actions**: https://docs.github.com/en/actions
- **LaTeX Packages**: https://www.ctan.org/

---

**Last Updated**: 2026-02-01
