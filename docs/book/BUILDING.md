# Building the Book

This book is built using **Quarto**, a professional scientific publishing system with native Mermaid diagram support.

## Why Quarto?

**Advantages over manual Pandoc + Mermaid conversion:**
- ✅ **Native Mermaid support**: No manual extraction/conversion needed
- ✅ **Automatic TOC**: Table of contents handled automatically
- ✅ **Professional layouts**: Built for books and scientific publishing
- ✅ **Simple workflow**: One command to build
- ✅ **No shell escaping issues**: No echo/printf/heredoc problems
- ✅ **Better error messages**: Clear feedback when things go wrong

## Building Locally

### Prerequisites

Install Quarto:
```bash
# macOS
brew install quarto

# Linux
curl -LO https://quarto.org/download/latest/quarto-linux-amd64.deb
sudo dpkg -i quarto-linux-amd64.deb

# Windows
# Download installer from https://quarto.org/docs/get-started/
```

### Build the PDF

```bash
cd docs/book
quarto render --to pdf
```

The PDF will be created at: `docs/book/_book/Building-AI-Orchestrated-Precision-Oncology-Systems.pdf`

### Preview HTML version

```bash
cd docs/book
quarto preview
```

Opens in browser at http://localhost:4200

## GitHub Actions

The book builds automatically on:
- **Push to main** (when book files change): `.github/workflows/quarto-book.yml`
- **Version tags** (v*.*.* ): `.github/workflows/quarto-release.yml`

## Configuration

Book structure is defined in `_quarto.yml`:
- Chapter order
- Part divisions
- Appendices
- PDF formatting options

## Troubleshooting

### Mermaid diagrams not rendering
Quarto handles this automatically. If diagrams don't appear:
```bash
quarto install chromium
```

### PDF not building
Check TinyTeX is installed:
```bash
quarto install tinytex
```

### Image paths not working
Images are resolved from the `docs/book/` directory. Use relative paths:
```markdown
![Caption](images/screenshots/filename.png)
```

## Old Workflows (Disabled)

Previous pandoc-based workflows have been disabled:
- `build-book-pdf.yml.disabled` - Manual Mermaid conversion approach
- `release-book.yml.disabled` - Release with manual conversion

These are kept for reference but no longer run (archive).
