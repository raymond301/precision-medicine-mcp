# Building the Book

This book is built using **Quarto**, a professional scientific publishing system with native Mermaid diagram support.

---

## Why Quarto?

Advantages over the previous manual Pandoc + Mermaid approach:

- **Native Mermaid support**: No manual extraction or conversion needed
- **Automatic TOC**: Table of contents handled automatically
- **Professional layouts**: Built for books and scientific publishing
- **Simple workflow**: One command to build
- **No shell escaping issues**: No echo/printf/heredoc problems
- **Better error messages**: Clear feedback when things go wrong

---

## Building Locally

### Prerequisites

**Quarto** (macOS):
```bash
brew install quarto
```

**Quarto** (Linux):
```bash
curl -LO https://quarto.org/download/latest/quarto-linux-amd64.deb
sudo dpkg -i quarto-linux-amd64.deb
```

**TinyTeX**:
```bash
quarto install tinytex
```

### Build the PDF

```bash
cd docs/book
quarto render --to pdf
```

Output: `docs/book/_book/Building-AI-Orchestrated-Precision-Oncology-Systems.pdf`

### Preview HTML version

```bash
cd docs/book
quarto preview
```

Opens in browser at http://localhost:4200

---

## GitHub Actions (CI)

Two workflows handle automated builds:

### On Push to Main — `quarto-book.yml`

Triggers when any of these change on `main`:
- `docs/book/**/*.md`
- `docs/book/images/**`
- `docs/book/_quarto.yml`
- `.github/workflows/quarto-book.yml`

Produces a PDF artifact retained for 90 days. Download it from the **Actions** tab on GitHub.

### On Release Tag — `quarto-release.yml`

Triggers on version tags matching `v*.*.*`. Builds the PDF, creates a GitHub Release, and attaches it.

```bash
git tag v1.0.0 -m "Book v1.0.0 - First complete draft"
git push origin v1.0.0
```

Semantic versioning:
- **Major** (v2.0.0): New chapters or breaking content changes
- **Minor** (v1.1.0): Significant new sections
- **Patch** (v1.0.1): Typo and link fixes

### What the CI does (both workflows)

1. Checkout repository
2. Install Chrome via `browser-actions/setup-chrome@v1` — needed for Mermaid → PNG rendering
3. Install Quarto 1.7.13
4. Install TinyTeX
5. Set `QUARTO_CHROME` env var to the Chrome path
6. Render book to PDF with xelatex
7. Upload PDF as artifact (or attach to a GitHub Release)

---

## Configuration

Book structure and PDF formatting live in `_quarto.yml`:
- Chapter order and part divisions
- Appendices
- PDF engine (xelatex), fonts, margins
- Mermaid rendering (PNG via Chrome)
- Hyperref options (breaklinks, pdfborder)

The LaTeX preamble (`header.tex`) adds:
- URL line-breaking via custom `UrlBreaks` (uses the `url` package already loaded by hyperref)
- Sans-serif default font
- Graphics support

---

## Troubleshooting

### Mermaid diagrams not rendering

Quarto needs Chrome/Chromium to convert Mermaid → PNG. In CI, Chrome is installed via `browser-actions/setup-chrome@v1`. Locally:
```bash
# Option A: let Quarto manage it
quarto install chromium

# Option B: point Quarto at your existing Chrome
export QUARTO_CHROME="$(which google-chrome)"
quarto render --to pdf
```

### LaTeX errors during render

Common pitfalls in markdown source:
- Special characters (`_`, `$`, `#`) outside code blocks need escaping or backtick-wrapping
- Image paths must be relative to `docs/book/` — e.g. `images/screenshots/filename.png`
- Internal links must use the correct depth (see below)

### Relative link depth reference

From any file in `docs/book/`:

| Prefix | Resolves to | Examples |
|--------|------------|----------|
| `./` | `docs/book/` | `./companion-notebooks/`, other chapters, appendices |
| `../` | `docs/` | `../architecture/`, `../for-developers/`, `../deployment/` |
| `../../` | repo root | `../../servers/`, `../../data/`, `../../infrastructure/` |

---

## Disabled Workflows (Archive)

The previous pandoc-based workflows are kept for reference but no longer run:
- `build-book-pdf.yml.disabled` — manual Mermaid CLI conversion + Pandoc
- `release-book.yml.disabled` — release with manual conversion
