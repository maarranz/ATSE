# Applied Time Series Econometrics (ATSE)

Course materials by Miguel A. Arranz (Universidad Carlos III de Madrid):
self-study chapters, slides, R and Python Labs, and capstone projects.

**Contact:** <maarranz@eco.uc3m.es> · **ORCID:** [0000-0002-5951-2284](https://orcid.org/0000-0002-5951-2284)

**Website:** <https://maarranz.github.io/ATSE/>

## Contents

| Folder | What it holds |
|---|---|
| `chapters/` | Self-study chapters (Quarto `.qmd`), one folder per module |
| `slides/` | RevealJS slide decks (Quarto `.qmd`) |
| `typst/` | Chapter versions adapted for the PDF (Typst) output |
| `labs/r/` | R Labs (Quarto `.qmd`) |
| `labs/python/` | Python Labs (Jupyter `.ipynb`, run on Google Colab) |
| `capstones/` | Capstone research projects |
| `resources/` | Links to the course's Python and R packages and web apps |
| `setup/` | Setting up R, Python, Colab and the FRED API key |
| `data/` | Datasets used in the course |
| `docs/` | The rendered website (generated; do not edit by hand) |

## Running the materials

- **R:** open `ATSE.Rproj`, then run `renv::restore()` once to install the exact package versions used.
- **Python:** open a Lab from its "Open in Colab" badge; no local install needed.
  For local use: `conda env create -f environment.yml`.

## Licence

Course materials: [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/).
Code: MIT. See `LICENSE` for details.
