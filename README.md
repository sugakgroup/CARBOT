# CARBOT

**CARBOT** (*π-Conjugated Aromatic Rule-Based Organic Transformation*) is a unit-based molecule generator for **π-conjugated hydrocarbons**.

This repository contains the source code and research materials accompanying the CARBOT paper/preprint (see **Citation** below). For the version used in the manuscript, please use the corresponding tagged release (e.g., `v1.0.0`).

---

## Features

- Unit-based construction of π-conjugated hydrocarbon structures
- Research-focused codebase with scripts/notebooks for generation, benchmarking, and analysis
- Built on RDKit for cheminformatics utilities

---

## Repository layout (high level)

- `src/` — core implementation
- `production/` — scripts/workflows used for producing datasets/results
- `analysis/` — analysis notebooks and figure-generation utilities
- `data/` — input datasets and supporting files
- `evomol_result/` — example/benchmark outputs (EvoMol-related results)

---

## Requirements

Tested with:

- Python `3.12.12`
- RDKit `2025.09.2`

Other versions may work, but have not been tested.

---

## Installation

Clone the repository:

```bash
git clone https://github.com/sugakgroup/CARBOT.git
cd CARBOT
```

### (Recommended) Create a clean environment

Using conda (example):

```bash
conda create -n carbot python=3.12 -y
conda activate carbot
conda install -c conda-forge rdkit=2025.09.2 -y
```

---

## Usage

This repository is currently organized as research code rather than a packaged library.

Typical workflow:

1. Set up the Python/RDKit environment (see above).
2. Run scripts and notebooks under `production/` and `analysis/` to reproduce generation/benchmarking/plots.
3. Use or extend the modules under `src/` in your own workflows.

---

## Reproducibility

- For paper-aligned code, check out the tagged release used for the manuscript (e.g., `v1.0.0`).

---

## License

This project is released under the **MIT License** (see `LICENSE`).

---

## Citation

If you use CARBOT in academic work, please cite:

K. Suga*, H. Takahashi, K. Terayama, M. Sumita, S. Saito, **ChemRxiv** (2026). DOI: `10.26434/chemrxiv-2025-ldqfv-v2`

### BibTeX

```bibtex
@article{suga_carbot_2026,
  author  = {Suga, Kensuke and Takahashi, H. and Terayama, K. and Sumita, M. and Saito, S.},
  title   = {Defining the Complete Chemical Space of π-Conjugated Hydrocarbons by Unit-Based Construction},
  journal = {ChemRxiv},
  year    = {2026},
  doi     = {10.26434/chemrxiv-2025-ldqfv-v2}
}
```

---

## Contact

For questions or collaboration inquiries, please contact the author via e-mail (sugak24[at]chem.sci.osaka-u.ac.jp).
