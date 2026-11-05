# Installation

BULLKpy requires **Python 3.10 or newer**.

## From PyPI

```bash
pip install bullkpy
```

A plain install includes everything the tutorials need: UMAP, Leiden clustering, GSEA, survival
analysis and network layouts. The optional extras are:

| Extra | Adds |
| --- | --- |
| `survival-ml` | scikit-survival, for `bk.pl.time_dependent_roc()` |
| `notebook` | the Jupyter kernel and JupyterLab, to run the tutorials |

```bash
pip install "bullkpy[notebook]"
```

The extras `umap`, `leiden` and `network` from 0.2 still exist but do nothing, because those packages are now installed
by default. Calling a function whose optional package is missing raises an `ImportError` naming the package.

## Check the installation

```bash
python scripts/smoke_test.py     # from a clone of the repository
```

runs a short analysis on synthetic data (load, QC, PCA, UMAP, Leiden, differential expression, a plot, an `.h5ad`
round trip) and prints the installed version.

## Tutorial data

The tutorials use TCGA lung cancer data. `bk.datasets.tcga_lung()` downloads it once (about 60 MB) into
`~/.cache/bullkpy`; set the environment variable `BULLKPY_DATA` to use another folder.

## From source

```bash
git clone https://github.com/malumbres/BULLKpy.git
cd BULLKpy
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,docs,umap,leiden,network,notebook]"
```

## Verifying the install

```bash
python -c "import bullkpy as bk; print(bk.__version__)"
pytest            # from a source checkout
```

## Upgrading from 0.1.x

0.2.0 contains fixes for defects that silently produced wrong results, plus two
breaking changes. See the
[changelog](https://github.com/malumbres/BULLKpy/blob/main/CHANGELOG.md) before
re-running existing analyses.
