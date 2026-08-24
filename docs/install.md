# Installation

BULLKpy requires **Python 3.10 or newer**.

## From PyPI

```bash
pip install bullkpy
```

## Optional extras

Some functions rely on heavier dependencies, kept optional so the base install
stays light:

| Extra | Enables |
| --- | --- |
| `umap` | `bk.tl.umap()`, `bk.pl.umap()` |
| `leiden` | `bk.tl.cluster(method="leiden")` |
| `network` | network layouts in the GSEA leading-edge plots |
| `notebook` | Jupyter kernel and the tutorial dependencies |

```bash
pip install "bullkpy[umap,leiden]"
```

Calling a function whose extra is missing raises an `ImportError` naming the
package to install.

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
