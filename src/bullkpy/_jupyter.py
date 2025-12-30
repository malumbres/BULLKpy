# bullkpy/_jupyter.py
from __future__ import annotations


def in_jupyter() -> bool:
    """Return True if running inside a Jupyter notebook or lab."""
    try:
        from IPython import get_ipython

        ip = get_ipython()
        if ip is None:
            return False
        # ZMQInteractiveShell = Jupyter
        return ip.__class__.__name__ == "ZMQInteractiveShell"
    except Exception:
        return False


def apply_inline_dpi_defaults():
    """
    Apply Scanpy-like inline backend settings:
    - screen DPI ~100
    - save DPI ~300
    Safe to call multiple times.
    """
    try:
        from IPython import get_ipython

        ip = get_ipython()
        if ip is None:
            return

        # Only touch inline backend
        ip.run_line_magic("matplotlib", "inline")

        # Do NOT overwrite user custom rc if already set
        cfg = ip.config
        rc = getattr(cfg, "InlineBackend", {}).get("rc", {})

        rc.setdefault("figure.dpi", 100)
        rc.setdefault("savefig.dpi", 300)

        cfg.InlineBackend.rc = rc
        cfg.InlineBackend.figure_format = "png"

    except Exception:
        # Fail silently — never break imports
        pass