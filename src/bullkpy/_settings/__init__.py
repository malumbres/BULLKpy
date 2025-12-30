from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal, Any


# -------------------------------------------------------------------
# Verbosity (kept compatible with your existing repr + logging usage)
# -------------------------------------------------------------------
@dataclass
class Verbosity:
    """Logging verbosity as a standard logging level integer."""
    level: int = 20  # INFO


def _verbosity_from_int(v: int) -> int:
    """
    Scanpy-like convenience mapping:
      0 -> WARNING (30)
      1 -> INFO    (20)
      2 -> DEBUG   (10)
      3 -> DEBUG   (10)
    """
    v = int(v)
    if v <= 0:
        return 30
    if v == 1:
        return 20
    return 10


PlotTheme = Literal["default", "paper", "talk"]
Where = Literal["obs", "var"]


@dataclass
class Settings:
    # -----------------------------
    # Existing core settings
    # -----------------------------
    _verbosity: Verbosity = field(default_factory=Verbosity)
    figdir: Path = Path("figures")
    autoshow: bool = True
    autosave: bool = False

    # -----------------------------
    # Plot settings (new)
    # -----------------------------
    plot_theme: PlotTheme = "default"
    plot_palette: str = "Set1"
    plot_fontsize: float = 12.0
    plot_dpi: int = 150

    # “Scanpy-like” on-screen behavior
    scale_fonts_with_figsize: bool = True

    # -----------------------------
    # Central categorical colors (new)
    # -----------------------------
    # Cache: {(where, key): {category: "#RRGGBB"}}
    categorical_colors: dict[tuple[str, str], dict[str, str]] = field(default_factory=dict)

    # Optional per-column palette override: {(where, key): "Set2"}
    categorical_palette_overrides: dict[tuple[str, str], str] = field(default_factory=dict)

    # -----------------------------
    # Backwards-compatible verbosity property
    # -----------------------------
    @property
    def verbosity(self) -> Verbosity:
        return self._verbosity

    @verbosity.setter
    def verbosity(self, v: int | Verbosity) -> None:
        # allow bk.settings.verbosity = 3  (like scanpy)
        if isinstance(v, Verbosity):
            self._verbosity = v
        else:
            self._verbosity = Verbosity(level=_verbosity_from_int(int(v)))

    # -----------------------------
    # New helpers
    # -----------------------------
    def set_palette_for(self, key: str, *, where: Where = "obs", palette: str = "Set1") -> None:
        """Assign a palette override for a specific categorical column."""
        if where not in ("obs", "var"):
            raise ValueError("where must be 'obs' or 'var'")
        self.categorical_palette_overrides[(where, key)] = str(palette)

    def reset_colors(self) -> None:
        """Clear cached categorical colors (forces regeneration on next plot)."""
        self.categorical_colors.clear()


# Singleton instance (what bk.settings should point to)
settings = Settings()

__all__ = ["settings", "Settings", "Verbosity"]