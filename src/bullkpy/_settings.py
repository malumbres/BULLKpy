from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal


PlotTheme = Literal["default", "paper", "talk"]


@dataclass
class Settings:
    """
    Global BULLKpy settings (Scanpy-like).

    These settings control:
      - matplotlib / seaborn defaults
      - font scaling (screen-stable)
      - categorical color mapping (centralized & reusable)
    """

    # -----------------------------
    # Plot appearance
    # -----------------------------
    plot_theme: PlotTheme = "default"   # "default" | "paper" | "talk"
    plot_palette: str = "Set1"
    plot_fontsize: float = 12.0
    plot_dpi: int = 150

    # If True, scale fonts inversely with figure size (Scanpy-like behavior)
    # This is what fixes the "small fig huge font / big fig tiny font" issue
    scale_fonts_with_figsize: bool = True

    # -----------------------------
    # Categorical colors (central)
    # -----------------------------
    # Cached colors: {(where, key): {category: hex}}
    categorical_colors: dict[tuple[str, str], dict[str, str]] = field(
        default_factory=dict
    )

    # Optional per-column palette overrides
    # {(where, key): palette_name}
    categorical_palette_overrides: dict[tuple[str, str], str] = field(
        default_factory=dict
    )

    # -----------------------------
    # Methods
    # -----------------------------
    def set_palette_for(
        self,
        key: str,
        *,
        where: Literal["obs", "var"] = "obs",
        palette: str = "Set1",
    ) -> None:
        """
        Assign a palette to a specific categorical column.
        """
        self.categorical_palette_overrides[(where, key)] = palette

    def reset_colors(self) -> None:
        """
        Clear cached categorical colors (forces regeneration).
        """
        self.categorical_colors.clear()

    def reset(self) -> None:
        """
        Reset all settings to defaults.
        """
        self.plot_theme = "default"
        self.plot_palette = "Set1"
        self.plot_fontsize = 12.0
        self.plot_dpi = 150
        self.scale_fonts_with_figsize = True
        self.categorical_colors.clear()
        self.categorical_palette_overrides.clear()


# Singleton (Scanpy-like)
settings = Settings()