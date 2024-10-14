from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum
from typing import Literal, List
import numpy as np

# A colorset suitable for color blindness
class Colors(str, Enum):
    BLUE = "#4575B4"
    RED = "#D73027"
    LIGHTBLUE = "#91BFDB"
    ORANGE = "#FC8D59"
    PURPLE = "#7B2C65"
    YELLOW = "#FCCF53"
    DARKGRAY = "#555555"
    LIGHTGRAY = "#BBBBBB"


@dataclass(kw_only=True)
class _BasePlotterConfig(ABC):
    title: str = "1D Plot"
    xlabel: str = "X-axis"
    ylabel: str = "Y-axis"
    height: int = 500
    width: int = 500
    relative_intensity: bool = False
    show_legend: bool = True


# Abstract Class for Plotting
class _BasePlotter(ABC):
    def __init__(self, config: _BasePlotterConfig) -> None:
        self.config = config
        self.fig = None  # holds the figure object

    def updateConfig(self, **kwargs):
        for key, value in kwargs.items():
            if hasattr(self.config, key):
                setattr(self.config, key, value)
            else:
                raise ValueError(f"Invalid config setting: {key}")

    def _get_n_grayscale_colors(self, n: int) -> List[str]:
        """Returns n evenly spaced grayscale colors in hex format."""
        hex_list = []
        for v in np.linspace(50, 200, n):
            hex = "#"
            for _ in range(3):
                hex += f"{int(round(v)):02x}"
            hex_list.append(hex)
        return hex_list

    def plot(self, data, **kwargs):
        return self._plot(data, **kwargs)

    @abstractmethod
    def _plot(self, data, **kwargs):
        pass