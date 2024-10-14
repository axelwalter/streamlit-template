from dataclasses import dataclass
from typing import Literal, Union

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import plotly.graph_objects as go

from src.plotting.BasePlotter import Colors, _BasePlotter, _BasePlotterConfig


@dataclass(kw_only=True)
class MSExperimentPlotterConfig(_BasePlotterConfig):
    bin_peaks: Union[Literal["auto"], bool] = "auto"
    num_RT_bins: int = 50
    num_mz_bins: int = 50
    plot3D: bool = False
    title: str = "Peak Map"
    xlabel: str = "RT (s)"
    ylabel: str = "m/z"
    height: int = 500
    width: int = 750
    relative_intensity: bool = False
    show_legend: bool = True


class MSExperimentPlotter(_BasePlotter):
    def __init__(self, config: MSExperimentPlotterConfig, **kwargs) -> None:
        """
        Initialize the MSExperimentPlotter with a given configuration and optional parameters.

        Args:
            config (MSExperimentPlotterConfig): Configuration settings for the spectrum plotter.
            **kwargs: Additional keyword arguments for customization.
        """
        super().__init__(config=config, **kwargs)

    def _prepare_data(self, exp: pd.DataFrame) -> pd.DataFrame:
        """Prepares data for plotting based on configuration (binning, relative intensity, hover text)."""
        if self.config.bin_peaks == True or (
            exp.shape[0] > self.config.num_mz_bins * self.config.num_RT_bins
            and self.config.bin_peaks == "auto"
        ):
            exp["mz"] = pd.cut(exp["mz"], bins=self.config.num_mz_bins)
            exp["RT"] = pd.cut(exp["RT"], bins=self.config.num_RT_bins)

            # Group by x and y bins and calculate the mean intensity within each bin
            exp = (
                exp.groupby(["mz", "RT"], observed=True)
                .agg({"inty": "mean"})
                .reset_index()
            )
            exp["mz"] = exp["mz"].apply(lambda interval: interval.mid).astype(float)
            exp["RT"] = exp["RT"].apply(lambda interval: interval.mid).astype(float)
            exp = exp.fillna(0)
        else:
            self.config.bin_peaks = False

        if self.config.relative_intensity:
            exp["inty"] = exp["inty"] / max(exp["inty"]) * 100

        exp["hover_text"] = exp.apply(
            lambda x: f"m/z: {round(x['mz'], 6)}<br>RT: {round(x['RT'], 2)}<br>intensity: {int(x['inty'])}",
            axis=1,
        )

        return exp.sort_values("inty")

    def _plotMatplotlib3D(
        self,
        exp: pd.DataFrame,
    ) -> plt.Figure:
        """Plot 3D peak map with mz, RT and intensity dimensions. Colored peaks based on intensity."""
        fig = plt.figure(
            figsize=(self.config.width / 100, self.config.height / 100),
            layout="constrained",
        )
        ax = fig.add_subplot(111, projection="3d")

        if self.config.title:
            ax.set_title(self.config.title, fontsize=12, loc="left")
        ax.set_xlabel(
            self.config.ylabel,
            fontsize=9,
            labelpad=-2,
            color=Colors["DARKGRAY"],
            style="italic",
        )
        ax.set_ylabel(
            self.config.xlabel,
            fontsize=9,
            labelpad=-2,
            color=Colors["DARKGRAY"],
        )
        ax.set_zlabel("intensity", fontsize=10, color=Colors["DARKGRAY"], labelpad=-2)
        for axis in ("x", "y", "z"):
            ax.tick_params(axis=axis, labelsize=8, pad=-2, colors=Colors["DARKGRAY"])

        ax.set_box_aspect(aspect=None, zoom=0.88)
        ax.ticklabel_format(axis="z", style="sci", useMathText=True, scilimits=(0,0))
        ax.grid(color="#FF0000", linewidth=0.8)
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        ax.view_init(elev=25, azim=-45, roll=0)

        # Plot lines to the bottom with colored based on inty
        for i in range(len(exp)):
            ax.plot(
                [exp["RT"].iloc[i], exp["RT"].iloc[i]],
                [exp["inty"].iloc[i], 0],
                [exp["mz"].iloc[i], exp["mz"].iloc[i]],
                zdir="x",
                color=plt.cm.magma_r((exp["inty"].iloc[i] / exp["inty"].max())),
            )
        return fig

    def _plotPlotly2D(
        self,
        exp: pd.DataFrame,
    ) -> go.Figure:
        """Plot 2D peak map with mz and RT dimensions. Colored peaks based on intensity."""
        layout = go.Layout(
            title=dict(text=self.config.title),
            xaxis=dict(title=self.config.xlabel),
            yaxis=dict(title=self.config.ylabel),
            showlegend=self.config.show_legend,
            template="simple_white",
            dragmode="select",
            height=self.config.height,
            width=self.config.width,
        )
        fig = go.Figure(layout=layout)
        fig.add_trace(
            go.Scattergl(
                name="peaks",
                x=exp["RT"],
                y=exp["mz"],
                mode="markers",
                marker=dict(
                    color=exp["inty"].apply(lambda x: np.log(x)),
                    colorscale="sunset",
                    size=8,
                    symbol="square",
                    colorbar=(
                        dict(thickness=8, outlinewidth=0, tickformat=".0e")
                        if self.config.show_legend
                        else None
                    ),
                ),
                hovertext=exp["hover_text"] if not self.config.bin_peaks else None,
                hoverinfo="text",
                showlegend=False,
            )
        )
        return fig

    def _plot(
        self,
        exp: pd.DataFrame,
    ) -> go.Figure:
        """Prepares data and returns Plotly 2D plot or Matplotlib 3D plot."""
        exp = self._prepare_data(exp)
        if self.config.plot3D:
            return self._plotMatplotlib3D(exp)
        return self._plotPlotly2D(exp)

# ============================================================================= #
## FUNCTIONAL API ##
# ============================================================================= #


def plotMSExperiment(
    exp: pd.DataFrame,
    plot3D: bool = False,
    relative_intensity: bool = False,
    bin_peaks: Union[Literal["auto"], bool] = "auto",
    num_RT_bins: int = 50,
    num_mz_bins: int = 50,
    width: int = 750,
    height: int = 500,
    title: str = "Peak Map",
    xlabel: str = "RT (s)",
    ylabel: str = "m/z",
    show_legend: bool = False,
):
    """
    Plots a Spectrum from an MSSpectrum object

    Args:
        spectrum (pd.DataFrame): OpenMS MSSpectrum Object
        plot3D: (bool = False, optional): Plot peak map 3D with peaks colored based on intensity. Disables colorbar legend. Works with "MATPLOTLIB" engine only. Defaults to False.
        relative_intensity (bool, optional): If true, plot relative intensity values. Defaults to False.
        bin_peaks: (Union[Literal["auto"], bool], optional): Bin peaks to reduce complexity and improve plotting speed. Hovertext disabled if activated. If set to "auto" any MSExperiment with more then num_RT_bins x num_mz_bins peaks will be binned. Defaults to "auto".
        num_RT_bins: (int, optional): Number of bins in RT dimension. Defaults to 50.
        num_mz_bins: (int, optional): Number of bins in m/z dimension. Defaults to 50.
        width (int, optional): Width of plot. Defaults to 500px.
        height (int, optional): Height of plot. Defaults to 500px.
        title (str, optional): Plot title. Defaults to "Spectrum Plot".
        xlabel (str, optional): X-axis label. Defaults to "m/z".
        ylabel (str, optional): Y-axis label. Defaults to "intensity" or "ion mobility".
        show_legend (int, optional): Show legend. Defaults to False.

    Returns:
        Plot: The generated plot using the specified engine.
    """
    config = MSExperimentPlotterConfig(
        plot3D=plot3D,
        relative_intensity=relative_intensity,
        bin_peaks=bin_peaks,
        num_RT_bins=num_RT_bins,
        num_mz_bins=num_mz_bins,
        width=width,
        height=height,
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
        show_legend=show_legend,
    )
    plotter = MSExperimentPlotter(config)
    return plotter.plot(exp.copy())