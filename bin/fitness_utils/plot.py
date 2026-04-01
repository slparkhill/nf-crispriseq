"""Plotting utilities."""

from typing import List, Optional, Union

from carabiner.mpl import add_legend, grid
import pandas as pd
import seaborn as sns

def plot_guide_data(
    df: pd.DataFrame,
    x: str,
    y: str,
    panel_group: Union[str, List[str]],
    group_title: str = "{}",
    panel_size: float = 3.5,
    control_prefix: str = "controls:",
    spike_prefix: str = "spike",
    guide_name_column: str = "guide_name",
    replicate_column: str = "replicate",
    plot_lines: bool = False,
    logx: bool = False,
    logy: bool = False,
    highlight: Optional = None
):
    """Plots violin, strip-plot, and optionally connected lines side-by-side.

    Plotting connected lines takes a minute or two, so setting `plot_lines=False`
    allows for quick checks.
    """
    if highlight is None:
        highlight = {}
        aspect_ratio = 1.
    else:
        aspect_ratio = 1.2
        
    fig, axes = grid(
        ncol=3 if plot_lines else 2,
        nrow=len(set(df[panel_group])),
        panel_size=3.5,
        squeeze=False,
        aspect_ratio=aspect_ratio,
    )
    grouping = [guide_name_column] + ([replicate_column] if replicate_column in df else [])
    
    for axrow, (group_name, group_df) in zip(axes, df.groupby(panel_group)):
        sns.violinplot(
            data=group_df,
            x=x,
            y=y,
            ax=axrow[0],
            native_scale=True,
            # log_scale=(logx, logy),
        )
        sns.stripplot(
            data=group_df,
            x=x,
            y=y,
            ax=axrow[1],
            size=1.,
            native_scale=True,
            log_scale=(logx, logy),
            color="lightgrey",
        )
        sns.stripplot(
            data=group_df.query(f"{guide_name_column}.str.startswith(@control_prefix)"),
            x=x,
            y=y,
            ax=axrow[1],
            size=1.,
            color="dimgrey",
            native_scale=True,
            log_scale=(logx, logy),
        )
        sns.stripplot(
            data=group_df.query(f"{guide_name_column}.str.startswith(@spike_prefix)"),
            x=x,
            y=y,
            ax=axrow[1],
            size=3.,
            color="C1",
            native_scale=True,
            log_scale=(logx, logy),
        )
        if plot_lines:
            for guide_name, guide_df in group_df.groupby(grouping):
                group_names = dict(zip(grouping, guide_name))
                is_control = group_names[guide_name_column].startswith(control_prefix)
                is_spike = group_names[guide_name_column].startswith(spike_prefix)
                if is_control:
                    color = "dimgrey"
                elif is_spike:
                    color = "C1"
                else:
                    color = "lightgrey"
                axrow[2].plot(
                    x, y,
                    data=guide_df.sort_values(x),
                    c=color,
                    zorder=5 if is_control or is_spike else 0,
                    alpha=.3 if not is_spike else 1.,
                    label="_none",
                )
            axrow[2].set(
                xscale="log" if logx else "linear",
                yscale='log' if logy else "linear",
            )
        total_i = 0
        for key, vals in highlight.items():
            for v in vals:
                this_data = group_df.query(f"{key} == @v")
                sns.stripplot(
                    data=this_data,
                    x=x,
                    y=y,
                    ax=axrow[1],
                    size=3.,
                    color=f"C{1 + total_i}",
                    native_scale=True,
                    log_scale=(logx, logy),
                )
                if plot_lines:
                    this_data2 = this_data.groupby([guide_name_column, x])[[y]].mean().reset_index()
                    for i, (guide_name, guide_df) in enumerate(this_data2.groupby(guide_name_column)):
                        axrow[2].plot(
                            x, y,
                            data=guide_df.sort_values(x),
                            color=f"C{1 + total_i}",
                            zorder=6,
                            label=v if i == 0 else "_none",
                        )
                total_i += 1
        if len(highlight) > 0:
            add_legend(axrow[2])
        for ax in axrow:
            ax.axhline(1., color="lightgrey", zorder=-1)
            ax.set(
                title=group_title.format(group_name),
                xlabel=x,
                ylabel=y,
            )
    return fig, axes