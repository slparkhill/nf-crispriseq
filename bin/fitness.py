#!/usr/bin/env python

"""Calculate fitness from guide counts."""

from argparse import FileType, Namespace
from functools import wraps
import os
import sys

from carabiner import (
    print_err,
    clicommand, 
    CLIOption, 
    CLICommand, 
    CLIApp,
)

__version__ = "0.0.1"


def column_check(df, required) -> None:
    missing_columns = set(required) - set(df.columns)
    if len(missing_columns) > 0:
        raise ValueError(f"Missing required columns: {', '.join(sorted(missing_columns))}")
    else:
        return None


def row_check(f):

    @wraps(f)
    def _f(df, *args, **kwargs):

        initial_rows = df.shape[0] 
        print_err(f"[INFO] Rows before {f}: {initial_rows=}")

        df = f(df, *args, **kwargs)

        after_rows = df.shape[0]
        print_err(f"[INFO] Rows after {f}: {after_rows=}")

        assert initial_rows >= after_rows

        return df

    return _f


def get_intitial_timepoint_data(
    df,
    timepoint_column: str,
    initial_timepoint=None
):
    initial_timepoint = initial_timepoint or df.query(f"{timepoint_column} >= 0.")[timepoint_column].min()
    return df.query(f"{timepoint_column} == @initial_timepoint")


def filter_to_initial_guides(
    df, 
    guide_name_column: str,
    timepoint_column: str,
    count_column: str,
    min_count: int = 1,
    initial_timepoint=None
):
    import numpy as np

    total_guides = df[guide_name_column].nunique()
    print_err(f"[INFO] Number of unique guide names: {total_guides=}")
    data_t0 = get_intitial_timepoint_data(
        df,
        timepoint_column=timepoint_column,
        initial_timepoint=initial_timepoint,
    ).query(f"{count_column} >= @min_count")

    all_guides = df[guide_name_column].unique()
    sample_id_column = "sample_id"
    n_samples = df[sample_id_column].nunique()
    # most_samples = int(np.floor(n_samples * .9))
    # print_err(f"{most_samples=}")
    # present_guides = sorted(
    #     df.query(f"{count_column} >= @min_count")
    #     .groupby(guide_name_column)
    #     [[sample_id_column]]
    #     .nunique()
    #     .reset_index()
    #     .query(f"{sample_id_column} >= @most_samples")
    #     [guide_name_column]
    #     .unique()
    # )
    present_guides = sorted(data_t0[guide_name_column].unique())
    absent_guides = sorted(
        set(all_guides) 
        - set(present_guides)
    )
    print_err(
        f"[INFO] These {len(absent_guides)} guides absent at time-0:",
        ", ".join(absent_guides[:3]),
        "..." if len(absent_guides) > 3 else "",
    )
    print_err(f"[INFO] Number of guides retained: {len(present_guides)}")
    print_err(
        "[INFO] Percentage of guides retained:",
        f"{(100. * len(present_guides) / total_guides):.2f}",
    )

    return df.query(f"{guide_name_column}.isin(@present_guides)")


def get_count_columns(df, count_column=None, suffixes=("read_count", "umi_count")):
    return [
        col for col in df 
        if col.endswith(suffixes) or col == (count_column or False)
    ]


def generate_reference_columns(
    df,
    reference_guide_name: str,
    guide_name_column: str,
    timepoint_column: str,
    count_column: str,
    sample_id_column: str,
    initial_timepoint=None
):
    import pandas as pd
    count_cols = get_count_columns(df, count_column=count_column)
    def renamer(df, suffix):
        return {
            col: f"{col}_{suffix}" 
            for col in get_count_columns(df, count_column=count_column)
        }

    data_t0 = get_intitial_timepoint_data(
        df,
        timepoint_column=timepoint_column,
        initial_timepoint=initial_timepoint,
    )

    WT_counts = (
        df
        .query(f"{guide_name_column}.str.startswith(@reference_guide_name)")
    )
    print_err(f"[INFO] Number of reference guides retained: {WT_counts[guide_name_column].nunique()=}")

    WT_counts = (
        WT_counts
        .groupby(sample_id_column)
        [count_cols]
        .sum() 
        .rename(columns=renamer(df, "WT"))
        .reset_index()
    )
    print_err(f"[INFO] {WT_counts.shape=}")

    input_counts = (
        data_t0
        [["replicate", guide_name_column] + count_cols]
        .drop_duplicates()
        .rename(columns=renamer(df, "t0"))
    )
    print_err(f"[INFO] {input_counts.shape=}")

    input_counts_WT = (
        data_t0
        .query(f"{guide_name_column}.str.startswith(@reference_guide_name)")
        .groupby("replicate")
        [count_cols]
        .sum()
        .rename(columns=renamer(df, "WT_t0"))
        .reset_index()
    )
    print_err(f"[INFO] {input_counts_WT.shape=}")

    for merge, check in (
        ({"right": WT_counts, "how": "left"}, True),
        ({"right": input_counts, "how": "right"}, False),  # guides could be present at input but not later
        ({"right": input_counts_WT, "how": "left"}, True),
    ):
        df = (row_check(pd.merge) if check else pd.merge)(df, **merge)

    return df.assign(
        ratio_to_wt=lambda x: x[count_column].fillna(0.) / x[f"{count_column}_WT"],
        input_ratio_to_wt=lambda x: x[f"{count_column}_t0"] / x[f"{count_column}_WT_t0"],
        fold_change=lambda x: x["ratio_to_wt"] / x["input_ratio_to_wt"],
        rank_fold_change=lambda x: x.groupby(sample_id_column)["fold_change"].rank(pct=True),
    )


@clicommand(
    message='Calculating fitness with the following parameters',
    name="Fitness",
)
def _fitness(args: Namespace) -> None:

    from carabiner.mpl import figsaver, grid, scattergrid
    from carabiner.pd import print_err, read_table
    import numpy as np
    import pandas as pd

    from fitness_utils.plot import plot_guide_data
    from fitness_utils.models import fit_model

    if args.spike is None and not os.path.exists(args.growth.name):
        raise ValueError(f"Spike not used, and {args.growth} doesn't exist!")

    min_umis = 1
    reference_guide_name = args.reference
    guide_name_column = args.guide_column
    sample_id_column = args.sample_column
    count_column = args.count_column
    timepoint_column = args.timepoint_column

    required_columns = [
        guide_name_column,
        sample_id_column,
        count_column,
        timepoint_column,
    ]

    data = pd.concat(map(read_table, args.inputs), axis=0)
    conditions = read_table(args.conditions)
    
    data = row_check(pd.merge)(
        data,
        conditions,
        how="inner",
    )

    if args.plot:
        figsave = figsaver(
            format=args.plot_format,
            output_dir=args.plot,
        )
        grouping = [guide_name_column, timepoint_column, "replicate"]
        if args.concentration_column:
            grouping += [args.concentration_column]
        data_rep = (
            data
            .groupby(grouping)
            [[count_column]]
            .sum()
            .reset_index()
            .pivot(
                index=[guide_name_column, timepoint_column] + ([args.concentration_column] if args.concentration_column else []),
                columns="replicate",
                values=count_column,
            )
        )
        fig, axes = scattergrid(
            data_rep.reset_index(),
            grid_columns=data_rep.columns.to_list(),
            log=data_rep.columns.to_list(),
            grouping=timepoint_column if data_rep.reset_index()[timepoint_column].nunique() > 1 else None,
        )
        for i, axrow in enumerate(axes):
            for j, ax in enumerate(axrow):
                if i != j:
                    r = np.corrcoef(data_rep.iloc[:,[i,j]].dropna().values.T)[0,-1]
                    ax.plot(
                        ax.get_xlim(), 
                        ax.get_xlim(), 
                        color="lightgrey", 
                        zorder=-1,
                    )
                    ax.set(title=f"r = {r:.2f}")
        figsave(
            fig, 
            name=f"replicates-{count_column}", 
            df=data_rep.reset_index(),
        )

    column_check(data, required_columns)

    data = filter_to_initial_guides(
        data, 
        guide_name_column=guide_name_column,
        timepoint_column=timepoint_column,
        count_column=count_column,
        min_count=min_umis,
    )

    data = generate_reference_columns(
        data, 
        reference_guide_name=reference_guide_name,
        guide_name_column=guide_name_column,
        timepoint_column=timepoint_column,
        count_column=count_column,
        sample_id_column=sample_id_column,
    )

    if args.plot:
        figsave = figsaver(
            format=args.plot_format,
            output_dir=args.plot,
        )
        plots_to_make = {
            count_column: {"y": count_column},
            "fold_change": {"y": "fold_change"},
        }
        for name, plot_args in plots_to_make.items():
            fig, axes = plot_guide_data(
                data,
                x=timepoint_column,
                **plot_args,
                logy=True,
                panel_group=args.concentration_column,
                replicate_column="replicate",
                group_title="{} µM",
                plot_lines=True,
            )
            figsave(
                fig, 
                name=f"violin_{name}", 
                df=data,
            )

    if args.spike is not None:
        expansion_data = (
            data
            .query(f"{guide_name_column}.str.startswith(@args.spike)")
            [[sample_id_column, count_column, "fold_change"]]
        )
        expansion_data = (
            row_check(pd.merge)(
                expansion_data,
                conditions,
                how="inner",
            )
            .rename(columns={
                "fold_change": "fold_change_spike", 
                count_column: f"{count_column}_spike",
            })
        )
        ref_col, expansion_col = f"{count_column}_spike", "fold_change_spike"

    else:
        expansion_data = (
            row_check(pd.merge)(
                read_table(args.growth),
                conditions,
                how="inner",
            )
            .assign(
                growth_t0=lambda x: x.query(f"{timepoint_column} == {timepoint_column}.min()")[args.growth_column].mean(),
                fold_change_growth=lambda x: x["growth_t0"] / x[args.growth_column],
            )
        )
        ref_col, expansion_col = args.growth_column, "fold_change_growth"

    if args.plot:
        fig, axes = grid(ncol=2, aspect_ratio=1.2)

        for ax, y in zip(axes, (ref_col, expansion_col)):
            sc = ax.scatter(
                timepoint_column, y,
                c=args.concentration_column,
                data=expansion_data,
                norm="log",
                cmap="cividis"
            )
            if args.concentration_column is not None and args.concentration_column in expansion_data:
                fig.colorbar(sc, ax=ax, label="Concentration")
            ax.set(
                xlabel="Time",
                ylabel=y,
                yscale="log" if y.startswith("fold_change") else "linear",
            )
        figsave(
            fig, 
            name=f"expansion-{'spike' if args.spike else 'growth'}", 
            df=expansion_data,
        )

    # super rough linear regression for each guide
    expansion_data = expansion_data[[sample_id_column, ref_col, expansion_col]].drop_duplicates()
    print_err(expansion_data.sort_values(sample_id_column))
    regression_input = row_check(pd.merge)(
        data,
        expansion_data,
        how="left",
    )

    regression_input.to_csv(
        "regression-input.tsv",
        sep="\t",
        index=False,
    )

    if args.plot:
        fig, axes = plot_guide_data(
            regression_input,
            x=expansion_col,
            y="fold_change",
            logx=True,
            logy=True,
            panel_group=args.concentration_column,
            replicate_column="replicate",
            group_title="{} µM",
            plot_lines=True,
        )
        figsave(
            fig, 
            name=f"violin_expansion-{'spike' if args.spike else 'growth'}", 
            df=regression_input,
        )

    regression_results, predictions = fit_model(
        regression_input,
        y="fold_change",
        x=expansion_col,
        strain_column=guide_name_column,
        concentration_column=args.concentration_column,
        mode="fitness",
        model="linear_regression",
    )
    regression_results = (
        regression_results
        .rename(columns={
            f"gradient_{expansion_col}": f"relative_fitness_{expansion_col}",
        })
        .assign(
            is_control=lambda x: x[guide_name_column].str.startswith(args.controls),
            is_spike=lambda x: False if args.spike is None else x[guide_name_column].str.startswith(args.spike),
        )
    )

    regression_results.to_csv(
        args.output,
        sep="\t",
        index=False,
    )
    predictions.to_csv(
        "predictions.tsv",
        sep="\t",
        index=False,
    )

    if args.plot:
        fig, axes = grid(aspect_ratio=1.2)
        sc = axes.scatter(
            "predicted_fold_change", "fold_change",
            c=expansion_col,
            data=predictions,
            s=1.,
            cmap="cividis",
            norm="log",
        )
        fig.colorbar(sc, ax=axes, label=expansion_col)
        axes.set(
            xlabel="Modelled",
            ylabel="Observed",
            xscale="log",
            yscale="log",
        )
        axes.plot(axes.get_ylim(), axes.get_ylim(), c="lightgrey", zorder=-1)
        figsave(
            fig, 
            name=f"predictions-{'spike' if args.spike else 'growth'}", 
            df=predictions,
        )
    return None


def _aggregate_concentrations():

    from collections import defaultdict

    from sklearn.metrics import auc

    auc_results = defaultdict(list)
    for guide_name, guide_df in tqdm(  # wrap around any iterator for a progress bar
        regression_results
        .query("condition_iptg_conc > 0.")
        .groupby("guide_name")
    ):
        if guide_df.shape[0] > 1:
            guide_df = guide_df.sort_values("condition_iptg_conc")
            auc_results["guide_name"].append(guide_name)
            this_auc = auc(
                np.log10(guide_df["condition_iptg_conc"]), 
                guide_df["relative_fitness_OD_neg"],
            )
            null_auc = auc(
                np.log10(guide_df["condition_iptg_conc"]), 
                np.ones_like(guide_df["condition_iptg_conc"]), 
            )
            auc_results["auc_relative_fitness_OD_neg"].append(this_auc)
            auc_results["auc_relative_fitness_OD_neg_rel"].append(this_auc / null_auc)
        from scipy.stats import median_abs_deviation

    def zscore(x):
        return (x - x.median()) / median_abs_deviation(x)

    auc_results = (
        pd.DataFrame(auc_results)
        .assign(
            zscore=lambda x: zscore(x["auc_relative_fitness_OD_neg_rel"]),
            z_lt_3=lambda x: x["zscore"] < -3.,
        )
    )
    fig, axes = grid(ncol=2, aspect_ratio=1.1)
    for ax, _x, title in zip(axes, ("auc_relative_fitness_OD_neg_rel", "zscore"), ("AUC relative fitness (normalized)", "AUC Z-score")):
        ax.hist(
            _x,
            data=auc_results,
            bins=80,
        )
        
        ax.set(
            xlabel=title,
            ylabel="Frequency",
        )
    return None


def main():
    inputs = CLIOption(
        'inputs',
        type=str,
        nargs='*',
        help='Count table(s) in TSV or CSV format.',
    )
    conditions = CLIOption(
        '--conditions', '-c', 
        type=FileType('r'),
        required=True,
        help='Table mapping sample IDs to timepoints.',
    )
    reference = CLIOption(
        '--reference', '-r',
        type=str,
        required=True,
        help="Non-targeting guide name(s) to use for reference.",
    )
    count_column = CLIOption(
        '--count-column', '-a',
        type=str,
        default="umi_count",
        help='Column of count table corresponding to counts.',
    )
    guide_column = CLIOption(
        '--guide-column', '-u',
        type=str,
        default="guide_name",
        help='Column of count table corresponding to guide names.',
    )
    growth = CLIOption(
        '--growth', '-g',
        type=FileType('r'),
        help='Table mapping timepoints to a growth measurement.',
    )
    growth_column = CLIOption(
        '--growth-column', '-d',
        type=str,
        default="OD600",
        help='Column of growth table corresponding to growth measurement.',
    )
    concentration_column = CLIOption(
        '--concentration-column', '-k',
        type=str,
        default=None,
        help='Column of conditions table corresponding to concentration.',
    )
    timepoint_column = CLIOption(
        '--timepoint-column', '-t',
        type=str,
        default="timepoint",
        help='Column of conditions table corresponding to timepoint.',
    )
    sample_column = CLIOption(
        '--sample-column', '-s',
        type=str,
        default="sample_id",
        help='Column of conditions table corresponding to sample ID.',
    )
    spike = CLIOption(
        '--spike', '-x',
        type=str,
        default=None,
        help="Instead of growth, use spike-in guide name(s) matching this pattern.",
    )
    controls = CLIOption(
        '--controls', '-z',
        type=str,
        default="controls:",
        help="Control guide name(s) match this pattern.",
    )
    plots = CLIOption(
        '--plot',
        type=str,
        default=None,
        help='If provided, plots will be generated with this filename prefix.',
    )
    plot_format = CLIOption(
        '--plot-format',
        type=str,
        default="png",
        choices=["png", "PNG", "pdf", "PDF"],
        help='File format for plots.',
    )
    output = CLIOption(
        '--output', '-o', 
        type=FileType('w'),
        default=sys.stdout,
        help='Output file. Default: STDOUT',
    )
    formatting = CLIOption(
        '--format', '-f', 
        type=str,
        default='TSV',
        choices=['TSV', 'CSV', 'tsv', 'csv'],
        help='Format of files.',
    )

    fitness_command = CLICommand(
        "fitness",
        description="Calculate fitness from guide counts.",
        options=[
            inputs, 
            conditions,
            reference,
            count_column,
            guide_column,
            sample_column,
            timepoint_column,
            growth,
            growth_column,
            concentration_column,
            spike,
            controls,
            plots,
            output, 
            plot_format,
            formatting,
        ],
        main=_fitness,
    )

    app = CLIApp(
        "fitness", 
        version=__version__,
        description="Calculate fitness from guide counts.",
        commands=[fitness_command],
    )

    app.run()

    return None


if __name__ == "__main__":
    main()
