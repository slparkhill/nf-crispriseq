"""Functions for fitting models."""

from typing import Optional

from collections import defaultdict

from carabiner import print_err
import numpy as np
from numpy.typing import ArrayLike
import pandas as pd
from tqdm.auto import tqdm  # for easy progress bars


def linear_regression(
    df: pd.DataFrame,
    x: str,
    y: str,
    strain_column: str,
    concentration_column: Optional[str] = None,
    death_rate: bool = False
):
    from sklearn.linear_model import LinearRegression

    regression_results = defaultdict(list)

    grouping = [strain_column] + (
        [concentration_column] if concentration_column is not None else []
    ) + (["expt_id"] if "expt_id" in df else [])
    df_grouped = df.groupby(grouping)
    predictions = []
    for (group_name), guide_df in tqdm(
        df_grouped, 
        desc=f"Fitting models within {grouping} groups",
    ):
        guide_df = guide_df.loc[~np.any(guide_df[[x, y]].isna(), axis=1)]
        guide_df = guide_df.loc[np.all(guide_df[[x, y]] > 0, axis=1)]
        if guide_df.shape[0] > 2:  # avoid trivial spike = spike fit
            regressor = LinearRegression(
                fit_intercept=False, 
                # positive=not death_rate,  # fitness is always positive
            )
            for name, val in zip(grouping, group_name):
                regression_results[name].append(val)
            X = np.log(guide_df[[x]]).values
            Y = np.log(guide_df[[y]]).values

            result = regressor.fit(X, X - Y)
            predictions.append(guide_df.assign(**{
                f"predicted_{y}": np.exp(regressor.predict(X))
            }))
            regression_results[f"gradient_{x}"].append(
                result.coef_[0].item()
            )
            regression_results[f"rsq_{x}"].append(
                result.score(X, X - Y)
            )
        else:
            print_err(f"[WARN] Skipping {group_name}!")

    regression_results = pd.DataFrame(regression_results)
    predictions = pd.concat(predictions, axis=0)
    n_groups_start = len(df_grouped)
    n_groups_end = regression_results.shape[0]
    print_err(f"[INFO] {n_groups_start=}, {n_groups_end=}")

    return regression_results, predictions
    

def fit_model(
    df: pd.DataFrame,
    x: str,
    y: str,
    strain_column: str,
    concentration_column: Optional[str] = None,
    mode: str = "fitness",
    model: str = "linear_regression"
):
    if model == "linear_regression":
        model_f = linear_regression
    else:
        raise NotImplementedError(
            f"Model type `{model}` is not implemented."
        )

    model_results = model_f(
        df=df,
        x=x,
        y=y,
        strain_column=strain_column,
        concentration_column=concentration_column,
        death_rate=mode == "killing",
    )

    return model_results


def dose_response_log(
    x: ArrayLike, 
    log_ic50: ArrayLike, 
    log_top: ArrayLike = 1., 
    h: ArrayLike = 1.
) -> np.ndarray:
    return log_top + log_expit(-(np.log(x) - log_ic50) * h)


def fit_dose_response(
    df: pd.DataFrame,
    y: str,
    strain_column: str,
    concentration_column: str,
    plot: Optional[str] = None
):

    from scipy.optimize import curve_fit
    from scipy.special import log_expit

    EPS = np.finfo(float).eps

    param_names = ("log_ic50", "log_top", "h")
    dose_results = defaultdict(list)
    dose_predictions = defaultdict(list)
    for strain_name, strain_df in tqdm(  # wrap around any iterator for a progress bar
        df
        .query(f"{concentration_column} > 0.")
        .groupby(strain_column)
    ):
        dose_results[strain_column].append(strain_name)
        try:
            popt, pcov = curve_fit(
                dose_response_log, 
                guide_df[concentration_column], 
                np.log(np.clip(guide_df[y], EPS, np.inf)),
                maxfev=10000,
            )
        except:
            # make a dummy result
            n_params = len(param_names)
            popt = np.zeros((n_params,))
            pcov = np.inf * np.ones((n_params, n_params))
        for param_name, v, v_var in zip(param_names, popt, np.diag(pcov)):
            dose_results[param_name].append(v.item())
            dose_results[f"{param_name}_var"].append(v_var.item())

        new_rows = {
            strain_column: [strain_name] * strain_df.shape[0],
            concentration_column: strain_df[concentration_column].tolist(),
            y: strain_df[y].tolist(),
            "predicted": np.exp(dose_response_log(
                guide_df[concentration_column], 
                *popt,
            )).tolist()
        }
        for key, val in new_rows.items():
            dose_predictions[key] += val

    dose_results = pd.DataFrame(dose_results)
    dose_results = (
        dose_results
        .assign(
            top=lambda x: np.exp(x["log_top"]),
            ic50=lambda x: np.exp(x["log_ic50"]),
            ic50_z=lambda x: x["log_ic50"] / x["log_ic50_var"],
        )
        .sort_values("ic50_z")
    )
    dose_predictions = pd.DataFrame(dose_predictions)
    if plot:
        fig, ax = grid(aspect_ratio=1.1)
        sc = ax.scatter(
            "predicted", y,
            c=concentration_column,
            data=dose_predictions,
            s=1.,
            norm="log",
        )
        fig.colorbar(sc, ax=ax)
        ax.plot(
            ax.get_ylim(), 
            ax.get_ylim(), 
            color="lightgrey",
        )
        ax.set(
            xlabel=f"Modelled {y}",
            ylabel=f"Observed {y}",
        )

    return None
