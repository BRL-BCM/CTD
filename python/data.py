import math
import pandas as pd
import numpy as np


def impute_data(data, ref):

    """
    Impute missing values as lowest observed value in a reference population.

    Parameters
    ----------
    data : Normalized data with some missingness. Data matrix with features as rows, samples as columns.
    ref : Reference sample data with features as rows, samples as columns. Can include some missingness.

    Returns
    -------
    imputed_data : Imputed data.

    Examples
    --------

    """

    # Remove metabolites that are NA for all data/ref samples
    ref.dropna(axis=0, how='all', inplace=True)
    data.dropna(axis=0, how='all', inplace=True)

    # Match metabolites between data and ref matrices
    ref_ind_in_data = [i for i in data.index if i in ref.index]
    data_ind_in_ref = [i for i in ref.index if i in data.index]

    data = data.loc[ref_ind_in_data, :]
    ref = ref.loc[data_ind_in_ref, :]

    imputed_data = data.copy()

    for i, row in ref.iterrows():
        row = row.dropna()
        # Impute using uniform random variable, where a = 0.99*observed minimum, and b = observed minimum
        min_row = min(row)
        data_row = data.loc[[i]]
        cols = data_row.columns[data_row.isna().any()]

        if min_row < 0:
            min_row = -1 * min_row
            i_val = -1
        else:
            i_val = 1

        try:
            imputed_data.loc[i, cols] = i_val * np.random.uniform(low=0.99 * min_row, high=min_row, size=len(cols))
        except Exception:
            pass

    return imputed_data



def surrogate_profiles(data, ref_data, std=1):

    """
    Generate surrogate profiles

    Fill in a data matrix rank with surrogate profiles., when your data is low n, high p.

    Parameters
    ----------
    data : Data matrix with observations (e.g., patient samples) as columns, features (e.g., metabolites or genes) as rows.
    std : The level of variability (standard deviation) around each observed feature's z-score you want to add to
    generate the surrogate profiles.
    ref_data : Data matrix for healthy control "reference" samples, observations (e.g., patient samples) as columns,
    features (e.g., metabolites or genes) as rows.

    Returns
    -------
    data_mx_surr : Data matrix with added surrogate profiles.

    Examples
    --------

    """

    # Match metabolites between data and ref matrices
    ref_ind_in_data = [i for i in data.index if i in ref_data.index]
    data_ind_in_ref = [i for i in ref_data.index if i in data.index]

    data = data.loc[ref_ind_in_data, :]
    ref_data = ref_data.loc[data_ind_in_ref, :]

    rpt = math.ceil(len(data) / len(data.columns) / 2)
    num_surr = math.ceil(len(data) / 2)

    # Generate disease surrogates
    if num_surr > len(data.columns):
        d_surr = pd.DataFrame(np.full((len(data), len(data.columns) + rpt * len(data.columns)), np.NaN))
        c_col = len(data.columns)
        for pt in range(len(data.columns)):
            d_surr.loc[:, pt] = data.iloc[:, pt].values

            for rrpt in range(rpt):
                rr = np.random.normal(loc=0, scale=std, size=len(data))
                d_surr.loc[:, c_col] = data.iloc[:, pt].values + rr
                c_col += 1

        d_surr.rename(columns={i: col for i, col in enumerate(data.columns)}, inplace=True)
        d_surr.rename(columns={(i + len(data.columns)): f'disease_surr{i + 1}' for i in
                               range(len(d_surr.columns) - len(data.columns))}, inplace=True)

        d_surr.index = data.index
    else:
        d_surr = data

    # Generate control surrogates
    if num_surr > len(ref_data.columns):
        rpt = math.ceil(len(ref_data) / len(ref_data.columns) / 2)
        c_surr = pd.DataFrame(np.full((len(ref_data), len(ref_data.columns) + rpt * len(ref_data.columns)), np.NaN))
        c_col = len(ref_data.columns)
        for pt in range(len(ref_data.columns)):
            c_surr.loc[:, pt] = ref_data.iloc[:, pt].values

            for rrpt in range(rpt):
                rr = np.random.normal(loc=0, scale=std, size=len(data))
                c_surr.loc[:, c_col] = ref_data.iloc[:, pt].values + rr
                c_col += 1

        c_surr.rename(columns={i: col for i, col in enumerate(ref_data.columns)}, inplace=True)
        c_surr.rename(columns={(i + len(ref_data.columns)): f'control_surr{i + 1}' for i in
                               range(len(c_surr.columns) - len(ref_data.columns))}, inplace=True)

        c_surr.index = ref_data.index
    else:
        c_surr = ref_data

    # TODO: translate the rest

    return 0

