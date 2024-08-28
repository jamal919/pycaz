# -*- coding: utf-8 -*-

import pandas as pd
from pathlib import Path


def read_candhis(fname: str | Path, rename: dict = None) -> pd.DataFrame:
    """
    Read Candhis text data file.

    :param fname: Filename of the Candhis text data.
    :param rename: Dictionary to rename columns.
    :return: DataFrame
    """
    df = pd.read_csv(
        fname,
        sep=';',
        parse_dates=[0],
        na_values=999.999
    )
    df = df.set_index('DateHeure')

    if rename is not None:
        df = df.rename(columns=rename)

    return df
