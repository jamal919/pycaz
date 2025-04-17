# -*- coding: utf-8 -*-

import pandas as pd


def flat_date_parser(date_string):
    """
    Parsing MATLAB compatible flat date lists. In pd.read_csv(), the parse_dates should be a nested
    list of the columns containing [[y, m, d, H, M. S]]

    :param date_string: the datestring to be parsed
    :return: pd.Datetime
    """
    y, m, d, H, M, S = date_string.split()
    y, m, d, H, M, S = int(y), int(m), int(d), int(H), int(M), int(S)
    date_formatted = f'{y:4d}-{m:02d}-{d:02d} {H:02d}:{M:02d}:{S:02d}'
    return pd.to_datetime(date_formatted)
