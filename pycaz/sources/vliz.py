# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import requests
from tqdm import tqdm

def list_stations() -> pd.DataFrame:
    """
    List all available stations at ioc-sealevelmonitoring.org

    :return: DataFrame
    """
    stn_url = "https://www.ioc-sealevelmonitoring.org/list.php?operator=&showall=all&output=general"
    parsed_items = pd.read_html(stn_url)

    # the largest table contains the stations
    istn_table = np.argmax([len(item) for item in parsed_items])
    stn_df = parsed_items[istn_table]

    # it gives two name columns that needs cleaning
    stn_df = stn_df.drop([0]) # unnecessary column

    # this is the good columnnames
    colnames = stn_df.iloc[0, :].values
    icols = np.logical_not(stn_df.iloc[0, :].isna().values)
    colnames = colnames[icols]
    stn_df = stn_df.drop([1])
    stn_df = stn_df.iloc[:, icols]
    stn_df.columns = colnames

    return stn_df

def get_data(stnid) -> pd.DataFrame:
    """
    Scrap the data for stnid

    :param stnid: stnid from ioc-sealevelmonitoring.org, use list_stations()
    :return: pd.DataFrame of scraped data
    """
    data = pd.DataFrame({})

    date_current = pd.Timestamp('now').round('D')
    date_until = pd.Timestamp('2000-01-01')
    dt = pd.Timedelta('30d')

    chunks = np.ceil((date_current - date_until)/dt).astype(int)
    dates_to_process = pd.DataFrame({
        'enddate':pd.to_datetime(date_current - np.arange(0, chunks) * dt)
        })
    dates_to_process.loc[:, 'startdate'] = dates_to_process.loc[:, 'enddate'] - dt
    dates_to_process.loc[:, 'available'] = False
    dates_to_process.loc[:, 'count'] = np.nan

    for i, date in enumerate(tqdm(dates_to_process.loc[:, 'enddate'])):
        endtime = date.strftime("%Y-%m-%d")
        url = f'https://www.ioc-sealevelmonitoring.org/bgraph.php?code={stnid}&output=tab&period=30&endtime={endtime}'
        response = requests.get(url)
        if "NO DATA" in response.text:
            pass
        else:
            dates_to_process.loc[i, 'available'] = True
            df = pd.read_html(response.text, header=0)[0]
            df.loc[:, 'Time (UTC)'] = pd.to_datetime(df.loc[:, 'Time (UTC)'])
            dates_to_process.loc[i, 'count'] = len(df)
            data = pd.concat([data, df])

    data_sorted = data.loc[np.logical_not(data.loc[:, 'Time (UTC)'].duplicated()), :]
    data_sorted = data_sorted.set_index('Time (UTC)').sort_index()

    return data_sorted