# -*- coding: utf-8 -*-

import numpy as np
import xarray as xr
import pandas as pd
from glob import glob
import requests
import os


class GFS_0p25_Hourly():
    def __init__(self, localdir='./', dataurl='http://nomads.ncep.noaa.gov:80/dods/gfs_0p25_1hr'):
        """
        Initialize the GFS_0p25_Hourly class for querying and downloading locally to `localdir`

        :param localdir: Path to the download directory
        :param dataurl: Data URL, defaults to http://nomads.ncep.noaa.gov:80/dods/gfs_0p25_1hr
        """
        self.dataurl = dataurl
        self.localdir = localdir
        self.folders = list()
        self.forecasts = pd.DataFrame()
        self.downloaded = list()

    def query(self):
        self._gfs_available_days()
        self._gfs_available_forecasts()
        self._list_downloaded_cycles()
        print(f'{len(self.downloaded)} available, {len(self.forecasts)} downloadable')

    def download(self):
        pass

    def _gfs_available_days(self):
        response = requests.get(self.dataurl)
        initial = response.text.split('<hr>')[1].replace('<br>', '').replace('<b>', '').replace('</b>', '')
        for line in initial.split('\n')[1:-2]:
            url = line.split('"')[1]
            self.folders.append(os.path.basename(url))

    def _gfs_available_cycles(self, dayid):
        url = f'{self.dataurl}/{dayid}'
        response = requests.get(url)
        initial = response.text.split('<hr>')[1]
        initial = initial.replace('<br>', '').replace('<b>', '').replace('</b>', '')
        initial = initial.replace('\n&nbsp;', '').replace('&nbsp', '').split('\n')
        initial = initial[1:-2]

        eitems = 5  # items expected
        nitems = len(initial) // eitems
        items = np.arange(nitems)

        data = {
            'folder': [dayid for item in items],
            #'ids':[int(initial[eitems*item].strip().replace(':', '')) for item in items],
            'cycle': [initial[eitems * item + 1].split(':')[0] for item in items],
            'inittime': [pd.to_datetime(initial[eitems * item + 1].split('from ')[1].split(',')[0].replace('Z', ''),
                                        format='%H%d%b%Y') for item in items],
            'dltime': [pd.to_datetime(initial[eitems * item + 1].split('Z')[1][5:-4].replace(', downloaded ', ''),
                                      format='%Y%b %d %H:%M') for item in items],
            #'info':[initial[eitems*item+2].split('"')[1] for item in items],
            #'dds':[initial[eitems*item+3].split('"')[1] for item in items],
            #'das':[initial[eitems*item+4].split('"')[1] for item in items]
        }

        data['url'] = [f'{url}/{cycle}' for cycle in data['cycle']]

        data = pd.DataFrame(data)
        data['fid'] = data.apply(lambda x: f'{x.forecast}_{x.cycle}', axis=1)
        data = data.set_index('fid')
        return (data)

    def _list_downloaded_cycles(self):
        fnames = glob(os.path.join(self.localdir, 'gfs*.nc'))
        for fname in fnames:
            fid = os.path.basename(fname).replace('.nc', '')
            if fid not in self.downloaded:
                self.downloaded.append(fid)

    def _gfs_available_forecasts(self):
        for dayid in self.folders:
            self.forecasts = self.forecasts.append(self._gfs_available_cycles(dayid))

    def _get_data(self, folder, cycle, extent=[75, 102, 5, 30]):
        url = f'{self.dataurl}/{folder}/{cycle}'
        outfile = os.path.join(self.localdir, f'{folder}_{cycle}.nc')

        with xr.open_dataset(url) as ds:
            try:
                lon_select = ds['lon'].where(np.logical_and(ds.lon >= extent[0], ds.lon <= extent[1])).dropna(dim='lon')
                lat_select = ds['lat'].where(np.logical_and(ds.lat >= extent[2], ds.lat <= extent[3])).dropna(dim='lat')

                ds_out = xr.Dataset(
                    {
                        'prmsl': ds['prmslmsl'].sel(lat=lat_select, lon=lon_select),
                        'u10': ds['ugrd10m'].sel(lat=lat_select, lon=lon_select),
                        'v10': ds['vgrd10m'].sel(lat=lat_select, lon=lon_select),
                        'stmp': ds['tmp2m'].sel(lat=lat_select, lon=lon_select),
                        'spfh': ds['rh2m'].sel(lat=lat_select, lon=lon_select),
                        'dlwrf': ds['dlwrfsfc'].sel(lat=lat_select, lon=lon_select),
                        'dswrf': ds['dswrfsfc'].sel(lat=lat_select, lon=lon_select),
                        'prate': ds['pratesfc'].sel(lat=lat_select, lon=lon_select),
                    }
                )

                ds_out.to_netcdf(outfile)
            except OSError:
                os.remove(outfile)
            finally:
                ds_out.close()
                ds.close()
