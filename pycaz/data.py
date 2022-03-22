import numpy as np
import xarray as xr
import pandas as pd
from glob import glob
import requests
import os

class GFS_0p25_hourly():
    def __init__(self, localdir='./', dataurl='http://nomads.ncep.noaa.gov:80/dods/gfs_0p25_1hr'):
        self.dataurl = dataurl
        self.localdir = localdir
        self.folders = []
        self.forecasts = pd.DataFrame()
        self.downloaded = []
    
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

        eitems = 5 # items expected 
        nitems = len(initial)//eitems
        items = np.arange(nitems)

        data = {
            'folder':[dayid for item in items],
            #'ids':[int(initial[eitems*item].strip().replace(':', '')) for item in items],
            'cycle':[initial[eitems*item+1].split(':')[0] for item in items],
            'inittime':[pd.to_datetime(initial[eitems*item+1].split('from ')[1].split(',')[0].replace('Z', ''), format='%H%d%b%Y') for item in items],
            'dltime':[pd.to_datetime(initial[eitems*item+1].split('Z')[1][5:-4].replace(', downloaded ', ''), format='%Y%b %d %H:%M') for item in items],
            #'info':[initial[eitems*item+2].split('"')[1] for item in items],
            #'dds':[initial[eitems*item+3].split('"')[1] for item in items],
            #'das':[initial[eitems*item+4].split('"')[1] for item in items]
        }

        data['url'] = [f'{url}/{cycle}' for cycle in data['cycle']]

        data = pd.DataFrame(data)
        data['fid'] = data.apply(lambda x: f'{x.forecast}_{x.cycle}', axis=1)
        data = data.set_index('fid')
        return(data)
    
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
                lon_select = ds['lon'].where(np.logical_and(ds.lon>=extent[0], ds.lon<=extent[1])).dropna(dim='lon')
                lat_select = ds['lat'].where(np.logical_and(ds.lat>=extent[2], ds.lat<=extent[3])).dropna(dim='lat')

                ds_out = xr.Dataset(
                    {
                        'prmsl':ds['prmslmsl'].sel(lat=lat_select, lon=lon_select),
                        'u10':ds['ugrd10m'].sel(lat=lat_select, lon=lon_select),
                        'v10':ds['vgrd10m'].sel(lat=lat_select, lon=lon_select),
                        'stmp':ds['tmp2m'].sel(lat=lat_select, lon=lon_select),
                        'spfh':ds['rh2m'].sel(lat=lat_select, lon=lon_select),
                        'dlwrf':ds['dlwrfsfc'].sel(lat=lat_select, lon=lon_select),
                        'dswrf':ds['dswrfsfc'].sel(lat=lat_select, lon=lon_select),
                        'prate':ds['pratesfc'].sel(lat=lat_select, lon=lon_select),
                    }
                )

                ds_out.to_netcdf(outfile)
            except OSError:
                os.remove(outfile)
            finally:
                ds_out.close()
                ds.close()

import pandas as pd
import requests
from bs4 import BeautifulSoup
import re

class NOAA_HWRF:
    def __init__(self, localdir='./hwrf', year=2022, basin='NIO'):
        self.weburl = 'https://www.emc.ncep.noaa.gov/gc_wmb/vxt/HWRF/index.php'
        self.tcallurl = 'https://www.emc.ncep.noaa.gov/gc_wmb/vxt/HWRF/tcall.php'
        self.forecasturl = f'https://www.emc.ncep.noaa.gov/gc_wmb/vxt/HWRFForecast/RT{year}_NIO/'
        self.decks_url = 'https://www.emc.ncep.noaa.gov/gc_wmb/vxt/DECKS'
        self.localdir = localdir
        self.last_checked = pd.NaT
        self.basin = basin
        self.nactive_storm = 0
        self.active_storms = {}
        self.available_basin = {
            'NIO':'activeNorth Indian Ocean'
        }
        
        if not os.path.exists(self.localdir):
            os.mkdir(self.localdir)
            
    def query(self):
        self._query_active_storms()
        self._query_hwrf_cycles()
        self.last_checked = pd.to_datetime('now')
        return(self)
        
    def download(self):
        self._download_hwrf()
        self._download_besttrack()
        self._download_all_models()
            
    def _query_active_storms(self):
        r = requests.get(self.weburl)
        soup = BeautifulSoup(r.text, "html.parser")

        # if there is an active storm it will appear as a selection form
        storms = soup.find('td', {'name':'activeNorth Indian Ocean'}).find_all('form')

        for storm in storms:
            # figure out how to call them through tcall.php
            params = {}
            for input_tag in storm.find_all("input"):
                input_type = input_tag.attrs.get("type", "text")
                input_name = input_tag.attrs.get("name")
                input_value =input_tag.attrs.get("value", "")

                if input_type == 'submit':
                    storm_name = input_value
                if input_type == 'hidden':
                    params[input_name] = input_value

            self.active_storms[storm_name] = {'params':params}
        
        return(self)
        
    def _query_hwrf_cycles(self):
        for storm in self.active_storms:
            req_tcall = requests.post(self.tcallurl, data=self.active_storms[storm]['params'])
            soup_tcall = BeautifulSoup(req_tcall.text, 'html.parser')
            cycles = soup_tcall.find('select', {'name':'selectCycle'}).find_all('option')
            cycles = [cycle.text for cycle in cycles]
            self.active_storms[storm]['cycles'] = cycles
            print(f'{storm}: {len(cycles)} HWRF forecast cycles available')
            
        return(self)
    
    def _download_hwrf(self):
        for storm in self.active_storms:
            for cycle in self.active_storms[storm]['cycles']:
                stormfilename = f'{cycle.lower()}.trak.hwrf.atcfunix'
                stormfile = os.path.join(self.localdir, stormfilename)
                if os.path.exists(stormfile):
                    print(f'\t{cycle} : {stormfilename} exists!')
                else:
                    hwrf_atcf_url = f'{self.forecasturl}/{self.active_storms[storm]["params"]["selectStorm"]}/{cycle}/{stormfilename}'
                    req_atcf = requests.get(hwrf_atcf_url)
                    req_atcf.text

                    with open(stormfile, 'w') as f:
                        f.writelines(req_atcf.text)

                    print(f'\t{cycle} : {stormfilename} downloaded!')
                
    def _download_besttrack(self):
        prefix = 'b' # b for jtwc best track
        for storm in self.active_storms:
            lastcycle, lasttime = self.active_storms[storm]['cycles'][0].split('.')
            rs = re.search(r"\d+", lastcycle) # finds the range for the numbered part
            trackid = lastcycle[rs.start():rs.end()]
            lasttime = pd.to_datetime(lasttime, format='%Y%m%d%H')

            storm_name = f'{prefix}io{trackid}{lasttime.year}.dat'
            decks_storm = f'{self.decks_url}/{storm_name}'
            print(f'Downloading file - {storm_name}')

            decks_req = requests.get(decks_storm)
            with open(os.path.join(self.localdir, f'{storm_name}'), 'w') as f:
                f.writelines(decks_req.text)
    
    def _download_all_models(self):
        prefix = 'a' # a for analysis/model
        for storm in self.active_storms:
            lastcycle, lasttime = self.active_storms[storm]['cycles'][0].split('.')
            rs = re.search(r"\d+", lastcycle) # finds the range for the numbered part
            trackid = lastcycle[rs.start():rs.end()]
            lasttime = pd.to_datetime(lasttime, format='%Y%m%d%H')

            storm_name = f'{prefix}io{trackid}{lasttime.year}.dat'
            decks_storm = f'{self.decks_url}/{storm_name}'
            print(f'Downloading file - {storm_name}')

            decks_req = requests.get(decks_storm)
            with open(os.path.join(self.localdir, f'{storm_name}'), 'w') as f:
                f.writelines(decks_req.text)