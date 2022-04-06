#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import requests
from bs4 import BeautifulSoup
import re
import os

class NOAA_HWRF:
    def __init__(self, localdir:str='./hwrf', year:int=2022, basin:str='NIO'):
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
        storms = soup.find('td', {'name':self.available_basin[self.basin]}).find_all('form')

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