# -*- coding: utf-8 -*-

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Union, List, Dict

import pandas as pd
import requests
from bs4 import BeautifulSoup
from copy import deepcopy

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

BASIN_CODES = {
    'E': 'ep',
    'C': 'cp',
    'W': 'wp',
    'A': 'io',
    'B': 'io',
    'S': 'sh',
    'P': 'sh',
    'L': 'al'
}


@dataclass
class Basin:
    """
    Class for storing information of the basins
    """
    id: str
    name: str
    code1c: str
    code2c: str
    reldir: str


@dataclass
class Storm:
    """
    Class for storing information of a storm
    """
    basin: Basin
    year: int
    name: str
    cycle: str

    @property
    def id(self):
        return self.name[-3:-1]

    @property
    def code1c(self):
        return self.name[-1]

    @property
    def cycletime(self):
        return self.cycle.split('.')[1]

    @property
    def bstr(self):
        return BASIN_CODES[self.code1c]

    @property
    def model_year(self):
        if (self.bstr == 'sh') & (int(self.cycle.split('.')[1][4:]) >= 70000):
            # year change in southern ocean after July
            return self.year + 1
        else:
            return self.year

    def download(self, baseurl, savedir):
        savedir = Path(savedir)

        # Downloading atcfunix
        atcfunix_path = savedir / '{}.atcfunix'.format(self.cycletime)
        atcfunix_url = '/'.join([
            baseurl,
            self.basin.reldir,
            self.name,
            self.cycle,
            f'{self.cycle.lower()[-14:]}.hfsa.trak.atcfunix'])
        save_to_text(url=atcfunix_url, fname=atcfunix_path)

        # Downloading Deck files
        deck_baseurl = '/'.join([baseurl, 'gc_wmb', 'vxt', 'DECKS'])

        # Deck A file
        deck_a_path = savedir / '{}.adeck.dat'.format(self.cycletime)
        deck_a_url = '/'.join([deck_baseurl, f'a{self.bstr}{self.id}{self.model_year}.dat'])
        save_to_text(url=deck_a_url, fname=deck_a_path)

        # Downloading Deck B file
        deck_b_path = savedir / '{}.bdeck.dat'.format(self.cycletime)
        deck_b_url = '/'.join([deck_baseurl, f'b{self.bstr}{self.id}{self.model_year}.dat'])
        save_to_text(url=deck_b_url, fname=deck_b_path)

        return {
            self.name: {
                'atcfunix': atcfunix_path,
                'deck_a': deck_a_path,
                'deck_b': deck_b_path
            }
        }


def fetch_webpage(
        baseurl: str = 'https://www.emc.ncep.noaa.gov',
        directory: str = 'hurricane/HFSA',
        index: str = 'index.php') -> BeautifulSoup:
    """
    Get the HWRF index webpage and return a soup

    :param baseurl: URL to the NCEP/ENC webpage. Defaults to 'https://www.emc.ncep.noaa.gov'.
    :param directory: Path to HWRF forecasts. Defaults to 'hurricane/HFSA'.
    :param index: Name of the index page. Defaults to 'index.php'.
    :return: BeautifulSoup: A BeautifulSoup representation of the webpage
    """
    url = '/'.join([baseurl, directory, index])
    r = requests.get(url)
    soup = BeautifulSoup(r.text, "html.parser")
    return soup


def has_data(tagdata) -> bool:
    """
    Check if intended variables are available in code

    :param tagdata: tagged data string
    :return:
    """
    data = tagdata.string
    if data is None:
        return False

    if len(re.findall('var basins=*', data)) > 0:
        return True
    else:
        return False


def get_datastring(soup: BeautifulSoup, item=0) -> str:
    """
    Filter and keep only valid scripts tags with code inside

    :param soup: soup created by beautifulsoup4
    :param item: which item to parse
    :return: string of parsable js embedded in js
    """
    valid_data_list = list(filter(has_data, soup.findAll('script')))
    if len(valid_data_list) == 0:
        return ''  # empty string
    elif (len(valid_data_list) > 1) & (item < len(valid_data_list)):
        return valid_data_list[item].string
    else:
        return valid_data_list[0].string


def parse_basins(data: str, item=0) -> List[Basin]:
    """
    Parse basin list from js code as dictionary

    :param data: data as js code scraped from the webpage
    :param item: item index to be parsed as basin
    :return:
    """
    values = re.findall('var basins=(.*?);', data)
    basins_list = json.loads(values[item])

    basins = []

    for i, basin in enumerate(basins_list):
        basin_data = Basin(id=basin['@attributes']['id'],
                           name=basin['name'],
                           code1c=basin['code1c'],
                           code2c=basin['code2c'],
                           reldir=basin['reldir'])
        basins.append(basin_data)

    return basins


def parse_active_storms(data: str, basins: List[Basin]) -> List[Storm]:
    """
    Parse active storm list from js code as dictionary

    :param data: data as js code scraped from the webpage
    :param basins: list of basins parsed using parse_basins
    :return:
    """
    storms = []
    values = re.findall('var actstorm=(.*?);', data)
    actstorms = json.loads(values[0])

    values = re.findall('var actcycle=(.*?);', data)
    actcycles = json.loads(values[0])

    years = re.findall('var year=(.*?);', data)
    year = int(years[0])

    for basinid in actstorms:
        for storm, cycle in zip(actstorms[basinid], actcycles[basinid]):
            storm_data = Storm(
                basin=basins[int(basinid)],
                year=year,
                name=storm,
                cycle=cycle
            )
            storms.append(storm_data)

    return storms


def save_to_text(url: str, fname: Union[str, Path]):
    """Download a textfile from the web

    Args:
        url (str): url to be download from
        fname (Union[str, Path]): filename where the results are to be saved

    Returns:
        bool: if the file is downloaded or not
    """
    fname = Path(fname)
    res = requests.get(url)

    if res.status_code != 200:
        res.raise_for_status()  # Will only raise for 4xx codes
        raise RuntimeError(f"Request to {url} returned status code {res.status_code}")

    fname.parent.mkdir(parents=True, exist_ok=True)

    with open(fname, 'w') as f:
        f.writelines(res.text)


class HAFS:
    def __init__(self, baseurl: str = 'https://www.emc.ncep.noaa.gov', url_kw: Dict = None):
        """
        Download module for Hurricane Aanlysis and Forecasting System forecasts by NOAA

        :param baseurl: The main url of the emc, ie., https://www.emc.ncep.noaa.gov
        :param url_kw: A dictionary with two keys, `directory` and `index`.
            Default value: {'directory': 'hurricane/HFSA', 'index': 'index.php'}.
            Other known option for `directory` is 'hurricane/HFSB'
        """
        self.baseurl: str = baseurl
        if url_kw is None:
            url_kw = dict(directory='hurricane/HFSA', index='index.php')
        self.url_kw = url_kw

        self.soup = None
        self.datastring = None
        self.last_updated = None
        self.basins = None
        self.storms = None

        self.update()

    @property
    def active_storms(self):
        return self.storms

    def copy(self):
        """
        Create a deepcopy of self

        :return: deepcopy of self
        """
        return deepcopy(self)

    def update(self):
        """
        Fetch website and update the basins and storms by scrapting it.

        :return: None
        """
        try:
            self.soup = fetch_webpage(baseurl=self.baseurl, **self.url_kw)
        except Exception as e:
            print('Webpage could not be fetched, following exception was thrown: {}'.format(e))

        self.datastring = get_datastring(self.soup)
        self.last_updated = pd.to_datetime('now')
        self.basins = parse_basins(data=self.datastring, item=0)
        self.storms = parse_active_storms(data=self.datastring, basins=self.basins)

    def filter(self, basin: Literal['AL', 'EP', 'CP', 'WP', 'IO', 'SH']):
        """
        Filter by a single basin using 2c codes

        :param basin: Basin code, one of 'AL', 'EP', 'CP', 'WP', 'IO', 'SH'
        :return: Filtered object
        """
        filtered_storms = list(filter(lambda x: x.basin.code2c == basin, self.storms))
        obj = self.copy()
        obj.storms = filtered_storms
        return obj

    def __repr__(self):
        nstorm = len(self.storms)
        if nstorm == 0:
            return 'No Active storms'
        else:
            return f'{nstorm} active storms! \nUse `storms` property to see the list, or download using .download()'

    def download(self, savedir: Union[str, Path], cycle_dir: bool = True, basin_dir: bool = True,
                 storm_dir: bool = True) -> Dict:
        """
        Download the current list of storms

        :param savedir: Path to the save directory
        :param cycle_dir: If the named cycle directory to be created
        :param basin_dir: If the named basin directory to be created
        :param storm_dir: If the named storm directory to be created
        :return:
        """
        savedir = Path(savedir)
        if not savedir.exists():
            savedir.mkdir()

        all_downloaded_files = {}

        for storm in self.storms:
            _savedir = savedir
            if cycle_dir:
                _savedir = _savedir / storm.cycletime
                if not _savedir.exists():
                    _savedir.mkdir()

            if basin_dir:
                _savedir = _savedir / storm.basin.code2c
                if not _savedir.exists():
                    _savedir.mkdir()

            if storm_dir:
                _savedir = _savedir / storm.name
                if not _savedir.exists():
                    _savedir.mkdir()

            downloaded_files = storm.download(baseurl=self.baseurl, savedir=_savedir)
            all_downloaded_files.update(downloaded_files)
        return all_downloaded_files
