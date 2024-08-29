#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from datetime import datetime
from pycaz.cyclone.track import Record, Track
from pycaz.convert import knot2mps, hpa2pa, ntm2m, ft2m
from logging import getLogger

logger = getLogger(__name__)


def read_jtwc(fname, replace_zero_radial=True):
    '''
    A reader function for reading JTWC file.

    The coded fields not found on the tracks are tagged as np.nan
    '''
    # Wind codes not implemented
    WINDCODE_NOTIMPLEMENTED = {
        'NNS': 'north semicircle',
        'NES': 'northeast semicircle',
        'EES': 'east semicircle',
        'SES': 'southeast semicircle',
        'SSS': 'south semicircle',
        'SWS': 'southwest semicircle',
        'WWS': 'west semicircle',
        'NWS': 'northwest semicirlce',
        'NNQ': 'north north',
        'EEQ': 'east east',
        'SEQ': 'south east',
        'SSQ': 'south south',
        'SWQ': 'south west',
        'WWQ': 'west west',
        'NWQ': 'north west'}

    # Possible wind code
    WINDCODE_IMPLEMENTED = {
        'AAA': 'full circle',
        'NEQ': 'north east'
    }

    # Dummy variable for storage
    track = {}

    with open(fname) as f:
        records = f.readlines()
        logger.info(f'Reading {fname}')

    for linum, record in enumerate(records):
        logger.info(f'Processing record {linum + 1}')
        fields = record.split(',')

        # BASIN - basin, e.g. WP, , SH, CP, EP, AL
        basin = fields[0].strip()

        # CY - annual cyclone number: 1 through 99
        ncyclone = int(fields[1].strip())

        # YYYYMMDDHH - Warning Date-Time-Group: 0000010100 through 9999123123
        # Previous cyclone may have 2 digit year
        timestamp = fields[2].strip()

        # TECHNUM - objective technique sorting number: 00 - 99
        technum = fields[3].strip()

        # TECH - acronym for each objective technique or CARQ or WRNG
        # BEST for best track.
        tech = fields[4].strip()

        # TAU - forecast period: -24 through 120 hours
        # 0 for best-track
        # negative taus used for CARQ and WRNG records.
        tau = int(fields[5].strip())

        # LatN/S - Latitude (tenths of degrees) for the DTG
        # 0 through 900, N/S is the hemispheric index.
        lat = fields[6].strip()
        if lat[-1] == 'N':
            lat = float(lat[0:-1]) / 10
        elif lat[-1] == 'S':
            lat = float(lat[0:-1]) / 10 * -1

        # LonE/W - Longitude (tenths of degrees) for the DTG
        # 0 through 1800, E/W is the hemispheric index.
        lon = fields[7].strip()
        if lon[-1] == 'E':
            lon = float(lon[0:-1]) / 10
        elif lon[-1] == 'W':
            lon = float(lon[0:-1]) / 10 * -1
        else:
            raise Exception(f'Unknown direction code {lon[-1]}!')

        # VMAX - Maximum sustained wind speed in knots: 0 through 300
        try:
            vmax = float(fields[8].strip())
        except:
            vmax = np.nan
        finally:
            vmax = knot2mps(vmax)

        # MSLP - Minimum sea level pressure, 1 through 1100 MB. (hpa)
        try:
            mslp = float(fields[9].strip())
        except:
            mslp = np.nan
        finally:
            mslp = hpa2pa(mslp)

        # TY - Level of tc development:
        #    TD - tropical depression,
        #    TS - tropical storm,
        #    TY - typhoon,
        #    ST - super typhoon,
        #    TC - tropical cyclone,
        #    HU - hurricane,
        #    SD - subtropical depression,
        #    SS - subtropical storm,
        #    EX - extratropical systems,
        #    MD - monsoon depression,
        #    IN - inland,
        #    DS - dissipating,
        #    LO - low,
        #    WV - tropical wave,
        #    ET - extrapolated,
        #    XX - unknown.
        try:
            ty = fields[10].strip()
        except:
            ty = 'XX'

        # RAD - Wind intensity (kts) for the radii defined in this record
        # Typically 35, 50, 65 or 100.
        # Also 34, 50, 64
        try:
            radv = float(fields[11].strip())
            if radv == 0:
                radv = np.nan
        except:
            radv = np.nan
        finally:
            radv = knot2mps(radv)

        # WINDCODE - Radius code
        #    AAA - full circle
        #    NNS - north semicircle
        #    NES - northeast semicircle
        #    EES - east semicircle
        #    SES - southeast semicircle
        #    SSS - south semicircle
        #    SWS - southwest semicircle
        #    WWS - west semicircle
        #    NWS - northwest semicirlce
        #    QQQ - quadrant (NNQ, NEQ, EEQ, SEQ, SSQ, SWQ, WWQ, NWQ)
        try:
            windcode = fields[12].strip()
        except:
            windcode = 'AAA'

        if windcode in WINDCODE_NOTIMPLEMENTED:
            raise Exception(f'L{linum + 1} - Windcode {windcode} not implemented at line {linum + 1}')

        # RAD1 - If full circle, radius of specified wind intensity
        # If semicircle or quadrant, radius of specified wind intensity of circle portion specified in radius code. 
        # 0 - 1200 nm.
        try:
            rad1 = float(fields[13].strip())
            if rad1 == 0:
                rad1 = np.nan
        except:
            rad1 = np.nan
        finally:
            rad1 = ntm2m(rad1)

        # RAD2 - If full circle this field not used
        # If semicicle, radius (nm) of specified wind intensity for semicircle not specified in radius code, 
        # If quadrant, radius (nm) of specified wind intensity for 2nd quadrant 
        # counting clockwise from quadrant specified in radius code
        # 0 through 1200 nm.
        try:
            rad2 = float(fields[14].strip())
            if rad2 == 0:
                rad2 = np.nan
        except:
            rad2 = np.nan
        finally:
            rad2 = ntm2m(rad2)

        # RAD3 - If full circle or semicircle this field not used
        # If quadrant, radius (nm) of specified wind intensity for 3rd quadrant
        # counting clockwise from quadrant specified in radius code
        # 0 through 1200 nm.
        try:
            rad3 = float(fields[15].strip())
            if rad3 == 0:
                rad3 = np.nan
        except:
            rad3 = np.nan
        finally:
            rad3 = ntm2m(rad3)

        # RAD4 - If full circle or semicircle this field not used
        # If quadrant, radius (nm) of specified wind intensity for 4th quadrant
        # counting clockwise from quadrant specified in radius code
        # 0 through 1200 nm.
        try:
            rad4 = float(fields[16].strip())
            if rad4 == 0:
                rad4 = np.nan
        except:
            rad4 = np.nan
        finally:
            rad4 = ntm2m(rad4)

        # Setting all the radius value to same for full circle
        if windcode == 'AAA':
            # if full circle
            rad2 = rad1
            rad3 = rad1
            rad4 = rad1

        # RADP - pressure in millibars of the last closed isobar, 900 - 1050 mb (hpa)
        try:
            pres = float(fields[17].strip())
        except:
            pres = np.nan
        finally:
            pres = hpa2pa(pres)

        # RRP - radius of the last closed isobar in nm, 0 - 9999 nm.
        try:
            rpres = float(fields[18].strip())
        except:
            rpres = np.nan
        finally:
            rpres = ntm2m(rpres)

        # MRD - radius of max winds, 0 - 999 nm.
        try:
            rmax = float(fields[19].strip())
        except:
            rmax = np.nan
        finally:
            rmax = ntm2m(rmax)

        # GUSTS - gusts, 0 through 995 kts.
        try:
            gusts = float(fields[20].strip())
        except:
            gusts = np.nan
        finally:
            gusts = knot2mps(gusts)

        # EYE - eye diameter, 0 through 999 nm.
        try:
            eye = float(fields[21].strip())
        except:
            eye = np.nan
        finally:
            eye = ntm2m(eye)

        # SUBREGN - subregion code: W, A, B, S, P, C, E, L.
        #    A - Arabian Sea
        #    B - Bay of Bengal
        #    C - Central Pacific
        #    E - Eastern Pacific
        #    L - Atlantic
        #    P - South Pacific (135E - 120W)
        #    S - South (20E - 135E)
        #    W - Western Pacific
        try:
            subregion = fields[22]
        except:
            subregion = None

        # MAXSEAS - max seas: 0 through 999 ft.
        try:
            maxseas = float(fields[23].strip())
        except:
            maxseas = np.nan
        finally:
            maxseas = ft2m(maxseas)

        # INITIALS - Forecaster's initials, used for tau 0 WRNG, up to 3 chars.
        try:
            initials = fields[24]
        except:
            initials = None

        # DIR - storm direction, 0 - 359 degrees.
        try:
            stormdir = float(fields[25].strip())
        except:
            stormdir = np.nan

        # SPEED - storm speed, 0 - 999 kts.
        try:
            stormspeed = float(fields[26].strip())
        except:
            stormspeed = np.nan
        finally:
            stormspeed = knot2mps(stormspeed)

        # STORMNAME - literal storm name, NONAME or INVEST
        # TCcyx used pre-1999
        #    where: cy = Annual cyclone number 01 through 99
        #            x = Subregion code: W, A, B, S, P, C, E, L.
        #                A - Arabian Sea
        #                B - Bay of Bengal
        #                C - Central Pacific
        #                E - Eastern Pacific
        #                L - Atlantic
        #                P - South Pacific (135E - 120W)
        #                S - South (20E - 135E)
        #                W - Western Pacific
        try:
            stormname = fields[27]
        except:
            stormname = None

        # DEPTH - system depth, D-deep, M-medium, S-shallow, X-unknown
        try:
            depth = fields[28]
        except:
            depth = None

        # SEAS - Wave height for radii defined in SEAS1-SEAS4, 0-99 ft.
        try:
            seas = float(fields[29])
        except:
            seas = np.nan
        finally:
            seas = ft2m(seas)

        # SEASCODE - Radius code:
        #    AAA - full circle
        #    NNS - north semicircle
        #    NES - northeast semicircle
        #    EES - east semicircle
        #    SES - southeast semicircle
        #    SSS - south semicircle
        #    SWS - southwest semicircle
        #    WWS - west semicircle
        #    NWS - northwest semicircle
        #    QQQ - quadrant (NNQ, NEQ, EEQ, SEQ, SSQ, SWQ, WWQ, NWQ)
        try:
            seacode = fields[30]
        except:
            seacode = 'AAA'

        # SEAS1 - first quadrant seas radius as defined by SEASCODE, 0 through 999 nm.
        try:
            seas1 = int(fields[31])
        except:
            seas1 = np.nan
        finally:
            seas1 = ntm2m(seas1)

        # SEAS2 - second quadrant seas radius as defined by SEASCODE, 0 through 999 nm.
        try:
            seas2 = int(fields[32])
        except:
            seas2 = np.nan
        finally:
            seas2 = ntm2m(seas2)

        # SEAS3 - third quadrant seas radius as defined by SEASCODE, 0 through 999 nm.
        try:
            seas3 = int(fields[33])
        except:
            seas3 = np.nan
        finally:
            seas3 = ntm2m(seas3)

        # SEAS4 - fourth quadrant seas radius as defined by SEASCODE, 0 through 999 nm.
        try:
            seas4 = int(fields[34])
        except:
            seas4 = np.nan
        finally:
            seas4 = ntm2m(seas4)

        if seacode == 'AAA':
            # Full circle
            seas2 = seas1
            seas3 = seas1
            seas4 = seas1

        # Creating fields for building tracks
        vinfo = np.array([radv])
        radinfo = np.array([rad1, rad2, rad3, rad4])
        if replace_zero_radial:
            # replace the zero values in radinfo with the maximum distance
            radinfo_old = radinfo.copy()
            radinfo[np.isnan(radinfo)] = np.nanmax(radinfo)
            logger.info(f'radinfo is updated from {radinfo_old} to {radinfo} for line {linum + 1}')

        seas = np.array([seas])
        seainfo = np.array([seas1, seas2, seas3, seas4])
        if replace_zero_radial:
            # replace the zero values in radinfo with maximum distance
            seainfo_old = seainfo.copy()
            seainfo[np.isnan(seainfo)] = np.nanmax(seainfo)
            logger.info(f'seainfo is updated from {seainfo_old} to {seainfo} for line {linum + 1}')

        if timestamp not in track:
            # new record
            track[timestamp] = {
                'basin': basin,
                'ncyclone': ncyclone,
                'technum': technum,
                'tech': tech,
                'tau': tau,
                'lon': lon,
                'lat': lat,
                'vmax': vmax,
                'ty': ty,
                'mslp': mslp,
                'pres': pres,
                'rpres': rpres,
                'gusts': gusts,
                'rmax': rmax,
                'eye': eye,
                'subregion': subregion,
                'maxseas': maxseas,
                'initials': initials,
                'stormdir': stormdir,
                'stormspeed': stormspeed,
                'stormname': stormname,
                'depth': depth,
                'seas': seas,
                'windcode': windcode,
                'vinfo': vinfo,
                'radinfo': radinfo,
                'seacode': seacode,
                'seainfo': seainfo}
        else:
            if radv not in np.array(track[timestamp]['vinfo']):
                track[timestamp]['vinfo'] = np.append(track[timestamp]['vinfo'], vinfo)
                track[timestamp]['radinfo'] = np.vstack((track[timestamp]['radinfo'], radinfo))

            if seas not in np.array(track[timestamp]['seas']):
                track[timestamp]['seas'] = np.append(track[timestamp]['seas'], seas)
                track[timestamp]['seainfo'] = np.vstack((track[timestamp]['seainfo'], seainfo))

    # Create the track class
    records = np.array([])
    for record in track:
        record_obj = Record(
            timestamp=datetime.strptime(record, '%Y%m%d%H'),
            lon=track[record]['lon'],
            lat=track[record]['lat'],
            mslp=track[record]['mslp'],
            vmax=track[record]['vmax'],
            rmax=track[record]['rmax'],
            vinfo=track[record]['vinfo'],
            radinfo=track[record]['radinfo']
        )
        records = np.append(records, record_obj)

    return (Track(records=records))
