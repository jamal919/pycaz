# -*- coding: utf-8 -*-

import numpy as np
from scipy import optimize

from pycaz.convert import pa2mb, km2m
from pycaz.typing import ArrayLike


def coriolis(lat: ArrayLike) -> ArrayLike:
    """
    Calculate the coriolis parameter with latitude (lat)

    :param lat: Latitude from -90 to 90
    :return: Coriolis parameter
    """
    f = 2 * 7.292e-5 * np.sin(np.deg2rad(lat))
    return f


def calc_rmax_s02(mslp: ArrayLike) -> ArrayLike:
    """
    Calculates rmax as a function of central pressure as described in Silva 2002.

    Ref: Silva,  R., G. Georges, S. Paulo, B. Gustavo and B. G. D. Gabriel (2002).
    Oceanographic vulnerability to hurricanes on the Mexican coast.
    Pro-ceedings of 28th Conference on Coastal Engineering, 39-51.

    :param mslp: central pressure in Pascal
    :return: rmax in m
    """
    try:
        mslp = pa2mb(mslp)
    except:
        raise Exception(f'mslp is a required input for S02 method of rmax calculation')

    rmax = 0.4785 * mslp - 413  # In Km
    rmax = km2m(rmax)

    return rmax


def calc_rmax_w04(vmax: ArrayLike, lat: ArrayLike) -> ArrayLike:
    """
    Calculates rmax as a function of the maximum velocity as described in
    Willoughby and Rahn (2004) paper.

    vmax: float, m/s
    lat: float, m/s

    Ref: Willoughby and Rahn (2004) Parametric Representation of the Primary Hurricane
    Vortex. Part I: Observations andEvaluation of the Holland (1980) Model, Monthly
    Wearher Review


    :param vmax: float, m/s
    :param lat: latitude from -90 to 90
    :return: rmax in m
    """
    try:
        rmax = 51.6 * np.exp(-0.0223 * vmax + 0.0281 * lat)  # Km
    except:
        raise Exception(f'vmax and lat are required input for W04 method of rmax calculation')

    rmax = km2m(rmax)

    return rmax


def calc_rmax_e11(v, r, vmax, f, solver='fsolve', limit=(5000, 500000), step=100):
    """
    Calculate maximum radius using Emanuel 2011 model.

    Ref: Emanuel 2011, Lin and Chavas 2012, Krien et al. 2017

    :param v: float, known velocity information at the boundary layer, m/s
    :param r: float, known radial information corresponding to v, m
    :param vmax: float, known maximum velocity, m/s
    :param f: coriolis coefficient, taken at the center as it is small
    :param solver: solver type. 'fsolve': use scipy.optimize.fsolve; 'scan': use linear scanning
    :param limit: limit for scanning in scan solver and limit[0] is x0 for fsolve
    :param step: stepping for scan solver
    :return: rmax in m
    """

    rmax_solved = np.nan
    resfunc = lambda rmax: v - (2 * r * (vmax * rmax + 0.5 * f * rmax ** 2) / (rmax ** 2 + r ** 2) - f * r / 2)
    try:
        if solver == 'scan':
            rm_range = np.arange(start=limit[0], stop=limit[1] + step, step=step)
            for rm in rm_range:
                res = resfunc(rm)
                if res < 0:
                    break
            rmax_solved = rm
        elif solver == 'fsolve':
            rmax_solved = optimize.fsolve(func=resfunc, x0=limit[0])
    except:
        raise Exception('Solver failed')

    return rmax_solved


def calc_holland_B(vmax, pc, pn, rhoair, bmax=2.5, bmin=0.5):
    """
    Calculates the simplified version of Holland's B parameters excluding the term with coriolis.

    :param vmax: Velocity of maximum wind at boundary layer, m/s
    :param pc: Central pressure, Pa
    :param pn: environmental pressure, Pa
    :param rhoair: Density of air kg/m**3
    :param bmax: maximum limit of B
    :param bmin: minimum limit of B
    :return: Holland B parameter
    """
    B = (vmax ** 2) * rhoair * np.exp(1) / (pn - pc)

    if isinstance(B, float):
        if B < bmin:
            B = bmin
        elif B > bmax:
            B = bmax
        else:
            B = B
    else:
        B[B < bmin] = bmin
        B[B > bmax] = bmax

    return B


def calc_holland_B_full(vmax, rmax, pc, f, pn, rhoair, bmax=2.5, bmin=0.5):
    """
    Calculates the holland B paramter using the full expression.

    :param vmax: Velocity of maximum wind at boundary layer, m/s
    :param rmax: radius of maximum wind, m
    :param pc: central pressure, Pa
    :param f: coriolis
    :param pn: Environmental pressure, e.g., 101325pa
    :param rhoair: density of air kg/m**3
    :param bmax: maximum limit of B
    :param bmin: minimum limit of B
    :return: Holland B parameter
    """
    B = (vmax ** 2 * rhoair * np.exp(1) + f * vmax * rmax * np.exp(1) * rhoair) / (pn - pc)

    if isinstance(B, float):
        if B < bmin:
            B = bmin
        elif B > bmax:
            B = bmax
        else:
            B = B
    else:
        B[B < bmin] = bmin
        B[B > bmax] = bmax

    return B


def calc_rmax_h80(v, r, pc, B, f, pn, rhoair, solver='fsolve', limit=(5000, 100000), step=100):
    """
    Solve rmax with for a given r and v using Holland 1980 model.

    :param v: float, velocity in m/s
    :param r: float, radial distance corresponding to v, in m
    :param pc: float, central pressure, Pa
    :param B: float, Holland B parameter, calculable by calc_holland_B()
    :param f: float, coriolis parameter
    :param pn: Environmental pressure, e.g., 101325pa
    :param rhoair: density of air kg/m**3
    :param solver: solver to use - 'scan' or 'fsolve' (default)
    :param limit: limit for scan solver or starting point for fsolve
    :param step: scan solver stepping, in m
    :return: rmax in m
    """
    rmax_solved = np.nan
    resfunc = lambda rmax: v - (np.sqrt(
        B / rhoair * (rmax / r) ** B * (pn - pc) * np.exp(-(rmax / r) ** B) + (r * f / 2) ** 2) - r * f / 2)

    try:
        if solver == 'scan':
            rm_range = np.arange(start=limit[0], stop=limit[1] + step, step=step)
            for rm in rm_range:
                res = resfunc(rm)
                if res < 0:
                    break
            rmax_solved = rm
        elif solver == 'fsolve':
            rmax_solved = optimize.fsolve(func=resfunc, x0=limit[0])
    except:
        raise Exception('Solver failed')

    return rmax_solved


def calc_vcirc_j92(r, Rm, Vm):
    """
    Calculate circular wind speed according to Jelesnianski et al. 1992 model.
    It is also known as SLOSH: Sea, lake, and overland surges from hurricanes
    NOAA model.

    It is a simple model, and apprently when the coriolis is neglected (f -> 0)
    Emanuel 2011 model reduces to SLOSH model.

    :param r: float, radial distance (m)
    :param Rm: float, radius of maximum wind (m)
    :param Vm: float, maximum wind speed (m/s)
    :return: float, circular wind speed (m/s) at r distance
    """
    vcirc = (2 * r * Rm * Vm) / (Rm ** 2 + r ** 2)
    return vcirc


def calc_vcirc_h80(r, Rm, pc, B, pn, rhoair, f):
    """
    Calculate circular wind speed based on Holland 1980 model.

    :param r: float, radial distance (m) where the circular wind speed is calculated
    :param Rm: float, radius of maximum wind (m)
    :param pc: float, central pressure
    :param B: float, Holland B parameter, can be calculated using calc_holland_B_full()
    :param pn: float, nominal pressure outide of the storm
    :param rhoair: density of air kg/m**3
    :param f: coriolis parameter, can be calculated with coriolis() function
    :return: float, circular wind speed (m/s) at r distance
    """
    r = r + 1e-8  # Avoid divided by zero issue for r==0
    vcirc = np.sqrt(
        (Rm / r) ** B * (B / rhoair) * (pn - pc) * np.exp(-(Rm / r) ** B) + (r * f / 2) ** 2
    ) - (r * f / 2)
    return vcirc


def calc_vcirc_h80c(r, Rm, pc, B, pn, rhoair):
    """
    Holland 1980 model with cyclostropic approximation. In this approximation
    all the terms related to coriolis is neglected (f->0)

    :param r: float, radial distance (m) where the circular wind speed is calculated
    :param Rm: float, radius of maximum wind (m)
    :param pc: float, central pressure
    :param B: float, Holland B parameter, can be calculated using calc_holland_B_full()
    :param pn: float, nominal pressure outide of the storm
    :param rhoair: density of air km/m**3
    :return: float, circular wind speed (m/s) at r distance
    """
    r = r + 1e-8
    vcirc = np.sqrt((Rm / r) ** B * (B / rhoair) * (pn - pc) * np.exp(-(Rm / r) ** B))
    return vcirc


def calc_vcirc_w06(r, Rm, Vm, n=0.79, X=243000):
    """
    Calculate circular wind speed using Willoughby et al. 2006 model.

    :param r: float, radial distance in (m)
    :param Rm: float, radius of maximum wind (m)
    :param Vm: float, maximum wind speed (m/s)
    :param n: power of increasing profile, default n=0.79
    :param X: distance to which extend the profile, (m) default 243km (243000m)
    :return: float, circular wind speed (m/s) at r distance
    """
    if r <= Rm:
        vcirc = Vm * (r / Rm) ** n
    else:
        vcirc = Vm * np.exp(-(r - Rm) / X)

    return vcirc


def calc_vcirc_e04(r, Rm, Vm, b=0.25, m=1.6, n=0.9, R0=420000):
    """
    Calculate circular wind speed using Emanuel 2004 model.

    :param r: float, radial distance (m)
    :param Rm: float, radisu of maximum wind (m)
    :param Vm: float, velocity of maximum wind (m/s)
    :param b: constant, default 0.25
    :param m: constant, default 1.6
    :param n: constant, default 0.9
    :param R0: distance, (m), default 420000m (420km)
    :return: float, circular wind speed (m/s) at r distance
    """
    multiplier = ((1 - b) * (n + m) / (n + m * (r / Rm) ** (2 * (n + m)))) + (
            b * (1 + 2 * m) / (1 + 2 * m * (r / Rm) ** (2 * m + 1)))
    vcirc = Vm * ((R0 - r) / (R0 - Rm)) * (r / Rm) ** m * np.sqrt(multiplier)

    return vcirc


def calc_vcirc_e11(r, Rm, Vm, f):
    """
    Calculate circular wind speed using Emanuel 2011 model. The model equation is given in
    Lin and Chavas 2012.

    :param r: float, radial distance (m)
    :param Rm: float, radius of maximum wind (m)
    :param Vm: float, maximum wind speed (m/s)
    :param f: coriolis, can be calc by coriolis() function
    :return: float, circular wind speed (m/s) at r distance
    """
    vcirc = 2 * r * (Rm * Vm + 0.5 * f * Rm ** 2) / (Rm ** 2 + r ** 2) - r * f / 2

    return vcirc


def calc_vcirc_m16(r, Rm, Vm, n=0.6):
    """
    Calculate circular wind speed using Murty et al. 2016 model, developed for indian ocean.

    :param r: float, radial distance (m)
    :param Rm: float, radius of maximum wind (m)
    :param Vm: float, maximum wind speed (m/s)
    :param n: power, default is 0.6
    :return: float, circular wind speed (m/s) at r distance
    """
    vcirc = Vm * ((2 * r * Rm) / (Rm ** 2 + r ** 2)) ** n
    return vcirc


def calc_mslp_h80(r, Rm, pc, B, pn):
    """
    Calculate mean-sea-level pressure based using Holland 1980 model.

    :param r: float, radial distance (m) where the circular wind speed is calculated
    :param Rm: float, radius of maximum wind (m)
    :param pc: float, central pressure
    :param B: float, Holland B parameter, can be calculated using calc_holland_B_full()
    :param pn: float, nominal pressure outide of the storm
    :return: float, pressure at r distance
    """
    r = r + 1e-8
    mslp = (pn - pc) * np.exp(-(Rm / r) ** B) + pc
    return mslp
