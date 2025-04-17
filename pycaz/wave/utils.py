# -*- coding: utf-8 -*-

import numpy as np
from scipy import stats


def isnum(x):
    try:
        x = float(x)
        r = True
    except:
        r = False
    return r


def filter_sinc(x, n):
    """
    A sinc-based filter

    :param x: the array to be filtered
    :param n: number of samples
    :return: the filtered array
    """
    f = 1 / n
    # print(f'f = {f}')
    nf = 2 * np.fix(2 / f)
    o = np.arange(1, nf + 1) - nf / 2
    filt = 0 * o + 2 * f
    selouter = np.asarray(np.hstack((np.arange(np.floor(nf / 2) - 1), np.arange(np.floor(nf / 2), nf))), dtype=int)
    filt[selouter] = np.sin(2 * np.pi * f * o[selouter]) / o[selouter] / np.pi
    hamm = 0.54 - 0.46 * np.cos(2 * np.pi * np.arange(1, nf + 1) / nf)
    filt = filt * hamm
    xf = np.convolve(x, filt, mode='same')
    return xf


def wavenum(f, h):
    """
    Compute wavenumber (k), by solving the dispersion relation solved using a polynomial approximation.

    Reference: George Voulgaris, SUDO, 1992

    :param f: wave frequency; f=1/T
    :param h: water depth (in m)
    :return:
    """

    w = 2 * np.pi * f
    dum1 = (w ** 2) * h / 9.81
    dum2 = dum1 + (1.0 + 0.6522 * dum1 + 0.4622 * dum1 ** 2 + 0.0864 * dum1 ** 4 + 0.0675 * dum1 ** 5) ** (-1)
    dum3 = np.sqrt(9.81 * h * dum2 ** (-1)) / f
    y = 2 * np.pi * dum3 ** (-1)

    return y


def compute_wavelength(wave_period, d_mean):
    """
    Computes wavelength based on linear wave theory

    :param wave_period: Wave period in seconds
    :param d_mean: Mean water depth in meters
    :return: wavelength in meters
    """

    return (2 * np.pi) / wavenum(1 / wave_period, d_mean)


def tfm_nh_correction(psd, freqs, h_mean, zpt):
    """
    Non-hydrostatic correction for buried sensor using transfer function method.

    :param psd: Computed psd on water level computed using hydrostatic approximation
    :param freqs: Corresponding frequencties of the input psd
    :param h_mean: Mean water depth at the water level location
    :param zpt: Depth of sensor from the sea bottom (e.g., buried sensor)
    :return: Corrected psd
    """
    k = wavenum(freqs, h_mean)
    kp = np.cosh(k * zpt) / np.cosh(k * h_mean)
    psd_cor = psd / kp ** 2

    return psd_cor


def calculate_spectral_ci(data, nfft, ci):
    """
    Calculate spectral Confidence interval

    :param data: Spectral density output from psd()
    :param nfft: Number of samples including for 50% overlap
    :param ci: Target confidence interval, e.g., 0.95
    :return: (ci_upper, ci_lower)
    """
    data_var = np.var(data, ddof=1)
    data_dof = (round(len(data) / nfft) * 2 - 1) * 2
    alpha = 1 - ci

    ci_upper = data_dof * data_var / stats.chi2.ppf(alpha / 2, df=data_dof)
    ci_lower = data_dof * data_var / stats.chi2.ppf(1 - alpha / 2, df=data_dof)

    return ci_upper, ci_lower


def calculate_Cg(f, h, shallow=False):
    g = 9.81
    if shallow:
        Cg = np.sqrt(g * h)
    else:
        k = wavenum(f, h)
        Cp = np.sqrt(g * np.tanh(k * h) / k)
        Cg = Cp * ((k * h / np.sinh(2 * k * h)) + 1 / 2)

    return Cg


def compute_cp_cg(f: float | np.ndarray, h: float | np.ndarray) -> list:
    """
    Compute the group velocity (Cg) and for a specific frequency (f) and depth (h)
    :param f: Frequency (Hz) for which the velocities are computed
    :param h: Depth (m) for which the velocities are computed
    :return:
    """
    if np.any(h < 0):
        raise ValueError('Depth h has to be positive')

    g = 9.81
    k = wavenum(f, h)
    Cp = np.sqrt(g * np.tanh(k * h) / k)
    Cg = Cp * ((k * h / np.sinh(2 * k * h)) + 1 / 2)

    return Cp, Cg
