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
    # y=wavenum(f,h): FUNCTION for the calculation of the wavenumber.
    #                   The dispertion relation is solved using a
    #                   polynomial approximation.
    #                   f, wave frequency; f=1/T.
    #                   h, water depth (in m).
    #
    #       George Voulgaris, SUDO, 1992
    w = 2 * np.pi * f
    dum1 = (w ** 2) * h / 9.81
    dum2 = dum1 + (1.0 + 0.6522 * dum1 + 0.4622 * dum1 ** 2 + 0.0864 * dum1 ** 4 + 0.0675 * dum1 ** 5) ** (-1)
    dum3 = np.sqrt(9.81 * h * dum2 ** (-1)) / f
    y = 2 * np.pi * dum3 ** (-1)

    return y


def compute_wavelength(T, d):
    """
    Computes wavelength based on linear wave theory

    Parameters
    ----------
    T : integer, float
        Wave period in seconds.
    d : integer, float
        Water depth in meters.

    Returns
    -------
    Wavelength in meters.

    """

    return (2 * np.pi) / wavenum(1 / T, d)


def TFM_nh_correction(PSD, freqs, h_mean, Zpt):
    k = wavenum(freqs, h_mean)
    Kp = np.cosh(k * Zpt) / np.cosh(k * h_mean)
    PSD_cor = PSD / Kp ** 2

    return PSD_cor


def calculate_spectral_CI(data, Nfft, CI):
    data_var = np.var(data, ddof=1)
    data_dof = (round(len(data) / Nfft) * 2 - 1) * 2
    alpha = 1 - CI

    CI_upper = data_dof * data_var / stats.chi2.ppf(alpha / 2, df=data_dof)
    CI_lower = data_dof * data_var / stats.chi2.ppf(1 - alpha / 2, df=data_dof)

    return CI_upper, CI_lower


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


