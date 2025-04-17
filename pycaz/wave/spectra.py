# -*- coding: utf-8 -*-

import numpy as np

from pycaz.wave.utils import compute_cp_cg


def compute_bulk_params(freq, efth, freq_low=None, freq_high=None):
    """
    Compute Bulk parameter from a 1D Spectra using method of moements.

    The integratin is done using trapozoidal rule with np.trapz. The low and high frequency cutoff
    can be provided, however it only slices, and do not ensure exact boundary. For making sure an
    exact boundary, use an interpolation function before passing here.

    Current list of bulk parameters:
        - Hm0: Significant wave height
        - Tm02: Mean period
        - Tp: Peak period
        - DTp: Discrete peak period

    :param freq: The list of frequencies for which
    :param efth: Spectral energy density
    :param freq_low: Lower cutoff frequency (typically from the buoy)
    :param freq_high: Higher cutoff frequency (typically from the buoy)
    :return: Dictionary of Bulk parameters

    TODO: Make it more modular
    """
    # Cutoff processing
    if freq_low is None:
        freq_low = np.min(freq)

    if freq_high is None:
        freq_high = np.max(freq)

    # Selections
    sel_freq = np.logical_and(
        freq >= freq_low,
        freq <= freq_high
    )

    freq_ = freq[sel_freq]
    efth_ = efth[sel_freq]

    # Moments
    m0_ = np.trapz(efth_, freq_)
    m1_ = np.trapz(efth_ * freq_, freq_)
    m2_ = np.trapz(efth_ * freq_ ** 2, freq_)
    m_2_ = np.trapz(efth_ * freq_ ** -2, freq_)

    # Bulk parameters
    hm0_ = 4 * np.sqrt(m0_)
    tm02_ = np.sqrt(m0_ / m2_)
    if np.any(np.isnan(efth_)):
        dtp_ = np.nan
    else:
        dfp_ = freq_[np.argmax(efth_)]
        dtp_ = 1 / dfp_
    fpc_ = (m0_ ** 2) / (m_2_ * m1_)
    tpc_ = 1 / fpc_

    return {
        'Hm0': hm0_,
        'Tm02': tm02_,
        'Tpc': tpc_,
        'Tp': dtp_
    }


def correct_shoaling(
        freq: float | np.ndarray,
        efth: float | np.ndarray,
        depth: float | np.ndarray,
        target_depth: float | np.ndarray) -> float | np.ndarray:
    """
    Correct the energy with a simple shoaling correction based on conservation of energy flux between depths

    :param freq: Frequency (Hz) for which the energy (efth) is provided.
    :param efth: Spectral energy.
    :param depth: Current depth in (m) or consistent units
    :param target_depth: The target depth (m) for which the corrected energy (efth) is provided.
    :return: Corrected efth
    """
    _, cg_in = compute_cp_cg(freq, depth)
    _, cg_out = compute_cp_cg(freq, target_depth)
    corr_factor = cg_in / cg_out
    efth_target_depth = efth * corr_factor

    return efth_target_depth
