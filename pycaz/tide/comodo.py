# -*- coding: utf-8 -*-
import numpy as np
from .utide import utide_resolved_consts
from pycaz.typing import TimestampConvertibleTypes

COMODO_CONSTS = np.array([
    'Z0', 'M0', 'LP2', 'LP3', 'Sa_1', 'Sa', 'Ssa_1', 'Ssa', 'Sta',
    'MSm', 'Mm', 'Mm_2', 'Mm_1', 'MSf', 'MSf_nL', 'Mf_1', 'Mf_2', 'Mf',
    'MStm_1', 'MStm', 'Mtm', 'MSqm', 'Mqm', '3OK1', '2Q1', 'Sig1',
    'Q1', 'Ro1', 'O1', 'MS1', 'MP1', 'M1_1', 'M1_3', 'M1', 'M1_2',
    'Ki1', 'Pi1', 'P1', 'S1', 'K1', 'Psi1', 'Phi1', 'Tta1', 'J1',
    'SO1', 'OO1', 'KQ1', '2MN2S2', '2NS2', '3M2S2', 'ST1', 'OQ2', 'E2',
    'MNS2', 'MNuS2', 'ST2', 'ST3', '2MK2', '2N2', 'Mu2', '2MS2',
    'SNK2', 'A76', 'N2', 'Nu2', 'MSK2', 'OP2', 'M(SK)2', 'M2',
    'M(KS)2', 'MKS2', 'La2', 'L2', 'A79', 'NKM2', 'T2', 'S2', 'R2',
    'K2', 'MSN2', 'KJ2', '2SM2', 'SKM2', '2SN2', '2SMu2', 'MQ3',
    '2MK3', 'MO3', 'M3', 'SO3', 'MS3', 'MK3', 'SP3', 'S3', 'SK3', 'K3',
    '2MNS4', 'N4', '3MS4', 'MN4', 'MNu4', 'MA4', 'M4', '2MKS4', 'SN4',
    'ML4', 'NK4', 'MS4', 'MK4', '2MSN4', 'S4', 'SK4', '2MQ5', '2MO5',
    '2NK5', '3MS5', '3MP5', 'M5', '2MP5', '2MS5', '2MK5', 'NSK5',
    '3MQ5', 'MSP5', 'MSK5', 'S5', '3MNK6', '3MNS6', '4MK6', '3MNL6',
    '2MN6', '2MNu6', '3MSK6', 'M6', '3MKS6', 'MSN6', '2ML6', 'MSNu6',
    'MNK6', '2MS6', '3MLN6', '2MK6', 'MSL6', '3MSN6', '3MKN6', '2SM6',
    'MSK6', 'S6', 'S7', 'N8', '4MNS8', '2M2N8', '3MN8', '3MNu8', 'M8',
    '2MSN8', '3ML8', '3MS8', '3MS8', '3MK8', 'MSNK8', '2M2S8', '2M2S8',
    '2MSK8', '2M2K8', '2SKN8', 'S8', 'S9', 'M10', 'S10', 'S11'
    ])

def get_comodo_names(consts:list) -> list:
    """Get the names that exists in the comodo detidor

    Args:
        consts (list): list of constituent names

    Returns:
        list: list of constituents that exists in comodo detidor lists
    """
    comodo_names = list(set(COMODO_CONSTS).intersection(set(consts)))
    return comodo_names


def comodo_resolved_names(rnday: float, dt: float = 1.0,
                        tref: TimestampConvertibleTypes = "2020-01-01",
                        resolver=utide_resolved_consts) -> list:
    """ Get the list of comodo constituents

    The list is first computed by resolved_constituent() from pycaz.tide.utide module

    :param rnday: (float) number of days of data
    :param dt: (float, optional) timestep of the data in hour. Defaults to 1.0.
    :param tref: (TimestampConvertibleTypes) Reference time, default to "2020-01-01"
    :param resolver: (callable) A function that can resolves the list of constituents and returns a dataframe
    :return: list of comodo constituent names
    """
    consts = resolver(rnday=rnday, dt=dt, tref=tref)
    comodo_names = get_comodo_names(consts.name)

    return comodo_names

