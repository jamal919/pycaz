# -*- coding: utf-8 -*-

import numpy as np

# Ref 
# Sonel: https://www.sonel.org/Filters-for-the-daily-mean-sea.html
# Pugh: https://eprints.soton.ac.uk/19157/1/sea-level.pdf
filters = {
    "Average": {
        "Filter": np.array([
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
        ]),
        "Denom": 25
    },
    "Demerliac": {
        "Filter": np.array([
            1, 3, 8, 15, 21, 32, 45, 55, 72, 91, 105, 128, 153, 171, 200, 231, 253, 288, 325, 351,
            392, 435, 465, 512, 558, 586, 624, 658, 678, 704, 726, 738, 752, 762, 766,
            768,
            766, 762, 752, 738, 726, 704, 678, 658, 624, 586, 558, 512, 465, 435, 392,
            351, 325, 288, 253, 231, 200, 171, 153, 128, 105, 91, 72, 55, 45, 32, 21, 15, 8, 3, 1]),
        "Denom": 24576
    },
    "Doodson": {
        "Filter": np.array([
            1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 2, 0, 1, 1, 0, 2, 1, 1, 2,
            0,
            2, 1, 1, 2, 0, 1, 1, 0, 2, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1]),
        "Denom": 30
    },
    "Munk": {
        "Filter": np.array([
            13307, 30073, 47028, 60772, 72261, 85349, 101603, 122665,
            146225, 165525, 180727, 195518, 208050, 219260, 234033, 251492,
            278167, 300054, 314959, 325633, 338603, 354118, 370094, 386839,
            395287,
            386839, 370094, 354118, 338603, 325633, 314959, 300054, 278167,
            251492, 234033, 219260, 208050, 195518, 180727, 165525, 146225,
            122665, 101603, 85349, 72261, 60772, 47028, 30073, 13307
        ]),
        "Denom": 1e7
    },
    "Godin": {
        "Filter": np.array([
            1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66, 78, 91,
            105, 120, 136, 153, 171, 190, 210, 231, 253, 276, 300, 323, 344,
            363, 380, 395, 408, 419, 428, 435, 440, 443,
            444,
            443, 440, 435, 428, 419, 408, 395, 380, 363,
            344, 323, 300, 276, 253, 231, 210, 190, 171, 153, 136, 120, 105,
            91, 78, 66, 55, 45, 36, 28, 21, 15, 10, 6, 3, 1
        ]),
        "Denom": 14400
    }
}


def apply_filter(in_array, filter_name):
    """
    Apply `filter_name` tide filter to `in_array` hourly dataset.

    :param in_array: Input hourly timeseries where the tide filter to be applied.
    :param filter_name: Name of filter to be applied.
    :return: Output hourly timeseries.
    """
    temp_array = np.ones(len(in_array)) * np.nan
    filter_info = filters.get(filter_name)
    crop_len = len(filter_info["Filter"]) // 2
    filter_array = filter_info["Filter"]
    filter_denom = filter_info["Denom"]
    temp_array[crop_len:-1 * crop_len] = np.convolve(in_array, filter_array / filter_denom, mode="valid")
    return temp_array
