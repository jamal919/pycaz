# -*- encoding: utf-8 -*-
import numpy as np
import cmocean.cm as cmo
import matplotlib.colors as mcolors

SCALES = {
    "saffir": {
        "bounds": np.array([0, 10, 18, 33, 43, 50, 58, 70, 100]),
        "ticklabels": ['NA', 'TD', 'TS', '1', '2', '3', '4', '5']
    }
}


def get_storm_cmap(scale: dict | str = None, from_cmap=cmo.matter) -> dict:
    """Generate cmaps for storms based on scales

    Args:
        scale (dict | str, optional): Scale to use. Can be a dict(bounds, ticklabels), or scale names. Defaults to saffir.
        from_cmap (optional): Colormap to pull from. Defaults to cmo.matter.

    Raises:
        KeyError: If the scale name is not found in the defines storm scales
        ValueError: If the bounds and ticklabels (+1) lengths are equal

    Returns:
        dict: A dictionary containing the norm, cmap, ticklocs, and ticklabels
    """
    # Parameter parsing
    if scale is None:
        scale = "saffir"

    if isinstance(scale, str):
        try:
            scale_def = SCALES[scale]
        except KeyError:
            raise KeyError(f"{scale} scale not found, avilable scales - {SCALES.keys()}")

    if isinstance(scale, dict):
        scale_def = scale

    bounds = np.atleast_1d(scale_def["bounds"])
    ticklabels = scale_def["ticklabels"]

    # Bounds length should be len(ticklabels) + 1
    try:
        assert len(bounds) == len(ticklabels) + 1
    except AssertionError:
        raise ValueError("Bounds length needs to be len(ticklabels) + 1")

    # Colorscales for cyclone
    cmap = from_cmap
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmaplist[0] = (.5, .5, .5, 1.0)  # First one to gray

    # create the new map
    cmap = mcolors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
    ticklocs = bounds[1:] - np.diff(bounds) / 2
    norm = mcolors.BoundaryNorm(bounds, cmap.N)

    obj = {
        "norm": norm,
        "cmap": cmap,
        "ticklocs": ticklocs,
        "ticklabels": ticklabels
    }

    return obj
