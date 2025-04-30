# -*- coding: utf-8 -*-

def get_country_geometry(country_name: str):
    """
    Get the country geometry from natural_earth

    :param country_name: Name of the country to get geometries
    :return: Shapely geometry
    """
    import cartopy.io.shapereader as shpreader

    fn = shpreader.natural_earth(
        category="cultural",
        name="admin_0_countries",
        resolution="10m"
    )

    reader = shpreader.Reader(fn)
    countries = reader.records()
    country = next(countries)
    while country.attributes["NAME"] != country_name:
        country = next(countries)

    return country.geometry