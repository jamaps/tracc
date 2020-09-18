import numpy as np
import pandas as pd
import geopandas as gpd
import libpysal

def radius(spatial_data_file_path, id_field):
    """
    returns dataframe of the radius of the circle which has the same area of each input polygon
    """
    gdf = gpd.read_file(spatial_data_file_path)
    crs = gdf.crs
    gdf = gdf.to_crs({'proj':'cea'})
    gdf["area"] = gdf["geometry"].area / 10**6
    gdf["radius"] = (gdf["area"] / 3.14159) ** 0.5
    gdf = gdf[[id_field,"radius"]]
    return gdf

def area(spatial_data_file_path, id_field):
    """
    returns dataframe of the area of each input polygon
    """
    gdf = gpd.read_file(spatial_data_file_path)
    crs = gdf.crs
    gdf = gdf.to_crs({'proj':'cea'})
    gdf["area"] = gdf["geometry"].area / 10**6
    gdf = gdf[[id_field,"area"]]
    return gdf

def get_neighbours(spatial_data_file_path, weight_type, idVariable, param = None):

    gdf = gpd.read_file(spatial_data_file_path)

    if weight_type == "Queen":

        w = libpysal.weights.Queen.from_dataframe(gdf,ids = idVariable)

    if weight_type == "KNN":

        if param is None:
            kn = 5
        else:
            kn = param

        w = libpysal.weights.KNN.from_dataframe(gdf,ids = idVariable, k = kn)

    return w.neighbors
