#
# estimate daily range based on existing points
#
# For each individual and each day, enclose points inside the smallest possible rotated rectangle, then return the length of the diagonal as 
# an approximation of the covered area.
# This doesn't discount night roosts. The method is vulnerable to data scarcity.
#

import geopandas as gpd
import pandas as pd
import numpy as np
from shapely import MultiPoint, LineString

def range_rectangle(gdf:gpd.GeoDataFrame) -> np.ndarray:
    if gdf.shape[0] <= 1:
        return np.nan
    geom = MultiPoint(gpd.points_from_xy(gdf.geometry.x, gdf.geometry.y)).minimum_rotated_rectangle
    if type(geom) == LineString:
        return geom.length
    xx, yy = geom.exterior.coords.xy
    return np.sqrt((xx[2]-xx[0])**2 + (yy[2]-yy[0])**2)

gpd.read_parquet(snakemake.input[0]).groupby([pd.Grouper(level=0),pd.Grouper(level=1, freq='d')]).agg(
    cohort=pd.NamedAgg(column='cohort', aggfunc='min'),
    age=pd.NamedAgg(column='age', aggfunc='max'),
    battery=pd.NamedAgg(column='battery', aggfunc='median'),
    npoints=pd.NamedAgg(column='lat', aggfunc='count'),
    range=pd.NamedAgg(column='geometry', aggfunc=range_rectangle)
).to_parquet(snakemake.output[0])
