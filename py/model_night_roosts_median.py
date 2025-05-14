#
# estimate where each night's roost took place for each individual.
#
# use the median lat / long, which buffers some of the GPS uncertainty and smoothes the outright outliers. It also helps during data scarcity
# where only 1 fix/night might be recorded.
#
# limitations:
#  - genuine night movement (change of roosting location, early activity, ...) reduces the accuracy or might select intermediate locations.
#
# takes: GeoDataFrame
# returns: GeoDataFrame with median lat/long, geometry in the original GeoDataFrame CRS, and the count and std (in the unit of the CRS) of the points
# that made up the night. Timestamps will be raw dates. Night on 2024-01-01 is considered to start the 01 in the evening and stop on 02 morning. Sunpos
# is also defined for compatibility.

import pandas as pd
import geopandas as gpd
import numpy as np

def m_night_roosts_median(gdf:gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    gdf = gdf.loc[gdf.sunpos=='night'].copy()
    gdf.index = gdf.index.set_levels(gdf.index.levels[1].shift(-12, 'h'), level=1)
    std = gdf.get_coordinates().groupby([pd.Grouper(level='ind'),pd.Grouper(level='ts', freq='d')]).std()
    df = pd.concat([
        gdf.groupby([pd.Grouper(level='ind'),pd.Grouper(level='ts', freq='d')]) \
            .agg(lat=('lat','median'), long=('long','median'), sunpos=('sunpos','first'), n=('lat','count')),
        np.sqrt(std.x**2 + std.y**2).rename('std')
    ], axis=1)
    return gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.long, df.lat, crs='EPSG:4326')).to_crs(gdf.crs)

m_night_roosts_median(gpd.read_parquet(snakemake.input[0])) \
        .to_parquet(snakemake.output[0])
