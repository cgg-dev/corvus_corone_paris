import geopandas as gpd
import pandas as pd
from shapely import LineString
import folium
import colorcet as cc
from typing import Optional

def _traj(gdf:gpd.GeoDataFrame, m:Optional[folium.Map]=None, linewidth:int=2) -> folium.Map:
    return gdf.groupby(level=0).filter(lambda x: x.index.shape[0] > 1) \
        .groupby(level=0).agg({'geometry': lambda g: LineString(g.values)}) \
        .set_geometry('geometry') \
        .set_crs(gdf.crs) \
        .reset_index() \
        .explore(column='ind', cmap=cc.cm.glasbey, m=m, style_kwds={'weight': linewidth})

def _expl(gdf:gpd.GeoDataFrame, m:Optional[folium.Map]=None) -> folium.Map:
    cmap = {'day':'blue','night':'black','twilight':'orange'}
    m = _traj(gdf, m)
    return gdf.reset_index().explore(color=[cmap[x] for x in gdf.sunpos], m=m, legend=False)
