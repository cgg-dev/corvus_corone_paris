import geopandas as gpd
import pandas as pd
from shapely import LineString
import folium
import xyzservices.providers as xyz
import xyzservices
import colorcet as cc
import matplotlib as mpl

def _traj(gdf:gpd.GeoDataFrame,
        provider:xyzservices.TileProvider | str=None,
        cmap:mpl.colors.ListedColormap | mpl.colors.LinearSegmentedColormap | str=cc.cm.glasbey,
        m:folium.Map | None =None,
        linewidth:int=2) -> folium.Map:
    if not provider:
        provider = 'OpenStreetMap Mapnik'
    gdf_traj = gdf.groupby(level=0).filter(lambda x: x.index.shape[0] > 1) \
        .groupby(level=0).agg({'geometry': lambda g: LineString(g.values)}) \
        .set_geometry('geometry') \
        .set_crs(gdf.crs)
    if gdf_traj.shape[0] == 0:
        return None
    return gdf_traj.reset_index().explore(column=gdf.index.names[0], cmap=cmap, tiles=provider, m=m, style_kwds={'weight': linewidth})

def _expl(gdf:gpd.GeoDataFrame, provider:xyzservices.TileProvider | str | None=None, cmap=cc.cm.glasbey, m:folium.Map | None=None) -> folium.Map:
    sunpos_cmap = {'day': 'blue', 'night': 'black', 'twilight': 'orange'}
    m = _traj(gdf, provider, cmap, m)
    return gdf.reset_index().explore(color=[sunpos_cmap[x] for x in gdf.sunpos], m=m, legend=False)
