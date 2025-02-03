#
# clusterize positions with DBSCAN, return Polygons
#
# Only accounts for core samples - edge points that are not core samples will be ignored.
#
# parameters:
#   per_individual
#   dist_threshold
#   min_core_points 
#   min_days
#

import geopandas as gpd
import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN
from shapely import MultiPoint


def clusterize(gdf):
    # dbscan
    db = DBSCAN(eps=snakemake.params['dist_threshold'], min_samples=snakemake.params['min_days'])
    dbclusters = db.fit(gdf.geometry.get_coordinates())
    gdf = gdf.copy()
    gdf['cluster'] = dbclusters.labels_
    gdf['core'] = 0
    gdf.loc[gdf.iloc[dbclusters.core_sample_indices_].index, 'core'] = 1

    # build new GeoDataframe with polygons that enclose each cluster's core samples
    index, npoints, ndays, geom = [], [], [], []
    for cl, gdf_cluster in gdf.loc[gdf.cluster >= 0].groupby('cluster'):
        n = np.unique(gdf_cluster.index.get_level_values(1).date).shape[0]
        if n < snakemake.params['min_days']:
            continue
        gdf_core = gdf_cluster.loc[gdf_cluster.core > 0]
        if gdf_core.shape[0] < snakemake.params['min_core_points']:
            continue
        if snakemake.params['per_individual']:
            index.append((gdf_cluster.index.get_level_values(0)[0], cl))
        else:
            index.append(cl)
        npoints.append(gdf_cluster.shape[0])
        ndays.append(n)
        geom.append(MultiPoint(gpd.points_from_xy(gdf_core.geometry.x, gdf_core.geometry.y)).convex_hull)

    return None if not len(geom) \
        else gpd.GeoDataFrame(
            index=pd.MultiIndex.from_tuples(index, names=['ind', 'c']) if snakemake.params['per_individual'] else index,
            data={'npoints': npoints, 'ndays': ndays},
            geometry=geom,
            crs=gdf.crs)

gdf = gpd.read_parquet(snakemake.input[0])
out = gdf.groupby(level=0, group_keys=False).apply(clusterize) if snakemake.params['per_individual'] else clusterize(gdf)
out.to_parquet(snakemake.output[0])
