#
# take 'areas' geometries (Polygons) and determine the co-occurrence for each geometry - i.e. the number of distinct individuals
# intersecting with this space.
#
# input:
#   geodataframe with 'areas' geometries
#   population geodataframe
#
# parameters:
#   min_days    min distinct days of presence for an individual to be considered.
#   restrict_y  expects a 'year' column in the areas geodataframe and only consider said year for matching positions.
#

import geopandas as gpd
import pandas as pd
import numpy as np


def cooccurence(area):
    gdf_match = gdf_population
    if snakemake.params.get('restrict_y'):
        gdf_match = gdf_match.loc[gdf_match.index.get_level_values(1).year == area.year.iloc[0]]
    cooccur = area.sjoin(gdf_match, rsuffix=None) \
        .groupby(pd.Grouper('ind')) \
        .filter(lambda x: np.unique(x.ts.dt.date).shape[0] >= snakemake.params['min_days']) \
        .groupby(pd.Grouper('ind')).lat.count().count()
    area.insert(len(area.columns)-1, 'cooccur', cooccur)
    return area

areas = gpd.read_parquet(snakemake.input[0])
gdf_population = gpd.read_parquet(snakemake.input[1])
areas.groupby(level=list(range(areas.index.nlevels)), group_keys=False).apply(cooccurence).to_parquet(snakemake.output[0])
