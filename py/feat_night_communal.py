#
# add communal feature to night roosts, depending on whether other individuals' night roosts intersect within the same Autumn/Winter (Oct to Mar).
#
# only consider 1st-year juveniles from the Parisian cohorts.
#
# compare to all other Parisian individuals, except Rescues (G279...) as they mostly contribute noise.
#
# limitations:
# - density-dependent, both spatial and temporal: 2020 and 2021 don't benefit from the concurrency of previous cohort (data from 2020 has already dissipated
# when 2021 is equipped). The more distant the location, the poorer the data. Conversely, the capture site has overwhelming noise.
#
# - assumes territorial stability - but in reality, across the winter, some territories disappear, some shift, and others might tolerate some levels of invasion.
#
# - noise makes bordering territories collapse into the communal roosts.
#
# - some areas, esp. the nearby forests and large cemeteries, have a much looser density, eroding the matching.
#
# - singular invasions are rare but do happen, and currently are not filtered.
#

import pandas as pd
import geopandas as gpd
import numpy as np

idx = pd.IndexSlice

def feat_night_individual(gdf_juvs:gpd.GeoDataFrame, gdf_population:gpd.GeoDataFrame, radius:int=25.) -> gpd.GeoDataFrame:
    s_communal = pd.Series(0, index=gdf_juvs.index, name='communal')
    gdf_population = gdf_population.copy()
    gdf_population.geometry = gdf_population.geometry.buffer(radius)
    for ind, gdf_juv in gdf_juvs.groupby(level='ind'):
        match = gdf_juv.sjoin(gdf_population.loc[gdf_population.index.get_level_values(1)!=ind], lsuffix=None)
        s_communal.loc[match.index] = 1
    return s_communal

gdf_nights = gpd.read_parquet(snakemake.input[0])
gdf_nights = gdf_nights.loc[~gdf_nights.index.get_level_values(0).isin(['RESCUE','RURAL'])]

s_communal = pd.concat([
    feat_night_individual(gdf_nights.loc[idx[['2020'],:,:'2021-03-31'],:], gdf_nights.loc[idx[:,:,'2020-09-01':'2021-03-31'],:], snakemake.params['radius']),
    feat_night_individual(gdf_nights.loc[idx[['2021'],:,:'2022-03-31'],:], gdf_nights.loc[idx[:,:,'2021-09-01':'2022-03-31'],:], snakemake.params['radius']),
    feat_night_individual(gdf_nights.loc[idx[['2022'],:,:'2023-03-31'],:], gdf_nights.loc[idx[:,:,'2022-09-01':'2023-03-31'],:], snakemake.params['radius']),
    feat_night_individual(gdf_nights.loc[idx[['2024'],:,:'2025-03-31'],:], gdf_nights.loc[idx[:,:,'2024-09-01':'2025-03-31'],:], snakemake.params['radius'])
])

pd.concat([gdf_nights.loc[s_communal.index,:], s_communal], axis=1) \
        .to_parquet(snakemake.output[0])
