#
# Resample to 1 fix/h per individual using the nearest neighbor method, with a maximum limit (default 30m).
#
# This controls the amount of interpolation since the resampled data is used in dyads analysis (timestamps must match).
# 
# The limit parameter restricts the amount of interpolation to increments of 10 minutes. For example, limit=2 allows a maximum
# of 20 minutes of interpolation, which sets the temporal distance for a dyadic match at a maximum of 40 minutes when data
# is sparse (interpolation over 20m for both sides of the dyad).
# This works by ignoring fixes that fall outside of the limit boundary, i.e. for limit=2, timestamps between hh:20 and hh:40
# will be dropped.
#
# WARNING: limits lower than the default (3) will drop data. Proceed as you see fit.
#
# The dataset offers little in ways of synchronization and equipment can vary significantly with vendors. Transmitters are fairly
# consistent as long as they hold charge, and because timing is regular, entire sequences might be dropped if the fixes fall
# within the interpolation limit's blind spot - e.g. the transmitter captures fixes at hh:24 every hour, when limit=2.
#
# Timestamps that fall exactly at hh:30:00 and match 2 timesteps are rounded up to hh+1:00:00.
#

import pandas as pd
import geopandas as gpd
import numpy as np

idx = pd.IndexSlice


gdf = gpd.read_parquet(snakemake.input[0]).groupby(level=0).apply(lambda x: x.droplevel(0).resample('10 min').nearest(limit=snakemake.params['limit']).dropna(how='all'))
timesteps = pd.date_range(
    pd.to_datetime(gdf.index.get_level_values(1).min().date(), utc=True),
    pd.to_datetime(gdf.index.get_level_values(1).max().date(), utc=True) +pd.Timedelta('1d'),
    freq='1h')
timesteps = timesteps.intersection(gdf.index.get_level_values(1))
gdf.loc[idx[:,timesteps],:] \
        .drop_duplicates(keep='last') \
        .to_parquet(snakemake.output[0])
