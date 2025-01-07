#
# filter out non-physiological points based on excessive speed or distance.
#
# Note that this is very crude and only helps so far that it removes the most obvious inconsistencies.
#

import pandas as pd
import geopandas as gpd


# Remove points that have speed >= kmh *to* and *from* that location. Other situations (one side only >= kmh) are trickier to solve.
def filter_speed(gdf:gpd.GeoDataFrame, kmh:int=65, verbose:bool=False) -> gpd.GeoDataFrame:
    speeds = gdf.groupby(level=0, group_keys=False).apply(lambda x: x.distance(x.shift()) / (x.index.get_level_values(1).diff() / pd.Timedelta('1h'))).rename('speed')
    speeds *= 1e-3
    rm = speeds.loc[(speeds >= kmh) & (speeds.shift(-1) >= kmh)]
    if verbose:
        print(f"filtered {rm.shape[0]} points - speed >= {kmh}")
    return gdf.drop(rm.index)

# Remove points that are separated from the preceding and following points by x km. Speed filtering won't pick up on such inconsistencies when data is scarce, though is still the preferred method.
def filter_distance(gdf:gpd.GeoDataFrame, km:int=50, verbose:bool=False) -> gpd.GeoDataFrame:
    dist = gdf.groupby(level=0, group_keys=False).apply(lambda g: g.distance(g.shift())).rename('dist')
    dist *= 1e-3
    rm = dist.loc[(dist >= km) & (dist.shift(-1) >= km)]
    if verbose:
        print(f"filtered {rm.shape[0]} points - distance >= {km}")
    return gdf.drop(rm.index)

gdf = gpd.read_parquet(snakemake.input[0])
gdf = filter_speed(gdf, snakemake.config['filtering']['max_kmh'], True)
gdf = filter_distance(gdf, snakemake.config['filtering']['max_km'], True)

gdf.to_parquet(snakemake.output[0])

