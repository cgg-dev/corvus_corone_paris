#
# featurize raw movebank dataframe - geometrize and add sunpos.
#

import pandas as pd
import geopandas as gpd
from suncalc import get_times


# return a 'sunpos' Series with values 'twilight', 'day', 'night'.
# Daytime is defined conservatively, by enhancing the duration of twilight - i.e. daytime starts and ends when
# the *base* (lower part) of the sun touches the horizon. This is to better isolate the transitional aspect of
# twilight and provide a clearer picture of the animal's behavior. 
# Night is more intricate since 'actual' astronomical night doesn't exist or is skewed in summer for high latitudes.
# Nautical twilight is used instead, which is geographically safer but less 'strict'.
# This setup *will* fail nonetheless for locations that are very up north or very down south.
def feat_sunpos(df:pd.DataFrame) -> pd.Series:
    ts = df.index.get_level_values(1)
    suntimes = pd.DataFrame(get_times(ts, df.long, df.lat)).set_axis(df.index)
    odf = pd.Series(data='twilight', index=df.index, name='sunpos')
    odf.loc[(ts < suntimes.nautical_dawn) | (ts > suntimes.nautical_dusk)] = 'night'
    odf.loc[(ts > suntimes.sunrise_end) & (ts < suntimes.sunset_start)] = 'day'
    return odf

def geometrize(df:pd.DataFrame, epsg:int) -> gpd.GeoDataFrame:
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.long, df.lat, crs='EPSG:4326'))
    return gdf.to_crs(epsg=epsg)


df = pd.read_parquet(snakemake.input[0])
# warnings from suncalc are expected (because astronomical night doesn't exist in Paris around the summer solstice)
df['sunpos'] = feat_sunpos(df)
gdf = geometrize(df, snakemake.config['epsg'])
gdf.to_parquet(snakemake.output[0])
