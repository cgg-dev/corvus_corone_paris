#
# load raw movebank csv into usable pandas dataframe
#

import pandas as pd
import numpy as np


#
# load
#

colmap = {
	'individual-local-identifier': 'ind',
	'timestamp': 'ts',
	'location-lat': 'lat',
	'location-long': 'long',
	'battery-charge-percent': 'battery',
    'orn:transmission-protocol': 'protocol' # used for filtering
}
dtypes = {
	'location-lat': np.float64,
	'location-long': np.float64,
	'battery-charge-percent': np.float64,
    'orn:transmission-protocol': str
}

df = pd.read_csv(snakemake.input[0], dtype=dtypes, usecols=list(colmap.keys())).rename(colmap, axis=1)
# Some 2024 tags were configured with SMS transmission on top of GPRS, presumably as a fallback for the rural cohort.
# This only resulted in low-quality duplicates -> remove these points
df = df.loc[df.protocol != 'SMS'].drop(columns='protocol')
df['ts'] = pd.to_datetime(df.ts, utc=True)
df = df.set_index(['ind', 'ts']).sort_index()
# Remove the extraneous FRP identifier from individual names ; drop the 5 rescued individuals from the 2022 release experiment.
# This experiment didn't yield data - none of the released individuals survived past a couple days.
df = df.rename(index={x:x[0:4] for x in df.index.get_level_values(0)}, level=0).drop('[FRP', axis=0)
idf = pd.read_csv(snakemake.input[1], index_col='ind')


#
# filter
#

# Filter points that are marked as non-physiological in individuals.csv.
# These are from GPS tags that are known to be 'lost' but still exist in the dataset.
# In some cases these points have been removed upstream at the level of movebank.
def filter_losttags(df:pd.DataFrame, individuals:pd.DataFrame, verbose:bool=False) -> pd.DataFrame:
    jdf = df.join(individuals)
    rm = jdf.loc[jdf.index.get_level_values(1) > jdf.tag_lost]
    if verbose:
        print(f"filtered {rm.shape[0]} points - lost GPS tags")
    return df.drop(rm.index)

df = filter_losttags(df, idf, True)


#
# featurize
#

def feat_cohort(df:pd.DataFrame, idf:pd.DataFrame) -> pd.Series:
    return df.join(idf).cohort

def feat_age_at_point(df:pd.DataFrame, idf:pd.DataFrame, cutoff_md:str='01-01') -> pd.Series:
    m, d = cutoff_md.split('-')
    jdf = df.join(idf, rsuffix='_r')
    ts = jdf.index.get_level_values(1)
    age = ts.year - jdf.birthyear
    age.loc[ts.month < float(m)] -= 1
    age.loc[(ts.month == float(m)) & (ts.day < float(d))] -= 1
    return age

df['cohort'] = feat_cohort(df, idf)
df['age'] = feat_age_at_point(df, idf, snakemake.config['age_cutoff'])


#
# output
#

df.to_parquet(snakemake.output[0])
