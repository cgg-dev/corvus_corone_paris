configfile: "pipeline.yaml"

# pseudo rule
rule all:
    input:
        "data/processed/movement_daily_range.parquet",
        "data/geometries/natal_territories.gpkg"

#
# transform raw data
#
rule load:
    input:
        "data/raw/movebank.csv",
        "data/raw/individuals.csv"
    output:
        temp("data/processed/raw.parquet")
    script:
        "py/load.py"

rule featurize:
    input:
        "data/processed/raw.parquet"
    output:
        temp("data/processed/featurized.parquet")
    script:
        "py/featurize.py"

rule filter:
    input:
        "data/processed/featurized.parquet"
    output:
        "data/processed/clean.parquet"
    script:
        "py/filter.py"

rule model_make_night_roosts_median:
    input:
        "data/processed/clean.parquet"
    output:
        "data/processed/m_night_roosts_median.parquet"
    script:
        "py/model_night_roosts_median.py"

rule make_feat_night_communal:
    params:
        radius=config['feat_night_communal']['radius'],
    input:
        "data/processed/m_night_roosts_median.parquet"
    output:
        "data/processed/feat_night_communal.parquet"
    script:
        "py/feat_night_communal.py"

rule make_feat_1styear_communal_roosting:
    input:
        "data/processed/feat_night_communal.parquet"
    output:
        "data/processed/feat_1styear_communal_roosting.parquet"
    run:
        import pandas as pd
        import geopandas as gpd
        gpd.read_parquet(input[0]) \
            .groupby([pd.Grouper(level=0),pd.Grouper(level=1)]).apply(lambda x:x.communal.sum()/x.communal.count()) \
            .rename('1styear_communal_roosting') \
            .to_frame().to_parquet(output[0])

#
# extract natal territories candidates
#
rule extract_territorial_season:
    input:
        "data/processed/clean.parquet"
    output:
        "data/processed/subset_territorial_season.parquet"
    run:
        import geopandas as gpd
        gdf_clean = gpd.read_parquet(input[0])
        gdf_clean.loc[gdf_clean.index.get_level_values(1).month.isin(config['territorial_season_months'])].to_parquet(output[0])

# ts are shifted backward by 12h so nights span a single calendar date.
rule extract_cohorts_juvenile_nights:
    input:
        "data/processed/subset_territorial_season.parquet"
    output:
        "data/processed/subset_cohorts_juvenile_nights_timeshifted.parquet"
    run:
        import geopandas as gpd
        gdf_population = gpd.read_parquet(input[0])
        gdf_night = gdf_population.loc[gdf_population.cohort.isin(config['natal_territories']['cohorts'])]
        gdf_night = gdf_night.loc[gdf_night.age == 0.]
        gdf_night = gdf_night.loc[gdf_night.sunpos=='night']
        gdf_night.index = gdf_night.index.remove_unused_levels()
        gdf_night.index = gdf_night.index.set_levels(gdf_night.index.levels[1].shift(-12, 'h'), level=1)
        gdf_night.to_parquet(output[0])

rule clusterize_natal_territory_candidates:
    input:
        "data/processed/subset_cohorts_juvenile_nights_timeshifted.parquet"
    output:
        temp("data/processed/areas_natal_territory_candidates_raw.parquet")
    params:
        per_individual=1,
        dist_threshold=config['natal_territories']['clust_dist_threshold'],
        min_core_points=config['natal_territories']['clust_min_core_points'],
        min_days=config['natal_territories']['clust_min_days']
    script:
        "py/clusterize.py"

rule add_natal_territory_candidates_years:
    input:
        "data/processed/areas_natal_territory_candidates_raw.parquet",
        "data/processed/subset_cohorts_juvenile_nights_timeshifted.parquet"
    output:
        temp("data/processed/areas_natal_territory_candidates_years.parquet")
    run:
        import geopandas as gpd
        import pandas as pd
        years = gpd.read_parquet(input[1]).groupby(level=0).apply(lambda x: x.index.get_level_values(1).year[0]).rename('year')
        areas = gpd.read_parquet(input[0]).join(years, on='ind')
        cols = areas.columns.to_list()
        areas[cols[-1:] + cols[:-1]].to_parquet(output[0])

rule add_natal_territory_candidates_cooccurrence:
    input:
        "data/processed/areas_natal_territory_candidates_years.parquet",
        "data/processed/subset_territorial_season.parquet"
    output:
        temp("data/processed/areas_natal_territory_candidates_cooccur.parquet")
    params:
        min_days=config['natal_territories']['cooccur_min_days'],
        restrict_y=1
    script:
        "py/cooccurrence.py"

rule rank_natal_territory_candidates:
    input:
        "data/processed/areas_natal_territory_candidates_cooccur.parquet"
    output:
        "data/processed/areas_natal_territory_candidates.parquet"
    run:
        import geopandas as gpd
        import numpy as np
        areas = gpd.read_parquet(input[0])
        # because non-core samples from DBSCAN are removed, some clusters might fall below the ndays requirement for cooccurrence. Make sure cooccur is at least 1.
        areas.cooccur = np.clip(areas.cooccur, 1, None)
        # sort by descending ndays
        areas = areas.groupby(level=0, group_keys=False).apply(lambda x: x.sort_values('ndays', ascending=False))
        # rank based on co-occurrence
        ranks = areas.groupby(level=0, group_keys=False).apply(lambda x: x.cooccur.rank(method='first')).astype(int).rename('rank')
        areas.droplevel(1).set_index(ranks, append=True).sort_index().to_parquet(output[0])

rule package_natal_territory_candidates:
    input:
        "data/processed/areas_natal_territory_candidates.parquet"
    output:
        "data/geometries/natal_territories.gpkg"
    run:
        import geopandas as gpd
        gpd.read_parquet(input[0]).loc[(slice(None),1),:].droplevel(1).to_file(output[0], layer='natal_territories')

#
# derive movement metrics
#
rule estimate_daily_range:
    input:
        "data/processed/clean.parquet"
    output:
        "data/processed/movement_daily_range.parquet"
    params:
        daily_min_points=config['movement']['range_daily_min_points']
    script:
        "py/movement_range_estimate.py"
