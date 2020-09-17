import sys
sys.path.append("..")

import tracc
import pandas as pd
import numpy as np
import geopandas as gpd

from functools import reduce


dfo = pd.merge(
        pd.read_csv("test_data/boston/destination_employment_lehd.csv"),
        pd.read_csv("test_data/boston/destination_groceries_snap.csv"),
        how = "outer",
        left_on = "block_group_id",
        right_on = "GEOID"
)

dfo = tracc.supply(
    supply_df = dfo,
    columns = ["block_group_id","C000","snap"]
    )

out_data = []

for matrix in [["test_data/boston/transit_time_matrix_8am_30_06_2020.zip","3006"],["test_data/boston/transit_time_matrix_8am_29_02_2020.zip","2902"]]:

    dft = tracc.costs(
        pd.read_csv(matrix[0], compression='zip')
        )
    dft.data.time = dft.data.time / 60 # converting time from seconds to minutes
    dft.data.time = dft.data.time.round(1) # rounding to just one decimal place

    # removes existing intrazonal times
    dft.intrazonal(
        cost_column = "time",
        origin_column = "o_block",
        destination_column = "d_block",
        method = "constant",
        value = 0
    )
    dft.data = dft.data[dft.data["time"] > 0]

    dft.fill_missing_costs(
        where = "origin",
        cost_column = "time",
        origin_column = "o_block",
        destination_column = "d_block",
        spatial_file_path = "test_data/boston/block_group_poly.geojson",
        spatial_file_id = "GEOID"
    )

    dft.fill_missing_costs(
        where = "destination",
        cost_column = "time",
        origin_column = "o_block",
        destination_column = "d_block",
        spatial_file_path = "test_data/boston/block_group_poly.geojson",
        spatial_file_id = "GEOID"
    )

    dft.intrazonal(
        cost_column = "time",
        origin_column = "o_block",
        destination_column = "d_block",
        method = "radius",
        value = 1 / 0.0833333333,
        polygon_file = "test_data/boston/block_group_poly.geojson",
        polygon_id = "GEOID"
    )


    dft.impedence_calc(
        cost_column = "time",
        impedence_func = "cumulative",
        impedence_func_params = 45,
        output_col_name = "fCij_c45",
        prune_output = False
    )

    dft.impedence_calc(
        cost_column = "time",
        impedence_func = "cumulative",
        impedence_func_params = 30,
        output_col_name = "fCij_c30",
        prune_output = False
    )

    dft.impedence_calc(
        cost_column = "time",
        impedence_func = "cumulative",
        impedence_func_params = 60,
        output_col_name = "fCij_c60",
        prune_output = False
    )

    acc = tracc.accessibility(
            travelcosts_df = dft.data,
            supply_df = dfo.data,
            travelcosts_ids = ["o_block","d_block"],
            supply_ids = "block_group_id"
        )

    dfa = acc.potential(
            opportunity = "C000",
            impedence = "fCij_c45"
            )


    dfa.columns = ['o_block','A_C000_TfCijc45_' + matrix[1] + "_MD"]

    out_data.append(dfa)

    dfa = acc.potential(
            opportunity = "C000",
            impedence = "fCij_c60"
            )


    dfa.columns = ['o_block','A_C000_TfCijc60_' + matrix[1] + "_MD"]

    out_data.append(dfa)

    dfa = acc.potential(
            opportunity = "C000",
            impedence = "fCij_c30"
            )


    dfa.columns = ['o_block','A_C000_TfCijc30_' + matrix[1] + "_MD"]

    out_data.append(dfa)

    dfa = acc.mintravelcost(
            travelcost = "time",
            opportunity = "snap",
            min_n = 1
            )


    dfa.columns = ['A_snap_Tmintime1_' + matrix[1] + "_MD"]

    out_data.append(dfa)

    dfa = acc.mintravelcost(
            travelcost = "time",
            opportunity = "snap",
            min_n = 3
            )


    dfa.columns = ['A_snap_Tmintime3_' + matrix[1] + "_MD"]

    out_data.append(dfa)


dfoa = reduce(lambda x, y: pd.merge(x, y, on = 'o_block'), out_data)

print(dfoa.head())

dfoa.to_csv("test_data/boston/sample_accesibility_V5.csv", index = False)
