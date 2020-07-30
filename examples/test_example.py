import tracc
import pandas as pd
import numpy as np
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
    
    print(dfa)
    
    dfa.columns = ['o_block','A_C000_fCij_c45_' + matrix[1]]
    
    out_data.append(dfa)
    
    dfa = acc.potential(
            opportunity = "C000",
            impedence = "fCij_c60"
            )
    
    print(dfa)
    
    dfa.columns = ['o_block','A_C000_fCij_c60_' + matrix[1]]
    
    out_data.append(dfa)
    
    dfa = acc.potential(
            opportunity = "C000",
            impedence = "fCij_c30"
            )
    
    print(dfa)
    
    dfa.columns = ['o_block','A_C000_fCij_c30_' + matrix[1]]
    
    out_data.append(dfa)

    dfa = acc.mintravelcost(
            travelcost = "time",
            opportunity = "snap",
            min_n = 1
            )
    
    print(dfa)
    
    dfa.columns = ['A_mintravelcost_time_snap_1_' + matrix[1]]
    
    out_data.append(dfa)
    
    dfa = acc.mintravelcost(
            travelcost = "time",
            opportunity = "snap",
            min_n = 3
            )
    
    print(dfa)
    
    dfa.columns = ['A_mintravelcost_time_snap_3_' + matrix[1]]
    
    out_data.append(dfa)
    
    
dfoa = reduce(lambda x, y: pd.merge(x, y, on = 'o_block'), out_data)

print(dfoa.head())

dfoa.to_csv("test_data/boston/sample_accesibility_V2.csv")