# TO DO:
# add comments
# floating catchment area / competitive functions
# test all flow for larger dataset


import sys
sys.path.append("..")


import tracc
import pandas as pd
import numpy as np

import time


# # dfc = pd.read_csv("tests/Boston8AM.csv")
# # dfo = pd.read_csv("tests/destination_employment_lehd.csv")
#
# #
# # origins = list(dfc["OriginName"].unique())
# #
# # dfo = tracc.supply(dfo)
# #
# # dfo.weight(
# #     columns = ["CA01","CA02"],
# #     weights = [0.8,0.2],
# #     prune_output = False
# #     )
#
# #
# # print(dfc[dfc["OriginName"].isin([250250104051,250250701013,250259801011,250250701016])])
#
# # chunks = [origins[x:x+100] for x in range(0, len(origins), 100)]
# #
# #
# # test = []
# #
# # i = 0
# # for chunk in chunks:
# #
# #     dfcs = dfc[dfc["OriginName"].isin(chunk)]
# #
# #     dfcs = tracc.costs(dfcs,["OriginName","DestName","Total_Time","Total_Dist"])
# #
# #     # dfcs.generalized_cost(["Total_Time","Total_Dist"],[1,0.05,1],[2,1], prune_output = True)
# #
# #     dfcs.impedence_calc("Total_Time","cumulative",30)
# #
# #     acc = tracc.accessibility(
# #         travelcosts_df = dfcs.data,
# #         supply_df = dfo,
# #         travelcosts_ids = ["OriginName","DestName"],
# #         supply_ids = "block_group_id"
# #     )
# #
# #     test.append(acc.gravity("C000", "fCij"))
# #
# #
# #     i += 1
# #     print(i)
# #
# #
# #
# # test = pd.concat(test)
# #
# # print(test)
# #
# # test.to_csv("tests/output_test_acc.csv")
# #
# # print(dfo["C000"].sum())
#
#
#
#
#
# # print(tracc.decay.linear(30,[[1,1,1,1],[2,1,1,1,1]], upper_bound = 1))
#
#
# dfc = pd.read_csv("tests/test_data_2_costs.csv")
# dfo = pd.read_csv("tests/test_data_1_supply.csv")
# dfd = pd.read_csv("tests/test_data_1_demand.csv")
#
#
# dfcs = tracc.costs(dfc,["o","d","c","f"])
#
# dfcs.impedence_calc("c","linear",-0.02,output_col_name = "fCij_time")
# dfcs.impedence_calc("f","cumulative",1,output_col_name = "fCij_fare")
#
# print(dfcs.data)
#
# dfcs.impedence_combine(["fCij_time","fCij_fare"])
#
# dfcs.max_impedence(["o","d"], "fCij")
#
# print(dfcs.data)
#
#
# #
#
# acc = tracc.accessibility(
#         travelcosts_df = dfcs.data,
#         supply_df = dfo,
#         demand_df = dfd,
#         travelcosts_ids = ["o","d"],
#         supply_ids = "dd",
#         demand_ids = "o"
#     )
#
# g = acc.potential("oj","fCij")
#
# print(g)
#
# print(dfd)

# dfa = pd.read_csv("tests/output_test_acc.csv")
#
# print(dfa)
#
# dfp = pd.read_csv("tests/demographics.csv")
# dfp.fillna(0, inplace=True)
#
# print(dfp)
#
# summary = tracc.summary(
#     accessibility_df = dfa,
#     summary_vars = dfp,
#     accessibility_id = "OriginName",
#     summary_vars_id = "geoid")
# #
# print(summary.data)
#
# # for p in ["B01001_001E", "C000"]:
# #
# #     print(p)
# #     print(summary.weighted_mean("jobs_30min_auto",p))
# #     print(summary.weighted_Gini("jobs_30min_auto",p))
#
# test = summary.quantiles(
#     access_var ="jobs_30min_auto",
#     group_vars = ["B01001_001E", "C000"],
#     nbins = 10)
#
# print(test)
#
# #
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


dft = tracc.costs(
    pd.read_csv("test_data/boston/transit_time_matrix_8am_30_06_2020.zip", compression='zip')
    )
dft.data.time = dft.data.time / 60 # converting time from seconds to minutes
dft.data.time = dft.data.time.round(1) # rounding to just one decimal place

dft.intrazonal(
    cost_column = "time",
    origin_column = "o_block",
    destination_column = "d_block",
    method = "constant",
    value = 0
)

print(dft.data)

# dft.data = dft.data[dft.data["time"] > 0]
#
# dft.fill_missing_costs(
#     where = "origin",
#     cost_column = "time",
#     origin_column = "o_block",
#     destination_column = "d_block",
#     spatial_file_path = "test_data/boston/block_group_poly.geojson",
#     spatial_file_id = "GEOID"
# )
#
# dft.fill_missing_costs(
#     where = "destination",
#     cost_column = "time",
#     origin_column = "o_block",
#     destination_column = "d_block",
#     spatial_file_path = "test_data/boston/block_group_poly.geojson",
#     spatial_file_id = "GEOID"
# )

# print(dft.data[dft.data["o_block"] == '250277261005'])


# print(dft.data[dft.data["o_block"] == '250173632022'])
# print(dft.data[dft.data["o_block"] == '250173231002'])
# print(dft.data[dft.data["o_block"] == '250173632013'])

start_time = time.time()


dft.intrazonal(
    cost_column = "time",
    origin_column = "o_block",
    destination_column = "d_block",
    method = "radius",
    value = 1 / 0.0833333333,
    polygon_file = "test_data/boston/block_group_poly.geojson",
    polygon_id = "GEOID"
)

print(time.time() - start_time)


# print(dft.data)


# dft.impedence_calc(
#     cost_column = "time",
#     impedence_func = "cumulative",
#     impedence_func_params = 45,
#     output_col_name = "fCij_c45",
#     prune_output = False
# )
#
#
# acc = tracc.accessibility(
#         travelcosts_df = dft.data,
#         supply_df = dfo.data,
#         travelcosts_ids = ["o_block","d_block"],
#         supply_ids = "block_group_id"
#     )
#
# dfa = acc.potential(
#         opportunity = "C000",
#         impedence = "fCij_c45"
#         )
#
# print(dfa)
# dfa.to_csv("test_data_nofill.csv")
