import tracc
import pandas as pd
import numpy as np


class costs:

    def __init__(self,
        travelcosts_df,
        columns = None
        ):

        """
        Inputs data and prunes columns if desired
        """

        if columns is not None:
            self.data = travelcosts_df[columns]

        else:
            self.data = travelcosts_df


    def intrazonal(self,
        cost_column,
        origin_column,
        destination_column,
        method = "constant",
        value = 0,
        polygon_file = None,
        polygon_id = None
        ):
        """
        Computes and updates intrazonal travel cost in a travel costs matrix. The output will include a travel cost between any origin or destination location in the matrix to itself.

        Parameters
        ----------
        cost_column : column name for travel costs

        origin_column : column name for origin IDs

        destinationn_column : column name for origin IDs

        method : "constant" applies a single @value to all intrazonal travel costs. "radius" applies a cost which is proportional to the radius of a circle with the same area as its input polygon

        value : parameters for the method

        polygon_file : file path to an input spatial polygon (e.g. geojson) if needed (it is for method = "radius")

        polygon_id : ID field for the polygon_file needed for joining to the cost matrix
        """

        # making sure ID columns are strings for a merge later on
        self.data[origin_column] = self.data[origin_column].astype(str)
        self.data[destination_column] = self.data[destination_column].astype(str)

        # getting set of unique locations in the dataset
        locations = list(self.data[origin_column].unique()) + list(self.data[destination_column].unique())
        locations = list(set(locations))

        if method == "constant":

            new_times = [value] * len(locations)

            df = pd.DataFrame(
                list(zip(locations, locations, new_times)),
                columns =[origin_column, destination_column, cost_column + "_i"])

        elif method == "radius":

            from tracc.spatial import radius

            # compute based on the equivilant radius of each polygon
            df = radius(polygon_file,polygon_id)
            df[origin_column] = df[polygon_id]
            df[destination_column] = df[polygon_id]
            del df[polygon_id]
            df[cost_column + "_i"] = value * df["radius"]
            del df["radius"]

        else:
            raise Exception("Method can only be 'constant' or 'radius'")

        df[origin_column] = df[origin_column].astype(str)
        df[destination_column] = df[destination_column].astype(str)

        # join in the newly created intrazonal travel times
        self.data = pd.merge(self.data, df,  how='outer', left_on=[origin_column, destination_column], right_on = [origin_column, destination_column])

        # replace the older intrazonal travel times
        self.data[cost_column] = np.where((self.data[cost_column + "_i"] >= 0),self.data[cost_column + "_i"],self.data[cost_column])

        del self.data[cost_column + "_i"]




    def fill_missing_costs(
        self,
        cost_column,
        origin_column,
        destination_column,
        spatial_file_path,
        spatial_file_id,
        where = "origin",
        weight_type = "Queen"
        ):
        """
        Completes an OD matrix by filling locations that were missing from the original matrix, based on a neighbourhood spatial weights matrix. For example if a origin zone has no travel costs, it presumes its travel costs to destinations are the average of the same costs of its neighbouring zones.
        """

        from tracc.spatial import area

        # get list of zones which are missing from the input costs table
        dfz = area(spatial_file_path, spatial_file_id)
        dfz[spatial_file_id] = dfz[spatial_file_id].astype(str)
        self.data[origin_column] = self.data[origin_column].astype(str)
        li1 = list(self.data[origin_column].unique())
        li2 = list(dfz[spatial_file_id].unique())
        missing = [x for x in li2 if x not in li1]
        del li1,li2

        if len(missing) == 0:
            return None

        if where == "origin":

            # get neighbours for each missing zone
            from tracc.spatial import get_neighbours
            neighbours = get_neighbours(spatial_file_path, "Queen", spatial_file_id)

            new_times = []

            # for each zone, compute average travel times to other zones based on neighbours
            for location in missing:

                locneigh = neighbours[location]

                temp = self.data[self.data[origin_column].isin(locneigh)]

                temp = pd.DataFrame(temp.groupby([destination_column], as_index=False)[cost_column].mean())

                temp[origin_column] = location

                new_times.append(temp)

            # combine the outputs, and concat to the input times
            new_times = pd.concat(new_times)
            self.data = pd.concat([self.data, new_times])

        elif where == "destination":

            # get neighbours for each missing zone
            from tracc.spatial import get_neighbours
            neighbours = get_neighbours(spatial_file_path, "Queen", spatial_file_id)

            new_times = []

            # for each zone, compute average travel times from other zones based on neighbours
            for location in missing:

                locneigh = neighbours[location]

                temp = self.data[self.data[destination_column].isin(locneigh)]

                temp = pd.DataFrame(temp.groupby([origin_column], as_index=False)[cost_column].mean())

                temp[destination_column] = location

                new_times.append(temp)

            # combine the outputs, and concat to the input times
            new_times = pd.concat(new_times)
            self.data = pd.concat([self.data, new_times])

        else:

            raise Exception("Input paramater @where should either be 'origin' or 'destination'")




    def generalized_cost(
        self,
        columns,
        coefficients,
        exponents = None,
        prune_output = True,
        output_cost_name = "GC"
        ):

        """
        Computes generalized costs
        """

        # need to add a column check warning, and make the intercept = 0 if none is provided

        # set all exponents as 1 if none are inputted
        if exponents is None:
            exponents = [1] * len(columns)

        # compute the generalized cost value
        self.data[output_cost_name] = coefficients[len(coefficients) - 1]
        i = 0
        while i < len(columns):
            self.data[output_cost_name] = self.data[output_cost_name] + coefficients[i] * self.data[columns[i]] ** exponents[i]
            i += 1

        # delete initital cost columns if desired
        if prune_output is True:
            for col in list(set(columns)):
                del self.data[col]


    def impedence_calc(
        self,
        cost_column,
        impedence_func,
        impedence_func_params,
        prune_output = True,
        output_col_name = "fCij"
        ):

        """
        Measures impdence given input of travel cost and selected impedence funciton and parameters

        # To Do: add in more impdence function options
        """

        if impedence_func == "cumulative":
            self.data[output_col_name] = self.data[cost_column].apply(tracc.decay.cumulative,args = (impedence_func_params,))

        elif impedence_func == "linear":
            self.data[output_col_name] = self.data[cost_column].apply(tracc.decay.linear,args = (impedence_func_params,))

        elif impedence_func == "neg_exp":
            self.data[output_col_name] = self.data[cost_column].apply(tracc.decay.neg_exp,args = (impedence_func_params,))

        else:
            raise Exception("Please select an appropriate decay function")

        if prune_output is True:
            del self.data[cost_column]


    def impedence_combine(self,
        columns,
        how = "product",
        output_col_name = "fCij",
        prune_output = True
        ):

        """
        If there are multiple impedences, and we want to combine them into a single impedence value. This is similar to genearlized cost.

        For example, if we have an impedence value for transit travel time, and we also want to remove any trips based on a fare criteria, it can be applied in this way.
        """

        if how == "product":
            self.data[output_col_name] = 1
            i = 0
            while i < len(columns):
                self.data[output_col_name] = self.data[output_col_name] * self.data[columns[i]]
                i += 1

        elif how == "sum":
            self.data[output_col_name] = 0
            i = 0
            while i < len(columns):
                self.data[output_col_name] = self.data[output_col_name] + self.data[columns[i]]
                i += 1

        else:
            raise Exception('the input @how must be one of "product" or "sum"')



    def max_impedence(self,
        columns,
        imp_col_name = "fCij"
        ):
        """
        Reduces the cost table to only include rows with the maximum impedence value for the set of input columns.

        For example, if there 3 transit trips from i to j, each with a different computed generalized_cost resulting from different route choices, this function will return the row with the one resulting in the greatest impedence value (i.e. lowest generalized cost)
        """

        self.data = self.data.groupby(columns)[imp_col_name].max().reset_index()



class supply:

    def __init__(self,
        supply_df,
        columns = None
        ):
        """
        intitializing can include pruning the dataset to a list of @column names
        """

        if columns is not None:
            self.data = supply_df[columns]

        else:
            self.data = supply_df



    def weight(self,
        columns,
        weights,
        weight_col_name = "Oj",
        prune_output = True
        ):
        """
        Creating a value based on a weighted linear combination other values. Can be used to weight by destinations by their desirability.

        Parameters
        ----------------
        columns : columns in which to input into the weights function

        weights : linear multipliers, the same length as the weights

        weight_col_name : output column name

        prune_output : if True, delete all input columns used in the weight function
        """

        if len(columns) != len(weights):
            raise Exception("Please make sure columns and weights are lists of the same length")

        if len(columns) < 2:
            raise Exception("Can only weight opportunities if 2 or more are inputted")

        if sum(weights) < 0.999 or sum(weights) > 1.001:
            print("WARNING: the inputted weights do not sum to 1.")


        self.data[weight_col_name] = 0
        i = 0
        while i < len(columns):
            self.data[weight_col_name] = self.data[weight_col_name] + weights[i] * self.data[columns[i]]
            i += 1

        if prune_output is True:
            for col in list(set(columns)):
                del self.data[col]



class demand:

    def __init__(self,
        demand_df,
        columns = None
        ):
        """
        intitializing can include pruning the dataset to a list of @column names
        """

        if columns is not None:
            self.data = demand_df[columns]

        else:
            self.data = demand_df


    def weight(self,
        columns,
        weights,
        weight_col_name = "Pi",
        prune_output = True
        ):
        """
        Creating a value based on a weighted linear combination other values. Can be used to weight by population groups by their propensity to travel to certain activity types.

        Parameters
        ----------------
        columns : columns in which to input into the weights function

        weights : linear multipliers, the same length as the weights

        weight_col_name : output column name

        prune_output : if True, delete all input columns used in the weight function
        """

        if len(columns) != len(weights):
            raise Exception("Please make sure columns and weights are lists of the same length")

        if len(columns) < 2:
            raise Exception("Can only weight opportunities if 2 or more are inputted")

        if sum(weights) < 0.999 or sum(weights) > 1.001:
            print("WARNING: the inputted weights do not sum to 1.")

        self.data[weight_col_name] = 0
        i = 0
        while i < len(columns):
            self.data[weight_col_name] = self.data[weight_col_name] + weights[i] * self.data[columns[i]]
            i += 1

        if prune_output is True:
            for col in list(set(columns)):
                del self.data[col]



class accessibility:

    def __init__(self,
        travelcosts_df,
        supply_df,
        demand_df = None,
        travelcosts_ids = ["origin_id","destination_id"],
        supply_ids = "destination_id",
        demand_ids = None
        ):
        """
        Parameters
        ----------
        travelcosts_df : a pandas dataframe containing travel costs from a set of locations (e.g. orignis) to another set of locations (e.g. destinations). Data should be in a long table format:

        origin_id | destination_id | travel_cost_1 | travel_cost_2 (optional) | etc (optional)

        supply_df : a pandas dataframe containing the number of opportunities (e.g. supply), relational to the destination IDs in travelcosts_df

        demand_df : a pandas dataframe containing the number of agents competiting for opportunities (e.g. demand), relational to the origin IDs in travelcosts_df. This is optional since several accessibility measures do not account for demand

        travelcosts_ids : a two item list of the column names for the origin and destination IDs in the travelcosts_df table

        supply_ids : a single variable string for the destination ID in the supply_df table

        demand_ids : a single variable string for the origin ID in the demand_df table. This is optional since several accessibility measures do not account for demand

        """

        self.travelcosts_ids = travelcosts_ids
        self.supply_ids = supply_ids
        self.demand_ids = demand_ids

        if demand_df is None and supply_df is None:
            raise Exception("Please input a supply_df or a demand_df")

        # setting ID columns to strings to aid merging
        travelcosts_df[travelcosts_ids[0]] =         travelcosts_df[travelcosts_ids[0]].astype(str)
        travelcosts_df[travelcosts_ids[1]] =         travelcosts_df[travelcosts_ids[1]].astype(str)

        # join supply data to the travel costs
        if supply_df is not None and demand_df is None:
            supply_df[supply_ids] = supply_df[supply_ids].astype(str)
            self.data = pd.merge(
                travelcosts_df,
                supply_df,
                left_on=travelcosts_ids[1],
                right_on=self.supply_ids,
                how = 'left'
            )

        # join demand data as well, if inputted
        elif demand_df is not None and supply_df is None:
            demand_df[demand_ids] = demand_df[demand_ids].astype(str)
            self.data = pd.merge(
                travelcosts_df,
                demand_df,
                left_on=travelcosts_ids[0],
                right_on=self.demand_ids,
                how = 'left'
            )

        else:
            supply_df[supply_ids] = supply_df[supply_ids].astype(str)
            demand_df[demand_ids] = demand_df[demand_ids].astype(str)
            self.data = pd.merge(
                travelcosts_df,
                supply_df,
                left_on=travelcosts_ids[1],
                right_on=self.supply_ids,
                how = 'left'
            )
            self.data = pd.merge(
                self.data,
                demand_df,
                left_on=travelcosts_ids[0],
                right_on=self.demand_ids,
                how = 'left'
            )


    def potential(self, opportunity, impedence, output_col_name = None):
        """
        Measures potential accessibility to destinations

        Parameters
        ----------
        opportunity : a string indicating the column name for which opportunity we are measuring access to (e.g. jobs, grocery stores, etc.). This column should be in the supply_df dataframe

        impedence : column from the travel costs object to weight opportunities by

        output_col_name : a string for the column name of the output accessibility measure


        Output
        ----------
        A pandas dataframe with the first column with the IDs of the origin point (self.travelcosts_ids[0]), and the second column accessibility measures based on the input parameters.

        """

        # set the output name for the accessibility measure
        if output_col_name is None:
            A_col_name = "A_" + opportunity + "_" + impedence
        else:
            A_col_name = output_col_name

        # multiply the opportunity by the impedence
        self.data[A_col_name] = self.data[opportunity] * self.data[impedence]

        # sum by the origin locations
        Ai = self.data.groupby(self.travelcosts_ids[0])[[A_col_name]].sum().reset_index()

        del self.data[A_col_name]

        return Ai




    def passive(self, population, impedence, output_col_name = None):

        """
        Measures passive accessibility to destinations

        Parameters
        ----------
        population : a string indicating the column name for which population we are measuring access to (e.g. overall population, employed population, etc.). This column should be in the demand_df dataframe

        impedence : column from the travel costs object to weight opportunities by

        output_col_name : a string for the column name of the output accessibility measure


        Output
        ----------
        A pandas dataframe with the first column with the IDs of the origin point (self.travelcosts_ids[0]), and the second column accessibility measures based on the input parameters.

        """

        # set the output name for the accessibility measure
        if output_col_name is None:
            A_col_name = "A_" + population + "_" + impedence
        else:
            A_col_name = output_col_name

        # multiply the opportunity by the impedence
        self.data[A_col_name] = self.data[population] * self.data[impedence]

        # sum by the origin locations
        Ai = self.data.groupby(self.travelcosts_ids[1])[[A_col_name]].sum().reset_index()

        del self.data[A_col_name]

        return Ai




    def mintravelcost(self, travelcost, opportunity, min_n,  output_col_name = None):
        """
        Parameters
        ----------
        opportunity : a string indicating the column name for which opportunity we are measuring access to (e.g. jobs, grocery stores, etc.). This column should be in the supply_df dataframe

        travelcost : a string indicating the column name for which travel cost shall be used (e.g. travel time, monetary cost, etc.). This column should be in the travelcosts_df dataframe

        min_n : an int indicating the number of desired reachable opportunities (e.g. 1 library, 3 grocery stores, 10k jobs, etc.)

        output_col_name : a string for the column name of the output accessibility measure



        Output
        ---------
        A pandas dataframe with the first column with the IDs of the origin point (self.travelcosts_ids[0]), and the second column are the accessibility measures based on the input parameters.
        """

        # set the output name for the accessibility measure
        if output_col_name is None:
            A_col_name = "A_mintravelcost_" + str(travelcost) + "_" + str(opportunity) + "_" +  str(min_n)
        else:
            A_col_name = output_col_name

        # internal function of returning the min travel time for n opportunities
        def get_min(df, tc, o, n):
            df = df.sort_values(by=[tc], ascending=True)
            df["cumsum"] = df[o].cumsum()
            df = df[df["cumsum"] >= n]
            return df[travelcost].min()

        # generating the accessibility measure
        out =  pd.DataFrame(self.data.groupby(self.travelcosts_ids[0]).apply(get_min, tc = travelcost, o = opportunity, n = min_n))

        # setting the column name of the output
        out.columns = [A_col_name]

        return out



class summary:
    """
    Computing various summary statistics of accessibility, usually with respect to different population groups

    Some of these can be used to assess distributions and equity of transport networks.
    """

    def __init__(
        self,
        accessibility_df,
        summary_vars,
        accessibility_id = "id",
        summary_vars_id = "id"
        ):

        # join the data
        self.data = pd.merge(
            accessibility_df,
            summary_vars,
            left_on=accessibility_id,
            right_on=summary_vars_id,
            how = 'left'
        )

    def weighted_mean(self, access_var, group_var):

        return tracc.statistics.weighted_mean(self.data, access_var, group_var)

    def weighted_var(self, access_var, group_var):

        return tracc.statistics.weighted_var(self.data, access_var, group_var)

    def weighted_sd(self, access_var, group_var):

        return tracc.statistics.weighted_sd(self.data, access_var, group_var)

    def weighted_CV(self, access_var, group_var):

        return tracc.statistics.weighted_CV(self.data, access_var, group_var)

    def weighted_Gini(self, access_var, group_var):

        return tracc.statistics.weighted_Gini(self.data, access_var, group_var)

    def quantiles(self, access_var, group_vars, nbins = 10, result = "percent"):

        # assign each observation a bin, based on nbins
        dfq = pd.DataFrame( tracc.statistics.weighted_qcut(self.data[access_var], self.data[group_vars[0]], nbins))

        # create a specific name for the quantile column
        q_col_name = 'q' + str(nbins) + "_" + (group_vars[0])
        dfq.columns = [q_col_name]
        self.data = self.data.join(dfq, how='outer')

        # group by each bin, susmmarize
        dfq = self.data.groupby([q_col_name])[group_vars].sum()

        # return as counts or percent
        if result == "count":
            return dfq
        elif result == "percent":
            for var in group_vars:
                dfq[var] = dfq[var] / dfq[var].sum()
            return dfq
