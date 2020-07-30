import math
import numpy as np
import pandas as pd

class decay:

    def __init__(self):
        None 

    def cumulative(cost, threshold):

        """
        Parameters
        ----------
        cost | an integer reperesenting the travel cost between two locations

        threshold | an integer representing a cutoff point

        Output
        ----------
        returns 1 if the cost is less than or equal to the threshold, 0 if not
        """

        if cost <= threshold:
            return 1
        else:
            return 0

    def linear(cost, params, lower_bound = 0, upper_bound = 1):

        """
        A linear combination of inpupt paramaters as follows

        params = [[x, y, ..., z],[p, q, r, ..., s]]

        returns f(C) = p + q C ^ x + r C ^ y + ... +

        if params is input as a single value, m

        returns f(C) = m C + 1
        """

        if isinstance(params, list):

            f = [params[1][0]]
            i = 0
            while i < len(params[0]):
                f.append(params[1][i + 1] * (cost ** params[0][i]))
                i += 1
            f = sum(f)

        else:
            intercept = 1
            slope = float(params)
            f = intercept + cost * slope

        if f < lower_bound:
            f = lower_bound
        if f > upper_bound:
            f = upper_bound
        return f


    def exponential(cost, params, lower_bound = 0, upper_bound = 1):

        """
        returns f(C) = e ^ (b C ^ a),

        params can be = [b, a, e], or = [b, a], or b

        if e is not input, it is assumed to be the constant approximately equal to 2.71828

        if a is not input, it is assumed to be 1, i.e. a simple negative exponential function,

        setting a = 2, creates a modifeid Gaussian
        """

        if isinstance(params, list):
            if len(params) == 1:
                f = math.exp(params[0] * cost)

            elif len(params) == 2:
                f = math.exp(params[0] * cost ** params[1])
            else:
                f = params[2] ** (params[0] * cost ** params[1])

        else:
            f = math.exp(params * cost ** 2)

        if f < lower_bound:
            f = lower_bound
        if f > upper_bound:
            f = upper_bound
        return f


    def power(cost, params, lower_bound = 0, upper_bound = 1):

        """
        returns f(C) = p (q C + r) ^ b + s

        params = [b,p,q,r,s], or params = b

        if params = b is inputted, then returns f(C) = C ^ b
        """

        if isinstance(params, list):
            if len(params) == 1:
                f = cost ** params[0]
            else:
                f = params[1] * ((params[2] * cost + params[3]) ** params[0]) + params[4]
        else:
            f = cost ** params

        if f < lower_bound:
            f = lower_bound
        if f > upper_bound:
            f = upper_bound
        return f



class statistics:

    def __init__(self):
        None

    def weighted_mean(df, x, w):

        return np.average(df[x],weights=df[w])

    def weighted_var(df, x, w):

        wm = np.average(df[x],weights=df[w])

        return np.average((df[x] - wm)**2, weights=df[w])

    def weighted_sd(df, x, w):

        wm = np.average(df[x],weights=df[w])

        return math.sqrt(np.average((df[x] - wm)**2, weights=df[w]))

    def weighted_CV(df, x, w):

        wm = np.average(df[x],weights=df[w])

        sd = math.sqrt(np.average((df[x] - wm)**2, weights=df[w]))

        return sd / wm

    def weighted_Gini(df, x, w):

        # borrowed from https://stackoverflow.com/questions/48999542/more-efficient-weighted-gini-coefficient-in-python

        # The rest of the code requires numpy arrays.
        x = np.asarray(df[x])
        w = np.asarray(df[w])
        sorted_indices = np.argsort(x)
        sorted_x = x[sorted_indices]
        sorted_w = w[sorted_indices]
        # Force float dtype to avoid overflows
        cumw = np.cumsum(sorted_w, dtype=float)
        cumxw = np.cumsum(sorted_x * sorted_w, dtype=float)
        return (np.sum(cumxw[1:] * cumw[:-1] - cumxw[:-1] * cumw[1:]) / (cumxw[-1] * cumw[-1]))


    def weighted_qcut(values, weights, q, **kwargs):
        'Return weighted quantile cuts from a given series, values.'
        quantiles = np.linspace(0, 1, q + 1)
        order = weights.iloc[values.argsort()].cumsum()
        bins = pd.cut(order / order.iloc[-1], quantiles, labels=False, **kwargs)
        return bins.sort_index()
