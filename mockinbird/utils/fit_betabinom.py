from functools import lru_cache
import warnings
import collections

import numpy as np

from rpy2 import robjects as ro
from rpy2.robjects.packages import STAP
import rpy2.robjects.numpy2ri as rpyn

from scipy.special import betaln


FIT_RFUN_STR = """
options(warn=-1)
library(VGAM)

fit_betabinom <- function(n, k) {
    fit <- vglm(cbind(k, n-k) ~ 1, betabinomialff)
    return(coef(fit, matrix = TRUE))
}

fit_betabinom_w <- function(n, k, w) {
    fit <- vglm(cbind(k, n-k) ~ 1, betabinomialff, weights=w)
    return(coef(fit, matrix = TRUE))
}

eval_betabinom <- function(n, k, a, b) {
    p <- dbetabinom.ab(k, size=n, shape1=a, shape2=b)
    return(p)
}

sample_betabinom <- function(n, N, a, b) {
    sample <- rbetabinom.ab(N, size=n, shape1=a, shape2=b) 
}
"""

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    bbfit = STAP(FIT_RFUN_STR, 'bbfit')


def fit_betabinom_ab(n, k, weights=None):
    assert np.all((k >= 0) & (k <= n))
    n_r = ro.FloatVector(n)
    k_r = ro.FloatVector(k)
    if weights is not None:
        assert len(weights) == len(k)
        assert len(weights) == len(n)
        weights_r = ro.FloatVector(weights)
        result_r = bbfit.fit_betabinom_w(n_r, k_r, weights_r)
    else:
        result_r = bbfit.fit_betabinom(n_r, k_r)
    result = rpyn.ri2py(result_r)
    log_a, log_b = result.flatten()
    return np.exp(log_a), np.exp(log_b)


def sample_betabinom_ab(n, N, a, b):
    result_r = bbfit.sample_betabinom(n, N, a, b)
    result = rpyn.ri2py(result_r)
    return result


def lp_kn_wrapper(a, b, idv_fit_ab):
    def lp_kn(k, n):

        if not isinstance(n, collections.Iterable):
            n = [n]

        def map_fun(n):
            return idv_fit_ab.get(n, (a, b))

        a_list, b_list = list(zip(*list(map(map_fun, n))))
        alpha = np.array(a_list)
        beta = np.array(b_list)

        lp = (
            - betaln(1 + n - k, 1 + k)
            - np.log(n + 1)
            + betaln(alpha + k, beta + n - k)
            - betaln(alpha, beta)
        )
        assert np.all(lp <= 0)
        return lp
    return lp_kn
