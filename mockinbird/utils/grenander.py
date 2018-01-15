import numpy as np

from rpy2.robjects.packages import importr, STAP
from rpy2 import robjects as ro
import rpy2.robjects.numpy2ri as rpyn


PVAL_ECDF_STR = '''
options(warn=-1)

ecdf_pval <- function (x)
{
    # compute empirical CDF as usual
    x = sort(x)
    n = length(x)
    if (n < 1)
        stop("'x' must have 1 or more non-missing values")
    vals = sort(unique(x))
    F.raw = cumsum(tabulate(match(x, vals)))/n

    # if necessary add an atom at 1 to make it a proper CDF
    if (vals[length(vals)] != 1)
    {
       F.raw = c(F.raw, 1)
       vals = c(vals, 1)
    }

    # if necessary also add an atom at 0 with weight zero to get support [0,1]
    if (vals[1] != 0)
    {
       F.raw = c(0, F.raw)
       vals = c(0, vals)
    }

    rval <- approxfun(vals, F.raw,
        method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
    class(rval) = c("ecdf", "stepfun", class(rval))
    attr(rval, "call") <- sys.call()
    rval
}

library(fdrtool)
plot_ecdf_lcm <- function(pvals, outfile, title) {
    x_uniq <- sort(unique(pvals))
    pval_ecdf <- ecdf(pvals)
    y <- pval_ecdf(x_uniq)
    res <- fdrtool::gcmlcm(x_uniq, y, type="lcm")
    pdf(outfile)
    plot(pval_ecdf, main=title, xlab='p-value', ylab='cumul. density', lwd=2)
    lines(res$x.knots, res$y.knots, col='red', lwd=2)
    dev.off()
}
'''

ecdf_pkg = STAP(PVAL_ECDF_STR, 'ecdf_pkg')


def pval_grenander_fit(pvals):
    fdrtool = importr('fdrtool')
    pval_vec = ro.FloatVector(pvals)
    pval_ecdf = ecdf_pkg.ecdf_pval(pval_vec)
    gre_fit = fdrtool.grenander(pval_ecdf, type='decreasing')
    x_knots = rpyn.ri2py(gre_fit.rx2('x.knots'))
    f_knots = rpyn.ri2py(gre_fit.rx2('f.knots'))
    if len(f_knots) == 0:
        x_knots = np.array([0, 1])
        f_knots = np.array([1, 1])
    assert len(f_knots) == len(x_knots)
    return x_knots, f_knots


PVALUE_EPSILON = 1e-12


def discrete_pval_grenander_fit(pvals, pval_freqs, eta0=None, full_output=False):

    x_ecdf = pvals
    y_ecdf = pval_freqs.cumsum() / len(pval_freqs)

    x_ecdf = np.append(0, x_ecdf)
    y_ecdf = np.append(0, y_ecdf)

    if eta0 is not None:
        y_ecdf = np.minimum(y_ecdf, 1 - eta0 * (1 - x_ecdf))
        y_ecdf = np.maximum(y_ecdf, eta0 * x_ecdf)

    # calculate the least concave majorant (lcm)
    # naive (?) implemenation:
    # starting from (0,0) in ecdf, draw a line to the point p with the steepest slope connection
    # repeat from point p, to p'' accoring to above criteria. Stop when (1,1) in ecdf is reached

    x_knots = [x_ecdf[0]]
    y_knots = [y_ecdf[0]]

    pval_slopes = []

    i = 0
    max_i = len(x_ecdf)
    while i < max_i - 1:
        x_cur = x_ecdf[i]
        y_cur = y_ecdf[i]
        delta_x = x_ecdf[i + 1:] - x_cur
        delta_y = y_ecdf[i + 1:] - y_cur

        slopes = delta_y / delta_x
        max_slope_index = np.argmax(slopes)
        max_slope = slopes[max_slope_index]

        jump_distance = 1 + max_slope_index
        pval_slopes.extend([max_slope] * jump_distance)
        i += jump_distance

        x_knots.append(x_ecdf[i])
        y_knots.append(y_ecdf[i])

    grenander_slopes = np.array(pval_slopes)

    if full_output:
        return grenander_slopes, x_ecdf, y_ecdf, np.array(x_knots), np.array(y_knots)
    else:
        return grenander_slopes


def plot_cumul_density(pvals, plot_basename, title):
    pval_vec = ro.FloatVector(pvals)
    ecdf_pkg.plot_ecdf_lcm(pval_vec, plot_basename + '.pdf', title)
