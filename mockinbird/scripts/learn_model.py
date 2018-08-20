import argparse
import json
import os
import pickle
import logging


import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.misc import logsumexp
import pandas as pd

from mockinbird.utils.fit_betabinom import fit_betabinom_ab, sample_betabinom_ab


logger = logging.getLogger()


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('site_statistics')
    parser.add_argument('max_mixture_components', type=int, default=3)
    parser.add_argument('bam_stat_json')
    parser.add_argument('out_dir')
    parser.add_argument('--no_global_fit', action='store_true')
    parser.add_argument('--n_iterations', default=100, type=int)
    parser.add_argument('--dump_data', action='store_true')
    parser.add_argument('--max_coverage', type=float, default=np.inf)
    return parser


def plot_fit(x, x_weights, w, a, title, file_name, out_dir, max_x=None):

    def geom_distr(k, w, a):
        return np.sum(w * (a / (1 + a) ** (k + 1)))

    if max_x:
        bins = max_x
    else:
        max_x = int(x.max())
        bins = 100

    fit = [geom_distr(k, w, a) for k in range(max_x + 1)]
    fit = np.array(fit)
    fig, ax = plt.subplots()

    subset_weights = x_weights[x <= max_x]
    eta = subset_weights.sum() / x_weights.sum()
    values, edges = np.histogram(x[x <= max_x], bins=bins, weights=subset_weights,
                                 density=True)
    height = values * eta

    ax.bar(x=edges[:-1], width=np.diff(edges), height=height, align='edge', alpha=0.5)
    ax.plot(fit, color='red')
    ax.set_yscale("log", nonposy='clip')
    if max_x:
        ax.set_xlim((0, max_x))
    fig.suptitle(title, fontsize=18)
    fname = os.path.join(out_dir, file_name)
    plt.savefig(fname)


def geom_mm_fit_fast(x, x_counts, n_components, n_iter=100):
    # unfortunately this time I parameterized with p = (1-q). Shame on me.
    weights = [0.01, 0.99]
    x_mean = np.mean(x)
    q_mean = x_mean / (1 + x_mean)
    x_max = np.max(x)
    q_max = x_max / (1 + x_max)

    p_mean = 1 - q_mean
    p_max = 1 - q_max

    p_initial = [p_mean, p_mean + 0.05 * p_mean]
    for m in range(2, n_components):
        p_initial.append(p_max + 0.01 * m * p_max)
        weights = np.hstack((weights, [0.0001]))
        weights /= np.sum(weights)

    ps, weights = fast_geom_mm_fit(x, p_initial, weights, n_iter=n_iter, weights=x_counts)
    w_m = np.array(weights)
    EPSILON = 1e-30
    g_m = np.array([1 - p + p*EPSILON for p in ps])
    g_m = (1 - g_m) / g_m
    return w_m, g_m


def main():
    parser = create_parser()
    args = parser.parse_args()

    ch = logging.StreamHandler()
    logger.addHandler(ch)
    logger.setLevel(logging.INFO)

    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    with open(args.bam_stat_json) as bam_json:
        bam_stat = json.load(bam_json)

    mock_table = pd.read_table(args.site_statistics,
                               usecols=['n_mock', 'k_mock', 'count'])
    mock_table = mock_table.query('n_mock < %s' % args.max_coverage)
    mock_table = mock_table.groupby(['n_mock', 'k_mock']).agg({'count': sum}).reset_index()

    k_collapsed = mock_table.groupby('k_mock').agg({'count': sum}).reset_index()
    k_vals = k_collapsed.k_mock
    k_counts = k_collapsed['count']

    n_collapsed = mock_table.groupby('n_mock').agg({'count': sum}).reset_index()
    n_vals = n_collapsed.n_mock
    n_counts = n_collapsed['count']

    if args.dump_data:
        n_file = os.path.join(args.out_dir, 'n_dump.tab')
        pd.DataFrame({
            'n': n_vals,
            'count': n_counts,
        }).to_csv(n_file, sep='\t', index=False)

        k_file = os.path.join(args.out_dir, 'k_dump.tab')
        pd.DataFrame({
            'k': k_vals,
            'count': k_counts,
        }).to_csv(k_file, sep='\t', index=False)

        kn_file = os.path.join(args.out_dir, 'kn_dump.tab')
        mock_table.to_csv(kn_file, sep='\t', index=False)

    coverage = bam_stat['total_coverage']

    parameters = {}
    parameters['coverage'] = coverage
    parameters['n_mixture_components'] = args.max_mixture_components

    # fit p(k|z=0)
    logger.info('Fitting p(k)')
    pw_m, g_m = geom_mm_fit_fast(k_vals, k_counts, args.max_mixture_components, args.n_iterations)
    plot_fit(k_vals, k_counts, pw_m, g_m, 'p(k)', 'pk_fit', args.out_dir, max_x=25)
    plot_fit(k_vals, k_counts, pw_m, g_m, 'p(k)', 'pk_fit_long', args.out_dir, max_x=500)
    pg_m = coverage * g_m
    parameters['pk_params'] = pw_m, pg_m

    # fit p(n)
    logger.info('Fitting p(n)')
    w, g = geom_mm_fit_fast(n_vals, n_counts, args.max_mixture_components, args.n_iterations)
    plot_fit(n_vals, n_counts, w, g, 'p(n)', 'pn_fit', args.out_dir, max_x=400)
    plot_fit(n_vals, n_counts, w, g, 'p(n)', 'pn_fit_long', args.out_dir, max_x=5000)
    pg = coverage * g
    parameters['pn_params'] = w, pg
    small_n_params = {}

    logger.info('Fitting p(k|n)')
    INDIVIDUAL_FIT_N_THRESH = 5
    if not args.no_global_fit:
        weights = np.array(mock_table['count'], dtype=float)
        k_vals = np.array(mock_table.k_mock, dtype=float)
        n_vals = np.array(mock_table.n_mock, dtype=float)
        mask = (n_vals > INDIVIDUAL_FIT_N_THRESH)
        assert len(k_vals) == len(n_vals)
        assert len(weights) == len(k_vals)
        alpha, beta = fit_betabinom_ab(n_vals[mask], k_vals[mask], weights=weights[mask])

        # fit small n separately
        mask = (n_vals > 0) & (n_vals <= INDIVIDUAL_FIT_N_THRESH)
        for n, w_df in mock_table.loc[mask, :].groupby('n_mock'):
            n_vec = np.ones(len(w_df)) * n
            alpha_idv, beta_idv = fit_betabinom_ab(n_vec, w_df.k_mock, weights=w_df['count'])
            small_n_params[n] = (alpha_idv, beta_idv)

    else:
        alpha_estim = []
        beta_estim = []
        for n in 3, 4, 5:
            subtable_n = mock_table[mock_table.n_mock == n]
            weights = subtable_n['count']
            k_vals = subtable_n.k_mock
            n_vals = subtable_n.n_n_mock  # all == n
            a, b = fit_betabinom_ab(n_vals, k_vals, weights)
            alpha_estim.append(a)
            beta_estim.append(b)
        alpha = np.mean(alpha_estim)
        beta = np.mean(beta_estim)

    parameters['pkn_params'] = alpha, beta
    parameters['pkn_idv_params'] = small_n_params

    model_file = os.path.join(args.out_dir, 'model.pkl')
    with open(model_file, 'wb') as pkl:
        pickle.dump(parameters, pkl)


def fast_geom_mm_fit(x, p_init, pi_init, n_iter=250, weights=None):
    if weights is None:
        N = len(x)
    else:
        N = sum(weights)
    p = np.array(p_init)
    pi = np.array(pi_init)
    assert len(p) == len(pi)

    # calculate the responsibility for each distinct x value just once
    # and weight by the number of occurences
    x_agg = np.bincount(x + 1, weights)[1:]
    x_new = np.arange(len(x_agg)) + 1


    EPSILON=1e-30

    n, d = len(x_new), len(p)
    r_matrix = np.zeros((d, n))
    for i in range(n_iter):
        # calculate the responsibilities
        for k in range(d):
            print('pi:', pi, 'p:', p)
            r_matrix[k, :] = np.log(pi[k]) + (x_new - 1) * np.log((1 - p[k] + p[k]*EPSILON)) + np.log(p[k])
            #r_matrix[k, :] = pi[k] * (1 - p[k]) ** (x_new - 1) * p[k]
        col_logexp_sum = logsumexp(r_matrix, axis=0)
        r_matrix = np.exp(r_matrix - col_logexp_sum)

        # optimize the parameters
        for k in range(d):
            Nk = np.sum(r_matrix[k, :] * x_agg)
            p[k] = Nk / np.dot(x_agg * x_new, r_matrix[k, :])
            pi[k] = Nk / (N)

    assert np.isclose(np.sum(pi), 1)
    return p, pi


if __name__ == '__main__':
    main()
