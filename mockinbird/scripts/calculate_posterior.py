import argparse
import os
import sys
import json
from functools import lru_cache
import pickle
import math

import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
from scipy.special import binom

from mockinbird.utils.grenander import discrete_pval_grenander_fit
from mockinbird.utils.fit_betabinom import p_kn_wrapper


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('factor_mock_table')
    parser.add_argument('factor_mock_statistics')
    parser.add_argument('mock_model')
    parser.add_argument('--max_iter', type=int, default=100)
    parser.add_argument('outfile')
    parser.add_argument('--plot_dir')
    parser.add_argument('--max_k_mock', type=int, default=10)
    parser.add_argument('--estimate_null_fraction', action='store_true')
    parser.add_argument('--null_fraction', type=float)
    parser.add_argument('--bam_statistics_json')
    parser.add_argument('--posterior_table')
    parser.add_argument('--scale_pvalues', action='store_true',
                        help='scaling pvalues such that the highest p-value is 1')
    parser.add_argument('--debug', action='store_true')
    return parser


def pk_wrapper(N, g_m, w_m):
    @lru_cache(maxsize=2**15)
    def pk_inner(k):
        return np.sum(w_m * (g_m / N) / ((1 + g_m / N)**(k + 1)))
    return pk_inner


def pk_k_mock_wrapper(N, N_mock, g_m, w_m, pk):
    @lru_cache(maxsize=2**15)
    def pk_k_mock_inner(k, k_mock):
        N_bar = (N + N_mock) / 2
        res = 0
        res -= np.log(pk(k_mock))
        res += np.log(binom(k + k_mock, k))
        res += k_mock * np.log(N_mock / N_bar)
        res += k * np.log(N / N_bar)
        res += np.log(np.sum(w_m * (g_m / N_bar) / ((2 + g_m / N_bar)**(k + k_mock + 1))))
        return np.exp(res)
    return pk_k_mock_inner


def pval_k_wrapper(pk_k_mock):
    @lru_cache(maxsize=2**15)
    def pval_k_inner(k, k_mock):
        p_sum = 1
        for k_i in range(k):
            p_sum -= pk_k_mock(k_i, k_mock)
            p_sum = max(0, p_sum)
        return p_sum
    return pval_k_inner


def pval_nk_wrapper(pn_k):
    @lru_cache(maxsize=2**15)
    def pval_n_inner(n, k):
        p_sum = 0
        for n_i in range(k, n + 1):
            p = pn_k(n_i, k)
            assert 0 <= p <= 1, '(n=%s, k=%s): invalid probability: %s' % (n, k, p)
            p_sum += p
        assert 0 <= p_sum <= 1.0001, "invalid p-value: %s" % p_sum
        return min(p_sum, 1)
    return pval_n_inner


def norm_wrapper(p_kn, p_n):
    @lru_cache(maxsize=2**15)
    def p_k_inner(k):
        p_k = 0
        for n in range(k, 50000):
            p_k += p_kn(k, n) * p_n(n)
        return p_k
    return p_k_inner


def p_nk_wrapper(p_kn, p_n, norm):
    @lru_cache(maxsize=2**15)
    def p_nk_inner(n, k):
        return p_kn(k, n) * p_n(n) / norm(k)
    return p_nk_inner


def p_n_wrapper(w, a, N):
    a = a / N

    @lru_cache(maxsize=2**15)
    def geom_distr(x):
        return np.sum(w * (a / (1 + a) ** x))
    return geom_distr


def plot_fit(k_vals, w, a, title, file_name, out_dir):

    def geom_distr(k, w, a):
        return np.sum(w * (a / (1 + a) ** (k + 1)))

    fit = [geom_distr(k, w, a) for k in range(np.max(k_vals) + 1)]
    fit = np.array(fit)
    fig, ax = plt.subplots()
    ax.hist(k_vals, log=True, bins=50, normed=True, alpha=0.5)
    ax.plot(fit, color='red')
    fig.suptitle(title, fontsize=18)
    fname = os.path.join(out_dir, file_name)
    plt.savefig(fname)


def main():
    parser = create_parser()
    args = parser.parse_args()

    with open(args.mock_model, 'rb') as model_file:
        mock_model = pickle.load(model_file)

    with open(args.bam_statistics_json) as json_handle:
        bam_statistics = json.load(json_handle)

    if args.estimate_null_fraction:
        if not (args.posterior_table and args.bam_statistics_json):
            print('estimating the null fraction requires commandline arguments for both'
                  '--posterior_table and --bam_statistics_json', file=sys.stderr)
            sys.exit(1)

        avg_read_len = bam_statistics['total_coverage'] / bam_statistics['total_reads']

        true_fraction = 0
        total_sites = 0
        true_coverage = 0
        with open(args.posterior_table) as post_handle:
            post_handle.readline()
            for line in post_handle:
                line = line.split()
                posterior = float(line[-1])
                k = int(line[2])
                total_sites += 1
                if k > 0:
                    true_fraction += posterior
                    true_coverage += k * posterior * avg_read_len

        n0 = bam_statistics['total_coverage'] - true_coverage
    else:
        total_sites = 0
        null_sites = 0
        with open(args.factor_mock_table) as post_handle:
            post_handle.readline()
            for line in post_handle:
                line = line.split()
                k = int(line[2])
                if k == 0:
                    null_sites += 1
                total_sites += 1

        n0 = bam_statistics['total_coverage']

    if args.plot_dir:
        if not os.path.exists(args.plot_dir):
            os.makedirs(args.plot_dir)

    def plot_pval_bf(p_values, bayes_factor, plot_file, title):
        if not args.plot_dir:
            return
        fig, ax = plt.subplots()
        ax.plot(p_values, bayes_factor)
        out_file = os.path.join(args.plot_dir, plot_file)
        ax.set_title(title)
        fig.savefig(out_file)
        plt.close()

    def plot_pval_ecdf(x_ecdf, y_ecdf, x_lcm, y_lcm, plot_file, title):
        if not args.plot_dir:
            return
        out_file = os.path.join(args.plot_dir, plot_file)
        fig, ax = plt.subplots()
        ax.plot(x_ecdf, y_ecdf, 'o')
        ax.plot(x_lcm, y_lcm)
        ax.set_title(title)
        fig.savefig(out_file)
        plt.close()

    max_k_mock = args.max_k_mock
    mock_cov = mock_model['coverage']
    pw_m, pg_m = mock_model['pk_params']

    pn_w, pn_g = mock_model['pn_params']
    alpha, beta = mock_model['pkn_params']

    if args.null_fraction is not None:
        n0 = int(bam_statistics['total_coverage'] * args.null_fraction)
    else:
        # this is a hack that sets the null fraction to 1.
        # by testing on various factors we have seen that setting this value does not seem to
        # change the predictions a lot.
        # everything between 0.8-2 seems to yield reasonable results.
        n0 = bam_statistics['total_coverage']
    factor_cov = n0

    if not 0.2 <= factor_cov / mock_cov <= 5:
        print('WARNING: trying to scale the mock coverage more than a factor of 5. ',
              'Reported posterior probabilities may be inaccurate.')

    if args.debug:
        print('factor_cov:', factor_cov)
        print('mock_cov:', mock_cov)
        print('factor_cov pct: %.4f' % (factor_cov / mock_cov))

    pk = pk_wrapper(mock_cov, pg_m, pw_m)
    pk_k_mock = pk_k_mock_wrapper(factor_cov, mock_cov, pg_m, pw_m, pk)
    pval_k_k_mock = pval_k_wrapper(pk_k_mock)

    p_n = p_n_wrapper(pn_w, pn_g, factor_cov)
    p_kn = p_kn_wrapper(alpha, beta)
    norm = norm_wrapper(p_kn, p_n)
    pn_k = p_nk_wrapper(p_kn, p_n, norm)
    pval_nk = pval_nk_wrapper(pn_k)

    # calculate p-values
    table_statistics = pd.read_table(args.factor_mock_statistics)
    # sum out n_mock
    table_statistics = table_statistics.groupby(['k_factor', 'n_factor', 'k_mock']).agg(
        {'count': sum}
    ).reset_index()

    k_lump = 10
    table_statistics.loc[:, 'k_lump'] = table_statistics.k_factor
    table_statistics.loc[table_statistics.k_lump >= k_lump, 'k_lump'] = k_lump

    pvals_nk = []
    pvals_kk = []
    for idx, row in table_statistics.iterrows():
        pvals_kk.append(pval_k_k_mock(row.k_factor, row.k_mock))
        # we do not have to make predictions for sites without transitions
        if row.k_factor > 0:
            pvals_nk.append(pval_nk(row.n_factor, row.k_factor))
        else:
            pvals_nk.append(np.nan)

    table_statistics['pval_kk'] = pvals_kk
    table_statistics['pval_nk'] = pvals_nk

    # set max(pval(k|k')) to 1
    for k_mock, kmock_df in table_statistics.groupby('k_mock'):
        if k_mock > max_k_mock:
            continue
        pvals = kmock_df.pval_kk
        table_statistics.loc[kmock_df.index, 'pval_kk'] = pvals / pvals.max()

    # set max(pval(n|k)) to 1
    for k, k_df in table_statistics.groupby('k_lump'):
        if k < 1:
            # we do not have to make predictions for sites without transitions
            continue
        pvals = k_df.pval_nk
        table_statistics.loc[k_df.index, 'pval_nk'] = pvals / pvals.max()

    # fit discrete grenander density estimations
    eta0_pval_thresh = 0.9

    for k_mock, kmock_df in table_statistics.groupby('k_mock'):
        if k_mock > max_k_mock:
            continue

        agg_df = kmock_df.groupby('pval_kk').aggregate({'count': sum}).reset_index()
        pvals = agg_df.pval_kk
        freqs = agg_df['count']

        # get a more robust eta0 estimation
        grenander_slopes, ecdf_x, ecdf_y, x_knots, y_knots = discrete_pval_grenander_fit(
            pvals, freqs, full_output=True)
        eta0 = estimate_eta0(pvals.as_matrix(), grenander_slopes, eta0_pval_thresh)

        grenander_slopes, ecdf_x, ecdf_y, x_knots, y_knots = discrete_pval_grenander_fit(
            pvals, freqs, eta0=eta0, full_output=True)

        plot_pval_ecdf(ecdf_x, ecdf_y, x_knots, y_knots, 'ecdf_kk%s' % k_mock,
                       'k_mock=%s' % k_mock)

        qk_prime = 1 - min(grenander_slopes)
        f_knots_bf_kk = (grenander_slopes / (1 - qk_prime)) - 1

        plot_pval_bf(pvals, f_knots_bf_kk, 'pval_kk%s_bf.png' % k_mock,
                     'k_mock=%s' % k_mock)

        agg_df.loc[:, 'bf_kk'] = f_knots_bf_kk
        merged_df = (
            kmock_df.merge(agg_df, on='pval_kk', how='left')
            .set_index(kmock_df.index)
            .dropna()
        )
        table_statistics.loc[merged_df.index, 'bf_kk'] = merged_df.bf_kk

    qk_store = {}
    for k, k_df in table_statistics.groupby('k_lump'):
        if k < 1:
            continue

        agg_df = k_df.groupby('pval_nk').aggregate({'count': sum}).reset_index()
        pvals = agg_df.pval_nk
        freqs = agg_df['count']

        # get a more robust eta0 estimation
        grenander_slopes, ecdf_x, ecdf_y, x_knots, y_knots = discrete_pval_grenander_fit(
            pvals, freqs, full_output=True)
        eta0 = estimate_eta0(pvals.as_matrix(), grenander_slopes, eta0_pval_thresh)

        # redo the grenander fit, but this time contrained by the estimated eta0
        grenander_slopes, ecdf_x, ecdf_y, x_knots, y_knots = discrete_pval_grenander_fit(
            pvals, freqs, eta0=eta0, full_output=True)

        plot_pval_ecdf(ecdf_x, ecdf_y, x_knots, y_knots, 'ecdf_nk%s.png' % k, 'k=%s' % k)

        qk = 1 - np.min(grenander_slopes)
        qk_store[k] = qk
        f_knots_bf_nk = (grenander_slopes / (1 - qk)) - 1

        plot_pval_bf(pvals, f_knots_bf_nk, 'pval_nk%s_bf.png' % k, 'k=%s' % k)

        agg_df.loc[:, 'bf_nk'] = f_knots_bf_nk
        merged_df = (
            k_df.merge(agg_df, on='pval_nk', how='left')
            .set_index(k_df.index)
            .dropna()
        )
        table_statistics.loc[merged_df.index, 'bf_nk'] = merged_df.bf_nk

    table_statistics = table_statistics[
        (table_statistics.k_factor > 0) & (table_statistics.k_mock <= max_k_mock)
    ]

    predictions = {}
    for idx, row in table_statistics.iterrows():
        # all predictions are positive
        if qk_store[row.k_lump] < 1e-20:
            posterior = 1
        else:
            qk = qk_store[row.k_lump]
            BF = row.bf_kk * row.bf_nk * (1 - qk) / qk

            if math.isinf(BF):
                posterior = 1
            elif math.isnan(BF):
                posterior = 0    
            else:
                posterior = BF / (1 + BF)

        predictions[row.k_factor, row.n_factor, row.k_mock] = (
            row.bf_kk, row.bf_nk, row.pval_kk, row.pval_nk, posterior
        )

    with open(args.factor_mock_table) as table,\
            open(args.outfile, 'w') as outfile:
        header = next(table).split()
        header += ['bf_kk', 'bf_nk', 'pval_kk', 'pval_nk', 'posterior']
        print(*header, sep='\t', file=outfile)

        for line in table:
            chrom, pos, k_factor, n_factor, k_mock, strand = line.split()
            try:
                bf_kk, bf_nk, pval_kk, pval_nk, posterior = predictions[
                    int(k_factor), int(n_factor), int(k_mock)
                ]
                print(chrom, pos, k_factor, n_factor, k_mock, strand, bf_kk, bf_nk,
                      pval_kk, pval_nk, posterior, sep='\t', file=outfile)
            except KeyError:
                continue


def estimate_eta0(pvals, slopes, threshold):
    ins = np.searchsorted(pvals, threshold, side='left')
    pvals_rel = np.append(threshold, pvals[ins:])
    slopes_rel = slopes[ins:]
    integral = np.dot(np.diff(pvals_rel), slopes_rel)
    return integral / (1 - threshold)


if __name__ == '__main__':
    main()
