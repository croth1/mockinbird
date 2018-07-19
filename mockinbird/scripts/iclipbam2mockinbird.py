from argparse import ArgumentParser
import sys
from os import path
from multiprocessing import Pool
from tempfile import TemporaryDirectory
from collections import Counter, namedtuple
from types import SimpleNamespace
import subprocess

import pandas as pd


class Region(namedtuple('Region', ['seqid', 'start', 'end'])):
    __slots__ = ()

    def __str__(self):
        return '%s:%s-%s' % (self.seqid, self.start, self.end)


def create_parser():
    parser = ArgumentParser()
    parser.add_argument('factor_bam')
    parser.add_argument('mock_bam')
    parser.add_argument('genome_fasta')
    parser.add_argument('binding_sites')
    parser.add_argument('site_statistics')
    parser.add_argument('--tmp_dir')
    parser.add_argument('--n_processes', type=int)
    parser.add_argument('--window_size', type=int, default=5)
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    fasta_idx = args.genome_fasta + '.fai'
    if not path.exists(fasta_idx):
        print('index: %s not found' % fasta_idx)
        sys.exit(1)

    regions = generate_regions(fasta_idx)
    write_sites_header(args.binding_sites)

    statistics = Counter()

    with Pool(args.n_processes) as pool:
        jobs = []
        for region in regions:

            process_args = {
                'factor_bam': args.factor_bam,
                'mock_bam': args.mock_bam,
                'genome_fasta': args.genome_fasta,
                'region': region,
                'window_size': args.window_size,
                'tmp_dir': args.tmp_dir,
            }

            job = pool.apply_async(process_region, args=(region, SimpleNamespace(**process_args)))
            jobs.append(job)

        for job in jobs:
            stats, sites = job.get()
            statistics.update(stats)
            write_sites(sites, args.binding_sites)

    write_statistics(statistics, args.site_statistics)


def generate_regions(fasta_idx_file):
    regions = []
    with open(fasta_idx_file) as handle:
        for line in handle:
            seq_id, seq_len_str, *_ = line.split()
            seq_len = int(seq_len_str)
            region = Region(seq_id, 1, seq_len)
            regions.append(region)
    return regions


SITE_COLS = ['chrom', 'pos', 'k_factor', 'n_factor', 'k_mock', 'strand']


def write_sites_header(sites_file):
    with open(sites_file, 'w') as file:
        print(*SITE_COLS, sep='\t', file=file)


def process_region(region, args):
    with TemporaryDirectory(dir=args.tmp_dir) as tmp_dir:
        statistics_file = path.join(tmp_dir, region.seqid + '.stat')
        sites_file = path.join(tmp_dir, region.seqid + '.sites')

        cmd = [
            'iclipbam2mockinbird',
            args.factor_bam,
            args.mock_bam,
            args.genome_fasta,
            args.region,
            args.window_size,
            sites_file,
            statistics_file,
        ]

        subprocess.run([str(token) for token in cmd])

        stats = Counter()
        with open(statistics_file) as stat:
            for line in stat:
                n_str, k_str, n_mock_str, k_mock_str, count_str = line.split()
                stats[int(k_str), int(n_str), int(k_mock_str), int(n_mock_str)] = int(count_str)

        sites_df = pd.read_table(sites_file, names=['pos', 'strand', 'n_factor', 'k_factor', 'k_mock'])
        sites_df.loc[:, 'chrom'] = region.seqid
        return stats, sites_df


def write_sites(sites, sites_file):
    sites[SITE_COLS].to_csv(sites_file, mode='a', sep='\t', header=False, index=False)


def write_statistics(stats, statistics_file):
    with open(statistics_file, 'w') as statistics_file:
        print('k_factor', 'n_factor', 'k_mock', 'n_mock', 'count', sep='\t', file=statistics_file)
        for (k, n, k_mock, n_mock), count in stats.items():
            print(k, n, k_mock, n_mock, count, sep='\t', file=statistics_file)


if __name__ == '__main__':
    main()
