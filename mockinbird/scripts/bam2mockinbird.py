from argparse import ArgumentParser
import subprocess
from subprocess import Popen, PIPE
import sys
from os import path
from multiprocessing import Pool
from tempfile import TemporaryDirectory

from collections import Counter, defaultdict, namedtuple
from types import SimpleNamespace

from dataclasses import dataclass
import pandas as pd


@dataclass
class BufferAgg:
    n: int
    k: int
    n_mock: int
    k_mock: int


SITE_COLS = ['chrom', 'pos', 'k_factor', 'n_factor', 'k_mock', 'strand']


class Region(namedtuple('Region', ['seqid', 'start', 'end'])):
    __slots__ = ()

    def __str__(self):
        return '%s:%s-%s' % (self.seqid, self.start, self.end)


class FullSite(namedtuple('Site', ['seqid', 'pos', 'strand', 'n', 'k', 'n_mock', 'k_mock'])):
    __slots__ = ()

    @property
    def _int_start(self):
        return self.pos

    @property
    def _int_end(self):
        return self.pos


class PileupSite(namedtuple('Site', ['seqid', 'pos', 'strand', 'n', 'k'])):
    __slots__ = ()


def create_parser():
    parser = ArgumentParser()
    parser.add_argument('factor_bam')
    parser.add_argument('mock_bam')
    parser.add_argument('genome_fasta')
    parser.add_argument('binding_sites')
    parser.add_argument('site_statistics')
    parser.add_argument('--bed_file')
    parser.add_argument('--tmp_dir')
    parser.add_argument('--transition_from', choices=['A', 'C', 'G', 'T'], default='T')
    parser.add_argument('--transition_to', choices=['A', 'C', 'G', 'T'], default='C')
    parser.add_argument('--n_processes', type=int)
    parser.add_argument('--aggregate_bp', type=int, default=1)
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

    bed_regions = defaultdict(list)
    with TemporaryDirectory(dir=args.tmp_dir) as tmp_dir:
        seq_bed_template = path.join(tmp_dir, '%s_regions.bed')
        if args.bed_file:
            with open(args.bed_file) as bed:
                for line in bed:
                    toks = line.split()
                    bed_regions[toks[0]].append(toks)

            for chrom, chrom_regions in bed_regions.items():
                with open(seq_bed_template % chrom, 'w') as out:
                    for reg in chrom_regions:
                        print(*reg, sep='\t', file=out)
        else:
            for region in regions:
                with open(seq_bed_template % region.seqid, 'w') as out:
                    for strand in '+', '-':
                        print(region.seqid, region.start - 1, region.end, '.', '.', strand,
                              sep='\t', file=out)

        with Pool(args.n_processes) as pool:
            jobs = []
            for region in regions:
                process_args = {
                    'factor_bam': args.factor_bam,
                    'mock_bam': args.mock_bam,
                    'genome_fasta': args.genome_fasta,
                    'region': region,
                    'window_size': args.aggregate_bp,
                    'tmp_dir': args.tmp_dir,
                    'bed_file': seq_bed_template % region.seqid,
                    'ref_nucleotide': args.transition_from,
                    'mut_nucleotide': args.transition_to,
                }

                job = pool.apply_async(pileup_region, args=(SimpleNamespace(**process_args),))
                jobs.append(job)

            for job in jobs:
                sites, stats = job.get()
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


def write_sites_header(sites_file):
    with open(sites_file, 'w') as file:
        print('chrom', 'pos', 'k_factor', 'n_factor', 'k_mock', 'strand',
              sep='\t', file=file)


def pileup_region(args):
    region = args.region
    with TemporaryDirectory(dir=args.tmp_dir) as tmp_dir:
        sites_file = path.join(tmp_dir, region.seqid + '.sites')
        stats_file = path.join(tmp_dir, region.seqid + '.stats')

        cmd = [
            'pcbam2mockinbird',
            args.factor_bam,
            args.mock_bam,
            args.genome_fasta,
            args.bed_file,
            region,
            args.window_size,
            sites_file,
            stats_file,
            args.ref_nucleotide,
            args.mut_nucleotide,
        ]

        subprocess.run([str(token) for token in cmd], check=True)

        stats = Counter()
        with open(stats_file) as stats_handle:
            for line in stats_handle:
                n_factor, k_factor, n_mock, k_mock, count = line.split()
                stats[int(k_factor), int(n_factor), int(k_mock), int(n_mock)] += int(count)

        C_SITE_FIELDS = ['pos', 'strand', 'n_factor', 'k_factor', 'k_mock']
        site_df = pd.read_table(sites_file, header=None, names=C_SITE_FIELDS, index_col=False)
        site_df['chrom'] = region.seqid

    return site_df, stats


def write_site_header(site_file):
    with open(site_file, 'w') as file:
        print(*SITE_COLS, sep='\t', file=file)


def write_sites(sites, sites_file):
    sites.loc[:, SITE_COLS].to_csv(sites_file, mode='a', sep='\t', header=False, index=False)


def write_statistics(stats, statistics_file):
    with open(statistics_file, 'w') as statistics_file:
        print('k_factor', 'n_factor', 'k_mock', 'n_mock', 'count', sep='\t', file=statistics_file)
        for (k, n, k_mock, n_mock), count in stats.items():
            print(k, n, k_mock, n_mock, count, sep='\t', file=statistics_file)


if __name__ == '__main__':
    main()
