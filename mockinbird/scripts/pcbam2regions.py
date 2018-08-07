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
    parser.add_argument('genome_fasta')
    parser.add_argument('regions_bed')
    parser.add_argument('--tmp_dir')
    parser.add_argument('--n_processes', type=int)
    parser.add_argument('--window_size', type=int, default=7)
    parser.add_argument('--transition_from', choices=['A', 'C', 'G', 'T'], default='T')
    parser.add_argument('--transition_to', choices=['A', 'C', 'G', 'T'], default='C')
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    fasta_idx = args.genome_fasta + '.fai'
    if not path.exists(fasta_idx):
        print('index: %s not found' % fasta_idx)
        sys.exit(1)

    regions = generate_regions(fasta_idx)
    write_region_header(args.regions_bed)

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
            regions = job.get()
            write_regions(regions, args.regions_bed)


def generate_regions(fasta_idx_file):
    regions = []
    with open(fasta_idx_file) as handle:
        for line in handle:
            seq_id, seq_len_str, *_ = line.split()
            seq_len = int(seq_len_str)
            region = Region(seq_id, 1, seq_len)
            regions.append(region)
    return regions


BED_COLS = ['seqid', 'start', 'end', 'name', 'score', 'strand']


def write_region_header(region_file):
    with open(region_file, 'w') as file:
        # empty the file
        pass
        #print(*BED_COLS, sep='\t', file=file)


def process_region(region, args):
    with TemporaryDirectory(dir=args.tmp_dir) as tmp_dir:
        region_file = path.join(tmp_dir, region.seqid + '_regions.bed')

        cmd = [
            'pcbam2regions',
            args.factor_bam,
            args.genome_fasta,
            args.region,
            args.window_size,
            region_file,
            args.transition_from,
            args.transition_to,
        ]

        subprocess.run([str(token) for token in cmd])

        region_df = pd.read_table(region_file, names=['start', 'end', 'strand'])
        region_df.loc[:, 'chrom'] = region.seqid
        return region_df


def write_regions(sites, sites_file):
    sites['name'] = '.'
    sites['score'] = '.'
    sites[BED_COLS].to_csv(sites_file, mode='a', sep='\t', header=False, index=False)


if __name__ == '__main__':
    main()
