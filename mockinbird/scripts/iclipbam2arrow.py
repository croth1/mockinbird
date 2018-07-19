from argparse import ArgumentParser
from subprocess import Popen, PIPE
import sys
import os
from os import path
import pyarrow as pa
from multiprocessing import Pool

from collections import defaultdict, namedtuple

from mockinbird.utils import cy_helpers
from dataclasses import dataclass, field
from typing import List
import pandas as pd


@dataclass
class SiteBuffer:
    pos: List[int] = field(default_factory=list)
    n: List[int] = field(default_factory=list)
    k: List[int] = field(default_factory=list)
    n_mock: List[int] = field(default_factory=list)
    k_mock: List[int] = field(default_factory=list)


class Region(namedtuple('Region', ['seqid', 'start', 'end'])):
    __slots__ = ()

    def __str__(self):
        return '%s:%s-%s' % (self.seqid, self.start, self.end)


def create_parser():
    parser = ArgumentParser()
    parser.add_argument('protein_bam')
    parser.add_argument('mock_bam')
    parser.add_argument('genome_fasta')
    parser.add_argument('out_dir')
    parser.add_argument('--n_processes', type=int)

    return parser


def create_df(data):
    df = pd.DataFrame({
        'pos': data.pos,
        'n': data.n,
        'k': data.k,
        'n_mock': data.n_mock,
        'k_mock': data.k_mock,
    })
    return df


def write_df(df, out_file):
    batch = pa.RecordBatch.from_pandas(df)
    writer = pa.RecordBatchFileWriter(out_file, batch.schema)
    writer.write_batch(batch)
    writer.close()


def main():
    parser = create_parser()
    args = parser.parse_args()

    fasta_idx = args.genome_fasta + '.fai'
    if not path.exists(fasta_idx):
        print('index: %s not found' % fasta_idx)
        sys.exit(1)

    regions = generate_regions(fasta_idx)

    plus_dir = path.join(args.out_dir, 'plus')
    minus_dir = path.join(args.out_dir, 'minus')

    for strand_dir in plus_dir, minus_dir:
        os.makedirs(strand_dir, exist_ok=True)

    with Pool(args.n_processes) as pool:
        jobs = []
        for region in regions:
            pileup_args = (region, args)

            job = pool.apply_async(pileup_region, args=pileup_args)
            jobs.append((region.seqid, job))

        for seqid, job in jobs:
            fwd_sites, rev_sites = job.get()
            fwd_df = create_df(fwd_sites)
            rev_df = create_df(rev_sites)
            del fwd_sites
            del rev_sites

            out_plus = path.join(plus_dir, seqid)
            write_df(fwd_df, out_plus)
            del fwd_df

            out_minus = path.join(minus_dir, seqid)
            write_df(rev_df, out_minus)
            del rev_df


def generate_regions(fasta_idx_file):
    regions = []
    with open(fasta_idx_file) as handle:
        for line in handle:
            seq_id, seq_len_str, *_ = line.split()
            seq_len = int(seq_len_str)
            region = Region(seq_id, 1, seq_len)
            regions.append(region)
    return regions


def pileup_region(region, args):
    pileup_cmd = 'samtools mpileup -aa -C 0 -d 1000000000 -q 0 -Q 0 -f %r %r -r %s'

    extra_args = {
        'universal_newlines': True,
        'stdout': PIPE,
        'stderr': PIPE,
        'shell': True
    }

    with Popen([pileup_cmd % (args.genome_fasta, args.protein_bam, region)], **extra_args) as factor_proc,\
            Popen([pileup_cmd % (args.genome_fasta, args.mock_bam, region)], **extra_args) as mock_proc:

        fwd_buffer = SiteBuffer()
        rev_buffer = SiteBuffer()

        factor = factor_proc.stdout
        mock = mock_proc.stdout

        sites = []
        site_statistics = defaultdict(int)

        try:
            factor_line = next(factor)
            mock_line = next(mock)
        except StopIteration:
            return sites, site_statistics

        while factor_line:
            factor_sites = cy_helpers.extract_ieclip_cy(factor_line)
            mock_sites = cy_helpers.extract_ieclip_cy(mock_line)
            for factor_site, mock_site in zip(factor_sites, mock_sites):
                f_seqid, f_pos, f_strand, f_n, f_k = factor_site
                m_seqid, m_pos, m_strand, m_n, m_k = mock_site

                buf = fwd_buffer if f_strand == '+' else rev_buffer
                buf.pos.append(f_pos)
                buf.n.append(f_n)
                buf.k.append(f_k)
                buf.n_mock.append(m_n)
                buf.k_mock.append(m_k)

            try:
                factor_line = next(factor)
                mock_line = next(mock)
            except StopIteration:
                break

    return fwd_buffer, rev_buffer


if __name__ == '__main__':
    main()
