from argparse import ArgumentParser
from subprocess import Popen, PIPE
import sys
from os import path
from multiprocessing import Pool

from collections import Counter, namedtuple
from dataclasses import dataclass


@dataclass
class Site:
    pos: int
    k: int


@dataclass
class BedRegion:
    seqid: str
    start: int
    end: int
    strand: str


@dataclass
class BufferAgg:
    n_interesting_sites: int


class Region(namedtuple('Region', ['seqid', 'start', 'end'])):
    __slots__ = ()

    def __str__(self):
        return '%s:%s-%s' % (self.seqid, self.start, self.end)


def create_parser():
    parser = ArgumentParser()
    parser.add_argument('protein_bam')
    parser.add_argument('genome_fasta')
    parser.add_argument('regions_of_interest')
    parser.add_argument('--transition_from', choices=['A', 'C', 'G', 'T'], default='T')
    parser.add_argument('--transition_to', choices=['A', 'C', 'G', 'T'], default='C')
    parser.add_argument('--n_processes', type=int)
    parser.add_argument('--window_size', type=int, default=7)
    parser.add_argument('--min_k', type=int, default=1)

    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    fasta_idx = args.genome_fasta + '.fai'
    if not path.exists(fasta_idx):
        print('index: %s not found' % fasta_idx)
        sys.exit(1)

    regions = generate_regions(fasta_idx)

    if args.n_processes == 1:
        for region in regions:
            regions_of_interest = extract_regions_of_interest(region, args)
            write_regions(regions_of_interest, args.regions_of_interest)

    else:
        with Pool(args.n_processes) as pool:
            jobs = []
            for region in regions:
                pileup_args = (region, args)
                job = pool.apply_async(extract_regions_of_interest, args=pileup_args)
                jobs.append(job)

            for job in jobs:
                regions_of_interest = job.get()
                write_regions(regions_of_interest, args.regions_of_interest)


def generate_regions(fasta_idx_file):
    regions = []
    with open(fasta_idx_file) as handle:
        for line in handle:
            seq_id, seq_len_str, *_ = line.split()
            seq_len = int(seq_len_str)
            region = Region(seq_id, 1, seq_len)
            regions.append(region)
    return regions


def extract_regions_of_interest(region, args):
    pileup_cmd = 'samtools mpileup -aa -C 0 -d 1000000000 -q 0 -Q 0 -f %r %r -r %s'

    extra_args = {
        'universal_newlines': True,
        'stdout': PIPE,
        'stderr': PIPE,
        'shell': True
    }

    rev_map = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }

    plus_transition_from = args.transition_from
    minus_transition_from = rev_map[args.transition_from]
    plus_transition_to = args.transition_to
    minus_transition_to = rev_map[args.transition_to].lower()

    def extract_transition(line):
        chrom, pos, ref_nuc, cov, cov_str, qual_str = line.split()
        if ref_nuc == plus_transition_from:
            strand = '+'
        elif ref_nuc == minus_transition_from:
            strand = '-'
        else:
            return None

        if cov == '0':
            n = k = 0
        else:
            cov_symb = Counter(cov_str)
            if ref_nuc == plus_transition_from:
                k = cov_symb[plus_transition_to]
                n = k + cov_symb['.']
            elif ref_nuc == minus_transition_from:
                k = cov_symb[minus_transition_to]
                n = k + cov_symb[',']
        return (chrom, int(pos), strand, n, k)

    with Popen([pileup_cmd % (args.genome_fasta, args.protein_bam, region)], **extra_args) as factor_proc:

        agg_bp = args.window_size
        min_k = args.min_k

        if agg_bp > 1:
            fwd_buffer = []
            fwd_agg = BufferAgg(0)
            rev_buffer = []
            rev_agg = BufferAgg(0)

        factor = factor_proc.stdout

        regions_of_interest = []

        try:
            factor_line = next(factor)
        except StopIteration:
            return regions_of_interest

        while factor_line:
            factor_site = extract_transition(factor_line)
            if factor_site:
                seqid, pos, strand, n, k = factor_site
                site = Site(pos=pos, k=k)

                if strand == '+':
                    site_buffer = fwd_buffer
                    agg_buffer = fwd_agg
                else:
                    site_buffer = rev_buffer
                    agg_buffer = rev_agg

                site_buffer.append(site)
                agg_buffer.n_interesting_sites += (k > 0)

                # we cannot make predictions for the first and last agg_bp // 2 sites
                if len(site_buffer) == agg_bp:
                    current_site = site_buffer[len(site_buffer) // 2]
                    prev_site = site_buffer.pop(0)
                    agg_buffer.n_interesting_sites -= (prev_site.k > 0)

                    if current_site.k >= min_k:
                        bed_region = BedRegion(seqid, prev_site.pos, current_site.pos + 1, strand)
                        regions_of_interest.append(bed_region)

            try:
                factor_line = next(factor)
            except StopIteration:
                break

    return regions_of_interest


def write_regions(bed_regions, sites_file):
    with open(sites_file, 'a') as handle:
        for r in bed_regions:
            print(r.seqid, r.start, r.end, '.', '.', r.strand, sep='\t', file=handle)


if __name__ == '__main__':
    main()
