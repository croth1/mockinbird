from argparse import ArgumentParser
from subprocess import Popen, PIPE
import sys
from os import path
from multiprocessing import Pool

from collections import Counter, defaultdict, namedtuple


class Region(namedtuple('Region', ['seqid', 'start', 'end'])):
    __slots__ = ()
    def __str__(self):
        return '%s:%s-%s' % (self.seqid, self.start, self.end)


def create_parser():
    parser = ArgumentParser()
    parser.add_argument('protein_bam')
    parser.add_argument('mock_bam')
    parser.add_argument('genome_fasta')
    parser.add_argument('binding_sites')
    parser.add_argument('site_statistics')
    parser.add_argument('--transition_from', choices=['A', 'C', 'G', 'T'], default='T')
    parser.add_argument('--transition_to', choices=['A', 'C', 'G', 'T'], default='C')
    parser.add_argument('--n_processes', type=int)

    return parser


BP_PER_REGION = 10000000


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
            pileup_args = (region, args)

            job = pool.apply_async(pileup_region, args=pileup_args)
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
            seq_id, bp_str, *_ = line.split()
            bp = int(bp_str)
            for start in range(1, bp, BP_PER_REGION):
                region = Region(seq_id, start, start + BP_PER_REGION - 1)
                regions.append(region)
    return regions


def write_sites_header(sites_file):
    with open(sites_file, 'w') as file:
        print('chrom', 'pos', 'k_factor', 'n_factor', 'k_mock', 'strand',
              sep='\t', file=file)


def pileup_region(region, args):
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
            return chrom, pos, None, 0, 0

        if cov == '0':
            n, k = 0, 0
        else:
            cov_symb = Counter(cov_str)
            if ref_nuc == plus_transition_from:
                k = cov_symb[plus_transition_to]
                n = k + cov_symb['.']
            elif ref_nuc == minus_transition_from:
                k = cov_symb[minus_transition_to]
                n = k + cov_symb[',']
        return chrom, pos, strand, k, n

    with Popen([pileup_cmd % (args.genome_fasta, args.protein_bam, region)], **extra_args) as factor_proc,\
            Popen([pileup_cmd % (args.genome_fasta, args.mock_bam, region)], **extra_args) as mock_proc:

        factor = factor_proc.stdout
        mock = mock_proc.stdout

        factor_line = next(factor)
        mock_line = next(mock)

        sites = []
        site_statistics = defaultdict(int)

        while factor_line:
            chrom, pos, strand, k_factor, n_factor = extract_transition(factor_line)
            if strand:
                _, _, _, k_mock, n_mock = extract_transition(mock_line)
                site_statistics[k_factor, n_factor, k_mock, n_mock] += 1

                if k_factor > 0:
                    sites.append((chrom, pos, k_factor, n_factor, k_mock, strand))

            try:
                factor_line = next(factor)
                mock_line = next(mock)
            except StopIteration:
                break

    return sites, site_statistics


def write_sites(sites, sites_file):
    with open(sites_file, 'a') as handle:
        for site in sites:
            print(*site, sep='\t', file=handle)


def write_statistics(stats, statistics_file):
    with open(statistics_file, 'w') as statistics_file:
        print('k_factor', 'n_factor', 'k_mock', 'n_mock', 'count', sep='\t', file=statistics_file)
        for (k, n, k_mock, n_mock), count in stats.items():
            print(k, n, k_mock, n_mock, count, sep='\t', file=statistics_file)


if __name__ == '__main__':
    main()
