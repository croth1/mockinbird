from argparse import ArgumentParser
from subprocess import Popen, PIPE
import sys
from os import path
from multiprocessing import Pool

from collections import Counter, defaultdict, namedtuple

from mockinbird.utils import cy_helpers
from dataclasses import dataclass

@dataclass
class BufferAgg:
    n: int
    k: int
    n_mock: int
    k_mock: int


class Region(namedtuple('Region', ['seqid', 'start', 'end'])):
    __slots__ = ()

    def __str__(self):
        return '%s:%s-%s' % (self.seqid, self.start, self.end)


class FullSite(namedtuple('Site', ['seqid', 'pos', 'strand', 'n', 'k', 'n_mock', 'k_mock'])):
    __slots__ = ()


class PileupSite(namedtuple('Site', ['seqid', 'pos', 'strand', 'n', 'k'])):
    __slots__ = ()


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
    parser.add_argument('--ieclip_model', action='store_true')
    parser.add_argument('--aggregate_bp', type=int, default=1)
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
    write_sites_header(args.binding_sites)

    statistics = Counter()

    if args.n_processes == 1:
        for region in regions:
            sites, stats = pileup_region(region, args)
            statistics.update(stats)
            write_sites(sites, args.binding_sites)

    else:
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
            seq_id, seq_len_str, *_ = line.split()
            seq_len = int(seq_len_str)
            region = Region(seq_id, 1, seq_len)
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
        return (chrom, pos, strand, n, k)

    with Popen([pileup_cmd % (args.genome_fasta, args.protein_bam, region)], **extra_args) as factor_proc,\
            Popen([pileup_cmd % (args.genome_fasta, args.mock_bam, region)], **extra_args) as mock_proc:

        agg_bp = args.aggregate_bp
        min_k = args.min_k
        is_ieclip = args.ieclip_model

        if agg_bp > 1:
            fwd_buffer = []
            fwd_agg = BufferAgg(0, 0, 0, 0)
            rev_buffer = []
            rev_agg = BufferAgg(0, 0, 0, 0)

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
            extracted_sites = []
            # PARCLIP mode
            if not is_ieclip:
                factor_site = extract_transition(factor_line)
                # None is a site that cannot be mutated
                if factor_site:
                    mock_site = extract_transition(mock_line)
                    f_seqid, f_pos, f_strand, f_n, f_k = factor_site
                    m_seqid, m_pos, m_strand, m_n, m_k = mock_site
                    full_site = FullSite(
                        f_seqid,
                        f_pos,
                        f_strand,
                        f_n,
                        f_k,
                        m_n,
                        m_k,
                    )
                    extracted_sites.append(full_site)
            else:
                factor_sites = cy_helpers.extract_ieclip_cy(factor_line)
                mock_sites = cy_helpers.extract_ieclip_cy(mock_line)
                for factor_site, mock_site in zip(factor_sites, mock_sites):
                    f_seqid, f_pos, f_strand, f_n, f_k = factor_site
                    m_seqid, m_pos, m_strand, m_n, m_k = mock_site

                    full_site = FullSite(
                        f_seqid,
                        f_pos,
                        f_strand,
                        f_n,
                        f_k,
                        m_n,
                        m_k,
                    )
                    extracted_sites.append(full_site)

            output_sites = []
            for site in extracted_sites:
                if agg_bp > 1:
                    # aggregation logic comes here
                    if site.strand == '+':
                        site_buffer = fwd_buffer

                    if site.strand:
                        site_buffer = fwd_buffer
                        agg_buffer = fwd_agg
                    else:
                        site_buffer = rev_buffer
                        agg_buffer = rev_agg

                    # we cannot make predictions for the first and last agg_bp // 2 sites
                    if len(site_buffer) == agg_bp - 1:

                        agg_buffer.n += site.n
                        agg_buffer.k += site.k
                        agg_buffer.n_mock += site.n_mock
                        agg_buffer.k_mock += site.k_mock

                        current_pos = site_buffer[len(site_buffer) // 2]
                        site_buffer.append(site)
                        prev_site = site_buffer.pop(0)

                        output_site = FullSite(
                            current_pos.seqid,
                            current_pos.pos,
                            current_pos.strand,
                            n=agg_buffer.n,
                            k=agg_buffer.k,
                            n_mock=agg_buffer.n_mock,
                            k_mock=agg_buffer.k_mock,
                        )
                        output_sites.append(output_site)

                        agg_buffer.n -= prev_site.n
                        agg_buffer.k -= prev_site.k
                        agg_buffer.n_mock -= prev_site.n_mock
                        agg_buffer.k_mock -= prev_site.k_mock

                    else:
                        # buffer not full yet, we cannot aggregate
                        site_buffer.append(site)
                        agg_buffer.n += site.n
                        agg_buffer.k += site.k
                        agg_buffer.n_mock += site.n_mock
                        agg_buffer.k_mock += site.k_mock

                else:
                    output_sites.append(site)

            for site in output_sites:
                site_statistics[site.k, site.n, site.k_mock, site.n_mock] += 1

                if site.k >= min_k:
                    sites.append((
                        site.seqid,
                        site.pos,
                        site.k,
                        site.n,
                        site.k_mock,
                        site.strand
                    ))

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
