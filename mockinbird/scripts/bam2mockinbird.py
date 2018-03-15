from argparse import ArgumentParser
from subprocess import Popen, PIPE

from collections import Counter, defaultdict


def create_parser():
    parser = ArgumentParser()
    parser.add_argument('protein_bam')
    parser.add_argument('mock_bam')
    parser.add_argument('genome_fasta')
    parser.add_argument('binding_sites')
    parser.add_argument('site_statistics')
    parser.add_argument('--transition_from', choices=['A', 'C', 'G', 'T'], default='T')
    parser.add_argument('--transition_to', choices=['A', 'C', 'G', 'T'], default='C')

    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    pileup_cmd = 'samtools mpileup -aa -C 0 -d 1000000 -q 0 -Q 0 -f %r %r'

    extra_args = {
        'universal_newlines': True,
        'stdout': PIPE,
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

    with Popen([pileup_cmd % (args.genome_fasta, args.protein_bam)], **extra_args) as factor_proc,\
            Popen([pileup_cmd % (args.genome_fasta, args.mock_bam)], **extra_args) as mock_proc,\
            open(args.binding_sites, 'w') as sites_file:

        factor = factor_proc.stdout
        mock = mock_proc.stdout

        print('chrom', 'pos', 'k_factor', 'n_factor', 'k_mock', 'strand',
              sep='\t', file=sites_file)

        factor_line = next(factor)
        mock_line = next(mock)

        site_statistics = defaultdict(int)

        while factor_line:
            chrom, pos, strand, k_factor, n_factor = extract_transition(factor_line)
            if strand:
                _, _, _, k_mock, n_mock = extract_transition(mock_line)
                site_statistics[k_factor, n_factor, k_mock, n_mock] += 1

                if k_factor > 0:
                    print(chrom, pos, k_factor, n_factor, k_mock, strand,
                          sep='\t', file=sites_file)

            try:
                factor_line = next(factor)
                mock_line = next(mock)
            except StopIteration:
                break

        with open(args.site_statistics, 'w') as statistics_file:
            print('k_factor', 'n_factor', 'k_mock', 'n_mock', 'count', sep='\t', file=statistics_file)
            for (k, n, k_mock, n_mock), count in site_statistics.items():
                print(k, n, k_mock, n_mock, count, sep='\t', file=statistics_file)


if __name__ == '__main__':
    main()
