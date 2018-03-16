import argparse
import sys
import subprocess

import mockinbird.utils.argparse_helper as aph


def register_arguments(parser):
    parser.add_argument('bam_file', type=aph.file_r, help='input bam file')
    parser.add_argument('output_bam', type=aph.file_rw, help='output bam file')
    parser.add_argument('--seed', type=int, default=42, help='random seed')
    parser.add_argument('n_reads', type=int, help='target read number')


def create_parser():
    description = 'subsample a bam file roughly to a given number of reads'
    parser = argparse.ArgumentParser(prog='mockinbird-subsample-bam', description=description)
    register_arguments(parser)
    return parser


def run(args):
    cmd = 'samtools view %r | wc -l' % args.bam_file
    proc = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
    n_entries = int(proc.stdout)

    if args.n_reads > n_entries:
        msg = (
            'Number of requested reads exceeds the number of reads in the bam file.\n'
            'This script can only subsample. Exiting.'
        )
        print(msg, file=sys.stderr)
        sys.exit(0)

    fraction = args.n_reads / n_entries
    subsample_fmt = '%s.%s' % (args.seed, str(fraction)[2:])

    with open(args.output_bam, 'w') as out:
        subsample_cmd = 'samtools view -b -h -s %s %r' % (subsample_fmt, args.bam_file)
        subprocess.run(subsample_cmd, shell=True, stdout=out)


def run_cmdline():
    parser = create_parser()
    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    run_cmdline()
