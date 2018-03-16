import argparse
import mockinbird.utils.argparse_helper as aph
import os
import glob
import subprocess
import tempfile
from os import path


def register_arguments(parser):
    parser.add_argument('output_bam', type=aph.file_rw, help='generated mock bam file')
    parser.add_argument('bam_folder', type=aph.dir_rx, help='folder with sorted bam files')
    fraction_help = 'number of datasets covering a site for it to be called a background site'
    parser.add_argument('--ubiquity_fraction', type=aph.float_low0_high1, default=0.33,
                        help=fraction_help)
    # parser.add_argument('--seed', type=int, default=42, help='seed for random generators')
    tmp_help = 'parent folder for storing temporary files'
    parser.add_argument('--tmp_root', type=aph.dir_rwx, help=tmp_help)


def create_parser():
    description = 'create a mock experiment from a set of CLIP datasets'
    parser = argparse.ArgumentParser(prog='mockinbird-create-metamock', description=description)
    register_arguments(parser)
    return parser


def mergefile2bed(merge_file, bed_file):
    with open(merge_file) as merge_file, open(bed_file, 'w') as bed:
        for line in merge_file:
            chrom, start, end, strand = line.split()
            print(chrom, start, end, '.', '.', strand, sep='\t', file=bed)


def extract_mock_positions(coverage_file, outfile, min_coverage):
    with open(coverage_file) as coverage_file, open(outfile, 'w') as out:
        for line in coverage_file:
            toks = line.split()
            chrom = toks[0]
            start = int(toks[1]) + int(toks[6]) - 1
            coverage = int(toks[7])
            strand = toks[5]

            if coverage >= min_coverage:
                print(chrom, start, start + 1, '.', coverage, strand, sep='\t', file=out)


def run_command(tokens, out=None):
    print(' '.join(tokens))
    subprocess.run(' '.join(tokens), shell=True, stdout=out)


def run(args):
    with tempfile.TemporaryDirectory(dir=args.tmp_root) as tmpdir:
        bed_dir = path.join(tmpdir, 'bed')
        bam_dir = path.join(tmpdir, 'bam')
        os.mkdir(bed_dir)
        os.mkdir(bam_dir)

        bam_files = glob.glob(path.join(args.bam_folder, '*.bam'))
        n_experiments = len(bam_files)

        # create covered regions from each bam file
        for bam_file in bam_files:
            prefix, _ = path.splitext(path.basename(bam_file))
            merge_file = path.join(tmpdir, 'bed', prefix + '.iv')
            cmd = [
                'bedtools', 'bamtobed', '-i %r' % bam_file,
                '|', 'bedtools', 'sort',
                '|', 'bedtools', 'merge', '-s', '-c 6', '-o', 'distinct',
            ]

            with open(merge_file, 'w') as out:
                run_command(cmd, out)

            bed_file = path.join(tmpdir, 'bed', prefix + '.merged_bed')
            mergefile2bed(merge_file, bed_file)

        # create one file with all regions
        bed_files = glob.glob(path.join(bed_dir, '*.merged_bed'))

        cmd = ['cat'] + bed_files + [
            '|', 'bedtools', 'sort',
            '|', 'bedtools', 'merge', '-s', '-c 6', '-o', 'distinct',
        ]

        combined_iv = path.join(tmpdir, 'combined.iv')
        with open(combined_iv, 'w') as out:
            run_command(cmd, out)

        combined_bed = path.join(tmpdir, 'combined.bed')
        mergefile2bed(combined_iv, combined_bed)

        # calculate coverage for joint regions
        cmd = [
            'bedtools', 'coverage', '-a', '%r' % combined_bed, '-s', '-d', '-b'
        ] + bed_files

        coverage_bed = path.join(tmpdir, 'combined_coverage.bed')
        with open(coverage_bed, 'w') as out:
            run_command(cmd, out)

        # generate mock positions
        mock_cov_bed = path.join(tmpdir, 'mock_cov.bed')
        min_coverage = int(args.ubiquity_fraction * n_experiments)
        extract_mock_positions(coverage_bed, mock_cov_bed, min_coverage)

        # merge to a full mock bed file
        mock_iv = path.join(tmpdir, 'mock.iv')
        mock_bed = path.join(tmpdir, 'mock.bed')

        cmd = [
            'bedtools', 'sort', '-i %r' % mock_cov_bed,
            '|', 'bedtools', 'merge', '-s', '-c 6', '-o', 'distinct',
        ]
        with open(mock_iv, 'w') as out:
            run_command(cmd, out)
        mergefile2bed(mock_iv, mock_bed)

        # create the mock bam files
        for bam_file in bam_files:
            prefix, _ = path.splitext(path.basename(bam_file))
            out_file = path.join(bam_dir, prefix + '.mock.bam')
            cmd = [
                'samtools', 'view', '-b', '-h',
                '%r' % bam_file,
                '-L %r' % mock_bed,
                '-o %r' % out_file
            ]
            run_command(cmd)

        # create merged_mock
        mock_bams = glob.glob(path.join(bam_dir, '*.mock.bam'))
        cmd = ['samtools', 'merge', '%r' % args.output_bam] + mock_bams
        run_command(cmd)


def run_cmdline():
    parser = create_parser()
    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    run_cmdline()
