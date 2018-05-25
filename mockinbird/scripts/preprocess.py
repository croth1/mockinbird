import argparse
import os
import sys
import logging
from collections import OrderedDict


from mockinbird.utils import argparse_helper as aph
from mockinbird import LOG_DEFAULT_FORMAT, LOG_LEVEL_MAP
from mockinbird.utils import prepare_output_dir
from mockinbird.utils import config_validation as cv
from mockinbird.utils import module_utils as mu
from mockinbird.utils import pipeline as pl
from mockinbird import __version__

logger = logging.getLogger()
cur_dir = os.path.dirname(os.path.realpath(__file__))
scriptPath = os.path.join(cur_dir, 'utils')


def register_arguments(parser):
    parser.add_argument('parclip_fastq', help='path to PAR-CLIP fastq reads',
                        type=aph.file_r)
    outdir_help = 'output directory - will be created if it does not exist'
    parser.add_argument('output_dir', help=outdir_help, type=aph.dir_rwx_create)
    parser.add_argument('prefix', help='filename prefix for newly created files')
    parser.add_argument('config_file', help='path to preprocessing config file',
                        type=aph.file_r)
    parser.add_argument('--log_level', help='verbosity level of the logger',
                        choices=LOG_LEVEL_MAP.keys(), default='info')
    aph.add_version_arguments(parser)


def create_parser():
    description = 'run the PAR-CLIP preprocessing pipeline'
    parser = argparse.ArgumentParser(prog='mockinbird-preprocess', description=description)
    register_arguments(parser)
    return parser


def prepare_dir_or_die(dir_path):
    try:
        prepare_output_dir(dir_path)
    except ValueError as e:
        logger.error('Error while creating output directory: %s', dir_path)
        logger.error(e)
        sys.exit(1)


def run(args):

    inputfile = args.parclip_fastq
    outputdir = args.output_dir
    prefix = args.prefix

    prepare_dir_or_die(outputdir)

    # activate logging
    logging_file = os.path.join(outputdir, 'preprocess.log')
    logger = logging.getLogger()
    logger.setLevel(LOG_LEVEL_MAP[args.log_level])
    formatter = logging.Formatter(LOG_DEFAULT_FORMAT)

    console_handler = logging.StreamHandler(stream=sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    file_handler = logging.FileHandler(logging_file, mode='w')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    logger.info('mockinbird version: %s', __version__)
    logger.info('working directory: %s', os.getcwd())
    logger.info('started preprocessing via %r', ' '.join(sys.argv))

    config = mu.parse_yaml(args.config_file)

    def relpath_conv(file_path):
        return cv.rel_file_r_validator(file_path, args.config_file)

    def genomefasta_validator(file_path):
        return cv.rel_genome_validator(file_path, args.config_file)

    general_fmt = OrderedDict([
        ('adapter5prime', cv.Annot(str, converter=cv.dnastr_validator)),
        ('adapter3prime', cv.Annot(str, converter=cv.dnastr_validator)),
        ('genomefasta', cv.Annot(str, converter=genomefasta_validator)),
        ('rmTemp', cv.Annot(cv.boolean, default=True)),
        ('n_threads', cv.Annot(int, default=2)),
    ])

    reads_fmt = OrderedDict([
        ('bc_5prime', cv.Annot(int, default=0, converter=cv.nonneg_integer)),
        ('bc_3prime', cv.Annot(int, default=0, converter=cv.nonneg_integer)),
        ('min_len', cv.Annot(int, default=20, converter=cv.nonneg_integer)),
        ('reference_nucleotide', cv.Annot(str, default='T', converter=cv.dnanuc_validator)),
        ('mutation_nucleotide', cv.Annot(str, default='C', converter=cv.dnanuc_validator)),
    ])

    mandatory_sections = 'pipeline', 'general', 'reads'
    for section in mandatory_sections:
        if section not in config:
            logger.error('the config file does not define the mandatory section %s', section)
            sys.exit(1)

    try:
        general_raw = config['general']
        general_cfg = cv.validate_section(general_raw, general_fmt)
    except cv.ConfigError:
        logger.error('Error while parsing section %r', 'general')
        sys.exit(1)

    general_cfg['prefix'] = prefix
    general_cfg['output_dir'] = outputdir

    try:
        reads_raw = config['reads']
        reads_cfg = cv.validate_section(reads_raw, reads_fmt)
    except cv.ConfigError:
        logger.error('Error while parsing section %r', 'reads')
        sys.exit(1)

    initial_files = {'fastq': inputfile}

    if 'custom_files' in config:
        for fmt, path in config['custom_files'].items():
            if isinstance(path, str):
                try:
                    rel_path = relpath_conv(path)
                    initial_files[fmt] = rel_path
                except ValueError:
                    logger.warning('key %r: invalid file path %r', fmt, path)

    gencfg = {
        'reads': reads_cfg,
        'general': general_cfg,
    }
    pipeline = pl.Pipeline(initial_files=initial_files, general_cfg=gencfg,
                           cfg_path=args.config_file)
    mu.queue_pipeline(config, pipeline, def_lookup_path='mockinbird.utils.preprocess_modules')
    mu.run_pipeline(pipeline)

    if general_cfg['rmTemp']:
        pipeline.cleanup()

    logger.info('all done. See you soon!')


def main():
    parser = create_parser()
    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()
