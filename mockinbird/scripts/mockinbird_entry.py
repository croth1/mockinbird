import argparse
import sys

try:
    from mockinbird import __version__
except ImportError:
    __version__ = 'unknown'

from mockinbird.scripts.preprocess import register_arguments as preprocess_register_args
from mockinbird.scripts.preprocess import run as preprocess_main
from mockinbird.scripts.postprocess import register_arguments as postprocess_register_args
from mockinbird.scripts.postprocess import run as postprocess_main
from mockinbird.scripts.flip_read2 import register_arguments as flipread2_register_args
from mockinbird.scripts.flip_read2 import run as flipread2_main
from mockinbird.utils.cmd_tools import check_config as check_cfg
from mockinbird.utils.cmd_tools import create_metamock as metamock_cfg
from mockinbird.utils.cmd_tools import subsample_bam as subsample_cfg


def create_parser():
    parser = argparse.ArgumentParser(
        'mockinbird',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    subparser = parser.add_subparsers(title='subcommands')
    for mod_cls in StammpModule.__subclasses__():
        mod_cls(subparser)
    parser.add_argument('--version', action='version',
                        version=__version__)
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    if not sys.argv[1:]:
        args = parser.parse_args(['--help'])
    else:
        args = parser.parse_args()

    args.subcommand_func(args)


class StammpModule(object):

    subcommand = None
    aliases = []

    def __init__(self, parser, description=None, help=None):
        subcommand_parser = parser.add_parser(
            self.subcommand,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            description=description,
            help=help,
            aliases=self.aliases
        )

        subcommand_parser.set_defaults(subcommand_func=self)
        self.subcommand_parser = subcommand_parser

    def __call__(self, args):
        pass


class PreprocessModule(StammpModule):
    subcommand = 'preprocess'

    def __init__(self, parser):
        help_msg = 'run preprocessing pipeline'
        description = 'start preprocessing pipeline using a config script'
        super().__init__(parser, help=help_msg, description=description)
        scp = self.subcommand_parser
        preprocess_register_args(scp)

    def __call__(self, args):
        preprocess_main(args)


class PostprocessModule(StammpModule):
    subcommand = 'postprocess'

    def __init__(self, parser):
        help_msg = 'run postprocessing pipeline'
        description = 'start postprocessing pipeline using a config script'
        super().__init__(parser, help=help_msg, description=description)
        scp = self.subcommand_parser
        postprocess_register_args(scp)

    def __call__(self, args):
        postprocess_main(args)


class CheckConfigModule(StammpModule):
    subcommand = 'check-config'

    def __init__(self, parser):
        help_msg = 'evaluate jinja code and print out the config file'
        description = help_msg
        super().__init__(parser, help=help_msg, description=description)
        scp = self.subcommand_parser
        check_cfg.register_arguments(scp)

    def __call__(self, args):
        check_cfg.run(args)


class CreateMetaMockModule(StammpModule):
    subcommand = 'create-metamock'

    def __init__(self, parser):
        help_msg = 'create a mock experiment from a set of CLIP datasets'
        description = help_msg
        super().__init__(parser, help=help_msg, description=description)
        scp = self.subcommand_parser
        metamock_cfg.register_arguments(scp)

    def __call__(self, args):
        metamock_cfg.run(args)


class SubsampleBamModule(StammpModule):
    subcommand = 'subsample-bam'

    def __init__(self, parser):
        help_msg = 'subsample a bam file roughly to a given number of reads'
        description = help_msg
        super().__init__(parser, help=help_msg, description=description)
        scp = self.subcommand_parser
        subsample_cfg.register_arguments(scp)

    def __call__(self, args):
        subsample_cfg.run(args)


class FlipMateModule(StammpModule):
    subcommand = 'flip-mate'

    def __init__(self, parser):
        help_msg = 'flip strand of second read'
        description = (
            'flip the strand of the second read. Used for generating a normalizing pileup from '
            'a paired-end sequenced library'
        )
        super().__init__(parser, help=help_msg, description=description)
        scp = self.subcommand_parser
        flipread2_register_args(scp)

    def __call__(self, args):
        flipread2_main(args)


if __name__ == '__main__':
    main()
