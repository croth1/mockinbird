import argparse

from mockinbird.utils.argparse_helper import file_r

from jinja2 import Environment, PackageLoader
import jinja2

def register_arguments(parser):
    parser.add_argument('config_file', type=file_r, help='path to config file to be checked')


def create_parser():
    description = 'evaluate jinja code and print out the config file'
    parser = argparse.ArgumentParser(prog='mockinbird-check-config', description=description)
    register_arguments(parser)
    return parser


def run(args):
    env = Environment(loader=PackageLoader('mockinbird'))
    with open(args.config_file) as infile:
        try:
            data = env.from_string(infile.read())
            print(data.render())
        except jinja2.exceptions.TemplateSyntaxError as e:
            print('Error in line %s: %s' %( e.lineno, e))
        except jinja2.exceptions.TemplateRuntimeError as e:
            print('Error while evaluating the jinja code:', e)


def run_cmdline():
    parser = create_parser()
    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    run_cmdline()
