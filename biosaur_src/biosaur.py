from __future__ import division
import argparse
import logging
import json_logging
from . import bio
import pyfiglet
import json
import coloredlogs
from termcolor import colored

json_logging.ENABLE_JSON_LOGGING = True
json_logging.init_non_web()
coloredlogs.install()


def run():

    logging.info(u'Starting program with following params...')
    parser = argparse.ArgumentParser(
        description='Detection of peptide features',
        epilog='''

    Example usage
    -------------
    $ biosaur input.mzml -mass_accuracy 7 -min_length 3
    -------------
    ''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input_mzml_path', nargs='+', help='input MZML file')
    parser.add_argument(
        '-np',
        '--number_of_processes',
        help='Number of processes',
        default=0)

    parser.add_argument(
        '-cm',
        '--correlation_map',
        help='Add correlation map to final table',
        action='store_true')

    parser.add_argument(
        '-nm',
        '--negative_mode',
        help='Add negative mode option',
        action='store_true')

    parser.add_argument(
        '-ac',
        '--mass_accuracy',
        help='Mass accuracy',
        type=float,
        default=8)
    parser.add_argument(
        '-minc',
        '--min_charge',
        help='Minimum charge',
        type=int,
        default=1)
    parser.add_argument(
        '-maxc',
        '--max_charge',
        help='Maximum charge',
        type=int,
        default=6)
    parser.add_argument(
        '-minl',
        '--min_length',
        help='Minimum length',
        type=int,
        default=3)
    parser.add_argument(
        '-minlh',
        '--min_length_hill',
        help='Minimum length for hill',
        type=int,
        default=2)
    parser.add_argument(
        '-mini',
        '--min_intensity',
        type=float,
        help='Minimum intensity',
        default=1)
    parser.add_argument(
        '-hvf',
        '--hill_valley_factor',
        help='Hill Valley Factor',
        type=float,
        default=1.3)
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Enable debugging output')
    parser.add_argument(
        '-out',
        '--output_file',
        help='Output File')
    parser.add_argument(
        '-pxfp',
        '--pep_xml_file_path',
        default='0',
        help='Pepxml filepath')
    # parser.add_argument(
        # '-faims',
        # '--faims',
        # help='Use when mzML contain FAIMS data',
        # action='store_true')
    args = vars(parser.parse_args())
    log_args = json.dumps(args, indent=2, sort_keys=False)
    logging.basicConfig(format='%(levelname)9s: %(asctime)s %(message)s',
                        datefmt='[%H:%M:%S]',
                        level=[logging.INFO, logging.DEBUG][args['debug']])

    logging.info('Starting with args: %s', log_args)
    return bio.process_files(args)


if __name__ == '__main__':
    run()
