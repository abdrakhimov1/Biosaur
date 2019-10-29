from __future__ import division
import argparse
import logging
import pkg_resources
from . import bio

def run():
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
    parser.add_argument('-ac', '--mass_accuracy', help='Mass accuracy', default=8)
    parser.add_argument('-minc', '--min_charge', help='Minimum charge', default=1)
    parser.add_argument('-maxc', '--max_charge', help='Maximum charge', default=6)
    parser.add_argument('-minl', '--min_length', help='Minimum length', default=3)
    parser.add_argument('-mini', '--min_intensity', help='Minimum intensity', default=1)
    parser.add_argument('-hvf', '--hill_valley_factor', help='Hill Valley Factor', default=0.8)
    parser.add_argument('--debug', action='store_true', help='Enable debugging output')
    parser.add_argument('-out', '--output_file', help='Output File', default = 'test_custom.features.tsv')
    args = vars(parser.parse_args())
    
    logging.basicConfig(format='%(levelname)9s: %(asctime)s %(message)s',
            datefmt='[%H:%M:%S]', level=[logging.INFO, logging.DEBUG][args['debug']])
    logger = logging.getLogger(__name__)
    logger.debug('Starting with args: %s', args)
    return bio.process_files(args)


if __name__ == '__main__':
    run()
