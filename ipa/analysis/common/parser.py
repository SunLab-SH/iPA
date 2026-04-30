# coding = 'utf-8'

import logging
import argparse 
import os

# Dynamically calculate project root directory for portability
# Assuming this file is at: ipa/analysis/common/parser.py
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(CURRENT_DIR, '..', '..', '..'))
DATA_ROOT = os.path.join(PROJECT_ROOT, 'data')

LOG = logging.getLogger('main')
parser = argparse.ArgumentParser(description='ET data processing')

## load data - Using relative paths based on PROJECT_ROOT
parser.add_argument('--root-dir', type=str, default=PROJECT_ROOT, help='root for all files')
parser.add_argument('--data-root-dir', type=str, default=os.path.join(DATA_ROOT, 'cryoET'), help = 'root path for all ET data' )
parser.add_argument('--output-dir', type=str, default='results', help='output root for all visualization data')
parser.add_argument('--parameter-dir', type=str, default='parameters', help='root for all parameters')
parser.add_argument('--processed-data-dir', type=str, default=os.path.join(DATA_ROOT, 'processed'), help='output root for all visualization data' )

# Raw data directories (Relative to DATA_ROOT)
parser.add_argument('--fa-raw-data-dir', type=str, default=os.path.join(DATA_ROOT, 'folcon', 'FA'), help='raw data root facal adhesion data' )
parser.add_argument('--ca-raw-data-dir', type=str, default=os.path.join(DATA_ROOT, 'folcon', 'Ca'), help='raw data root calcium data' )

# processed data patterns
parser.add_argument('--actin-file', default="*Actin.xml")
parser.add_argument('--mt-file', default="*MT.xml")
parser.add_argument('--isg-file', default='*ISG.mrc')
parser.add_argument('--mito-file', default='*mito_filled.mrc')
parser.add_argument('--actin-modified-file', default="*filament_filled_points.json")
parser.add_argument('--mt-modified-file', default="*MT_filled_points.json")

# multi threading
parser.add_argument('--thread-num', default=4)



arg = parser.parse_args(args=[])