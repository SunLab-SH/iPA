# coding = 'utf-8'

import logging
import argparse 


LOG = logging.getLogger('main')
parser = argparse.ArgumentParser(description='ET data processing')
## load data
parser.add_argument('--root-dir', type=str, default='E:\\filament_project_weimin\\data', help='root for all files')
parser.add_argument('--data-root-dir', type=str, default='Segmentation\\raw_data', help = 'root path for all ET data' )
parser.add_argument('--output-dir', type=str, default=f'results', help='output root for all visualization data')
parser.add_argument('--parameter-dir', type=str, default='parameters', help='root for all parameters')
parser.add_argument('--processed-data-dir', type=str, default=f'Segmentation\\processed', help='output root for all visualization data' )
parser.add_argument('--fa-raw-data-dir', type=str, default=f'E:\\filament_project_weimin\\data\\folcon\\FA', help='raw data root facal adhesion data' )
parser.add_argument('--ca-raw-data-dir', type=str, default=f'E:\\filament_project_weimin\\data\\folcon\\Ca', help='raw data root calcium data' )

# processed data
parser.add_argument('--actin-file', default="*Actin.xml")
parser.add_argument('--mt-file', default="*MT.xml")
parser.add_argument('--isg-file', default='*ISG.mrc')
parser.add_argument('--mito-file', default='*mito_filled.mrc')
parser.add_argument('--actin-modified-file', default="*filament_filled_points.json")
parser.add_argument('--mt-modified-file', default="*MT_filled_points.json")
# parser.add_argument('--ca-modified-file', default="*focal_to_ca_distance.csv")
# multi threading
parser.add_argument('--thread-num',default=4)



arg = parser.parse_args(args=[])