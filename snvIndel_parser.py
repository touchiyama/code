import io
import os
import re
import sys
import glob
import json
import openpyxl
import datetime
import argparse
import pandas as pd
import seaborn as sns
import logging
from logging import Logger, getLogger, Formatter, StreamHandler, FileHandler


DATE = datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
LOGFILE = os.path.join(os.path.abspath('.'),
                       os.path.basename(sys.argv[0]).replace('.py', '') + '_' + DATE + '_log.txt')
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

#===================   praparation   ===================#
# 以下のclassをrun_vcftools_indelSnv.pyに移行させるか?

class make_json:
    def __init__(self, indir, legend, data, json_name, n, flag):
        self.indir = indir
        self.legend = legend
        self.data = data
        self.json_name = json_name
        self.n = n
        self.flag = flag

    def get_colorpalette(self):
        n_colors = len(self.n)
        palette = sns.color_palette('hls', n_colors) #colorpalette = 'hls'
        rgb = ['rgb({},{},{})'.format(*[x*256 for x in rgb]) for rgb in palette]
        return rgb

    def data2json(self):
        if os.path.exists(self.indir):
            logger.info(f'vcf_path: {self.indir}')

            colors = self.get_colorpalette()
            legend = {}
            snv = {}
            indel = {}

            #===================   get variants   ===================#
            for file in glob.glob(os.path.join(self.indir, '*log')):
                logger.info(f'vcf_file: {file}')
                file_name = file.split('/')
                ID = re.compile(r'(.*)_pass_(.*)\.log').search(file_name[-1]).group(1)
                vcf_type = re.compile(r'(.*)_pass_(.*)\.log').search(file_name[-1]).group(2)
                name = re.compile(r'(.*)_(\d+)').search(ID).group(1)
                if legend.get(name):
                    legend[name] += ',' + ID
                else:
                    legend[name] = ID

                with open(file, 'r') as rf:
                    for line in rf.readlines():
                        line = line.rstrip('\n')
                        if 'Sites' in line:
                            variants = re.compile(r'kept (\d+) out').search(line).group(1)
                            if vcf_type ==  'snv':
                                snv[ID] = int(variants)
                            else:
                                indel[ID] = int(variants)

            #===================   output   ===================#

            samples = ['snv', 'indel', 'snvIndel']

            for i in range(len(samples)):
                plot_data = pd.DataFrame()
                for j, (name, value) in enumerate(legend.items()):
                    X = []
                    Y = []
                    for ID in sorted(set([i for i in value.split(',')])):
                        x = re.compile(r'(.*)_(\d+)').search(ID).group(2)
                        X.append(int(x))
                        if samples[i] == 'snv':
                            Y.append(int(snv[ID]))
                            filename = self.json_name + '_snv'
                        elif samples[i] == 'indel':
                            Y.append(int(indel[ID]))
                            filename = self.json_name + '_indel'
                        else:
                            Y.append((int(snv[ID]) + int(indel[ID])))
                            filename = self.json_name + '_snvIndel'

                col = pd.Series([name, Y, X, 'solid', colors[j], 'none'])
                plot_data = plot_data.append(col, ignore_index=True)

                plot_data.columns = ['name', 'y', 'x', 'line_type', 'color', 'fill']

                outfile = filename + '.json'
                plot_data.to_json(outfile, orient='records', indent=4)
                logger.info(f'Output json: {outfile}')

        else:
            logger.error(f'vcf_path: {self.indir} not found..')


def parameters(__desc__):
    parser = argparse.ArgumentParser(
        description = __desc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter #show default
    )
    parser.add_argument('-i',
                        dest='indir',
                        help='input vcf_path (abusolute path) ',
                        type=str,
                        required=True
                        )
    parser.add_argument('-n',
                        dest='legend',
                        help='input sample name based on the input data header',
                        nargs='+',
                        type=str,
                        required=True
                        )
    parser.add_argument('-x1',
                        dest='X1',
                        help='tumor_data or x-axis based on the experiment',
                        nargs='+',
                        type=int,
                        required=True
                        )
    parser.add_argument('-x2',
                        dest='X2',
                        help='input normal data based on the experiment',
                        nargs='+',
                        type=int
                        )
    parser.add_argument('-j',
                        dest='json_name',
                        help='input json filename',
                        default=os.path.join(os.path.abspath('.'), os.path.basename(sys.argv[0]).replace('.py', ''))
                        )
    parser.add_argument('-log',
                        dest='loglevel',
                        help='choose log level (default: INFO)',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        default='INFO'
                        )
    args = parser.parse_args()

    return args

def show_param(args: argparse.Namespace):
    """入力として与えられたコマンドライン引数を表示

    Args:
        args (argparse.Namespace): コマンドライン引数パース結果
    """
    logger.info("--------------------------------------------------------------")
    logger.info(f"program: {__file__}")
    cla = vars(args)
    for key in sorted(cla.keys()):
        logger.info(f"{key}: {cla[key]}")
    logger.info("--------------------------------------------------------------")
    logger.info("")

def main(args):

    logger.info('')
    logger.info('Output log..')

    # set up parameters
    show_param(args)

    indir = args.indir
    legend = args.legend
    tumor_data = args.X1
    normal_data = args.X2
    json_name = args.json_name

    # run process
    # must define how to use this program again!!!!!
    if normal_data is None:
        n = len(legend)
        data = tumor_data
        flag = 0
    else:
        tumorData = []
        normalData = []
        for td in sorted(tumor_data, key=int):
            for nd in sorted(normal_data, key=int):
                data = str(td) + '_' + str(nd)
                tumorData.append(data)
                data = str(nd) + '_' + str(td)
                normalData.append(data)

        n = len(tumorData)
        data = [tumorData, normalData]
        flag = 1

    make_json(indir, legend, data, json_name, n, flag).data2json()

    logger.info('Done!')
    logger.info('')


if __name__ == '__main__':
    __version__ = '1.0'
    __desciption__ = 'Some useful program commands are:'

    parm = parameters(__desciption__)

    logger = getLogger(__name__)
    logger.setLevel(parm.loglevel)

    FORMAT = '%(levelname)s:[%(asctime)s] %(message)s'
    #dt_fmt = '%Y-%m-%d %H:%M:%S'
    formatter = Formatter(FORMAT)

    stream_handler = StreamHandler()
    stream_handler.setLevel(parm.loglevel)
    stream_handler.setFormatter(formatter)

    file_handler = FileHandler(filename=LOGFILE, mode='w', encoding='utf-8')
    file_handler.setLevel(parm.loglevel)
    file_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    main(parm)

#
# END
#