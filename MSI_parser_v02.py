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

class make_table:
    def __init__(self, lotID, indir, outdir):
        self.lotID = lotID
        self.indir = indir
        self.outdir = outdir

    def get_variants(self):
        if os.path.exists(self.indir):
            logger.info(f'vcf_path: {self.indir}')

            snv = {}
            indel = {}
            for file in glob.glob(os.path.join(self.indir, '*log')):
                logger.info(f'vcf_file: {file}')
                file_name = file.split('/')
                ID = re.compile(r'(.*)_pass_(.*)\.log').search(file_name[-1]).group(1)
                vcf_type = re.compile(r'(.*)_pass_(.*)\.log').search(file_name[-1]).group(2)
                with open(file, 'r') as rf:
                    for line in rf.readlines():
                        line = line.rstrip('\n')
                        if 'Sites' in line:
                            variants = re.compile(r'kept (\d+) out').search(line).group(1)
                            if vcf_type ==  'snv':
                                snv[ID] = int(variants)
                            else:
                                indel[ID] = int(variants)
        else:
            logger.error(f'vcf_path: {self.indir} not found..')

        return snv, indel

    def output(self):
        snv, indel = self.get_variants()

        tmp_file = os.path.join(os.path.abspath('.'), 'tmp.txt')
        with open(tmp_file, 'w') as wf:
            wf.write('ID\tSNV\tINDEL\tSNV+INDEL\n')
            for ID in sorted(snv.keys()):
                sum = snv[ID] + indel[ID]
                wf.write(f'{ID}\t{snv[ID]}\t{indel[ID]}\t{sum}\n')

        xlsxname = self.lotID + '_snvIndel_summary.xlsx'
        xlsxfile = os.path.join(self.outdir, xlsxname)
        logger.info(f'result file: {xlsxfile}')

        wb = openpyxl.Workbook()
        sheet = wb.active
        sheet.title = 'result_' + DATE

        #tmp.txtの内容を張り出す
        with open(tmp_file, 'r') as rf:
            for line in rf.readlines():
                line = line.rstrip('\n')

        wb.save(xlsxfile)

        os.system(f'rm {tmp_file}')


class make_json:
    def __init__(self, infile, legend, data, json_name, n, flag):
        self.infile = infile
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

        colors = self.get_colorpalette()
        y = {}
        x = {}
        color = {}

        #===================   data processing   ===================#
        for i in range(len(self.infile)):
            if os.path.exists(self.infile[i]):
                logger.info(f'Input file: {self.infile[i]} found.')

                df = pd.read_excel(self.infile[i], engine='openpyxl')

                if self.flag == 0:
                    # for TMB
                    ID = df['Unnamed: 0']
                    index = [f'{ID[i]}' for i in range(len(ID))]
                    df = df.sort_index(axis=1)
                    df.index = index
                    response_variable = 'TMB (Mutations/Mb)' # indicate TMB value

                else:
                    # for MSI, snv+indel
                    ID = df['ID']
                    columns = [f'{ID[i]}' for i in range(len(df.index))]
                    df = df.sort_index(axis=0).transpose()
                    df.columns = columns
                    response_variable = '%'  # indicate MSI value

                index = df.index.get_loc(response_variable)

                for j in range(len(self.legend)):
                    if df.columns.str.contains(self.legend[j]).sum():
                        logger.info(f'header: {self.legend[j]} found in {self.infile[i]}')
                        logger.info('#======================================#')
                        logger.info(' Data processing.. ')
                        logger.info('#======================================#')
                        if self.flag == 0:
                            for tumor_data in self.data:
                                header = self.legend[j] + '_' + str(tumor_data) + '_10'
                                columns = df.columns.get_loc(header)

                                if i == 0:
                                    name = self.legend[j]
                                else:
                                    name = self.legend[j] + '_OBF'

                                if y.get(name):
                                    y[name] += ',' + str(df.iloc[index, columns])
                                    x[name] += ',' + str(tumor_data)
                                else:
                                    y[name] = str(df.iloc[index, columns])
                                    x[name] = str(tumor_data)

                                color[self.legend[j]] = colors[j]

                        else:
                            for k, (tumor_data, normal_data) in enumerate(zip(self.data[0], self.data[1])):
                                # tumor
                                tumor_header = self.legend[j] + '_' + tumor_data
                                tumor_columns = df.columns.get_loc(tumor_header)

                                tumor_value = re.compile(r'(.*)_(.*)').search(tumor_data).group(1)
                                normal_value = re.compile(r'(.*)_(.*)').search(tumor_data).group(2)

                                name = 'tumor,' + self.legend[j] + '_' + normal_value

                                if y.get(name):
                                    y[name] += ',' + str(df.iloc[index, tumor_columns])
                                    x[name] += ',' + str(tumor_value)
                                else:
                                    y[name] = str(df.iloc[index, tumor_columns])
                                    x[name] = str(tumor_value)
                                    color[name] = None

                                # normal
                                normal_header = self.legend[j] + '_' + normal_data
                                normal_columns = df.columns.get_loc(normal_header)

                                name = 'normal,' + self.legend[j] + '_' + tumor_value

                                if y.get(name):
                                    y[name] += ',' + str(df.iloc[index, normal_columns])
                                    x[name] += ',' + str(normal_value)
                                else:
                                    y[name] = str(df.iloc[index, normal_columns])
                                    x[name] = str(normal_value)
                                    color[name] = None

                    else:
                        logger.error(f'header: {self.legend[j]} not found in {self.infile[i]}')
                        logger.error('Pls cheak parameter \'-n\'.')
                        sys.exit()
            else:
                logger.error(f'Input file: {self.infile[i]} not found..')
                logger.error('Pls cheak parameter \'-i1 or -i2\'.')
                sys.exit()

        #===================   output   ===================#
        if self.flag == 0:
            plot_data = pd.DataFrame()
            for name in y.keys():
                Y = [float(i) for i in y[name].split(',')]
                X = [int(i) for i in x[name].split(',')]
                if 'OBF' in name:
                    name = re.compile(r'(.*)_OBF').search(name).group(1)
                    col = pd.Series([name, Y, X, 'dash', color[name], 'none'])
                else:
                    col = pd.Series([name, Y, X, 'solid', color[name], 'none'])

                plot_data = plot_data.append(col, ignore_index=True)

            plot_data.columns = ['name', 'y', 'x', 'line_type', 'color', 'fill']

            outfile = self.json_name + '.json'
            plot_data.to_json(outfile, orient='records', indent=4)
            logger.info(f'Output json: {outfile}')

        else:
            samples = ['tumor', 'normal']
            for ii in range(len(samples)):
                plot_data = pd.DataFrame()
                for name in y.keys():
                    if samples[ii] in name:
                        rename = re.compile(r'(.*),(.*)').search(name).group(2)
                        Y = [float(i) for i in y[name].split(',')]
                        X = [int(i) for i in x[name].split(',')]
                        col = pd.Series([rename, Y, X, 'solid', color[name], 'none'])
                        plot_data = plot_data.append(col, ignore_index=True)

                plot_data.columns = ['name', 'y', 'x', 'line_type', 'color', 'fill']

                outfile = self.json_name + '_' + samples[ii] + '.json'
                plot_data.to_json(outfile, orient='records', indent=4)
                logger.info(f'Output json: {outfile}')


def parameters(__desc__):
    parser = argparse.ArgumentParser(
        description = __desc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter #show default
    )
    parser.add_argument('-i1',
                        dest='infile1',
                        help='input normal TMB file ',
                        type=str,
                        required=True
                        )
    parser.add_argument('-i2',
                        dest='infile2',
                        help='input orientation bias filtered TMB file',
                        type=str,
                        default=None
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

    infile1 = args.infile1
    infile2 = args.infile2
    legend = args.legend
    tumor_data = args.X1
    normal_data = args.X2
    json_name = args.json_name

    # run process
    # must define how to use this program again!!!!!

    if infile2 is None:
        infile = [infile1]
    else:
        infile = [infile1, infile2]


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

    make_json(infile, legend, data, json_name, n, flag).data2json()

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