import io
import os
import re
import sys
import glob
import openpyxl
import datetime
import argparse
import pandas as pd
import seaborn as sns
import logging
from logging import Logger, getLogger, Formatter, StreamHandler, FileHandler

DATE = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
LOGFILE = os.path.join(os.path.abspath("."),
                       os.path.basename(sys.argv[0]).replace(".py", "") + "_" + DATE + "_log.txt")
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def get_colorpalette(n_colors):
    """legendとRBG値を対応させる関数

    Returns:
        rgb (list): RBG値が格納された配列
    """
    palette = sns.color_palette('hls', n_colors) #colorpalette = 'hls'
    rgb = ['rgb({},{},{})'.format(*[x*256 for x in rgb]) for rgb in palette]
    return rgb

class data2json:
    def __init__(self, indir, legend, gene, data, outdir, json_name):
        self.indir = indir
        self.legend = legend
        self.gene = gene
        self.data = data
        self.outdir = outdir
        self.json_name = json_name

    def read_cnvTable(self):
        """CNV variant summary xlsxからのデータを抽出する関数

        Returns:
            copy_num　(dict): サンプルに対応したCNV数
            group (dict): サンプルに対応したID群
        """
        logger.info('#--------------------------------------------------------------')
        logger.info('# data processing ...')
        logger.info('#--------------------------------------------------------------')

        markerGene = {}
        for gene in self.gene:
            markerGene[gene] = '-'

        copy_num = {}
        group = {}
        data_list = ','.join(self.data)

        for legend in self.legend:
            files = os.path.join(self.indir, '*CNV*' + legend + '*xlsx')

            for file in sorted(glob.glob(files)):
                path, df = input_xlsx(file)
                name = path.split('/')
                #A07_CNV25_TN_20_cnv_variant_summary_1.xlsx
                CNV_data = re.compile(r'(.*)_(.*)_(.*)_(.*)_cnv_(.*)').search(name[-1]).group(2) # will change
                sample = re.compile(r'(.*)_(.*)_(.*)_(.*)_cnv_(.*)').search(name[-1]).group(3)
                data = re.compile(r'(.*)_(.*)_(.*)_(.*)_cnv_(.*)').search(name[-1]).group(4)

                if data in data_list:
                    logger.info(f'xlsx file: {file} found.')
                    sample_data = sample + '_' + data

                    for gene in self.gene:
                        groupID = gene + '_' + sample
                        ID = gene + '_' + sample_data + ',' + CNV_data
                        if group.get(groupID):
                            group[groupID] += '#' + ID
                        else:
                            group[groupID] = ID

                    for i in range(2, len(df)):
                        Chr = str(df.loc[i, 'Chr'])
                        if (Chr == 'X') | (Chr == 'Y'):
                            pass
                        else:
                            if df.loc[i, 'FILTER'] == 'PASS':
                                if markerGene.get(df.loc[i, 'Gene name']):
                                    ID = df.loc[i, 'Gene name'] + '_' + sample_data + ',' + CNV_data
                                    copy_num[ID] = df.loc[i, 'CN']

        return copy_num, group

    def prep_output(self):
        """出力への準備する関数

        Returns:
            jsonfile (dict): jsonファイルに記載するサンプル群
            y (dict): サンプルに対応したy
            x (dict): サンプルに対応したx
        """
        copy_num, group = self.read_cnvTable()

        logger.info('#--------------------------------------------------------------')
        logger.info('# preparation for output ...')
        logger.info('#--------------------------------------------------------------')

        #============= preparation for jsonfile ================#
        sample = {}
        jsonfile = {}
        for name, IDs in sorted(group.items(), key=lambda x:x[0]):
            sample_groups = []
            for ID in IDs.split('#'):
                #============== convert N.A to 2 ================#
                if copy_num.get(ID) is None:
                    copy_num[ID] = 2

                #================== regrouping ==================#
                sample_group = re.compile(r'(.*),(.*)').search(ID).group(1)
                cnv_data = re.compile(r'(.*),CNV(\d+)').search(ID).group(2)
                if sample.get(sample_group):
                    sample[sample_group] += ',' + cnv_data
                else:
                    sample[sample_group] = cnv_data

                sample_groups.append(sample_group)

            jsonfile[name] = sorted(list(set(sample_groups)))

        #============= preparation for y and x ==============#
        y = {}
        x = {}
        for i, j in sample.items():
            for dt in sorted(j.split(','), key=int):
                gene = re.compile(r'(.*)_(.*)_(\d+)').search(i).group(1)
                samp = re.compile(r'(.*)_(.*)_(\d+)').search(i).group(2)
                name = gene + '_' + samp
                y_id = i + ',CNV' + dt
                x_id = gene + '_CNV' + dt
                if y.get(i):
                    y[i] += ',' + str(copy_num[y_id])
                    x[i] += ',' + str(known_copynum(x_id))
                else:
                    y[i] = str(copy_num[y_id])
                    x[i] = str(known_copynum(x_id))

        return jsonfile, y, x

    def output_json(self):
        """jsonファイル出力関数
        """
        jsonfile, y, x = self.prep_output()

        logger.info('#--------------------------------------------------------------')
        logger.info('# output ...')
        logger.info('#--------------------------------------------------------------')

        #===================   output   ===================#
        for name, IDs in sorted(jsonfile.items()):
            plot_data = pd.DataFrame()
            for ID in IDs:
                Y = [int(i) for i in y[ID].split(',')]
                X = [float(i) for i in x[ID].split(',')]
                col = pd.Series([ID, Y, X, 'dash', None, 'none'])
                plot_data = plot_data.append(col, ignore_index=True)

            plot_data.columns = ['name', 'y', 'x', 'line_type', 'color', 'fill']

            outfile = os.path.join(self.outdir, self.json_name + '_' + name + '_linearity' + '.json')
            plot_data.to_json(outfile, orient='records', indent=4)
            logger.info(f'output json: {outfile}')

def input_xlsx(infile):
    """xlsxデータ読み込み関数

    Returns:
        path (str): 拡張子を除いたファイル名
        df (DataFrame): excelファイルの情報
    """
    logger.info('#--------------------------------------------------------------')
    logger.info('# read input file ...')
    logger.info('#--------------------------------------------------------------')

    if len(glob.glob(infile)) == 0:
        logger.error(f'xlsx file: {infile} not found..')
        logger.error('Pls cheak parameter \'-i \'.')
        sys.exit()

    path, ext = os.path.splitext(infile)

    #logger.info(f'xlsx file: {infile} found.')
    df = pd.read_excel(infile, engine='openpyxl', sheet_name=0)
    columns = df.loc[1:2].values[0]
    df.columns = columns

    return path, df

def known_copynum(arg):
    """変数対応ハッシュテーブル（逐次変更）

    Args:
        args (str): 与えらえた数値

    Returns:
        var[arg] (dict): 変数に対応した値
    """
    var = {
            'MYCN_CNV0'   : 2.0,
            'MYCN_CNV25'  : 3.625,
            'MYCN_CNV50'  : 5.25,
            'MYCN_CNV75'  : 6.875,
            'MYCN_CNV100' : 8.5,
            'MET_CNV0'    : 2.0,
            'MET_CNV25'   : 2.625,
            'MET_CNV50'   : 3.25,
            'MET_CNV75'   : 3.875,
            'MET_CNV100'  : 4.5
          }

    return var[arg]

def parameters(__desc__):
    """入力されたコマンドライン引数をパースする関数

    Args:
        __desc__ (str): usage文

    Returns:
        args (argparse.Namespace): コマンドライン引数パース結果
    """
    parser = argparse.ArgumentParser(
        description = __desc__,
        #formatter_class=argparse.ArgumentDefaultsHelpFormatter #show default
    )
    parser.add_argument('-i',
                        dest='indir',
                        help='input cnv summary variant path (absolute path) ',
                        type=str,
                        required=True
                        )
    parser.add_argument('-l',
                        dest='legend',
                        help='input the following sample names' + '\n' +
                        '(TN: Tumor-Normal, TO: Tumor-Only, Both: TN+TO) (default: Both)' ,
                        choices = ['TN', 'TO', 'Both'],
                        type=str,
                        default = 'Both'
                        )
    parser.add_argument('-g',
                        dest='gene',
                        help='input gene name' + '\n' +
                        '(default: ' + 'MYCN, ' + 'MET' + ')',
                        nargs='*',
                        default=['MYCN', 'MET']
                        )
    parser.add_argument('-x1',
                        dest='X1',
                        help='input tumor_data or x-axis based on the experiment',
                        nargs='+',
                        type=int,
                        default=None,
                        required=True
                        )
    parser.add_argument('-x2',
                        dest='X2',
                        help='input normal data based on the experiment',
                        nargs='+',
                        type=int,
                        default=None
                        )
    parser.add_argument('-o',
                        dest='out_dir',
                        help='input output path (absolution path)',
                        type=str,
                        default = os.path.abspath('.')
                        )
    parser.add_argument('-n',
                        dest='json_name',
                        help='input json filename (default: ' + os.path.join(os.path.abspath('.'), os.path.basename(sys.argv[0]).replace('.py', '') + '.json') + ')' \
                            + '(You should include LotID)',
                        default=None
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
    logger.info("#--------------------------------------------------------------")
    logger.info(f"program: {__file__}")
    cla = vars(args)
    for key in sorted(cla.keys()):
        logger.info(f"{key}: {cla[key]}")
    logger.info("#--------------------------------------------------------------")
    logger.info("")

def main(args):

    logger.info('')
    logger.info('Output log..')

    #=================== setUp parameter ===================#
    show_param(args)

    indir = args.indir
    legend = args.legend
    gene = args.gene
    tumor_data = args.X1
    normal_data = args.X2
    outdir = args.out_dir
    json_name = args.json_name

    #===================   praparation   ===================#
    if legend == 'Both':
        leg_list = ['TN', 'TO']
    else:
        leg_list = [legend]

    if gene is None:
        gene = ['MYCN', 'MET']

    if normal_data is None:
        Data = list(map(str, sorted(tumor_data, key=int)))
    else:
        Data = []
        for x1 in sorted(tumor_data, key=int):
            for x2 in sorted(normal_data, key=int):
                data = str(x1) + '_' + str(x2)
                Data.append(data)
    

    #===================   run process   ===================#
    data2json(indir, leg_list, gene, Data, outdir, json_name).output_json()

    logger.info('Done!')
    logger.info('')

    #===================   move file     ===================#
    filename = json_name + '_' + DATE + '_' + os.path.basename(sys.argv[0]).replace('.py', '') + '_log.txt' 
    logfile = os.path.join(outdir, filename)
    os.system(f'mv {LOGFILE} {logfile}')

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