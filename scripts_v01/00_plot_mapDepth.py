import io
import re
import os
import sys
import json
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

class make_json:
    def __init__(self, infile, y, legend, res_var, data, json_name):
        self.infile = infile
        self.y = y
        self.legend = legend
        self.res_var = res_var
        self.data = data
        self.json_name = json_name

    def get_colorpalette(self):
        """legendとRBG値を対応させる関数

        Returns:
            rgb (list): RBG値が格納された配列
        """
        n_colors = len(self.legend)
        palette = sns.color_palette('hls', n_colors) #colorpalette = 'hls'
        rgb = ['rgb({},{},{})'.format(*[x*256 for x in rgb]) for rgb in palette]
        return rgb

    def data2json(self):
        """目的データをjson形式にする関数
        """
        logger.info('#--------------------------------------------------------------')
        logger.info('# output ...')
        logger.info('#--------------------------------------------------------------')
        #colors = self.get_colorpalette()
        y = {}
        df = input_xlsx(self.infile, self.legend)
        x = df.index.to_list()

        if self.y == '1':
            all(df, self.legend, self.data, self.res_var, y, x, self.json_name)
        elif self.y =='2':
            robustness(df, self.legend, self.data, self.res_var, y, x, self.json_name)

def input_xlsx(infile, legends):
    """xlsxデータ読み込み関数

    Returns:
        df (DataFrame): excelファイルの情報
    """
    logger.info('#--------------------------------------------------------------')
    logger.info('# read input file ...')
    logger.info('#--------------------------------------------------------------')

    if os.path.exists(infile) is False:
        logger.error(f'input file: {infile} not found..')
        logger.error('Pls cheak parameter \'-i \'.')
        sys.exit()

    logger.info(f'input file: {infile} found.')
    df = pd.read_excel(infile, engine='openpyxl', sheet_name=0)
    item = [f'{df.iloc[:, 0][i]}' for i in range(len(df.iloc[:, 0]))]

    flag = 0
    for legend in legends:
        if df.columns.str.contains(legend).any():
            flag = 1

    if flag == 1:
        df = df.sort_index(axis=1)
        df.index = item
    else:
        df = df.sort_index(axis=0).transpose()
        df.columns = item

    return df

def all(df, legend, data, res_var, y, x, json_name):
    """総当たりのデプス分布グラフ対応jsonファイル作成関数

    Args:
        df (DataFrame): excelファイルの情報
        legend (list): legend名が格納された配列
        data (list): データ値が格納された配列
        res_var (str): 応答変数
        y (dict): legendに対応したyの値
        x (list): xの値
        json_name (str): 出力するjsonファイル名
    """
    sample = ['tumor', 'normal']
    for ii in range(len(sample)):
        for i in range(len(legend)):
            for k, (tumor_data, normal_data) in enumerate(zip(data[0], data[1])):
                tumor_value = re.compile(r'(\d+)_(\d+)').search(tumor_data).group(1)
                normal_value = re.compile(r'(\d+)_(\d+)').search(tumor_data).group(2)
                if ii == 0:
                    y = extract_data(df, legend[i], tumor_data, normal_value, res_var, y)
                else:
                    y = extract_data(df, legend[i], normal_data, tumor_value, res_var, y)

            plot_data = pd.DataFrame()
            for name, yy in sorted(y.items(), key=lambda x:x[0]):
                col = pd.Series([name, list(map(float, yy)), list(map(int, x)), 'solid', None, 'tozeroy']) # fill='tozeroy'
                plot_data = plot_data.append(col, ignore_index=True)

            plot_data.columns = ['name', 'y', 'x', 'line_type', 'color', 'fill']
            outfile = json_name + '_' + legend[i] + '_' + sample[ii] + '.json'
            plot_data.to_json(outfile, orient='records', indent=4)
            logger.info(f'output json: {outfile}')

def robustness(df, legend, data, res_var, y, x, json_name):
    """室内再現性(Rep:reproducibility)または、
       頑健性(Rob:robustness)デプス分布グラフ対応jsonファイル作成関数

    Args:
        df (DataFrame): excelファイルの情報
        legend (list): legend名が格納された配列
        data (list): データ値が格納された配列
        res_var (str): 応答変数
        y (dict): legendに対応したyの値
        x (list): xの値
        json_name (str): 出力するjsonファイル名
    """
    for i in range(len(legend)):
        for j in range(len(data)):
            if re.fullmatch(r'(\d+)_(\d+)', data[j]):
                data_value = re.compile(r'(\d+)_(\d+)').search(data[j]).group(1)
            else:
                data_value = data[j]

            y = extract_data(df, legend[i], data[j], data_value, res_var, y)

        plot_data = pd.DataFrame()
        for name, yy in sorted(y.items(), key=lambda x:x[0]):
            col = pd.Series([name, list(map(float, yy)), list(map(int, x)), 'solid', None, 'tozeroy']) # fill='tozeroy'
            plot_data = plot_data.append(col, ignore_index=True)

        plot_data.columns = ['name', 'y', 'x', 'line_type', 'color', 'fill']
        outfile = json_name + '_' + legend[i] + '_RepOrRob' + '.json'
        plot_data.to_json(outfile, orient='records', indent=4)
        logger.info(f'output json: {outfile}')

def extract_data(df, legend, data, data_value, res_var, y):
    """目的データを抽出する関数

    Args:
        df (DataFrame): excelファイルの情報
        legend (list): legend名が格納された配列
        data (list): データ値が格納された配列
        data_value (str): データ値
        res_var (str): 応答変数
        y (dict): legendに対応したyの値

    Returns:
        y (dict): legendに対応したyの値
    """
    logger.info('#--------------------------------------------------------------')
    logger.info('# data processing ...')
    logger.info('#--------------------------------------------------------------')

    name = legend + '_' + data_value
    header = legend + '_' + data + res_var
    if header in df.columns:
        logger.info(f'data: {header} found.')
        y[name] = df.loc[:, header].to_list()
    else:
        logger.error(f'data: {header} not found..')
        logger.error('Pls cheak parameter \'-x1 or -x2\'.')
        sys.exit()

    return y

def parameters(__desc__):
    """入力されたコマンドライン引数をパースする関数

    Args:
        __desc__ (str): usage文

    Returns:
        args (argparse.Namespace): コマンドライン引数パース結果
    """
    parser = argparse.ArgumentParser(
        description = __desc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter #show default
    )
    parser.add_argument('-i',
                        dest='infile',
                        help='input file',
                        type=str,
                        required=True
                        )
    parser.add_argument('-l',
                        dest='legend',
                        help='input sample name based on the input data header',
                        nargs='+',
                        type=str,
                        required=True
                        )
    parser.add_argument('-y',
                        dest='y',
                        help='choose the following number' + '\n' +
                        '(1: all, 2: robustness)',
                        choices = ['1', '2'],
                        type=str,
                        required=True
                        )
    parser.add_argument('-x1',
                        dest='X1',
                        help='input tumor_data or x-axis based on the experiment',
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
    parser.add_argument('-o',
                        dest='out_dir',
                        help='input output path (absolution path)',
                        type=str,
                        default = os.path.abspath('.')
                        )
    parser.add_argument('-n',
                        dest='json_name',
                        help='input json filename',
                        default=os.path.join(os.path.abspath('.'), os.path.basename(sys.argv[0]).replace('.py', '') + '.json')
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

    infile = args.infile
    legend = sorted(args.legend)
    y = args.y
    tumor_data = args.X1
    normal_data = args.X2
    outdir = args.out_dir
    json_name = args.json_name
    res_var = '_%_Distribution'

    #===================   praparation   ===================#
    if normal_data is None:
        Data = list(map(str, sorted(tumor_data, key=int)))
    else:
        Data = []
        if y == '1':
            for td in sorted(tumor_data, key=int):
                for nd in sorted(normal_data, key=int):
                    data = str(td) + '_' + str(nd)
                    Data.append(data)
        else:
            tumorData = []
            normalData = []
            for td in sorted(tumor_data, key=int):
                for nd in sorted(normal_data, key=int):
                    data = str(td) + '_' + str(nd)
                    tumorData.append(data)
                    data = str(nd) + '_' + str(td)
                    normalData.append(data)

            Data = [tumorData, normalData]

    #===================   run process   ===================#
    jsonName = os.path.join(outdir, json_name + '_distDepth')
    make_json(infile, y, legend, res_var, Data, jsonName).data2json()

    #===================  move log file  ===================#
    logfile = jsonName + '_' + DATE + '_' + os.path.basename(sys.argv[0]).replace('.py', '') + '_log.txt'
    os.system(f'mv {LOGFILE} {logfile}')

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