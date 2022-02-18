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

class data2json:
    def __init__(self, infile, legend, res_var, data, json_name, color, n):
        self.infile = infile
        self.legend = legend
        self.res_var = res_var
        self.data = data
        self.json_name = json_name
        self.color = color
        self.n = n

    def get_colorpalette(self, color, n_colors):
        """legendとRBG値を対応させるための関数

        Returns:
            rgb (list): RBG値が格納された配列
        """
        palette = sns.color_palette(color, n_colors) #colorpalette = 'hls'
        rgb = ['rgb({},{},{})'.format(*[x*256 for x in rgb]) for rgb in palette]
        return rgb

    def data_processing(self):
        """データを抽出する関数

        Returns:
            y(dict): legendに対応したy
            x(dict): legendに対応したx
            color(dict): legendに対応したRBG値
        """
        y = {}
        x = {}
        color = {}
        colors = self.get_colorpalette(self.color, self.n)
        #===================   data processing   ===================#
        df = input_xlsx(self.infile, self.legend)
        if self.res_var in df.index:
            index = df.index.get_loc(self.res_var)
            logger.error(f'index: {self.res_var} found.')
        else:
            logger.error(f'input file {self.infile} was incorrect..')
            logger.error('Pls cheak parameter \'-i\'.')
            sys.exit()

        for j in range(len(self.legend)):
            if df.columns.str.contains(self.legend[j]).any() is False:
                logger.error(f'legend: {self.legend[j]} not found..')
                logger.error('Pls cheak parameter \'-n\'.')
                sys.exit()

            logger.info(f'legend: {self.legend[j]} found.')

            logger.info('#--------------------------------------------------------------')
            logger.info('# data processing ...')
            logger.info('#--------------------------------------------------------------')

            for k, (tumor_data, normal_data) in enumerate(zip(self.data[0], self.data[1])):
                tumor_value = re.compile(r'(\d+)_(\d+)').search(tumor_data).group(1)
                normal_value = re.compile(r'(\d+)_(\d+)').search(tumor_data).group(2)
                #===================   tumor data   ===================#
                y, x, color = extract_data(
                    df, index, self.legend[j], tumor_data, normal_value,
                    'tumor', y, x, color, colors[k]
                )
                #===================  normal data   ===================#
                y, x, color = extract_data(
                    df, index, self.legend[j], normal_data, tumor_value,
                    'normal', y, x, color, colors[k]
                )

        return y, x, color

    def output_json(self):
        """jsonファイル出力関数
        """
        y, x, color = self.data_processing()

        logger.info('#--------------------------------------------------------------')
        logger.info('# output ...')
        logger.info('#--------------------------------------------------------------')

        #===================   output   ===================#
        sample = ['tumor', 'normal']
        for ii in range(len(sample)):
            plot_data = pd.DataFrame()
            for name in sorted(y.keys()):
                if sample[ii] in name:
                    rename = re.compile(r'(.*),(.*)').search(name).group(2)
                    Y = [float(i) for i in y[name].split(',')]
                    X = [int(i) for i in x[name].split(',')]
                    col = pd.Series([rename, Y, X, 'solid', color[name], 'none'])

                plot_data = plot_data.append(col, ignore_index=True)

            plot_data.columns = ['name', 'y', 'x', 'line_type', 'color', 'fill']

            type = re.compile(r'(.*)_(.*)_(\d+)').search(rename).group(2)
            outfile = self.json_name + '_' + type + '_' + sample[ii] + '.json'
            plot_data.to_json(outfile, orient='records', indent=4)
            logger.info(f'output json: {outfile}')

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

def extract_data(df, index, legend, data, data_value, sample, y, x, color, color_type):
    """目的データを抽出する関数

    Args:
        df (DataFrame): excelファイルの情報
        index (int): データフレーム座標上で応答変数に対応するindex
        legend (list): legend名が格納された配列
        data (list): データ値が格納された配列
        data_value (str): データ値
        sample (str): サンプル名
        y (dict): legendに対応したyの値
        x (dict): legendに対応したxの値
        color (dict): legendに対応したRBG値

    Returns:
        y (dict): legendに対応したy
        x (dict): legendに対応したx
        color (dict): legendに対応したRBG値
    """
    header = legend + '_' + data
    if header in df.columns:
        logger.info(f'data: {header} found.')
        columns = df.columns.get_loc(header)

        name = sample + ',' + legend + '_' + data_value

        if y.get(name):
            y[name] += ',' + str(df.iloc[index, columns])
            x[name] += ',' + str(data_value)
        else:
            y[name] = str(df.iloc[index, columns])
            x[name] = str(data_value)
            color[name] = color_type
    else:
        logger.error(f'data: {header} not found..')
        logger.error('Pls cheak parameter \'-x1 or -x2\'.')
        sys.exit()

    return y, x, color

def display_color():
    color_list = [
        'Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r',
        'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r','Dark2', 'Dark2_r', 'GnBu', 'GnBu_r',
        'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r',
        'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r',
        'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r',
        'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy','RdGy_r',
        'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r',
        'Set1', 'Set1_r','Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r',
        'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r',
        'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r',
        'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r',
        'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'crest', 'crest_r',
        'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'flare', 'flare_r', 'gist_earth', 'gist_earth_r',
        'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r','gist_rainbow', 'gist_rainbow_r',
        'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r',
        'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'icefire', 'icefire_r',
        'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'mako', 'mako_r',
        'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r',
        'prism', 'prism_r', 'rainbow', 'rainbow_r', 'rocket', 'rocket_r', 'seismic', 'seismic_r',
        'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r',
        'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'turbo', 'turbo_r',
        'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'vlag', 'vlag_r',
        'winter', 'winter_r'
    ]
    return color_list


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
                        help='input summary xlsx file',
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
                        '(1: mapping-related, 2: TMB, 3: MSI, 4: SNV+INDEL, 5: sensitivity, 6: CNV, 7: CNV_Corr)' ,
                        choices = ['1', '2', '3', '4', '5', '6', '7'],
                        type=str,
                        required=True
                        )
    parser.add_argument('-x1',
                        dest='X1',
                        help='input tumor_data based on the experiment',
                        nargs='+',
                        type=int,
                        required=True
                        )
    parser.add_argument('-x2',
                        dest='X2',
                        help='input normal data based on the experiment',
                        nargs='+',
                        type=int,
                        required=True
                        )
    parser.add_argument('-o',
                        dest='out_dir',
                        help='input output path (absolution path)',
                        type=str,
                        default = os.path.abspath('.')
                        )
    parser.add_argument('-n',
                        dest='json_name',
                        help='input json filename  (You should include LotID)',
                        default=os.path.join(os.path.abspath('.'), os.path.basename(sys.argv[0]).replace('.py', ''))
                        )
    parser.add_argument('-c',
                        dest='color',
                        help='choose the following colors',
                        choices=display_color(),
                        type=str,
                        default = 'hls'
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

def response_variable(args):
    """変数対応ハッシュテーブル（逐次変更）

    Args:
        args (str): 与えらえた数値

    Returns:
        var[args] (dict): 変数に対応した値
    """
    var = { '1': [['Mapped reads (%)',
                   'Number of duplicate marked reads (%)',
                   'Insert length: mean',
                   'Aligned reads in target region (%)',
                   'Average alignment coverage over target region'],
                  ['mapping_rate',
                   'duplicate_reads',
                   'insert_length',
                   'On-target',
                   'depth']],
            '2': [['TMB (Mutations/Mb)'], ['TMB']],
            '3': [['%'], ['MSI']],
            '4': [['SNV', 'INDEL', 'SNV+INDEL'], ['SNV', 'INDEL', 'SnvIndel']],
            '5': [['sensitivity(%)'], ['Se']],
            '6': [['Total CNV (PASS)'], ['CNV']],
            '7': [['corr(int)', 'corr(float)'], ['CnvCorr_int', 'CnvCorr_flt']]
          }

    return var[args]

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
    color = args.color
    res_var = response_variable(y)

    #===================   praparation   ===================#
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

    #===================   run process   ===================#
    for i, rv in enumerate(res_var[0]):
        jsonName = os.path.join(outdir, json_name + '_' + res_var[1][i])
        data2json(infile, legend, rv, data, jsonName, color, n).output_json()

    logger.info('Done!')
    logger.info('')

    #===================  move log file  ===================#
    logfile = jsonName + '_' + DATE + '_' + os.path.basename(sys.argv[0]).replace('.py', '') + '_log.txt'
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