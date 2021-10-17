import io
import os
import re
import sys
import glob
import datetime
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from sklearn.linear_model import LinearRegression
from logging import Logger, getLogger, Formatter, StreamHandler, FileHandler

DATE = datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
LOGFILE = os.path.join(os.path.abspath('.'),
                       os.path.basename(sys.argv[0]).replace('.py', '') + '_' + DATE + '_log.txt')
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

class plot_corr_with_OncoScan:
    def __init__(self, indir, legend, data, outdir, file_name, y_name, x_name,
                 x_title=None, y_title=None, main_title=None, r2_show=True):

        self.indir = indir
        self.legend = legend
        self.data = data
        self.outdir = outdir
        self.filename = file_name
        self.y_name = y_name
        self.x_name = x_name
        self.x_title = x_title
        self.y_title = y_title
        self.main_title = main_title
        self.r2_show = r2_show

    def make_scatterPlot(self):
        """複数の散布図を同時に出力する関数
        """
        logger.info('#--------------------------------------------------------------')
        logger.info('# make multiple scatter plots ...')
        logger.info('#--------------------------------------------------------------')

        n_row = len(self.legend)
        n_col = len(self.data)

        colors = get_colorpalette(n_col)

        if self.x_title is None:
            self.x_title = 'X'
        if self.y_title is None:
            self.y_title = 'Y'
        if self.main_title is None:
            self.main_title = os.path.basename(sys.argv[0]).replace('.py', '') + '_' + DATE

        titles = []
        for i in range(n_row):
            for j in range(n_col):
                name = self.legend[i] + '_' + self.data[j]
                titles.append(name)

        fig = make_subplots(rows=n_row,
                            cols=n_col,
                            subplot_titles=titles,
                            horizontal_spacing=0.07,  # will change
                            vertical_spacing=0.15     # will change
                            )

        for i in range(1, n_col+1):
            for j in range(1, n_row+1):
                ii = i-1
                jj = j-1
                x, y = data_extract(self.indir, self.legend[jj], self.data[ii], self.y_name, self.x_name)
                fig.add_trace(go.Scatter(x=x, y=y, mode='markers', line=dict(color=colors[ii])),
                              row=j,
                              col=i)
                fig.add_trace(go.Scatter(x=x, y=calc_lm(x, y).predict(x.reshape(-1, 1)),
                              line_width=2.0,
                              line=dict(dash='dot', color='gray')),
                              row=j,
                              col=i)
                fig.update_layout(plot_bgcolor='white',
                                  height=1500,
                                  width=1500)
                fig.update_xaxes(title=self.x_title,
                                 showline=True,
                                 linewidth=1,
                                 linecolor='black',
                                 mirror=True,
                                 ticks='inside',
                                 row=j,
                                 col=i)
                fig.update_yaxes(title=self.y_title,
                                 showline=True,
                                 linewidth=1,
                                 linecolor='black',
                                 mirror=True,
                                 ticks='inside',
                                 row=j,
                                 col=i)
                if self.r2_show:
                    r2 = calc_lm(x, y).score(x.reshape(-1, 1), y)
                    text = '*' + 'R^2=' + str(round(r2, 3))
                    fig.add_annotation(x=0.13, y=max(y)+0.05,
                                       text=text,
                                       font=dict(size=8),
                                       showarrow=False,
                                       arrowhead=1,
                                       row=j,
                                       col=i)

        fig.for_each_xaxis(lambda axis: axis.title.update(font=dict(color='black', size=10)))
        fig.for_each_yaxis(lambda axis: axis.title.update(font=dict(color='black', size=10)))

        fig.update_layout(title=dict(text=self.main_title,
                          x=0.5,
                          xanchor='center'),
                          showlegend=False)
        fig.update_annotations(font=dict(size=10))

        pngfile = os.path.join(self.outdir, self.filename + '.png')
        logger.info(f'png file: {pngfile}')
        fig.write_image(pngfile, format='png')
        htmlfile = os.path.join(self.outdir, self.filename + '.html')
        logger.info(f'html file: {htmlfile}')
        fig.write_html(htmlfile)

def get_colorpalette(n_colors):
    """legendとRBG値を対応させる関数

    Returns:
        rgb (list): RBG値が格納された配列
    """
    palette = sns.color_palette('hls', n_colors) #colorpalette = 'hls'
    rgb = ['rgb({},{},{})'.format(*[x*256 for x in rgb]) for rgb in palette]
    return rgb

def calc_lm(x, y):
    """決定係数を求める関数

    Args:
        x (ndarray): xの値
        y (ndarray): yの値

    Returns:
        [type]: yとxの決定係数
    """
    return LinearRegression().fit(x.reshape(-1, 1), y)

def data_extract(indir, legend, data, y_name, x_name):
    """ファイルからyとxの値を抽出する関数

    Args:
        indir (str): DragenとOncoScanのCN数間の相関係数が記載された表ファイルのあるディレクトリ
        legend (str): legend名
        data (str): データ値
        y_name (str): yに対応するヘッダー名
        x_name (str): xに対応するヘッダー名

    Returns:
        y (ndarray): yの値
        x (ndarray): xの値
    """
    logger.info('#--------------------------------------------------------------')
    logger.info('# read input file ...')
    logger.info('#--------------------------------------------------------------')

    file = os.path.join(indir, legend + '_' + data + '_cmpCnv_table.xlsx')
    if len(glob.glob(file)) == 0:
        logger.info(f'input file: {file} not found..')
        sys.exit()

    logger.info(f'input file: {file}')

    logger.info('#--------------------------------------------------------------')
    logger.info('# data processing ...')
    logger.info('#--------------------------------------------------------------')

    df = pd.read_excel(file, engine='openpyxl', sheet_name=0)
    x = df.loc[:, x_name].to_numpy()
    y = df.loc[:, y_name].to_numpy()

    return x, y

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
                        help='input cnv correlation summary path (absolute path) ',
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
                        dest='file_name',
                        help='input png filename (default: ' + os.path.join(os.path.abspath('.'), os.path.basename(sys.argv[0]).replace('.py', '') + '.png') + ')' \
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
    tumor_data = args.X1
    normal_data = args.X2
    outdir = args.out_dir
    file_name = args.file_name

    #===================   praparation   ===================#
    if normal_data is None:
        Data = list(map(str, sorted(tumor_data, key=int)))
    else:
        Data = []
        for x1 in sorted(tumor_data, key=int):
            for x2 in sorted(normal_data, key=int):
                data = str(x1) + '_' + str(x2)
                Data.append(data)

    #===================   run process   ===================#
    y_name = 'OncoScan_CN'
    x_name = 'Dragen_CN'
    y_title = 'OncoScan_CN'
    x_title = 'Dragen_CN'
    main_title = 'corr_with_OncoScan(float)'
    filename = file_name + '_cnvCorr_float'
    pcwO = plot_corr_with_OncoScan(indir, legend, Data, outdir, filename, y_name, x_name,
                                   x_title, y_title, main_title, r2_show=True)
    pcwO.make_scatterPlot()

    y_name = 'OncoScan_CN(int)'
    x_name = 'Dragen_CN'
    y_title = 'OncoScan_CN(int)'
    x_title = 'Dragen_CN'
    main_title = 'corr_with_OncoScan(int)'
    filename = file_name + '_cnvCorr_int'
    pcwO = plot_corr_with_OncoScan(indir, legend, Data, outdir, filename, y_name, x_name,
                                   x_title, y_title, main_title, r2_show=True)
    pcwO.make_scatterPlot()

    logger.info('Done!')
    logger.info('')

    #===================   move file     ===================#
    filename = file_name + '_' + DATE + '_' + os.path.basename(sys.argv[0]).replace('.py', '') + '_log.txt'
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