import io
import os
import re
import sys
import glob
import gzip
import openpyxl
import datetime
import argparse
import pandas as pd
import logging
from logging import Logger, getLogger, Formatter, StreamHandler, FileHandler

DATE = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
LOGFILE = os.path.join(os.path.abspath("."),
                       os.path.basename(sys.argv[0]).replace(".py", "") + "_" + DATE + "_log.txt")
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

class data2venn:
    def __init__(self, indir, legend, data, outdir, png_name, y):
        self.indir = indir
        self.legend = legend
        self.data = data
        self.outdir = outdir
        self.png_name = png_name
        self.y = y

    def extract_data(self):
        """データを抽出する関数

        Returns:
            txt_dir (str): ベン図作成用いるファイルの格納ディレクトリ
        """
        logger.info('#--------------------------------------------------------------')
        logger.info('# data processing ...')
        logger.info('#--------------------------------------------------------------')

        txt_dir = os.path.join(self.indir, 'venn')
        os.makedirs(txt_dir, exist_ok=True)

        if self.y == '1':
            vcf_format(txt_dir, self.indir, self.legend, self.data)
        elif self.y == '2':
            xlsx_format(txt_dir, self.indir, self.legend, self.data)

        return txt_dir

    def run_Rscript(self):
        """vennパッケージを実行する関数
        """
        txt_dir = self.extract_data()
        Rscript = make_Rscript(self.png_name, self.legend, self.data, txt_dir, self.outdir)
        R = '/NGSWORK/NGS/DOCUMENT/SOP/GA/SOP-42_Validation/program/wes_validation/environment/miniconda/envs/py38/bin/Rscript'
        logger.info('#--------------------------------------------------------------')
        logger.info('# run Rscript (venn) ...')
        logger.info('#--------------------------------------------------------------')
        os.system(f'{R} {Rscript}')
        logger.info(f'{R} {Rscript}')

def input_xlsx(infile):
    """xlsxデータ読み込み関数

    Returns:
        path (str): 拡張子を除くファイル名
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

def xlsx_format(txt_dir, indir, legends, data):
    """xlsxファイルからのデータを抽出する関数

    Args:
        txt_dir (str): ベン図作成用いるファイルの出力ディレクトリ
        indir (str): cnv variants summary xlsxが格納されたディレクトリ
        legends (list): legend名が格納された配列
        data (list): データ値が格納された配列
    """
    logger.info('#--------------------------------------------------------------')
    logger.info('# xlsx format ...')
    logger.info('#--------------------------------------------------------------')
    data_list = ',' .join(data)
    for legend in legends:
        files = os.path.join(indir, legend + '*.xlsx')

        for file in sorted(glob.glob(files)):
            path, df = input_xlsx(file)
            name = path.split('/')

            data_value = re.compile(r'[A-Z]\d+_(.*)_[A-Z]{2}_(.*)_cnv_(.*)').search(name[-1]).group(2)

            if data_value in data_list:
                logger.info(f'xlsx file: {file} found.')
                rename = re.compile(r'(.*)_cnv_(.*)').search(name[-1]).group(1)
                txt_file = os.path.join(txt_dir, rename + '.txt')
                with open(txt_file, 'w') as wf:
                    for i in range(2, len(df)):
                        Chr = str(df.loc[i, 'Chr'])
                        if (Chr == 'X') | (Chr == 'Y'):
                            pass
                        else:
                            if df.loc[i, 'View in variant(full) or gene(split)']=='split':
                                if df.loc[i, 'FILTER'] == 'PASS':
                                    if (int(df.loc[i, 'CN']) > 2 ) | (int(df.loc[i, 'CN']) < 2):
                                        info = df.loc[i, 'SV Type'] + ':' + df.loc[i, 'Gene name']
                                        wf.write(f'{info}\n')

def vcffile_type(infile):
    """vcf fileの種類を見分ける関数

    Returns:
        path (str): 拡張子を除くファイル名
        fh (TextIO | TextIOWrapper): ファイルヘッダー
    """
    logger.info('#--------------------------------------------------------------')
    logger.info('# read input file ...')
    logger.info('#--------------------------------------------------------------')

    if os.path.exists(infile) is False:
        logger.error(f'vcf file: {infile} not found..')
        sys.exit()

    logger.info(f'vcf file: {infile} found.')
    path, ext = os.path.splitext(infile)

    if ext == '.gz':
        fh = gzip.open(infile, 'rt')
        name = path.split('/')
        path = re.compile(r'(.*)\.(.*)').search(name[-1]).group(1)
    elif ext == '.vcf':
        fh = open(infile, 'r')
    else:
        logger.info(f'pls input vcf file! --> your input file {infile}')
        sys.exit()

    return path, fh

def vcf_format(txt_dir, indir, legends, data):
    """vcfファイルからのデータを抽出する関数

    Args:
        txt_dir (str): ベン図作成用いるファイルの出力ディレクトリ
        indir (str): cnv variants summary xlsxが格納されたディレクトリ
        legends (list): legend名が格納された配列
        data (list): データ値が格納された配列
    """
    logger.info('#--------------------------------------------------------------')
    logger.info('# vcf format ...')
    logger.info('#--------------------------------------------------------------')
    data_list = ',' .join(data)
    for legend in legends:
        files = os.path.join(indir, legend + '*.vcf.gz')

        for file in sorted(glob.glob(files)):
            path, fh = vcffile_type(file)
            name = path.split('/')
            data_value = re.compile(r'[A-Z]\d+_(.*)_[A-Z]{2}_(.*)').search(name[-1]).group(2)

            if data_value in data_list:
                txt_file = os.path.join(txt_dir, name[-1] + '.txt')
                with open(txt_file, 'w') as wf:
                    with fh as vcf:
                        for line in vcf.readlines():
                            line = line.rstrip('\n')
                            tmp = line.split()
                            if '#' in tmp[0]:
                                pass
                            else:
                                if tmp[6] == 'PASS':
                                    info = tmp[0] + ':' + tmp[1] + ':' + tmp[3] + ':' + tmp[4]
                                    wf.write(f'{info}\n')

def make_Rscript(out_f, legends, data, indir, outdir):
    """vennパッケージRscriptの作成

    Args:
        out_f ([type]): ベン図出力ファイル名
        legends (list): legend名が格納された配列
        data (list): データ値が格納された配列
        indir (str): ベン図作成用いるファイルの格納ディレクトリ
        outdir (str): vennパッケージRscriptの出力ディレクトリ

    Returns:
        Rscript (str): vennパッケージRscript名
    """
    logger.info('#--------------------------------------------------------------')
    logger.info('# make Rscript (venn) ...')
    logger.info('#--------------------------------------------------------------')

    Rscript = os.path.join(outdir, os.path.basename(sys.argv[0]).replace(".py", "") + '.R')
    pngfile = os.path.join(outdir, out_f + '.png')

    with open(Rscript, 'w') as wf:
        wf.write('library(venn)\n')
        wf.write('\n')

        names = []
        for legend in legends:
            for dt in data:
                file = os.path.join(indir, legend + '_' + dt + '.txt')
                if os.path.exists(file) is False:
                    logger.info(f'input file: {file} not found..')
                    sys.exit()
                wf.write(f'{legend} <- scan(\"{file}\", what=character(), sep=\"\\n\", blank.lines.skip=F)\n')
                ele = legend + '=' + legend
                names.append(ele)

        content = ','.join(names)
        wf.write('\n')
        wf.write(f'data <- list({content})\n')
        wf.write('\n')
        wf.write(f'png(\"{pngfile}\")\n')
        wf.write(f'venn(data, ilab=TRUE, zcolor=\"style\", ilcs=1.5, sncs=1.5, box=FALSE)\n')
        wf.write(f'dev.off()\n')

        return Rscript

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
                        help='input vcf path or cnv variants summary file path (absolute path) ',
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
                        '(1: SNV+INDEL, 2: CNV)' ,
                        choices = ['1', '2'],
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
                        dest='fig_name',
                        help='input xlsx filename(default: ' + os.path.join(os.path.abspath('.'), os.path.basename(sys.argv[0]).replace('.py', '') + '.png') + ')' \
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
    legend = sorted(args.legend)
    y = args.y
    tumor_data = args.X1
    normal_data = args.X2
    outdir = args.out_dir
    fig_name = args.fig_name

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
    data2venn(indir, legend, Data, outdir, fig_name, y).run_Rscript()

    logger.info('Done!')
    logger.info('')

    #===================   move file     ===================#
    filename = fig_name + '_' + DATE + '_' + os.path.basename(sys.argv[0]).replace('.py', '') + '_log.txt' 
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