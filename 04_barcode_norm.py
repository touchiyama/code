import os
import re
import sys 
import platform
import getpass
import argparse
import openpyxl
import datetime
import pandas as pd
import numpy as np
from pathlib import Path
from natsort import natsorted
import matplotlib.pyplot as plt
from logging import Logger, getLogger, Formatter, StreamHandler, FileHandler
import warnings

warnings.filterwarnings('ignore')

DATE = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
LOGFILE = os.path.join(
    os.path.abspath('.'),
    os.path.basename(sys.argv[0]).replace('.py', '') + '_' + DATE + '.log'
)

def Ascan_with_BC(lotID, outdir, prefix, Ascan_file, Norm_file, samples):
    logger.info('')
    logger.info('Start Ascan_with_BC..')

    df_ascan = pd.read_excel(Ascan_file)
    df_tmp = df_ascan.groupby(['pRC', 'AA'])['pAAV'].agg(','.join).reset_index()
    df_tmp.columns = ['Ascan', 'AA_seq', 'BC_ID']

    df_fc_norm = pd.read_excel(Norm_file, sheet_name='FC_Norm')
    norm_column = [col for col in df_fc_norm.columns if 'Norm' in col]
    df_norm = df_fc_norm[['BC_ID'] + norm_column]
    df_info = df_ascan.loc[:, ['pRC', 'pAAV']].rename(columns={'pRC': 'Ascan', 'pAAV': 'BC_ID'})
    df_merge = pd.merge(df_info, df_norm, how='left', on='BC_ID')
    df_mean = df_merge.groupby(['Ascan']).agg('mean').reset_index()
    df_ascan_merge = pd.merge(df_tmp, df_mean, how='left', on='Ascan')

    for i, name in enumerate(samples):
        df_melt = pd.melt(
            df_ascan_merge[[f'{name}#101_Norm', f'{name}#102_Norm']],
            value_vars=['1', '2']
        )
        ave = df_melt.mean()
        std = df_melt.std(ddof=1)
        df = pd.concat([df_ascan_merge.loc[:, 'Ascan'], ave, std], axis=1)
        df.columns = ['Ascan', f'{name}_Mean', f'{name}_Se']
        if i == 0:
            df_new = df
        else:
            df_new = pd.merge(df_new, df, how='left', on='Ascan')

    df_summary = pd.merge(df_ascan_merge, df_new, how='left', on='Ascan')

    if prefix:
        fname = lotID + '_' + prefix + '_Ascan_Summary.xlsx'
    else:
        fname = lotID + '_Ascan_Summary.xlsx'

    outfile = os.path.join(outdir, fname)
    df_summary.to_excel(outfile, index=False)

    return df_summary



def parameters(__desc__):
    """入力されたコマンドライン引数をパースする関数
    Args:
        __desc__ (str): usege文

    Returns:
        args (argparse.Namespace): コマンドライン引数パース結果
    """
    parser = argparse.ArgumentParser(
        description = __desc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter #show default
    )
    parser.add_argument(
        '-l',
        dest='lotID',
        help='input LotID',
        type=str,
        required=True
    )
    parser.add_argument(
        '-x',
        dest='xlsx_file',
        help='input normaliazion table (XLSX format)',
        type=str,
        required=True
    )
    parser.add_argument(
        '-q',
        dest='qPCR_file',
        help='input qPCR data each sample for normalization (XLSX format)',
        type=str,
        required=True
    )
    parser.add_argument(
        '-t',
        dest='target_sample',
        help='input the target samples',
        nargs='*',
        type=str,
        required=True
    )
    parser.add_argument(
        '-p',
        dest='prefix',
        help='input prefix (etc. sample group name..)',
        type=str,
        default=None
    )
    parser.add_argument(
        '-o',
        dest='out_dir',
        help='input output directory path',
        type=str,
        default=os.path.join(os.path.abspath('.'), '02_fold_change_with_aa')
    )
    parser.add_argument(
        '-log',
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
    logger.info('')
    logger.info('#--------------------------------------------------------------')
    logger.info(f"python: {platform.python_version()}")
    logger.info(f'program: {__file__}')
    logger.info(f"user: {getpass.getuser()}")
    logger.info(f"current directory: {Path.cwd()}")
    cla = vars(args)
    for key in sorted(cla.keys()):
        logger.info(f"{key}: {cla[key]}")
    logger.info("#--------------------------------------------------------------")

def main(args):
    logger.info('')
    logger.info('Output log...')

    # setUp parameter ---
    show_param(args)

    lotID = args.lotID
    xlsx_file = args.xlsx_file
    prefix = args.prefix
    bg_sample = args.bg_sample
    target_sample = args.target_sample
    scaling = args.scaling
    barcode = args.barcode
    qPCR_file = args.qPCR_file
    cluster = args.cluster
    identity = args.identity
    length = args.length
    wordsize = args.word_size
    outdir = args.out_dir

    target_sample_sort = sorted(target_sample)

    samplename_num_only = 0
    for sample in target_sample:
        if sample.isnumeric():
            samplename_num_only += 1

    if samplename_num_only == len(target_sample):
        target_sample_sort = natsorted(target_sample)

    os.makedirs(outdir, exist_ok=True)

    # run process ---
    if barcode:
        make_table_for_Ascan_with_BC(
            lotID, outdir, prefix, xlsx_file,
            qPCR_file, bg_sample, target_sample
        )

    logger.info('')
    logger.info('Done!')
    logger.info('')

    # move file ---
    # os.system(f'mv {LOGFILE} {outdir}')

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