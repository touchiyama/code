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
    """_summary_

    Args:
        lotID (_type_): _description_
        outdir (_type_): _description_
        prefix (_type_): _description_
        Ascan_file (_type_): _description_
        Norm_file (_type_): _description_
        samples (_type_): _description_

    Returns:
        _type_: _description_
    """
    logger.info('')
    logger.info('Start Ascan_with_BC..')

    # norm table ---
    df_fc_norm = pd.read_excel(Norm_file, sheet_name='FC_Norm')

    norm_column = []
    norm_group = {}
    for col in df_fc_norm.columns:
        for sample in samples:
            if (sample in col) and ('Norm' in col):
                norm_column.append(col)
                if norm_group.get(sample):
                    norm_group[sample] += '//' + col
                else:
                    norm_group[sample] = col
    flag = 1
    for col_name in norm_column:
        if not col_name in df_fc_norm.columns:
            flag = 1
            logger.error(f'No found sample name: {col_name}')

    if flag == 1:
        logger.error('Please make sure the input sample names')
        sys.exit()

    df_norm = df_fc_norm[['BC_ID'] + norm_column]

    # ascan table ---
    df_ascan = pd.read_excel(Ascan_file)
    df_tmp = df_ascan.groupby(['pRC', 'AA'])['pAAV'].agg(','.join).reset_index()
    df_tmp.columns = ['Ascan', 'AA_seq', 'BC_ID']
    df_info = df_ascan.loc[:, ['pRC', 'pAAV']].rename(columns={'pRC': 'Ascan', 'pAAV': 'BC_ID'})

    # merge norm and ascan table ---
    df_merge = pd.merge(df_info, df_norm, how='left', on='BC_ID')
    df_mean = df_merge.groupby(['Ascan']).agg('mean').reset_index()
    df_ascan_merge = pd.merge(df_tmp, df_mean, how='left', on='Ascan')

    # summarize stats ---
    for i, name in enumerate(norm_group.keys()):
        groups = norm_group[name].split('//')
        col_vars = [group for group in groups]
        df_melt = pd.melt(
            df_ascan_merge[col_vars],
            value_vars=col_vars
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

def Tab_norm(df_summary, samples):
    """_summary_

    Args:
        df_summary (_type_): _description_
        samples (_type_): _description_

    Returns:
        _type_: _description_
    """
    logger.info('')
    logger.info('Start Tab normalization..')

    # replicate samples grouping ---
    norm_group = {}
    for col in df_summary.columns:
        for sample in samples:
            if (sample in col) and ('Norm' in col):
                if norm_group.get(sample):
                    norm_group[sample] += '//' + col
                else:
                    norm_group[sample] = col

    tmp = []
    for col in norm_group.keys():
        cols = [i for i in norm_group[col].split('//')]
        tmp.append(cols)
    df_group = pd.DataFrame(tmp)

    # prepare for making table ---
    tab_ave = {}
    tab_std = {}
    Tab_ave = [f'{sample}_Tab_Mean' for sample in samples]
    Tab_std = [f'{sample}_Tab_Se' for sample in samples]
    for sample, ave, std in zip(samples, Tab_ave, Tab_std):
        tab_ave[sample] = ave
        tab_std[sample] = std

    # make Tab table ---
    for i in range(df_summary.index.size):
        df_rep = pd.DataFrame()
        df_merge = pd.DataFrame()
        for j in range(df_group.columns.size):
            col_list = df_group[j].to_list()
            df_val = df_summary.loc[i, col_list] / df_summary.loc[i, col_list].sum()
            index = {}
            for col_name, com_name in zip(col_list, samples):
                index[col_name] = com_name
            df_prep = df_val.rename(index=index)
            if j == 0:
                df_rep = df_val
                df_merge = df_prep
            else:
                df_rep = pd.concat([df_rep, df_val])
                df_merge = pd.concat([df_merge, df_prep], axis=1)
        # summarize stats ---
        ave = df_merge.mean(axis=1).rename(tab_ave)
        std = df_merge.std(ddof=1, axis=1).rename(tab_std)
        df_reps = pd.DataFrame(df_rep).T
        df_stats = pd.DataFrame(pd.concat([ave, std]), columns=[i]).T
        df_info = pd.DataFrame(df_summary.loc[i, ['Ascan', 'AA_seq']]).T
        df_line = pd.concat([df_info, df_reps, df_stats], axis=1)
        if i == 0:
            df_Tab = df_line
        else:
            df_Tab = pd.concat([df_Tab, df_line])

    logger.info('Finish Tab normalization..')

    return df_Tab

def Tab_plot(lotID, df_Tab, outdir, prefix):
    """_summary_

    Args:
        lotID (_type_): _description_
        df_Tab (_type_): _description_
        outdir (_type_): _description_
        prefix (_type_): _description_
    """
    logger.info('')
    logger.info('Start Tab_mean visualization..')

    Tab_ave = [col for col in df_Tab.columns if '_Tab_Mean' in col]
    Tab_std = [col for col in df_Tab.columns if '_Tab_Se' in col]

    for i in range(df_Tab.index.size):
        seqID = df_Tab.loc[i, 'Ascan']
        AA_seq = df_Tab.loc[i, 'AA_seq']
        x = [name.replace('_Tab_Mean', '').replace('cDNA#', '') for name in Tab_ave]
        y = df_Tab.loc[i, Tab_ave].tolist()
        z = df_Tab.loc[i, Tab_std].tolist()

        fig, ax = plt.subplots() # figsize=(8, 3)
        error_bar_set = dict(lw=1, capthick=0.5, capsize=2)
        ax.bar(x, y, color='gray', edgecolor='gray', yerr=z , error_kw=error_bar_set, align='center')
        #ax.axvspan(x[0],  x[10], color='yellow', alpha=0.1)
        plt.ylim(0, 1.0)
        plt.xticks(rotation=90)
        #plt.axhline(y=1, color='red', linestyle='--')
        plt.title(f'{seqID} ({AA_seq})')
        plt.ylabel('Tab_mean')

        if prefix:
            fname = lotID + '_' + prefix + '_' + seqID + '_Tab_mean.png'
        else:
            fname = lotID + '_' + seqID + '_Tab_mean.png'

        out_dir = os.path.join(outdir, 'Tab')
        os.makedirs(out_dir, exist_ok=True)
        outfile = os.path.join(out_dir, fname)
        plt.savefig(outfile, bbox_inches='tight', bbox_inches='tight', dpi=300)

    logger.info('Finish Tab_mean visualization..')

def Vab_norm(df_summary, samples):
    logger.info('')
    logger.info('Start Vab normalization..')

    # replicate samples grouping ---
    norm_group = {}
    for col in df_summary.columns:
        for sample in samples:
            if (sample in col) and ('Norm' in col):
                if norm_group.get(sample):
                    norm_group[sample] += '//' + col
                else:
                    norm_group[sample] = col

    # make Vab table ---
    for i, col in enumerate(norm_group.keys()):
        cols = [i for i in norm_group[col].split('//')]
        df_rep = pd.DataFrame()
        for j, col_name in enumerate(cols):
            df_val = df_summary.loc[:, col_name] / df_summary.loc[:, col_name].sum()
            Vab_name = col_name.replace('_Norm', '_Vab')
            df_val = df_val.rename(Vab_name)
            if j == 0:
                df_rep = df_val
            else:
                df_rep = pd.concat([df_rep, df_val], axis=1)
        # summarize stats ---
        df_ave = pd.DataFrame(df_rep.mean(axis=1), columns=[f'{col}_Vab_Mean'])
        df_std = pd.DataFrame(df_rep.std(ddof=1, axis=1), columns=[f'{col}_Vab_Se'])
        df_info = df_summary.loc[:, ['Ascan', 'AA_seq']]
        df_line = pd.concat([df_info, df_rep, df_ave, df_std], axis=1)
        if i == 0:
            df_Vab = df_line
        else:
            df_Vab = pd.merge(df_Vab, df_line, how='left', on=['Ascan', 'AA_seq'])

    logger.info('Finish Vab normalization..')

    return df_Vab

def Vab_plot(lotID, df_Vab, outdir, prefix):
    """_summary_

    Args:
        lotID (_type_): _description_
        df_Tab (_type_): _description_
        outdir (_type_): _description_
        prefix (_type_): _description_
    """
    logger.info('')
    logger.info('Start Vab_mean visualization..')

    Vab_ave = [col for col in df_Vab.columns if '_Vab_Mean' in col]
    Vab_std = [col for col in df_Vab.columns if '_Vab_Se' in col]

    for ave_col, std_col in zip(Vab_ave, Vab_std):
        df_sort = df_Vab.sort_values(ave_col, ascending=False).reset_index()
        x = df_sort['Ascan'].tolist()
        y = df_sort[ave_col].values.tolist()
        z = df_sort[std_col].values.tolist()

        fig, ax = plt.subplots(figsize=(12, 4))
        error_bar_set = dict(lw=1, capthick=0.5, capsize=2)
        bar = ax.bar(x, y, color='white', edgecolor='black', yerr=z , error_kw=error_bar_set)

        idx = x.index('CereAAV')
        bar[idx].set_color('red')

        idx = x.index('mi-342')
        bar[idx].set_color('blue')

        plt.xticks(rotation=90)
        plt.title(f'{ave_col}')
        plt.ylabel('Vab_mean')

        tissue = re.search(f'.*#(.*)_Vab_Mean', ave_col).group(1)
        if prefix:
            fname = lotID + '_' + prefix + '_' + tissue + '_Vab_mean.png'
        else:
            fname = lotID + '_' + tissue + '_Tab_mean.png'

        out_dir = os.path.join(outdir, 'Vab')
        os.makedirs(out_dir, exist_ok=True)
        outfile = os.path.join(out_dir, fname)
        plt.savefig(outfile, bbox_inches='tight', dpi=300)

        logger.info('Finish Vab_mean visualization..')

def output_Vab_Tab(lotID, df_Vab, df_Tab, outdir, prefix):
    """_summary_

    Args:
        lotID (_type_): _description_
        df_Vab (_type_): _description_
        df_Tab (_type_): _description_
        outdir (_type_): _description_
        prefix (_type_): _description_
    """
    logger.info('')
    logger.info('Start output_Vab_Tab..')
    if prefix:
        fname = lotID + '_' + prefix + '_Ascan_Vab_Tab.xlsx'
    else:
        fname = lotID + '_Ascan_Vab_Tab.xlsx'

    outfile = os.path.join(outdir, fname)
    with pd.ExcelWriter(outfile) as writer:
        df_Vab.to_excel(writer, sheet_name='Vab', index=False)
        df_Tab.to_excel(writer, sheet_name='Tab', index=False)

    logger.info('Finish output_Vab_Tab..')

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
        default=os.path.join(os.path.abspath('.'), os.path.basename(sys.argv[0]).replace('.py', ''))
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
    Ascan_file = args.Ascan_file
    Norm_file = args.Norm_file
    samples = args.samples
    prefix = args.prefix
    outdir = args.out_dir

    os.makedirs(outdir, exist_ok=True)

    # run process ---
    df_summary = Ascan_with_BC(lotID, outdir, prefix, Ascan_file, Norm_file, samples)
    df_Tab = Tab_norm(df_summary, samples)
    Tab_plot(lotID, df_Tab, outdir, prefix)
    df_Vab = Vab_norm(df_summary, samples)
    Vab_plot(lotID, df_Vab, outdir, prefix)
    output_Vab_Tab(lotID, df_Vab, df_Tab, outdir, prefix)

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