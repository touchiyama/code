import os
import re
import sys
import platform
import getpass
import argparse
import openpyxl
import datetime
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from natsort import natsorted
from logging import Logger, getLogger, Formatter, StreamHandler, FileHandler
import warnings

warnings.filterwarnings('ignore')

DATE = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
LOGFILE = os.path.join(
    os.path.abspath('.'),
    os.path.basename(sys.argv[0]).replace('.py', '') + '_' + DATE + '.log'
)
CD_HIT_PATH = '/NGSWORK/TEST/TEST_UCHIYAMA/tools/cdhit/cd-hit'

def translation(nucl_seq):
    codon = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G' 
    }
    aa = []
    for i in range(0, int(len(nucl_seq)/3)):
        triple = nucl_seq[3*i:3*i+3]
        if codon.get(triple):
            cdn = codon[triple]
        else:
            cdn = '-'
        aa.append(cdn)
    aa_seq = ''.join(aa)

    return aa_seq

def make_table_for_AA(lotID, outdir, prefix, xlsx_file, bg_sample, target_sample, scaling):
    """_summary_

    Args:
        xlsx_file (_type_): _description_
        input_sample (_type_): _description_
    """
    logger.info('')
    logger.info('Start making the frequency table for AA..')

    df_count = pd.read_excel(xlsx_file)

    # cheak sample name ---
    flag = 0
    bg_sample = bg_sample + '_Freq'
    target_sample = [f'{i}_Freq' for i in target_sample]
    legends = target_sample + [bg_sample]

    for legend in legends:
        if not legend in df_count.columns:
            flag = 1
            logger.error(f'No found sample name: {legend}')

    if flag == 1:
        logger.error('Please make sure the input sample names')
        sys.exit()

    # sum counts all over the samples ---
    target_column = [col for col in df_count.columns if 'Count' in col]
    df_seq_samp_sum = pd.concat(
        [df_count['Sequence'], df_count[target_column].sum(axis=1)],
        axis=1
    )
    seq_samp_sum = {}
    for i in range(len(df_seq_samp_sum)):
        sequence = df_seq_samp_sum.loc[i, 'Sequence']
        seq_samp_sum[sequence] = df_seq_samp_sum.loc[i, 0]

    # groupby the same AA seq ---
    df_AA_seq_join = df_count[['AA_seq', 'Seq_ID', 'Sequence']].groupby('AA_seq', sort=False).agg(','.join).reset_index()
    df_AA_seq_cnt = df_count[['AA_seq', 'Seq_ID']].groupby('AA_seq', sort=False).count().reset_index()
    df_AA_seq_cnt.columns = ['AA_seq', 'Seq_count']
    df_AA_seq_merge = pd.merge(df_AA_seq_join, df_AA_seq_cnt, how='left', on='AA_seq')

    # add clusterID to dataframe ---
    df_cluster_ID = pd.DataFrame(df_AA_seq_merge.index, columns=['cluster_ID'])
    df_AA_seq_merge = pd.concat([df_cluster_ID, df_AA_seq_merge], axis=1)

    # select representative nucl sequence ---
    for i in range(df_AA_seq_merge.index.size):
        sequences = df_AA_seq_merge.loc[i, 'Sequence'].split(',')
        sum_cnts = [seq_samp_sum[seq] for seq in sequences]
        max_index = sum_cnts.index(max(sum_cnts))
        df_AA_seq_merge.loc[i, 'Representative_Seq'] = sequences[max_index]

    # sum index ---
    df_AA_seq_sum = df_count.groupby('AA_seq', sort=False).sum().reset_index()
    df_AA_seq = pd.merge(df_AA_seq_merge, df_AA_seq_sum, how='left', on='AA_seq')

    # fold_change ---
    df_AA_fc = df_AA_seq[['cluster_ID', 'AA_seq', 'Seq_ID', 'Sequence', 'Seq_count', 'Representative_Seq']]
    df_target = df_AA_seq[target_sample]
    df_bg = df_AA_seq[bg_sample]

    # prevent ZeroDivisionError to add one to both numerator and denominator ---
    for sample in target_sample:
        name = re.search(r'(.*)_Freq', sample).group(1)
        FC_name = name + '_FC'
        if scaling == 'log10':
            df_AA_fc[FC_name] = np.log10((df_target[sample] + 1) / (df_bg + 1))
        elif scaling == 'log2':
            df_AA_fc[FC_name] = np.log2((df_target[sample] + 1) / (df_bg + 1))
        else:
            df_AA_fc[FC_name] = (df_target[sample] + 1) / (df_bg + 1)

    if prefix:
        fname = lotID + '_' + prefix + '_Consensus_frequency.xlsx'
    else:
        fname = lotID + '_Consensus_frequency.xlsx'

    outfile = os.path.join(outdir, fname)
    with pd.ExcelWriter(outfile) as writer:
        df_AA_seq.to_excel(writer, sheet_name='Cnt_Freq', index=False)
        df_AA_fc.to_excel(writer, sheet_name='FoldChange', index=False)

    logger.info('Finish making the frequency table for AA..')

def make_table_for_BC(lotID, outdir, prefix, xlsx_file,
                              qPCR_file, bg_sample, target_sample):
    """_summary_

    Args:
        lotID (_type_): _description_
        outdir (_type_): _description_
        prefix (_type_): _description_
        xlsx_file (_type_): _description_
        bg_sample (_type_): _description_
        target_sample (_type_): _description_
    """

    logger.info('')
    logger.info('Start making the frequency table for Ascan..')

    df_count = pd.read_excel(xlsx_file)

    # cheak sample name ---
    flag = 0
    bg_sample = bg_sample + '_Freq'
    target_sample = [f'{i}_Freq' for i in target_sample]
    legends = target_sample + [bg_sample]

    for legend in legends:
        if not legend in df_count.columns:
            flag = 1
            logger.error(f'No found sample name: {legend}')

    if flag == 1:
        logger.error('Please make sure the input sample names')
        sys.exit()

    # qPCR ---
    df_qPCR = pd.read_excel(qPCR_file)

    qPCR = {}
    for i in range(len(df_qPCR)):
        Ti = df_qPCR.loc[i, 'Tissue']
        ID = df_qPCR.loc[i, 'Animal No.']
        # DNA ---
        name = 'DNA#' + Ti + '#' + str(ID)
        qPCR[name] = df_qPCR.loc[i, 'Mean']
        # cDNA ---
        name = 'cDNA#' + Ti + '#' + str(ID)
        qPCR[name] = df_qPCR.loc[i, 'Mean']

    # fold_change and Normalization ---
    df_norm = df_count[['BC_ID', 'Sequence']]
    df_target = df_count[target_sample]
    df_bg = df_count[bg_sample]

    for sample in target_sample:
        name = re.search(r'(.*)_Freq', sample).group(1)
        FC_name = name + '_FC'
        Norm_name = name + '_Norm'
        qPCR_value = qPCR[sample]
        FC_value = (df_target[sample] + 1) / (df_bg + 1)
        Norm_value = FC_value * qPCR_value
        df_norm[FC_name] = FC_value
        df_norm[Norm_name] = Norm_value

    if prefix:
        fname = lotID + '_' + prefix + '_Barcode_Summary.xlsx'
    else:
        fname = lotID + '_Barcode_Summary.xlsx'

    outfile = os.path.join(outdir, fname)
    with pd.ExcelWriter(outfile) as writer:
        df_count.to_excel(writer, sheet_name='Cnt_Freq', index=False)
        df_norm.to_excel(writer, sheet_name='FC_Norm', index=False)

    logger.info('Finish making the frequency table for AA..')

def make_table_AA_clustering(lotID, outdir, prefix, xlsx_file, bg_sample,
                             target_sample, scaling, identity, length, wordsize):

    logger.info('')
    logger.info('Start making the frequency table for AA by clustering..')

    df_count = pd.read_excel(xlsx_file)

    # cheak sample name ---
    flag = 0
    bg_sample = bg_sample + '_Freq'
    target_sample = [f'{i}_Freq' for i in target_sample]
    legends = target_sample + [bg_sample]

    for legend in legends:
        if not legend in df_count.columns:
            flag = 1
            logger.error(f'No found sample name: {legend}')

    if flag == 1:
        logger.error('Please make sure the input sample names')
        sys.exit()

    # sum counts all over the samples ---
    target_column = [col for col in df_count.columns if 'Count' in col]
    df_seq_samp_sum = pd.concat(
        [df_count['Sequence'], df_count[target_column].sum(axis=1)],
        axis=1
    )
    seq_samp_sum = {}
    for i in range(len(df_seq_samp_sum)):
        sequence = df_seq_samp_sum.loc[i, 'Sequence']
        seq_samp_sum[sequence] = df_seq_samp_sum.loc[i, 0]

    # clustering ---
    out_dir = os.path.join(outdir, 'cluster')
    os.makedirs(out_dir, exist_ok=True)
    outfa = os.path.join(out_dir, 'AA_seq.fa')
    with open(outfa, 'w') as fa:
        for i in range(len(df_count)):
            header = df_count.loc[i, 'Seq_ID']
            AA_seq = df_count.loc[i, 'AA_seq']
            fa.write(f'>{header}\n')
            fa.write(f'{AA_seq}\n')

    clus_fa = os.path.join(out_dir, 'AA_seq.nr.fa')
    logger.info(f'cd-hit -i {outfa} -o {clus_fa} -n {wordsize} -l {length-1} -c {identity} -bak 1')
    cmd = [
        f'{CD_HIT_PATH} -i {outfa} -o {clus_fa} -n {str(wordsize)} -l {str(length-1)} -c {str(identity)} -bak 1'
    ]
    cdhit_log = os.path.join(out_dir, 'cd-hit.log')
    try:
        subprocess.run(
            cmd, check=True, shell=True, stdout=open(cdhit_log, 'a'),
            stderr=subprocess.PIPE, universal_newlines=True
        )
    except subprocess.CalledProcessError as e:
        logger.error(f'Failed: {e}')
        logger.error('Exit..')
        sys.exit()

    # parse cd-hit.bak.clstr file ---
    clus_bak = os.path.join(out_dir, 'AA_seq.nr.fa.bak.clstr')
    df_clstr = pd.read_csv(clus_bak, sep='\t', header=None).sort_values(0, ascending=True).reset_index(drop=True)
    df_clstr_seq = df_clstr[1].str.split(' ', expand=True)[1]
    df_seqID = df_clstr_seq.str.extract('>(.*)...', expand=True)
    df_clstr_mod = pd.concat([df_clstr[0], df_seqID], axis=1)
    df_clstr_mod.columns = ['cluster_ID', 'Seq_ID']
    df_merge = pd.merge(df_clstr_mod, df_count,  how='left', on='Seq_ID')

    # groupby the same cluster ID ---
    df_AA_seq_join = df_merge[['cluster_ID', 'AA_seq', 'Seq_ID', 'Sequence']].groupby('cluster_ID', sort=False).agg(','.join).reset_index()
    df_AA_seq_cnt = df_merge[['cluster_ID', 'AA_seq']].groupby('cluster_ID', sort=False).count().reset_index()
    df_AA_seq_cnt.columns = ['cluster_ID', 'AA_seq_count']
    df_AA_seq_merge = pd.merge(df_AA_seq_join, df_AA_seq_cnt, how='left', on='cluster_ID')

    # select representative nucl sequence ---
    for i in range(df_AA_seq_merge.index.size):
        sequences = df_AA_seq_merge.loc[i, 'Sequence'].split(',')
        sum_cnts = [seq_samp_sum[seq] for seq in sequences]
        max_index = sum_cnts.index(max(sum_cnts))
        df_AA_seq_merge.loc[i, 'Representative_Seq'] = sequences[max_index]
        df_AA_seq_merge.loc[i, 'Representative_AA_Seq'] = translation(sequences[max_index])

    # sum index ---
    df_AA_seq_sum = df_merge.groupby('cluster_ID', sort=False).sum().reset_index()
    df_AA_seq = pd.merge(df_AA_seq_merge, df_AA_seq_sum, how='left', on='cluster_ID')

    # fold_change ---
    df_AA_fc = df_AA_seq[
        ['cluster_ID', 'AA_seq', 'Seq_ID', 'Sequence',
         'AA_seq_count', 'Representative_Seq', 'Representative_AA_Seq']
    ]
    df_target = df_AA_seq[target_sample]
    df_bg = df_AA_seq[bg_sample]

    # prevent ZeroDivisionError to add one to both numerator and denominator ---
    for sample in target_sample:
        name = re.search(r'(.*)_Freq', sample).group(1)
        FC_name = name + '_FC'
        if scaling == 'log10':
            df_AA_fc[FC_name] = np.log10((df_target[sample] + 1) / (df_bg + 1))
        elif scaling == 'log2':
            df_AA_fc[FC_name] = np.log2((df_target[sample] + 1) / (df_bg + 1))
        else:
            df_AA_fc[FC_name] = (df_target[sample] + 1) / (df_bg + 1)

    if prefix:
        fname = lotID + '_' + prefix + '_Consensus_frequency_clustering.xlsx'
    else:
        fname = lotID + '_Consensus_frequency_clustering.xlsx'

    outfile = os.path.join(outdir, fname)
    with pd.ExcelWriter(outfile) as writer:
        df_AA_seq.to_excel(writer, sheet_name='Cnt_Freq', index=False)
        df_AA_fc.to_excel(writer, sheet_name='FoldChange', index=False)

    logger.info('Finish making the frequency table for AA by clustering..')

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
        help='input count and frequecy table (XLSX format)',
        type=str,
        required=True
    )
    parser.add_argument(
        '-b',
        dest='barcode',
        help='Turn on sequences with barcode analysis ',
        action='store_true'
    )
    parser.add_argument(
        '-q',
        dest='qPCR_file',
        help='input qPCR data each sample for normalization (XLSX format)',
        type=str,
        required=True
    )
    parser.add_argument(
        '-g',
        dest='bg_sample',
        help='input a background sample',
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
        '-s',
        dest='scaling',
        help='choose scaling method',
        choices = ['log10', 'log2', None],
        type=str,
        default=None
    )
    parser.add_argument(
        '-c',
        dest='cluster',
        help='Turn on amino acid sequence clustering analysis using CD-HIT',
        action='store_true'
    )
    parser.add_argument(
        '-i',
        dest='identity',
        help='sequence identity threshold',
        type=float,
        default=0.8
    )
    parser.add_argument(
        '-n',
        dest='length',
        help='length of throw_away_sequences',
        type=int,
        default=10
    )
    parser.add_argument(
        '-w',
        dest='word_size',
        help='word_length, see user\'s guide of CD-HIT for choosing it',
        type=int,
        default=5
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
        make_table_for_BC(
            lotID, outdir, prefix, xlsx_file,
            qPCR_file, bg_sample, target_sample
        )
    elif cluster:
        make_table_for_AA(
            lotID, outdir, prefix, xlsx_file,
            bg_sample, target_sample_sort, scaling
        )
        make_table_AA_clustering(
            lotID, outdir, prefix, xlsx_file, bg_sample,
            target_sample, scaling, identity, length, wordsize
        )
    else:
        make_table_for_AA(
            lotID, outdir, prefix, xlsx_file,
            bg_sample, target_sample_sort, scaling
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