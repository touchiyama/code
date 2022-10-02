#!/usr/bin/python
# -*- coding: utf-8 -*-
"""RNA-Seq_checker_vXX.py: cheaker for outputs from DRAGEN RNA-Seq pipeline

 | Author: uchiyamat, Shogo Satoyama
 | Version: 1.0
 | © 2022/09/28 Takara Bio Inc.

 * A supplementary script in DRAGEN RNA-Seq pipeline

Todo:
    * DRAGEN RNA-Seq pipeline

Examples:
        ::
        $python3 RNA-Seq_checker_v01.py -y XXXX.yaml -x XXXX_fastq_file.xlsx
"""
from distutils.log import error
import os
from glob import iglob
from pathlib import Path
from random import sample
import re
from select import select
import sys
import time
import datetime
import subprocess
import argparse
from logging import getLogger, StreamHandler, Formatter, FileHandler, basicConfig
from typing import no_type_check, List, Set, Dict, Tuple, Optional, Union
from pathlib import Path
import yaml
import platform
import getpass
from packaging import version
import glob
import openpyxl
import pandas as pd
import filecmp
import warnings
warnings.simplefilter('ignore')

# SOP(pipeline path) ---
__file__ = '/NGSWORK/NGS/DOCUMENT/SOP/GA/SOP-33_Dragen/RNA/pipeline'
base_dir = Path(__file__)

# log file ---
DATE = datetime.datetime.today().strftime('%Y_%m_%d_%H_%M_%S') # 作業日の取得(YYYY_MM_DD_hh_mm_ss形式)
LOGFILE = os.path.join(
    __file__,
    os.path.basename(sys.argv[0]).replace('.py', '') + '_' + DATE + '.log'
)

# function ---
def get_fastqFile_xlsx(xlsx_file, tkr_format='_R1.fastq.gz'):
    """fastq_file_listに記載された依頼書記載名称、サンプルID、ライブラリーID、fastqファイル名R1とR2を獲得
    Args:
        xlsx_file (str): fastq file list名（xlsx）
        tkr_format (str): takara仕様のfastq file名（illumina出力とは異なることに留意）
    Returns:
        xlsx_data (dir): takaraIDとリンクする依頼書記載名称、ライブラリーID、ライブラリー名
    """
    logger.info("")
    logger.info('Reading FASTQ file list ...')
    logger.info(f'Input: {xlsx_file}')
    if not os.path.exists(xlsx_file):
        logger.error(f'No found :( , {xlsx_file}')
        sys.exit()

    df = pd.read_excel(xlsx_file)
    col_name = "依頼書記載名称"
    col_id = "サンプルID"
    col_lib = "ライブラリーID"
    col_r1 = "fastqファイル名称(Read1)"
    col_r2 = "fastqファイル名称(Read2)"
    xlsx_data = {}
    for i in range(len(df)):
        sample_name = str(df.loc[i, col_name]).replace(' ', '')
        sample_id = str(df.loc[i, col_id])
        lib_id = str(df.loc[i, col_lib])
        r1 = str(df.loc[i, col_r1])
        lib_name = re.compile(r'(.*)'+ tkr_format).search(r1).group(1)
        xlsx_data[sample_id] = sample_name + '//' + lib_id + '//' + lib_name

    return xlsx_data

def read_yaml(yaml_file: str):
    """yamlファイルを読み込み辞書型に変換
    Args:
        yaml_file (str): yamlファイルへのパス
    Returns:
        yamlファイルをパースし、値を格納した辞書型オブジェクト
    """
    if not os.path.exists(yaml_file):
        logger.error(f'No found :( , {yaml_file}')
        sys.exit()

    with open(yaml_file, "r", encoding="utf-8") as F:
        if version.parse(yaml.__version__) >= version.parse("5.1"):
            dat = yaml.load(F, Loader=yaml.FullLoader)
        else:
            dat = yaml.load(F)
    return dat

def refs(base_dir):
    ref_dir = base_dir.joinpath("dat", "reference")
    yamls = iglob(f"{ref_dir.joinpath('*.yaml')}")
    refs = []
    for yaml in yamls:
        refs.append(Path(yaml).stem)
    refs = sorted(refs)

    return refs

def refs_yaml(base_dir):
    ref_hash = {}
    fasta = {}
    annotation = {}
    annotation_txt = {}
    for ref in refs(base_dir):
        ref_dir = base_dir.joinpath("dat", "reference")
        yaml_file = os.path.join(ref_dir, ref + '.yaml')
        config = read_yaml(yaml_file)
        ref_hash[ref] = config['ref_hash']
        fasta[ref] = config['fasta']['file']
        annotation[ref] = config['annotation']['file']
        annotation_txt[ref] = config['annotation_txt']

    return ref_hash, fasta, annotation, annotation_txt

class cheak_expression_report:
    def __init__(self, report_dir, xlsx_file, mapping_dir, lot):
        self.report_file = os.path.join(report_dir, '11-analysis.tex')
        self.xlsx_file = xlsx_file
        self.mapping_file = os.path.join(mapping_dir, f'{lot}_mapping_info.xlsx')

    def get_mapping_info(self):
        logger.info("")
        logger.info('Reading maping_info file ...')
        logger.info(f'Input: {self.mapping_file}')
        if not os.path.exists(self.mapping_file):
            logger.error(f'No found :( , {self.mapping_file}')
            return None
        else:
            df = pd.read_excel(self.mapping_file)
            col_id = 'Sample ID'
            col_total_reads = 'Total input reads'
            col_map_reads = 'Mapped reads'
            col_map_readp = 'Mapped reads (%)'
            col_Covered_bases = 'Covered bases (%)'
            col_Total_aligns = 'Total alignments'
            map_data = {}
            for i in range(len(df)):
                sample_id = str(df.loc[i, col_id])
                total_reads = str(int(df.loc[i, col_total_reads]))
                map_reads = str(int(df.loc[i, col_map_reads]))
                map_readp =  '{:.2f}'.format(df.loc[i, col_map_readp])
                covered_bases = '{:.2f}'.format(df.loc[i, col_Covered_bases])
                total_aligns = str(int(df.loc[i, col_Total_aligns]))
                map_data[sample_id] = total_reads + '//' + map_reads + '//' + map_readp + '//' + covered_bases + '//' + total_aligns

        return map_data

    def result(self):
        logger.info('')
        logger.info('<--- Cheak expression_report --->')
        xlsx_data = get_fastqFile_xlsx(self.xlsx_file, tkr_format='_R1.fastq.gz')
        map_data  = self.get_mapping_info()
        if map_data is None:
            return False
        else:
            NG_flag = False
            NG_list = []
            flag = 0
            with open(self.report_file, 'r', encoding='utf-8') as rf:
                for line in rf.readlines():
                    line = line.rstrip('\n')
                    if 'caption{使用したシーケンスデータ}' in line:
                        flag = 1
                        cnt = 0
                    elif 'caption{マッピング結果}' in line:
                        flag = 1
                        cnt = 0
                    else:
                        if 'endhead' in line:
                            if flag == 1:
                                flag = 2
                        elif 'hline' in line:
                            if flag == 2:
                                flag = 0
                        else:
                            if flag == 2:
                                line = line.replace(' ', '').replace('\\', '').split('&')
                                if len(line) <= 3:
                                    cnt += 1
                                    if cnt == 1:
                                        logger.info('')
                                        logger.info(f'解析対象のサンプルIDとリード配列ファイル(FASTQ形式)の対応表の確認')
                                    if len(line) < 3:
                                        NG_flag = True
                                        NG_list.append('解析対象のサンプルIDとリード配列ファイル(FASTQ形式)の対応表')
                                        logger.error(f'{line[0]} was BAD (table column)')
                                    else:
                                        if xlsx_data.get(line[0]):
                                            names = xlsx_data[line[0]].split('//')
                                            if (names[0] == line[1]) & (names[2] in line[2]):
                                                logger.info(f'{line[0]} -> GOOD :)')
                                            else:
                                                NG_flag = True
                                                NG_list.append('解析対象のサンプルIDとリード配列ファイル(FASTQ形式)の対応表')
                                                if names[0] != line[1]:
                                                    logger.error(f'{line[0]} -> [依頼書記載名称] was BAD (tex: {line[0]}, fastq_file_list: {names[1]})')
                                                else:
                                                    logger.error(f'{line[0]} -> [ファイル名] was BAD (tex: {line[2]}, fastq_file_list: {names[2]})')
                                        else:
                                            NG_flag = True
                                            NG_list.append('解析対象のサンプルIDとリード配列ファイル(FASTQ形式)の対応表')
                                            logger.error(f'No found :( , {line[0]} ')
                                else:
                                    cnt += 1
                                    if cnt == 1:
                                        logger.info('')
                                        logger.info(f'リード配列をゲノム配列にマッピングした結果の統計値の確認')
                                        # サンプルID & Total reads & Mapped reads & Mapped reads (\%) & Covered bases (\%) & Total alignments \\
                                    if len(line) < 6:
                                        NG_flag = True
                                        NG_list.append('解析対象のサンプルIDとリード配列ファイル(FASTQ形式)の対応表')
                                        logger.error(f'{line[0]} was BAD (table column)')
                                    else:
                                        if map_data.get(line[0]):
                                            tex_join = '//'.join(line[1:])
                                            tex_data = tex_join.replace(',', '')
                                            if tex_data == map_data[line[0]]:
                                                logger.info(f'{line[0]} -> GOOD :)')
                                            else:
                                                NG_flag = True
                                                NG_list.append('リード配列をゲノム配列にマッピングした結果の統計値')
                                                title = [
                                                    '[Total reads]', '[Mapped reads]', '[Mapped reads%]',
                                                    '[Covered bases%]', '[Total alignments]'
                                                ]
                                                texs = [i for i in tex_data.split('//')]
                                                maps = [i for i in map_data[line[0]].split('//')]
                                                for i in range(len(texs)):
                                                    if texs[i] != maps[i]:
                                                        logger.error(f'{line[0]} -> {title[i]} was BAD (tex: {texs[i]}, map_info: {maps[i]})')
                                        else:
                                            NG_flag = True
                                            NG_list.append('リード配列をゲノム配列にマッピングした結果の統計値')
                                            logger.error(f'No found :( , {line[0]} ')

            logger.info('>--- RESULT (Cheak expression_report) ---<')
            if NG_flag:
                logger.error(f"DEG_report looked BAD..")
                logger.error("### NG items ###")
                logger.error('Pls cheak the following items: ')
                for i, NG in enumerate(list(set(NG_list)), 1):
                    logger.error(f'({i}) {NG}')
                return False
            else:
                logger.info('Everything was GOOD :)')
                return True

class cheak_expression_appendix:
    def __init__(self, report_dir, xlsx_file, mapping_dir, expression_dir, lot):
        self.report_file = os.path.join(report_dir, '12-apendix.tex')
        self.xlsx_file = xlsx_file
        self.mapping_file = os.path.join(mapping_dir, f'{lot}_mapping_info.xlsx')
        self.expression_dir = expression_dir
        self.lot = lot

    def cheak_mapping_info(self):
        logger.info('')
        logger.info(f'<--- Cheak {self.lot}_map_info_appendix --->')
        logger.info(f'- {self.mapping_file}')
        if not os.path.exists(self.mapping_file):
            logger.error(f'No found :( , {self.mapping_file}')
            return True
        else:
            df = pd.read_excel(self.mapping_file)
            map_info = []
            for col in df.columns.to_list():
                map_info.append(col.replace(' ', ''))

            flag = 0
            tex = []
            with open(self.report_file, 'r', encoding='utf-8') as rf:
                for line in rf.readlines():
                    line = line.rstrip('\n')
                    if 'caption{使用したシーケンスデータ}' in line:
                        table_name = line
                        flag = 1
                        cnt = 0
                    else:
                        if 'endhead' in line:
                            if flag == 1:
                                flag = 2
                        elif 'hline' in line:
                            if flag == 2:
                                flag = 0
                        else:
                            if flag == 2:
                                line = line.replace(' ', '').replace('\\', '').split('&')
                                cnt += 1
                                if cnt == 1:
                                    logger.info('')
                                    logger.info(f'{table_name}')
                                if '#' in line[0]:
                                    pass
                                else:
                                    tex.append(line[0])

        logger.info('>--- RESULT (Cheak map_info_appendix) ---<')
        if map_info == tex:
            logger.info('Everything was GOOD :)')
            return True
        else:
            diff_list = list(set(tex) ^ set(map_info))
            logger.error('map_info_appendix looked BAD..')
            logger.error("### NG items ###")
            logger.error('Pls cheak the following items: ')
            for i, NG in enumerate(list(set(diff_list)), 1):
                if len(tex) > len(map_info):
                    logger.error(f'{i}[{NG}] was BAD (Deleted item from tex)')
                else:
                    logger.error(f'{i}[{NG}] was BAD (Add item to tex)')
            return False

    def cheak_expression(self):
        logger.info("")
        logger.info("<--- Cheak expression_appendix --->")
        logger.info(f'- lotID: {self.lot}')
        logger.info(f'- expression_dir: {self.expression_dir}')

        NG_flag = False
        NG_list = []
        files = [
            os.path.join(self.expression_dir, f'{self.lot}_gene_expression.xlsx'),
            os.path.join(self.expression_dir, f'{self.lot}_transcript_expression.xlsx'),
        ]
        for file in files:
            if os.path.exists(file):
                pass
            else:
                logger.error(f"> RESULT: Not found :( , {os.path.basename(file)}")
                NG_flag = True

        if NG_flag == False:
            exp_col = []
            for file in files:
                df = pd.read_excel(file)
                for col in df.columns.to_list():
                    if 'TPM' in col:
                        exp_col.append('TPM')
                    elif 'count' in col:
                        exp_col.append('count')
                    else:
                        exp_col.append(col.replace(' ', ''))
            exp = list(set(exp_col))

            flag = 0
            tex_col = []
            with open(self.report_file, 'r', encoding='utf-8') as rf:
                for line in rf.readlines():
                    line = line.rstrip('\n')
                    if 'caption{使用したシーケンスデータ}' in line:
                        table_name = line
                        flag = 1
                        cnt = 0
                    elif 'caption{使用したシーケンスデータ}' in line:
                        table_name = line
                        flag = 1
                        cnt = 0
                    else:
                        if 'endhead' in line:
                            if flag == 1:
                                flag = 2
                        elif 'hline' in line:
                            if flag == 2:
                                flag = 0
                        else:
                            if flag == 2:
                                line = line.replace(' ', '').replace('\\', '').split('&')
                                cnt += 1
                                if cnt == 1:
                                    logger.info('')
                                    logger.info(f'{table_name}')
                                if '#' in line[0]:
                                    pass
                                else:
                                    if 'TPM' in col:
                                        tex_col.append('TPM')
                                    elif 'count' in col:
                                        tex_col.append('count')
                                    else:
                                        tex_col.append(line[0])
            tex = list(set(tex_col))

        logger.info('>--- RESULT (Cheak expression_appendix) ---<')
        if exp == tex:
            logger.info('Everything was GOOD :)')
            return True
        else:
            diff_list = list(set(tex) ^ set(exp))
            logger.error('expression_appendix looked BAD..')
            logger.error("### NG items ###")
            logger.error('Pls cheak the following items: ')
            for i, NG in enumerate(list(set(diff_list)), 1):
                if len(tex) > len(exp):
                    logger.error(f'{i}[{NG}] was BAD (Deleted item from tex)')
                else:
                    logger.error(f'{i}[{NG}] was BAD (Add item to tex)')
            return False

class cheak_fusion_report:
    def __init__(self, report_dir, xlsx_file, mapping_dir, fusion_subdir, lot):
        self.report_file = os.path.join(report_dir, '11-analysis.tex')
        self.xlsx_file = xlsx_file
        self.mapping_file = os.path.join(mapping_dir, f'{lot}_mapping_info.xlsx')
        self.fusion_file = os.path.join(fusion_subdir, f'{lot}_fusion_info.txt')

    def get_mapping_info(self):
        logger.info("")
        logger.info('Read maping_info file ...')
        logger.info(f'Input: {self.mapping_file}')
        if not os.path.exists(self.mapping_file):
            logger.error(f'No found :( , {self.mapping_file}')
            return None
        else:
            df = pd.read_excel(self.mapping_file)
            col_id = 'Sample ID'
            col_total_reads = 'Total input reads'
            col_map_reads = 'Mapped reads'
            col_map_readp = 'Mapped reads (%)'
            col_Covered_bases = 'Covered bases (%)'
            col_Total_aligns = 'Total alignments'
            map_data = {}
            for i in range(len(df)):
                sample_id = str(df.loc[i, col_id])
                total_reads = str(int(df.loc[i, col_total_reads]))
                map_reads = str(int(df.loc[i, col_map_reads]))
                map_readp =  '{:.2f}'.format(df.loc[i, col_map_readp])
                covered_bases = '{:.2f}'.format(df.loc[i, col_Covered_bases])
                total_aligns = str(int(df.loc[i, col_Total_aligns]))
                map_data[sample_id] = total_reads + '//' + map_reads + '//' + map_readp + '//' + covered_bases + '//' + total_aligns

            return map_data

    def get_fusion_info(self):
        logger.info('')
        logger.info('Reading fusion_info file ...')
        logger.info(f'Input: {self.fusion_file}')
        if not os.path.exists(self.fusion_file):
            logger.error(f'No found :( , {self.fusion_file}')
            sys.exit()

        df = pd.read_csv(self.fusion_file, sep='\t')
        col_id = 'Sample ID'
        col_total = 'Total'
        col_intra = 'Intra'
        col_inter = 'Inter'
        fusion_data = {}
        for i in range(len(df)):
            sample_id = str(df.loc[i, col_id])
            total = str(df.loc[i, col_total])
            intra = str(df.loc[i, col_intra])
            inter = str(df.loc[i, col_inter])
            fusion_data[sample_id] = total + '//' + intra + '//' + inter

        return fusion_data

    def result(self):
        logger.info('')
        logger.info('<--- Cheak fusion_report --->')
        xlsx_data = get_fastqFile_xlsx(self.xlsx_file, tkr_format='_R1.fastq.gz')
        map_data  = self.get_mapping_info()
        fusion_data = self.get_fusion_info()
        if map_data is None:
            return False
        else:
            NG_flag = False
            NG_list = []
            flag = 0
            with open(self.report_file, 'r', encoding='utf-8') as rf:
                for line in rf.readlines():
                    line = line.rstrip('\n')
                    if 'caption{使用したシーケンスデータ}' in line:
                        flag = 1
                        cnt = 0
                    elif 'caption{マッピング結果}' in line:
                        flag = 1
                        cnt = 0
                    elif 'caption{融合遺伝子検出結果}' in line:
                        flag = 1
                        cnt = 0
                    else:
                        if 'endhead' in line:
                            if flag == 1:
                                flag = 2
                        elif 'hline' in line:
                            if flag == 2:
                                flag = 0
                        else:
                            if flag == 2:
                                line = line.replace(' ', '').replace('\\', '').split('&')
                                if len(line) <= 3:
                                    cnt += 1
                                    if cnt == 1:
                                        logger.info('')
                                        logger.info(f'解析対象のサンプルIDとリード配列ファイル(FASTQ形式)の対応表の確認')
                                    if len(line) < 3:
                                        NG_flag = True
                                        NG_list.append('解析対象のサンプルIDとリード配列ファイル(FASTQ形式)の対応表')
                                        logger.error(f'{line[0]} was BAD (table column)')
                                    else:
                                        if xlsx_data.get(line[0]):
                                            names = xlsx_data[line[0]].split('//')
                                            if (names[0] == line[1]) & (names[2] in line[2]):
                                                logger.info(f'{line[0]} -> GOOD :)')
                                            else:
                                                NG_flag = True
                                                NG_list.append('解析対象のサンプルIDとリード配列ファイル(FASTQ形式)の対応表')
                                                if names[0] != line[1]:
                                                    logger.info(f'{line[0]} -> [依頼書記載名称] was BAD (tex: {line[0]}, fastq_file_list: {names[1]})')
                                                else:
                                                    logger.info(f'{line[0]} -> [ファイル名] was BAD (tex: {line[2]}, fastq_file_list: {names[2]})')
                                        else:
                                            NG_flag = True
                                            NG_list.append('解析対象のサンプルIDとリード配列ファイル(FASTQ形式)の対応表')
                                            logger.info(f'No found :( , {line[0]} ')

                                elif len(line) <= 4:
                                    cnt += 1
                                    if cnt == 1:
                                        logger.info('')
                                        logger.info('各サンプルで検出された融合遺伝子数の確認')
                                    if len(line) < 4:
                                        NG_flag = True
                                        NG_list.append('各サンプルで検出された融合遺伝子数の対応表')
                                        logger.error(f'{line[0]} was BAD (table column)')
                                    else:
                                        if fusion_data.get(line[0]):
                                            tex_join = '//'.join(line[1:])
                                            tex_data = tex_join.replace(',', '')
                                            if tex_data == fusion_data[line[0]]:
                                                logger.info(f'{line[0]} -> GOOD :)')
                                            else:
                                                NG_flag = True
                                                NG_list.append('各サンプルで検出された融合遺伝子数の対応表')
                                                title = ['[Total]', '[Intra]', '[Inter]']
                                                texs = [i for i in tex_data.split('//')]
                                                fusions = [i for i in fusion_data[line[0]].split('//')]
                                                for i in range(len(texs)):
                                                    if texs[i] != fusions[i]:
                                                        logger.error(f'{line[0]} -> {title[i]} was BAD (tex: {texs[i]}, fusion_info: {maps[i]})')
                                        else:
                                            NG_flag = True
                                            NG_list.append('各サンプルで検出された融合遺伝子数の対応表')
                                            logger.error(f'No found :( , {line[0]} ')

                                else:
                                    cnt += 1
                                    if cnt == 1:
                                        logger.info('')
                                        logger.info(f'リード配列をゲノム配列にマッピングした結果の統計値の確認')
                                        # サンプルID & Total reads & Mapped reads & Mapped reads (\%) & Covered bases (\%) & Total alignments \\
                                    if len(line) < 6:
                                        NG_flag = True
                                        NG_list.append('解析対象のサンプルIDとリード配列ファイル(FASTQ形式)の対応表')
                                        logger.error(f'{line[0]} was BAD (table column)')
                                    else:
                                        if map_data.get(line[0]):
                                            tex_join = '//'.join(line[1:])
                                            tex_data = tex_join.replace(',', '')
                                            if tex_data == map_data[line[0]]:
                                                logger.info(f'{line[0]} -> GOOD :)')
                                            else:
                                                NG_flag = True
                                                NG_list.append('リード配列をゲノム配列にマッピングした結果の統計値')
                                                title = [
                                                    '[Total reads]', '[Mapped reads]', '[Mapped reads%]',
                                                    '[Covered bases%]', '[Total alignments]'
                                                ]
                                                texs = [i for i in tex_data.split('//')]
                                                maps = [i for i in map_data[line[0]].split('//')]
                                                for i in range(len(texs)):
                                                    if texs[i] != maps[i]:
                                                        logger.info(f'{line[0]} -> {title[i]} was BAD (tex: {texs[i]}, map_info: {maps[i]})')
                                        else:
                                            NG_flag = True
                                            NG_list.append('リード配列をゲノム配列にマッピングした結果の統計値')
                                            logger.info(f'No found :( , {line[0]} ')

            logger.info('>--- RESULT (Cheak fusion_report) ---<')
            if NG_flag:
                logger.error(f"FUSION_report looked BAD..")
                logger.error("### NG items ###")
                logger.error('Pls cheak the following items: ')
                for i, NG in enumerate(list(set(NG_list)), 1):
                    logger.error(f'({i}) {NG}')
                return False
            else:
                logger.info('Everything was GOOD :)')
                return True

class cheak_fusion_appendix:
    def __init__(self, report_dir, xlsx_file, fusion_dir, lot):
        self.report_file = os.path.join(report_dir, '12-apendix.tex')
        self.xlsx_file = xlsx_file
        self.fusion_file = os.path.join(fusion_dir, f'{lot}_fusion_info.txt')
        self.lot = lot

    def cheak_fusion(self):
        logger.info("")
        logger.info("<--- Cheak fusion_appendix --->")
        logger.info(f'- lotID: {self.lot}')
        logger.info(f'- fusion_dir: {self.expression_dir}')

        NG_flag = False
        NG_list = []
        files = [
            os.path.join(self.expression_dir, f'{self.lot}_gene_expression.xlsx'),
            os.path.join(self.expression_dir, f'{self.lot}_transcript_expression.xlsx'),
        ]
        for file in files:
            if os.path.exists(file):
                pass
            else:
                logger.error(f"> RESULT: Not found :( , {os.path.basename(file)}")
                NG_flag = True

        if NG_flag == False:
            exp_col = []
            for file in files:
                df = pd.read_excel(file)
                for col in df.columns.to_list():
                    if 'TPM' in col:
                        exp_col.append('TPM')
                    elif 'count' in col:
                        exp_col.append('count')
                    else:
                        exp_col.append(col.replace(' ', ''))
            exp = list(set(exp_col))

            flag = 0
            tex_col = []
            with open(self.report_file, 'r', encoding='utf-8') as rf:
                for line in rf.readlines():
                    line = line.rstrip('\n')
                    if 'caption{使用したシーケンスデータ}' in line:
                        table_name = line
                        flag = 1
                        cnt = 0
                    elif 'caption{使用したシーケンスデータ}' in line:
                        table_name = line
                        flag = 1
                        cnt = 0
                    else:
                        if 'endhead' in line:
                            if flag == 1:
                                flag = 2
                        elif 'hline' in line:
                            if flag == 2:
                                flag = 0
                        else:
                            if flag == 2:
                                line = line.replace(' ', '').replace('\\', '').split('&')
                                cnt += 1
                                if cnt == 1:
                                    logger.info('')
                                    logger.info(f'{table_name}')
                                if '#' in line[0]:
                                    pass
                                else:
                                    if 'TPM' in col:
                                        tex_col.append('TPM')
                                    elif 'count' in col:
                                        tex_col.append('count')
                                    else:
                                        tex_col.append(line[0])
            tex = list(set(tex_col))

        logger.info('>--- RESULT (Cheak fusion_appendix) ---<')
        if exp == tex:
            logger.info('Everything was GOOD :)')
            return True
        else:
            diff_list = list(set(tex) ^ set(exp))
            logger.error('map_info_appendix looked BAD..')
            logger.error("### NG items ###")
            logger.error('Pls cheak the following items: ')
            for i, NG in enumerate(list(set(diff_list)), 1):
                if len(tex) > len(exp):
                    logger.error(f'{i}[{NG}] was BAD (Deleted item from tex)')
                else:
                    logger.error(f'{i}[{NG}] was BAD (Add item to tex)')
            return False

class display_cheak:
    def __init__(self, yaml_file, xlsx_file, qc):
        """コンストラクタの定義
        Args:
            yaml_file (str): YAMLファイルパス
        """
        self.config = read_yaml(yaml_file)
        self.xlsx_file = xlsx_file
        self.qc = qc

    def response(self, sign):
        """標準入力に対して応答する関数
        Args:
            sign (str): 標準入力結果
        Returns:
            受け取った標準入力の結果を返す
        """
        print(f'Input: {sign}')
        if sign == 'y':
            print('OK, continue this process..')
        elif sign == 'n':
            print('Pls Retry..')
            print('==============================================================')
            print('\n')
            os.system(f'rm -rf {LOGFILE}')
            sys.exit()
        else:
            print('Pls input \"y\" or \"n\"')
            sign = input('Are all items correct? [y/n] ').strip()
            self.response(sign)

        return sign

    def cheak_args(self):
        """入力された引数を表示させる関数
        """
        lotID = self.config['lot_id']
        expression = self.config['expression']
        fusion = self.config['fusion']
        try:
            lib_type = self.config['library_type']
        except KeyError:
            lib_type = None
        bam_nouhin = self.config['bam_nouhin']
        seq_data = os.path.dirname(self.config['samples'][0]['libs'][0]['fastqs'][0][0])
        samples = len([ sample['id'] for sample in self.config['samples']])
        references = refs(base_dir)
        ref_hash, _, _, _ = refs_yaml(base_dir)
        annotation = self.config['gtf']
        dup_mark = self.config['dup_mark']
        dragen_pe = self.config['dragen_pe']

        input_list = []
        print('\n')
        print('==============================================================')
        input_list.append('==============================================================')
        print('以下の内容がDRAGEN_RNA-Seq解析作業指示書の記載内容と一致していることを確認して下さい。:')
        input_list.append('Pls cheak the following items using instruction manual:')
        print('\n')
        print(f'lotID: \"{lotID}\"')
        input_list.append(f'lotID: \"{lotID}\"')
        if (expression == True) & (fusion == False):
            print(f'analysis_type: \"expression only\"')
            input_list.append(f'analysis_type: \"expression only\"')
        elif (expression == True) & (fusion == True):
            print(f'analysis_type: \"expression + fusion\"')
            input_list.append(f'analysis_type: \"expression + fusion\"')
        print(f'library_type: \"{lib_type}\"')
        input_list.append(f'library_type: \"{lib_type}\"')
        print(f'bam_nouhin: \"{bam_nouhin}\"')
        input_list.append(f'bam_nouhin: \"{bam_nouhin}\"')
        print(f'sequence data: \"{seq_data}\"')
        input_list.append(f'Sequence data: \"{seq_data}\"')
        print(f'FASTQ file lists: \"{self.xlsx_file}\"')
        input_list.append(f'FASTQ file lists: \"{self.xlsx_file}\"')
        print(f'num of samples: \"{samples}\"')
        input_list.append(f'num of samples: \"{samples}\"')
        genome = None
        for ref in references:
            if ref_hash[ref] == self.config['ref_hash']:
                genome = ref
                print(f'reference: \"{genome}\"')
        if genome is None:
            print(f'reference: None or Unregistered')
        input_list.append(f'reference: \"{genome}\"')
        print(f'annotation: \"{annotation}\"')
        input_list.append(f'annotation: \"{annotation}\"')
        print(f'dup_mark: \"{dup_mark}\"')
        input_list.append(f'dup_mark: \"{dup_mark}\"')
        print(f'dargen version: \"{dragen_pe}\"')
        input_list.append(f'dargen version: \"{dragen_pe}\"')
        if self.qc:
            if fusion == True:
                fusion_bam = self.config['fusion_bam']
                print(f'fusion_bam: \"{fusion_bam}\"')
                input_list.append(f'fusion_bam: \"{fusion_bam}\"')
            dragen_memory = self.config['dragen_memory']
            print(f'dragen_memory: \"{dragen_memory}\"')
            input_list.append(f'dragen_memory: \"{dragen_memory}\"')
        print('\n')
        sign = input("Are all items correct? [y/n] ").strip()
        input_list.append(f'Are all items correct? [y/n] {sign}')
        self.response(sign)
        print('==============================================================')
        input_list.append('==============================================================')
        print('\n')

        return lotID, genome, input_list

class file_operation:
    def mkDir(self, out_dir, work_dir, lotID):
        logger.info('')
        logger.info('---- Log OutPut Directory ----')
        if out_dir is None:
            out_dir = os.path.join(work_dir, lotID, '000_DOC') #change
            #out_dir = os.path.abspath('.')
            if os.path.exists(out_dir):
                logger.info(f'dir_name: {out_dir}')
            else:
                logger.info(f'dir_name: {out_dir}')
                os.makedirs(out_dir, exist_ok=True)
                logger.info(f'mkdir -p {out_dir}')
        logger.info('')

        return out_dir

class cheak_yaml(file_operation):
    def __init__(self, lotID, genome, yaml_file, xlsx_file, out_dir, work_dir):
        self.lotID = lotID
        #self.analysis_type = analysis_type
        self.genome = genome
        #self.library_type = library_type
        self.yaml_file = yaml_file
        self.xlsx_file = xlsx_file
        self.out_dir = out_dir
        self.work_dir = work_dir

    def find_file(self):
        yaml_flag = 0
        xlsx_flag = 0
        if self.yaml_file is not None:
            yaml_flag = 1
            yaml_file = self.yaml_file
            if os.path.exists(yaml_file):
                logger.info(f'YAML File: {yaml_file}')
            else:
                logger.error(f'YAML File for {self.lotID} not found!')
                sys.exit()

        if self.xlsx_file is not None:
            xlsx_flag = 1
            xlsx_file = self.xlsx_file
            if os.path.exists(xlsx_file):
                logger.info(f'FASTQ File List: {xlsx_file}')
            else:
                logger.error(f'FASTQ File List for {self.lotID} not found!')
                sys.exit()

        if (yaml_flag == 0) | (xlsx_flag == 0):
            cur_dir = os.path.join(self.work_dir, self.lotID)
            yaml_name = self.lotID + '.yaml'
            xlsx_name = self.lotID + '_fastq_file_list.xlsx'

            for dirPath, dirList, files in os.walk(cur_dir):
                for f in files:
                    if self.yaml_file is None:
                        if '200_INFO' in dirPath :
                            if f == yaml_name:
                                yaml_flag = 1
                                yaml_file = os.path.join(dirPath,f)
                                logger.info(f'YAML File: {yaml_file}')

                    if self.xlsx_file is None:
                        if f == xlsx_name:
                            if '100_COMMON' in dirPath:
                                xlsx_flag = 1
                                xlsx_file = os.path.join(dirPath,f)
                                logger.info(f'FASTQ File List: {xlsx_file}')

            if yaml_flag == 0:
                logger.error(f'YAML File for {self.lotID} not found!')
                sys.exit()

            if xlsx_flag == 0:
                logger.error(f'FASTQ File List for {self.lotID} not found!')
                sys.exit()

        return yaml_file, xlsx_file

    def results(self):
        yaml_file, xlsx_file = self.find_file()
        # read FASTQ File List ---
        logger.info("")
        logger.info('Reading FASTQ file list ...')
        xlsx_df = pd.read_excel(xlsx_file)
        R1_list = xlsx_df['fastqファイル名称(Read1)'].to_list()
        R2_list = xlsx_df['fastqファイル名称(Read2)'].to_list()
        # read YAML file ---
        logger.info("")
        logger.info('Reading YAML file ...')
        logger.info(f'Input: {yaml_file}')
        config = read_yaml(yaml_file)

        # YAMLファイルの確認（lot_id）---
        logger.info("")
        logger.info('<--- Cheak YAML file --->')

        OK_flag = 0
        NG_list = []
        WARN_list = []

        logger.info('@ lotID: ')
        lotID = config['lot_id']
        logger.info(f'- pipeline yaml: {lotID}')
        logger.info(f'- your input: {self.lotID}')
        if lotID == self.lotID:
            logger.info('> RESULT: OK :) ')
        else:
            OK_flag = 1
            logger.error(f'> RESULT: NG :( ')
            NG_list.append('lotID')
        logger.info("")

        # YAMLファイルの確認（work_dirの存在）---
        logger.info('@ work_dir: ')
        workDir = config['work_dir']
        logger.info(f'- work_dir: {workDir}')
        if workDir:
            if os.path.exists(workDir):
                logger.info(f'> RESULT: Found :)')
            else:
                OK_flag = 1
                logger.error(f'> RESULT: Not found :(')
                NG_list.append('work_dir')
        else:
            OK_flag = 1
            logger.error(f'> RESULT: Not found :(')
            NG_list.append('work_dir')
        logger.info("")

        # YAMLファイルの確認（report_dirの存在）---
        logger.info('@ report_dir: ')
        reportDir = config['report_dir']
        logger.info(f'- report_dir: {reportDir}')
        if reportDir:
            if os.path.exists(reportDir):
                logger.info(f'> RESULT: Found :)')
            else:
                OK_flag = 1
                logger.info(f'> RESULT: Not Found :(')
                NG_list.append('report_dir')
        else:
            OK_flag = 1
            logger.info(f'> RESULT: Found :(')
            NG_list.append('report_dir')
        logger.info("")

        # YAMLファイルの確認（nouhin_dirの存在）---
        logger.info('@ nouhin_dir: ')
        nouhinDir = config['nouhin_dir']
        logger.info(f'- nouhin_dir: {nouhinDir}')
        if nouhinDir:
            if os.path.exists(nouhinDir):
                logger.info(f'> RESULT: Found :)')
            else:
                OK_flag = 1
                logger.info(f'> RESULT: Not Found :(')
                NG_list.append('nouhin_dir')
        else:
            OK_flag = 1
            logger.info(f'> RESULT: Not Found :(')
            NG_list.append('nouhin_dir')
        logger.info("")

        # YAMLファイルの確認（bam_nouhin）---
        logger.info('@ bam_nouhin: ')
        bamNouhin = config['bam_nouhin']
        logger.info(f'- pipeline yaml: {bamNouhin}')
        logger.info('- default: True')
        if config['bam_nouhin']:
            logger.info(f'> RESULT: Marched default setting :)')
        else:
            logger.warning(f'> RESULT: Mismarched default setting ..')
            logger.warning(f'> Please cheak bam_nouhin setting!')
            WARN_list.append('bam_nouhin')
        logger.info("")

        # YAMLファイルの確認（sampleの存在）---
        OK_cnt = 0
        NG_cnt = 0
        NG_faq = {}
        logger.info('@ samples: ')
        for sample in config['samples']:
            sample_id = sample['id']
            #name = sample['name']
            map_flag = sample['map']
            logger.info(f'- sample_id: {sample_id}')
            for lib in sample['libs']:
                for fq in lib['fastqs']:
                    _, R1 = os.path.split(fq[0])
                    _, R2 = os.path.split(fq[1])
                    if (R1 in R1_list) & (R2 in R2_list):
                        OK_cnt += 1
                        logger.info(f'- R1: {fq[0]}')
                        logger.info(f'- R2: {fq[1]}')
                        logger.info(f'> RESULT: Found in FASTQ file :)')
                    else:
                        NG_cnt += 1
                        OK_flag = 1
                        NG_faq[fq[0]] = '-'
                        NG_faq[fq[1]] = '-'
                        logger.error(f'- R1:{fq[0]}')
                        logger.error(f'- R2:{fq[1]}')
                        logger.error(f'> RESULT: No found in FASTQ file :(')
                        if NG_cnt == 1:
                            NG_list.append('samples(fastqs)')

            if map_flag:
                pass
            else:
                sample_prev = sample['prev']
                logger.warning(f'- prev:')
                if sample_prev is None:
                    logger.warning(f'> RESULT: No Setting, is it OK?')
                    WARN_list.append('samples(prev)')
                else:
                    logger.info(f'> RESULT: {sample_prev}')
        logger.info("")

        logger.info('@ num of FASTQ files:')
        num_xlsx = int(len(xlsx_df))
        logger.info(f'- # of FASTQ_file_list: {num_xlsx}')
        logger.info(f'- # of FASTQ_files in yaml: {OK_cnt}')
        if num_xlsx == OK_cnt:
            logger.info(f'> RESULT: Matched count :) ')
        else:
            OK_flag = 1
            NG_list.append('num_of_FASTQ_files')
            logger.error(f'> RESULT: Mismatched count :( ')
            for fq in NG_faq.keys():
                logger.error(f'No found in yaml :( , {fq}')
        logger.info("")

        # YAMLファイルの確認（参照配列系） ---
        ref_ha, ref_fa, ref_anno, ref_anno_txt = refs_yaml(base_dir)
        logger.info(f'@ input ref genome: {self.genome}')
        if ref_fa.get(self.genome) is None:
            if self.genome is None:
                logger.error(f'No found in SOP-33 Dragen Pipeline :(')
            else:
                logger.error(f'No found {self.genome} in SOP-33 Dragen Pipeline :(')
            logger.error(f'Pls cheak your input ..')
            logger.error(f'Pls cheak the following list of reference or annotation files ..')
            fasta = config['fasta']['file']
            refhash = config['ref_hash']
            gtf = config['gtf']
            anno_txt = config['annotation_txt']
            logger.error(f'- fasta: {fasta}')
            logger.error(f'- ref_hash: {refhash}')
            logger.error(f'- gtf: {gtf}')
            logger.error(f'- annotation_txt: {anno_txt}')
            OK_flag = 1
            NG_list.append('reference_genome')
        else:
            logger.info('@ fasta:')
            fasta = config['fasta']['file']
            logger.info(f'- fasta_file(yaml): {fasta}')
            logger.info(f'- fasta_file(SOP): {ref_fa[self.genome]}')
            if fasta == ref_fa[self.genome]:
                logger.info(f'> RESULT: Matched :) ')
            else:
                OK_flag = 1
                logger.error('> RESULT: Mismatched :( ')
                NG_list.append('fasta(reference genome)')
            logger.info("")

            logger.info('@ ref_hash:')
            refhash = config['ref_hash']
            logger.info(f'- ref_hash(yaml): {refhash}')
            logger.info(f'- ref_hash(SOP): {ref_ha[self.genome]}')
            if refhash == ref_ha[self.genome]:
                logger.info(f'> RESULT: Matched :) ')
            else:
                OK_flag = 1
                logger.error('> RESULT: Mismatched :( ')
                NG_list.append('ref_hash(reference genome)')
            logger.info("")

            logger.info('@ gtf:')
            gtf = config['gtf']
            logger.info(f'- gtf(yaml): {gtf}')
            logger.info(f'- gtf(SOP): {ref_anno[self.genome]}')
            if gtf == ref_anno[self.genome]:
                logger.info(f'> RESULT: Matched :) ')
            else:
                OK_flag = 1
                logger.error('> RESULT: Mismatched :( ')
                NG_list.append('gtf(reference genome)')
            logger.info("")

            try:
                logger.info('@ annotation_txt:')
                anno_txt = config['annotation_txt']
                logger.info(f'- annotation_txt(yaml): {anno_txt}')
                logger.info(f'- annotation_txt(SOP): {ref_anno_txt[self.genome]}')
                if gtf == ref_anno[self.genome]:
                    logger.info(f'> RESULT: Matched :) ')
                else:
                    OK_flag = 1
                    logger.error('> RESULT: Mismatched :( ')
                    NG_list.append('annotation_txt(reference genome)')
                logger.info("")
            except KeyError:
                pass

        # YAMLファイルの確認（fusion:True -> lib） ---
        if config['fusion'] :
            if 'GRCh38' in self.genome:
                logger.info('@ fusion_lib:')
                if config['fusion_lib'] is None:
                    logger.error('> RESULT: No found :( ')
                    OK_flag = 1
                    NG_list.append('fusion_lib')
                else:
                    logger.info('> RESULT: Found :) ')
                    logger.info('- items')
                    for i, j in config['fusion_lib'].items():
                        logger.info(f'{i}: {j}')
                logger.info("")

        logger.info('>--- RESULT (Cheak YAML file) ---<')
        if OK_flag == 0:
            logger.info('Everything was GOOD :)')
            logger.info('')
            return True
        else:
            logger.error('Your analysis looked BAD.. Pls cheak YAML file.')
            logger.error('')
            if len(NG_list) != 0:
                logger.error('### NG items ###')
                logger.error('Pls cheak the following items: ')
                for i, NG in enumerate(NG_list, 1):
                    logger.error(f'({i}) {NG}')
            if len(WARN_list) != 0:
                logger.error('### WARNING items ###')
                logger.error('Pls cheak the following items: ')
                for i, WARN in enumerate(WARN_list, 1):
                    logger.warning(f'({i}) {WARN}')
            logger.info('')
            return False

def check_log(work_dir: str, bam_nouhin):
    """作業フォルダ以下にあるログファイルにエラーがないか確認する関数

    Args:
        work_dir (str): 作業フォルダへのパス
        bam_nouhin: 解析変数対策
    Returns:
        bool: エラーが生じていなければTrue、エラーが検出されていればFalse
    """
    logger.info("")
    logger.info(f'<--- Cheak log files --->')
    logger.info(f'- work_dir: {work_dir}')

    error_flag = False
    for current_dir, dirs, files in os.walk(work_dir):
        for file in sorted(files):
            file_path = os.path.join(current_dir, file)
            if re.search("log", file_path):
                if bam_nouhin:
                    fh = open(file_path, "r", encoding='utf-8')
                else:
                    fh = open(file_path, "r", encoding='CP932')
                with fh as f:
                    if re.search("error", f.read(), flags=re.IGNORECASE):
                        logger.error(f"ERROR file: {file_path}")
                        error_flag = True

    if error_flag:
        logger.error(f"> RESULT: NG :(")
        logger.info("")
        return False
    else:
        logger.info(f"> RESULT: OK :)")
        logger.info("")
        return True

def check_sequence(config: Dict, sequence_dir: str, fastq_file_list: str, qc: bool):
    """sequence フォルダの確認
    Args:
        sequence_dir (str): sequenceフォルダへのパス
        fastq_file_list (str): [Lot]_fastq_file_list.xlsxへのパス

    Returns:
        bool: ディレクトリにfastq_file_list記載のfastqファイルがあればTrue、なければFalse
    """
    logger.info("")
    logger.info("<--- Cheak sequence dir --->")
    logger.info(f'- sequence_dir: {sequence_dir}')
    logger.info(f'- fastq_file_list: {fastq_file_list}')

    df = pd.read_excel(fastq_file_list)
    fastq_files = list(glob.glob(os.path.join(sequence_dir, "*.gz")))

    logger.info('@ num of FASTQ files:')
    logger.info(f'- num of FASTQ files: {len(fastq_files)}')
    logger.info(f'- num of FASTQ files in FASTQ_file_list: {len(df)*2}')

    NG_flag = False
    WARN_flag = False
    NG_list = []
    WARN_list = []

    if len(df)*2 == len(fastq_files):
        logger.info(f'> RESULT: Matched :)')
    else:
        logger.error('> RESULT: Mismatched :( ')
        read1_fastqs = df.loc[:, "fastqファイル名称(Read1)"].values.tolist()
        read2_fastqs = df.loc[:, "fastqファイル名称(Read2)"].values.tolist()
        #logger.debug(read1_fastqs)
        #logger.debug(read2_fastqs)
        NG_flag = True
        NG_list.append('num_of_fastq_file')

        for fastq in fastq_files:
            if os.path.basename(fastq) in read1_fastqs or os.path.basename(fastq) in read2_fastqs:
                pass
            else:
                logger.error(f"> RESULT: Not found {os.path.basename(fastq)} in {os.path.basename(fastq_file_list)}")
                NG_flag = True

    logger.info('@ fastq_file_list_file:')
    sequence_fastq_file_list = os.path.join(sequence_dir, os.path.basename(fastq_file_list))
    logger.info(f'- Your input: {fastq_file_list}')
    flag = 0
    seq_dir_df = None
    try:
        logger.info(f'- 900_NOUHIN: {sequence_fastq_file_list}')
        input_df = df
        seq_dir_df = pd.read_excel(sequence_fastq_file_list)
    except FileNotFoundError:
        nouhin_dir = config['nouhin_dir']
        logger.warning(f'> RESULT: No {sequence_fastq_file_list} found in {nouhin_dir}')
        WARN_flag = True
        WARN_list.append('Existence_of_fastq_file_list_file')

        datasheet_dir = os.path.join(nouhin_dir, 'datasheet')
        if not os.path.exists(datasheet_dir):
            logger.warning(f'> RESULT: No {datasheet_dir} found in {nouhin_dir}')
            WARN_flag = True
            WARN_list.append('Existence_of_datasheet_dir')
            seq_dir_df = None
        else:
            sequence_fastq_file_list = glob.glob(os.path.join(datasheet_dir, '*_datasheet.xlsx'))[0]
            logger.info(f'- 900_NOUNIN: {sequence_fastq_file_list}')
            input_df = df.loc[:, ['依頼書記載名称', 'サンプルID', 'fastqファイル名称(Read1)', 'fastqファイル名称(Read2)']]
            datasheet_df = pd.read_excel(sequence_fastq_file_list, sheet_name='sequence stats')
            seq_dir_df = datasheet_df.loc[:, ['依頼書記載名称', 'サンプルID', 'fastqファイル名（Read1）', 'fastqファイル名（Read2）']]

    if seq_dir_df:
        for i in range(len(input_df)):
            for j in range(len(input_df.columns)):
                if str(input_df.iloc[i, j]) == str(seq_dir_df.iloc[i, j]):
                    pass
                else:
                    flag = 1
                    logger.error(f'> RESULT: Mismatched :( (Your input: {input_df.iloc[i, j]}, 900_NOUHIN: {seq_dir_df.iloc[i, j]})')
                    NG_flag = True
                    NG_list.append(f'Mismatched_fastq_file_list_file')
        if flag == 0:
            logger.info(f'> RESULT: Matched :)')

    if qc:
        flag = 0
        logger.info("")
        if not os.path.exists(sequence_fastq_file_list.replace(".xlsx", ".txt")):
            logger.info(f"<--- Cheak fastq_file_list.txt (QC service) --->")
            logger.error(f'> RESULT: No found :(')
            NG_flag = True
            NG_list.append(f'Existence_of_fastq_file_list_file(QC service)')
        else:
            fastq_file_list_txt = sequence_fastq_file_list.replace(".xlsx", ".txt")
            fname = os.path.basename(fastq_file_list_txt)
            logger.info(f"<--- Cheak {fname} (QC service) --->")
            logger.info('@ fastq_file_list_file:')
            logger.info(f'- Your input: {fastq_file_list}')
            logger.info(f'- 900_NOUHIN: {fastq_file_list_txt}')
            txt_df = pd.read_csv(fastq_file_list_txt, sep='\t')
            for i in range(len(df)):
                for j in range(len(df.columns)):
                    if str(df.iloc[i, j]) == str(txt_df.iloc[i, j]):
                        pass
                    else:
                        flag = 1
                        logger.error(f'> RESULT: Mismatched :( (Your input: {df.iloc[i, j]}, 900_NOUHIN: {txt_df.iloc[i, j]})')
                        NG_flag = True
                        NG_list.append(f'Mismatched_fastq_file_list_file(QC service)')
            if flag == 0:
                logger.info(f'> RESULT: Matched :)')

    logger.info(">--- RESULT (Cheak sequence dir) ---<")
    if NG_flag:
        logger.error(f"sequence dir looked BAD..")
        logger.error("### NG items ###")
        logger.error('Pls cheak the following items: ')
        for i, NG in enumerate(list(set(NG_list)), 1):
            logger.error(f'({i}) {NG}')
        return False
    elif WARN_flag:
        logger.error('Your analysis looked BAD?')
        logger.error('### WARNING items ###')
        logger.error('Pls cheak the following items: ')
        for i, WARN in enumerate(WARN_list, 1):
            logger.warning(f'({i}) {WARN}')
        return False
    else:
        logger.info('Everything was GOOD :)')
        return True

def check_map_info(dir_path: str, fastq_file_list: str, lot: str, qc: bool):
    """マッピング結果サマリー情報の確認
    Args:
        dir_path (str): [Lot]_map_info.xlsxがあるディレクトリまでのパス
        fastq_file_list (str): [Lot]_fastq_file_list.xlsxへのパス
        lot (str): Lot ID
        qc (bool): QCサービス付き案件に対応
    Returns:
        bool: [Lot]_map_info.xlsxに問題がない場合はTrue、ファイルがない、sample idの行がない場合はFalse
    """
    logger.info("")
    logger.info(f"<--- Cheak {lot}_mapping_info.xlsx --->")
    logger.info(f'- lotID: {lot}')
    logger.info(f'- map_info_dir: {dir_path}')
    logger.info(f'- fastq_file_list: {fastq_file_list}')

    NG_flag = False
    NG_list = []
    if not os.path.exists(os.path.join(dir_path, f"{lot}_mapping_info.xlsx")):
        logger.error(f'> RESULT: Not found :( , {lot}_mapping_info.txt')
        NG_flag = True
        NG_list.append(f'Existence_of_{lot}_mapping_info.xlsx')
    else:
        map_info_path = os.path.join(dir_path, f"{lot}_mapping_info.xlsx")
        df = pd.read_excel(fastq_file_list)
        map_info_df = pd.read_excel(map_info_path)
        sample_ids = df["サンプルID"].unique().tolist()
        rows_num_fq = df.shape[0]
        rows_num_map = map_info_df.shape[0]
        logger.info('@ num of sample_id:')
        logger.info(f'- num of sample_id in {os.path.basename(fastq_file_list)}: {rows_num_fq}')
        logger.info(f'- num of sample_id in {lot}_mapping_info.xlsx: {rows_num_map}')
        if len(sample_ids) == len(map_info_df):
            logger.info('> RESULT: Matched :) ')
        else:
            logger.error('> RESULT: Mismatched :( ')
            NG_flag = True
            NG_list.append(f'sample_id_of_{lot}_mapping_info.xlsx')
            for sample_id in sample_ids:
                if sample_id in map_info_df["Sample ID"].values.tolist():
                    pass
                else:
                    logger.error(f"> RESULT: Not found :( , {sample_id}")

    if qc:
        map_info_txt_path = os.path.join(dir_path, f"{lot}_mapping_info.txt")
        logger.info("")
        logger.info(f"<--- Cheak {lot}_mapping_info.txt (QC service) --->")
        logger.info(f"{dir_path}")
        logger.info(f"{fastq_file_list}")
        logger.info(f"{lot}")
        if not os.path.exists(map_info_txt_path):
            logger.error(f"> RESULT: Not found :( , {lot}_mapping_info.txt")
            NG_flag = True
            NG_list.append(f'Existence_of_{lot}_mapping_info(QC service)')

        map_info_txt_df = pd.read_csv(map_info_txt_path, sep="\t")
        rows_num_map = map_info_txt_df.shape[0]
        logger.info('@ num of sample_id:')
        logger.info(f'- num of sample_id in {os.path.basename(fastq_file_list)}: {rows_num_fq}')
        logger.info(f'- num of sample_id in {lot}_mapping_info.txt: {rows_num_map}')
        if len(sample_ids) == len(map_info_txt_df):
            logger.info('> RESULT: Matched :) ')
        else:
            logger.error('> RESULT: Mismatched :( ')
            NG_flag = True
            NG_list.append(f'sample_id of {lot}_mapping_info(QC service)')
            for sample_id in sample_ids:
                if sample_id in map_info_txt_df["Sample ID"].values.tolist():
                    pass
                else:
                    logger.error(f"> RESULT: Not found :( , {sample_id}")

    logger.info(f">--- RESULT ({lot}_mapping_info) ---<")
    if NG_flag:
        logger.error(f"{lot}_mapping_info looked BAD..")
        logger.error("### NG items ###")
        logger.error('Pls cheak the following items: ')
        for i, NG in enumerate(NG_list, 1):
            logger.error(f'({i}) {NG}')
        return False
    else:
        logger.info('Everything was GOOD :)')
        return True

def check_expression(work_dir: str, expression_dir: str, fastq_file_list: str, lot: str):
    """expressionフォルダの確認
    Args:
        expression_dir (str): expressionフォルダへのパス
        fastq_file_list (str): [Lot]_fastq_file_list.xlsxへのパス
        lot (str): Lot ID
    Returns:
        bool: 各ファイルに[sample id].TPM、[sample id].countが含まれていればTrue、含まれていなければFalse
    """
    logger.info("")
    logger.info("<--- Cheak expression_dir --->")
    logger.info(f'- lotID: {lot}')
    logger.info(f'- expression_dir: {expression_dir}')
    logger.info(f'- fastq_file_list: {fastq_file_list}')
    NG_flag = False
    NG_list = []
    # ファイルの存在確認
    files = [
        os.path.join(expression_dir, f"{lot}_gene_expression.xlsx"),
        os.path.join(expression_dir, f"{lot}_transcript_expression.xlsx"),
        os.path.join(os.path.join(expression_dir, "supplement"), f"{lot}_gene_expression.txt"),
        os.path.join(os.path.join(expression_dir, "supplement"), f"{lot}_transcript_expression.txt")
    ]
    for file in files:
        if os.path.exists(file):
            pass
        else:
            logger.error(f"> RESULT: Not found :( , {os.path.basename(file)}")
            NG_list.append(f'Existence_of_{lot}_expression.[txt|xlsx]')
            NG_flag = True

    if NG_flag == False:
        # sample idの確認 (全検体の発現値が記載されているか)
        df = pd.read_excel(fastq_file_list)
        sample_ids = df["サンプルID"].unique().tolist()
        rows_num = {}
        for file in files:
            logger.info(f"Cheaking {file} ...")
            name = os.path.basename(file)
            file_df = pd.read_excel(file) if file.endswith(".xlsx") else pd.read_csv(file, sep="\t")
            rows_num[name] = file_df.shape[0]

            # アノテーション情報が付与されているか
            df_anno = file_df.loc[:, 'Entrez Gene ID':'Biological Process Name']
            if df_anno.empty:
                NG_flag = True
                NG_list.append('annotation')

            for sample_id in sample_ids:
                if f"{sample_id}.TPM" in file_df.columns and f"{sample_id}.count" in file_df.columns:
                    pass
                else:
                    logger.error(f"> RESULT: No found in file :( , {sample_id}")
                    NG_list.append(f'Row_num_of_{lot}_expression.[txt|xlsx]')
                    NG_flag = True

        # ファイルの行数がdat/basic_annotation_[gene|transcript].txtと一致するか
        dat_files = [
            os.path.join(os.path.join(work_dir, "dat"), "basic_annotation_gene.txt"),
            os.path.join(os.path.join(work_dir, "dat"), "basic_annotation_transcript.txt")
        ]
        for file in dat_files:
            df = pd.read_csv(file, sep="\t")
            name = os.path.basename(file)
            for fname in rows_num.keys():
                if ('gene' in name) & ('gene' in fname) :
                    col_num = df.shape[0]
                    logger.info(f'Col num of {name}: {col_num}')
                    logger.info(f'Col num of {fname}: {rows_num[fname]}')
                    if df.shape[0] == rows_num[fname]:
                        logger.info('> RESULT: Matched :) ')
                    else:
                        logger.error('> RESULT: Mismatched :( ')
                        NG_list.append(f'Col_num_of_{lot}_gene_expression.[txt|xlsx]')
                        NG_flag = True
                elif ('transcript' in name) & ('transcript' in fname) :
                    col_num = df.shape[0]
                    logger.info(f'Col num of {name}: {col_num}')
                    logger.info(f'Col num of {fname}: {rows_num[fname]}')
                    if df.shape[0] == rows_num[fname]:
                        logger.info('> RESULT: Matched :)')
                    else:
                        logger.error('> RESULT: Mismatched :(')
                        NG_list.append(f'Col_num_of_{lot}_transcript_expression.[txt|xlsx]')
                        NG_flag = True

    logger.info(f">--- RESULT (Cheak expression_dir) ---<")
    if NG_flag:
        logger.error("expression_dir looked BAD.. ")
        logger.error("### NG items ###")
        logger.error('Pls cheak the following items: ')
        for i, NG in enumerate(NG_list, 1):
            logger.error(f'({i}) {NG}')
        return False
    else:
        logger.info('Everything was GOOD :)')
        return True

def check_fusion(config: Dict, fusion_dir: str, fastq_file_list: str, qc: bool):
    """fusionフォルダの確認
    Args:
        config (dict): configファイル
        fusion_dir (str): fusionフォルダへのパス
        fastq_file_list (str): [Lot]_fastq_file_list.xlsxへのパス
        qc (bool): QCサービス付き案件に対応
    Returns:
        bool: ディレクトリにsample idのファイルがあればTrue、なければFalse
    """
    logger.info("")
    logger.info(f"<--- Cheak fusion_dir --->")
    logger.info(f'- fision_dir: {fusion_dir}')
    logger.info(f'- fastq_file_list: {fastq_file_list}')

    fusion_files = list(glob.glob(os.path.join(fusion_dir, '*.xlsx')))
    sample_ids = [sample['id'] for sample in config['samples']]

    logger.info('@ num of fusion files:')
    logger.info(f'- num of FASTQ files in FASTQ_file_list: {len(sample_ids)}')
    logger.info(f'- num of fusion files: {len(fusion_files)}')

    NG_flag = False
    NG_list = []
    if len(sample_ids) == len(fusion_files):
        logger.info('> RESULT: Matched :) ')
    else:
        logger.error('RESULT: Mismatched :( ')
        NG_flag = True
        NG_list.append(f'Existence_of_fusion.xlsx')
        for sample_id in sample_ids:
            if not os.path.exists(os.path.join(fusion_dir, f"{sample_id}_fusion.txt")):
                logger.error(f"> RESULT: Not found :( , {sample_id}_fusion.xlsx")

    if qc:
        logger.info("")
        logger.info(f"<--- Cheak fusion.txt (QC service) --->")
        fusion_txt_files = list(glob.glob(os.path.join(fusion_dir, "*.txt")))
        if len(sample_ids) == len(fusion_txt_files):
            logger.info('> RESULT: Matched :) ')
        else:
            logger.error('RESULT: Mismatched :( ')
            NG_flag = True
            NG_list.append(f'Existence_of_fusion.txt(QC service)')
            for sample_id in sample_ids:
                if not os.path.exists(os.path.join(fusion_dir, f"{sample_id}_fusion.txt")):
                    logger.error(f"> RESULT: Not found :( , {sample_id}_fusion.txt")

    logger.info(f">--- RESULT (Cheak fusion_dir) ---<")
    if NG_flag:
        logger.error(f"fusion_dir looked BAD..")
        logger.error("### NG items ###")
        logger.error('Pls cheak the following items: ')
        for i, NG in enumerate(NG_list, 1):
            logger.error(f'({i}) {NG}')
        return False
    else:
        logger.info('Everything was GOOD :)')
        return True

def check_reference(reference_dir: str, config: Dict)-> bool:
    """referenceフォルダの確認
    Args:
        reference_dir (str): referenceフォルダへのパス
        config (Dict): dictionary of yaml file
    Returns:
        bool: ディレクトリにyamlファイル記載の参照配列があればTrue、なければFalse
    """
    logger.info("")
    logger.info(f"<--- Cheak reference_dir --->")
    logger.info(f'- reference_dir: {reference_dir}')
    NG_flag = False
    NG_list = []
    if not os.path.exists(reference_dir):
        logger.error("> RESULT: No reference_dir found :(")
        NG_flag = True
        NG_list.append('Existence_of_reference_dir')
    else:
        files = sorted(os.listdir(reference_dir))
        NG_flag = False
        NG_list = []
        if len(files) == 0:
            logger.error("> RESULT: No found :(")
            NG_list.append('Existence_of_references')
        logger.info("")

        logger.info(f"@ List of reference or annotation files:")
        for file in files:
            logger.info(f'- {file}')
        logger.info("")

        fasta = os.path.join(
            reference_dir,
            os.path.basename(config["fasta"]["file"]) if not config["fasta"]["file"].endswith(".gz") else os.path.basename(config["fasta"]["file"].replace(".gz", ""))
        )
        if not os.path.exists(fasta):
            logger.error("> RESULT: No fasta found :(")
            NG_flag = True
            NG_list.append('Existence_of_fasta')

        fai = os.path.join(
            reference_dir,
            os.path.basename(config["fasta"]["file"])+".fai" if not config["fasta"]["file"].endswith(".gz") else os.path.basename(config["fasta"]["file"].replace(".gz", ""))+".fai"
        )
        if not os.path.exists(fai):
            logger.error("> RESULT: No fa.fai found :(")
            NG_flag = True
            NG_list.append('Existence_of_fa.fai')

        gene_def = os.path.join(
            reference_dir,
            os.path.basename(config["gtf"]) if config["gtf"].endswith(".gz") else os.path.basename(config["gtf"])+".gz"
        )
        if not os.path.exists(gene_def):
            logger.error("> RESULT: No gtf found :(")
            NG_flag = True
            NG_list.append('Existence_of_gtf')

        gene_def_tbi = os.path.join(
            reference_dir,
            os.path.basename(config["gtf"])+".tbi" if config["gtf"].endswith(".gz") else os.path.basename(config["gtf"])+".gz.tbi"
        )
        if not os.path.exists(gene_def_tbi):
            logger.error("> RESULT: No gtf.gz.tbi found :(")
            NG_flag = True
            NG_list.append('Existence_of_gtf.gz.tbi')

    logger.info(f">--- RESULT (Cheak reference_dir) ---<")
    if NG_flag:
        logger.error("reference_dir looked BAD..")
        logger.error("### NG items ###")
        logger.error('Pls cheak the following items: ')
        for i, NG in enumerate(NG_list, 1):
            logger.error(f'({i}) {NG}')
        return False
    else:
        logger.info('Everything was GOOD :)')
        return True

def check_document(document_dir: str, bam_nouhin: bool)-> None:
    """documentフォルダの確認する関数
    Args:
        document_dir (str): documentフォルダへのパス
    Returns:
        bool: ディレクトリに補足資料があればTrue、なければFalse
    """
    logger.info("")
    logger.info(f"<--- Cheak document_dir --->")
    logger.info(f'- document_dir: {document_dir}')
    NG_flag = False
    NG_list = []
    if not os.path.exists(document_dir):
        logger.error("> RESULT: No document_dir found :(")
        NG_flag = True
        NG_list.append('Existence_of_document_dir')
    else:
        files = sorted(os.listdir(document_dir))
        if len(files) == 0:
            logger.error("> RESULT: No found :(")
            NG_flag = False
            NG_list.append('Existence_of_pdf_files')
        logger.info("")

        logger.info(f"@ List of pdf files:")
        for file in files:
            logger.info(f'- {file}')
        logger.info("")

        if not os.path.exists(os.path.join(document_dir, "補足資料.pdf")):
            logger.error("> RESULT: No 補足資料.pdf found :(")
            NG_flag = True
            NG_list.append('Existence_of_補足資料.pdf')

        if bam_nouhin:
            if not os.path.exists(os.path.join(document_dir, "RNA-Seqデータ解析ガイド.pdf")):
                logger.error("> RESULT: No RNA-Seqデータ解析ガイド.pdf found :(")
                NG_flag = True
                NG_list.append('Existence_of_RNA-Seqデータ解析ガイド.pdf')

    logger.info(f">--- RESULT (Cheak document_dir) ---<")
    if NG_flag:
        logger.error("document_dir looked BAD..")
        logger.error("### NG items ###")
        logger.error('Pls cheak the following items: ')
        for i, NG in enumerate(NG_list, 1):
            logger.error(f'({i}) {NG}')
        return False
    else:
        logger.info('Everything was GOOD :)')
        return True

def check_fastqc(fastqc_dir: str, fastq_file_list: str)-> bool:
    """fastqcフォルダを確認する関数(QC service)
    Args:
        fastqc_dir (str): フォルダへのパス
        fastq_file_list (str): [Lot]_fastq_file_list.xlsxへのパス
    Returns:
        bool: fastqファイルに対応するhtmlがあればTrue、なければFalse
    """

    logger.info("")
    logger.info("<--- Cheak fastqc_dir --->")
    logger.info(f'- fastqc_dir: {fastqc_dir}')
    logger.info(f'- fastq_file_list: {fastq_file_list}')

    df = pd.read_excel(fastq_file_list)
    fastqc_files = list(glob.glob(os.path.join(fastqc_dir, "*.html")))

    logger.info('@ num of FastQC files:')
    logger.info(f'- num of FastQC files: {len(fastqc_files)}')
    logger.info(f'- num of FASTQ files in FASTQ_file_list: {len(df)*2}')

    NG_flag = False
    NG_list = []
    if len(df)*2 == len(fastqc_files):
        logger.info(f'> RESULT: Matched :)')
    else:
        logger.error('> RESULT: Mismatched :( ')
        NG_flag = True
        NG_list.append('num_of_FastQC_files')
        for i in range(len(df)):
            for col_name in ["fastqファイル名称(Read1)", "fastqファイル名称(Read2)"]:
                fastqc_result = os.path.join(
                    fastqc_dir,
                    os.path.join(fastqc_dir, df.loc[i, col_name].replace(".fastq.gz","_fastqc.html"))
                )
                if not os.path.exists(fastqc_result):
                    logger.error(f"> RESULT: Not found :( , {fastqc_result}")

    logger.info(f">--- RESULT (Cheak fastqc_dir) ---<")
    if NG_flag:
        logger.error("fastqc_dir looked BAD..")
        logger.error("### NG items ###")
        logger.error('Pls cheak the following items: ')
        for i, NG in enumerate(NG_list, 1):
            logger.error(f'({i}) {NG}')
        return False
    else:
        logger.info('Everything was GOOD :)')
        return True

def check_mapping(mapping_dir: str, fastq_file_list: str, bam_nouhin: bool)-> None:
    """mappingフォルダの確認
    Args:
        mapping_dir (str): mappingフォルダへのパス
        fastq_file_list (str): [Lot]_fastq_file_list.xlsxへのパス
    Returns:
        bool: ディレクトリにsample idに対応するbamがあればTrue、なければFalse
    """
    if bam_nouhin:
        logger.info("")
        logger.info("<--- Cheak mapping_dir --->")
        logger.info(f'- mapping_dir: {mapping_dir}')
        logger.info(f'- fastq_file_list: {fastq_file_list}')

        df = pd.read_excel(fastq_file_list)
        sample_ids = df["サンプルID"].unique().tolist()
        bam_files = list(glob.glob(os.path.join(mapping_dir, "*.bam")))
        bai_files = list(glob.glob(os.path.join(mapping_dir, "*.bai")))

        logger.info('@ num of mapping files:')
        logger.info(f'- num of FASTQ files in FASTQ_file_list: {len(df)}')
        logger.info(f'- num of bam files: {len(bam_files)}')
        NG_flag = False
        NG_list = []
        if len(sample_ids) == len(bam_files):
            logger.info(f'> RESULT: Matched :)')
        else:
            logger.error(f'> RESULT: Mismatched :)')
            NG_flag = True
            NG_list.append('num_of_bam_files')
            for sample_id in sample_ids:
                bam = os.path.join(mapping_dir, f"{sample_id}.bam")
                if not os.path.exists(bam):
                    logger.error(f"> RESULT: Not found :( , {bam}")

        logger.info(f'- num of bai files: {len(bai_files)}')
        if len(sample_ids) == len(bai_files):
            logger.info(f'> RESULT: Matched :)')
        else:
            logger.error(f'> RESULT: Mismatched :)')
            NG_flag = True
            NG_list.append('num_of_bai_files')
            for sample_id in sample_ids:
                bai = os.path.join(mapping_dir, f"{sample_id}.bam.bai")
                if not os.path.exists(bai):
                    logger.error(f"> RESULT: Not found :( , {bai}")

        logger.info(f">--- RESULT (Cheak mapping_dir) ---<")
        if NG_flag:
            logger.error("mapping_dir looked BAD..")
            logger.error("### NG items ###")
            logger.error('Pls cheak the following items: ')
            for i, NG in enumerate(NG_list, 1):
                logger.error(f'({i}) {NG}')
            return False
        else:
            logger.info('Everything was GOOD :)')
            return True

def show_param(args: argparse.Namespace):
    """入力として与えられたコマンドライン引数を表示
    Args:
        args (argparse.Namespace): コマンドライン引数パース結果
    """
    logger.info("==============================================================")
    logger.info(f"python: {platform.python_version()}")
    scripts = os.path.join(__file__, os.path.basename(sys.argv[0]))
    logger.info(f'program: {scripts}')
    logger.info(f"user: {getpass.getuser()}")
    logger.info(f"current directory: {Path.cwd()}")
    logger.info("--------------------------------------------------------------")
    logger.info("                    Command Line Arguments                    ")
    for i, (key, value) in enumerate(vars(args).items(), 1):
        logger.info(f'args_{i} ({key}) : {value}')
    logger.info("==============================================================")
    logger.info("")

def parameters(__desc__):
    parser = argparse.ArgumentParser(
        description = __desc__,
        #formatter_class=argparse.ArgumentDefaultsHelpFormatter #default
    )
    abspath = os.path.abspath('.')
    parser.add_argument(
        '-y',
        dest='yaml_file',
        help='input yaml file (default: ' + \
        '/NGSWORK/PROJECT/[lotID]/100_COMMON/200_INFO/[lotID].yaml)',
        default=None
    )
    parser.add_argument(
        '-f',
        dest='xlsx_file',
        help='input FASTQ file list (default: ' + \
        '/NGSWORK/PROJECT/[lotID]/100_COMMON/110_GAII/[lotID]_fastq_file_list.xlsx)',
        default=None
    )
    parser.add_argument(
        '-qc',
        dest='qc',
        action='store_true',
        default=False,
        help='option for QC service.'
    )
    parser.add_argument(
        '-o',
        dest='out_dir',
        help='input result directory (default:'+abspath+')',
        default=None
    )
    parser.add_argument(
        '-pd',
        dest='project_dir',
        help='input project directory (default: ' + \
        '/NGSWORK/PROJECT)', #/NGSWORK/PROJECT
        default='/NGSWORK/PROJECT'
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

def main(args):
    # setup parameters ---
    yaml_file = args.yaml_file
    xlsx_file= args.xlsx_file
    out_dir = args.out_dir
    project_dir = args.project_dir
    qc = args.qc
    config = read_yaml(yaml_file)

    # cheak Command Line Arguments ---
    if not os.path.exists(yaml_file):
        logger.error('')
        logger.error(f'No found yaml file ({yaml_file})..')
        logger.error('Pls cheak your input.')
        logger.error('')
        sys.exit()
    if not os.path.exists(xlsx_file):
        logger.error('')
        logger.error(f'No found fastq_file_list ({xlsx_file})..')
        logger.error('Pls cheak your input.')
        logger.error('')
        sys.exit()

    # show parameters ---
    show_param(args)

    # cheak instruction manual ---
    lotID, genome, input_list = display_cheak(yaml_file, xlsx_file, qc).cheak_args()

    logger.info('')
    for term in input_list:
        logger.info(f'{term}')

    # cheak yaml ---
    out_dir =  file_operation().mkDir(out_dir, project_dir, lotID)
    yaml_flag = cheak_yaml(lotID, genome, yaml_file, xlsx_file, out_dir, project_dir).results()

    # cheak log ---
    error_flag = check_log(config['work_dir'], config['bam_nouhin'])

    # cheak outputs in 900_NOUHIN ---
    nouhin_lot_dir = config['nouhin_dir']
    if not os.path.exists(os.path.join(nouhin_lot_dir, 'sequence')):
        logger.info("")
        logger.info("<--- Cheak sequence dir --->")
        logger.warning(f'> RESULT: No found :(')
        sequence_flag = False
    else:
        sequence_flag = check_sequence(config, os.path.join(nouhin_lot_dir, 'sequence'), xlsx_file, qc)
    if qc:
        if not os.path.exists(os.path.join(nouhin_lot_dir, 'fastqc')):
            logger.info("<--- Cheak fastqc_dir --->")
            logger.warning(f'> RESULT: No found :(')
            fastqc_flag = False
        else:
            fastqc_flag = check_fastqc(os.path.join(nouhin_lot_dir, 'fastqc'), xlsx_file)

    if config['bam_nouhin']:
        map_info_flag = check_map_info(os.path.join(nouhin_lot_dir, 'mapping'), xlsx_file, lotID, qc)
    else:
        if config['expression']:
            map_info_flag = check_map_info(os.path.join(nouhin_lot_dir, 'expression'), xlsx_file, lotID, qc)
        elif config['fusion']:
            map_info_flag = check_map_info(os.path.join(nouhin_lot_dir, 'fusion'), xlsx_file, lotID, qc)

    mapping_flag = check_mapping(os.path.join(nouhin_lot_dir, 'mapping'), xlsx_file, config['bam_nouhin'])

    if (config['expression']) & (config['fusion']):
        expression_flag = check_expression(config['work_dir'], os.path.join(nouhin_lot_dir, 'expression'), xlsx_file, lotID)
        fusion_flag = check_fusion(config, os.path.join(nouhin_lot_dir, 'fusion'), xlsx_file, qc)
    else:
        if config['expression']:
            expression_flag = check_expression(config['work_dir'], os.path.join(nouhin_lot_dir, 'expression'), xlsx_file, lotID)
        elif config['fusion']:
            fusion_flag = check_fusion(config, os.path.join(nouhin_lot_dir, 'fusion'), xlsx_file, qc)

    reference_flag = check_reference(os.path.join(nouhin_lot_dir, 'reference'), config)
    document_flag = check_document(os.path.join(nouhin_lot_dir, 'document'), config['bam_nouhin'])

    if config['expression'] & config['fusion']:
        if os.path.exists(os.path.join(nouhin_lot_dir, 'mapping')):
            mapping_dir = os.path.join(nouhin_lot_dir, 'mapping')
        else:
            mapping_dir = os.path.join(config['work_dir'], config['dragen_subdir'])
        report_flag = cheak_fusion_report(
            config['report_dir'], xlsx_file, mapping_dir,
            os.path.join(config['work_dir'], config['fusion_subdir']), lotID
        ).result()
    else:
        if config['expression']:
            if config['bam_nouhin']:
                if os.path.exists(os.path.join(nouhin_lot_dir, 'mapping')):
                    mapping_dir = os.path.join(nouhin_lot_dir, 'mapping')
                else:
                    mapping_dir = os.path.join(config['work_dir'], config['dragen_subdir'])
            else:
                if os.path.exists(os.path.join(nouhin_lot_dir, 'expression')):
                    mapping_dir = os.path.join(nouhin_lot_dir, 'expression')
                else:
                    mapping_dir = os.path.join(config['work_dir'], config['dragen_subdir'])
            report_flag = cheak_expression_report(config['report_dir'], xlsx_file, mapping_dir, lotID).result()
        elif config['fusion']:
            if config['bam_nouhin']:
                if os.path.exists(os.path.join(nouhin_lot_dir, config['dragen_subdir'])):
                    mapping_dir = os.path.join(nouhin_lot_dir, config['dragen_subdir'])
                else:
                    mapping_dir = os.path.join(config['work_dir'], config['dragen_subdir'])
            else:
                if os.path.exists(os.path.join(nouhin_lot_dir, 'fusion')):
                    mapping_dir = os.path.join(nouhin_lot_dir, 'fusion')
                else:
                    mapping_dir = os.path.join(config['work_dir'], config['dragen_subdir'])
            report_flag = cheak_fusion_report(
                config['report_dir'], xlsx_file, mapping_dir,
                os.path.join(config['work_dir'], config['fusion_subdir']), lotID
            ).result()

    logger.info('')
    logger.info('==============================================================')
    logger.info('                      Result Summary                          ')
    logger.info('==============================================================')
    logger.info('@ Cheak yaml file')
    logger.info('> RESULT: GOOD :)') if yaml_flag else logger.info('> RESULT: BAD :( Pls cheak yaml result! (WARNING or ERROR?)')
    logger.info('')
    logger.info('@ Cheak log files')
    logger.info('> RESULT: GOOD :)') if error_flag else logger.info('> RESULT: BAD :( Pls cheak log result!')
    logger.info('')
    logger.info('@ Cheak outputs')
    logger.info('> sequence_dir: GOOD :)') if sequence_flag else logger.info('> sequence_dir: BAD :( Pls cheak sequence_dir result! (WARNING or ERROR?)')
    if qc:
        logger.info('> fastqc_dir: GOOD :)') if fastqc_flag else logger.info('> fastqc_dir: BAD :( Pls cheak fastqc_dir result!')

    if config['bam_nouhin'] is True:
        logger.info('> mapping_dir: GOOD :)') if mapping_flag else logger.info('> mapping_dir: BAD :( Pls cheak mapping_dir result!')

    logger.info('> map_info: GOOD :)') if map_info_flag else logger.info('> map_info: BAD :( Pls cheak map_info result!')

    if (config['expression']) & (config['fusion']):
        logger.info('> expression_dir: GOOD :)') if expression_flag else logger.info('> expression_dir: BAD :( Pls cheak expression_dir result!')
        logger.info('> fusion_dir: GOOD :)') if fusion_flag else logger.info('> fusion_dir: BAD :( Pls cheak fusion_dir result!')
    else:
        if config['expression']:
            logger.info('> expression_dir: GOOD :)') if expression_flag else logger.info('> expression_dir: BAD :( Pls cheak expression_dir result!')
        elif config['fusion']:
            logger.info('> fusion_dir: GOOD :)') if fusion_flag else logger.info('> fusion_dir: BAD :( Pls cheak fusion_dir result!')

    logger.info('> reference_dir: GOOD :)') if reference_flag else logger.info('> reference_dir: BAD :( Pls cheak reference_dir result!')
    logger.info('> document_dir: GOOD :)') if document_flag else logger.info('> document_dir: BAD :( Pls cheak document_dir result!')
    logger.info('> report: GOOD :)') if report_flag else logger.info('> report: BAD :( Pls cheak report result!')
    logger.info('')
    logger.info('Done!')
    logger.info('')
    os.system(f'mv {LOGFILE} {out_dir}')

if __name__ == "__main__":
    __version__ = '1.0' # プログラムのバージョン
    __desc__ = 'Cheak tool for outputs from DRAGEN RNA-Seq pipeline'
    parser = argparse.ArgumentParser(
        description=__desc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parm = parameters(__desc__)

    logger = getLogger(__name__)
    logger.setLevel(parm.loglevel)
    FORMAT = "%(levelname)s: %(message)s"
    #FORMAT = '%(asctime)s %(name)s %(levelname)s: %(message)s'
    dt_fmt = '%Y-%m-%d %H:%M:%S'
    formatter = Formatter(FORMAT,dt_fmt)

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