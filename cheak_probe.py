import os
import re
import sys
import argparse
import datetime
import openpyxl
import gzip
import logging
from logging import getLogger, StreamHandler, FileHandler, Formatter

LOGTIME = datetime.date.today().strftime('%y%m%d')

def parameters(desc):
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('-ib',
                        dest='input_bed',
                        help='input your probe bed file'
                        )
    parser.add_argument('-rb',
                        dest='ref_bed',
                        help='input reference probe bed file'
                        )
    parser.add_argument('-rg',
                        dest='ref_gtf',
                        help='input reference gtf file'
                        )
    parser.add_argument('-l',
                        dest='log_level',
                        help='choose log level (default: INFO)',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        default='INFO'
                        )

    args = parser.parse_args()

    return args

class cheak_probe_location:
    def __init__(self, ref_gtf, ref_bed, input_bed):
        self.ref_gtf = ref_gtf
        self.ref_bed = ref_bed
        self.input_bed = input_bed

    def get_CDS(self):
        """[summary]

        Returns:
            [type]: [description]
        """
        path, ext = os.path.splitext(self.ref_gtf)

        if ext == '.gz':
            fh = gzip.open(self.ref_gtf, 'rt')
            logger.info(f'reference gtf file: {self.ref_gtf}')
        elif ext == '.gtf':
            fh = open(self.ref_gtf, 'r')
            logger.info(f'reference gtf file: {self.ref_gtf}')
        else:
            logger.error(f'pls input gtf file! --> your input file {self.ref_gtf}')
            sys.exit()

        CDS = {}
        with fh as gtf:
            for line in gtf.readlines():
                line = line.rstrip('\n')
                tmp=line.split('\t')
                if tmp[2] == 'CDS':
                    ele = tmp[8].split('; ')
                    gene_id = re.compile(r'gene_id\s\"(.*)\"').search(ele[0]).group(1)
                    trans_id = re.compile(r'transcript_id\s\"(.*)\"').search(ele[1]).group(1)
                    ID = gene_id + ',' + trans_id
                    locus = tmp[0] + ':' + tmp[3] + '-' + tmp[4] # exon番号を追加する可能性あり
                    if CDS.get(ID):
                        CDS[ID] += ',' + locus
                    else:
                        CDS[ID] = locus

        return CDS

    def info_IDT_exome_probe(self):
        """[summary]

        Args:
            ref_bed ([type]): [description]

        Returns:
            [type]: [description]
        """
        logger.info(f"IDT_exome_probe: {self.ref_bed}")

        IDT = {}
        with open(self.ref_bed, 'r') as rb:
            for line in rb.readlines():
                line = line.rstrip('\n')
                tmp = line.split()
                gene = re.compile(r'(.*)_(\d+)').search(tmp[3]).group(1)
                locus = tmp[0] + ':' + tmp[1] + '-' + tmp[2]
                if IDT.get(gene):
                    IDT[gene] += locus
                else:
                    IDT[gene] = locus

        return IDT

    def result(self):
        """[summary]
        """
        CDS = self.CDS
        IDT = self.IDT
        MIS_list = []
        BAD_list = []

        logger.info(f'input your probe bed file:{self.input_bed}')

        # 結果の比較部（結論までの過程）
        tmp_file1 = './tmp1.txt'
        with open(tmp_file1, 'w') as wf:
            wf.write('quary_info\tquary_probe\tCDS_region\tIDT_probe\tcheak(OK/NG)\n')
            with open(self.input_bed, 'r') as ib:
                for line in ib.readlines():
                    line = line.rstrip('\n')
                    tmp = line.split()
                    # ファイルの中身に合わせて、データの処理の仕方を変更　
                    # ->　作成されたプローブの位置が何番目のexon上にあるかの情報の有無を確認
                    # 有 -> get_CDSの変数locusにexon番号を追加、
                    # 無 -> しらみ潰しに確認して(splitで崩して、for分を回す)、該当するCDS領域が1つに収束したらOK
                    gene = re.compile(r'(.*)_(\d+)').search(tmp[3]).group(1)
                    quary_ID = gene + ',' + tmp[5]
                    wf.write('\t')
                    if CDS.get(quary_ID):
                        quary_chr = tmp[0]
                        quary_start = tmp[1]
                        quary_end = tmp[2]

                        """ファイルの中身に合わせて、データの処理の仕方を変更
                        if :
                            wf.write('OK\n')
                        else:
                            logger.error(f'{quary_info} not observed on target CDS().. Pls cheak probe location.')
                            wf.write('NG\n')
                        """

                    else:
                        logger.error(f'{quary_info} not found in annotation file.. Pls cheak probe location.')
                        wf.write('NG\n')
                        BAD_list.append(quary_info)

        # 結果のまとめの文章（結論）
        tmp_file2 = './tmp2.txt'
        with open(tmp_file2, 'w') as wf:
            wf.write('RESULTS')

        result_file = './tmp.txt'
        os.system(f'cat {tmp_file2} {tmp_file1} > {result_file}')
        os.system(f'rm -f {tmp_file2} {tmp_file1}')

        # excelに結果を書き込む
        xlsx_file = './cheak_probe' + LOGTIME + '_result.xlsx'
        logger.info(f'result file:{xlsx_file}')
        wb = openpyxl.Workbook()
        sheet = wb.active
        sheet.title = 'result_' + LOGTIME
        #tmp.txtの内容を張り出す
        with open(result_file, 'r') rf:
            for line in ib.readlines():
                line = line.rstrip('\n')

        wb.save(xlsx_file)

def main(args):
    # setup parameters

    # run process


if __name__ == '__main__':
    __version__ = '1.0'
    __desciption__ = 'Some useful program commands are:'

    parm = parameters(__desciption__)

    logger = getLogger('__main__')
    logger.setLevel(parm.loglevel)

    handler_format = Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    stream_handler = StreamHandler()
    stream_handler.setLevel(logging.DEBUG)
    stream_handler.setFormatter(handler_format)

    logfile = './cheak_probe' + LOGTIME + '_log.txt'
    file_handler = FileHandler(filename=logfile, encoding='utf-8')
    file_handler.setLevel(logging.loglevel)
    file_handler.setFormatter(handler_format)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    main(args)