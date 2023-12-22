import os
import sys
import glob
import gzip
import datetime
import platform
import getpass
import argparse
import numpy as np
import pandas as pd
import anndata
import hdf5plugin
from anndata import (
    AnnData,
    read_csv,
    read_text,
    read_excel,
    read_mtx,
    read_loom,
    read_hdf,
)
import scanpy as sc
from anndata import read as read_h5ad
from io import BytesIO
from pathlib import Path
from scipy.io import mmwrite
from scipy.sparse import csr_matrix, coo_matrix
from logging import Logger, getLogger, Formatter, StreamHandler, FileHandler
import warnings

warnings.filterwarnings('ignore')

DATE = datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
LOGFILE = os.path.join(
    os.path.abspath('.'),
    os.path.basename(sys.argv[0]).replace('.py', '') + '_' + DATE + '.log'
)

def read_gtf_info(gtf_file, feature):
    """_summary_

    Args:
        gtf_file (_type_): _description_
        feature (_type_): _description_

    Returns:
        _type_: _description_
    """

    keys = ['gene_id', 'intron_id', 'gene_name', 'gene_biotype']
    if feature == 'transcript':
        keys = ['transcript_id', 'transcript_name', 'transcript_biotype']
    gtfInfo = {key:[] for key in keys}

    try:
        with open(gtf_file, 'r') as inp:
            for line in inp:
                prev = 'NA'
                line = line.rstrip().split('\t')
                if line[0].startswith('#'):
                      continue
                elif line[2] == feature:
                    annot = line[8].split(';') 
                    id = ''
                    name = ''
                    type = ''
                    for i in range(0, len(annot)):
                        entry = annot[i]
                        if feature == 'transcript':    
                            if 'transcript_id' in entry: id = entry.split('"')[1]
                            if 'transcript_name' in entry: name = entry.split('"')[1]
                            if 'transcript_biotype' in entry: type = entry.split('"')[1]
                        else:
                            if 'gene_id' in entry: id = entry.split('"')[1]
                            if 'gene_name' in entry: name = entry.split('"')[1]
                            if 'gene_biotype' in entry: type = entry.split('"')[1]
                        
                    if id == prev:
                        continue
                    prev = id
                    
                    if id == '':
                        continue
                    
                    append_data = {'gene_id': id, 'intron_id': id + '-I', 'gene_name': name, 'gene_biotype': type}
                    if feature == 'transcript':
                        append_data = {'transcript_id': id, 'transcript_name': name, 'transcript_biotype': type}
                    for key in keys:
                        gtfInfo[key].append(append_data[key])
        
        df_gtf = pd.DataFrame(gtfInfo)

    except EnvironmentError as err:
        print(f'Unable to open files: {gtf_file}')
    
    return df_gtf

"""
def make_DataFrame(gtf_file, featurefile, tx2genefile):
    #DataFrame for all genes and transcripts

    #Args:
    #    gtf_file (str): GTF file
    #    featurefile (str): File with corresponding spliced and unspliced gene IDs
    #    tx2genefile (str): File containing a mapping of transcripts to genes

    #Returns:
    #    DataFrame: retern gene and transcript list
    
    # read gtf ---
    df_gtf_gene = read_gtf_info(gtf_file, 'gene')
    df_gtf_tx = read_gtf_info(gtf_file, 'transcript')

    # gene ---
    df_feat = pd.read_csv(featurefile, sep='\t')
    df_gene = pd.merge(df_feat, df_gtf_gene, left_on='spliced', right_on='gene_id', how='left')
    df_gene['intron_name'] = df_gene['gene_name'] + '-I'
    df_gene = df_gene.drop(columns=['gene_id'])

    # tx ---
    df_tx_all = pd.read_csv(tx2genefile, sep='\t', header=None)
    df_tx_spl = df_tx_all[0][~df_tx_all[0].str.contains('-I')].reset_index(drop=True).rename('tx')
    df_tx_spliced = pd.merge(df_tx_spl, df_gtf_tx, left_on='tx', right_on='transcript_id', how='left')
    df_tx_spliced = df_tx_spliced.drop(columns=['transcript_id'])

    ### the following line should be changed to transcript_name and transcript_id ---
    # df_gene_uns = df_tx_all[1][df_tx_all[1].str.contains('-I')].reset_index(drop=True).rename('unspliced_gene')
    # df_gene_intron = pd.merge(df_gene_uns, df_gene[['intron', 'gene_name']], left_on='unspliced_gene', right_on='intron', how='left')
    # df_gene_intron = df_gene_intron.drop(columns=['intron'])

    # df_tx_ins = df_tx_all[df_tx_all[1].str.contains('-I')].reset_index(drop=True)
    # df_tx_ins.columns = ['tx', 'unspliced_gene']
    # df_tx_unspliced = pd.merge(df_tx_ins, df_gene_intron, on='unspliced_gene', how='left')
    # df_tx_unspliced = df_tx_unspliced.drop(columns=['unspliced_gene'])
    # df_tx_unspliced['gene_name'] = df_tx_unspliced['gene_name'] + '-I'
    # df_tx_unspliced = df_tx_unspliced.rename(columns={'gene_name': 'transcript_name', 'gene_biotype': 'transcript_biotype'})
    ###

    df_tx = pd.concat([df_tx_spliced, df_tx_unspliced], axis=0, ignore_index=True)

    return df_gene, df_tx
"""
    
def count_genes(df_gene, salmon_dir, outdir, prefix='', csv=True):
    """Count matrix for genes and reteined introns
    salmon_dir = '/Users/TomoyaUchiyama/WORKSPACE/salmon'
    
    Args:
        df_gene (DataFrame): gene list including retained introns
        
        salmon_dir (str): Path to the directory that has separate directories for each barcode, 
                        each containing an quant.genes.sf file 
        
        outdir (str): Path to output directory
        
        prefix (str, optional): Any prefix before `matrix.mtx`, `genes.tsv` and `barcodes.tsv`. For instance,
                            if the files are named `patientA_matrix.mtx`, `patientA_genes.tsv` and
                            `patientA_barcodes.tsv` the prefix is `patientA_`. Defaults to ''. Defaults to ''.

    Returns:
        NumPy : gene and reteined intron matrix
    """
    outdir = os.path.join(outdir, 'output')
    os.makedirs(outdir, exist_ok=True)

    genefile = os.path.join(outdir, f'{prefix}genes.tsv')
    intronfile = os.path.join(outdir, f'{prefix}introns.tsv')
    df_gene[['spliced', 'gene_name', 'gene_biotype']].to_csv(genefile, sep='\t', index=False)
    df_gene[['intron', 'intron_name', 'gene_biotype']].to_csv(intronfile, sep='\t', index=False)
    
    gene_sfs = sorted(glob.glob(os.path.join(salmon_dir, '*/quant.genes.sf')))
    bcfile = os.path.join(outdir, f'{prefix}barcodes.tsv')

    mtx_gene = np.empty((len(gene_sfs), len(df_gene['spliced'])))
    mtx_intron = np.empty((len(gene_sfs), len(df_gene['intron'])))
    with open(bcfile, 'w') as fh:
        for i, gene_sf in enumerate(gene_sfs):
            bc = os.path.dirname(gene_sf).split('/')[-1]
            fh.write(f'{bc}\n')

            df_qnt = pd.read_csv(gene_sf, sep='\t')
            df_spliced_genes = df_qnt[~df_qnt['Name'].str.contains('-I')].reset_index(drop=True)
            df_unspliced_genes = df_qnt[df_qnt['Name'].str.contains('-')].reset_index(drop=True)    # // modify HNF1A-AS1

            df_gene_cnt = pd.concat([
                df_gene['spliced'],
                df_spliced_genes[['Name', 'NumReads']].set_index('Name').
                reindex(df_gene['spliced'].values).reset_index(drop=True)
            ], axis=1).fillna(0)['NumReads'].to_numpy()
            mtx_gene[i] = df_gene_cnt

            df_intron_cnt = pd.concat([
                df_gene['intron'],
                df_unspliced_genes[['Name', 'NumReads']].set_index('Name').
                reindex(df_gene['intron'].values).reset_index(drop=True)
            ], axis=1).fillna(0)['NumReads'].to_numpy()
            mtx_intron[i] = df_intron_cnt
    
    if csv:
        barcode = pd.read_csv(bcfile, sep='\t', header=None)
        gene_cnts = pd.DataFrame(mtx_gene.T, columns=barcode[0])
        genes = pd.read_csv(genefile, sep='\t')
        GeneID = genes['gene_id'] + '_' + genes['gene_name']
        GeneID.name = 'GeneID'
        mtx_gene = pd.concat([GeneID, gene_cnts], axis=1)

        barcode = pd.read_csv(bcfile, sep='\t', header=None)
        intron_cnts = pd.DataFrame(mtx_intron.T, columns=barcode[0])
        mtx_intron = pd.concat([GeneID, intron_cnts], axis=1)

    return mtx_gene, mtx_intron

"""
outfile = '/wgbs_global1/davinci/tuchiyama/A549/output_combine_decoy/genes_mtx.csv'
mtx_gene.to_csv(outfile, index=False)

outfile = '/wgbs_global1/davinci/tuchiyama/A549/output_combine_decoy/introns_mtx.csv'
mtx_intron.to_csv(outfile, index=False)
"""

def count_txs(df_tx_all, salmon_dir, outdir, prefix=''):
    """_summary_
    # file = '/Users/TomoyaUchiyama/WORKSPACE/salmon/Homo_sapiens.GRCh38.94.expanded.tx2gene.tsv'
    Args:
        tx2gene (DaraFrame): transcript list
        
        salmon_dir (str): Path to the directory that has separate directories for each barcode, 
                        each containing an quant.sf file 
        
        outdir (str): Path to output directory

        prefix (str, optional): Any prefix before `matrix.mtx`, `genes.tsv` and `barcodes.tsv`. For instance,
                            if the files are named `patientA_matrix.mtx`, `patientA_genes.tsv` and
                            `patientA_barcodes.tsv` the prefix is `patientA_`. Defaults to ''.

    Returns:
        NumPy : transcript matrix
    """
    df_tx = df_tx_all[~df_tx_all['tx'].str.contains('-I')].reset_index(drop=True)

    genefile = os.path.join(outdir, f'{prefix}transcripts.tsv')
    df_tx.to_csv(genefile, sep='\t', index=False, header=None)

    tx_sfs = sorted(glob.glob(os.path.join(salmon_dir, '*/quant.sf')))
    mtx_tx = np.empty((len(tx_sfs), len(df_tx)))
    for i, tx_sf in enumerate(tx_sfs):
        df_qnt = pd.read_csv(tx_sf, sep='\t')
        df_spliced_txs = df_qnt[~df_qnt['Name'].str.contains('-I')].reset_index(drop=True)
        
        df_tx_cnt = pd.concat([
            df_tx['tx'],
            df_spliced_txs[['Name', 'NumReads']].set_index('Name').
            reindex(df_tx['tx'].values).reset_index(drop=True)
        ], axis=1).fillna(0)['NumReads'].to_numpy()
        mtx_tx[i] = df_tx_cnt

    return mtx_tx

def output_CooMatrix(counts, outgzfile):
    """Output a sparse matrix in Coordinate format

    Args:
        counts (NumPy): count matrix
        outgzfile (str): output file (e.g. /AAA/BBB/[genes|introns|transcripts].mtx.gz)
    """
    target = BytesIO()
    mmwrite(target, coo_matrix(counts)) # transpose
    with gzip.open(outgzfile, 'wt') as fh:
        fh.write(target.getvalue().decode('latin1'))

def createAnnData(
    mtxfile,
    cellgroup_file,
    path,
    var_names = 'gene_symbols',
    make_unique = True,
    prefix=''
):
    """_summary_
    var_names = 'gene_symbols'
    make_unique = True
    prefix = ''
    path = '/Users/TomoyaUchiyama/WORKSPACE/data/filtered_gene_bc_matrices/test'

    Args:
        mtxfile (str): target matrix file (e.g. [genes|introns|transcripts].mtx.gz)
        path (str): Path to directory for `.tsv` files,
        var_names (str, optional): The variables index. Defaults to 'gene_symbols'.
        make_unique (bool, optional): make sure gene names are unique. Defaults to True.
        prefix (str, optional): Any prefix before `matrix.mtx`, `genes.tsv` and `barcodes.tsv`. For instance,
        if the files are named `patientA_matrix.mtx`, `patientA_genes.tsv` and
        `patientA_barcodes.tsv` the prefix is `patientA_`. Defaults to ''.

    Returns:
        AnnData: Annotated data matrix (stores observations (samples) of variables/features in the rows of a matrix)
    """
    adata = read_mtx(mtxfile) # convert coo_matrix to csr_matrix during reading mtx file
    path = Path(path)
    genes = pd.read_csv(path / f'{prefix}genes.tsv', sep='\t')
    var_names = genes['gene_name'].values # var_names = genes[1].values
    if make_unique:
        var_names = anndata.utils.make_index_unique(pd.Index(var_names))
    adata.var_names = var_names
    adata.var['gene_ids'] = genes['spliced'].values # adata.var['gene_ids'] = genes[0].values
    adata.obs_names = pd.read_csv(path / f'{prefix}barcodes.tsv', header=None)[0].values

    df_barcode = pd.read_csv(path / f'{prefix}barcodes.tsv', header=None)
    adata.obs_names = df_barcode[0].values
    if cellgroup_file:
        df_cellgroup = pd.read_csv(cellgroup_file)
        df_merge = pd.merge(df_barcode, df_cellgroup, left_on=0, right_on='Barcode', how='left')['BC1'] # should be changed column name 
        adata.obs['clusters'] = pd.Categorical(df_merge) 

    return adata

def AnnSplicedscVelo(gene_adata, intron_adata):
    """Add the spliced and unspliced layers to the created Anndata object 
    for RNA velocity estimation using scVelo

    Args:
        gene_adata (AnnData): Annotated data matrix for genes
        intron_adata (AnnData): Annotated data matrix fir retained introns

    Returns:
        AnnData: Annotated data matrix with the spliced and unspliced layers
    """
    gene_adata.layers['spliced'] = gene_adata.X
    gene_adata.layers['unspliced'] = intron_adata.X
    return gene_adata


def ProcessingByScanpy(adata):
    """_summary_

    Args:
        adata (_type_): _description_
    """
    #sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
    #sc.logging.print_header()
    #sc.settings.set_figure_params(dpi=80, facecolor='white')
    #sc.pl.highest_expr_genes(adata, n_top=25, )

    # get embeddings ----
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.tsne(adata)
    sc.tl.umap(adata, n_components = 2)

def parameters(__desc__):
    """Perse arguments
    Args:
        __desc__ (str): usage

    Returns:
        args (argparse.Namespace): retern the input arguments  
    """
    parser = argparse.ArgumentParser(
        description = __desc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter #show default
    )
    parser.add_argument(
        '-s',
        dest='salomn_dir',
        help='input salmon directory path',
        type=str,
        required=True
    )
    parser.add_argument(
        '-g',
        dest='gtf_file',
        help='gtf_file',
        type=str,
        required=True
    )
    parser.add_argument(
        '-f',
        dest='feature_file',
        help='input feature file',
        type=str,
        required=True
    )
    parser.add_argument(
        '-t',
        dest='tx2gene_file',
        help='input the tx2gene file',
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
    """display the input argument

    Args:
        args (argparse.Namespace): the input argument
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
    Norm_file = args.BC_Norm_file
    samples = args.target_tissues
    prefix = args.prefix
    outdir = args.out_dir

    os.makedirs(outdir, exist_ok=True)
        
    # run process ---
   

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