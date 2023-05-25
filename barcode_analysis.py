import os
import re
import sys
import glob
import pandas as pd
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt

# intialization
file = 'A:/PROJECT/PG4796/300_INFO/00_data/PG4796_BC_info.xlsx'
bc_df = pd.read_excel(file)
BC_ID = pd.Series(['BC_'+ str(i) for i in bc_df.loc[:, 'BC No.']])
cols =  ['BC_ID', 'Sequence']
cnt_index = pd.concat([BC_ID, bc_df.loc[:, 'target_seq']], axis=1)
cnt_index.set_axis(cols, axis='columns', inplace=True)

# join count value data (cDNA)
insert_cnt = 'A:/PROJECT/PG4796/300_INFO/03_BCcnt/cDNA/*.txt'
insert_cnts =  sorted(glob.glob(insert_cnt))

for i, txt in enumerate(insert_cnts):
    sample_id = os.path.splitext(os.path.basename(txt))[0]
    df = pd.read_csv(txt, header=None, delim_whitespace=True)
    df.columns = [f'{sample_id}_Count', 'Sequence']
    df_mod = df.loc[:, ['Sequence', f'{sample_id}_Count']]
    if i == 0:
        df_new = pd.merge(cnt_index, df_mod, how='left', on='Sequence')
    else:
        df_new = pd.merge(df_new, df_mod, how='left', on='Sequence')
cdna_cnt = df_new.fillna(0)

# cDNA_freqの算出
for i in cdna_cnt.columns[2:].to_list():
    name = re.compile(r'(.*)_Count').search(i).group(1)
    freq = name + "_Freq"
    cdna_cnt[freq] = cdna_cnt[i] / cdna_cnt[i].sum()

# qPCR
file = 'A:/PROJECT/PG4796/300_INFO/00_data/PG4796_normalizationGbeta.xlsx'
qPCR_df = pd.read_excel(file, sheet_name='Gbeta')
qPCR = {}
for i in range(len(qPCR_df)):
    Ti = qPCR_df.loc[i, 'Tissue']
    ID = qPCR_df.loc[i, 'Animal No.']
    name = 'cDNA#' + Ti + '#' + str(ID)
    qPCR[name] = qPCR_df.loc[i, 'Mean']

# FCとFDRを算出 (cDNA) qPCRの結果に合わせ、cortexを含めない
cdna_c_list = [
    #'Library#Input_1#101_Count',
    #'Library#Input_2#102_Count',
    'cDNA#Brain#101_Count',
    'cDNA#Brain#102_Count',
    #'cDNA#Brainstem.Pons#101_Count',
    #'cDNA#Brainstem.Pons#102_Count',
    #'cDNA#Cerebellum#101_Count',
    #'cDNA#Cerebellum#102_Count',
    #'cDNA#Cortex#101_Count',
    #'cDNA#Cortex#102_Count',
    'cDNA#Heart#101_Count',
    'cDNA#Heart#102_Count',
    #'cDNA#Hippocampus#101_Count',
    #'cDNA#Hippocampus#102_Count',
    'cDNA#Liver#101_Count',
    'cDNA#Liver#102_Count',
    'cDNA#Lung#101_Count',
    'cDNA#Lung#102_Count',
    #'cDNA#Midbrain#101_Count',
    #'cDNA#Midbrain#102_Count',
    'cDNA#Muscle#101_Count',
    'cDNA#Muscle#102_Count',
    #'cDNA#Olfactory_bulb#101_Count',
    #'cDNA#Olfactory_bulb#102_Count',
    'cDNA#Spinal_cord#101_Count',
    'cDNA#Spinal_cord#102_Count',
    #'cDNA#Striatum#101_Count',
    #'cDNA#Striatum#102_Count',
    #'cDNA#Thalamus.Hypothalamus#101_Count',
    #'cDNA#Thalamus.Hypothalamus#102_Count',
]

# FCとFDRを算出 (cDNA)
before_cnt_1 = 'Library#Input_1#101_Count'
before_cnt_2 = 'Library#Input_2#102_Count'

before_freq_1 = 'Library#Input_1#101_Freq'
before_freq_2 = 'Library#Input_2#102_Freq'


for name in cdna_c_list:
    a = re.compile(r'(.*)_Count').search(name).group(1)
    name_freq = a + '_Freq'
    if '101' in name:
        print(f'{name_freq} vs {before_freq_1}\n')
        FC_name = a + '_FC'
        after_freq = cdna_cnt[name_freq]
        before_freq = cdna_cnt[before_freq_1]
        cdna_cnt[FC_name] = after_freq / before_freq

        print(f'{name} vs {before_cnt_1}\n')
        FDR_name = a + '_FDR'
        after_cnt = cdna_cnt[name]
        before_cnt = cdna_cnt[before_cnt_1]

        print(f'{name} normalization\n')
        Norm_name = a + '_Norm'
        qPCR_value = qPCR[a]
        cdna_cnt[Norm_name] = (after_freq / before_freq) * qPCR_value

    elif '102' in name:
        print(f'{name_freq} vs {before_freq_2}\n')
        FC_name =  a + '_FC'
        after_freq = cdna_cnt[name_freq]
        before_freq = cdna_cnt[before_freq_2]
        cdna_cnt[FC_name] = after_freq / before_freq

        print(f'{name} vs {before_cnt_2}\n')
        after_cnt = cdna_cnt[name]
        before_cnt = cdna_cnt[before_cnt_2]
        FDR_name = a + '_FDR'

        print(f'{name} normalization\n')
        Norm_name = a + '_Norm'
        qPCR_value = qPCR[a]
        cdna_cnt[Norm_name] = (after_freq / before_freq) * qPCR_value

cDNA_norm_list = [
    'BC_ID', 'Sequence',
    'cDNA#Brain#101_Norm',
    'cDNA#Brain#102_Norm',
    'cDNA#Spinal_cord#101_Norm',
    'cDNA#Spinal_cord#102_Norm',
    'cDNA#Heart#101_Norm',
    'cDNA#Heart#102_Norm',
    'cDNA#Liver#101_Norm',
    'cDNA#Liver#102_Norm',
    'cDNA#Lung#101_Norm',
    'cDNA#Lung#102_Norm',
    'cDNA#Muscle#101_Norm',
    'cDNA#Muscle#102_Norm',
]
cDNA_norm = cdna_cnt.loc[:, cDNA_norm_list]

outfile = 'A:/PROJECT/PG4796/900_NOUHIN/report_v03_230331/PG4796_cDNA_Barcode_Summary.xlsx'
with pd.ExcelWriter(outfile) as writer:
    cDNA_norm.to_excel(writer, sheet_name='cDNA_norm', index=False)
    cdna_cnt.to_excel(writer, sheet_name='cDNA_stats', index=False)


file = 'A:/PROJECT/PG4796/300_INFO/00_data/PG4796_Ascan_combination.xlsx'
ascan_df = pd.read_excel(file)

df_tmp = ascan_df.groupby(['pRC', 'AA'])['pAAV'].agg(','.join).reset_index()
df_tmp.columns = ['Ascan', 'AA_seq', 'BC_ID']

df_dna_norm = cDNA_norm.drop(columns='Sequence')
df_info = ascan_df.loc[:, ['pRC', 'pAAV']].rename(columns={'pRC': 'Ascan', 'pAAV': 'BC_ID'})
df_norm = pd.merge(df_info, df_dna_norm, how='left', on='BC_ID')
df_mean = df_norm.groupby(['Ascan']).agg('mean').reset_index()

df_dna_merge = pd.merge(df_tmp, df_mean, how='left', on='Ascan')

dna_c_list = [
    'cDNA#Brain',
    'cDNA#Spinal_cord',
    'cDNA#Heart',
    'cDNA#Liver',
    'cDNA#Lung',
    'cDNA#Muscle',
]

for i, name in enumerate(dna_c_list):
    rep1 = df_dna_merge.loc[:, f'{name}#101_Norm']
    rep2 = df_dna_merge.loc[:, f'{name}#102_Norm']
    ave = (rep1 + rep2) / 2
    total = ((((rep1 - ave) **2) + ((rep2 - ave) **2))) / (2 - 1)
    std = np.sqrt(total) / np.sqrt(2)
    df = pd.concat([df_dna_merge.loc[:, 'Ascan'], ave, std], axis=1)
    df.columns = ['Ascan', f'{name}_Mean', f'{name}_Se']
    if i == 0:
        df_new = df
    else:
        df_new = pd.merge(df_new, df, how='left', on='Ascan')

df_dna_final = pd.merge(df_dna_merge, df_new, how='left', on='Ascan')

outfile = 'A:/PROJECT/PG4796/900_NOUHIN/report_v03_230331/PG4796_cDNA_Ascan_Summary.xlsx'
with pd.ExcelWriter(outfile) as writer:
    df_dna_final.to_excel(writer, sheet_name='cDNA', index=False)


# Tab/Vab ---
file = 'A:/PROJECT/PG4796/900_NOUHIN/report_v03_230331/PG4796_cDNA_Ascan_Summary.xlsx'
df_dna = pd.read_excel(file, sheet_name='cDNA')

dna_1 = [
    'cDNA#Brain#101_Norm', 'cDNA#Spinal_cord#101_Norm', 'cDNA#Heart#101_Norm', 
    'cDNA#Liver#101_Norm', 'cDNA#Lung#101_Norm', 'cDNA#Muscle#101_Norm'
]
dna_2 = [
    'cDNA#Brain#102_Norm', 'cDNA#Spinal_cord#102_Norm', 'cDNA#Heart#102_Norm', 
    'cDNA#Liver#102_Norm', 'cDNA#Lung#102_Norm', 'cDNA#Muscle#102_Norm'
]
df_dna_1 = df_dna.loc[:, ['Ascan', 'AA_seq'] + dna_1]
df_dna_2 = df_dna.loc[:, ['Ascan', 'AA_seq'] + dna_2]

dna_12 = [
    'Ascan', 'AA_seq', 'cDNA#Brain', 'cDNA#Spinal_cord', 'cDNA#Heart',
    'cDNA#Liver', 'cDNA#Lung', 'cDNA#Muscle'
]
df_dna_1.columns = dna_12
df_dna_2.columns = dna_12

# Tab ---
df_Tab = pd.DataFrame()
for i in range(len(df_dna)):
    rep1_Tab = df_dna_1.iloc[i, 2:] / df_dna_1.iloc[i, 2:].sum()
    rep2_Tab = df_dna_2.iloc[i, 2:] / df_dna_2.iloc[i, 2:].sum()
    ave = (rep1_Tab + rep2_Tab) / 2
    total = ((((rep1_Tab - ave) **2) + ((rep2_Tab - ave) **2))) / 2
    std = np.sqrt(total.values.tolist())
    info = df_dna_1.iloc[i, :2].tolist()
    col = rep1_Tab.tolist() + rep2_Tab.tolist() + ave.tolist() + std.tolist()
    df_Tab = df_Tab.append(pd.DataFrame(info + col).T, ignore_index=True)

dna_12_mean = [
    'cDNA#Brain_Tab_Mean', 'cDNA#Spinal_cord_Tab_Mean', 'cDNA#Heart_Tab_Mean',
    'cDNA#Liver_Tab_Mean', 'cDNA#Lung_Tab_Mean', 'cDNA#Muscle_Tab_Mean'
]
dna_12_sd = [
    'cDNA#Brain_Tab_Sd', 'cDNA#Spinal_cord_Tab_Sd', 'cDNA#Heart_Tab_Sd',
    'cDNA#Liver_Tab_Sd', 'cDNA#Lung_Tab_Sd', 'cDNA#Muscle_Tab_Sd'
]
df_Tab.columns =  ['Ascan', 'AA_seq'] + dna_1 + dna_2 + dna_12_mean + dna_12_sd

df_Tab_12 = ['Ascan', 'AA_seq'] + dna_12_mean + dna_12_sd

# 組織で順位付け ---
for i in range(len(df_Tab)):
    seqID = df_Tab.loc[i, 'Ascan']
    AA_seq = df_Tab.loc[i, 'AA_seq']
    x = [name.replace('_Tab_Mean', '').replace('cDNA#', '') for name in dna_12_mean]
    y = df_Tab.loc[i, dna_12_mean].tolist()
    z = df_Tab.loc[i, dna_12_sd].tolist()

    fig, ax = plt.subplots() # figsize=(8, 3)
    error_bar_set = dict(lw=1, capthick=0.5, capsize=2)
    ax.bar(x, y, color='gray', edgecolor='gray', yerr=z , error_kw=error_bar_set, align='center')
    #ax.axvspan(x[0],  x[10], color='yellow', alpha=0.1)
    plt.ylim(0, 1.0)
    plt.xticks(rotation=90)
    #plt.axhline(y=1, color='red', linestyle='--')
    plt.title(f'{seqID} ({AA_seq})')
    plt.ylabel('Tab_mean')
    fname = 'PG4796_' + seqID + '_Tab_mean.png'
    plt.savefig('A:/PROJECT/PG4796/900_NOUHIN/report_v03_230331/Tab_cDNA/' + fname, bbox_inches='tight', dpi=300)
    plt.show()


# Vab ---
info_list = []
for name in ['_101_Vab', '_102_Vab', '_Vab_Mean', '_Vab_Sd']:
    for ID in df_dna['Ascan']:
        rename = ID + name
        info_list.append(rename)
info = pd.DataFrame(info_list, columns=['Name'])

for i in range(2, len(df_dna_1.columns)):
    rep1_Vab = df_dna_1.iloc[:, i] / df_dna_1.iloc[:, i].sum()
    rep2_Vab = df_dna_2.iloc[:, i] / df_dna_2.iloc[:, i].sum()
    ave = (rep1_Vab + rep2_Vab) / 2
    total = ((((rep1_Vab - ave) **2) + ((rep2_Vab - ave) **2))) / 2
    std = np.sqrt(total.values.tolist())
    col = rep1_Vab.tolist() + rep2_Vab.tolist() + ave.tolist() + std.tolist()
    df_tissue = pd.DataFrame(col, columns=[df_dna_1.columns[i]])
    df_tmp = pd.concat([info, df_tissue], axis=1)
    if i == 2:
        df_Vab = df_tmp
    else:
        df_Vab = pd.merge(df_Vab, df_tmp, how='left', on='Name')

mean_list = []
for ID in df_dna['Ascan']:
    rename = ID + '_Vab_Mean'
    mean_list.append(rename)

sd_list = []
for ID in df_dna['Ascan']:
    rename = ID + '_Vab_Sd'
    sd_list.append(rename)

mean_idx = [df_Vab[df_Vab['Name'] == name].index[0] for name in mean_list]
sd_idx = [df_Vab[df_Vab['Name'] == name].index[0] for name in sd_list]
df_Vab_12 = df_Vab.T

all_tissues = [
    'cDNA#Brain', 'cDNA#Spinal_cord', 'cDNA#Heart',
    'cDNA#Liver', 'cDNA#Lung', 'cDNA#Muscle'
]

for tissue in all_tissues:
    df_Vab_mean = df_Vab_12.loc[tissue, mean_idx].reset_index(drop=True)
    df_Vab_sd = df_Vab_12.loc[tissue, sd_idx].reset_index(drop=True)
    df_stats = pd.concat([df_dna['Ascan'], df_Vab_mean, df_Vab_sd], axis=1)
    df_stats.columns = ['Ascan', f'{tissue}_mean', f'{tissue}_sd']
    df_sort = df_stats.sort_values(f'{tissue}_mean', ascending=False).reset_index()
    x = df_sort['Ascan'].tolist()
    y = df_sort[f'{tissue}_mean'].values.tolist()
    z = df_sort[f'{tissue}_sd'].values.tolist()

    fig, ax = plt.subplots(figsize=(12, 4))
    error_bar_set = dict(lw=1, capthick=0.5, capsize=2)
    bar = ax.bar(x, y, color='white', edgecolor='black', yerr=z , error_kw=error_bar_set)

    idx = x.index('CereAAV')
    bar[idx].set_color('red')

    idx = x.index('mi-342')
    bar[idx].set_color('blue')

    plt.xticks(rotation=90)
    #plt.axhline(y=1, color='red', linestyle='--')
    #ax.legend(
    #    bbox_to_anchor=(1.05, 1),
    #    loc='upper right',
    #    #borderaxespad=0,
    #    fontsize=8
    #)
    #plt.ylim(0, 1)
    title = tissue.replace('cDNA#', '')
    plt.title(f'{title}')
    plt.ylabel('Vab_mean')
    fname = 'PG4796_' + title + '_Vab_mean.png'
    plt.savefig('A:/PROJECT/PG4796/900_NOUHIN/report_v03_230331/Vab_cDNA/' + fname, bbox_inches='tight', dpi=300)
    plt.show()

# output ---
outfile = 'A:/PROJECT/PG4796/900_NOUHIN/report_v03_230331/PG4796_cDNA_Ascan_Vab_Tab.xlsx'
with pd.ExcelWriter(outfile) as writer:
    df_Vab.to_excel(writer, sheet_name='Vab', index=False)
    df_Tab.to_excel(writer, sheet_name='Tab', index=False)
