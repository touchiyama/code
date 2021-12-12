# %%
import os
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# %%
file = '/Users/tomoyauchiyama/code/Consensus_frequency.xlsx'
df = pd.read_excel(file, header=None)
data = df.iloc[2:, 15:int(df.columns.size)+1]
data.columns = df.iloc[1, 3:15]

print('# max ---------')
print(data.max())
print('# min ---------')
print(data.min())

# %%
# PG4373_[13/14/15/16/17/18]_freq
groupA = df.iloc[2:, 15:21]
groupA_sum = groupA.sum(axis=1)
groupA_con = pd.DataFrame([df.iloc[2:, 0], groupA_sum]).transpose()
groupA_con.columns = ['seqID', 'freq']
groupA_sort = groupA_con.sort_values(by='freq')

# PG4373_[19/20/21/22/23/24]_freq
groupB = df.iloc[2:, 21:int(df.columns.size)+1]
groupB_sum = groupB.sum(axis=1)
groupB_con = pd.DataFrame([df.iloc[2:, 0], groupB_sum]).transpose()
groupB_con.columns = ['seqID', 'freq']
groupB_sort = groupB_con.sort_values(by='freq')

# %%
fig = make_subplots()

fig.add_trace(
    go.Scatter(
        x=groupA_sort['seqID'],
        y=groupA_sort['freq'],
        mode='markers',
        marker=dict(
            color=None,
            size=2,
            line=dict(color='black', width=1)
        )
    )
)
"""
fig.add_trace(
    go.Scatter(
        x=groupB_sort['seqID'],
        y=groupB_sort['freq'],
        mode='markers',
        marker=dict(
            color=None,
            size=2,
            line=dict(color='black', width=1)
        )
    )
)
"""
fig.add_hline(
    y=groupA_sort['freq'].mean(),
    line_width=1.0,
    line_dash='dash',
    line_color='red'
)
fig.update_layout(
    plot_bgcolor='white'
    #height=800,
    #width=900
)
fig.update_xaxes(
    title='SeqID',
    showline=True,
    linewidth=1,
    linecolor='black',
    mirror=True,
    ticks='inside'
    #range=(min_x, max_x),
    #dtick=5
)
fig.update_yaxes(
    title='Freq',
    showline=True,
    linewidth=1,
    linecolor='black',
    mirror=True,
    ticks='inside'
    #range=(min_y, max_y),
    #dtick=5
)

fig.show()

outdir = '/Users/tomoyauchiyama/code'
htmlfile = os.path.join(
    outdir,
    'groupA_seq_freq.html'
)
fig.write_html(htmlfile)

# %%
print(df.index.size*0.05)
# %%
# 基本統計量 ----------
print('# groupA ----------')
print(groupA_sort['freq'].describe())
print('# groupB ----------')
print(groupB_sort['freq'].describe())

# %%
import matplotlib.pyplot as plt

plt.boxplot([groupA_sort['freq'], groupB_sort['freq']], labels=['groupA', 'groupB'])
plt.ylim([0, 0.05])
plt.show()

# %%
print(groupA_sort)
# %%
print(groupB_sort)
# %%
seq_data = df.iloc[2:, [0, 1]]
seq_data.columns = df.iloc[1, [0, 1]]
seq_data = seq_data.reset_index(drop=True)
print(seq_data.head())

# %%
groupA_maxID = groupA_sort[groupA_sort['freq'] == groupA_sort['freq'].max()]['seqID'].values[0]
groupA_seq = seq_data[seq_data['Header_ID'] == groupA_maxID]['Sequence'].values[0]
groupB_maxID = groupB_sort[groupB_sort['freq'] == groupB_sort['freq'].max()]['seqID'].values[0]
groupB_seq = seq_data[seq_data['Header_ID'] == groupB_maxID]['Sequence'].values[0]

# make group_A db -----------
fasta_file = os.path.join(
    os.path.abspath('.'),
    'groupA_maxFreqSeq_db.fa'
)

with open(fasta_file, 'w') as wf:
    wf.write(f'>{groupA_maxID}\n')
    wf.write(f'{groupA_seq}\n')


# make group_B db -----------
fasta_file = os.path.join(
    os.path.abspath('.'),
    'groupB_maxFreqSeq_db.fa'
)

with open(fasta_file, 'w') as wf:
    wf.write(f'>{groupB_maxID}\n')
    wf.write(f'{groupB_seq}\n')

# %%
# make group_AB quary -----------
fasta_file = os.path.join(
    os.path.abspath('.'),
    'groupAB_FreqSeq_quary.fa'
)

with open(fasta_file, 'w') as wf:
    for i in range(len(seq_data)):
        SeqID = seq_data.loc[i, 'Header_ID']
        Seq = seq_data.loc[i, 'Sequence']
        wf.write(f'>{SeqID}\n')
        wf.write(f'{Seq}\n')

# %%
# data processing

infile = '/Users/tomoyauchiyama/code/Consensus_frequency/result_groupA_2.txt'
outfile = '/Users/tomoyauchiyama/code/Consensus_frequency/result_groupA_maxscore.txt'
nrA = {}
groupA_score = {}
with open(outfile, 'w') as wf:
    with open(infile, 'r') as rf:
        for line in rf.readlines():
            line = line.rstrip('\n')
            tmp = line.split()
            ID = tmp[0] + ',' + tmp[1]
            if nrA.get(ID):
                pass
            else:
                nrA[ID] = '-'
                groupA_score[tmp[0]] = tmp[-1]
                wf.write(f'{line}\n')

infile = '/Users/tomoyauchiyama/code/Consensus_frequency/result_groupB_2.txt'
outfile = '/Users/tomoyauchiyama/code/Consensus_frequency/result_groupB_maxscore.txt'
nrB = {}
groupB_score = {}
with open(outfile, 'w') as wf:
    with open(infile, 'r') as rf:
        for line in rf.readlines():
            line = line.rstrip('\n')
            tmp = line.split()
            ID = tmp[0] + ',' + tmp[1]
            if nrB.get(ID):
                pass
            else:
                nrB[ID] = '-'
                groupB_score[tmp[0]] = tmp[-1]
                wf.write(f'{line}\n')

# %%
scoreA = pd.DataFrame(groupA_score.values(), index=groupA_score.keys()).reset_index(drop=True)
scoreA = scoreA.astype(float)
scoreA.columns = ['Score']
groupA_con = groupA_con.reset_index(drop=True)

# %%
scoreB = pd.DataFrame(groupB_score.values(), index=groupB_score.keys()).reset_index(drop=True)
scoreB = scoreB.astype(float)
scoreB.columns = ['Score']
groupB_con = groupB_con.reset_index(drop=True)

# %%
from sklearn.linear_model import LinearRegression

def calc_lm(x, y):
    return LinearRegression().fit(x.reshape(-1, 1), y)

# %%

fig = make_subplots()

fig.add_trace(
    go.Scatter(
        x=scoreA['Score'] / scoreA['Score'].max(),
        y=groupA_con['freq'] / groupA_con['freq'].max(),
        mode='markers',
        marker=dict(
            color=None,
            size=2,
            line=dict(color='black', width=1)
        )
    )
)
"""
fig.add_trace(
    go.Scatter(
        x=scoreB['Score'] / scoreB['Score'].max(),
        y=groupB_con['freq'] / groupB_con['freq'].max(),
        mode='markers',
        marker=dict(
            color=None,
            size=2,
            line=dict(color='black', width=1)
        )
    )
)
"""
"""
fig.add_hline(
    y=groupA_sort['freq'].mean(),
    line_width=1.0,
    line_dash='dash',
    line_color='red'
)
"""
x = scoreA['Score'].to_numpy() / scoreA['Score'].max()
y = groupA_con['freq'].to_numpy() / groupA_con['freq'].max()
fig.add_trace(
    go.Scatter(
        x=x,
        y=calc_lm(x, y).predict(x.reshape(-1, 1)),
        line_width=1.0,
        line=dict(dash='dot', color='gray')
    )
)
fig.update_layout(
    plot_bgcolor='white',
    height=800,
    width=900
)
fig.update_xaxes(
    title='Similarity Score',
    showline=True,
    linewidth=1,
    linecolor='black',
    mirror=True,
    ticks='inside'
    #range=(0, scoreA['Score'].max()),
    #dtick=5
)
fig.update_yaxes(
    title='Freq Score',
    showline=True,
    linewidth=1,
    linecolor='black',
    mirror=True,
    ticks='inside'
    #range=(groupA_con['freq'].min(), groupA_con['freq'].max()),
    #dtick=5
)

fig.show()


outdir = '/Users/tomoyauchiyama/code/Consensus_frequency'
htmlfile = os.path.join(
    outdir,
    'groupA_Freq_Simiarity_Score.html'
)
fig.write_html(htmlfile)

# %%
# groupA correlation --------------------------------------------
s1 = scoreA['Score'] / scoreA['Score'].max()
s2 = groupA_con['freq'].astype(float) / groupA_con['freq'].max()

res = s1.corr(s2, method='pearson')
#methed: 'pearson', 'spearman', 'kendall'
print(f'groupA correlation: {res}')

# %%
# groupB correlation --------------------------------------------
s1 = scoreB['Score'] / scoreB['Score'].max()
s2 = groupB_con['freq'].astype(float) / groupB_con['freq'].max()

res = s1.corr(s2, method='pearson')
#methed: 'pearson', 'spearman', 'kendall'
print(f'groupB correlation: {res}')

# %%
"""[summary]
correlationの計算に関しては、
各グループで、Freq scoreとSimilarity Scoreの分布を見て判断する

"""


# %%
# fastq feature -----------------------------------
import re
import gzip

file = '/Users/tomoyauchiyama/Downloads/K4BR8_PG4536_01B2003_H1_L001_R1.fastq.gz'
outf = './fw.fa'

flag = 0
with open(outf, 'w') as wf:
    with gzip.open(file, 'rt') as fq:
        for line in fq.readlines():
            line = line.rstrip('\n')
            if re.compile(r'^@').search(line):
                flag = 1
            else:
                if flag == 1:
                    forward = line[:46]
                    wf.write(f'{forward}\n')
                    flag = 0

# %%
file = './fw.fa'
forward = 'GCAGTATGGTTCTGTATCTACCAACCTCCAGAGAGGCCAGAGAGGC'
outf = './mis_fw.fa'

cnt = 0
with open(outf, 'w') as wf:
    with open(file, 'r') as rf:
        for line in rf.readlines():
            line = line.rstrip('\n')
            if line != forward:
                cnt += 1
                id = '>seq_' + str(cnt)
                wf.write(f'{id}\n')
                wf.write(f'{line}\n')
# %%
file = '/Users/tomoyauchiyama/Downloads/K4BR8_PG4536_01B2003_H1_L001_R1.fastq.gz'
outf = './01_R1_fastq_feat/rv.fa'

flag = 0
with open(outf, 'w') as wf:
    with gzip.open(file, 'rt') as fq:
        for line in fq.readlines():
            line = line.rstrip('\n')
            if re.compile(r'^@').search(line):
                flag = 1
            else:
                if flag == 1:
                    reverse = line[68:114]
                    wf.write(f'{reverse}\n')
                    flag = 0

# %%
file = './01_R1_fastq_feat/rv.fa'
reverse = 'CCCAGGCGGCCACCGCAGATGTCAACACACAAGGCGTTCTTCCAGG'
outf = './01_R1_fastq_feat/mis_rv.fa'

cnt = 0
with open(outf, 'w') as wf:
    with open(file, 'r') as rf:
        for line in rf.readlines():
            line = line.rstrip('\n')
            if line != reverse:
                cnt += 1
                id = '>seq_' + str(cnt)
                wf.write(f'{id}\n')
                wf.write(f'{line}\n')

# %%

file = '/Users/tomoyauchiyama/code/Consensus_frequency/01_R1_fastq_feat/K4BR8_PG4536_01B2003_H1_L001_R1.fastq.gz'
outf = './01_R1_fastq_feat/insert_seq.fa'

cnt = 0
with open(outf, 'w') as wf:
   with gzip.open(file, 'rt') as fq:
        for line in fq.readlines():
            line = line.rstrip('\n')
            if re.compile(r'^@').search(line):
                flag = 1
                cnt += 1
            else:
                if flag == 1:
                    id = '>seq_' + str(cnt)
                    wf.write(f'{id}\n')
                    seq = line[46:67]
                    wf.write(f'{seq}\n')
                    flag = 0

# %%
# meme parse ----------------------------------
import pandas as pd

meme_res = '/Users/tomoyauchiyama/code/Consensus_frequency/01_R1_fastq_feat/insert_freq/meme.txt'

#ALPHABET= ACGT
#letter-probability matrix: alength= 4 w= 21 nsites= 61416 E= 5.5e-672
flag = 0
cnt = 0
nucl_freq = pd.DataFrame()
with open(meme_res, 'r') as meme:
    for line in meme.readlines():
        line = line.rstrip('\n')
        if 'letter-probability' in line:
            flag = 1
            tmp = line.split()
            N = int(tmp[5])
        else:
            if flag == 1:
                cnt += 1
                if cnt == N:
                    flag = 0
                tmp = line.split()
                freq = pd.Series([float(i) for i in tmp])
                nucl_freq = nucl_freq.append(freq, ignore_index=True)

nucl_freq.columns = ['A', 'C', 'G', 'T']
nucl_freq['pos'] = nucl_freq.index + 1

# %%
print(nucl_freq)
# %%
print(nucl_freq.plot.bar())

# %%
print(type(nucl_freq.loc[i, nucl_freq.columns[4]]))

# %%
# nucl frequency visialization ---------------------
import plotly.graph_objects as go
from plotly.subplots import make_subplots

fig = make_subplots()
colors = ['red', 'blue', '#f1c40f', 'green'] # ['A', 'C', 'G', 'T']

for j in range(4):
    fig.add_trace(
        go.Bar(
            x=nucl_freq[nucl_freq.columns[4]],
            y=nucl_freq[nucl_freq.columns[j]],
            name=nucl_freq.columns[j],
            marker=dict(color=colors[j])
        )
    )

fig.update_layout(
    plot_bgcolor='white'
    #height=800,
    #width=900
)
fig.update_xaxes(
    title='Position',
    showline=True,
    linewidth=1,
    linecolor='black',
    mirror=True,
    ticks='inside',
    #range=(min_x, max_x),
    dtick=1
)
fig.update_yaxes(
    title='Frequncy',
    showline=True,
    linewidth=1,
    linecolor='black',
    mirror=True,
    ticks='inside',
    range=(0, 1),
    dtick=0.1
)

fig.show()

outdir = '/Users/tomoyauchiyama/code/Consensus_frequency/01_R1_fastq_feat'
htmlfile = os.path.join(
    outdir,
    'insert_nucl_freq.html'
)
fig.write_html(htmlfile)

# %%
