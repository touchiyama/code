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