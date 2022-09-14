# %%
from operator import index
import re
import os
import sys
import json
import pandas as pd
import numpy as np
import cmdstanpy as stan
import logomaker as lm
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
# %%
target_seq = '/Users/tomoyauchiyama/code/CNN/test/enrich.fa'
nucl = pd.DataFrame()

with open(target_seq, 'rt') as fq:
    for line in fq.readlines():
        line = line.rstrip('\n')
        if re.compile(r'^>').search(line):
                flag = 1
        else:
            if flag == 1:
                insert = pd.Series([str(i) for i in line])
                nucl = nucl.append(insert, ignore_index=True)
                flag = 0

# %%
col = [str(i) for i in range(1, nucl.columns.size + 1)]
nucl.columns = col

# %%
# 事前確率　--------------------------------------------------------
total = nucl.index.size
pos_freq = pd.DataFrame()
for pos in col:
    col_list = list(nucl[pos])
    A = col_list.count('A') / total
    T = col_list.count('T') / total
    G = col_list.count('G') / total
    C = col_list.count('C') / total
    #N = col_list.count('N') / total
    freq = pd.Series([A, T, G, C])
    pos_freq = pos_freq.append(freq, ignore_index=True)

# %%
pos_freq.columns = ['A', 'T', 'G', 'C']
pos_freq.index += 1
print(pos_freq)

# %%
pos_freq = pd.read_excel('/Users/tomoyauchiyama/code/CNN/test/PR2688_C_60nt.xlsx')
pos_freq.index += 1
print(pos_freq)

# %%
# logomekerでシーケンスロゴを作成-------------------------------------
# https://logomaker.readthedocs.io/en/latest/examples.html
# pip install logomaker

color_scheme = {
    'T' : [0, 0.5, 0], # green
    'A' : [1, 0, 0], # red
    'G' : [1, 0.65, 0], # yellow
    'C' :[0, 0, 1], # blue
    'N': 'gray'
}
logo = lm.Logo(
    pos_freq,
    baseline_width=0.1,
    vpad=0.08,
    fade_probabilities=False,
    stack_order='small_on_top',
    color_scheme=color_scheme
    #font_name='Luxi Mono',
    #color_scheme='class',
    #font_name='Rosewood Std',
    #figsize=(30, 2.5)
)

logo.style_spines(spines=['left', 'right'], visible=False)

# style using Axes methods
logo.ax.set_xticks(np.arange(1, len(pos_freq)+1, step=1))
logo.ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
logo.ax.set_xlabel('Position')
logo.ax.set_ylabel('Probability')
plt.show()

# %%
# calc PWM score
def pwm(cnt, total, bg=0.25):
    p = (cnt + (np.sqrt(total) * 0.25)) / (total + (4 * (np.sqrt(total) * 0.25)))
    pwm_score = np.log2(p/bg)

    return pwm_score

# %%
col

# %%
pos_cnt = pd.DataFrame()
for pos in col:
    col_list = list(nucl[pos])
    A = col_list.count('A')
    T = col_list.count('T')
    G = col_list.count('G')
    C = col_list.count('C')
    freq = pd.Series([A, T, G, C])
    pos_cnt = pos_cnt.append(freq, ignore_index=True)

cnt = pos_cnt.to_numpy().astype(int)
print(cnt)

# %%
total = nucl.index.size
pwm_score = pwm(cnt, total)
print(pwm_score)

# %%
#
# modeling
#

# 事前確率　--------------------------------------------------------
total = nucl.index.size
pos_refreq = pd.DataFrame()
for pos in col:
    col_list = list(nucl[pos])
    A = col_list.count('A') / total
    T = col_list.count('T') / total
    G = col_list.count('G') / total
    C = col_list.count('C') / total
    #N = col_list.count('N') / total
    freq = pd.Series([A, T, G, C])
    pos_refreq = pos_refreq.append(freq, ignore_index=True)

# %%
"""
# smoothing for sampling convergence -----------------------------
# usage or example
# https://pypi.org/project/loess/#id11

#========== preparation ==========#
d_melt = pd.melt(
    pos_refreq.reset_index(),
    id_vars='index',
    ignore_index=True
)
d_melt.columns = ('i', 'j', 'Y')

print(d_melt)
"""

# %%
"""
#==========  smoothing  ==========#
from loess.loess_2d import loess_2d

x = d_melt['i'].to_numpy()
y = d_melt['j'].astype('int64').to_numpy()
z = d_melt['Y'].to_numpy()
#z = x + y

loess_res, _ = loess_2d(x, y, z, frac=0.2)

i, j = cnt.shape
print(loess_res.reshape(i, j))
"""
# %%
# json file ------------------------------------------------------
pre_freq = pos_refreq.transpose().to_numpy()
base, pos = pre_freq.shape

observed_data = dict(
    Ni=pos,
    Nt=base,
    Y=pre_freq.tolist()
)

json_name = 'PrefreqPositionMatrix' # have to change

json_file = os.path.join(
    '/Users/tomoyauchiyama/code',
    'STAN',
    json_name + '.json'
)

with open(json_file, 'w') as wf:
    json.dump(observed_data, wf, indent=2)

print(json_file)

# %%
# stan file ------------------------------------------------------
code = """data {
    int Ni;
    int Nt;
    real Y[Nt, Ni];
}
parameters {
    real r[Nt, Ni];
    real<lower=0> s_r;
    real<lower=0> s_Y;
}
model {
    for (i in 1:Nt)
        for (j in 3:Ni)
            target += normal_lpdf(
                r[i, j] | 2*r[i, j-1] - r[i, j-2],
                s_r
            );

    for (i in 1:Nt)
        for (j in 1:Ni)
            Y[i, j] ~ normal(r[i, j], s_Y);
}
generated quantities {
    matrix[Nt, Ni] Y_model;
    for (i in 1:Nt)
        for (j in 1:Ni)
            Y_model[i, j] = normal_rng(r[i, j], s_Y);
}
"""

stan_name = 'motif_model' # have to change

stan_file = os.path.join(
    '/Users/tomoyauchiyama/code',
    'STAN',
    stan_name + '.stan'
)

with open(stan_file, 'w') as wf:
    wf.writelines(code)

print(stan_file)
# %%
# compile stan file ----------------------------------------------
sm = stan.CmdStanModel(stan_file=stan_file)

# %%
print(sm.name)
print(sm.stan_file)
print(sm.exe_file)
print(sm.code())

# %%
"""
# initilazation --------------------------------------------------
init = dict(
    r=loess_res.reshape(i, j),
    s_r=1,
    s_Y=1
)

print(init)
"""

# %%
# MCMC sampling (NUTS method) ------------------------------------
fit = sm.sample(
    data=json_file,
    chains=3,
    iter_sampling=2500,
    iter_warmup=500,
    thin=5,
    seed=1234
    # inits=init
)

# %%
# MCMC sampling summary ------------------------------------------

# summaries of
# the total joint log-probability density lp__
# plus all model parameters
# and quantities of interest in a pandas.DataFrame

# type -> pandas (basic statistic)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
print(fit.summary().loc[:, ['Mean', 'StdDev', '5%', '50%', '95%', 'N_Eff', 'R_hat']])

# %%
# MCMC sampling analyzation --------------------------------------

# analyze the per-draw sampler parameters across all chains
# looking for potential problems which indicate that
# the sample isn’t a representative sample from the posterior
print(fit.diagnose())

# %%
# parameter data extraction --------------------------------------
# extract parmetar 'Y_model'
parm_Y = fit.stan_variable('Y_model')

# %%
print(parm_Y.shape)
# -> output: (1500, 4, 21)

# %%
s_r = fit.stan_variable('s_r')
print(s_r.mean())

# %%
# calc ave of parm_r -------------------
Y_mean = np.median(parm_Y, axis=0)

# realize array shape ------------------
print(Y_mean.shape)
# output: (4, 21)

# %%
# constract DataFrame ------------------
# -> 1024 base x (10 data x [p10, p25, p50, p75, p90])
probs = (10, 25, 50, 75, 90)
column = pd.Series([f'p{p}' for p in probs])
prob_Y = np.percentile(parm_Y, probs, axis=0).transpose()

est_d = pd.DataFrame()
for i in range(len(prob_Y)):
    val = i + 1
    for id, prob_y in enumerate(prob_Y[i], 1):
        df1 = pd.Series([id, val])
        df2 = pd.Series(prob_y)
        col = pd.concat([df1, df2])
        est_d = est_d.append(col, ignore_index=True)

est_d.columns = pd.concat([pd.Series(['ID', 'pos']), column])
est_d = est_d.sort_values(['ID', 'pos']).reset_index()

pd.set_option('display.max_rows', None)
print(est_d)

# %%
# prepare for data visualization ---------------------------------
p10 = {}
p25 = {}
p50 = {}
p75 = {}
p90 = {}

for i in range(len(est_d)):
    ID = int(est_d.loc[i, 'ID'])
    p10_val = str(est_d.loc[i, 'p10'])
    p25_val = str(est_d.loc[i, 'p25'])
    p50_val = str(est_d.loc[i, 'p50'])
    p75_val = str(est_d.loc[i, 'p75'])
    p90_val = str(est_d.loc[i, 'p90'])
    if p10.get(ID):
        p10[ID] += ',' + p10_val
        p25[ID] += ',' + p25_val
        p50[ID] += ',' + p50_val
        p75[ID] += ',' + p75_val
        p90[ID] += ',' + p90_val
    else:
        p10[ID] = p10_val
        p25[ID] = p25_val
        p50[ID] = p50_val
        p75[ID] = p75_val
        p90[ID] = p90_val

pred_df = pd.DataFrame()
for ID in sorted(est_d['ID'].unique()):
    xx = np.linspace(1, 1025, 1024).astype(int).tolist()
    p10y = [float(i) for i in p10[ID].split(',')]
    p25y = [float(i) for i in p25[ID].split(',')]
    p50y = [float(i) for i in p50[ID].split(',')]
    p75y = [float(i) for i in p75[ID].split(',')]
    p90y = [float(i) for i in p90[ID].split(',')]
    col = pd.Series([ID, xx, p10y, p25y, p50y, p75y, p90y])
    pred_df = pred_df.append(col, ignore_index=True)

pred_df.columns = ['ID', 'pos',
                   'p10', 'p25', 'p50', 'p75', 'p90']

print(pred_df)

# %%
import seaborn as sns

def get_colorpalette(n_colors):
    """legendとRBG値を対応させる関数

    Returns:
        rgb (list): RBG値が格納された配列
    """
    palette = sns.hls_palette(n_colors+1, l=0.6, s=1)
    #palette = sns.color_palette('plasma', n_colors+1)
    rgb = ['rgb({},{},{})'.format(*[x*256 for x in rgb]) for rgb in palette]
    return rgb

# %%
# 各位置におけるATGCの出現頻度(テストデータ)の可視化 ----------------

titles = ['A', 'T', 'G', 'C']
probs = (10, 25, 50, 75, 90)
pxx = [f'p{p}' for p in probs]

n = len(titles)
n_color = n
colors = get_colorpalette(n_color)

n_col = int(n/2)
if n_col > 5:
    n_col = 5
    n_row = int(n/n_col) + 1
else:
    if (n%2) == 0:
        n_row = int(n/n_col)
    else:
        n_row = int(n/n_col) + 1


# %%
outdir = '/Users/tomoyauchiyama/code/CNN/test'
outfile = 'motif_model_with_stan.html'

y_title = 'Frequency'
x_title = 'Position'
main_title = 'Freq of ATGC per postion (*σ(r) = ' + str(round(s_r.mean(), 2)) + ')'

fig = make_subplots(
    rows=n_row,
    cols=n_col,
    subplot_titles=titles,
    y_title=y_title,
    x_title=x_title,
    horizontal_spacing=0.04, # will change
    vertical_spacing=0.07    # will change
)

ii = -1
min_x = 1
max_x = pos
min_y = 0
max_y = 1
for i in range(1, n_row+1):
    for j in range(1, n_col+1):
        ii += 1
        y = pos_refreq.transpose().to_numpy()[ii, :]
        x = [i+1 for i in range(max_x)]
        pre_y = pred_df.loc[ii, 'p50']
        pre_x = x

        # fill area between trace0 and trace1
        for k in range(2):
            fig.add_trace(
                go.Scatter(
                    x=pre_x,
                    y=pred_df.loc[ii, pxx[k]],
                    fill=None,
                    mode='lines',
                    line=dict(color='black', width=0.2),
                        opacity=0.2*(k+1)
                ),
                row=i,
                col=j
            )
            fig.add_trace(
                go.Scatter(
                    x=pre_x,
                    y=pred_df.loc[ii, pxx[-k-1]],
                    fill='tonexty',
                    mode='lines',
                    line=dict(color='black', width=0.2),
                    opacity=0.2*(k+1)
                ),
                row=i,
                col=j
            )
        fig.add_trace(
            go.Scatter(
                x=pre_x,
                y=pre_y,
                mode='lines',
                line=dict(
                    color=colors[ii],
                    width=2.0
                )
            ),
            row=i,
            col=j
        )
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                mode='markers',
                marker=dict(
                    color=colors[ii],
                    size=5.0,
                    line=dict(
                        color='black',
                        width=2.0
                    )
                )
            ),
            row=i,
            col=j
        )
        """
        fig.add_trace(
            go.Bar(
                x=x,
                y=y,
                marker_color=colors[ii],
                opacity=1.0
            ),
            row=i,
            col=j
        )
        """
        fig.update_layout(
            plot_bgcolor='white'
            #height=800,
            #width=900
        )
        fig.update_xaxes(
            showline=True,
            linewidth=1,
            linecolor='black',
            mirror=True,
            ticks='inside',
            range=(min_x, max_x),
            row=i,
            col=j
        )
        fig.update_yaxes(
            showline=True,
            linewidth=1,
            linecolor='black',
            mirror=True,
            ticks='inside',
            range=(min_y, max_y),
            row=i,
            col=j
        )

fig.for_each_xaxis(
    lambda axis: axis.title.update(
        font=dict(
            color='black',
            size=15
        )
    )
)
fig.for_each_yaxis(
    lambda axis: axis.title.update(
        font=dict(
            color='black',
            size=15
        )
    )
)
fig.update_annotations(
    font=dict(size=12),
)
fig.update_layout(
    title=dict(
        text=main_title,
        x=0.5,
        xanchor='center'
    ),
    showlegend=False
)
fig.show()

htmlfile = os.path.join(outdir, outfile)
fig.write_html(htmlfile)


# %%
#========== visialize obserbed data ==========#
fig=plt.figure()
ax=fig.add_subplot(111)
x = np.arange(1, pos+1, 1)
legend =['A', 'T', 'G', 'C']
cols = ['red', 'green', 'orange', 'blue']

for i in range(4):
    plt.plot(
        x,
        pos_refreq.transpose().to_numpy()[i, :],
        linestyle='--',
        color=cols[i]
    )
    plt.scatter(
        x,
        pos_refreq.transpose().to_numpy()[i, :],
        label=legend[i],
        color=cols[i]
    )
    plt.plot(
        x,
        pred_df.loc[i, 'p50'],
        color=cols[i]
    )

ax.set_xlabel('Position')
ax.set_ylabel('Freqency')
ax.legend(loc="upper left")
ax.set_ylim(0, 1)
plt.legend(
    bbox_to_anchor=(1.05, 1),
    loc='upper left',
    borderaxespad=0,
    fontsize=10
)
#ax.legend(loc="upper right")
plt.xticks(x)
plt.show()

# %%
# 50%予測区間(平均値)で、各位置でATGCの出現頻度の総和が1になることを確認

tmp = []
for i in range(4):
    tmp.append(pred_df.loc[i, 'p50'])
arr = np.array(tmp)
np.sum(arr.transpose(), axis=1)

"""
array([0.990883 , 1.010306 , 1.009148 ,
       1.004706 , 1.009713 , 1.002827 ,
       1.004017 , 1.0040835, 1.006968 ,
       0.98814  , 0.9915585, 1.007901 ,
       0.9963495, 0.998081 , 0.988662 ,
       0.9996325, 0.9829535, 1.004559 ,
       1.0164875, 0.9981655, 1.0183   ])
"""

# %%
import numpy as np
import matplotlib.pyplot as plt

mu1, sigma1 = 100, 15
mu2, sigma2 = 70, 6
x1 = mu1 + sigma1 * np.random.randn(10000)
x2 = mu2 + sigma2 * np.random.randn(10000)

fig=plt.figure()
ax=fig.add_subplot(111)

plt.hist(x1, bins=50, density=True, rwidth=0.8, color='red', alpha=0.5)
plt.hist(x2, bins=50, density=True, rwidth=0.8, color='blue', alpha=0.5)
ax.set_title('sixth histogram $\mu1=100,\ \sigma1=15,\ \mu2=50,\ \sigma2=4$')
ax.set_xlabel('x')
ax.set_ylabel('freq')
ax.set_ylim(0,0.1)
plt.show()
# %%
# ライブラリー結果のモデリング --------------------------------
import pandas as pd

file = '/Users/tomoyauchiyama/code/PC0174_PC0131_LIB_summmary_211110.xlsx'
lib_df = pd.read_excel(file).dropna(subset=['Depth_80Gb'])

# %%
pd.set_option('display.max_columns', None)
print(lib_df)

# %%
lib_df['d_Depth'] = lib_df['Depth_80Gb'] - lib_df['Depth_40Gb']

# %%
df1 = lib_df[lib_df['batch']=='PC0174_A'].loc[:, ['batch','name', 'input', 'frag_total', 'Depth_40Gb', 'Depth_80Gb', 'd_Depth']]
df2 = lib_df[lib_df['batch']=='PC0174_B'].loc[:, ['batch', 'name', 'input', 'frag_total', 'Depth_40Gb', 'Depth_80Gb', 'd_Depth']]
df3 = lib_df[lib_df['batch']=='PC0174_C'].loc[:, ['batch', 'name', 'input', 'frag_total', 'Depth_40Gb', 'Depth_80Gb', 'd_Depth']]
data = df1.append(df2).append(df3)

print(data)

# %%
df1 = lib_df[lib_df['bs']=='BS0401_210929'].loc[:, ['Depth_40Gb', 'Depth_80Gb']]
df2 = lib_df[lib_df['batch']=='PC0174_B'].loc[:, ['Depth_40Gb', 'Depth_80Gb']]
df3 = lib_df[lib_df['batch']=='PC0174_C'].loc[:, ['Depth_40Gb', 'Depth_80Gb']]
data = df1.append(df2).append(df3)
df_melt = pd.melt(data)

print(df_melt)
# %%

# データ分布の可視化 --------------------------------
import matplotlib.pyplot as plt
import seaborn as sns

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

sns.boxplot(
    x='variable',
    y='value',
    data=df_melt,
    showfliers=False,
    width=0.4,
    ax=ax
)
sns.stripplot(
    x='variable',
    y='value',
    data=df_melt,
    jitter=False,
    color='black',
    ax=ax,
    size=3
)
plt.plot(
    ['Depth_40Gb', 'Depth_80Gb'],
    [data['Depth_40Gb'].to_list(), data['Depth_80Gb'].to_list()],
    color = 'black',
    linewidth = 1.0,
    linestyle = '--'
)
ax.set_xlim(-0.5, 1.5)
ax.set_ylim(0, data['Depth_80Gb'].max()+100)

# %%
# 対応のある2群の検定（Wilcoxonの符号付順位和検定）-----
"""
WilcoxonResult(statistic=0.0, pvalue=1.9073486328125e-06)
Depthが80Gbの時、有意に上昇している
"""
import numpy as np
from scipy import stats

A = data['Depth_40Gb'].to_numpy()
B = data['Depth_80Gb'].to_numpy()
res = stats.wilcoxon(A, B)

print(res)

# y軸:オッズ比ーx軸:seqID

# %%
# 2群の比率の検定（両側検定）-> BH法でFDR算出 ------------------
import statsmodels.stats.multitest as multi
from statsmodels.stats.proportion import proportions_ztest
import scipy.stats as st

# %%
x1_total = 12762
x2_total = 11291

x1 = 2 + 2 * np.random.randn(500).astype(int)
x1 = np.abs(x1)
x2 = 5 + 2 * np.random.randn(500).astype(int)
x2 = np.abs(x2)

print(x1)
print(x2)

# %%
# フィッシャー正確確率検定 ---------------------------
pval = []
for i in range(len(x1)):
    """
    p = proportions_ztest(
        [x1[i], x2[i]],
        [x1_total, x2_total],
        alternative='two-sided'
    )
    """
    x1u = x1_total - int(x1[i])
    x2u = x2_total - int(x2[i])
    p = st.fisher_exact(
        [[int(x1[i]), int(x2[i])], [x1u, x2u]],
        alternative='two-sided'
    )
    pval.append(p[1])

# %%
# カイ2乗分布による比率の差の検定(イェーツの補正あり) ---
# フィッシャーの正確確率検定では、
# サンプルサイズが増加するとp値の算出に時間がかかる
pval = []
for i in range(len(x1)):
    if (x1[i] == 0) | (x2[i] == 0):
        pval.append(1)
    else:
        x1u = x1_total - x1[i]
        x2u = x2_total - x2[i]
        _, p, _, _ = st.chi2_contingency(
            np.array([[x1[i], x2[i]], [x1u, x2u]]),
            correction=True
        )
        pval.append(p)

# %%
print(pval)

# %%
fdr = multi.multipletests(
    pval,
    alpha=0.05,
    method='fdr_bh'
)
x1x2_fdr = fdr[1]

# %%
print(x1x2_fdr)

# %%
"""test
p = proportions_ztest(
    [1, 4],
    [23777220, 20266770],
    alternative='two-sided'
)
print(p[1])

# %%
p = st.fisher_exact(
    [[1, 4], [23777220-1, 20266770-4]]
)
print(p)
print(p[1])

# %%
x2, pval, dof, e = st.chi2_contingency(
    [[1, 4], [23777220-1, 20266770-4]],
    correction=True
)
print(x2, pval, dof, e)
"""

# %%
x1_de = x1_total - x1
x2_de = x2_total - x2
x1_ratio = x1 / x1_de
x2_ratio = x2 / x2_de
x1x2_odds = x1_ratio / x2_ratio
# %%
x1x2_odds

# %%
# DataFrameの作成
id = ['Seq_'+ str(i+1) for i in range(len(x1))]
df = pd.DataFrame([id, x1, x2, x1x2_odds, x1x2_fdr]).transpose()
df.columns = ['SeqID', 'A', 'C', 'Odds-ratio', 'FDR']

#%%
df

# %%
df_sort = df.sort_values(by='Odds-ratio', ascending=False)
print(df_sort)
# %%
df_sort = df.sort_values(by='FDR', ascending=True)
print(df_sort)

# %%
# x軸:SeqID、y1軸:-log10(FDR)、y2軸:オッズ比とした時の2軸プロット -----
import matplotlib.pyplot as plt

fig=plt.figure()

ax1 = fig.subplots()
ax2 = ax1.twinx()

x = [i+1 for i in range(len(x1x2_odds))]
y1 = (-1) * np.log10(x1x2_fdr)
y2 = x1x2_odds

ax1.scatter(x, y1, color='red', s=1.2, label='-log10(FDR)')
ax2.scatter(x, y2, color='black', s=1.2, label='Odds-ratio')

ax1.set_ylim(0, y1.max()+2)
ax1.axhline((-1) * np.log10(0.05), color='gray', linestyle='--')
ax2.set_ylim(0, x1x2_odds.max()+1)

ax1.set_xlabel('SeqID')
ax1.set_ylabel('-log10(FDR)')
ax2.set_ylabel('Odds-ratio')

ax1.legend(
    bbox_to_anchor=(1.1, 1),
    loc='upper left',
    borderaxespad=0,
    fontsize=10
)
ax2.legend(
    bbox_to_anchor=(1.1, 0.90),
    loc='upper left',
    borderaxespad=0,
    fontsize=10
)
plt.show()

# %%
# x軸: オッズ比、y軸: -log10(FDR)とした時の散布図 ------------------
fig=plt.figure()
ax = fig.subplots()

x = x1x2_odds
y = (-1) * np.log10(x1x2_fdr)
ax.scatter(x, y, color='black', s=8)

ax.set_xlim(0, x1x2_odds.max()+1)
ax.set_ylim(0, y.max()+2)
ax.axvline(1, color='gray', linestyle='--')
ax.axhline((-1) * np.log10(0.05), color='gray', linestyle='--')

ax.set_xlabel('Odds-ratio')
ax.set_ylabel('-log10(FDR)')

plt.show()

# %%
(-1) * np.log10(0.05)

# %%
file = '/Users/tomoyauchiyama/code/CNN/test/TFBS_len.txt'
cnt = {}
with open(file, 'r') as rf:
    for line in rf.readlines():
        line = line.rstrip('\n')
        if line == 'V1':
            pass
        else:
            id = int(line)
            if cnt.get(id):
                cnt[id] += 1
            else:
                cnt[id] = 1

# %%
df = pd.read_csv(file)
df.describe()['V1']['min']

# %%
df.mean()

# %%
fig=plt.figure()
ax=fig.add_subplot(111)

#plt.hist(df, bins=20, rwidth=0.8, color='red', alpha=0.5, range=(0, 40), align='left')
x = sorted(cnt.keys())
y = [cnt[i] for i in x]
plt.bar(x, y, align='center', color='orange')
ax.set_title('Distribution of nr TFBS length (JASPAR_vertebrates)')
ax.set_xlabel('Length(nt)')
ax.set_ylabel('Count')
#plt.vlines(min(), 0, 300, colors='red', linestyle='-.', linewidth=2)
plt.vlines(df.mean(), 0, 300, colors='gray', linestyle='--', linewidth=2, label='mean_len=12.48 nt')
ax.set_xlim(0, max(x)+5)
ax.set_ylim(0, max(y)+20)
ax.legend(
    #bbox_to_anchor=(1.05, 1),
    loc='upper right',
    #borderaxespad=0,
    fontsize=10
)
plt.show()

# %%

def k_mer(k, nucl):
    if k == 1:
        k_mer()
        #return tmp
    else:
        tmp = []
        for i in nucl:
            for base in ['A', 'T', 'G', 'C']:
                seq = i + base
                tmp.append(seq)
        new_nucl = tmp
        k_mer(k-1, new_nucl)


# %%

nucl = ['A', 'T', 'G', 'C']
res = k_mer(4, nucl)

# %%
print(res)

# %%
for i in ['A', 'T', 'G', 'C']:
    seq = k_mer(4, i)
    print(seq)
# %%
print(k_mer(4))
# %%
ATGC = ['A', 'T', 'G', 'C']
tmp = []
for i in ATGC:
    for j in ATGC:
        seq = i + j
        tmp.append(seq)


# %%
# 5/2の作業 -----------------------------------------

def sameNucl_pos(k, ref):
    if k >= 0:
        init = 'N' * k
        tar = init + ref
    elif k < 0:
        cut = (-1) * k
        tar = ref[cut:]

    tmp = []
    for i in range(len(ref)):
        if len(tar) > i:
            if ref[i] == tar[i]:
                pos = i+1
                tmp.append(pos)

    if len(tmp) == 0:
        tmp = -1

    return tmp


# %%
# 1塩基のずれ ---------------
ref = 'ACTCTTCTGGT'
pos = sameNucl_pos(1, ref)
print(pos)
for i in pos:
    print(ref[i-1])

# %%
ref = 'ACTCTTCTGGT'
pos = sameNucl_pos(-1, ref)
print(pos)
for i in pos:
    print(ref[i-1])

# %%
# 2塩基のずれ ---------------
ref = 'ACTCTTCTGGT'
pos = sameNucl_pos(2, ref)
print(pos)
for i in pos:
    print(ref[i-1])

# %%
# 3塩基のずれ ---------------
ref = 'ACTCTTCTGGT'
pos = sameNucl_pos(3, ref)
print(pos)
for i in pos:
    print(ref[i-1])

# %%
"""test
# サンプルCで各塩基位置で出現頻度が最大値をとる塩基から構成されるDNA断片 -------------------
pos_freq = pd.read_excel('/Users/tomoyauchiyama/code/CNN/test/PR2688_C_60nt.xlsx')
maxnucl = pos_freq.idxmax(axis=1).tolist()
maxfrag = ''.join(maxnucl)

print(len(maxfrag))
print(maxfrag)

# 構築したDNA断片上で、5'UTR領域の開始位置を探索（一時的に推定)
ref = 'ACTCTTCTGGT'
five_utr = maxfrag.find(ref)
print(five_utr)
"""

"""test
# %%
df = pd.read_excel('/Users/tomoyauchiyama/code/CNN/test/PR2688_C_60nt.xlsx') * 10000
cnt = df.to_numpy().astype(int)
pwm_score = pwm(cnt, 10000)
pwm_df = pd.DataFrame(pwm_score, columns=df.columns)
print(pwm_df[49:])

# %%
ref = 'ACTCTTCTGGT'
length = len(ref)
start = 49
tar_pwm = pwm_df[start:start+length].reset_index(drop=True)
tar_min = tar_pwm.min(axis=1)
tar_max = tar_pwm.max(axis=1)
score = sum([tar_pwm.loc[i, nucl] for i, nucl in enumerate(ref)])
score_norm = (score - sum(tar_min)) / (sum(tar_max) - sum(tar_min))
print(score)
print(sum(tar_min))
print(sum(tar_max))
print(sum(tar_max) - sum(tar_min))
print(score_norm)
"""

# %%
# 全サンプルにおける5'UTR領域の開始位置の一時的な決定
import re
import glob

indir = '/Users/tomoyauchiyama/code/CNN/test'
ref = 'ACTCTTCTGGT'
five_utr = {}
for xlsx in sorted(glob.glob(os.path.join(indir, '*.xlsx'))):
    file_name = os.path.splitext(os.path.basename(xlsx))[0]
    name = re.compile(r'.*_(.*)_60nt').search(file_name).group(1)
    pos_freq = pd.read_excel(xlsx)
    maxnucl = pos_freq.idxmax(axis=1).tolist()
    maxfrag = ''.join(maxnucl)
    pos = maxfrag.find(ref)
    if pos != -1:
        five_utr[name] = pos
    else:
        flist = [m.start() for m in re.finditer('ACT', maxfrag)]
        if len(flist) == 0:
            five_utr[name] = -1
        else:
            for i in flist:
                if five_utr.get(name):
                    five_utr[name] += ',' + str(i)
                else:
                    five_utr[name] = str(i)

# %%
# calc PWM score
def pwm(cnt, total, bg=0.25):
    p = (cnt + (np.sqrt(total) * 0.25)) / (total + (4 * (np.sqrt(total) * 0.25)))
    pwm_score = np.log2(p/bg)

    return pwm_score

# %%
# pwmスコアの算出 ----------------------------
def pwmScore(pwm_df, start, length, ref):
    tar_pwm = pwm_df[start:start+length].reset_index(drop=True)
    ele = [nucl + '#' + str(tar_pwm.loc[i, nucl]) for i, nucl in enumerate(ref)]
    total = float(sum([tar_pwm.loc[i, nucl] for i, nucl in enumerate(ref)]))
    return ele, total

"""
スコアを出現頻度の対数によって定義した場合、スコアの和が配列の対数尤度に相当するため、
確率モデルに基づく推定と見なすことができる。
"""
# %%
# pwmスコアの正規化 ----------------------------
def pwmScore_norm(pwm_df, start, length, ref):
    tar_pwm = pwm_df[start:start+length].reset_index(drop=True)
    tar_min = tar_pwm.min(axis=1)
    tar_max = tar_pwm.max(axis=1)
    score = float(sum([tar_pwm.loc[i, nucl] for i, nucl in enumerate(ref)]))
    score_norm = (score - sum(tar_min)) / (sum(tar_max) - sum(tar_min))

    return score_norm

"""
pwm scoreの正規化(0~1)を行うことで、解釈をしやすくなると思ったが、複雑となったので却下。
pwm_norm(nuclFrag, i) = (pwm_score(nuclFrag) - min(i)) / (max(i) - min(i))
"""

# %%
def pwm_pattern(k, start, ref, file, total):
    pos_cnt = pd.read_excel(file) * total
    cnt = pos_cnt.to_numpy().astype(int)
    pwm_score = pwm(cnt, total)
    pwm_df = pd.DataFrame(pwm_score, columns=df.columns)

    score = {}
    for i in range(k+1):
        if i == 0:
            s = start
            #length = len(ref)
            length = len(ref) - 1
            ref = ref[:length]
            score[s] = pwmScore(pwm_df, s, length, ref)
            #score[s] = pwmScore_norm(pwm_df, s, length, ref)
        else:
            s_m = start - i
            length_m = len(ref)
            score[s_m] = pwmScore(pwm_df, s_m, length_m, ref)
            #score[s_m] = pwmScore_norm(pwm_df, s_m, length_m, ref)
            s_p = start + i
            #length_p = len(ref) - i
            length_p = len(ref) - (i+1)
            ref_p = ref[:length_p]
            score[s_p] = pwmScore(pwm_df, s_p, length_p, ref_p)
            #score[s_p] = pwmScore_norm(pwm_df, s_p, length_p, ref_p)

    return score

# %%
ref = 'ACTCTTCTGG'
start = 50
length = len(ref)
total = 21044437
file = '/Users/tomoyauchiyama/code/CNN/test/PR2688_C_60nt.xlsx'
pos_cnt = pd.read_excel(file) * total
cnt = pos_cnt.to_numpy().astype(int)
pwm_score = pwm(cnt, total)
pwm_df = pd.DataFrame(pwm_score, columns=pos_cnt.columns)
ele,  score = pwmScore(pwm_df, start, length, ref)

# %%
print(ele)
print(score)

# %%
# {'C': 49, 'D': '50', 'E': '50', 'F': 49, 'G': 49, 'H': 49}
ref = 'ACTCTTCTGGT'
file = '/Users/tomoyauchiyama/code/CNN/test/PR2688_C_60nt.xlsx'
total = 21044437
score = pwm_pattern(3, 49, ref, file, total)

# %%
file = '/Users/tomoyauchiyama/code/CNN/test/PR2688_D_60nt.xlsx'
total = 12028692
score = pwm_pattern(3, 50, ref, file, total)

# %%
file = '/Users/tomoyauchiyama/code/CNN/test/PR2688_E_60nt.xlsx'
total = 15283034
score = pwm_pattern(3, 50, ref, file, total)

# %%
file = '/Users/tomoyauchiyama/code/CNN/test/PR2688_F_60nt.xlsx'
total = 20323468
score = pwm_pattern(3, 49, ref, file, total)

# %%
file = '/Users/tomoyauchiyama/code/CNN/test/PR2688_G_60nt.xlsx'
total = 19745701
score = pwm_pattern(3, 49, ref, file, total)

# %%
file = '/Users/tomoyauchiyama/code/CNN/test/PR2688_H_60nt.xlsx'
total = 15582685
score = pwm_pattern(3, 49, ref, file, total)

# %%
# +側で、ライブラリー構造のATGC比率の調査領域を+3以上伸ばす。
print(sorted(score.items()))

# %%
fig=plt.figure()
ax = fig.subplots()

x = [i+1 for i in sorted(score.keys())]
y = [score[i-1] for i in x]
ax.scatter(x, y, color='black')
ax.plot(x, y, color='black', linestyle='--')

ax.set_xlim(min(x), max(x))
ax.set_ylim(min(y)-2, max(y)+2)
ax.axhline(0, color='gray', linestyle='--')
ax.axvline(50+1, color='red', linestyle='--', label='Predicted 5\'UTR start site')
#ax.axhline((-1) * np.log10(0.05), color='gray', linestyle='--')
ax.legend(
    bbox_to_anchor=(1.05, 1),
    loc='upper left',
    borderaxespad=0,
    fontsize=10
)
ax.set_xlabel('5\'UTR start position from 5\'index')
ax.set_ylabel('Score')

plt.show()

# %%
# 5/6の作業 -----------------------------------------
# randomにk-merをつくる -------
import itertools

bases = ['A', 'T', 'G', 'C']

# %%
k = 4
kmer_list = [''.join(p) for p in itertools.product(bases, repeat=k)]
print(len(kmer_list))
print(kmer_list)

# %%
# randomに作成された16塩基長の配列を100000個作る ------
# 正規分布に従うのか？
import random

def random_dna_sequence(length):
    return ''.join(random.choice('ACTG') for _ in range(length))

# %%
outfile = '/Users/tomoyauchiyama/code/CNN/test/jellyfish/random_seq.fa'
#N = 10000000
N = 100000
with open(outfile, 'w') as wf:
    for i in range(N):
        id = '>seq_' + str(i+1)
        seq = random_dna_sequence(16)
        wf.write(f'{id}\n')
        wf.write(f'{seq}\n')

# %%
# 作成した16塩基長の配列がランダム組成になっているかを確認 --------------
target_seq = '/Users/tomoyauchiyama/code/CNN/test/jellyfish/random_seq.fa'
nucl = pd.DataFrame()
flag = 0
with open(target_seq, 'rt') as fq:
    for line in fq.readlines():
        line = line.rstrip('\n')
        if re.compile(r'^>').search(line):
                flag = 1
        else:
            if flag == 1:
                insert = pd.Series([str(i) for i in line])
                nucl = nucl.append(insert, ignore_index=True)
                flag = 0

# %%
col = [str(i) for i in range(1, nucl.columns.size + 1)]
nucl.columns = col

# %%
# 事前確率　--------------------------------------------------------
total = nucl.index.size
pos_freq = pd.DataFrame()
for pos in col:
    col_list = list(nucl[pos])
    A = col_list.count('A') / total
    T = col_list.count('T') / total
    G = col_list.count('G') / total
    C = col_list.count('C') / total
    #N = col_list.count('N') / total
    freq = pd.Series([A, T, G, C])
    pos_freq = pos_freq.append(freq, ignore_index=True)

# %%
pos_freq.columns = ['A', 'T', 'G', 'C']
pos_freq.index += 1
print(pos_freq)

# %%
# logomekerでシーケンスロゴを作成-------------------------------------
color_scheme = {
    'T' : [0, 0.5, 0], # green
    'A' : [1, 0, 0], # red
    'G' : [1, 0.65, 0], # yellow
    'C' :[0, 0, 1], # blue
    'N': 'gray'
}
logo = lm.Logo(
    pos_freq,
    baseline_width=0.1,
    vpad=0.08,
    fade_probabilities=False,
    stack_order='small_on_top',
    color_scheme=color_scheme
    #font_name='Luxi Mono',
    #color_scheme='class',
    #font_name='Rosewood Std',
    #figsize=(30, 2.5)
)

logo.style_spines(spines=['left', 'right'], visible=False)

# style using Axes methods
logo.ax.set_xticks(np.arange(1, len(pos_freq)+1, step=2))
logo.ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
logo.ax.set_xlabel('Position')
logo.ax.set_ylabel('Probability')
plt.show()

# %%
hist_file = '/Users/tomoyauchiyama/code/CNN/test/B/histogram_12.txt'
hist_df = pd.read_csv(hist_file, delimiter=' ', header=None)
hist_file = '/Users/tomoyauchiyama/code/CNN/test/H/histo_12_random.txt'
hist_df2 = pd.read_csv(hist_file, delimiter=' ', header=None)
hist_file = '/Users/tomoyauchiyama/code/CNN/test/H/histo_12.txt'
hist_df3 = pd.read_csv(hist_file, delimiter=' ', header=None)

# %%
"""
hist_file = '/Users/tomoyauchiyama/code/CNN/test/H_hist/histo_7.txt'
hist_df3 = pd.read_csv(hist_file, delimiter=' ', header=None)
hist_file = '/Users/tomoyauchiyama/code/CNN/test/H_hist/histo_8.txt'
hist_df4 = pd.read_csv(hist_file, delimiter=' ', header=None)
hist_file = '/Users/tomoyauchiyama/code/CNN/test/H_hist/histo_9.txt'
hist_df5 = pd.read_csv(hist_file, delimiter=' ', header=None)
hist_file = '/Users/tomoyauchiyama/code/CNN/test/H_hist/histo_10.txt'
hist_df6 = pd.read_csv(hist_file, delimiter=' ', header=None)
hist_file = '/Users/tomoyauchiyama/code/CNN/test/H_hist/histo_11.txt'
hist_df7 = pd.read_csv(hist_file, delimiter=' ', header=None)
hist_file = '/Users/tomoyauchiyama/code/CNN/test/H_hist/histo_12.txt'
hist_df8 = pd.read_csv(hist_file, delimiter=' ', header=None)
"""

# %%

fig=plt.figure()
ax=fig.add_subplot(111)

plt.bar(hist_df[0], hist_df[1], align='center', width=1.0, label='k=12(B)')
#plt.scatter(hist_df[0], hist_df[1], label='k=5')
#plt.bar(hist_df2[0], hist_df2[1], align='center', width=1.0, label='k=12(random)')
#plt.scatter(hist_df2[0], hist_df2[1], label='k=6')
plt.bar(hist_df3[0], hist_df3[1], alpha=0.5, align='center', width=1.0, label='k=12(H)')

#plt.bar(hist_df4[0], hist_df4[1], alpha=0.5, align='center', width=1.0, label='k=8')
#plt.bar(hist_df5[0], hist_df5[1], alpha=0.5, align='center', width=1.0, label='k=9')
#plt.bar(hist_df6[0], hist_df6[1], alpha=0.5, align='center', width=1.0, label='k=10')
#plt.bar(hist_df7[0], hist_df7[1], alpha=0.5, align='center', width=1.0, label='k=11')
#plt.bar(hist_df8[0], hist_df8[1], alpha=0.5, align='center', width=1.0, label='k=12')

#ax.set_xlim(-500, 14000)
plt.xlabel('# of occurance at k-mer')
plt.ylabel('# of unique k-mers')
ax.legend(
    #bbox_to_anchor=(1.05, 1),
    loc='upper right',
    #borderaxespad=0,
    fontsize=8
)
#plt.ylim(0, higher_frequency)
plt.show()

# %%
# 16塩基内に局在するk-mer配列の開始位置の同定 -------------------------
km_file = '/Users/tomoyauchiyama/code/CNN/test/jellyfish/k6.txt'
kmer_df = pd.read_csv(km_file, delimiter=' ', header=None)
kmer_list = kmer_df[0]

outfile = '/Users/tomoyauchiyama/code/CNN/test/jellyfish/k6.fa'
with open(outfile, 'w') as wf:
    for i in range(len(kmer_list)):
        id = '>seq_' + str(i+1)
        seq = kmer_list[i]
        wf.write(f'{id}\n')
        wf.write(f'{seq}\n')


# %%
def kmer_start(kmer, dna):
    tmp = []
    for m in re.finditer(kmer, dna):
        start = str(m.start() + 1)
        tmp.append(start)

    return ','.join(tmp)

# %%
fa_file = '/Users/tomoyauchiyama/code/CNN/test/jellyfish/random_seq.fa'
kmer_pos = {}
for kmer in kmer_list:
    with open(fa_file, 'r') as rf:
        for line in rf.readlines():
            line = line.rstrip('\n')
            if re.compile(r'^>').search(line):
                pass
            else:
                if kmer in line:
                    pos = kmer_start(kmer, line)
                    if kmer_pos.get(kmer):
                        kmer_pos[kmer] += ',' + pos
                    else:
                        kmer_pos[kmer] = pos

# %%
print(len(kmer_pos))
print(kmer_pos)

# %%
# 上記のスクリプトで実行すると、17分かかるのでseqkitで実行することにした。
# usage: https://bioinf.shenwei.me/seqkit/usage/#locate
# conda install -c bioconda seqkit
# jellyfish dump mer_counts.jf -c > k6.txt
# 16塩基内に局在するk-mer配列の開始位置の同定 -------------------------
"""
km_file = '/Users/tomoyauchiyama/code/CNN/test/jellyfish/k6.txt'
kmer_df = pd.read_csv(km_file, delimiter=' ', header=None)
kmer_list = kmer_df[0]

outfile = '/Users/tomoyauchiyama/code/CNN/test/jellyfish/k6.fa'
with open(outfile, 'w') as wf:
    for i in range(len(kmer_list)):
        id = '>seq_' + str(i+1)
        seq = kmer_list[i]
        wf.write(f'{id}\n')
        wf.write(f'{seq}\n')
"""
# time cat random_seq.fa | seqkit locate -P -f k6.fa -j 30 > k6_pos.txt
# user 16m59.860s (default thread:4) - 自作のスクリプト実行と変わらないが、thread数をあげると実行速度が上がるだろう
# jellyfishによるcount結果とseqkitによるヒット数の結果は一致している
# %%
# seqkit locate 出力ファイルの処理 ------------------------------------
infile = '/Users/tomoyauchiyama/code/CNN/test/jellyfish/k6_pos.txt'
"""
seqID   patternName     pattern strand  start   end     matched
seq_1   seq_451 ACTAAG  +       5       10      ACTAAG
"""
kmerCnt = {}
with open(infile, 'r') as rf:
    for line in rf.readlines():
        line = line.rstrip('\n')
        if 'seqID' in line:
            pass
        else:
            tmp = line.split()
            id = tmp[2] + ',' + str(tmp[4])
            if kmerCnt.get(id):
                kmerCnt[id] += 1
            else:
                kmerCnt[id] = 1
            """
            if kmerStart.get(kmer):
                kmerStart[kmer] += ',' + pos
            else:
                kmerStart[kmer] = pos
            """

# %%
"""test
print(kmerCnt['AAAAGG,2'])

pos= ['AAAAA,19', 'AAAAA,1', 'AAAAA,10', 'AAAAA,5', 'AAAAA,8']
print(sorted(pos, key=lambda x: int(re.search(r'.*,(\d+)', x).group(1))))
"""

# %%
# グラフによる可視化の準備 ------------------------------
y = {}
x = {}
for id in sorted(kmerCnt.keys(), key=lambda x: int(re.search(r'.*,(\d+)', x).group(1))):
    kmer, xval = id.split(',')
    if y.get(kmer):
        y[kmer] += ',' + str(kmerCnt[id])
        x[kmer] += ',' + str(xval)
    else:
        y[kmer] = str(kmerCnt[id])
        x[kmer] = str(xval)

kmer_df = pd.DataFrame()
for kmer in sorted(y.keys()):
    df1 = pd.Series([kmer, [int(i) for i in x[kmer].split(',')]])
    df2 = pd.Series([[int(i) for i in y[kmer].split(',')]])
    col = pd.concat([df1, df2])
    kmer_df = kmer_df.append(col, ignore_index=True)

kmer_df.columns = ['kmer', 'pos', 'count']

# %%
"""test
kmer_df
kmer_df[kmer_df['kmer'] == 'AAAAGG']
"""
# %%
# あるk-mer配列におけるターゲット配列上の出現位置とその回数 ------
# 例として、最も出現回数が多かったAAAAGGの場合を示す

seq = 'TGCGCC'

fig=plt.figure()
ax=fig.add_subplot(111)

idx = kmer_df.index[kmer_df['kmer'] == seq][0]
x_val = kmer_df.loc[idx, 'pos']
y_val = kmer_df.loc[idx, 'count']
plt.bar(x_val, y_val, align='center', color='orange')
ax.set_title(f'{seq}')
ax.set_xlabel('Position')
ax.set_ylabel('Count')
#plt.vlines(min(), 0, 300, colors='red', linestyle='-.', linewidth=2)
#plt.vlines(df.mean(), 0, 300, colors='gray', linestyle='--', linewidth=2, label='mean_len=12.48 nt')
ax.set_xlim(0, max(x_val)+1)
ax.set_ylim(0, max(y_val)+10)
"""
ax.legend(
    #bbox_to_anchor=(1.05, 1),
    loc='upper right',
    #borderaxespad=0,
    fontsize=10
)
"""
plt.xticks(x_val)
plt.show()
# random性が確認できた

# %%
# 目的配列の切り出し
df = pd.read_excel('/Users/tomoyauchiyama/code/CNN/test/BG_k10_top100.xlsx')
df['Seq'].head(3)

# %%
nucl = pd.DataFrame()
for seq in df['Seq']:
    insert = pd.Series([str(i) for i in seq])
    nucl = nucl.append(insert, ignore_index=True)

col = [str(i) for i in range(1, nucl.columns.size + 1)]
nucl.columns = col

total = nucl.index.size
pos_freq = pd.DataFrame()
for pos in col:
    col_list = list(nucl[pos])
    A = col_list.count('A') / total
    T = col_list.count('T') / total
    G = col_list.count('G') / total
    C = col_list.count('C') / total
    #N = col_list.count('N') / total
    freq = pd.Series([A, T, G, C])
    pos_freq = pos_freq.append(freq, ignore_index=True)

pos_freq.columns = ['A', 'T', 'G', 'C']
pos_freq.index += 1
print(pos_freq)

# %%
color_scheme = {
    'T' : [0, 0.5, 0], # green
    'A' : [1, 0, 0], # red
    'G' : [1, 0.65, 0], # yellow
    'C' :[0, 0, 1], # blue
    'N': 'gray'
}
logo = lm.Logo(
    pos_freq,
    baseline_width=0.1,
    vpad=0.08,
    fade_probabilities=False,
    stack_order='small_on_top',
    color_scheme=color_scheme
    #font_name='Luxi Mono',
    #color_scheme='class',
    #font_name='Rosewood Std',
    #figsize=(30, 2.5)
)

logo.style_spines(spines=['left', 'right'], visible=False)

# style using Axes methods
logo.ax.set_xticks(np.arange(1, len(pos_freq)+1, step=1))
logo.ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
logo.ax.set_xlabel('Position')
logo.ax.set_ylabel('Frequency')
plt.show()


# %%
# 8/29 -------------------------------------
# ToDo: p値の補正をbofferoniであることを確認すること

import os
import re
import sys
import pandas as pd
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# %%
# read data ---
infile = '/Users/tomoyauchiyama/code/PG4699_target_ratio.txt'
mesh2D = pd.read_csv(infile, sep='\t')
# %%
mesh2D.head(5)

# %%
# rep1 freq data ---
rep1_a = mesh2D.loc[:, 'PG4699_01_a_Freq':'PG4699_20_a_Freq']
rep1_b = mesh2D.loc[:, 'PG4699_43_a_Freq']
rep1_freq = pd.concat([rep1_a, rep1_b], axis=1)
rep1_freq.head(5)
# %%
rep1_freq.shape
# %%
# paired data ---
# paired-sample(controlのカウント値は1以上)
rep1_freq_m = rep1_freq[rep1_freq['PG4699_43_a_Freq'] != 0]
rep1_freq_m.head(5)

# %%
# rep2 freq data ---
rep2_a = mesh2D.loc[:, 'PG4699_21_a_Freq':'PG4699_42_a_Freq']
rep2_b = mesh2D.loc[:, 'PG4699_44_a_Freq']
rep2_freq = pd.concat([rep2_a, rep2_b], axis=1)
rep2_freq.head(5)
# %%
rep2_freq.shape
# %%
# paired data ---
# paired-sample(controlのカウント値は1以上)
rep2_freq_m = rep2_freq[rep2_freq['PG4699_44_a_Freq'] != 0]
rep2_freq_m.head(5)
# %%
rep2_freq_m.shape

# %%
# rep1 and rep2 FC ---
FC_col_1 = []
FC_col_2 = []
for cols in mesh2D.columns:
    if '43_a_FC' in cols:
        FC_col_1.append(cols)
    elif '44_a_FC' in cols:
        FC_col_2.append(cols)
# %%
rep1 = pd.DataFrame()
for cols in FC_col_1:
    name = re.compile(r'(.*)_43_a_FC').search(cols).group(1)
    ID = name + '_a'
    rep1[ID] = mesh2D.loc[rep1_freq_m.index.tolist(), cols]
# %%
rep1.shape
# %%
rep2 = pd.DataFrame()
for cols in FC_col_2:
    name = re.compile(r'(.*)_44_a_FC').search(cols).group(1)
    ID = name + '_a'
    rep2[ID] = mesh2D.loc[rep2_freq_m.index.tolist(), cols]
# %%
rep2.shape

# %%
# zscore ---
# サンプル間で比較するため、サンプル間でスケールを合わせる
z_rep1 = st.zscore(rep1, axis=1)
z_rep2 = st.zscore(rep2, axis=1)

# %%
# 比較対象とはならない欠損値がある行、配列は調査外とする
z_rep1_m = z_rep1.dropna()
z_rep2_m = z_rep2.dropna()

# %%
# make heatmap ---
y = z_rep1_m.index.values
x = z_rep1_m.columns.values

fig = go.Figure(
    data=go.Heatmap(
        colorbar=dict(
            title='z-score'
        ),
        z=z_rep1_m,
        x=x,
        y=y,
        colorscale='Viridis'
    )
)

fig.update_xaxes(title='Sample(normalized)')
fig.update_yaxes(title='Seq_ID')

fig.update_layout(
    title='heatmap',
    xaxis_nticks=36)

fig.show()
htmlfile = os.path.join(
    os.path.abspath('.'),
    'PG4699_rep1_heatmap.html'
)
fig.write_html(htmlfile)

# %%
# make heatmap ---
y = z_rep2_m.index.values
x = z_rep2_m.columns.values

fig = go.Figure(
    data=go.Heatmap(
        colorbar=dict(
            title='z-score'
        ),
        z=z_rep2_m,
        x=x,
        y=y,
        colorscale='Viridis'
    )
)

fig.update_xaxes(title='Sample(normalized)')
fig.update_yaxes(title='Seq_ID')

fig.update_layout(
    title='heatmap',
    xaxis_nticks=36)

fig.show()
htmlfile = os.path.join(
    os.path.abspath('.'),
    'PG4699_rep2_heatmap.html'
)
fig.write_html(htmlfile)

# %%
# サンプル間のz-scoreを用いて、サンプル間の類似度を確認(クラスタリング)
import scipy.spatial.distance as distance
from scipy.cluster.hierarchy import dendrogram, ward

#ward法で分類
linkage_array = ward(z_rep1_m.T)

ax = plt.figure(figsize=(20,10)).gca()
dendrogram(linkage_array)
bounds = ax.get_xbound()

plt.xlabel("sample index",fontsize=10)
plt.ylabel("Cluster distance",fontsize=10)

# %%
# サンプルの散布図 ---
def scatter_with_dist(df, x_label, y_label, outfile):
    #file = r'C:\Users\T22178\Desktop\test/BG_k9.xlsx'
    #df = pd.read_excel(file)

    #file = r'C:\Users\T22178\Desktop\test/BG_k9_top100.xlsx'
    #df_top = pd.read_excel(file)
    #top100_xy = df[df['Seq'].isin(df_top['Seq'].to_list())]

    left, width = 0.1, 0.65 #散布図の左余白と、幅を定義
    bottom, height = 0.1, 0.65 #散布図の下余白と、高さを定義
    spacing = 0.05 #散布図とヒストグラムの間隔
    fig_scatter = [left, bottom, width, height] #散布図の定義
    fig_histx = [left, bottom + height + spacing, width, 0.2] #上のヒストグラムの定義
    fig_histy = [left + width + spacing, bottom, 0.2, height] #右のヒストグラムの定義

    fig = plt.figure(figsize=(8, 8)) #描写領域を正方形で生成

    ax = fig.add_axes(fig_scatter) #図1として散布図領域を生成
    ax_histx = fig.add_axes(fig_histx, sharex=ax) #上のヒストグラム領域を生成
    ax_histy = fig.add_axes(fig_histy, sharey=ax) #右のヒストグラム領域を生成
    ax_histx.tick_params(axis="x", labelbottom=False) #x軸ラベル無し
    ax_histy.tick_params(axis="y", labelleft=False) #y軸ラベル無し

    x = df[x_label].to_numpy()
    y = df[y_label].to_numpy()
    xymax = max(np.max(x), np.max(y)) + 1 #xとyの最大値を定義
    binwidth = 0.05 #binの値幅を定義。
    bins = np.arange(0, xymax, binwidth) #binの個数を定義

    ax.set_xlim([0, xymax])
    ax.set_ylim([0, xymax])
    ax.scatter(x, y, color='gray', s=2.0)
    #ax.scatter(top100_xy[x_label], top100_xy[y_label], color='red', s=2.0, label='top100')
    ax.plot([0, xymax], [0, xymax], color='darkturquoise', linestyle='--', label='FC=1')
    ax.plot(bins, 1.5*bins, color='dodgerblue', linestyle='--', label='FC=1.5')
    ax.plot(bins, 2.0*bins, color='blue', linestyle='--', label='FC=2')

    #r = df.loc[:, x_label:y_label].corr(method='pearson').iloc[0, 1]
    #lm = LinearRegression().fit(x.reshape(-1, 1), y)
    #ax.plot(x, lm.predict(x.reshape(-1, 1)), color='purple', linestyle='--', label='r='+str(round(r, 2)))

    ax_histx.hist(x, bins=bins, color='gray') #xについてのヒストグラム作成
    ax_histy.hist(y, bins=bins, orientation='horizontal', color='gray') #yについてのヒストグラム作成

    #ax.set_xlabel(x_label)
    #ax.set_ylabel(y_label)
    ax.set_xlabel(f'log10(({x_label})x10000)')
    ax.set_ylabel(f'log10(({y_label})x10000)')

    ax.legend(
        bbox_to_anchor=(1.3, 1.4),
        loc='upper right',
        #borderaxespad=0,
        fontsize=8
    )

    plt.tick_params(labelsize=8)
    plt.show()
    plt.savefig(outfile)

# %%
#paired-sample(controlのカウント値は1以上)　*生データのスケールを合わせた場合
outdir = '/Users/tomoyauchiyama/code/CNN/test/PG4699'
df_10000 = rep1_freq_m  * 10000
df_10000_m = df_10000.replace(0, 1)
df_log10 = np.log10(df_10000_m)
x_label = 'PG4699_43_a_Freq'
y_label = 'PG4699_01_a_Freq'
name1 = re.compile(r'(.*)_a').search(y_label).group(1)
name2 = re.compile(r'PG.*_(.*_a)_Freq').search(x_label).group(1)
name = name1 + '_' + name2 + '_scatter.png'
outfile = os.path.join(outdir, name)
scatter_with_dist(df_log10, x_label, y_label, outfile)

# %%
#paired-sample(controlのカウント値は1以上)　*生データの場合
outdir = '/Users/tomoyauchiyama/code/CNN/test/PG4699'
df_10000 = rep1_freq_m
for cols in df.columns.values:
    x_label = 'PG4699_43_a_Freq'
    if cols != x_label:
        y_label = cols
        name1 = re.compile(r'(.*)_a').search(y_label).group(1)
        name2 = re.compile(r'PG.*_(.*_a)_Freq').search(x_label).group(1)
        name = name1 + '_' + name2 + '_scatter.png'
        outfile = os.path.join(outdir, name)
        scatter_with_dist(df, x_label, y_label, outfile)

# %%
#paired-sample(controlのカウント値は1以上)
df = rep2_freq_m
for cols in df.columns.values:
    x_label = 'PG4699_44_a_Freq'
    if cols != x_label:
        y_label = cols
        name1 = re.compile(r'(.*)_a').search(y_label).group(1)
        name2 = re.compile(r'PG.*_(.*_a)_Freq').search(x_label).group(1)
        name = name1 + '_' + name2 + '_scatter.png'
        outfile = os.path.join(outdir, name)
        scatter_with_dist(df, x_label, y_label, outfile)

# %%
# 3D plot ----
rep1_fcfdr = mesh2D.loc[
    rep1_freq_m.index.tolist(),
    'PG4699_01_43_a_FC':'PG4699_20_43_a_FDR'
].reset_index(drop=True)
rep1_fcfdr['SeqID'] = pd.Series([f'Seq_'+ str(i) for i in rep1_freq_m.index.tolist()])
log10FDR = (-1) * np.log10(rep1_fcfdr['PG4699_01_43_a_FDR'].replace(0, 1e-300))
data = pd.concat([rep1_fcfdr['SeqID'], rep1_fcfdr['PG4699_01_43_a_FC'], log10FDR], axis=1)

# %%
# ToDo: p値の補正をbofferoni法であることを確認すること

outdir = '/Users/tomoyauchiyama/code/CNN/test/PG4699'

mesh_size = 1.5
#xrange = np.arange(0, mesh2D.index.size, mesh_size)
#yrange = np.arange(0, mesh2D.columns.size, mesh_size)

fig = px.scatter_3d(
    data,
    y='SeqID',
    x='PG4699_01_43_a_FDR',
    z='PG4699_01_43_a_FC'
    #color='batch'
)
fig.update_layout(
    showlegend=True,
    legend=dict(
        x=-0.1,
        xanchor='left',
        y=1,
        yanchor = 'auto'
    )
)
fig.update_traces(
    marker=dict(
        size=1.0,
        line=dict(width=1.0, color='DarkSlateGrey')
        #color='white'
    ),
    selector=dict(mode='markers')
)
"""
px0 = np.linspace(0, max(X), 5)
px1 = np.linspace(0, max(Y), 5)
px0, px1 = np.meshgrid(px0, px1)
z = w[0]*px0 + w[1] * px1 + w[2]
fig.add_traces(
    go.Surface(
        y=px0,
        x=px1,
        z=z,
        colorscale='Viridis',
        colorbar=dict(title='d_depth')
    )
)

fig.update_layout(
    title='display 3D Surface Plots',
    xaxis_nticks=36
)
fig.update_layout(showlegend=True)

f = str(round(w[0], 2)) + 'x+' + str(round(w[1], 2)) + 'y' + str(round(w[2], 2))
text = '*f(x,y) =' + f
fig.add_annotation(
    x=-0.1,
    y=0.05,
    text=text,
    font=dict(size=8),
    showarrow=False,
    arrowhead=1
)

text = '*SD = ' + str(round(np.sqrt(mse), 2))
fig.add_annotation(
    x=-0.1,
    y=0,
    text=text,
    font=dict(size=8),
    showarrow=False,
    arrowhead=1,
)
"""
fig.show()
htmlfile = os.path.join(
    outdir,
    'PG4699_01_43_a_3D_Plots.html'
)
fig.write_html(htmlfile)

# %%
# x軸:SeqID、y軸:オッズ比とした時のプロット -----
fig = plt.figure()
ax1 = fig.subplots()

x = rep1_freq_m.index
y = rep1_fcfdr['PG4699_01_43_a_FC']

ax1.scatter(x, y, color='black', s=1.2, label='Fold-change')
ax1.set_ylim(0, y.max()+1)

ax1.set_xlabel('SeqID')
ax1.set_ylabel('Fold-change')

plt.show()

# %%
# cheak ----
# FC=281のサンプルを確認
data[data['SeqID'] == 'Seq_1345002']
rep1_freq_m.iloc[86011, :]
rep1_freq_m.describe()

# %%
outdir = '/Users/tomoyauchiyama/code/CNN/test/PG4699'

main_title = 'FC'
x = rep1_freq_m.index
min_x = x.min()
max_x = x.max()
y = rep1_fcfdr['PG4699_01_43_a_FC']
min_y = y.min()
max_y = y.max() + 10


fig = make_subplots()
fig.add_trace(
    go.Scatter(
        x=x,
        y=y,
        mode='markers',
        marker=dict(
            color='gray',
            size=2.0,
            line=dict(
                color='black',
                width=0.7
            )
        )
    )
)
fig.update_layout(
    plot_bgcolor='white'
    #height=800,
    #width=900
)
fig.update_xaxes(
    title = 'Seq_ID',
    showline=True,
    linewidth=1,
    linecolor='black',
    mirror=True,
    ticks='inside',
    range=(min_x, max_x)
)
fig.update_yaxes(
    title = 'Fold-change',
    showline=True,
    linewidth=1,
    linecolor='black',
    mirror=True,
    ticks='inside',
    range=(min_y, max_y)
)

fig.for_each_xaxis(
    lambda axis: axis.title.update(
        font=dict(
            color='black',
            size=15
        )
    )
)
fig.for_each_yaxis(
    lambda axis: axis.title.update(
        font=dict(
            color='black',
            size=15
        )
    )
)
fig.update_annotations(
    font=dict(size=12),
)
fig.update_layout(
    title=dict(
        text=main_title,
        x=0.5,
        xanchor='center'
    ),
    showlegend=False
)
#fig.show()
htmlfile = os.path.join(
    outdir,
    'PG4699_01_43_a_FC_plot.html'
)
fig.write_html(htmlfile)

# %%
# make heatmap ---
fig = make_subplots(
    rows=1,
    cols=2,
    y_title='Seq_ID',
    x_title='Sample(normalized)',
    column_widths=[30, 1],
    horizontal_spacing=0.0001,
    shared_yaxes=True
)


y = rep1_filter['Header_ID']
x = rep1_filter.iloc[:, 1:len(rep1_filter.columns) - 1].columns.values
z = rep1_filter.iloc[:, 1:len(rep1_filter.columns) - 1]

fig.add_trace(
    go.Heatmap(
        colorbar=dict(
            title='z-score'
        ),
        z=z,
        x=x,
        y=y,
        colorscale='Viridis'
    ),
    row=1,
    col=1
)

fig.add_trace(
    go.Heatmap(
        z=rep1_filter.loc[:, ['cluster']],
        x=rep1_filter.loc[:, ['cluster']].columns.values,
        colorscale='rainbow',
        showscale=False
    ),
    row=1,
    col=2
)

fig.update_layout(
    #plot_bgcolor='white'
    height=800,
    width=700,
    font=dict(
        size=5
    )
)

fig.update_annotations(
   font=dict(size=10),
)

fig.update_yaxes(autorange='reversed')

fig.update_layout(
    title=dict(
        text='Brain_High_heatmap_rep1',
        x=0.5,
        xanchor='center'
    ),
)


fig.show()
htmlfile = os.path.join(
    '/Users/tomoyauchiyama/code/CNN/test/PG4699',
    'PG4699_Brain_High_heatmap_rep1.html'
)
fig.write_html(htmlfile)