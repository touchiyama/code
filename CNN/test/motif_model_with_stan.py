# %%
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
target_seq = '/Users/tomoyauchiyama/code/Consensus_frequency/12B2106/insert.fa'
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
    fade_probabilities=True,
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