# %%
import h5py
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# %%

"""目的
全長131072baseからなる配列に対して、128baseごとにChIPによるカバレッジ値が記載されたデータ
を1次元の状態空間モデルで表現する。
"""

"""h5ファイル構造
・HDF5ファイルは，階層的にデータを格納することができ，
　行列やテンソルデータをそれぞれの位置で名前付きで格納できる
・with構文直下で、ファイルの中身を操作できる。
・keyでデータを取得
"""

infile = '/Users/tomoyauchiyama/code/CNN/seq.h5'

with h5py.File('seq.h5', 'r') as hf:
    for key in hf.keys():
        print(key, hf[key].shape, hf[key].dtype)

"""output
target_labels (10,) |S29
test_in (500, 131072, 4) bool
test_out (500, 1024, 10) float16
train_in (5000, 131072, 4) bool
train_out (5000, 1024, 10) float16
valid_in (500, 131072, 4) bool
valid_out (500, 1024, 10) float16

・*_outに、10種類のChIP-seqの結果の1024bp毎のカバレッジ値が
　格納(5000セット)
"""

with h5py.File(infile, 'r') as h5:
    #for key in h5.keys(
    print(h5['test_out'][0, :])

"""1セット目の10種のChIP-seqデータ
[[1.417    0.3848   0.121    ... 1.691    0.       5.47    ]
 [2.098    0.4763   0.811    ... 2.877    0.       2.102   ]
 [1.057    0.       0.6177   ... 1.28     0.       1.48    ]
 ...
 [2.781    4.715    0.       ... 0.03041  1.3      1.687   ]
 [0.013245 7.195    0.       ... 1.354    1.794    0.9497  ]
 [0.2449   5.56     0.       ... 0.4526   1.958    1.149   ]]

・縦方向（行）が塩基の位置、横方向（列）が、各種ChIP-seqデータ
　を表している
"""

# %%
with h5py.File(infile, 'r') as h5:
    #for key in h5.keys(
    print(h5['test_out'][:,:,:].min())

# %%
import matplotlib.pyplot as plt
print(plt.rcParams["figure.figsize"])

# %%
# 3種のChIP-seqデータ(テストデータ)の可視化 ----------------------
fig = make_subplots()
K = 3
with h5py.File(infile, 'r') as h5:
    for i in range(K):
        cov = h5['test_out'][0, :, i]
        x = [i+1 for i in range(len(cov))]
        y = cov
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                name='ChIP_' + str(i+1),
                mode='lines',
                line_width=0.5,
                fill='tozeroy'
                #mode='markers',
                #marker=dict(
                #    size=5,
                #    line=dict(color='black', width=1)
                #)
            )
        )

    fig.update_layout(
        plot_bgcolor='white'
    )
    min_x = min(x)
    max_x = max(x)
    fig.update_xaxes(
        title='Position',
        showline=True,
        linewidth=1,
        linecolor='black',
        mirror=True,
        ticks='inside',
        range=(min_x, max_x)
        #dtick=1.0
    )
    min_y = h5['test_out'][:,:,:].min()
    max_y = h5['test_out'][:,:,:].max()
    fig.update_yaxes(
        title='Covarage',
        showline=True,
        linewidth=1,
        linecolor='black',
        mirror=True,
        ticks='inside',
        range=(min_y, max_y)
        #dtick=1.0
    )
    fig.for_each_xaxis(
        lambda axis: axis.title.update(
            font=dict(
                color='black',
                size=10
            )
        )
    )
    fig.for_each_yaxis(
        lambda axis: axis.title.update(
            font=dict(
                color='black',
                size=10
            )
        )
    )
    fig.update_layout(
        title=dict(
            text='test_data',
            x=0.5,
            xanchor='center'
        ),
        showlegend=True
    )

    fig.update_annotations(
        font=dict(size=10)
    )

    fig.show()

# %%

# ChIP-seqデータのカバレッジ値を連続値から離散値に変更
# ChIP-seqデータは0以上のカウントデータなので、ポアソン分布に従うと仮定する
# 連続値ならば、ガンマ分布でいけるだろうか。。。

with h5py.File(infile, 'r') as h5:
    test_out = h5['test_out'][:, :, :]
    train_out = h5['train_out'][:, :, :]
    valid_out = h5['valid_out'][:, :, :]

# %%
# data processing
test_num, test_pos, test_type = test_out.shape
train_num, train_pos, train_type = train_out.shape
valid_num, valid_pos, valid_type = valid_out.shape

#test_sum = np.round(np.log2(np.sum(test_out, axis=0))).astype(int)
#train_sum = np.round(np.log2(np.sum(train_out, axis=0))).astype(int)
#valid_sum = np.round(np.log2(np.sum(valid_out, axis=0))).astype(int)

# 一番最初のデータを使う
test_data = np.round(test_out[0, :, :].astype(int)).transpose()
train_data = np.round(train_out[0, :, :].astype(int)).transpose()
valid_data = np.round(valid_out[0, :, :].astype(int)).transpose()


# %%
# 250番目のデータを使う
test_data = np.round(test_out[249, :, :].astype(int)).transpose()
train_data = np.round(train_out[249, :, :].astype(int)).transpose()
valid_data = np.round(valid_out[249, :, :].astype(int)).transpose()

# %%
print(test_data.max())
print(train_data.max())
print(valid_data.max())

# %%
print(train_data[0, :].max())

# %%
#--------------------------------------
# make json file
# 各10種のCHIP-seqデータに対応する
# -> 個々のパラメータを推定するような操作
#--------------------------------------
import os
import json

type, N = test_data.shape

for i in range(type):
    observed_data = dict(
        N=N,
        Y=test_data[i, :].tolist()
    )

    json_name = 'ChIP_' + str(i + 1) # have to change

    json_file = os.path.join(
        '/Users/tomoyauchiyama/code',
        'STAN',
        'CHIP_data',
        json_name + '.json'
    )

    with open(json_file, 'w') as wf:
        json.dump(observed_data, wf, indent=2)

    print(json_file)

# %%
#--------------------------------------
# make json file
# 全10種のCHIP-seqに対応
# -> 10種に共通するパラメータを推定する操作
#--------------------------------------
import os
import json

type, N = test_data.shape

observed_data = dict(
    Ni=N,
    Nt=type,
    Y=test_data.tolist()
)

json_name = 'ChIP_merge_250' # have to change

json_file = os.path.join(
    '/Users/tomoyauchiyama/code',
    'STAN',
    'CHIP_data',
    json_name + '.json'
)

with open(json_file, 'w') as wf:
    json.dump(observed_data, wf, indent=2)

print(json_file)

# %%
# stan file ---------------------------------
# -> 個々のパラメータを推定する方法

code = """data {
    int N;
    int Y[N];
}
parameters {
    vector[N] r;
    real<lower=0> s_r;
}
model {
    target += normal_lpdf(r[3:N] | 2*r[2:(N-1)] - r[1:(N-2)], s_r);
    Y ~ poisson_log(r);
}
generated quantities {
  vector[N] Y_mean;
  Y_mean = exp(r);
}
"""

stan_name = 'CHIP_cov_model' # have to change

stan_file = os.path.join(
    '/Users/tomoyauchiyama/code',
    'STAN',
    stan_name + '.stan'
)

with open(stan_file, 'w') as wf:
    wf.writelines(code)

print(stan_file)

# %%
# stan file ---------------------------------
# -> 10種に共通するパラメータを推定する方法
code = """data {
    int Ni;
    int Nt;
    int Y[Nt, Ni];
}
parameters {
    matrix[Nt, Ni] r;
    real<lower=0> s_r;
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
            Y[i, j] ~ poisson_log(r[i, j]);
}
generated quantities {
    matrix[Nt, Ni] Y_mean;
    for (i in 1:Nt)
        for (j in 1:Ni)
            Y_mean[i, j] = exp(r[i, j]);
}
"""

stan_name = 'CHIP_covCom_model' # have to change

stan_file = os.path.join(
    '/Users/tomoyauchiyama/code',
    'STAN',
    stan_name + '.stan'
)

with open(stan_file, 'w') as wf:
    wf.writelines(code)

print(stan_file)

# %%
#--------------------------------------
# compile stan file
#--------------------------------------
import cmdstanpy as stan

sm = stan.CmdStanModel(stan_file=stan_file)

# %%
print(sm.name)
print(sm.stan_file)
print(sm.exe_file)
print(sm.code())

# %%
#--------------------------------------
# MCMC sampling (NUTS method)
#--------------------------------------

fit = sm.sample(
    data=json_file,
    chains=3,
    iter_sampling=2500,
    iter_warmup=500,
    thin=5,
    seed=1234
)

"""
1番目のデータの実行時間: 102m 17.8s
250番目のデータの実行時間: 100m 44.5s
"""

# %%
#--------------------------------------
# MCMC sampling summary
#--------------------------------------
import pandas as pd
# summaries of
# the total joint log-probability density lp__
# plus all model parameters
# and quantities of interest in a pandas.DataFrame

# type -> pandas (basic statistic)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
print(fit.summary().loc[:, ['Mean', 'StdDev', '5%', '50%', '95%', 'N_Eff', 'R_hat']])

# %%
#--------------------------------------
# MCMC sampling analyzation
#--------------------------------------

# analyze the per-draw sampler parameters across all chains
# looking for potential problems which indicate that
# the sample isn’t a representative sample from the posterior

print(fit.diagnose())

"""
Processing csv files:
/var/folders/15/d1_kk_5j0cj54kx3_k_76tcm0000gn/T/tmpt0fzz5f_/CHIP_covCom_model-202201161824-1-b7k93_on.csv,
/var/folders/15/d1_kk_5j0cj54kx3_k_76tcm0000gn/T/tmpt0fzz5f_/CHIP_covCom_model-202201161824-2-umt4dauu.csv,
/var/folders/15/d1_kk_5j0cj54kx3_k_76tcm0000gn/T/tmpt0fzz5f_/CHIP_covCom_model-202201161824-3-spjzw7jb.csv

Checking sampler transitions treedepth.
1500 of 1500 (1e+02%) transitions hit the maximum treedepth limit of 10, or 2^10 leapfrog steps.
Trajectories that are prematurely terminated due to this limit will result in slow exploration.
For optimal performance, increase this limit.

Checking sampler transitions for divergences.
No divergent transitions found.

Checking E-BFMI - sampler transitions HMC potential energy.
E-BFMI satisfactory.

Effective sample size satisfactory.

Split R-hat values satisfactory all parameters.

Processing complete.

# -> 確認項目に満たした結果で、OK
# 特にR-hat値が1.0に収束していたのでOK
"""

# %%
#--------------------------------------
# parameter visualization
#--------------------------------------
import arviz

arviz.plot_trace(fit)

"""
実行時間：102m 17.8s
"""

# %%
print(fit.summary())

# %%
#--------------------------------------
# parameter data extraction
#--------------------------------------
# extract parmetar 'Y_mean' -----------
parm_Y = fit.stan_variable('Y_mean')

# %%
print(parm_Y.shape)
# -> output: (1500, 10, 1024)

# %%
s_r = fit.stan_variable('s_r')
print(s_r.mean())

# %%
# calc ave of parm_r -------------------
Y_mean = np.median(parm_Y, axis=0)

# realize array shape ------------------
print(Y_mean.shape)
# output: (10, 1024)

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
#--------------------------------------
# statistic modeling visualization
#--------------------------------------

# prepare for data visualization ------
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
    palette = sns.hls_palette(n_colors, l=0.6, s=1)
    #palette = sns.color_palette('hls', n_colors)
    #colorpalette = 'hls'
    rgb = ['rgb({},{},{})'.format(*[x*256 for x in rgb]) for rgb in palette]
    return rgb

# %%
# 3種のChIP-seqデータ(テストデータ)の可視化 ----------------

y_title = 'Covarage'
x_title = 'Position'
main_title = 'CHIP-seq data (*σ(r) = ' + str(round(s_r.mean(), 2)) + ')'

titles = []
n = 10
for i in range(n):
    ID = int(i + 1)
    name = 'CHIP_' + str(ID)
    titles.append(name)

colors = get_colorpalette(n + 1)
probs = (10, 25, 50, 75, 90)
pxx = [f'p{p}' for p in probs]
n_row = 2
n_col = 5

fig = make_subplots(
    rows=n_row,
    cols=n_col,
    subplot_titles=titles,
    horizontal_spacing=0.05, # will change
    vertical_spacing=0.15    # will change
)

ii = -1
with h5py.File(infile, 'r') as h5:
    min_x = 1
    max_x = 1024
    min_y = h5['test_out'][:,:,:].min()
    max_y = h5['test_out'][:,:,:].max()

    for i in range(1, n_row+1):
        for j in range(1, n_col+1):
            ii += 1

            cov = h5['test_out'][249, :, ii]
            y = cov
            x = [i+1 for i in range(len(cov))]
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
                        width=0.5
                    )
                ),
                row=i,
                col=j
            )
            """
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode='markers',
                    marker=dict(
                        color=colors[ii],
                        size=0.9,
                        line=dict(
                            color='black',
                            width=0.4
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
            fig.update_layout(
                plot_bgcolor='white'
                #height=800,
                #width=900
            )
            fig.update_xaxes(
                title=x_title,
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
                title=y_title,
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
                size=10
            )
        )
    )
    fig.for_each_yaxis(
        lambda axis: axis.title.update(
            font=dict(
                color='black',
                size=10
            )
        )
    )
    fig.update_layout(
        title=dict(
            text=main_title,
            x=0.5,
            xanchor='center'
        ),
        showlegend=False
    )
    fig.update_annotations(
        font=dict(size=10)
    )
    fig.show()
    htmlfile = os.path.join(
    '/Users/tomoyauchiyama/code/CNN/',
    'CHIP-seq_data_250.html'
    )
    fig.write_html(htmlfile)

# %%
"""
転写因子結合部位の同定、motif配列の探索
・転写因子結合部位の同定で、どのピークを選択するかを決める
① backgroundを調べる。実測値からbackgroundを引いたcovarageを
　chip-peakと定義するか？
② または、実測値の分布と予測値の散布図を描いて、
　covrage上位5%の大きさ（両方の箱ひげを描く）を持つピークを探索対象とするか？
-> ②をやってみる
"""
# valid_out

y_title = 'Covarage'
x_title = 'Position'
main_title = 'CHIP-seq validation data'

titles = []
n = 9
for i in range(n):
    ID = int(i + 1)
    name = 'CHIP_' + str(ID)
    titles.append(name)

colors = get_colorpalette(n + 1)
probs = (10, 25, 50, 75, 90)
pxx = [f'p{p}' for p in probs]
n_row = 2
n_col = 5

fig = make_subplots(
    rows=n_row,
    cols=n_col,
    subplot_titles=titles,
    horizontal_spacing=0.05, # will change
    vertical_spacing=0.15    # will change
)

ii = -1
with h5py.File(infile, 'r') as h5:
    min_x = 1
    max_x = 1024
    min_y = h5['valid_out'][:,:,:].min()
    max_y = h5['valid_out'][:,:,:].max()

    for i in range(1, n_row+1):
        for j in range(1, n_col+1):
            ii += 1

            cov = h5['valid_out'][:, :, ii].mean(axis=0)
            y = cov
            x = [i+1 for i in range(len(cov))]
            fig.add_trace(
                go.Bar(
                    x=x,
                    y=y,
                    marker_color=colors[ii]
                    #opacity=
                ),
                row=i,
                col=j
            )
            fig.update_layout(
                plot_bgcolor='white'
                #height=800,
                #width=900
            )
            fig.update_xaxes(
                title=x_title,
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
                title=y_title,
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
                size=10
            )
        )
    )
    fig.for_each_yaxis(
        lambda axis: axis.title.update(
            font=dict(
                color='black',
                size=10
            )
        )
    )
    fig.update_layout(
        title=dict(
            text=main_title,
            x=0.5,
            xanchor='center'
        ),
        showlegend=False
    )
    fig.update_annotations(
        font=dict(size=10)
    )
    fig.show()
    htmlfile = os.path.join(
    '/Users/tomoyauchiyama/code/CNN/',
    'CHIP-seq_valdata_1.html'
    )
    fig.write_html(htmlfile)

# %%
"""
転写因子結合部位の同定、motif配列の探索
・転写因子結合部位の同定で、どのピークを選択するかを決める
① backgroundを調べる。実測値からbackgroundを引いたcovarageを
　chip-peakと定義するか？
② または、実測値の分布と予測値の散布図を描いて、
　covrage上位5%の大きさ（両方の箱ひげを描く）を持つピークを探索対象とするか？
-> ②をやってみる
"""
import matplotlib.pyplot as plt

pre = pred_df.loc[0, 'p50']
with h5py.File(infile, 'r') as h5:
    cov = h5['test_out'][249, :, 0]

print(pd.Series(pre).quantile(0.95))
print(pd.Series(cov).quantile(0.95))

plt.scatter(pre, cov, c='k', s=7)
plt.axhline(y=pd.Series(cov).quantile(0.95))
plt.axvline(x=pd.Series(pre).quantile(0.95))
plt.show()

# %%

data = []
for i in range(10):
    with h5py.File(infile, 'r') as h5:
        pre = pred_df.loc[i, 'p50']
        cov = h5['test_out'][249, :, i]

        pre_95 = pd.Series(pre).quantile(0.95)
        cov_95 = pd.Series(cov).quantile(0.95)

        tmp = []
        for j in range(len(pre)):
            if (pre_95 < pre[j]) & (cov_95 < cov[j]):
                tmp.append(j)

        data.append(np.array(tmp, dtype=int))

data = np.array(data)

# %%

out_dir = '/Users/tomoyauchiyama/code/CNN/fragment'
nucl = ['A', 'T', 'C', 'G']
for i in range(10):
    fasta = os.path.join(out_dir, 'CHIP_' + str(i + 1) + '.fa')
    print(fasta)
    with open(fasta, 'w') as wf:
        with h5py.File(infile, 'r') as h5:
            for index in data[i]:
                tmp = []
                id = '>pos_' + str(index + 1)
                start = index * 128
                end = (index + 1) * 128
                seq = h5['test_in'][249, start:end, :]
                for j in range(len(seq)):
                    tmp.append(nucl[seq[j].argmax()])
                fragnemt = ''.join(tmp)
                wf.write(f'{id}\n')
                wf.write(f'{fragnemt}\n')

# %%
import plotly.express as px

# %%
fig = px.colors.sequential.swatches()
fig.show()

# %%
print(px.colors.named_colorscales())

# %%

print(px.colors.sequential.Reds)

# %%
import seaborn as sns

# %%
color_list = [
    'Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r',
    'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r','Dark2', 'Dark2_r', 'GnBu', 'GnBu_r',
    'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r',
    'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r',
    'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r',
    'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy','RdGy_r',
    'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r',
    'Set1', 'Set1_r','Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r',
    'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r',
    'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r',
    'bone', 'bone_r', 'brg', 'brg_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r',
    'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'crest', 'crest_r',
    'cubehelix', 'cubehelix_r', 'flag', 'flag_r', 'flare', 'flare_r', 'gist_earth', 'gist_earth_r',
    'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r','gist_rainbow', 'gist_rainbow_r',
    'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r',
    'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'icefire', 'icefire_r',
    'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'mako', 'mako_r',
    'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r',
    'prism', 'prism_r', 'rainbow', 'rainbow_r', 'rocket', 'rocket_r', 'seismic', 'seismic_r',
    'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r',
    'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'turbo', 'turbo_r',
    'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'viridis', 'viridis_r', 'vlag', 'vlag_r',
    'winter', 'winter_r'
]

# %%
print(sns.color_palette('hls'))

# %%
print(len(color_list))
# %%
pos = ['B1', 'B2', 'B3', 'B2', 'A1', 'A2', 'A3', 'D1', 'D2', 'C1', 'C2']

print(sorted(pos, key=lambda x: (int(x[1:]), x[0])))
# %%
import pandas as pd

# %%
sales_list1=[
    ["B1", "iPhoneX", 85000, 1],
    ["B2", "iPhoneXX", 260000, 2],
    ["B3", "iPhoneKK", 37000, 1],
    ["B2", "iPhoneXX-D", 130000, 1],
    ["A1", "iPhoneCX", 130000, 1],
    ["A2", "iPhoneYY-D", 130000, 1],
    ["A2", "iPhoneYY", 130000, 1],
    ["A3", "iPhoneSX", 130000, 1],
    ["D1", "iPhoneTAX", 130000, 1],
    ["D2", "iPhoneKKK", 130000, 1],
    ["C1", "iPhoneOOP", 130000, 1],
    ["C2", "iPhonePPA", 130000, 1],
]
columns1 = ["ID", "Name", "Amount (JPY)", "Qty"]
df = pd.DataFrame(data=sales_list1, columns=columns1)
print(df)

# %%
# Nameでsort
# -> MassArray WSの座標位置決定の準備
df_s = df.sort_values('Name')
print(df_s)

# %%
# DataFrameで文字と数字を含む文字列の数字でsortする方法
# -> MassArray WSの座標位置決定
df_sort = df_s.sort_values(
    by='ID',
    key=lambda col: col.map(lambda x: (int(x[1:]), x[0]))
)
print(df_sort)

# %%
# 末尾に「-D」を含まないIDとNameの抽出
# -> MassArray WSの座標位置決定後に行い、biomec_samplelistのdest wellの所に貼り付ける。
print(df_sort[~df_sort['Name'].str.endswith('-D')][['ID','Name']])

# %%

import re
import glob
import pandas as pd
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots


# %%

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
tsv_path = '/Users/tomoyauchiyama/code/CNN/DIST'
outdir = '/Users/tomoyauchiyama/code/CNN/DIST'
outfile = 'DIST.html'

titles = []
for tsv in sorted(glob.glob(os.path.join(tsv_path, '*tsv'))):
    name = re.compile(r'DIST/(.*)_dist\.tsv').search(tsv).group(1)
    titles.append(name)

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

main_title = 'Distribution of length'
y_title = 'Count'
x_title = 'Position'

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
for i in range(1, n_row+1):
    for j in range(1, n_col+1):
        ii += 1
        if n-1 < ii:
            pass
        else:
            file = os.path.join(
                tsv_path,
                titles[ii] + '_dist.tsv'
            )
            df = pd.read_table(file)
            x = df['V1']
            y = df['count']
            min_x = x.min()
            max_x = x.max()
            min_y = 0
            max_y = y.max() + 50
            fig.add_trace(
                go.Bar(
                    x=x,
                    y=y,
                    marker=dict(
                        color=colors[ii],
                        line=dict(
                            color='black',
                            width=0.8
                        )
                    )
                ),
                row=i,
                col=j
            )
            fig.update_layout(
                plot_bgcolor='white'
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

file = '/Users/tomoyauchiyama/code/PR2688/oRNAment/tmp/*PWM'
outdir = '/Users/tomoyauchiyama/code/PR2688/oRNAment'
for pwm in glob.glob(file):
    name = os.path.splitext(os.path.basename(pwm))[0]
    outfile = os.path.join(outdir, name + '.PWM')
    with open(outfile, 'w') as wf:
        with open(pwm, 'r') as rf:
            for line in rf.readlines():
                line = line.rstrip('\n')
                if 'PO' in line:
                    wf.write(f'>{name}\n')
                else:
                    tmp = line.split()
                    row = '\t'.join(tmp[1:])
                    wf.write(f'{row}\n')
# %%
file = '/Users/tomoyauchiyama/code/PR2688/oRNAment/*PWM'
outdir = '/Users/tomoyauchiyama/code/PR2688/oRNAment'
cmd = '/usr/local/Cellar/meme/5.1.0/libexec/libexec/meme-5.1.0/chen2meme'
for pwm in glob.glob(file):
    name = os.path.splitext(os.path.basename(pwm))[0]
    outfile = os.path.join(outdir, name + '.meme')
    print(f'{cmd} {pwm} > {outfile}')
    os.system(f'{cmd} {pwm} > {outfile}')
# %%

import statsmodels.stats.multitest as multi
from statsmodels.stats.proportion import proportions_ztest

pval = proportions_ztest([5, 4], [100, 100], alternative='two-sided')
pval_corr = multi.multipletests(pval[1], alpha=0.05, method='fdr_bh')

# %%
pval_corr[1]
# %%
