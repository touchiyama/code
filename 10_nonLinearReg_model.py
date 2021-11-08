# %%

#
# 線形ガウス基底回帰モデル構築
#

import os
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.optimize import minimize
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# %%
#--------------------------------------
# 関数の定義
#--------------------------------------

def model(w, x):
    """線形ガウス基底回帰モデル

    Args:
        w ([type]): [description]
        x ([type]): [description]

    Returns:
        [type]: [description]
    """
    #return w[0]-w[1]*np.exp(-w[2]*x)
    return w[0]*(1-np.exp(-w[1]*x))

def mse_model(w, x, y):
    """平均2乗誤差

    Args:
        w ([type]): [description]
        x ([type]): [description]
        t ([type]): [description]

    Returns:
        [type]: [description]
    """
    m = model(w, x)
    return np.mean((m-y)**2)

def fit_model(w_init, x, t):
    """パラメータ最適化

    Args:
        w_init ([type]): [description]
        x ([type]): [description]
        t ([type]): [description]

    Returns:
        [type]: [description]
    """
    res = minimize(mse_model, w_init, args=(x, t), method="powell")
    return res.x

# %%
#--------------------------------------------
# 欠損値を含まないデータを縦長のデータフレームに加工
#--------------------------------------------

#======  read row data without N.A  ======#
normal_df = pd.read_csv('/Users/tomoyauchiyama/RStanBook/chap08/input/data-conc-2.txt')

print(normal_df)

# %%
#======  convert wide DataFreme to long DataFrame  ======#
"""remove row including N.A
   substitute int type of number to add Time value
   add TimeID column linked to Time value after sorting DataFrame by 'Time'
   sort DataFrame by 'PersonID' and 'TimeID'
   reset index
"""

normal_l_df = pd.melt(normal_df, id_vars='PersonID', ignore_index=True).dropna()
normal_l_df.columns = ('PersonID', 'Time', 'Y')
normal_l_df['Time'] = normal_l_df.loc[:, 'Time'].str.extract(r'Time(\d+)').astype(int)
normal_l_df_sort = normal_l_df.sort_values('Time')
normal_l_df_sort['TimeID'] = normal_l_df_sort['Time'].factorize()[0] + 1
normal_l_df_sort = normal_l_df_sort.sort_values(['PersonID', 'Time']).reset_index()

print(len(normal_l_df_sort))
print(normal_l_df_sort)

# %%
#======  prepare for data visualization ======#
X = {}
Y = {}
for i in range(len(normal_l_df_sort)):
    ID = int(normal_l_df_sort.loc[i, 'PersonID'])
    x_val = str(normal_l_df_sort.loc[i, 'Time'])
    y_val = str(normal_l_df_sort.loc[i, 'Y'])
    if X.get(ID):
        X[ID] += ',' + x_val
        Y[ID] += ',' + y_val
    else:
        X[ID] = x_val
        Y[ID] = y_val

proc_n_df = pd.DataFrame()
for ID in sorted(X.keys()):
    xx = [int(i) for i in X[ID].split(',')]
    yy = [float(i) for i in Y[ID].split(',')]
    col = pd.Series([ID, xx, yy])
    proc_n_df = proc_n_df.append(col, ignore_index=True)

proc_n_df.columns = ['PersonID', 'Time', 'Y']

print(proc_n_df)


# %%
#--------------------------------------
# main関数
#--------------------------------------

x_title = 'Time(hour)'
y_title = 'Conc(mg/ml)'

fig = make_subplots()

max_x = normal_l_df_sort.loc[:, 'Time'].unique().max() + 3
min_x = -0.2
max_y = normal_df.describe().loc['max', :].max() + 3
min_y = 0


y = np.array(proc_n_df.loc[2, 'Y'])
x = np.array(proc_n_df.loc[2, 'Time'])

w_init=[1, 1]
w = fit_model(w_init, x, y)
mse = mse_model(w, x, y)

# %%
pre_x = np.linspace(0, max_x, 60)
pre_y = model(w, pre_x)

fig.add_trace(
    go.Scatter(
        x=pre_x,
        y=pre_y,
        mode='lines',
        line=dict(color=None, width=2)
    )
)
fig.add_trace(
    go.Scatter(
        x=x,
        y=y,
        mode='markers',
        marker=dict(
            color=None,
            size=5,
            line=dict(color='black', width=1)
        )
    )
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
    dtick=5
)
fig.update_yaxes(
    title=y_title,
    showline=True,
    linewidth=1,
    linecolor='black',
    mirror=True,
    ticks='inside',
    range=(min_y, max_y),
    dtick=5
)

# %%
#print('w0={0:.1f}, w1={1:.1f}, w2={2:.1f}'.format(w[0], w[1], w[2]))
print('w0={0:.1f}, w1={1:.1f}'.format(w[0], w[1]))
print('SD={0:.2f} mg/ml'.format(np.sqrt(mse)))

# %%
#--------------------------------------
# 欠損値を含むデータフレームの構築
#--------------------------------------

df_na = pd.read_csv('/Users/tomoyauchiyama/RStanBook/chap09/input/data-conc-2-NA-wide.txt')

print(df_na)
# %%
df_na_long = pd.melt(df_na, id_vars='PersonID', ignore_index=True).dropna()
df_na_long.columns = ('PersonID', 'Time', 'Y')
df_na_long['Time'] = df_na_long.loc[:, 'Time'].str.extract(r'Time(\d+)').astype(int)
df_na_long_sort = df_na_long.sort_values('Time')
df_na_long_sort['TimeID'] = df_na_long_sort['Time'].factorize()[0] + 1
df_na_long_sort = df_na_long_sort.sort_values(['PersonID', 'Time']).reset_index()

print(len(df_na_long_sort))
print(df_na_long_sort)

# %%
#======  prepare for data visualization ======#
X = {}
Y = {}
for i in range(len(df_na_long_sort)):
    ID = int(df_na_long_sort.loc[i, 'PersonID'])
    x_val = str(df_na_long_sort.loc[i, 'Time'])
    y_val = str(df_na_long_sort.loc[i, 'Y'])
    if X.get(ID):
        X[ID] += ',' + x_val
        Y[ID] += ',' + y_val
    else:
        X[ID] = x_val
        Y[ID] = y_val

proc_na_df = pd.DataFrame()
for ID in sorted(X.keys()):
    xx = [int(i) for i in X[ID].split(',')]
    yy = [float(i) for i in Y[ID].split(',')]
    col = pd.Series([ID, xx, yy])
    proc_na_df = proc_na_df.append(col, ignore_index=True)

proc_na_df.columns = ['PersonID', 'Time', 'Y']

print(proc_na_df)

# %%
#--------------------------------------
# main関数
#--------------------------------------

x_title = 'Time(hour)'
y_title = 'Conc(mg/ml)'

max_x = df_na_long_sort.loc[:, 'Time'].unique().max() + 3
min_x = -0.2
max_y = df_na.describe().loc['max', :].max() + 3
min_y = 0


y = np.array(proc_na_df.loc[1, 'Y'])
x = np.array(proc_na_df.loc[1, 'Time'])

w_init=[1, 1]
w = fit_model(w_init, x, y)
mse = mse_model(w, x, y)

# %%
pre_x = np.linspace(0, max_x, 60)
pre_y = model(w, pre_x)

fig = make_subplots()

fig.add_trace(
    go.Scatter(
        x=pre_x,
        y=pre_y,
        mode='lines',
        line=dict(color=None, width=2)
    )
)
fig.add_trace(
    go.Scatter(
        x=x,
        y=y,
        mode='markers',
        marker=dict(
            color=None,
            size=5,
            line=dict(color='black', width=1)
        )
    )
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
    dtick=5
)
fig.update_yaxes(
    title=y_title,
    showline=True,
    linewidth=1,
    linecolor='black',
    mirror=True,
    ticks='inside',
    range=(min_y, max_y),
    dtick=5
)

# %%
#print('w0={0:.1f}, w1={1:.1f}, w2={2:.1f}'.format(w[0], w[1], w[2]))
print('w0={0:.1f}, w1={1:.1f}'.format(w[0], w[1]))
print('SD={0:.2f} mg/ml'.format(np.sqrt(mse)))

# %%
#--------------------------------------
# 教科書の第5章の例題に従うデータセットの作成
#--------------------------------------

np.random.seed(seed=1)            # 乱数を固定
X_min = 4                         # X の下限（表示用）
X_max = 30                        # X の上限（表示用）
X_n = 16                          # データの個数
X = 5 + 25 * np.random.rand(X_n)  # X の生成
Prm_c = [170, 108, 0.2]           # 生成パラメータ
T = Prm_c[0] - Prm_c[1] * np.exp(-Prm_c[2] * X) + 4 * np.random.randn(X_n)
T_max = T.max() + 10
T_min = 0

# %%
print(type(T))

# %%
x_title = 'X'
y_title = 'Y'

fig = make_subplots()

fig.add_trace(
    go.Scatter(
        x=X,
        y=T,
        mode='markers',
        marker=dict(
            color='cornflowerblue',
            size=5,
            line=dict(color='black', width=1)
        )
    )
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
    range=(0, X_max),
    dtick=5
)
fig.update_yaxes(
    title=y_title,
    showline=True,
    linewidth=1,
    linecolor='black',
    mirror=True,
    ticks='inside',
    range=(T_min, T_max),
    dtick=20
)

# %%
#--------------------------------------
# make stan file
#--------------------------------------
code = """data {
    int N;
    vector[N] X;
    vector[N] T;
    int new_X;
    real X_new[new_X];
}
parameters {
    real m_a;
    real m_b;
    real m_c;
    real <lower=120, upper=190> a;
    real <lower=0, upper=120> b;
    real <lower=0, upper=1> c;
    real<lower=0> s_a;
    real<lower=0> s_b;
    real<lower=0> s_c;
    real <lower=0> s_T;
}
transformed parameters {
    vector[N] mu;
    mu = a-b*exp(-c*X);
}
model {
    a ~ normal(m_a, s_a);
    b ~ normal(m_b, s_b);
    c ~ normal(m_c, s_c);
    to_vector(T) ~ normal(to_vector(mu), s_T);
}
generated quantities {
    real T_new[new_X];
    for (i in 1:new_X){
        T_new[i] = normal_rng(a-b*exp(-c*X_new[i]), s_T);
    }
}
"""

stan_name = 'nonLinearReg_model' # have to change

stan_file = os.path.join(os.path.abspath('.'), 'STAN', stan_name + '.stan')

with open(stan_file, 'w') as wf:
    wf.writelines(code)

print(stan_file)

# %%
#--------------------------------------
# make json file
#--------------------------------------
import json

observed_data = dict(
    N=X_n,
    X=X.tolist(),
    T=T.tolist(),
    new_X=60,
    X_new=np.linspace(0, X_max, 60).tolist()
)

json_name = 'randam_data' # have to change

json_file = os.path.join(os.path.abspath('.'), 'STAN', json_name + '.json')

with open(json_file, 'w') as wf:
    json.dump(observed_data, wf, indent=2)

print(json_file)

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

# %%
#--------------------------------------
# MCMC sampling summary
#--------------------------------------

# summaries of
# the total joint log-probability density lp__
# plus all model parameters
# and quantities of interest in a pandas.DataFrame

# type -> pandas (basic statistic)
#pd.set_option('display.max_columns', None)
print(fit.summary().loc[:, ['Mean', 'StdDev', '5%', '50%', '95%', 'N_Eff', 'R_hat']])

# %%
#--------------------------------------
# MCMC sampling analyzation
#--------------------------------------

# analyze the per-draw sampler parameters across all chains
# looking for potential problems which indicate that
# the sample isn’t a representative sample from the posterior

print(fit.diagnose())
# %%
#--------------------------------------
# parameter visualization
#--------------------------------------
import arviz

arviz.plot_trace(fit)

# %%
#--------------------------------------------
# 欠損値を含まない場合のモデリング結果の可視化
#--------------------------------------------
#--------------------------------------
# parameter data extraction
#--------------------------------------

#====== extract parmetar 'new_Y' ======#
parm_T = fit.stan_variable('T_new')

#====== realize array shape  ======#
print(parm_T.shape)

# %%
# %%
probs = (10, 25, 50, 75, 90)
column = [f'p{p}' for p in probs]
X_new = np.linspace(0, X_max, 60)

probs = (10, 25, 50, 75, 90)
columns = [f'p{p}' for p in probs]

d_est = pd.DataFrame(
    np.percentile(parm_T, probs, axis=0).transpose(),
    columns=columns
)
d_est['x'] = pd.Series(X_new)

print(d_est)

# %%

fig = make_subplots()

y = T
x = X
pre_y = d_est[columns[2]]
pre_x = d_est['x']

# fill area between trace0 and trace1

for k in range(2):
    fig.add_trace(
        go.Scatter(
            x=d_est['x'],
            y=d_est[columns[k]],
            fill=None,
            mode='lines',
            line=dict(color='black', width=0.4),
            opacity=0.2*(k+1)
        )
    )
    fig.add_trace(
        go.Scatter(
            x=d_est['x'],
            y=d_est[columns[-k-1]],
            fill='tonexty',
            mode='lines',
            line=dict(color='black', width=0.4),
            opacity=0.2*(k+1)
        )
    )
fig.add_trace(
    go.Scatter(
        x=pre_x,
        y=pre_y,
        mode='lines',
        line=dict(color=None, width=2)
    )
)
fig.add_trace(
    go.Scatter(
        x=x,
        y=y,
        mode='markers',
        marker=dict(
            color='cornflowerblue',
            size=5,
            line=dict(color='black', width=1)
        )
    )
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
    range=(0, X_max),
    dtick=5
)
fig.update_yaxes(
    title=y_title,
    showline=True,
    linewidth=1,
    linecolor='black',
    mirror=True,
    ticks='inside',
    range=(T_min, T_max),
    dtick=20
)
fig.update_layout(showlegend=False)

