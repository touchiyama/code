# %%
import os
import json
import arviz
from matplotlib import colors
import matplotlib
import numpy as np
import pandas as pd
import cmdstanpy as stan
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy.stats.morestats import yeojohnson_normplot

# %%
"""
384 plateの位置の影響を考える

(1) 観測データにさらにノイズを上乗せする。
上乗せするノイズの大きさをどれほどにすると、
位置の影響は局面でなく、平面になるか。
→ ランダム効果項をいじる

(2) 位置の影響がある時、
s_Yに情報のある事前分布を設定し、位置の影響を推定する必要がある。
→ 一様分布から、正規分布に変更、パラメータはどうするか

(3) stanコードはベクトル化に変更
"""
# %%
# read input data -------------------------------------

mesh2D = pd.read_csv(
    '/Users/tomoyauchiyama/RStanBook/chap12/input/data-2Dmesh.txt',
    header=None
)

mesh2D_design = pd.read_csv(
    '/Users/tomoyauchiyama/RStanBook/chap12/input/data-2Dmesh-design.txt',
    header=None
)

# %%
print(mesh2D.describe().loc['std', :].mean())
print(mesh2D.describe().loc['mean', :].mean())

# %%
# make mesh2D heatmap --------------------------------------

y = mesh2D.index.values
x = mesh2D.columns.values

fig = go.Figure(
    data=go.Heatmap(
        colorbar=dict(
            title='2Dmesh_value'
        ),
        z=mesh2D,
        x=x,
        y=y,
        colorscale='Viridis'
    )
)

fig.update_xaxes(title='Plate column')
fig.update_yaxes(title='Plate row')

fig.update_layout(
    title='384 plate (2Dmesh) ',
    xaxis_nticks=36
)

fig.show()

# %%
# make mesh2D_design heatmap -------------------------------

y = mesh2D_design.index.values
x = mesh2D_design.columns.values

fig = go.Figure(
    data=go.Heatmap(
        colorbar=dict(
            title='2Dmesh_value'
        ),
        z=mesh2D_design,
        x=x,
        y=y,
        colorscale='Viridis'
    )
)

fig.update_xaxes(title='Plate column')
fig.update_yaxes(title='Plate row')

fig.update_layout(
    title='384 plate (2Dmesh_design) ',
    xaxis_nticks=36
)

fig.show()


# %%
# mesh2D data processing -------------------------------
d_melt = pd.melt(
    mesh2D.reset_index(),
    id_vars='index',
    ignore_index=True
)

d_melt.columns = ('Plate Row', 'Plate Column', 'Y')


# %%
# 3D Surface Plots with plotly -------------------------------
import plotly.express as px

mesh_size = 1.5
xrange = np.arange(0, mesh2D.index.size, mesh_size)
yrange = np.arange(0, mesh2D.columns.size, mesh_size)

fig = px.scatter_3d(
    d_melt,
    x='Plate Column',
    y='Plate Row',
    z='Y'
)
fig.update_traces(
    marker=dict(
        size=1.0,
        line=dict(width=0.5, color='DarkSlateGrey'),
        color='white'
    ),
    selector=dict(mode='markers')
)
fig.add_traces(
    go.Surface(
        z=mesh2D,
        colorscale='Viridis',
        colorbar=dict(title='2Dmesh_value')
    )
)

fig.update_layout(
    title='display 3D Surface Plots',
    xaxis_nticks=36
)

fig.show()

"""
観測値データが記載されたプレートでは、
ヒートマップ同様に、位置によって数値にムラがある
→ プレートの縁の方は、値が小さい
"""

# %%
# mesh2D_design data processing -------------------------------
dd_melt = pd.melt(
    mesh2D_design.reset_index(),
    id_vars='index',
    ignore_index=True
)

dd_melt.columns = ('Plate Row', 'Plate Column', 'Y')

# %%
# 3D Surface Plots with plotly --------------------------------

mesh_size = 1.5
xrange = np.arange(0, mesh2D_design.index.size, mesh_size)
yrange = np.arange(0, mesh2D_design.columns.size, mesh_size)

fig = px.scatter_3d(
    dd_melt,
    x='Plate Column',
    y='Plate Row',
    z='Y'
)
fig.update_traces(
    marker=dict(
        size=1.0,
        line=dict(width=0.5, color='DarkSlateGrey'),
        color='white'
    ),
    selector=dict(mode='markers')
)
fig.add_traces(
    go.Surface(
        z=mesh2D_design,
        #colorscale='Viridis',
        colorbar=dict(title='2Dmesh_value')
    )
)

fig.update_layout(
    title='display 3D Surface Plots(2Dmesh_design)',
    xaxis_nticks=36
)

fig.show()

"""
処理を加えたプレートでは
→ 位置によらずランダムに数値が割り振られていることから、
　ランダムに処理されていることが分かる
"""


# %%
# prepere for sampling -------------------------------
# smoothing

from loess.loess_2d import loess_2d

x = d_melt['Plate Row'].astype('int64').to_numpy()
y = d_melt['Plate Column'].astype('int64').to_numpy()
z = d_melt['Y'].to_numpy()

loess_res, _ = loess_2d(x, y, z, frac=0.2)

i = mesh2D.index.size
j = mesh2D.columns.size
int_r = loess_res.reshape(i, j)

print(int_r)

# %%
# modeling -------------------------------------------

#--------------------------------------
# make stan file
#--------------------------------------
code = """data {
    int I;
    int J;
    matrix[I, J] Y;
    int T;
    int<lower=1, upper=T> TID[I, J];
}
parameters {
    matrix[I, J] r;
    real<lower=0> s_r;
    vector[T] beta;
    real<lower=0> s_beta;
    real<lower=0> s_Y;
}
model {
    target += normal_lpdf(
        to_vector(r[1:I, 3:J]) | to_vector(2*r[1:I, 2:(J-1)] - r[1:I, 1:(J-2)]),
        s_r
    );
    target += normal_lpdf(
        to_vector(r[3:I, 1:J]) | to_vector(2*r[2:(I-1), 1:J] - r[1:(I-2), 1:J]),
        s_r
    );
    beta ~ student_t(6, 0, s_beta);
    to_vector(Y') ~ normal(to_vector(r') + beta[to_array_1d(TID)], s_Y);
}
"""

stan_name = 'mesh2D_model_2' # have to change

stan_file = os.path.join(
    os.path.abspath('.'),
    'STAN',
    stan_name + '.stan'
)

with open(stan_file, 'w') as wf:
    wf.writelines(code)

print(stan_file)

# %%
#--------------------------------------
# make json file
#--------------------------------------
i = mesh2D.index.size
j = mesh2D.columns.size

observed_data = dict(
    I=i,
    J=j,
    Y=mesh2D.to_numpy().reshape(i, j).tolist(),
    T=int(mesh2D_design.values.max()),
    TID=mesh2D_design.to_numpy().reshape(i, j).tolist()
)

json_name = 'mesh2D_data' # have to change

json_file = os.path.join(
    os.path.abspath('.'),
    'STAN',
    json_name + '.json'
)

with open(json_file, 'w') as wf:
    json.dump(observed_data, wf, indent=2)

print(json_file)


# %%
#--------------------------------------
# compile stan file
#--------------------------------------

sm = stan.CmdStanModel(stan_file=stan_file)

# %%
print(sm.name)
print(sm.stan_file)
print(sm.exe_file)
print(sm.code())

# %%
#--------------------------------------
# initilazation
#--------------------------------------

#print(mesh2D.describe().loc['std', :].mean())
#print(mesh2D.describe().loc['mean', :].mean())

init = dict(
    r=int_r,
    s_r=1,
    beta=np.random.normal(0, 1, mesh2D_design.values.max()),
    s_beta=1,
    s_Y=1
)

print(init)

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
    seed=1234,
    inits=init
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

# %%
#--------------------------------------
# parameter visualization
#--------------------------------------
arviz.plot_trace(fit)

# %%
#--------------------------------------
# parameter data extraction
#--------------------------------------

#====== extract parmetar 'r' ======#
parm_r = fit.stan_variable('r')

#====== realize array shape  ======#
print(parm_r.shape)
# output: (800, 16, 24)
# -> (16x24 array) x 800

#====== realize variant type ======#
print(type(parm_r))
# output: <class 'numpy.ndarray'>

#====== calc ave of parm_r   ======#
r_mean = np.median(parm_r, axis=0)
print(r_mean)

#====== realize array shape  ======#
print(r_mean.shape)
# output: (16, 24)
# -> (16x24 array) x 1


# %%
# plot 3D graph with plotly  ======#


I, J = r_mean.shape
mesh_size = 1.5
xrange = np.arange(0, I, mesh_size)
yrange = np.arange(0, J, mesh_size)

fig = px.scatter_3d(d_melt, x='Plate Column', y='Plate Row', z='Y')
fig.update_traces(
    marker=dict(
        size=1.0,
        line=dict(width=0.5, color='DarkSlateGrey'),
        color='white'
    ),
    selector=dict(mode='markers')
)
fig.add_traces(
    go.Surface(
        z=r_mean,
        colorscale='Viridis',
        colorbar=dict(title='pred_Y_value')
    )
)
fig.update_layout(
    title='statistic modeling (2D_mesh)',
    xaxis_nticks=36
)

fig.show()

# %%
"""
処理ごとの数値をまとめる
→ 処理の違いによってばらつきがある
→ データ分布を確認した上で、乱数を生成する
"""

# %%
# mesh2D data processing -------------------------------
data_2D = pd.melt(
    mesh2D.reset_index(),
    id_vars='index',
    ignore_index=True
)

data_2D.columns = ('Plate Row', 'Plate Column', 'Y')

# %%
# mesh2D data processing -------------------------------
data_2DD = pd.melt(
    mesh2D_design.reset_index(),
    id_vars='index',
    ignore_index=True
)

data_2DD.columns = ('Plate Row', 'Plate Column', 'design_Y')

# %%
data_2D['design_Y'] = data_2DD['design_Y']

# %%
print(data_2D)
# %%
design_data = {}
for i in range(len(data_2D)):
    ID = data_2D.loc[i, 'design_Y']
    if design_data.get(ID):
        design_data[ID] += ',' + str(data_2D.loc[i, 'Y'])
    else:
        design_data[ID] = str(data_2D.loc[i, 'Y'])

print(design_data)
# %%
probs = (10, 25, 50, 75, 90)
column = pd.Series([f'p{p}' for p in probs])

proc_2D_df = pd.DataFrame()
for ID in sorted(design_data.keys()):
    y = [float(i) for i in design_data[ID].split(',')]
    df1 = pd.Series([ID, y])
    prob_y =np.percentile(np.array(y), probs)
    df2 = pd.Series(prob_y)
    col = pd.concat([df1, df2])
    proc_2D_df = proc_2D_df.append(col, ignore_index=True)

proc_2D_df.columns = pd.concat([pd.Series(['ID', 'Y']), column])

print(proc_2D_df)

# %%
# %%
"""観測値が正規分布に従うと仮定して乱数を作成
　　
"""
#scale = mesh2D.describe().loc['std', :].mean()

scale = 0.5
z = np.random.normal(
    loc=proc_2D_df['p50'].to_numpy(),
    scale=scale
)
print(len(z))
print(z)

# %%

fig = make_subplots()

X_min = proc_2D_df['ID'].min()
X_max = proc_2D_df['ID'].max()
y_min = data_2D['Y'].min() - 2
y_max = data_2D['Y'].max() + 2

fig.add_trace(
    go.Scatter(
        x=proc_2D_df['ID'],
        y=proc_2D_df['p50'],
        mode='markers',
        marker=dict(
            color='cornflowerblue',
            size=5,
            line=dict(color='black', width=1)
        ),
        error_y = dict(
            type='percent',
            thickness=0.5,
            color='black',
            visible=True
        )
    )
)
fig.add_trace(
    go.Scatter(
        x=proc_2D_df['ID'],
        y=z,
        mode='lines+markers',
        line=dict(color='black', width=2),
        marker=dict(
            color='pink',
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
    title='ID',
    showline=True,
    linewidth=1,
    linecolor='black',
    mirror=True,
    ticks='inside',
    range=(X_min-1, X_max+2),
    dtick=4
)
fig.update_yaxes(
    title='Y',
    showline=True,
    linewidth=1,
    linecolor='black',
    mirror=True,
    ticks='inside',
    range=(y_min, y_max)
    #dtick=20
)
fig.update_layout(showlegend=False)


# %%
"""
元の観測データにノイズを加える
この時、正規分布に従うと仮定
"""
# modeling -------------------------------------------

#--------------------------------------
# make stan file
#--------------------------------------
code = """data {
    int I;
    int J;
    matrix[I, J] Y;
    int T;
    int<lower=1, upper=T> TID[I, J];
}
parameters {
    matrix[I, J] r;
    real<lower=0> s_r;
    vector[T] beta;
    real<lower=0> s_beta;
}
model {
    target += normal_lpdf(
        to_vector(r[1:I, 3:J]) | to_vector(2*r[1:I, 2:(J-1)] - r[1:I, 1:(J-2)]),
        s_r
    );
    target += normal_lpdf(
        to_vector(r[3:I, 1:J]) | to_vector(2*r[2:(I-1), 1:J] - r[1:(I-2), 1:J]),
        s_r
    );
    beta ~ student_t(6, 0, s_beta);
    to_vector(Y') ~ normal(to_vector(r') + beta[to_array_1d(TID)], 0.1);
}
"""

stan_name = 'mesh2D_model_2' # have to change

stan_file = os.path.join(
    os.path.abspath('.'),
    'STAN',
    stan_name + '.stan'
)

with open(stan_file, 'w') as wf:
    wf.writelines(code)

print(stan_file)

# %%
#--------------------------------------
# make json file
#--------------------------------------

i = mesh2D.index.size
j = mesh2D.columns.size

np.random.seed(1234)

scale = 0.1
Y = np.random.normal(
    loc=d_melt['Y'].to_numpy(),
    scale=scale
)

observed_data = dict(
    I=i,
    J=j,
    Y=Y.reshape(i, j).tolist(),
    T=int(mesh2D_design.values.max()),
    TID=mesh2D_design.to_numpy().reshape(i, j).tolist()
)

json_name = 'mesh2D_data_2' # have to change

json_file = os.path.join(
    os.path.abspath('.'),
    'STAN',
    json_name + '.json'
)

with open(json_file, 'w') as wf:
    json.dump(observed_data, wf, indent=2)

print(json_file)

# %%
# prepere for sampling -------------------------------
# smoothing

from loess.loess_2d import loess_2d

x = d_melt['Plate Row'].astype('int64').to_numpy()
y = d_melt['Plate Column'].astype('int64').to_numpy()
z = Y

loess_res, _ = loess_2d(x, y, z, frac=0.2)

i = mesh2D.index.size
j = mesh2D.columns.size
int_r = loess_res.reshape(i, j)

print(int_r)


# %%
#--------------------------------------
# compile stan file
#--------------------------------------

sm = stan.CmdStanModel(stan_file=stan_file)

# %%
print(sm.name)
print(sm.stan_file)
print(sm.exe_file)
print(sm.code())

# %%
#--------------------------------------
# initilazation
#--------------------------------------

#print(mesh2D.describe().loc['std', :].mean())
#print(mesh2D.describe().loc['mean', :].mean())

init = dict(
    r=int_r,
    s_r=1,
    beta=np.random.normal(0, 0.1, mesh2D_design.values.max()),
    s_beta=1
)

print(init)

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
    seed=1234,
    inits=init
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
pd.set_option('display.max_rows', None)
print(fit.summary().loc[:, ['Mean', 'StdDev', '5%', '50%', '95%', 'N_Eff', 'R_hat']])

# %%
#--------------------------------------
# parameter data extraction
#--------------------------------------

#====== extract parmetar 'r' ======#
parm_r = fit.stan_variable('r')

#====== realize array shape  ======#
print(parm_r.shape)
# output: (800, 16, 24)
# -> (16x24 array) x 800

#====== realize variant type ======#
print(type(parm_r))
# output: <class 'numpy.ndarray'>

#====== calc ave of parm_r   ======#
r_mean = np.median(parm_r, axis=0)
print(r_mean)

#====== realize array shape  ======#
print(r_mean.shape)
# output: (16, 24)
# -> (16x24 array) x 1

# %%


# %%
# plot 3D graph with plotly  ======#

I, J = r_mean.shape
mesh_size = 1.5
xrange = np.arange(0, I, mesh_size)
yrange = np.arange(0, J, mesh_size)

fig = px.scatter_3d(
    d_melt,
    x='Plate Column',
    y='Plate Row',
    z='Y'
)
fig.update_traces(
    marker=dict(
        size=1.0,
        line=dict(width=0.5, color='DarkSlateGrey'),
        color='white'
    ),
    selector=dict(mode='markers')
)
fig.add_traces(
    go.Surface(
        z=r_mean,
        colorscale='Viridis',
        colorbar=dict(title='pred_Y_value')
    )
)
fig.update_layout(
    title='statistic modeling (2D_mesh)',
    xaxis_nticks=36
)

fig.show()

# %%
d_melt['Y_rand'] = pd.Series(Y)

# %%
print(d_melt)
# %%
# 3D Surface Plots with plotly --------------------------------

mesh_size = 1.5
xrange = np.arange(0, mesh2D_design.index.size, mesh_size)
yrange = np.arange(0, mesh2D_design.columns.size, mesh_size)

fig = px.scatter_3d(
    d_melt,
    x='Plate Column',
    y='Plate Row',
    z='Y_rand'
)
fig.update_traces(
    marker=dict(
        size=1.0,
        line=dict(width=0.5, color='DarkSlateGrey'),
        color='white'
    ),
    selector=dict(mode='markers')
)
fig.add_traces(
    go.Surface(
        z=mesh2D,
        #colorscale='Viridis',
        colorbar=dict(title='2Dmesh_value')
    )
)

fig.update_layout(
    title='display 3D Surface Plots(2Dmesh)',
    xaxis_nticks=36
)

fig.show()

# %%
print(d_melt)

# %%
# -------------------------------------
# 解析解によるモデリング
# -------------------------------------

# 解析解--------------------------------
def fit_plane(x0, x1, t):
    c_tx0 = np.mean(t * x0) - np.mean(t) * np.mean(x0)
    c_tx1 = np.mean(t * x1) - np.mean(t) * np.mean(x1)
    c_x0x1 = np.mean(x0 * x1) - np.mean(x0) * np.mean(x1)
    v_x0 = np.var(x0)
    v_x1 = np.var(x1)
    w0 = (c_tx1 * c_x0x1 - v_x1 * c_tx0) / (c_x0x1**2 - v_x0 * v_x1)
    w1 = (c_tx0 * c_x0x1 - v_x0 * c_tx1) / (c_x0x1**2 - v_x0 * v_x1)
    w2 = -w0 * np.mean(x0) - w1 * np.mean(x1) + np.mean(t)
    return np.array([w0, w1, w2])

# %%

x = d_melt['Plate Row'].astype('int64').to_numpy()
y = d_melt['Plate Column'].astype('int64').to_numpy()
z = d_melt['Y'].to_numpy()
w = fit_plane(x, y, z)

# %%
def mse_plane(x0, x1, t, w):
    y = w[0] * x0 + w[1] * x1 + w[2]
    mse = np.mean((y - t)**2)
    return mse

mse = mse_plane(x, y, z, w)
# %%
print(mse)
# %%
import matplotlib.pyplot as plt

def show_data2(ax, x0, x1, t):
    ax.plot(
        x0,
        x1,
        t,
        'o',
        color='cornflowerblue',
        markeredgecolor='black',
        markersize=6
    )
    ax.view_init(elev=35, azim=-75)

def show_plane(ax, w):
    px0 = np.linspace(0, 16, 5)
    px1 = np.linspace(0, 24, 5)
    px0, px1 = np.meshgrid(px0, px1)
    y = w[0]*px0 + w[1] * px1 + w[2]
    ax.plot_surface(
        px0,
        px1,
        y,
        rstride=1,
        stride=1,
        alpha=0.3,
        color='blue',
        edgecolor='black'
    )

plt.figure(figsize=(6, 5))
ax = plt.subplot(1,1,1,projection='3d')
show_plane(ax, w)
show_data2(ax, x, y, z)
plt.show()

# %%
# ライブラリー結果のモデリング --------------------------------

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

# %%
print(data)

# %%
# -------------------------------------
# 解析解によるモデリング
# -------------------------------------

# 解析解--------------------------------
def fit_plane(x0, x1, t):
    c_tx0 = np.mean(t * x0) - np.mean(t) * np.mean(x0)
    c_tx1 = np.mean(t * x1) - np.mean(t) * np.mean(x1)
    c_x0x1 = np.mean(x0 * x1) - np.mean(x0) * np.mean(x1)
    v_x0 = np.var(x0)
    v_x1 = np.var(x1)
    w0 = (c_tx1 * c_x0x1 - v_x1 * c_tx0) / (c_x0x1**2 - v_x0 * v_x1)
    w1 = (c_tx0 * c_x0x1 - v_x0 * c_tx1) / (c_x0x1**2 - v_x0 * v_x1)
    w2 = -w0 * np.mean(x0) - w1 * np.mean(x1) + np.mean(t)
    return np.array([w0, w1, w2])

# %%
# 変数--------------------------------

X = df1['input'].astype('int64').to_numpy()
Y = df1['frag_total'].astype('int64').to_numpy()
Z = df1['d_Depth'].to_numpy()
w = fit_plane(X, Y, Z)

print(w)

# %%
def mse_plane(x0, x1, t, w):
    y = w[0] * x0 + w[1] * x1 + w[2]
    mse = np.mean((y - t)**2)
    return mse

mse = mse_plane(X, Y, Z, w)

print(mse)
print("SD={0:.2f}".format(np.sqrt(mse)))

# %%
import matplotlib.pyplot as plt

def show_data2(ax, x0, x1, t):
    ax.plot(
        x0,
        x1,
        t,
        'o',
        color='cornflowerblue',
        markeredgecolor='black',
        markersize=6
    )
    ax.view_init(elev=35, azim=-75)

def show_plane(ax, w):
    px0 = np.linspace(0, max(X), 5)
    px1 = np.linspace(0, max(Y), 5)
    px0, px1 = np.meshgrid(px0, px1)
    y = w[0]*px0 + w[1] * px1 + w[2]
    ax.plot_surface(
        px0,
        px1,
        y,
        rstride=1,
        cstride=1,
        alpha=0.3,
        color='blue',
        edgecolor='black'
    )

plt.figure(figsize=(6, 5))
ax = plt.subplot(1,1,1,projection='3d')
show_plane(ax, w)
show_data2(ax, X, Y, Z)
plt.show()

# %%
# 3D Surface Plots with plotly -------------------------------
import plotly.express as px

mesh_size = 1.5
#xrange = np.arange(0, mesh2D.index.size, mesh_size)
#yrange = np.arange(0, mesh2D.columns.size, mesh_size)

fig = px.scatter_3d(
    data,
    x='frag_total',
    y='input',
    z='d_Depth',
    color='batch'
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
        size=3.0,
        line=dict(width=1.0, color='DarkSlateGrey')
        #color='white'
    ),
    selector=dict(mode='markers')
)
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

htmlfile = os.path.join(
    os.path.abspath('.'),
    '3D_Surface_Plots.html'
)
fig.write_html(htmlfile)

# %%

file = '/Users/tomoyauchiyama/code/PC0174_PC0131_LIB_summmary_211110.xlsx'
lib_df = pd.read_excel(file).dropna(subset=['Depth_80Gb'])

print(lib_df)
# %%

df4 = lib_df[lib_df['batch']=='PC0174_G'].loc[:, ['batch', 'name', 'input', 'frag_total', 'adapter', 'Depth_40Gb', 'Depth_80Gb', 'd_Depth']]
print(df4)

# %%
df4['d_Depth'] = lib_df['Depth_80Gb'] - lib_df['Depth_40Gb']

# %%
mesh_size = 1.5
#xrange = np.arange(0, mesh2D.index.size, mesh_size)
#yrange = np.arange(0, mesh2D.columns.size, mesh_size)

fig = px.scatter_3d(
    df4,
    x='frag_total',
    y='adapter',
    z='d_Depth',
    color='batch'
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
        size=3.0,
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
"""

fig.update_layout(
    title='display 3D Surface Plots',
    xaxis_nticks=36
)
fig.update_layout(showlegend=True)

"""
f = str(round(w[0], 2)) + 'x+' + str(round(w[1], 2)) + 'y' + str(round(w[2], 2))
text = '*f(x,y) =' + f
fig.add_annotation(
    x=-0.1,
    y=0.05,
    text=text,
    font=dict(size=8),
    showarrow=False,
    arrowhead=1,
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

# %%

"""
apapter濃度とフラグメント数の関係（tumorかつqPCRで解析された検体を対象）
(1) input量が25と100の時の比較を行う (ある程度サンプル数が多いから)
(2) input量が100の時におけるサンプル間の比較を行う
この時、傾きは同じなのだろうか？
"""
import numpy as np

af = lib_df[lib_df['qubit_qpcr']=='qPCR'].loc[:, ['batch', 'name', 'input', 'frag_total', 'adapter', 'Depth_40Gb', 'Depth_80Gb', 'd_Depth']]
print(af)

# %%
# ---------------------------------------------
# (1) input量が25と100の時の比較を行う
# ---------------------------------------------

# input量が25ngの時 ----------------------------

a25 = af[af['input']==25].loc[:, 'adapter']
f25 = af[af['input']==25].loc[:, 'frag_total']

plt.scatter(a25, f25)
plt.xlim(0, 900)
plt.ylim(0, 25)
plt.xlabel('adapter')
plt.ylabel('frag_total')

"""全体的にバラついており、安定していないと考えられる
"""

# %%
# ---------------------------------------------
# (1) input量が25と100の時の比較を行う
# ----------------------------------------------

# input量が100ngの時 ----------------------------
a100 = af[af['input']==100].loc[:, 'adapter']
f100 = af[af['input']==100].loc[:, 'frag_total']

plt.scatter(a100, f100)
#plt.scatter(np.log10(a100), np.log10(f100))
plt.xlim(0, 250)
plt.ylim(0, 200)
plt.xlabel('adapter')
plt.ylabel('frag_total')

"""非線形？な概形をしている、2次の曲線が描けそう
"""

# %%
# ---------------------------------------------
# (1) input量が100の時におけるサンプル間の比較を行う
# ---------------------------------------------

af100 = af[af['input']==100]
print(af100)

# %%
# PC0174_Gの時 ----------------------------

a_G = af100[af100['batch']=='PC0174_G'].loc[:, 'adapter']
f_G = af100[af100['batch']=='PC0174_G'].loc[:, 'frag_total']

# %%
# PC0174_ABCの時 ----------------------------
A_100 = af100[af100['batch']=='PC0174_A']
B_100 = af100[af100['batch']=='PC0174_B']
C_100 = af100[af100['batch']=='PC0174_C']

a_A = A_100.loc[:, 'adapter']
f_A = A_100.loc[:, 'frag_total']
a_B = B_100.loc[:, 'adapter']
f_B = B_100.loc[:, 'frag_total']
a_C = C_100.loc[:, 'adapter']
f_C = C_100.loc[:, 'frag_total']

# %%
# PC0131の時 ----------------------------
PC0131 = af100[af100['batch']=='PC0131_adapter']
a_PC0131 = PC0131.loc[:, 'adapter']
f_PC0131 = PC0131.loc[:, 'frag_total']

# %%

plt.scatter(a_G, f_G, color='gray', label='PC0174_G')
plt.scatter(a_A, f_A, color='green', label='PC0174_A')
plt.scatter(a_B, f_B, color='purple', label='PC0174_B')
plt.scatter(a_C, f_C, color='orange', label='PC0174_C')
plt.scatter(a_PC0131, f_PC0131, color='red', label='PC0131_adapter')
plt.xlim(0, 250)
plt.ylim(0, 200)
plt.xlabel('adapter')
plt.ylabel('frag_total')
plt.legend(
    bbox_to_anchor=(1.05, 1),
    loc='upper left',
    borderaxespad=0,
    fontsize=10
)
plt.show()

# %%
mesh_size = 1.5
#xrange = np.arange(0, mesh2D.index.size, mesh_size)
#yrange = np.arange(0, mesh2D.columns.size, mesh_size)

fig = px.scatter_3d(
    af100,
    x='frag_total',
    y='adapter',
    z='d_Depth',
    color='batch',
    range_x=(0, 250),
    range_y=(0, 250),
    range_z=(0, 250)
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
        size=3.0,
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
"""
fig.update_layout(
    title='display 3D Surface Plots',
    xaxis_nticks=36
)
fig.update_layout(showlegend=True)

"""
f = str(round(w[0], 2)) + 'x+' + str(round(w[1], 2)) + 'y' + str(round(w[2], 2))
text = '*f(x,y) =' + f
fig.add_annotation(
    x=-0.1,
    y=0.05,
    text=text,
    font=dict(size=8),
    showarrow=False,
    arrowhead=1,
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

"""adapterの結合ができていないのか？
　　赤点と青点になぜ差が生じているのか？
　　非階層クラスタリングで分類できるか
"""

# %%
#---------------------------------------------------
# scikit-learnを用いて線形モデルと非線形特徴量について考える
#---------------------------------------------------

# sklearn test -------------------------------------
from sklearn.datasets import make_blobs

# 4種のクラスターを作成
X, y = make_blobs(centers=4, random_state=8)

# %%
plt.scatter(
    X[:,0],
    X[:,1],
    c=y,
    cmap='rainbow',
    marker='o',
    edgecolor='gray'
)
plt.xlabel('Feature 1')
plt.ylabel('Feature 2')
plt.show()

# %%
print(y)

# %%
import mglearn

# 2種のクラスター結果を可視化
y = y % 2
mglearn.discrete_scatter(X[:, 0], X[:, 1], y)

plt.xlabel('Feature 1')
plt.ylabel('Feature 2')
plt.show()

# %%
# Xに2番目の特徴量の要素を２乗したものを追加
X_new = np.hstack([X, X[:, 1:] ** 2])

# X_newにyを追加
X_new = np.insert(X_new, 3, y, axis=1)

# %%
fig = px.scatter_3d(
    X_new,
    x=X_new[:, 0],
    y=X_new[:, 1],
    z=X_new[:, 2],
    color=X_new[:, 3],
    range_x=(X_new[:, 0].min()-5, X_new[:, 0].max()+5),
    range_y=(X_new[:, 1].min()-5, X_new[:, 1].max()+5),
    range_z=(-20, X_new[:,2].max()+5)
)
fig.update_traces(
    marker=dict(
        size=3.0,
        line=dict(width=1.0, color='DarkSlateGrey')
        #color='white'
    ),
    selector=dict(mode='markers')
)
fig.update_layout(showlegend=False)
fig.show()

# %%
from sklearn.svm import LinearSVC

X_new = np.hstack([X, X[:, 1:] ** 2])

#LinearSVCのサンプリング試行回数はdefaultで1000回
linear_svm_3d = LinearSVC(max_iter=10000).fit(X_new, y)

# %%
coef = linear_svm_3d.coef_.ravel()
intercept = linear_svm_3d.intercept_

print(f'coef: {coef}')
print(f'intercept: {intercept}')
print(f'class: {linear_svm_3d.classes_}')

# %%
mesh_size = 1.5
#xrange = np.arange(0, mesh2D.index.size, mesh_size)
#yrange = np.arange(0, mesh2D.columns.size, mesh_size)

# X_newにyを追加
X_new = np.insert(X_new, 3, y, axis=1)

fig = px.scatter_3d(
    X_new,
    x=X_new[:, 0],
    y=X_new[:, 1],
    z=X_new[:, 2],
    color=X_new[:, 3],
    range_x=(X_new[:, 0].min()-5, X_new[:, 0].max()+5),
    range_y=(X_new[:, 1].min()-5, X_new[:, 1].max()+5),
    range_z=(-20, X_new[:,2].max()+5)
)
fig.update_traces(
    marker=dict(
        size=2.0,
        line=dict(width=1.0, color='DarkSlateGrey')
        #color='white'
    ),
    selector=dict(mode='markers')
)

fig.update_layout(
    legend=dict(
        x=-0.1,
        xanchor='left',
        y=1,
        yanchor = 'auto'
    )
)
fig.update_layout(showlegend=False)

px0 = np.linspace(X_new[:, 0].min()-2, X_new[:, 0].max()+2, 50)
px1 = np.linspace(X_new[:, 1].min()-2, X_new[:, 1].max()+2, 50)
px0, px1 = np.meshgrid(px0, px1)
z = -(px0 * coef[0] + px1 * coef[1] + intercept) / coef[2]

fig.add_traces(
    go.Surface(
        x=px0,
        y=px1,
        z=z
        #colorscale='Viridis'
        #colorbar=dict(title='d_depth')
    )
)

fig.update_layout(
    title='display 3D Surface Plots',
    xaxis_nticks=36
)
fig.update_layout(showlegend=True)

"""
f = str(round(w[0], 2)) + 'x+' + str(round(w[1], 2)) + 'y' + str(round(w[2], 2))
text = '*f(x,y) =' + f
fig.add_annotation(
    x=-0.1,
    y=0.05,
    text=text,
    font=dict(size=8),
    showarrow=False,
    arrowhead=1,
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
# %%

ZZ = px1 ** 2

#dec = linear_svm_3d.decision_function(
#    np.c_[px0.ravel(), px1.ravel(), ZZ.ravel()]
#)
dec = linear_svm_3d.predict(
    np.c_[px0.ravel(), px1.ravel(), ZZ.ravel()]
)
plt.contourf(
    px0,
    px1,
    dec.reshape(px0.shape),
    labels=[dec.min(), 0, dec.max()],
    cmap=mglearn.cm2,
    alpha=0.5
)
mglearn.discrete_scatter(X[:, 0], X[:, 1], y)

# %%
X = af100.loc[:, ['frag_total', 'adapter', 'd_Depth']].to_numpy()

# %%
# kmeansを実装 -------------------------------------
from sklearn.cluster import KMeans

y_KMeans = KMeans(n_clusters=3, random_state=0).fit_predict(X)

# %%
print(y_KMeans)
print(type(y_KMeans))

af100_kmeans = af100.reset_index(drop=True)
af100_kmeans['kmeans'] = pd.Series(y_KMeans)

print(af100_kmeans)

# %%
# 凝集型クラスタリングを実装 -------------------------------------
from sklearn.cluster import AgglomerativeClustering

y_agg = AgglomerativeClustering(n_clusters=3).fit_predict(X)

# %%
#LinearSVCのサンプリング試行回数はdefaultで1000回
linear_svm_3d = LinearSVC(max_iter=10000).fit(X, y_agg)

# %%
coef = linear_svm_3d.coef_
intercept = linear_svm_3d.intercept_

print(f'coef: {coef}')
print(f'intercept: {intercept}')
print(f'class: {linear_svm_3d.classes_}')

# %%
print(y_agg)
print(type(y_agg))

af100_agg = af100.reset_index(drop=True)
af100_agg['agg'] = pd.Series(y_agg)

print(af100_agg)

# %%
mesh_size = 1.5
#xrange = np.arange(0, mesh2D.index.size, mesh_size)
#yrange = np.arange(0, mesh2D.columns.size, mesh_size)

fig = px.scatter_3d(
    af100_agg,
    x='frag_total',
    y='adapter',
    z='d_Depth',
    color='agg',
    symbol='agg',
    range_x=(0, 250),
    range_y=(0, 250),
    range_z=(0, 250)
)
fig.update_traces(marker_coloraxis=None)
fig.update_layout(
    showlegend=True,
    legend=dict(
        x=-0.1,
        xanchor='left',
        y=1,
        yanchor='auto'
    )
)
fig.update_traces(
    marker=dict(
        size=3.0,
        line=dict(width=1.0, color='DarkSlateGrey')
        #color='white'
    ),
    selector=dict(mode='markers')
)
"""
px0 = np.linspace(0, af100_agg['frag_total'].max(), 100)
px1 = np.linspace(0, af100_agg['adapter'].max(), 100)
px0, px1 = np.meshgrid(px0, px1)
z = -(px0 * coef[0] + px1 * coef[1] + intercept) / coef[2]

fig.add_traces(
    go.Surface(
        x=px0,
        y=px1,
        z=z,
        colorscale='Viridis',
        colorbar=dict(title='d_depth')
    )
)
"""
fig.update_layout(
    title='display 3D Surface Plots (AgglomerativeClustering)',
    xaxis_nticks=36
)

"""
f = str(round(w[0], 2)) + 'x+' + str(round(w[1], 2)) + 'y' + str(round(w[2], 2))
text = '*f(x,y) =' + f
fig.add_annotation(
    x=-0.1,
    y=0.05,
    text=text,
    font=dict(size=8),
    showarrow=False,
    arrowhead=1,
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

# %%
#階層的クラスタリングはscipy
import scipy.spatial.distance as distance
from scipy.cluster.hierarchy import dendrogram, ward

#ward法で分類
linkage_array = ward(X)

ax = plt.figure(figsize=(20,10)).gca()
dendrogram(linkage_array)
bounds = ax.get_xbound()

plt.xlabel("sample index",fontsize=10)
plt.ylabel("Cluster distance",fontsize=10)

# %%
"""凝集型クラスタリングを用いると説明が付きそう
　　->　凝集型クラスタリングで分類されたクラスを1-of-K 符号化法で表現して
　　　　ロジスティック回帰を用いたクラス分類モデルを実装してみる
"""

# %%
# クラスタ数が3の時 --------------------------------------------
y_agg = AgglomerativeClustering(n_clusters=3).fit_predict(X)

# %%
linear_svm_3d = LinearSVC(max_iter=10000).fit(X, y_agg)

# %%
coef = linear_svm_3d.coef_
intercept = linear_svm_3d.intercept_

print(f'coef: {coef}')
print(f'intercept: {intercept}')
print(f'class: {linear_svm_3d.classes_}')

# %%

# data visualization ------------------------------------------

mesh_size = 1.5
#xrange = np.arange(0, mesh2D.index.size, mesh_size)
#yrange = np.arange(0, mesh2D.columns.size, mesh_size)

fig = px.scatter_3d(
    af100_agg,
    x='frag_total',
    y='adapter',
    z='d_Depth',
    color='agg',
    symbol='agg',
    range_x=(0, 250),
    range_y=(0, 250),
    range_z=(0, 250)
)
fig.update_traces(marker_coloraxis=None)
fig.update_layout(
    showlegend=True,
    legend=dict(
        x=-0.1,
        xanchor='left',
        y=1,
        yanchor='auto'
    )
)
fig.update_traces(
    marker=dict(
        size=3.0,
        line=dict(width=1.0, color='DarkSlateGrey')
        #color='white'
    ),
    selector=dict(mode='markers')
)

px0 = np.linspace(0, af100_agg['frag_total'].max(), 100)
px1 = np.linspace(0, af100_agg['adapter'].max(), 100)
px0, px1 = np.meshgrid(px0, px1)

i = 0
for w, b in zip(coef, intercept):
    i += 1
    z = -(px0 * w[0] + px1 * w[1] + b) / w[2]

    if i == 1:
        fig.add_traces(
            go.Surface(
                x=px0,
                y=px1,
                z=z,
                colorscale='Viridis',
                colorbar=dict(title='d_depth')
            )
        )
    else:
        fig.add_traces(
            go.Surface(
                x=px0,
                y=px1,
                z=z,
                colorscale='Viridis',
                showscale=False
            )
        )

fig.update_layout(
    title='display 3D Surface Plots (AgglomerativeClustering)',
    xaxis_nticks=36
)

"""
f = str(round(w[0], 2)) + 'x+' + str(round(w[1], 2)) + 'y' + str(round(w[2], 2))
text = '*f(x,y) =' + f
fig.add_annotation(
    x=-0.1,
    y=0.05,
    text=text,
    font=dict(size=8),
    showarrow=False,
    arrowhead=1,
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
    os.path.abspath('.'),
    '3D_Linear_SVM.html'
)
fig.write_html(htmlfile)

"""
平面モデルによるカーネル法を用いたSVMでは、
3平面間の共通部に対する解釈が難しい
-> クラスター群とは関係ないのない領域の存在が気になる
"""

# %%
"""凝集型クラスタを1-of-K 符号化法に変換
"""
X = af100.loc[:, ['frag_total', 'adapter', 'd_Depth']].to_numpy()
X_new = np.insert(X, 3, y_agg, axis=1)
print(X_new)

# %%
N = len(X)
oneOfk = np.zeros((N, 3))

for i in range(len(y_agg)):
    oneOfk[i, y_agg[i]] = 1

# %%
print(oneOfk)

# %%
# ----------------------------------------
# 関数の定義
# ----------------------------------------
# 尤度関数から交差エントロピー誤差を算出
# 勾配効果法によるパラメータ推定

"""境界線の個数（N(N-1)/2）
・2クラス分類 -> 1個
・3クラス分類 -> 3個
"""

# %%
"""ロジスティック回帰式の定義
"""
def logistic_2class(x0, x1, w):
    """2クラス分類の2次元のロジスティック回帰モデル

    Args:
        x0 ([type]): [description]
        x1 ([type]): [description]
        w ([type]): [description]

    Returns:
        [type]: [description]
    """
    z = w[0] * x0 + w[1] * x1 + w[2]
    return 1 / (1 + np.exp(-z))


def logistic_3class(x0, x1, w):
    """3クラス以上分類の2次元入力のロジスティック回帰モデル

    Args:
        x0 ([type]): [description]
        x1 ([type]): [description]
        w ([type]): [description]
        D ([type]): [description]
        n_class ([type]): [description]

    Returns:
        [type]: [description]
    """

    D = 2
    n_class = 3
    w = w.reshape((int(D+1), n_class))
    N = int(len(x1))
    y = np.zeros((N, n_class))
    for k in range(n_class):
        y[:, k] = np.exp(w[k, 0] * x0 + w[k, 1] * x1 + w[k, 2] * 1)

    ak = y.T / np.sum(y, axis=1)

    return ak.T

# %%
"""交差エントロピー誤差の算出
"""
def cee_logistic_2class(x0, x1, w, t):
    """2クラス分類2次元入力のロジスティック回帰モデルの交差エントロピー誤差

    Args:
        x0 ([type]): [description]
        x1 ([type]): [description]
        w ([type]): [description]
        t ([type]): [description]

    Returns:
        [type]: [description]
    """
    y = logistic_2class(x0, x1, w)
    n_data = len(y)

    cee = 0
    for i in range(n_data):
        cee += t[i, 0] * np.log(y[i]) + (1 - t[i, 0]) * np.log(1 - y[i])

    E = (-1) * (cee / n_data)

    return E


def cee_logistic_3class(w, x, t):
    """3クラス以上分類2次元入力のロジスティック回帰モデルの交差エントロピー誤差

    Args:
        x0 ([type]): [description]
        x1 ([type]): [description]
        w ([type]): [description]
        t ([type]): [description]
        n_class ([type]): [description]

    Returns:
        [type]: [description]
    """
    y = logistic_3class(x[:, 0], x[:, 1], w)
    n_data, n_class = y.shape

    cee = 0
    for i in range(n_data):
        for j in range(n_class):
            cee += t[i, j] * np.log(y[i, j])

    E = (-1) * (cee / n_data)

    return E

# %%
"""交差エントロピー誤差の微分の算出
パラメータ数(={(n_class)x(n_class-1)/2}x(説明変数+1))分の偏微分の値を得る
"""

def dcee_logistic_2class(x, w, t):
    """2クラス分類2次元入力のロジスティック回帰モデルの交差エントロピー誤差の微分

    Args:
        x0 ([type]): [description]
        x1 ([type]): [description]
        w ([type]): [description]
        t ([type]): [description]
        n_class ([type]): [description]

    Returns:
        [type]: [description]
    """
    D = 2
    n_class = 3
    n_parm = (n_class * (n_class - 1) / 2) * (D + 1)
    y = logistic_2class(x[:, 0], x[:, 1], w)
    n_data = len(y)

    dE = np.zeros(n_parm)
    for i in range(n_data):
        item = y[i] - t[i, 0]
        dE[0] += item * x[i, 0]
        dE[1] += item * x[i, 1]
        dE[2] += item * 1
    dE = dE / n_data

    return dE

def dcee_logistic_3class(w, x, t):
    D = 2
    n_class = 3
    y = logistic_3class(x[:, 0], x[:, 1], w)
    n_data, n_class = y.shape

    dE = np.zeros((n_class, int(D+1)))  # ( クラスの数 K) x (x の次元 D+1)
    for i in range(n_data):
        for j in range(n_class):
            dE[j, :] += (y[i, j] - t[i, j]) * np.r_[x[i, :], 1]
    dE = dE / n_data

    return dE.reshape(-1)

# %%
"""共役勾配法によるパラメータ推定
パラメータ引数の指定に注意する
-> cee_logistic_3classの引数の種類と関数で引数を与える時の位置など
"""
from scipy.optimize import minimize

def fit_logistic_2class(w_init, x, t):
    res = minimize(
        cee_logistic_2class,
        w_init,
        args=(x, t),
        jac=dcee_logistic_2class,
        method='CG'
    )

    return res.x

def fit_logistic_3class(w_init, x, t):
    res = minimize(
        cee_logistic_3class,
        w_init,
        args=(x, t),
        jac=dcee_logistic_3class,
        method='CG'
    )

    return res.x


# %%
# testデータ生成 --------------------------------
np.random.seed(seed=1)  # 乱数を固定
N = 100  # データの数
K = 3  # 分布の数
T3 = np.zeros((N, 3), dtype=np.uint8)
T2 = np.zeros((N, 2), dtype=np.uint8)
X = np.zeros((N, 2))
X_range0 = [-3, 3]  # X0 の範囲 , 表示用
X_range1 = [-3, 3]  # X1 の範囲 , 表示用
Mu = np.array([[-.5, -.5], [.5, 1.0], [1, -.5]])  # 分布の中心
Sig = np.array([[.7, .7], [.8, .3], [.3, .8]])  # 分布の分散
Pi = np.array([0.4, 0.8, 1])  # (A) 各分布への割合 0.4 0.8 1
for n in range(N):
    wk = np.random.rand()
    for k in range(K): # (B)
        if wk < Pi[k]:
            T3[n, k] = 1
            break
    for k in range(2):
        X[n, k] = (np.random.randn() * Sig[T3[n, :] == 1, k]
                   + Mu[T3[n, :] == 1, k])
T2[:, 0] = T3[:, 0]
T2[:, 1] = T3[:, 1] | T3[:, 2]

# %%
print(X)
print(T3)

# %%
# test------
W = np.array([1, 2, 3, 4 ,5, 6, 7, 8, 9])
y = logistic_3class(X[:, 0], X[:, 1], W, 3)
print(np.round(y, 3))

# %%
# test------
W = np.array([1, 2, 3, 4 ,5, 6, 7, 8, 9])
cee = cee_logistic_3class(X[:, 0], X[:, 1], W, T3, 3)
print(cee)

# %%
# test------
W = np.array([1, 2, 3, 4 ,5, 6, 7, 8, 9])
dcee = dcee_logistic_3class(X[:, 0], X[:, 1], W, T3, 3)
print(dcee)


# %%
# 実行 --------------------------------
x = np.vstack((X[:, 0], X[:, 1])).T
t = oneOfk
w_init = np.zeros((3, 3))

w = fit_logistic_3class(w_init, x, t)
cee = cee_logistic_3class(w, x, t)

print(np.round(w.reshape((3, 3)), 2))
print('CEE = {0:.2f}'.format(cee))

# %%
# データ処理 ---------------------------
px0 = np.linspace(0, af100_agg['frag_total'].max(), 100)
px1 = np.linspace(0, af100_agg['adapter'].max(), 100)
px0, px1 = np.meshgrid(px0, px1)

# %%
N = int(len(px0))
a = np.zeros((N, N))
y = []
for k in range(3):
    a = np.exp(w[k, 0] * px0 + w[k, 1] * px1 + w[k, 2] * 1)
    y.append(a)

# %%
print(np.array(y).shape)
#-> (3, 100, 100)
# %%
yk = np.array(y)
u = np.sum(yk, axis=0)


# %%
# data visualization ------------------------------------------

mesh_size = 1.5
#xrange = np.arange(0, mesh2D.index.size, mesh_size)
#yrange = np.arange(0, mesh2D.columns.size, mesh_size)

fig = px.scatter_3d(
    af100_agg,
    x='frag_total',
    y='adapter',
    z='d_Depth',
    color='agg',
    symbol='agg',
    range_x=(0, 250),
    range_y=(0, 250),
    range_z=(0, 250)
)
fig.update_traces(marker_coloraxis=None)
fig.update_layout(
    showlegend=True,
    legend=dict(
        x=-0.1,
        xanchor='left',
        y=1,
        yanchor='auto'
    )
)
fig.update_traces(
    marker=dict(
        size=3.0,
        line=dict(width=1.0, color='DarkSlateGrey')
        #color='white'
    ),
    selector=dict(mode='markers')
)

px0 = np.linspace(0, af100_agg['frag_total'].max(), 100)
px1 = np.linspace(0, af100_agg['adapter'].max(), 100)
px0, px1 = np.meshgrid(px0, px1)


for i in range(3):
    z = (yk[i] / u) * af100_agg['d_Depth'].max()

    if i == 0:
        fig.add_traces(
            go.Surface(
                x=px0,
                y=px1,
                z=z,
                colorscale='Viridis',
                colorbar=dict(title='d_depth')
            )
        )
    else:
        fig.add_traces(
            go.Surface(
                x=px0,
                y=px1,
                z=z,
                colorscale='Viridis',
                showscale=False
            )
        )

fig.update_layout(
    title='display 3D Surface Plots (AgglomerativeClustering)',
    xaxis_nticks=36
)

"""
f = str(round(w[0], 2)) + 'x+' + str(round(w[1], 2)) + 'y' + str(round(w[2], 2))
text = '*f(x,y) =' + f
fig.add_annotation(
    x=-0.1,
    y=0.05,
    text=text,
    font=dict(size=8),
    showarrow=False,
    arrowhead=1,
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
    os.path.abspath('.'),
    '3D_logistic_SVM.html'
)
fig.write_html(htmlfile)
# %%
