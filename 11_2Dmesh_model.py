# %%
import os
import json
import arviz
import numpy as np
import pandas as pd
import cmdstanpy as stan
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

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
    ax.plot(x0, x1, t, 'o',
            color='cornflowerblue', markeredgecolor='black',
            markersize=6)
    ax.view_init(elev=35, azim=-75)

def show_plane(ax, w):
    px0 = np.linspace(0, 16, 5)
    px1 = np.linspace(0, 24, 5)
    px0, px1 = np.meshgrid(px0, px1)
    y = w[0]*px0 + w[1] * px1 + w[2]
    ax.plot_surface(px0, px1, y, rstride=1, cstride=1, alpha=0.3,
                    color='blue', edgecolor='black')

plt.figure(figsize=(6, 5))
ax = plt.subplot(1,1,1,projection='3d')
show_plane(ax, w)
show_data2(ax, x, y, z)
plt.show()

# %%
