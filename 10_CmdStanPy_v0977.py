# %%
import cmdstanpy
cmdstanpy.install_cmdstan()

# %%
import os
import cmdstanpy as stan

#========== TEST1 ==========#
bernoulli_stan = os.path.join(stan.cmdstan_path(), 'examples', 'bernoulli', 'bernoulli.stan')
bernoulli_model = stan.CmdStanModel(stan_file=bernoulli_stan)

# %%
print(bernoulli_model.name)
print(bernoulli_model.stan_file)
print(bernoulli_model.exe_file)
print(bernoulli_model.code())
# %%
bernoulli_data = os.path.join(stan.cmdstan_path(), 'examples', 'bernoulli', 'bernoulli.data.json')
print(bernoulli_data)
# %%
bern_fit = bernoulli_model.sample(data=bernoulli_data)
# %%
print(bern_fit.summary())
# %%
print(bern_fit)
# %%
print(bern_fit.diagnose())

# %%
import arviz
import warnings

warnings.warn('ignore')

arviz.plot_trace(bern_fit)
# %%
results = arviz.summary(bern_fit)
print(results)
print(bern_fit.summary())

# %%
#========== TEST2 ==========#
schools_code = """data {
    int<lower=0> J;
    real y[J];
    real<lower=0> sigma[J];
}
parameters {
    real mu;
    real<lower=0> tau;
    vector[J] eta;
}
transformed parameters {
    vector[J] theta = mu + tau * eta;
}
model {
    target += normal_lpdf(eta | 0, 1);
    target += normal_lpdf(y | theta, sigma);
}
"""
# %%
sc = os.path.join(stan.cmdstan_path(), 'examples', 'schools', 'schools_code.stan')
sm = stan.CmdStanModel(stan_file=sc)

# %%
print(sm.name)
print(sm.stan_file)
print(sm.exe_file)
print(sm.code())


# %%
schools_data = os.path.join(stan.cmdstan_path(), 'examples', 'schools', 'schools_data.json')
fit = sm.sample(data=schools_data)

# %%
print(fit.summary())
# %%
print(fit.diagnose())
# %%
arviz.plot_trace(fit)
# %%
#--------------------------------------
# plot test using arviz
#--------------------------------------

# https://arviz-devs.github.io/arviz/api/plots.html
arviz.plot_energy(fit)

# %%

arviz.plot_mcse(fit)
# %%
arviz.plot_pair(fit)
# %%
arviz.plot_violin(fit)


# %%
#--------------------------------------
# make stan file
#--------------------------------------
code = """data {
    int<lower=0> J;
    real y[J];
    real<lower=0> sigma[J];
}
parameters {
    real mu;
    real<lower=0> tau;
    vector[J] eta;
}
transformed parameters {
    vector[J] theta = mu + tau * eta;
}
model {
    target += normal_lpdf(eta | 0, 1);
    target += normal_lpdf(y | theta, sigma);
}
"""

stan_name = 'schools_code' # have to change

stan_file = os.path.join(os.path.abspath('.'), 'STAN', stan_name + '.stan')

with open(stan_file, 'w') as wf:
    wf.writelines(code)


# %%
import pyper
#========== read data ==========#
r = pyper.R(use_pandas="True")
r('load(\'/Users/tomoyauchiyama/statisticModel/kubobook_2012/spatial/Y.RData\')')
data = r.get('Y') # class 'numpy.ndarray'

#%%
print(data)
# %%
#--------------------------------------
# make json file
#--------------------------------------
import json

observed_data = dict(
    N=len(data),
    Y=data.tolist()
)

json_name = 'schools_data' # have to change

json_file = os.path.join(os.path.abspath('.'), 'STAN', json_name + '.json')

with open(json_file, 'w') as wf:
    json.dump(observed_data, wf, indent=2)

# %%
import pandas as pd
import numpy as np

#--------------------------------------
# read input data
#--------------------------------------

mesh2D = pd.read_csv('/Users/tomoyauchiyama/RStanBook/chap12/input/data-2Dmesh.txt', header=None)
mesh2D_design = pd.read_csv('/Users/tomoyauchiyama/RStanBook/chap12/input/data-2Dmesh-design.txt', header=None)

#m = np.loadtxt('/Users/tomoyauchiyama/RStanBook/chap12/input/data-2Dmesh.txt', delimiter=',')
#md = np.loadtxt('/Users/tomoyauchiyama/RStanBook/chap12/input/data-2Dmesh-design.txt', delimiter=',')

# %%
import plotly.graph_objects as go
from plotly.subplots import make_subplots

#--------------------------------------
# make mesh2D heatmap
#--------------------------------------

y = mesh2D.index.values
x = mesh2D.columns.values

fig = go.Figure(data=go.Heatmap(
        colorbar=dict(
            title='2Dmesh_value'
        ),
        z=mesh2D,
        x=x,
        y=y,
        colorscale='Viridis'))

fig.update_xaxes(title='Plate column')
fig.update_yaxes(title='Plate row')

fig.update_layout(
    title='384 plate (2Dmesh) ',
    xaxis_nticks=36)

fig.show()

# %%
#--------------------------------------
# make mesh2D_design heatmap
#--------------------------------------

y = mesh2D_design.index.values
x = mesh2D_design.columns.values

fig = go.Figure(data=go.Heatmap(
        colorbar=dict(
            title='2Dmesh_d_value'
        ),
        z=mesh2D_design,
        x=x,
        y=y,
        colorscale='Viridis'))

fig.update_xaxes(title='Plate column')
fig.update_yaxes(title='Plate row')

fig.update_layout(
    title='384 plate (2Dmesh_design)',
    xaxis_nticks=36)

fig.show()

# %%
#========== preparation ==========#
d_melt = pd.melt(mesh2D.reset_index(), id_vars='index', ignore_index=True)
d_melt.columns = ('i', 'j', 'Y')

print(d_melt)

# %%
#--------------------------------------
# smoothing for sampling convergence
#--------------------------------------

# usage or example
# https://pypi.org/project/loess/#id11

#==========  smoothing  ==========#
from loess.loess_2d import loess_2d

x = d_melt['i'].to_numpy()
y = d_melt['j'].astype('int64').to_numpy()
z = d_melt['Y'].to_numpy()
#z = x + y

loess_res, _ = loess_2d(x, y, z, frac=0.2)

# %%

#--------------------------------------
# modeling
#--------------------------------------

i = int(mesh2D.index.size)
j = int(mesh2D.columns.size)

print(loess_res.reshape(i, j))

#--------------------------------------
# make stan file
#--------------------------------------
code = """data {
    int I;
    int J;
    real Y[I, J];
    int T;
    int<lower=0, upper=T> TID[I, J];
}
parameters {
    real r[I, J];
    real<lower=0> s_r;
    vector[T] beta;
    real<lower=0> s_beta;
    real<lower=0> s_Y;
}
model {
    for (i in 1:I)
        for (j in 3:J)
            target += normal_lpdf(r[i, j] | 2*r[i, j-1] - r[i, j-2], s_r);

    for (i in 3:I)
        for (j in 1:J)
            target += normal_lpdf(r[i, j] | 2*r[i-1, j] - r[i-2, j], s_r);

    beta ~ student_t(6, 0, s_beta);
    for (i in 1:I)
        for (j in 1:J)
            Y[i, j] ~ normal(r[i, j] + beta[TID[i, j]], s_Y);
}
"""

stan_name = 'mesh2D_model' # have to change

stan_file = os.path.join(os.path.abspath('.'), 'STAN', stan_name + '.stan')

with open(stan_file, 'w') as wf:
    wf.writelines(code)

print(stan_file)

# %%
#--------------------------------------
# make json file
#--------------------------------------
import json

i = int(mesh2D.index.size)
j = int(mesh2D.columns.size)

observed_data = dict(
    I=i,
    J=j,
    Y=mesh2D.to_numpy().reshape(i, j).tolist(),
    T=int(mesh2D_design.values.max()),
    TID=mesh2D_design.to_numpy().reshape(i, j).tolist()
)

json_name = 'mesh2D_data' # have to change

json_file = os.path.join(os.path.abspath('.'), 'STAN', json_name + '.json')

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

init = dict(
    r=loess_res.reshape(i, j),
    s_r=1,
    beta=np.random.normal(0, 0.1, mesh2D_design.values.max()),
    s_beta=1,
    s_Y=1
)

print(init)

# %%
#--------------------------------------
# MCMC sampling (NUTS method)
#--------------------------------------

fit = sm.sample(data=json_file, iter_warmup=500, thin=5, seed=1234, inits=init)

# %%
#--------------------------------------
# MCMC sampling summary
#--------------------------------------

# summaries of
# the total joint log-probability density lp__
# plus all model parameters
# and quantities of interest in a pandas.DataFrame

# type -> pandas (basic statistic)
print(fit.summary())
print(fit.summary().loc['s_Y', 'Mean'])

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
#--------------------------------------
# statistic modeling visualization
#--------------------------------------

#======  matplotlib method  ======#
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

%matplotlib inline

r_mean = np.median(parm_r, axis=0)
I, J = r_mean.shape
ii, jj = np.mgrid[:I, :J]

ax = Axes3D(plt.gcf())
ax.plot_wireframe(ii, jj, r_mean, color='k')
ax.plot_surface(ii, jj, r_mean, color='k', alpha=0.2)
plt.setp(ax, xlabel='Plate Row', ylabel='Plate Column', zlabel='r')

plt.show()

# %%
d_melt = pd.melt(mesh2D.reset_index(), id_vars='index', ignore_index=True)
d_melt.columns = ('Plate Row', 'Plate Column', 'Y')

# %%
#======  plot 3D graph with plotly  ======#
import plotly.express as px
import plotly.graph_objects as go

I, J = r_mean.shape
mesh_size = 1.5
xrange = np.arange(0, I, mesh_size)
yrange = np.arange(0, J, mesh_size)

fig = px.scatter_3d(d_melt, x='Plate Column', y='Plate Row', z='Y')
fig.update_traces(marker=dict(size=1.0,
                              line=dict(width=0.5,color='DarkSlateGrey'),
                              color='white'),
                  selector=dict(mode='markers'))
fig.add_traces(go.Surface(z=r_mean, colorscale='Viridis',
                          colorbar=dict(title='pred_Y_value')))

fig.update_layout(
    title='statistic modeling (2D_mesh)',
    xaxis_nticks=36)

fig.show()

# %%
#
# 患者別の薬効濃度を予測する非線形型の階層ベイズモデルを構築する
#
"""構築後、線形ガウス基底回帰モデルと比較する
"""

import pandas as pd

cons_data = pd.read_csv('/Users/tomoyauchiyama/RStanBook/chap08/input/data-conc-2.txt')

print(cons_data)

# %%
print(cons_data.describe())


# %%

"""
xx= cons_data.columns.str.extract(r'Time(\d+)')[0].to_list()[1:]
x = [int(i) for i in xx]
max_x = max(x)
print(max_x)
min_y = cons_data.iloc[1, 1:].max()
print(min_y)

n = len(cons_data.index.values)
for ii in range(n):
    if (ii+1) % n_row == 0:
        i = n_row
        j = n_row
    else:
        i = (ii+1) % n_row
        j = (ii+1) % n_row
    print(i, j)
"""

# %%
import seaborn as sns
import plotly.graph_objects as go
from plotly.subplots import make_subplots

#======  data visualization  ======#
def get_colorpalette(n_colors):
    """legendとRBG値を対応させる関数

    Returns:
        rgb (list): RBG値が格納された配列
    """
    palette = sns.color_palette('hls', n_colors) #colorpalette = 'hls'
    rgb = ['rgb({},{},{})'.format(*[x*256 for x in rgb]) for rgb in palette]
    return rgb


n = len(cons_data.index.values)
colors = get_colorpalette(n)

n_row = 4
n_col = 4

y_title = 'Conc(mg/ml)'
x_title = 'Time(hour)'
main_title = 'effect time each person'

titles = []
for i in range(n):
    ID = int(i + 1)
    name = 'PersonID_' + str(ID)
    titles.append(name)

xx= cons_data.columns.str.extract(r'Time(\d+)')[0].to_list()[1:]
x = [int(i) for i in xx]

max_x = max(x) + 3
min_x = -0.2
max_y = cons_data.describe().loc['max', :].max() + 3
min_y = 0

fig = make_subplots(rows=n_row,
                    cols=n_col,
                    subplot_titles=titles,
                    horizontal_spacing=0.1  # will change
                    #vertical_spacing=0.15     # will change
                    )

ii = -1
for i in range(1, n_row+1):
    for j in range(1, n_col+1):
        ii += 1

        y = cons_data.iloc[ii, 1:]

        fig.add_trace(go.Scatter(x=x, y=y, mode='lines+markers', line=dict(color=colors[ii])),
                    row=i,
                    col=j)
        fig.update_layout(plot_bgcolor='white',
                        height=800,
                        width=900)
        fig.update_xaxes(title=x_title,
                        showline=True,
                        linewidth=1,
                        linecolor='black',
                        mirror=True,
                        ticks='inside',
                        range=(min_x, max_x),
                        dtick=5,
                        row=i,
                        col=j)
        fig.update_yaxes(title=y_title,
                        showline=True,
                        linewidth=1,
                        linecolor='black',
                        mirror=True,
                        ticks='inside',
                        range=(min_y, max_y),
                        dtick=5,
                        row=i,
                        col=j)

fig.for_each_xaxis(lambda axis: axis.title.update(font=dict(color='black', size=10)))
fig.for_each_yaxis(lambda axis: axis.title.update(font=dict(color='black', size=10)))

fig.update_layout(title=dict(text=main_title,
                             x=0.5,
                             xanchor='center'),
                  showlegend=False)

fig.update_annotations(font=dict(size=10))

# %%
#--------------------------------------
# make stan file
#--------------------------------------
code = """data {
    int Np;
    int Nt;
    vector[Nt] T;
    matrix[Np, Nt] Y;
    int new_T;
    real new_Time[new_T];
}
parameters {
    real m_a;
    real m_b;
    vector[Np] log_a;
    vector[Np] log_b;
    real<lower=0> s_a;
    real<lower=0> s_b;
    real <lower=0> s_Y;
}
transformed parameters {
    vector[Np] a;
    vector[Np] b;
    matrix[Np, Nt] mu;
    a = exp(log_a);
    b = exp(log_b);
    for (j in 1:Nt)
        for (i in 1:Np)
            mu[i, j] = a[i]*(1-exp(-b[i]*T[j]));
}
model {
    log_a ~ normal(m_a, s_a);
    log_b ~ normal(m_b, s_b);
    to_vector(Y) ~ normal(to_vector(mu), s_Y);
}
generated quantities {
    real new_Y[Np, new_T];
    for (i in 1:Np)
        for (j in 1:new_T)
            new_Y[i, j] = normal_rng(a[i]*(1-exp(-b[i]*new_Time[j])), s_Y);
}
"""

stan_name = 'nonlinear_model_4' # have to change

stan_file = os.path.join(os.path.abspath('.'), 'STAN', stan_name + '.stan')

with open(stan_file, 'w') as wf:
    wf.writelines(code)

print(stan_file)

# %%
#--------------------------------------
# make json file
#--------------------------------------
import numpy as np
import json

col = cons_data.columns.str.match(r'Time(\d+)')
Time = cons_data.columns[col].str.extract(r'Time(\d+)').astype(int)[0].to_list()
Np = cons_data.index.size
Nt = len(Time)
Y = cons_data.loc[:, ['Time1', 'Time2', 'Time4', 'Time8', 'Time12', 'Time24']].to_numpy().reshape(Np, Nt).tolist()
new_T = 60
new_Time = np.linspace(0, max(Time), new_T).tolist()


observed_data = dict(
    Np=cons_data.index.size,
    Nt=len(Time),
    T=Time,
    Y=Y,
    new_T=new_T,
    new_Time=new_Time
)

json_name = 'nonlinear_data' # have to change

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

fit = sm.sample(data=json_file, chains=3, iter_sampling=2500, iter_warmup=500, thin=5, seed=1234)

# %%
#--------------------------------------
# MCMC sampling summary
#--------------------------------------

# summaries of
# the total joint log-probability density lp__
# plus all model parameters
# and quantities of interest in a pandas.DataFrame

# type -> pandas (basic statistic)
pd.set_option('display.max_columns', None)
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
#
# 欠損値を含む患者別の薬効濃度を予測する非線形型の階層ベイズモデルを構築する
#

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
#--------------------------------------------
# 欠損値を含むデータを縦長のデータフレームに加工
#--------------------------------------------

"""
#======  read row data including N.A  ======#
long_df = pd.read_csv('/Users/tomoyauchiyama/RStanBook/chap09/input/data-conc-2-NA-long.txt')

print(len(long_df))
print(long_df)
"""

#======  read row data including N.A  ======#
wide_df = pd.read_csv('/Users/tomoyauchiyama/RStanBook/chap09/input/data-conc-2-NA-wide.txt')

print(wide_df)

# %%
#======  convert wide DataFreme to long DataFrame  ======#
"""remove row including N.A
   substitute int type of number to add Time value
   add TimeID column linked to Time value after sorting DataFrame by 'Time'
   sort DataFrame by 'PersonID' and 'TimeID'
   reset index
"""

long_df = pd.melt(wide_df, id_vars='PersonID', ignore_index=True).dropna()
long_df.columns = ('PersonID', 'Time', 'Y')
long_df['Time'] = long_df.loc[:, 'Time'].str.extract(r'Time(\d+)').astype(int)
long_df_sort = long_df.sort_values('Time')
long_df_sort['TimeID'] = long_df_sort['Time'].factorize()[0] + 1
long_df_sort = long_df_sort.sort_values(['PersonID', 'Time']).reset_index()

print(len(long_df_sort))
print(long_df_sort)

# %%
#======  prepare for data visualization ======#
X = {}
Y = {}
for i in range(len(long_df_sort)):
    ID = int(long_df_sort.loc[i, 'PersonID'])
    x_val = str(long_df_sort.loc[i, 'Time'])
    y_val = str(long_df_sort.loc[i, 'Y'])
    if X.get(ID):
        X[ID] += ',' + x_val
        Y[ID] += ',' + y_val
    else:
        X[ID] = x_val
        Y[ID] = y_val

proc_df = pd.DataFrame()
for ID in sorted(X.keys()):
    xx = [int(i) for i in X[ID].split(',')]
    yy = [float(i) for i in Y[ID].split(',')]
    col = pd.Series([ID, xx, yy])
    proc_df = proc_df.append(col, ignore_index=True)

proc_df.columns = ['PersonID', 'Time', 'Y']

print(proc_df)

# %%
#--------------------------------------------
# 欠損値を含むデータの可視化とモデリング
#--------------------------------------------

#======  data visualization  ======#
def get_colorpalette(n_colors):
    """legendとRBG値を対応させる関数

    Returns:
        rgb (list): RBG値が格納された配列
    """
    palette = sns.color_palette('hls', n_colors) #colorpalette = 'hls'
    rgb = ['rgb({},{},{})'.format(*[x*256 for x in rgb]) for rgb in palette]
    return rgb


n = len(sorted(long_df.loc[:, 'PersonID'].unique().tolist()))
colors = get_colorpalette(n)

n_row = 4
n_col = 4

y_title = 'Conc(mg/ml)'
x_title = 'Time(hour)'
main_title = 'effect time each person(including N.A)'

titles = []
for i in range(n):
    ID = int(i + 1)
    name = 'PersonID_' + str(ID)
    titles.append(name)

max_x = long_df_sort.loc[:, 'Time'].unique().max() + 3
min_x = -0.2
max_y = wide_df.describe().loc['max', :].max() + 3
min_y = 0

fig = make_subplots(rows=n_row,
                    cols=n_col,
                    subplot_titles=titles,
                    horizontal_spacing=0.1  # will change
                    #vertical_spacing=0.15     # will change
                    )

ii = -1
for i in range(1, n_row+1):
    for j in range(1, n_col+1):
        ii += 1

        y = proc_df.loc[ii, 'Y']
        x = proc_df.loc[ii, 'Time']

        fig.add_trace(go.Scatter(x=x, y=y, mode='lines+markers', line=dict(color=colors[ii])),
                    row=i,
                    col=j)
        fig.update_layout(plot_bgcolor='white',
                        height=800,
                        width=900)
        fig.update_xaxes(title=x_title,
                        showline=True,
                        linewidth=1,
                        linecolor='black',
                        mirror=True,
                        ticks='inside',
                        range=(min_x, max_x),
                        dtick=5,
                        row=i,
                        col=j)
        fig.update_yaxes(title=y_title,
                        showline=True,
                        linewidth=1,
                        linecolor='black',
                        mirror=True,
                        ticks='inside',
                        range=(min_y, max_y),
                        dtick=5,
                        row=i,
                        col=j)

fig.for_each_xaxis(lambda axis: axis.title.update(font=dict(color='black', size=10)))
fig.for_each_yaxis(lambda axis: axis.title.update(font=dict(color='black', size=10)))

fig.update_layout(title=dict(text=main_title,
                             x=0.5,
                             xanchor='center'),
                  showlegend=False)

fig.update_annotations(font=dict(size=10))

# %%
#--------------------------------------
# make stan file
#--------------------------------------
code = """data {
    int Na;
    int Np;
    int Nt;
    vector[Nt] T;
    int<lower=1, upper=Np> PersonID[Na];
    int<lower=1, upper=Nt> TimeID[Na];
    vector[Na] Y;
    int new_T;
    real new_Time[new_T];
}
parameters {
    real m_a;
    real m_b;
    vector[Np] log_a;
    vector[Np] log_b;
    real<lower=0> s_a;
    real<lower=0> s_b;
    real <lower=0> s_Y;
}
transformed parameters {
    vector[Np] a;
    vector[Np] b;
    matrix[Np, Nt] mu;
    a = exp(log_a);
    b = exp(log_b);
    for (j in 1:Nt)
        for (i in 1:Np)
            mu[i, j] = a[i]*(1-exp(-b[i]*T[j]));
}
model {
    log_a ~ normal(m_a, s_a);
    log_b ~ normal(m_b, s_b);
    for (i in 1:Na){
        Y[i] ~ normal(mu[PersonID[i], TimeID[i]], s_Y);
    }
}
generated quantities {
    real new_Y[Np, new_T];
    for (i in 1:Np)
        for (j in 1:new_T)
            new_Y[i, j] = normal_rng(a[i]*(1-exp(-b[i]*new_Time[j])), s_Y);
}
"""

stan_name = 'nonlinear_model_5' # have to change

stan_file = os.path.join(os.path.abspath('.'), 'STAN', stan_name + '.stan')

with open(stan_file, 'w') as wf:
    wf.writelines(code)

print(stan_file)


# %%
#--------------------------------------
# make json file
#--------------------------------------
import json

Na = len(long_df_sort)
Time = sorted(long_df_sort.loc[:, 'Time'].unique().tolist())
Np = len(long_df_sort.loc[:, 'PersonID'].unique())
Nt = len(Time)
PersonID = long_df_sort.loc[:, 'PersonID'].tolist()
TimeID = long_df_sort.loc[:, 'TimeID'].tolist()
Y = long_df_sort.loc[:, 'Y'].tolist()
new_T = 60
new_Time = np.linspace(0, max(Time), new_T).tolist()


observed_data = dict(
    Na=Na,
    Np=Np,
    Nt=Nt,
    T=Time,
    PersonID=PersonID,
    TimeID=TimeID,
    Y=Y,
    new_T=new_T,
    new_Time=new_Time
)

json_name = 'nonlinear_NA_data' # have to change

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

fit_na = sm.sample(data=json_file, chains=3, iter_sampling=2500, iter_warmup=500, thin=5, seed=1234)

# %%
#--------------------------------------
# MCMC sampling summary
#--------------------------------------

# summaries of
# the total joint log-probability density lp__
# plus all model parameters
# and quantities of interest in a pandas.DataFrame

# type -> pandas (basic statistic)
pd.set_option('display.max_columns', None)
print(fit_na.summary().loc[:, ['Mean', 'StdDev', '5%', '50%', '95%', 'N_Eff', 'R_hat']])

# %%
#--------------------------------------
# MCMC sampling analyzation
#--------------------------------------

# analyze the per-draw sampler parameters across all chains
# looking for potential problems which indicate that
# the sample isn’t a representative sample from the posterior

print(fit_na.diagnose())
# %%
#--------------------------------------
# parameter visualization
#--------------------------------------
import arviz

arviz.plot_trace(fit_na)

# %%
#--------------------------------------------
# 欠損値を含まない場合のモデリング結果の可視化
#--------------------------------------------
#--------------------------------------
# parameter data extraction
#--------------------------------------

#====== extract parmetar 'new_Y' ======#
parm_Y = fit.stan_variable('new_Y')

#====== realize array shape  ======#
print(parm_Y.shape)
# output: (1500, 16, 60)
# -> (16 people x 60 bin of time) x 1500 (= 500 samplings x 3 chains)

# %%
#====== realize variant type ======#
print(type(parm_Y))
# output: <class 'numpy.ndarray'>

#====== calc ave of parm_r   ======#
Y_mean = np.median(parm_Y, axis=0)

#====== realize array shape  ======#
print(Y_mean.shape)
# output: (16, 60)
# -> (16x60 array) x 1

# %%
probs = (10, 25, 50, 75, 90)
column = pd.Series([f'p{p}' for p in probs])
Time = sorted(long_df_sort.loc[:, 'Time'].unique().tolist())
new_T = 60
new_Time = np.linspace(0, max(Time), new_T)

#======  convert long DataFrame  ======#
# -> 60 bins of time x (16 persons x [p10, p25, p50, p75, p90])
prob_Y = np.percentile(parm_Y, probs, axis=0).transpose()

est_d = pd.DataFrame()
for i, time_id in enumerate(prob_Y):
    time_val = new_Time[i]
    for person_id, prob_y in enumerate(time_id, 1):
        df1 = pd.Series([person_id, time_val])
        df2 = pd.Series(prob_y)
        col = pd.concat([df1, df2])
        est_d = est_d.append(col, ignore_index=True)

est_d.columns = pd.concat([pd.Series(['PersonID', 'Time']), column])
est_d = est_d.sort_values(['PersonID', 'Time']).reset_index()

pd.set_option('display.max_rows', None)
print(est_d)

# %%
#--------------------------------------
# statistic modeling visualization
#--------------------------------------

#=== prepare for data visualization ===#
p10 = {}
p25 = {}
p50 = {}
p75 = {}
p90 = {}

for i in range(len(est_d)):
    ID = int(est_d.loc[i, 'PersonID'])
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
for ID in sorted(X.keys()):
    xx = np.linspace(0, max(Time), new_T).tolist()
    p10y = [float(i) for i in p10[ID].split(',')]
    p25y = [float(i) for i in p25[ID].split(',')]
    p50y = [float(i) for i in p50[ID].split(',')]
    p75y = [float(i) for i in p75[ID].split(',')]
    p90y = [float(i) for i in p90[ID].split(',')]
    col = pd.Series([ID, xx, p10y, p25y, p50y, p75y, p90y])
    pred_df = pred_df.append(col, ignore_index=True)

pred_df.columns = ['PersonID', 'Time',
                   'p10', 'p25', 'p50', 'p75', 'p90']

print(pred_df)

# %%
#======  data visualization with generated quantities  ======#

"""欠損値を含まない場合の可視化
"""

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


n = len(sorted(long_df.loc[:, 'PersonID'].unique().tolist()))
colors = get_colorpalette(n)

n_row = 4
n_col = 4

y_title = 'Conc(mg/ml)'
x_title = 'Time(hour)'
main_title = 'effect time each person'

titles = []
for i in range(n):
    ID = int(i + 1)
    name = 'PersonID_' + str(ID)
    titles.append(name)

max_x = normal_l_df_sort.loc[:, 'Time'].unique().max() + 3
min_x = -0.2
max_y = normal_df.describe().loc['max', :].max() + 3
min_y = 0

probs = (10, 25, 50, 75, 90)
pxx = [f'p{p}' for p in probs]

fig = make_subplots(rows=n_row,
                    cols=n_col,
                    subplot_titles=titles,
                    horizontal_spacing=0.1  # will change
                    #vertical_spacing=0.15     # will change
                    )

ii = -1
for i in range(1, n_row+1):
    for j in range(1, n_col+1):
        ii += 1
        y = proc_n_df.loc[ii, 'Y']
        x = proc_n_df.loc[ii, 'Time']
        pre_y = pred_df.loc[ii, 'p50']
        pre_x = pred_df.loc[ii, 'Time']

        # fill area between trace0 and trace1
        for k in range(2):
            fig.add_trace(go.Scatter(
                x=pre_x,
                y=pred_df.loc[ii, pxx[k]],
                fill=None,
                mode='lines',
                line=dict(color='black', width=0.4),
                opacity=0.2*(k+1)
                ),
                row=i,
                col=j)
            fig.add_trace(go.Scatter(
                x=pre_x,
                y=pred_df.loc[ii, pxx[-k-1]],
                fill='tonexty',
                mode='lines',
                line=dict(color='black', width=0.4),
                opacity=0.2*(k+1)
                ),
                row=i,
                col=j)
        fig.add_trace(go.Scatter(
            x=pre_x,
            y=pre_y,
            mode='lines',
            line=dict(color=colors[ii], width=2)
            ),
            row=i,
            col=j)
        fig.add_trace(go.Scatter(
            x=x,
            y=y,
            mode='markers',
            marker=dict(
                color=colors[ii],
                size=5,
                line=dict(color='black', width=1)
                )
            ),
            row=i,
            col=j)
        fig.update_layout(plot_bgcolor='white',
                        height=800,
                        width=900)
        fig.update_xaxes(title=x_title,
                        showline=True,
                        linewidth=1,
                        linecolor='black',
                        mirror=True,
                        ticks='inside',
                        range=(min_x, max_x),
                        dtick=5,
                        row=i,
                        col=j)
        fig.update_yaxes(title=y_title,
                        showline=True,
                        linewidth=1,
                        linecolor='black',
                        mirror=True,
                        ticks='inside',
                        range=(min_y, max_y),
                        dtick=5,
                        row=i,
                        col=j)

fig.for_each_xaxis(lambda axis: axis.title.update(
    font=dict(color='black', size=10)))
fig.for_each_yaxis(lambda axis: axis.title.update(
    font=dict(color='black', size=10)))

fig.update_layout(
    title=dict(
        text=main_title,
        x=0.5,
        xanchor='center'
    ),
    showlegend=False)

fig.update_annotations(font=dict(size=10))

# %%
#--------------------------------------------
# 欠損値を含む場合のモデリング結果の可視化
#--------------------------------------------
#--------------------------------------
# parameter data extraction
#--------------------------------------

#====== extract parmetar 'new_Y' ======#
parm_Y = fit_na.stan_variable('new_Y')

#====== realize array shape  ======#
print(parm_Y.shape)
# output: (1500, 16, 60)
# -> (16 people x 60 bin of time) x 1500 (= 500 samplings x 3 chains)

# %%
#====== realize variant type ======#
print(type(parm_Y))
# output: <class 'numpy.ndarray'>

#====== calc ave of parm_r   ======#
Y_mean = np.median(parm_Y, axis=0)

#====== realize array shape  ======#
print(Y_mean.shape)
# output: (16, 60)
# -> (16x60 array) x 1

# %%
probs = (10, 25, 50, 75, 90)
column = pd.Series([f'p{p}' for p in probs])
Time = sorted(long_df_sort.loc[:, 'Time'].unique().tolist())
new_T = 60
new_Time = np.linspace(0, max(Time), new_T)

#======  convert long DataFrame  ======#
# -> 60 bins of time x (16 persons x [p10, p25, p50, p75, p90])
prob_Y = np.percentile(parm_Y, probs, axis=0).transpose()

est_d = pd.DataFrame()
for i, time_id in enumerate(prob_Y):
    time_val = new_Time[i]
    for person_id, prob_y in enumerate(time_id, 1):
        df1 = pd.Series([person_id, time_val])
        df2 = pd.Series(prob_y)
        col = pd.concat([df1, df2])
        est_d = est_d.append(col, ignore_index=True)

est_d.columns = pd.concat([pd.Series(['PersonID', 'Time']), column])
est_d = est_d.sort_values(['PersonID', 'Time']).reset_index()

pd.set_option('display.max_rows', None)
print(est_d)

# %%
#--------------------------------------
# statistic modeling visualization
#--------------------------------------

#=== prepare for data visualization ===#
p10 = {}
p25 = {}
p50 = {}
p75 = {}
p90 = {}

for i in range(len(est_d)):
    ID = int(est_d.loc[i, 'PersonID'])
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
for ID in sorted(X.keys()):
    xx = np.linspace(0, max(Time), new_T).tolist()
    p10y = [float(i) for i in p10[ID].split(',')]
    p25y = [float(i) for i in p25[ID].split(',')]
    p50y = [float(i) for i in p50[ID].split(',')]
    p75y = [float(i) for i in p75[ID].split(',')]
    p90y = [float(i) for i in p90[ID].split(',')]
    col = pd.Series([ID, xx, p10y, p25y, p50y, p75y, p90y])
    pred_df = pred_df.append(col, ignore_index=True)

pred_df.columns = ['PersonID', 'Time',
                   'p10', 'p25', 'p50', 'p75', 'p90']

print(pred_df)

# %%
#======  data visualization with generated quantities  ======#

"""欠損値を含む場合の可視化
"""

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


n = len(sorted(long_df.loc[:, 'PersonID'].unique().tolist()))
colors = get_colorpalette(n)

n_row = 4
n_col = 4

y_title = 'Conc(mg/ml)'
x_title = 'Time(hour)'
main_title = 'effect time each person without N.A'

titles = []
for i in range(n):
    ID = int(i + 1)
    name = 'PersonID_' + str(ID)
    titles.append(name)

max_x = long_df_sort.loc[:, 'Time'].unique().max() + 3
min_x = -0.2
max_y = wide_df.describe().loc['max', :].max() + 3
min_y = 0

probs = (10, 25, 50, 75, 90)
pxx = [f'p{p}' for p in probs]

fig = make_subplots(rows=n_row,
                    cols=n_col,
                    subplot_titles=titles,
                    horizontal_spacing=0.1  # will change
                    #vertical_spacing=0.15     # will change
                    )

ii = -1
for i in range(1, n_row+1):
    for j in range(1, n_col+1):
        ii += 1
        y = proc_df.loc[ii, 'Y']
        x = proc_df.loc[ii, 'Time']
        pre_y = pred_df.loc[ii, 'p50']
        pre_x = pred_df.loc[ii, 'Time']

        # fill area between trace0 and trace1
        for k in range(2):
            fig.add_trace(go.Scatter(
                x=pre_x,
                y=pred_df.loc[ii, pxx[k]],
                fill=None,
                mode='lines',
                line=dict(color='black', width=0.4),
                opacity=0.2*(k+1)
                ),
                row=i,
                col=j)
            fig.add_trace(go.Scatter(
                x=pre_x,
                y=pred_df.loc[ii, pxx[-k-1]],
                fill='tonexty',
                mode='lines',
                line=dict(color='black', width=0.4),
                opacity=0.2*(k+1)
                ),
                row=i,
                col=j)
        fig.add_trace(go.Scatter(
            x=pre_x,
            y=pre_y,
            mode='lines',
            line=dict(color=colors[ii], width=2)
            ),
            row=i,
            col=j)
        fig.add_trace(go.Scatter(
            x=x,
            y=y,
            mode='markers',
            marker=dict(
                color=colors[ii],
                size=5,
                line=dict(color='black', width=1)
                )
            ),
            row=i,
            col=j)
        fig.update_layout(plot_bgcolor='white',
                        height=800,
                        width=900)
        fig.update_xaxes(title=x_title,
                        showline=True,
                        linewidth=1,
                        linecolor='black',
                        mirror=True,
                        ticks='inside',
                        range=(min_x, max_x),
                        dtick=5,
                        row=i,
                        col=j)
        fig.update_yaxes(title=y_title,
                        showline=True,
                        linewidth=1,
                        linecolor='black',
                        mirror=True,
                        ticks='inside',
                        range=(min_y, max_y),
                        dtick=5,
                        row=i,
                        col=j)

fig.for_each_xaxis(lambda axis: axis.title.update(
    font=dict(color='black', size=10)))
fig.for_each_yaxis(lambda axis: axis.title.update(
    font=dict(color='black', size=10)))

fig.update_layout(
    title=dict(
        text=main_title,
        x=0.5,
        xanchor='center'
    ),
    showlegend=False)

fig.update_annotations(font=dict(size=10))
