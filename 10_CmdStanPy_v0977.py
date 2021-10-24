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
