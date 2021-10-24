# %%

import stan_jupyter as stan
import pyper
import pandas as pd
import numpy as np
import collections
import matplotlib.pyplot as plt
from matplotlib.figure import figaspect
%matplotlib inline

# %%

#========== read data ==========#
r = pyper.R(use_pandas="True")
r('load(\'/Users/tomoyauchiyama/statisticModel/kubobook_2012/spatial/Y.RData\')')
data = r.get('Y') # class 'numpy.ndarray'

#df.columns = df.columns.str.strip() # DataFrameのカラム名に挿入された空白文字を削除
#print(df.columns.values) #DataFrameのカラム名を列挙

# %%

#========== realize data info ==========#
print(len(data))
print(type(data))

print(pd.Series(data).describe())

print(data)

# %%

#========== visialize obserbed data ==========#
fig=plt.figure()
ax=fig.add_subplot(111)

x = np.arange(0, len(data), 1)
plt.plot(x, data, linestyle='--')
plt.scatter(x, data, label='observed data')

ax.set_xlabel('position[i]')
ax.set_ylabel('# of seeds')
ax.legend(loc="upper left")
ax.set_ylim(min(data)-2, max(data)+3)
ax.legend(loc="upper right")
plt.show()

# %%

# ========== statics modeling ========== #
# 離散値、閾値は0~上限なし
# 空間差、個体差を取り入れた空間状態モデルを考える
# 1階差分のトレンド項と2階差分のトレンド項でモデルを構築してみる
# 平均種子数λ[i] = β + r[i]
# 観測数Y[i] ~ Possion(λ[i])
CARS_model_1 = '''
    data {
        int<lower=0> N;
        int Y[N]; // should not be vector
    }
    parameters {
        vector[N] r;
        real beta;
        real<lower=0> s_r; // define uniform distribution
    }
    transformed parameters {
        real lamda[N];
        for (i in 1:N)
            lamda[i] = exp(beta + r[i]);
    }
    model {
        target += normal_lpdf(r[2:N] | r[1:(N-1)], s_r);
        Y ~ poisson(lamda);
    }
    generated quantities {
        vector[N] Y_mean;
        Y_mean = exp(beta + r);
    }
'''

observed_data = {
    'N': len(data),
    'Y': data
}

# %%
print(observed_data)
# %%
# ========== model compile ========== #
posterior = stan.build(CARS_model_1, data=observed_data)

# %%
# ========== MCMC sampling ========== #
fit = posterior.sample(num_chains=4, num_samples=1000)
# %%
# ========== visialize parameters ========== #
import arviz

arviz.plot_trace(fit)
plt.show()

# %%
# ========== summary ========== #
results = arviz.summary(fit)

print(results)
# %%
CARS_model_2 = '''
    data {
        int<lower=0> N;
        int Y[N]; // should not be vector
    }
    parameters {
        vector[N] r;
        real<lower=0> s_r; // define uniform distribution
    }
    model {
        target += normal_lpdf(r[2:N] | r[1:(N-1)], s_r);
        Y ~ poisson_log(r);
    }
    generated quantities {
        vector[N] Y_mean;
        Y_mean = exp(r);
    }
'''

observed_data = {
    'N': len(data),
    'Y': data
}

# %%
# ========== model compile ========== #
posterior = stan.build(CARS_model_2, data=observed_data)
# %%
# ========== MCMC sampling ========== #
fit = posterior.sample(num_chains=4, num_samples=1000)
# %%
# ========== visialize parameters ========== #
arviz.plot_trace(fit)
plt.show()
# %%
results = arviz.summary(fit)

print(results)
# %%
CARS_model_3 = '''
    data {
        int<lower=0> N;
        int Y[N]; // should not be vector
    }
    parameters {
        vector[N] r;
        real beta;
        real<lower=0> s_r; // define uniform distribution
    }
    model {
        target += normal_lpdf(r[2:N] | r[1:(N-1)], s_r);
        Y ~ poisson_log(beta + r);
    }
    generated quantities {
        vector[N] Y_mean;
        Y_mean = exp(beta + r);
    }
'''

observed_data = {
    'N': len(data),
    'Y': data
}
# %%
# ========== model compile ========== #
posterior = stan.build(CARS_model_3, data=observed_data)

# %%
# ========== MCMC sampling ========== #
fit = posterior.sample(num_chains=4, num_samples=1000)
# %%
# ========== visialize parameters ========== #
arviz.plot_trace(fit)
plt.show()

# %%
results = arviz.summary(fit)

print(results)
# %%
CARS_model_4 = '''
    data {
        int<lower=0> N;
        int Y[N]; // should not be vector
    }
    parameters {
        vector[N] r;
        real beta;
        real<lower=0> s_beta;
        real<lower=0> s_r; // define uniform distribution
    }
    model {
        beta ~ normal(0, s_beta);
        target += normal_lpdf(r[2:N] | r[1:(N-1)], s_r);
        Y ~ poisson_log(beta + r);
    }
    generated quantities {
        vector[N] Y_mean;
        Y_mean = exp(beta + r);
    }
'''

observed_data = {
    'N': len(data),
    'Y': data
}

# %%
# ========== model compile ========== #
posterior = stan.build(CARS_model_4, data=observed_data)

# %%
# ========== MCMC sampling ========== #
fit = posterior.sample(num_chains=4, num_samples=1000)

# %%
# ========== visialize parameters ========== #
arviz.plot_trace(fit)
plt.show()

# %%
results = arviz.summary(fit)

print(results)

# %%
# ========== modeling summary========== #
# model1 (r_hat: 3.06),
# model2 (r_hat: 1.0),
# model3 (r_hat: 2.45),
# model4 (r_hat: 2.35)

# %%
CARS_model_2 = '''
    data {
        int<lower=0> N;
        int Y[N]; // should not be vector
    }
    parameters {
        vector[N] r;
        real<lower=0> s_r; // define uniform distribution
    }
    model {
        target += normal_lpdf(r[2:N] | r[1:(N-1)], s_r);
        Y ~ poisson_log(r);
    }
    generated quantities {
        vector[N] Y_mean;
        Y_mean = exp(r);
    }
'''

observed_data = {
    'N': len(data),
    'Y': data
}

# %%
# ========== model compile ========== #
posterior = stan.build(CARS_model_2, data=observed_data)
# %%
# ========== MCMC sampling ========== #
fit = posterior.sample(num_chains=4, num_samples=1000)
# %%
# ========== visialize parameters ========== #
import arviz

arviz.plot_trace(fit)
plt.show()
# %%
results = arviz.summary(fit)

print(results)

# %%
# ========== preparation for visializing statitics modeling ========== #
probs = (10, 25, 50, 75, 90)
columns = [f'p{p}' for p in probs]

d_est = pd.DataFrame(np.percentile(fit['Y_mean'], probs, axis=1).transpose(), columns=columns)
d_est['x'] = d_est.index + 1

print(d_est)

# %%
# ========== visialize modeling ========== #
fig=plt.figure()
ax=fig.add_subplot(111)

x = np.arange(0, len(data), 1)
plt.plot(x, data, linestyle='--', alpha=0.9)
plt.scatter(x, data, label='observed data', alpha=0.9)

ax.set_xlabel('position[i]')
ax.set_ylabel('# of seeds')
ax.legend(loc="upper left")
ax.set_ylim(min(data)-2, max(data)+3)
ax.legend(loc="upper right")

ax.plot('x', columns[2], data=d_est, color='k')
for j in range(2):
    ax.fill_between('x', columns[j], columns[-j-1],
                    data=d_est, color='k', alpha=0.2*(j+1))

plt.show()

# %%
CARS2_model = '''
    data {
        int<lower=0> N;
        int Y[N]; // should not be vector
    }
    parameters {
        vector[N] r;
        real<lower=0> s_r; // define uniform distribution
    }
    model {
        target += normal_lpdf(r[3:N] | 2*r[2:(N-1)] - r[1:(N-2)], s_r);
        Y ~ poisson_log(r);
    }
    generated quantities {
        vector[N] Y_mean;
        Y_mean = exp(r);
    }
'''

observed_data = {
    'N': len(data),
    'Y': data
}

# %%
print(observed_data)

#%%

# ========== model compile ========== #
posterior = stan.build(CARS2_model, data=observed_data)
# %%
# ========== MCMC sampling ========== #
fit = posterior.sample(num_chains=4, num_samples=1000)
# %%
# ========== visialize parameters ========== #
arviz.plot_trace(fit)
plt.show()
# %%
results = arviz.summary(fit)

print(results)
# %%
# %%
# ========== preparation for visializing statitics modeling ========== #
probs = (10, 25, 50, 75, 90)
columns = [f'p{p}' for p in probs]

d_est = pd.DataFrame(np.percentile(fit['Y_mean'], probs, axis=1).transpose(), columns=columns)
d_est['x'] = d_est.index + 1

print(d_est)

# %%
# ========== visialize modeling ========== #
fig=plt.figure()
ax=fig.add_subplot(111)

x = np.arange(0, len(data), 1)
plt.plot(x, data, linestyle='--', alpha=0.9)
plt.scatter(x, data, label='observed data', alpha=0.9)

ax.set_xlabel('position[i]')
ax.set_ylabel('# of seeds')
ax.legend(loc="upper left")
ax.set_ylim(min(data)-2, max(data)+3)
ax.legend(loc="upper right")

ax.plot('x', columns[2], data=d_est, color='k')
for j in range(2):
    ax.fill_between('x', columns[j], columns[-j-1],
                    data=d_est, color='k', alpha=0.2*(j+1))

plt.show()

# %%
#
# 2次元の空間状態モデルを考える
#

import stan_jupyter as stan
import arviz
import numpy as np
import pandas as pd
from loess.loess_2d import loess_2d
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots
%matplotlib inline


# %%
# ========== read data ========== #
mesh2D = pd.read_csv('/Users/tomoyauchiyama/RStanBook/chap12/input/data-2Dmesh.txt', header=None)

m = np.loadtxt('/Users/tomoyauchiyama/RStanBook/chap12/input/data-2Dmesh.txt', delimiter=',')
print(m)

# %%
# ========== make heatmap ========== #
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
    title='384 plate',
    xaxis_nticks=36)

fig.show()
# %%
# ========== preparation ========== #
d_melt = pd.melt(mesh2D.reset_index(), id_vars='index', ignore_index=True)
d_melt.columns = ('i', 'j', 'Y')

print(d_melt['i'])

# %%
# ========== smoothing ========== #
x = d_melt['i'].to_numpy()
y = d_melt['j'].astype('int64').to_numpy()
z = d_melt['Y'].to_numpy()
#z = x + y

loess_res, _ = loess_2d(x, y, z, frac=0.2)

# %%
print(loess_res)
print(type(loess_res))

# %%

# ========== make smoothed data matrix ========== #
I = mesh2D.index.size
J = mesh2D.columns.size
smoothed = loess_res.reshape(I, J)
print(I)
print(J)
print(smoothed)
print(type(smoothed))

# %%
mesh2D_design = pd.read_csv('/Users/tomoyauchiyama/RStanBook/chap12/input/data-2Dmesh-design.txt', header=None)
dd_melt = pd.melt(mesh2D_design.reset_index(), id_vars='index', ignore_index=True)
dd_melt.columns = ('i', 'j', 'TID')

md = np.loadtxt('/Users/tomoyauchiyama/RStanBook/chap12/input/data-2Dmesh-design.txt', delimiter=',')
print(md)
print(mesh2D_design.to_numpy().reshape(I, J))


# %%
# ========== modeling ========== #
mesh2D_model = '''
    data {
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
'''
# %%

observed_data = {
    'I': int(mesh2D.index.size),
    'J': int(mesh2D.columns.size),
    'Y': mesh2D.to_numpy().reshape(I, J),
    'T': int(mesh2D_design.values.max()),
    'TID': mesh2D_design.to_numpy().reshape(I, J)
}

# %%
print(observed_data)

# %%
observed_data = dict(
    I=mesh2D.index.size,
    J=mesh2D.columns.size,
    Y=m,
    T=int(mesh2D_design.values.max()),
    TID=md
)

# %%
print(observed_data)

# %%
# ========== model compile ========== #
posterior = stan.build(mesh2D_model, data=observed_data)

# %%
# ========== MCMC sampling ========== #
init = {
    'r': smoothed,
    's_r': 1,
    'beta': np.random.normal(0, 0.1, mesh2D_design.values.max()),
    's_beta': 1,
    's_Y': 1
}

fit = posterior.sample(num_chains=4, num_samples=1000, init=init)

# %%
import pyper

r = pyper.R()
r('install.packages(\'reshape2\')')

r('library(\'reshape2\')')
r('d <- as.matrix(read.csv(\'/Users/tomoyauchiyama/RStanBook/chap12/input/data-2Dmesh.txt\', header=F))')
r('I <- nrow(d)')
r('J <- ncol(d)')
r('rownames <- 1:I')
r('colnames <- 1:J')
r('d_melt <- reshape2::melt(d)')
r('colnames(d_melt) <- c(\'i\', \'j\', \'Y\')')

r('d_melt$j = as.numeric(d_melt$j)')

r('loess_res <- loess(Y ~ i + j, data=d_melt, span=0.1)')
r('smoothed <- matrix(loess_res$fitted, nrow=I, ncol=J)')

data = r.get('smoothed') # class 'numpy.ndarray'
# %%
print(data)



# %%
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from sklearn.svm import SVR

df = px.data.iris()

#Xからの写像を作っておく

margin = 0
X = df[['sepal_width', 'sepal_length']]
y = df['petal_width']
model = SVR(C=1.)
model.fit(X, y)

#細かい点(メッシュ,グリッド)を発生

mesh_size = .02
x_min, x_max = X.sepal_width.min() - margin, X.sepal_width.max() + margin
y_min, y_max = X.sepal_length.min() - margin, X.sepal_length.max() + margin
xrange = np.arange(x_min, x_max, mesh_size)
yrange = np.arange(y_min, y_max, mesh_size)
xx, yy = np.meshgrid(xrange, yrange)

#メッシュのすべての点について予測

pred = model.predict(np.c_[xx.ravel(), yy.ravel()])
pred = pred.reshape(xx.shape)

#元の点をplotしてから、x1,x2によるgrid面をzによって押し上げる
#面は全ての点をつなぐsurfaceをつかって描く

fig = px.scatter_3d(df, x='sepal_width', y='sepal_length', z='petal_width')
fig.update_traces(marker=dict(size=1))
fig.add_traces(go.Surface(x=xrange, y=yrange, z=pred, name='pred_surface'))
fig.show()

