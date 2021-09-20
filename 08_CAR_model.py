# %%

#
# instrinstic CAR modelの実装
#

# 解析データの用意
from numpy.core.fromnumeric import shape
import pyper
import pandas as pd



# %%

# 隣接情報と重み情報を作成
import numpy as np

N = len(df)

adj = []
for i in range(N):
    if i == 0:
        adj.append([i+1])
    elif i == N-1:
        adj.append([i])
    else:
        adj.append([i-1, i+1])

# リスト型からnumpy.array型に変換
adj = np.array(adj)
print(adj)

# %%

weight = []
for i in range(len(adj)):
    w = [1] * len(adj[i])
    weight.append(w)

# リスト型からnumpy.array型に変換
weight = np.array(weight)
print(weight)

# %%

# 隣接情報と重み情報を行列化
amat2 = np.zeros((N, N))
wmat2 = np.zeros((N, N))
for i, a in enumerate(adj):
    amat2[i, a] = 1
    wmat2[i, a] = weight[i]

# %%

print(amat2)
print(wmat2)

# %%

import pymc3 as pm
import theano.tensor as tt
from pymc3.distributions import continuous, distribution

# CAR2関数の読み込み
class CAR2(distribution.Continuous):
    """
    Conditional Autoregressive (CAR) distribution

    Parameters
    ----------
    a : adjacency matrix
    w : weight matrix
    tau : precision at each location
    """

    def __init__(self, w, a, tau, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.a = a = tt.as_tensor_variable(a)
        self.w = w = tt.as_tensor_variable(w)
        self.tau = tau * tt.sum(w, axis=1)
        self.mode = 0.0

    def logp(self, x):
        tau = self.tau
        w = self.w
        a = self.a

        mu_w = tt.sum(x * a, axis=1) / tt.sum(w, axis=1)
        return tt.sum(continuous.Normal.dist(mu=mu_w, tau=tau).logp(x))

# %%

with pm.Model() as CAR_pymc3:
    s = pm.Uniform('s', lower=0, upper=100)
    beta = pm.Normal('beta', mu=0, sd=100)

    tau = 1/(s*s)
    r = CAR2('r', w=wmat2, a=amat2, tau=tau, shape=N)
    link = np.exp(beta + r[np.arange(50)])
    lamda = pm.Deterministic('lamda', np.exp(link))
    y = pm.Poisson(name='y', mu=lamda, observed=df)

# %%

import pymc3 as pm

with pm.Model() as CAR_model:
    rZero = pm.Normal('rZero', mu=0, sd=1)
    sigma = pm.Uniform('s', lower=0, upper=10000)

    r = [0] * N
    r[0] = pm.Normal('r0', mu=rZero, sd=sigma)
    for i in range(1, N):
        r[i] = pm.Normal('r' + str(i), mu=r[i-1], sd=sigma)

    for i in range(0, N):
        #lamda = pm.Deterministic('lamda', np.exp(r[i]))
        y = pm.Poisson(name='y', mu=np.exp(r[i]), observed=df)
    #rr = pm.Normal('rr', mu=r, sd=sigma, observed=df)

# %%

print(r + str(np.arange(50)))

#pm.model_to_graphviz(CAR_model)
# %%

with CAR_model:
    draws=2000
    p_start = pm.find_MAP()
    steps = pm.NUTS()
    #steps=mc.HamiltonianMC()
    tune = 500
    njobs = 3
    random_seed = 10

    trace = pm.sample(draws=draws, start=p_start, steps=steps, tune=tune, cores=1)
# %%
pm.summary(trace)
# %%

np.random.seed(100)
r = np.array([trace['r' + str(i)] for i in range(N)])
lamda = np.exp(r.mean(axis=1))
pred_y = np.random.poisson(lam=lamda)

print(pred_y)

# %%


