# %%
import pyper
import stan_jupyter as stan
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
