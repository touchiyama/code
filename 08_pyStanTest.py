# %%

#
# pystanを用いて、統計モデリングを行ってみる
#

import stan
import stan_jupyter as stan
import pandas as pd
import numpy as np
import collections
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score

# %%

#
# scikit-learnを用いて、サンプルデータの生成を行う
#

def generate_sampleData(num,seed):
    YY = [] # 応答変数
    XX = [] # 説明変数

    XX_num = 8
    intercept = 0.2
    weight = [0.2, 0.3, 0.5, -0.4, 0.1, 0.2, 0.5, -0.3] # 各特徴量の重み

    np.random.seed(seed=seed)
    for i in range(num):
        X = [np.random.rand() for j in range(XX_num)]
        noise = [np.random.normal(0, 0.1) for j in range(XX_num)]
        Y = sum([intercept + X[j]*weight[j] + noise[j] for j in range(XX_num)])

        YY.append(Y)
        XX.append(X)

    df = pd.DataFrame(np.c_[YY, XX],
                      columns = ['Y', 'X0', 'X1', 'X2', 'X3',
                                 'X4', 'X5', 'X6', 'X7']
                     )

    return df

data = generate_sampleData(1000, 0)

X = data.drop('Y', axis=1)
Y = data['Y']

X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=0)

# %%
print(X_train)

# %%
print(X_test)

# %%

#
# ベイズ線形回帰をモデリング　
#

sample_code = """
    data {
        int<lower=0> N;
        int<lower=0> D;
        matrix[N, D] X;
        vector[N] y;
        int<lower=0> N_new;
        matrix[N_new, D] X_new;
    }
    parameters {
        real w0;
        vector[D] w;
        real<lower=0> sigma;
    }
    model {
        for (i in 1:N)
            y[N] ~ normal(w0 + dot_product(X[i], w), sigma);
    }
    generated quantities {
        vector[N_new] y_new;
        for (i in 1:N)
            y_new[i] = normal_rng(w0 + dot_product(X[i], w), sigma);
    }
"""

sample_data = {
    'N': X_train.shape[0],
    'D': X_train.shape[1],
    'X': X_train,
    'y': Y_train,
    'N_new': X_test.shape[0],
    'X_new': X_test
}

#%%

posterior = stan.build(sample_code, data=sample_data)

# %%

schools_code = """
    data {
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

schools_data = {'J': 8,
                'y': [28, 8, -3, 7, -1, 1, 18, 12],
                'sigma': [15, 10, 16, 11, 9, 11, 10, 18]}

# %%
print(schools_data)

# %%

sm = pystan.StanModel(model_code=schools_code)
#posterior = stan.build(schools_code, data=schools_data)

#%%

fit = sm.sampling(data=schools_data, iter=1000, chains=4)
#fit = posterior.sample(num_chains=4, num_samples=1000)

# %%
eta = fit['eta']
print(eta)
df = fit.to_frame()
print(df)



# %%
import json

#schools_data =  pd.DataFrame({ "J": 8,
#                  "y": [28, 8, -3, 7, -1, 1, 18, 12],
#                  "sigma": [15, 10, 16, 11, 9, 11, 10, 18]
#                }).to_json()

#json_open = open('qitta_json.json', 'r')
json_str = '''
{
    "test":"json",
    "test2":"json2"
}
'''


json_load = json.load(json_str)

print(json_load)

# %%
import json

json_str = '''
{
    "test":"json",
    "test2":"json2"
}
'''

json_dict = json.loads(schools_data)
print(f'{schools_data}')


# %%

#
# Pystanで総計モデリングを行う
#

#
# データの読み込み
#
data = pd.read_csv("/Users/tomoyauchiyama/statisticModel/kubobook_2012/hbm/data7a.csv")
#print(data.describe())

seeds = np.arange(min(data.y), max(data.y)+1, 1)
c = collections.Counter(data.y)

cnt = []
for seed in seeds:
    cnt.append(c[seed])

fig=plt.figure()
ax=fig.add_subplot(111)

ax.scatter(seeds, cnt)
ax.plot(seeds, cnt, linestyle="--", label="observed data")

ax.set_xlabel("# of seed")
ax.set_ylabel("# of plant")
ax.legend(loc="upper left")
ax.set_ylim(0, 25)
ax.legend(loc="upper left")



# %%

#
# 階層ベイズモデルの構築
#

import stan_jupyter as stan

GLMM_model = """
    data {
        int<lower=0> N;
        int y[N];  // observed data
    }
    parameters {
        real beta;
        real s;
        vector[N] r;
    }
    transformed parameters {
        real p[N];
        for (i in 1:N)
            p[i] = inv_logit(beta + r[i]);
    }
    model{
        beta ~ uniform(0, 100);

        s ~ uniform(0, 100);

        for (i in 1:N)
            r[i] ~ normal(0, s);

        y ~  binomial(8, p);
    }
"""

#
# モデリングで用いられる変数にサンプルデータを代入
#
sample_data = {
    'N': len(data.y),
    'y': data.y.tolist()
}

# %%
print(sample_data)

# %%
#
# modelの読み込み
#
posterior = stan.build(GLMM_model, data=sample_data)

# %%
#
# MCMCサンプリングの実行
#

fit = posterior.sample(num_chains=4, num_samples=1000)

# %%
print(fit.to_frame())

# %%

import arviz

arviz.plot_trace(fit)
plt.show()

# %%
results = fit.to_frame()

for i in results.columns:
    print(i)

# %%
results = arviz.summary(fit)



# %%
import stan_jupyter as stan
import pandas as pd
import numpy as np
import collections
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score

# %%
def generate_sampleData(num,seed):
    YY = [] # 応答変数
    XX = [] # 説明変数

    XX_num = 8
    intercept = 0.2
    weight = [0.2, 0.3, 0.5, -0.4, 0.1, 0.2, 0.5, -0.3] # 各特徴量の重み

    np.random.seed(seed=seed)
    for i in range(num):
        X = [np.random.rand() for j in range(XX_num)]
        noise = [np.random.normal(0, 0.1) for j in range(XX_num)]
        Y = sum([intercept + X[j]*weight[j] + noise[j] for j in range(XX_num)])

        YY.append(Y)
        XX.append(X)

    df = pd.DataFrame(np.c_[YY, XX],
                      columns = ['Y', 'X0', 'X1', 'X2', 'X3',
                                 'X4', 'X5', 'X6', 'X7']
                     )

    return df

data = generate_sampleData(1000, 0)

X = data.drop('Y', axis=1)
Y = data['Y']

X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=0)
# %%
#
# ベイズ線形回帰をモデリング　
#

sample_code = """
    data {
        int<lower=0> N;
        int<lower=0> D;
        matrix[N, D] X;
        vector[N] y;
        int<lower=0> N_new;
        matrix[N_new, D] X_new;
    }
    parameters {
        real w0;
        vector[D] w;
        real<lower=0> sigma;
    }
    model {
        for (i in 1:N)
            y[N] ~ normal(w0 + dot_product(X[i], w), sigma);
    }
    generated quantities {
        vector[N_new] y_new;
        for (i in 1:N)
            y_new[i] = normal_rng(w0 + dot_product(X[i], w), sigma);
    }
"""

sample_data = {
    'N': X_train.shape[0],
    'D': X_train.shape[1],
    'X': X_train,
    'y': Y_train,
    'N_new': X_test.shape[0],
    'X_new': X_test
}
# %%
posterior = stan.build(sample_code, data=sample_data)

# %%
print(sample_data)

# %%
import stan_jupyter as stan

schools_code = """
    data {
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

schools_data = {'J': 8,
                'y': [28, 8, -3, 7, -1, 1, 18, 12],
                'sigma': [15, 10, 16, 11, 9, 11, 10, 18]
                }

# %%
# モデルのコンパイル
posterior = stan.build(schools_code, data=schools_data)

#%%
# MCMCサンプリングの実行
fit = posterior.sample(num_chains=4, num_samples=1000)

# %%
# 出力を確認
eta = fit['eta']
print(eta)
df = fit.to_frame()
print(df)

