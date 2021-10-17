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
