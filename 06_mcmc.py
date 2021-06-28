# %%

#RのglmmMLというパッケージでGLMMの最尤推定を行う
#Pythonには、それに対応したモジュールがないっぽい
#pyperモジュールでRを読み込む
from numpy.core.fromnumeric import size
import pyper

r=pyper.R()
r("library(\"glmmML\")")
r("data <- read.csv(\"/Users/tomoyauchiyama/statisticModel/kubobook_2012/glmm/data.csv\")")
r("head <- head(data)")
print(r.get("head"))
r("summary <- summary(data)")
print(r.get("summary"))

r("""
glmm <- glmmML(cbind(y, N-y) ~ x, data = data, family = binomial, cluster = id)
""")
glmm=r.get("glmm")
print(type(glmm)) #dictionary型で表現される
for key, value in glmm.items():
    print(key, value)

"""
boot 0
converged True
coefficients [-4.18990083  1.0048429 ]     => b1=-4.18990083, b2=1.0048429
coef.sd [0.87768931 0.20746726]
sigma 2.40818441259159 　　　　　　　　　　　　=> sの最尤推定値
sigma.sd 0.220207121126685                 => sの推定値のばらつき
variance [[ 7.70338527e-01 -5.69100843e-04 -1.45622415e-03]
 [-5.69100843e-04  4.30426647e-02  1.18130781e-02]
 [-1.45622415e-03  1.18130781e-02  4.84911762e-02]]
aic 275.41422406129
bootP None
deviance 269.41422406129
df.residual 97
cluster.null.deviance 513.836264179973
cluster.null.df 98
posterior.modes [-1.33717127  0.19525782  0.96712293  2.00676514  0.19525782 -1.33717127
 -1.33717127  3.55622124  0.19525782  2.96238793 -1.33717127 -1.33717127
 -1.33717127 -1.33717127 -1.33717127  0.19525782  1.52627594  2.46931607
  1.52627594  3.55622124  0.06883161 -0.64844365 -1.95068189  2.05046367
 -0.64844365  0.06883161 -1.95068189 -1.95068189 -1.95068189  0.60820467
  1.54625283  1.54625283 -1.95068189 -1.95068189 -1.95068189  3.63912978
  1.08201629  0.06883161 -1.95068189  3.63912978  1.14187599 -0.83434762
 -2.63806929  0.62433193  1.14187599 -0.83434762  2.88426326 -0.83434762
  1.80066545  1.80066545  2.88426326 -2.63806929  0.15699358  1.80066545
 -0.31172889 -1.50661638 -0.31172889 -0.83434762  2.88426326 -2.63806929
  0.93825437 -0.76806875 -1.74167413 -1.23319731 -1.23319731  0.93825437
  2.17690982  0.93825437 -3.37847072  0.93825437  2.17690982  2.17690982
  2.17690982  0.23721234  0.93825437 -0.76806875 -0.29615681 -2.37646571
  2.17690982 -1.74167413 -3.25586934 -3.25586934 -4.15753739  1.53551116
  1.53551116  1.53551116  1.53551116  1.53551116  1.53551116  1.53551116
  1.53551116 -4.15753739 -0.6628297  -2.15590396 -2.15590396  0.08927529
 -1.21489937  0.08927529 -2.65252652  1.53551116]
prior gaussian
terms ['~' 'cbind(y, N - y)' 'x']
info 0
call ['glmmML' 'cbind(y, N - y) ~ x' 'binomial' 'data' 'id']
"""
#Rでの出力（）
#Call:  glmmML(formula = cbind(y, N - y) ~ x, family = binomial, data = data,      cluster = id)

#             coef se(coef)      z Pr(>|z|)
#(Intercept) -4.190   0.8777 -4.774 1.81e-06
#x            1.005   0.2075  4.843 1.28e-06

#Scale parameter in mixing distribution:  2.408 gaussian => sの最尤推定値
#Std. Error:                              0.2202         => sの推定値のばらつき

#       LR p-value for H_0: sigma = 0:  2.136e-55

#Residual deviance: 269.4 on 97 degrees of freedom 	AIC: 275.4

print(glmm["coefficients"][0]) #np.array型

# %%
import numpy as np
import matplotlib.pyplot as plt

#glmmによる予測モデル ランダム効果項を生成する正規分布の平均値は０に近い値をとるのか?
def logistic(x):
    z=glmm["coefficients"][0]+glmm["coefficients"][1]*x
    return 1/(1+np.exp((-1)*z))

df=r.get("data")
plt.scatter(df.iloc[0:,2],df.iloc[0:,1],marker=".")
plt.plot(df.iloc[0:,2],logistic(df.iloc[0:,2])*8)
plt.xlabel("# of seed")
plt.ylabel("# of leaf")
plt.show()

# %%
import scipy.stats as st

#真の葉数依存性モデル
def logistic_rand(x):
    effect=st.norm.pdf(x,3.81,glmm["sigma"]) #生存種子数の平均値=3.81,正規分布パラメータs=2.408
    z=glmm["coefficients"][0]+glmm["coefficients"][1]*x+effect
    return 1/(1+np.exp((-1)*z))

#滑らかな曲線を描く
x=np.arange(df.iloc[0:,2].min(),df.iloc[0:,2].max(),(df.iloc[0:,2].max()-df.iloc[0:,2].min())/1000)
plt.scatter(df.iloc[0:,2],df.iloc[0:,1],s=20, alpha=0.3)
plt.plot(x,logistic(x)*8,color="black",label="glmm model")
plt.plot(x,logistic_rand(x)*8,color="red",label="real model")
plt.xlabel("# of seed")
plt.ylabel("# of leaf")
plt.legend(loc="upper left")
plt.show()

# %%
import collections

x4=df.iloc[0:,1][df.iloc[0:,2]==4].tolist()
c=collections.Counter(x4)
y,cnt=map(list, zip(*c.most_common()))

def logistic_rand(x):
    effect=st.norm.pdf(x,4,glmm["sigma"]) #生存種子数の平均値=3.81,正規分布パラメータs=2.408
    z=glmm["coefficients"][0]+glmm["coefficients"][1]*x+effect
    return 1/(1+np.exp((-1)*z))

glmm_binom_pmf=[(st.binom.pmf(n=8, p=logistic(4), k=i))*sum(cnt) for i in cnt]
real_binom_pmf=[(st.binom.pmf(n=8, p=logistic_rand(4), k=i))*sum(cnt) for i in cnt]

fig=plt.figure()
ax=fig.add_subplot(111)
ax.scatter(y,cnt,label="observed data")
ax.scatter(y,real_binom_pmf,label="real model",marker="x")
ax.scatter(y,glmm_binom_pmf,label="glmm model",marker="^")
ax.set_xlabel("# of seed")
ax.set_ylabel("# of plant")
ax.set_ylim(0,6)

# %%

"""
#rand_norm=[np.random.normal(i,0,int(glmm["sigma"])) for i in df.iloc[0:,1]]
for r in rand_norm:
    z=glmm["coefficients"][0]+glmm["coefficients"][1]*x+r
    q=1
"""
# %%
