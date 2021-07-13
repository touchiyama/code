# %%

#pymc3でベイス統計モデリングを行う
import pymc3 as mc
import pandas as pd
import numpy as np
from pymc3 import data
from pymc3.distributions.discrete import Poisson

"""
pymc3.sampling.sample(draws=500, step=None, init='auto', n_init=200000,
                    start=None, trace=None, chain=0, njobs=1, tune=500,
                    nuts_kwargs=None, step_kwargs=None, progressbar=True, model=None,
                    random_seed=-1, live_plot=False, discard_tuned_samples=True,
                    live_plot_kwargs=None, **kwargs
                    )
"""

"""
Pymc3の使い方の説明
https://pymc3-testing.readthedocs.io/en/rtd-docs/api/inference.html#step-methods
基本的な引数について
・draws: 最初の引数で、サンプル数を指定
・steps: MCMCサンプリングアルゴリズムを指定する
        (Metropolis, HMC, NUTS, default: NUTS）
・start: パラメータの初期値を与える
        (その時のパラメータの推定法(点推定): MAP推定法 => find_MAP(method='Powell’)）
・tune:　先頭からサンプル数を指定
        (サンプルの最初は、ランダムに選ばれた初期値の影響を大きく受けるため、
        捨てた方がよい => startsで初期値を設定したらどうなるか？）
・njobs(chain): njobsが２以上ならば、chainがnjobに応じて、その数分立ち上がる。
・random_seed: njobsが２以上ならば、リストが作成される
        (njobsに対応して、random_seedも設定)
"""


# %%
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt

np.random.seed(10) #乱数の固定

n = 50
p1 = st.bernoulli.rvs(0.50, size=n)
p2 = st.bernoulli.rvs(0.55, size=n)

h1 = sum(p1)
h2 = sum(p2)

print(f"likelihood_p1(h1/n)={h1/n}")
print(f"likelohood_p2(h1/n)={h2/n}")

prob=np.arange(0, 1, 0.01)
ll_1=[sum(st.bernoulli.logpmf(p1,p=i)) for i in prob]
ll_2=[sum(st.bernoulli.logpmf(p2,p=i)) for i in prob]

plt.plot(prob, ll_1)
plt.plot(prob, ll_2)

plt.axvline(x=prob[np.argmax(ll_1)],linestyle="--")
plt.axvline(x=prob[np.argmax(ll_2)],linestyle="--")

print(f"Loglikelihood_p1(logsum)={prob[np.argmax(ll_1)]}")
print(f"Loglikelohood_p2(logsum)={prob[np.argmax(ll_2)]}")

# %%

#pymc3によるモデリング
# パラメータpの事前分布は一様分布に従うと仮定
# ベルヌーイ分布の定義

import pymc3 as mc

with mc.Model() as model:
    p1 = mc.Uniform("p1", lower=0.0, upper=1.0)
    y1 = mc.Binomial("y1", n=n, p=p1, observed=h1)
    p2 = mc.Uniform("p2", lower=0.0, upper=1.0)
    y2 = mc.Binomial("y2", n=n, p=p2, observed=h2)
    delta = mc.Deterministic('delta', p2-p1)

#pymc3は環境を正しく構築しないと働かない
#バージョンが違うだけで、Cのコンパイラがが動かない

# %%

#モデルの可視化 => graphvizのライブラリーが必要
#anacondaでインストールする必要がある。
#conda install python-graphviz
mc.model_to_graphviz(model)

# %%

#MCMCサンプリングの実行
#'find_MAP should not be used to initialize the NUTS sampler,
# simply call pymc3.sample() and
# it will automatically initialize NUTS in a better way.'
# find_MapはNUTSの初期値の決定で使うことを推奨していない
# pm.sampleは自動的にNUTSの初期値を探している
with model:
        p_start = mc.find_MAP(method="Powell")
        trace = mc.sample(2000, tune=500, step=mc.NUTS())

# %%
#結果の出力
mc.summary(trace)

# %%
#MCMCサンプリングの結果を可視化
mc.traceplot(trace)

#事後分布の表示
mc.plot_posterior(trace)

# %%
import scipy.stats as st
import matplotlib.pyplot as plt
import pandas as pd

df=pd.read_csv("/Users/tomoyauchiyama/statisticModel/kubobook_2012/gamma/data06.csv")
print(df.head())

plt.scatter(df.x, df.y)

#花の質量yについてガンマ分布でモデルを作る
#花の質量y（応答変数）、葉の質量x（説明変数）とすると、
#花の質量は、正の値をとり、連続値をとるのでガンマ分布に従うと仮定できる
#ここで、y=A*x**bに従うと仮定すると
#この時、A＝exp(a)とおくと、y=exp(a)*x**b →　y=exp(a+blogx)と表せる
#これより、logy=a+blogx

#scipy.statsのppf (Percent point function) パーセント点関数で
#ガンマ分布 => 平均=shape/rate、scale=1/rate、分散=shape/rate**2
#ppfで必須のパラメーター：shapeとscare
#ここで、平均=exp(a+b*logx)と表されるので、rate=shape/平均=shape/exp(a+b*logx)と表せる
#これより、scale=1/rate=exp(a+b*logx)と表せる

#ガンマ分布のパラメータを
#scipy.statsのfitメソッドで求めると,shape=1.9172130853026834
shape, loc, scale=st.gamma.fit(df.y)
print(shape, loc, scale)

# %%
import pymc3 as mc
import numpy as np

#ガンマ分布のパラメータをpymc3で求める
with mc.Model() as model:
        a=mc.Normal("a", mu=0, sd=100)
        b=mc.Normal("b", mu=0, sd=100)

        shape=mc.Uniform("shape", lower=0.001, upper=1000)
        #Deterministicメソッドは
        #関連する他の変数から（確定的に求められる）変数を導出するために用いる
        rate=mc.Deterministic("rate", shape/np.exp(a+b*np.log(df.x)))

        y=mc.Gamma("y", alpha=shape, beta=rate, observed=df.y)

#構築したモデルが不適ならば、以下のような警告が出る
#UserWarning: The variable specified for alpha has negative support for Gamma,
#likely making it unsuitable for this parameter
#ガンマ分布のパラメータshapeは常に正の値をとる
mc.model_to_graphviz(model)

# %%

with model:
        draws=2000
        p_start=mc.find_MAP()
        steps=mc.NUTS()
        #steps=mc.HamiltonianMC()
        tune=500
        njobs=3
        random_seed=10

        trace=mc.sample(draws=draws, start=p_start, steps=steps, tune=tune)

# %%
mc.summary(trace)
# %%
mc.traceplot(trace)

# %%

print(trace["shape"].mean())
print(trace["a"].mean())
print(trace["b"].mean())
print(trace["rate"].mean())

# %%

#pymc3で推定したパラメータを用いて、統計モデルの可視化
def pred_range(per):
    pred=[]
    for x in df.x:
        a = trace["shape"].mean()
        rate = a/np.exp(trace["a"].mean()+trace["b"].mean()*(np.log(x)))
        scale=1/rate
        lower, upper=st.gamma.ppf([per,1-per], a=a, scale=scale)
        pred.append((lower,upper))

    return pred

plt.scatter(df.x,df.y,label="weight of flower")
x=np.arange(df.x.min(),df.x.max(),(df.x.max()-df.x.min())/100)
def gamma_model(x):
    return(np.exp(trace["a"].mean()+trace["b"].mean()*(np.log(x))))

plt.plot(x,gamma_model(x),color="black",label="mean")

#plt.plot(x,gamma_model2(x),color="gray",linestyle="--",label="2.5%")
#plt.plot(x,gamma_model3(x),color="gray",linestyle="--",label="97.5%")

y1, y2 = map(list, zip(*pred_range(0.25)))
plt.fill_between(x=df.x, y1=y1, y2=y2, alpha=0.2, color="blue")
y3, y4 = map(list, zip(*pred_range(0.05)))
plt.fill_between(x=df.x, y1=y3, y2=y4, alpha=0.1, color="green")

plt.xlabel("weight of leaf")
plt.ylabel("weight of flower")
plt.legend(loc="upper left")
plt.show()

# %%
import pyper
import pandas as pd

r = pyper.R(use_pandas="True")
r("load(\"/Users/tomoyauchiyama/statisticModel/kubobook_2012/gibbs/d.RData\")")
df = r.get("d")
df.columns = df.columns.str.strip() # DataFrameのカラム名に挿入された空白文字を削除
print(df.columns.values) #DataFrameのカラム名を列挙

# %%
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf

poisson_fit = smf.glm("y~x", data=df, family=sm.families.Poisson()).fit()
print(poisson_fit.summary())

df["predict"] = poisson_fit.predict()
#print(df.head())

plt.scatter(df.x, df.y)
plt.plot(df.x, df.predict, linestyle="--", color="gray")
plt.ylabel("# of seeds")
plt.xlabel("size")
plt.show()

# %%
import pymc3 as mc
import numpy as np

#種子数が植物の体サイズに依存する統計モデルの構築を考える
#poisson分布に従うモデル => 平均種子数logλ=β1+β2*x

#overflowを起こして、計算できない
with mc.Model() as model:
        beta1 = mc.Normal("beta1", mu=0, sd=100)
        beta2 = mc.Normal("beta2", mu=0, sd=100)
        #link = mc.Deterministic("link", beta1+beta2*df.x)
        xx = np.array(df.x, dtype="float64")
        yy = np.array(df.y, dtype="float64")
        link = beta1+beta2*xx

        y = mc.Poisson("y", mu=np.exp(link), observed=yy)


# %%
mc.model_to_graphviz(model)
# %%
with model:
    # 101個目から3個置きでサンプルを取得するチェインを3つ作る
    # NUTSではburnとthinが効いていない？
    trace = mc.sample(1600, step=mc.HamiltonianMC(), burn=100, thin=100, njobs=3, random_seed=0)[::3]
# %%
mc.summary(trace)
# %%
mc.traceplot(trace)
# %%
with model:
        trace = mc.sample(1600, step=mc.Metropolis(), burn=100, thin=100, njobs=3, random_seed=0)[::3]
# %%
mc.summary(trace)
# %%
mc.traceplot(trace)

#MCMCサンプリングの手法
#NUTS > HamiltonianMC > Metropolis

# %%

#MCMCサンプリング過程で得られた全パラメータβ1とβ2を用いて、モデリングを可視化

def MCMC_model(b1, b2, x):
        return(np.exp(b1+b2*x))

def poisson_model(x):
    return(np.exp(trace["beta1"].mean()+trace["beta2"].mean()*x))

plt.scatter(df.x, df.y)

for b1, b2 in zip(trace["beta1"], trace["beta2"]):
        x = np.arange(min(df.x), max(df.x), (max(df.x)-min(df.x))/100)
        plt.plot(x, MCMC_model(b1, b2, x), linestyle="--", color="gray", alpha=0.03)

x = np.arange(min(df.x), max(df.x), (max(df.x)-min(df.x))/100)
plt.plot(x, poisson_model(x), color="red", label="mean")
plt.ylabel("# of seeds")
plt.xlabel("size")
plt.legend()
plt.show()

# %%
import collections

data7a = pd.read_csv("/Users/tomoyauchiyama/statisticModel/kubobook_2012/hbm/data7a.csv")

c=collections.Counter(data7a.y)
cnt = []
yy = np.arange(min(data7a.y), max(data7a.y)+1, 1)
for y in yy:
        cnt.append(c[y])

#y,cnt=map(list, zip(*c.most_common()))

fig=plt.figure()
ax=fig.add_subplot(111)
ax.scatter(yy,cnt,label="observed data")
ax.plot(yy, cnt, linestyle="--")
ax.set_xlabel("# of seed")
ax.set_ylabel("# of plant")
ax.set_ylim(0, max(cnt)+5)
ax.legend(loc="upper left")

#pymc3で個体差を組み込んだ二項分布モデルを作る


# %%
