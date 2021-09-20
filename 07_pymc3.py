# %%

#pymc3でベイス統計モデリングを行う
from matplotlib import colors
from numpy.core.fromnumeric import mean, size
from numpy.lib.function_base import average
import pymc3 as mc
import pandas as pd
import numpy as np

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
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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

# %%

#pymc3で個体差を組み込んだ二項分布モデルを作る
#二項分布の確率 q = 1 / (1 + exp(-z))
#z = beta + r

#log(q/1-q) = z
#logit(q) = log(q/1-q) を指し、logit(q) = z と表せる。
import pymc3 as mc

with mc.Model() as binom_MC_model:
        beta = mc.Normal("beta", mu=0, sd=100)

        #全100個体分の正規分布を求める
        s = mc.Uniform("s", lower=0, upper=10000)
        r = mc.Normal("r", mu=0, sd=s, shape=len(data7a.y)) #shapeの意味を調べる

        z = beta + r
        q = mc.Deterministic("q", 1/(1+np.exp((-1)*z)))
        y = mc.Binomial("y", n=8, p=q, observed=data7a.y)

# %%
mc.model_to_graphviz(binom_MC_model)

# %%
with binom_MC_model:
    # NUTSで101個目からサンプルを取得するチェインを3つ作る
    trace = mc.sample(1500, step=mc.NUTS(), tune=100, njobs=3, random_seed=0)

# %%

mc.summary(trace, varnames=['beta', 's'])

# %%

mc.traceplot(trace, varnames=['beta', 's', 'r'])

# %%

#全ペア(1500x3=4500個）の{beta,s}で、
#100個体分の線形予測子z=bata+rを求める

#まず、平均0、各標準偏差s（全1500x3=4500個）から100個体分の個体差rを正規分布で生成する
#全ペアの{beta,r}で線形予測子を求め、それに対応する生存確率qを求める
#生存確率q=1/(1+exp(-z))、ロジスティック回帰を求める

#100個の各生存確率に対応した二項分布を求める

#イメージ：[(1)[1,2,..,100], (2)[1,2,..,100],..., (4500)[1,2,..,100]]

prob = []
for beta, ri in zip(trace['beta'], trace['r']):
        #ri = np.random.normal(loc=0, scale=sd)#xの範囲がわからない
        zi = beta + ri
        qi = 1/(1+np.exp((-1)*zi))
        prob.append(qi)

y_pred = []
for pi in np.array(prob):
        yi = np.random.binomial(n=8, p=pi) #二項分布に従う乱数を発生
        y_pred.append(yi)

#4500サンプル・100個体分の生存確率の平均値を計算する
#print(len(np.array(y))) <- 4500サンプル
#print(len(np.array(y[0]))) <- 100個体


# %%

plt.hist(np.array(y[0]))
plt.hist(np.array(y[1]))

#上記2サンプルのヒストグラムの分布を確認すると、
#パラメータの事後分布に基づいて統計モデルを可視化する方法をうまく探せそう

# %%

fig=plt.figure()
ax=fig.add_subplot(111)

yy = np.arange(min(data7a.y), max(data7a.y)+1, 1)

for yi in y_pred:
        c=collections.Counter(yi)
        cnt = []
        for y in yy:
                cnt.append(c[y])

        ax.scatter(yy, cnt, color="gray", alpha=0.03)
        ax.plot(yy, cnt, linestyle="--", color="gray", alpha=0.03)


data7a = pd.read_csv("/Users/tomoyauchiyama/statisticModel/kubobook_2012/hbm/data7a.csv")
c=collections.Counter(data7a.y)
cnt = []
for y in yy:
        cnt.append(c[y])

ax.scatter(yy, cnt, label="observed data")
ax.plot(yy, cnt, linestyle="--")

cnt = []
for yi in y_pred:
        c += collections.Counter(yi)

for y in yy:
        ave=c[y]/(4500)
        cnt.append(ave)

ax.scatter(yy, cnt, label="predict data")
ax.plot(yy, cnt, linestyle="--")

ax.set_xlabel("# of seed")
ax.set_ylabel("# of plant")
ax.legend(loc="upper left")
ax.set_ylim(0, 35)
ax.legend(loc="upper left")

#%%

fig=plt.figure()
ax=fig.add_subplot(111)

yy = np.arange(min(data7a.y), max(data7a.y)+1, 1)

data7a = pd.read_csv("/Users/tomoyauchiyama/statisticModel/kubobook_2012/hbm/data7a.csv")
c=collections.Counter(data7a.y)
cnt = []
for y in yy:
        cnt.append(c[y])

ax.scatter(yy, cnt, label="observed data")
ax.plot(yy, cnt, linestyle="--")

cnt = []
for yi in y_pred:
        c += collections.Counter(yi)

for y in yy:
        ave=c[y]/(4500)
        cnt.append(ave)

ax.scatter(yy, cnt, label="predict data")
ax.plot(yy, cnt, linestyle="--")

ax.set_xlabel("# of seed")
ax.set_ylabel("# of plant")
ax.legend(loc="upper left")
ax.set_ylim(0, 30)
ax.legend(loc="upper left")

# %%

#pymc3で個体差と場所差を組み込んだポアソン分布モデルを作る（階層ベイズモデル）
import matplotlib.pyplot as plt
import matplotlib.collections as mc
import pandas as pd

d1 = pd.read_csv("/Users/tomoyauchiyama/statisticModel/kubobook_2012/hbm/nested/d1.csv")
#print(d1.head())
#print(d1.tail())
#print(d1.y)

#データフレームで、複数条件で行を抽出
#plt.scatter(d1.id[(d1.pot=="A") & (d1.f=="C")], d1.y[(d1.pot=="A") & (d1.f=="C")])

#大文字のアルファベット26文字全て羅列する
#print([chr(ord("A")+i) for i in range(26)])

fig=plt.figure()
ax=fig.add_subplot(111)

for i in range(10):
        pot = chr(ord("A")+i)
        if i <= 4:
                label=pot + "-C"
        else:
                label=pot + "-T"
        ax.scatter(d1.id[d1.pot==pot], d1.y[d1.pot==pot], label=label)

"""
line_1 = [(1, d1.y[d1.f=="C"].mean()), (50, d1.y[d1.f=="C"].mean())]
line_2 = [(51, d1.y[d1.f=="T"].mean()), (100, d1.y[d1.f=="T"].mean())]
line_1_2 = mc.LineCollection([line_1, line_2],
                              colors=["#1f77b4", "#ff7f0e"],
                              linestyle=["--", "--"],
                              label=["C", "T"])
ax.add_collection(line_1_2)
"""

x1=np.arange(1, 51, 1)
y1=[d1.y[d1.f=="C"].mean() for i in range(50)]
ax.plot(x1, y1, color="#1f77b4", linestyle="--", label="C_mean")

x2=np.arange(51, 101, 1)
y2=[d1.y[d1.f=="T"].mean() for i in range(50)]
ax.plot(x2, y2, color="#ff7f0e", linestyle="--", label="T_mean")

ax.set_xlabel("plant ID")
ax.set_ylabel("# of seed")
ax.set_ylim(0, max(d1.y)+5)
ax.set_xlim(0, max(d1.id)+35)
ax.legend(loc="upper right")


# %%

fig=plt.figure()
ax=fig.add_subplot(111)

label = [chr(ord("A")+i) for i in range(10)]
bp = ax.boxplot([d1.y[d1.pot=="A"], d1.y[d1.pot=="B"], d1.y[d1.pot=="C"],
                 d1.y[d1.pot=="D"], d1.y[d1.pot=="E"], d1.y[d1.pot=="F"],
                 d1.y[d1.pot=="G"], d1.y[d1.pot=="H"], d1.y[d1.pot=="I"],
                 d1.y[d1.pot=="J"]], labels=label, patch_artist=True,
                 medianprops=dict(color='black', linewidth=1)
                )

# 細かい設定をできるようにする
# patch_artist=True

color = []
for i in range(10):
        if i <= 4:
                color.append("#1f77b4")
        else:
                color.append("#ff7f0e")

for b, c in zip(bp['boxes'], color):
        b.set(color='black', linewidth=1)  # boxの外枠の色
        b.set_facecolor(c) # boxの色

ax.set_xlabel("pot")
ax.set_ylabel("# of seed")
ax.set_ylim(0, max(d1.y)+5)


# %%

#d1データフレームの加工
#d1.f=="C"の時、0
#d1.f=="T"の時、1

d1.f[d1.f == 'C'] = int(0)
d1.f[d1.f == 'T'] = int(1)
d1.f = d1.f.astype(int)

# %%

#d1データフレームの加工
#A -> 1, B -> 2, C -> 3
alphabet = {}
for i in range(65, 65+10):
        alphabet[chr(i)] = i - 65


for i in range(len(d1)):
        d1.iloc[i, 1] = alphabet[d1.iloc[i, 1]]

d1.pot = d1.pot.astype(int)

# %%


print(d1.dtypes)

print(d1.head())
print(d1.tail())
print(d1.iloc[49,:])
print(d1.iloc[50,:])

# %%

#再帰関数使用時に発生したエラーの対策
#RecursionError: maximum recursion depth exceeded in comparison
#再帰上限の変更（デフォルト上限:10000）

import sys

print(sys.getrecursionlimit())
sys.setrecursionlimit(100000) # RecursionError対策

# 再帰上限を上げたら計算できない状況になった　→ pystanでは実行できるのだろうか？
# pymc3とpystanで比較検討を行ってみる　

# %%
print(sys.getrecursionlimit())

# %%
print(len(d1.id))

# %%
import pymc3 as mc

#logλ=b1+b2x+ri+rj

with mc.Model() as pois_MC_model:
        beta1 = mc.Normal("beta1", mu=0, sd=100)
        beta2 = mc.Normal("beta2", mu=0, sd=100)

        s = mc.Uniform("s", lower=0, upper=10000)
        r = mc.Normal("r", mu=0, sd=s, shape=len(d1.id))

        sp = mc.Uniform("sp", lower=0, upper=10000)
        rp = mc.Normal("rp", mu=0, sd=sp, shape=len(d1.pot.unique()))

        link = beta1 + beta2 * d1.f + r + rp[d1.pot]
        lamda = mc.Deterministic("lamda", np.exp(link))
        y = mc.Poisson("y", mu=lamda, observed=d1.y)

#%%
print(lamda)

# %%

mc.model_to_graphviz(pois_MC_model)


# %%

with pois_MC_model:
        draws=2000
        p_start=mc.find_MAP()
        steps=mc.NUTS()
        #steps=mc.HamiltonianMC()
        tune=500
        njobs=3
        random_seed=10

        trace=mc.sample(draws=draws, start=p_start, steps=steps, tune=tune, cores=1)

# %%

mc.summary(trace, varnames=['rp'])

# %%

import numpy as np

lamda_sum = np.zeros(len(d1.y))

for li in trace['lamda']:
        lamda_sum +=li

lamda_mean = lamda_sum / len(trace['lamda'])

print(lamda_mean)

# %%

mc.summary(trace, varnames=['lamda'])

# %%

d1 = pd.read_csv("/Users/tomoyauchiyama/statisticModel/kubobook_2012/hbm/nested/d1.csv")

np.random.seed(100)
id = np.arange(min(d1.id), max(d1.id)+1)
pred_y = np.random.poisson(lam=lamda_mean)

fig=plt.figure()
ax=fig.add_subplot(111)

for i in range(10):
        pot = chr(ord("A")+i)
        if i <= 4:
                label=pot + "-C"
        else:
                label=pot + "-T"
        ax.scatter(d1.id[d1.pot==pot], d1.y[d1.pot==pot], label=label)

x1=np.arange(1, 51, 1)
y1=[d1.y[d1.f=="C"].mean() for i in range(50)]
ax.plot(x1, y1, color="#1f77b4", linestyle="--", label="C_mean")

x2=np.arange(51, 101, 1)
y2=[d1.y[d1.f=="T"].mean() for i in range(50)]
ax.plot(x2, y2, color="#ff7f0e", linestyle="--", label="T_mean")

ax.scatter(id, pred_y, label="observed data", color="gray", alpha=0.5)
ax.plot(id, pred_y, linestyle="--", color="gray", alpha=0.5)

ax.set_xlabel("# of seed")
ax.set_ylabel("plant ID")
#ax.set_ylim(0, max(y_pred)+5)
ax.set_ylim(0, max(d1.y)+5)
ax.set_xlim(0, max(d1.id)+2)
ax.legend(bbox_to_anchor=(1.4, 1), loc="upper right")
#ax.legend(loc="upper left")


# %%
plt.scatter(d1.y, pred_y)

# %%
print(y_pred)


label1 = [chr(ord("A")+i) + 'o' for i in range(10)]
label2 = [chr(ord("A")+i) + 'p' for i in range(10)]
label = label1 + label2
sort_label = sorted(label)

print(label)

print(sort_label)

# %%

#
# 新規データフレームの構築
#

# データフレームの初期化
d2 = pd.DataFrame()
for i in range(len(pred_y)):
        idx = i % 10
        label = 'p_' + chr(ord('A') + idx)
        col = pd.Series([i+1, label, pred_y[i]])
        d2 = d2.append(col, ignore_index=True)

d2.columns = ['id', 'pot', 'y']
d2.id = d2.id.astype(int)
d2.pot = d2.pot.astype(str)
d2.y = d2.y.astype(int)

# %%

#
# 2つのデータフレームの連結
#

# 連結する前に配列の形を整えた
d3 = pd.concat([d1.iloc[:, [0, 1, 3]], d2])
print(d3)

# %%

fig=plt.figure()
ax=fig.add_subplot(111)

label1 = [chr(ord("A")+i) + 'o' for i in range(10)]
label2 = [chr(ord("A")+i) + 'p' for i in range(10)]
label = label1 + label2
sort_label = sorted(label)

bp = ax.boxplot([d3.y[d3.pot=="A"], d3.y[d3.pot=="p_A"], d3.y[d3.pot=="B"],
                 d3.y[d3.pot=="p_B"], d3.y[d3.pot=="C"], d3.y[d3.pot=="p_C"],
                 d3.y[d3.pot=="D"], d3.y[d3.pot=="p_D"], d3.y[d3.pot=="E"],
                 d3.y[d3.pot=="p_E"], d3.y[d3.pot=="F"], d3.y[d3.pot=="p_F"],
                 d3.y[d3.pot=="G"], d3.y[d3.pot=="p_G"], d3.y[d3.pot=="H"],
                 d3.y[d3.pot=="p_H"], d3.y[d3.pot=="I"], d3.y[d3.pot=="p_I"],
                 d3.y[d3.pot=="J"], d3.y[d3.pot=="p_J"]],labels=sort_label, patch_artist=True, medianprops=dict(color='black', linewidth=1)
                )

# 細かい設定をできるようにする
# patch_artist=True

color = []
for i in range(len(label)):
        if i < int(len(label)/2):
                color.append("#1f77b4")
        else:
                color.append("#ff7f0e")

for b, c in zip(bp['boxes'], color):
        b.set(color='black', linewidth=1)  # boxの外枠の色
        b.set_facecolor(c) # boxの色


ax.set_xlabel("pot")
ax.set_ylabel("# of seed")
ax.set_ylim(0, max(d1.y)+5)

# 総括
# observed dataとpredict dataの比較について
# 個体差間のプロットでは、predict dataがobserved dataに
# よく当てはまっていたことがわかった（過学習？）
# 場所差間のプロットでは、両者の分布に差が見られるように思えた


# %%
