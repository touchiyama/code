# %%

#2つの統計モデル（単純と複雑モデル）の逸脱度の差に着目した尤度比検定を行う
#一定モデル（単純モデル）とxモデル（複雑モデル）の逸脱度の差を考えた時に
#検定統計量である逸脱度の差から単純モデルに比べて複雑モデルの方が
#より良い予測力を有するかどうかを調べる　（帰無仮説：単純モデル、対立仮説：複雑モデルとする）

import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf

df=pd.read_csv("/Users/tomoyauchiyama/statisticModel/kubobook_2012/poisson/data3a.csv")

fit_x=smf.glm("y~x",data=df,family=sm.families.Poisson()).fit()
fit_cons=smf.glm("y~1",data=df,family=sm.families.Poisson()).fit()
dd=fit_cons.deviance-fit_x.deviance
print(dd)

# %%

#汎用性のあるパラメトリックブートストラップ(BP)法
#サンプルサイズが小さい場合でも、使える
def get_dd12(df, n_bp):
    dd12=[]
    #以下の操作を"n_bp"分繰り返す
    for i in range(n_bp):
        #観測データ（種子数y）に基づくポアソン分布で乱数を個体数分、作成
        df["y_rnd"]=np.random.poisson(df.y.mean(),len(df.y))
        #作成した乱数データに統計モデルを当てはめる
        fit_x=smf.glm("y_rnd~x",data=df,family=sm.families.Poisson()).fit()
        fit_cons=smf.glm("y_rnd~1",data=df,family=sm.families.Poisson()).fit()
        #2つのモデルの逸脱度の差を算出
        dd=fit_cons.deviance-fit_x.deviance
        dd12.append(dd)
    #リスト型からpandasのseries型に変換
    return pd.Series(dd12)

dd12=get_dd12(df, 1000)

print(dd12.describe())

# %%

#BP法による結果→→→近似計算ではない
import matplotlib.pyplot as plt

#ポアソン分布に従う乱数から得られた逸脱度の差における確率分布
plt.hist(dd12, 100)
plt.axvline(x=dd, linestyle='--', color="red")

#乱数は固定していないので、出力値は変化する
p_dd=((dd12 > dd).sum()/len(dd12))
print(p_dd)

#P=0.05となる逸脱度の差
p_05=dd12.quantile(0.95)
plt.axvline(x=p_05, linestyle='--', color="gray")
plt.xlabel("Δdeviance(cons_model-x_model)")
print(p_05)

plt.show()

# %%

#カイ2乗分布を用いた近似計算法
#カイ2乗分布とは、
#確率変数が互いに独立で、標準正規分布N(0,1)に従う時、自由度kのカイ2乗が従う確率分布
#→分布は自由度kに依存する

#カイ2乗分布を用いる場合、
#母分散の区間推定、適応度の検定、独立性の検定
#今回は母分散の区間推定のことを指す?

#サンプルサイズが大きい場合に、有効な近似計算である。
#→ n=100の時には、近似的に得られたP値は正確でない可能性がある

#R→anova(fit_x,fit_cons, test="Chipq")を実装
#＝ばらつきの一種である逸脱度を調べている
#statsmodelsのanova関数はGLMの結果に対応していない

#一定モデルの自由度k=1、xモデルの自由度k=2より、
#自由度1(2-1)のカイ2乗分布で近似する。

fit_cons=smf.glm("y~1",data=df,family=sm.families.Poisson()).fit()
fit_x=smf.glm("y~x",data=df, family=sm.families.Poisson()).fit()
dd12=fit_cons.deviance-fit_x.deviance

#自由度1のカイ2乗分布に従う乱数＝ばらつきの一種である逸脱度の差を近似？
chi=pd.Series(np.random.chisquare(1,1000))
plt.hist(chi, 100)
plt.axvline(x=dd12, linestyle='--', color="red")

p_dd=((chi > dd12).sum())/1000
print(p_dd)

p_05=chi.quantile(0.95)
plt.axvline(x=p_05, linestyle='--', color="gray")
plt.xlabel("Δdeviance(cons_model-x_model)")
print(p_05)

plt.show()

#目的するデータの分布が、ポアソン分布でなく、等分散正正規分布の場合、
#小標本の場合の統計検定量の確率分布が利用できる
#その場合、カイ2乗分布に従う近似を使わない　→　平均の差：t分布、分散比の検定：F分布が用いられる

#%%
