# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as st

df=pd.read_csv("/Users/tomoyauchiyama/statisticModel/kubobook_2012/poisson/data3a.csv")
print(df.head())
print("\n")
print(df.x) #xの列のみを表示
print("\n")
print(df.info()) #データの型を把握
print("\n")
print(df.describe()) #基本統計量
print("\n")
print(df.f.value_counts()) #サンプル群のn数

# %%
plt.scatter(df.x[df.f=="T"],df.y[df.f=="T"],label="C",marker="o",color="pink")
plt.scatter(df.x[df.f=="C"],df.y[df.f=="C"],label="T",marker="x",color="green")

plt.legend(loc="upper left")
plt.xlabel("size")
plt.ylabel("# of seed")
plt.show()

# %%
plt.boxplot([df.y[df.f=="C"], df.y[df.f=="T"]], labels=["C", "T"])
plt.ylabel("# of seed")
plt.show()

plt.boxplot([df.x[df.f=="C"], df.x[df.f=="T"]], labels=["C", "T"])
plt.ylabel("size")
plt.show()

# %%

#ポアソン回帰の統計モデル
import scipy.stats as st
import statsmodels.api as sm
import statsmodels.formula.api as smf

fit=smf.glm(formula="y~x",data=df,family=sm.families.Poisson())
print(fit.fit().summary())

#ポアソン回帰モデルで可能な最小逸脱度
# <- full model：全ての観測データに対応したパラメータを当てはめたモデル

#最大対数尤度
full_model=sum(st.poisson.logpmf(i,i) for i in df.y)
print(full_model)
#逸脱度
f=-2*full_model
print(f)

#一定モデルのパラメータの推定
fit=smf.glm(formula="y~1",data=df,family=sm.families.Poisson())
print(fit.fit().summary())

#AIC <- -2{(最大対数尤度)-(最尤推定したパラメータの数:k)}
#一定モデルの最大対数尤度（k=1）

#xモデルの対数対数尤度（k=2）

#モデルを複雑化（パラメータの数を増やす）をするだけで、
#観測データへの当てはまりのよさ最大対数尤度は改善されるので、
#AICでモデル選択を行う

#AICは統計モデルの予測のよさを示す平均対数尤度の推定値で、
#最大対数尤度のバイアス補正によって評価される

#予測力の高い統計モデルは、サンプルサイズに依存する
#データ量が少ない時は、真のモデルよりもパラメータの数に少ない方が
#良い予測ができる

# %%
