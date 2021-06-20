# %%
from matplotlib import colors
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

#xモデルのパラメータ推定
fit_x=smf.glm(formula="y~x",data=df,family=sm.families.Poisson())
print(fit_x.fit().summary())

#ポアソン回帰モデルで可能な最小逸脱度
#-> full model：全ての観測データに対応したパラメータを当てはめたモデル
#最大対数尤度
full_model=np.sum(st.poisson.logpmf(i,i) for i in df.y)
print(f"AIC_full: {full_model}")
#逸脱度
f=-2*full_model
print(f"Deviance: {f}")

#一定モデルのパラメータの推定
fit_cons=smf.glm(formula="y~1",data=df,family=sm.families.Poisson())
print(fit_cons.fit().summary())


#AIC <- -2{(最大対数尤度)-(最尤推定したパラメータの数:k)}
#一定モデルの最大対数尤度（k=1）
AIC_cons=-2*(fit_cons.fit().llf-1)
#xモデルの対数対数尤度（k=2）
AIC_x=-2*(fit_x.fit().llf-2)

print(f"AIC_cons: {AIC_cons}, AIC_x: {AIC_x}")

#モデルを複雑化（パラメータの数を増やす）をするだけで、
#観測データへの当てはまりのよさ最大対数尤度は改善されるので、
#AICでモデル選択を行う

#AICは統計モデルの予測のよさを示す平均対数尤度の推定値で、
#最大対数尤度のバイアス補正によって評価される

#予測力の高い統計モデルは、サンプルサイズに依存する
#データ量が少ない時は、真のモデルよりもパラメータの数に少ない方が
#良い予測ができる


# %%

#一定モデルの最大対数尤度と平均対数尤度の関係（バイアスの平均値を求めるため）
x=np.arange(1,200,1)
max_logll=[]
ave_logll=[]
for i in x:
    #真のモデルをλ=8のポアソン分布モデルから生成される50個のデータと仮定して、乱数を発生
    data_200=[np.random.poisson(8,50) for j in x]
    sum_logL=0
    for set in data_200:
        max_logL=np.sum(st.poisson.logpmf(set,np.exp(fit_cons.fit().params["Intercept"])))
        sum_logL=sum_logL+max_logL
        if i == 1:
            max_logll.append(max_logL)
    ave=sum_logL/len(data_200)
    ave_logll.append(ave)

plt.scatter(x,max_logll,color="green",marker="o",label="logL*")
plt.scatter(x,ave_logll,color="red",marker="x",label="E(logL)")
plt.legend(loc="upper right")
plt.ylabel("Log-Likelihood")
plt.xlabel("Data sets")
plt.show()

# %%

#一定モデルのバイアス
#バイアスb=logL*-E(logL)
#b=max_logll-ave_logll
b=np.array(max_logll)-np.array(ave_logll)
plt.boxplot(b)
print(b.mean())

# %%

#一定モデルの最大対数尤度と平均対数尤度の関係（バイアスの平均値を求めるため）
max_logll_cons=max_logll
ave_logll_cons=ave_logll

#xモデル
x=np.arange(1,200,1)
max_logll=[]
ave_logll=[]
for i in x:
    data_200=[np.random.poisson(8,50) for j in x]
    sum_logL=0
    for set in data_200:
        linear=fit_x.fit().params["Intercept"]+(fit_x.fit().params["x"]*set)
        max_logL=np.sum(st.poisson.logpmf(set,np.exp(linear)))
        sum_logL=sum_logL+max_logL
        if i == 1:
            max_logll.append(max_logL)
    ave=sum_logL/len(data_200)
    ave_logll.append(ave)

plt.scatter(x,max_logll,color="green",marker="o",label="logL*")
plt.scatter(x,ave_logll,color="red",marker="x",label="E(logL)")
plt.legend(loc="upper right")
plt.ylabel("Log-Likelihood")
plt.xlabel("Data sets")
plt.show()

max_logll_x=max_logll
ave_logll_x=ave_logll

# %%

#実際のデータ解析では、真の統計モデルが不明なので、E(logL)の値は未知
#バイアスb=logL*-E(logL)より、E(logL)=logL*-b
#bの平均値とlogL*がわかれば、E(logL)の推定値が得られそう
#バイアス補正という

#xモデル(k=2)と一定モデル(k=1)の最大対数尤度と平均対数尤度の比較
diff_max=np.array(max_logll_x)-np.array(max_logll_cons)
print(diff_max.mean())
diff_ave=np.array(ave_logll_x)-np.array(ave_logll_cons)
print(diff_ave.mean())
diff_b=diff_max-diff_ave

plt.boxplot(diff_b)
print(diff_b.mean())

# %%

#test
data_200=[np.random.poisson(8,50) for j in x]
for set in data_200:
    linear=fit_x.fit().params["Intercept"]+(fit_x.fit().params["x"]*set)
    max_logL=np.sum(st.poisson.logpmf(set,np.exp(linear)))
    print(max_logL)
# %%
