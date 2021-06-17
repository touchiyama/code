# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels import formula

df=pd.read_csv("./03_poisson.csv")
print(df.head())
# %%
#種子数の箱ひげ図
plt.boxplot([df.y[df.f=="C"],df.y[df.f=="T"]],labels=["C","T"])
plt.ylabel("# of seed")
plt.show()

#体サイズの箱ひげ図
plt.boxplot([df.x[df.f=="C"],df.x[df.f=="T"]],labels=["C","T"])
plt.ylabel("size")
plt.show()

# %%
#説明変数を体サイズx、応答変数を種子数yとした時のポアソン回帰モデル
#種子数はカウントサイズで、離散値である。種子数の域値は0~∞
#y=exp(b1+b2x) 
#logy=b1+b2x
x=np.arange(0,4,0.01)
y1=np.exp(-2-0.8*x)
y2=np.exp(-1+0.4*x)

plt.plot(x,y1,label="{-2,-0.8}")
plt.plot(x,y2,label="{-1,0.4}")
plt.legend()
plt.xlabel("size")
plt.ylabel("mean of seed")
plt.show()

# %%
#pythonの有名なglmモジュール →　scikit-learn、statsmodel
#statmodel→Rのglmに似ている
import statsmodels.api as sm
import statsmodels.formula.api as smf

#モデルの当てはめと当てはまりの良さを評価
#fomulaの指定→応答変数~説明変数
fit1=smf.glm(formula='y~x',data=df,family=sm.families.Poisson())
fit2=smf.glm(formula='y~x',data=df,family=sm.families.Gamma())

print(fit1.fit().summary()) #パラメータ推定の結果もRと同様だった

#Rの出力結果
"""
> summary(glm)
Call:
glm(formula = y ~ x, family = poisson, data = data)
Deviance Residuals:
    Min       1Q   Median       3Q      Max
-2.3679  -0.7348  -0.1775   0.6987   2.3760
Coefficients:
            Estimate Std. Error z value Pr(>|z|)
(Intercept)  1.29172    0.36369   3.552 0.000383 ***
x            0.07566    0.03560   2.125 0.033580 *
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
(Dispersion parameter for poisson family taken to be 1)

    Null deviance: 89.507  on 99  degrees of freedom
Residual deviance: 84.993  on 98  degrees of freedom
AIC: 474.77
Number of Fisher Scoring iterations: 4
"""

print("\n")
print(fit2.fit().summary())

# %%
#ポアソン回帰モデルによる予測
"""
上記のstatmodelモジュールのfourmulaによる
切片と説明変数xの最尤推定値は1.2917、0.0757
Log-Likelihood: -235.39
Deviance: 84.993 (k=1)
これより、y=exp(1.2917+0.0757x)

"""

plt.scatter(df.x[df.f=="C"],df.y[df.f=="C"],label="C", color="pink")
plt.scatter(df.x[df.f=="T"],df.y[df.f=="T"],label="T",color="green")
x=np.arange(min(df.x),max(df.x),(max(df.x)-min(df.x))/100)
y=np.exp(1.2917+0.0757*x)
plt.plot(x,y)
plt.xlabel("size")
plt.ylabel("mean of seed")
plt.legend()
plt.show()

# %%
#説明変数が数量型＋因子型
#Log-Likelihood: -235.29
#Deviance: 84.808 (k=2)
#y=exp(1.2631+0.0801x-0.0320f) Tをダミー変数

fit3=smf.glm(formula="y~x+f",data=df,family=sm.families.Poisson())
print(fit3.fit().summary())

# %%
