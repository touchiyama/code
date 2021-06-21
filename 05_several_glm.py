# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df=pd.read_csv("/Users/tomoyauchiyama/statisticModel/kubobook_2012/binomial/data4a.csv")
print(df.head())

plt.scatter(df.x[df.f=="C"],df.y[df.f=="C"],label="C", color="pink")
plt.scatter(df.x[df.f=="T"],df.y[df.f=="T"],label="T",color="green")
plt.legend(loc="upper left")
plt.ylabel("# of seed")
plt.xlabel("size")
plt.show()

# %%
import scipy.stats as st

#二項分布で表現するカウントデータ
x=np.arange(0,9,1)

#N=8、q=0.1の時、y=iとなる確率
p1=[st.binom.pmf(i,8,0.1) for i in x]
p2=[st.binom.pmf(i,8,0.3) for i in x]
p3=[st.binom.pmf(i,8,0.8) for i in x]

plt.plot(x,p1,label="q=0.1",marker="o",linestyle='--')
plt.plot(x,p2,label="q=0.3",marker="x",linestyle='--')
plt.plot(x,p3,label="q=0.8",marker="^",linestyle='--')
plt.legend(loc="upper right")
plt.ylabel("p(y|8,q)")
plt.xlabel("y")
plt.show()

# %%
#二項分布を用いたGLMの1つであるロジスティック回帰について
#ロジスティック回帰の確率分布は二項分布
#変数zは線形予測子

def logistic(z):
    return 1/(1+np.exp((-1)*z))

x=np.arange(-6,6,0.1)
plt.plot(x,logistic(x))
# %%

#生存確率q=1/1+exp(-z)
#式変形：log(q/1-q)=z=b1+b2x+b3f
#パラメータ推定

import statsmodels.api as sm
import statsmodels.formula.api as smf

fit=smf.glm("y+I(N-y)~x+f",data=df,family=sm.families.Binomial()).fit()
print(fit.summary())

# %%
#肥料を与えた場合（f=T）
plt.scatter(df.x[df.f=="T"],df.y[df.f=="T"],label="T")
xx=np.arange(df.x[df.f=="T"].min(),df.x[df.f=="T"].max(),0.1)
#predictメソッドを用いてパラメータの最尤推定値(fit)を基にモデル式を予測
yy=fit.predict(exog=pd.DataFrame({"x":xx,"f":["T"]*len(xx)}))

plt.plot(xx,yy*8)

plt.show()

# %%
#pandasで、xとfについてデータフレームを構築
print(pd.DataFrame({"x":xx, "f":["T"]*len(xx)}))

# %%

plt.scatter(df.x[df.f=="T"],df.y[df.f=="T"],label="T")
def logistic(x):
    fit=smf.glm("y+I(N-y)~x+f",data=df,family=sm.families.Binomial()).fit()
    z=fit.params["Intercept"]+fit.params["x"]*x+fit.params["f[T.T]"]
    return 1/(1+np.exp((-1)*z))

x=np.arange(df.x[df.f=="T"].min(),df.x[df.f=="T"].max(),0.1)
plt.plot(x,logistic(x)*8)
plt.legend()
plt.show()
# %%
plt.scatter(df.x[df.f=="C"],df.y[df.f=="C"],label="C")
def logistic(x):
    fit=smf.glm("y+I(N-y)~x+f",data=df,family=sm.families.Binomial()).fit()
    z=fit.params["Intercept"]+fit.params["x"]*x
    return 1/(1+np.exp((-1)*z))

x=np.arange(df.x[df.f=="T"].min(),df.x[df.f=="T"].max(),0.1)
plt.plot(x,logistic(x)*8)
plt.legend()
plt.show()

# %%

def logistic(x,fit,f):
    z=fit.params["Intercept"]+fit.params["x"]*x+fit.params["f[T.T]"]*f
    return 1/(1+np.exp((-1)*z))

plt.scatter(df.x[df.f=="C"],df.y[df.f=="C"],label="C",marker="o",color="pink")
plt.scatter(df.x[df.f=="T"],df.y[df.f=="T"],label="T",marker="x",color="green")

x=np.arange(df.x[df.f=="T"].min(),df.x[df.f=="T"].max(),0.1)
fit=smf.glm("y+I(N-y)~x+f",data=df,family=sm.families.Binomial()).fit()
C=logistic(x,fit,0)*8
T=logistic(x,fit,1)*8

plt.plot(x,C,color="pink")
plt.plot(x,T,color="green")
plt.legend(loc="upper left")
plt.ylabel("# of seed")
plt.xlabel("size")
plt.show()

# %%

#二項分布を用いた5つのモデルのAIC
#1)一定モデル（定数）log(q)=β1
#2)xモデル log(q)=β1+β2x
#3)fモデル log(q)=β1+β3f
#4)x+fモデル log(q)=β1+β2x+β3f
#5)x*fモデル（相互作用項）log(q)==β1+β2x+β3f+β4xf

formula=["y+I(N-y)~1","y+I(N-y)~x","y+I(N-y)~f","y+I(N-y)~x+f","y+I(N-y)~x*f"]
for i in formula:
    fit=smf.glm(i,data=df,family=sm.families.Binomial()).fit().aic
    print(f"{i} AIC: {fit}")

"""
y+I(N-y)~1 AIC: 644.4093416623778
y+I(N-y)~x AIC: 364.34544328371595
y+I(N-y)~f AIC: 637.7597534566678
y+I(N-y)~x+f AIC: 272.21112928522336
y+I(N-y)~x*f AIC: 273.6105967259739
"""

# %%

#x*fモデル（相互作用項）log(q)==β1+β2x+β3f+β4xf
fit=smf.glm(i,data=df,family=sm.families.Binomial()).fit()
print(fit.summary())

# %%

def logistic_1(x,fit,f):
    z=fit.params["Intercept"]+fit.params["x"]*x+fit.params["f[T.T]"]*f
    return 1/(1+np.exp((-1)*z))

def logistic_2(x,fit,f):
    z=fit.params["Intercept"]+fit.params["x"]*x+fit.params["f[T.T]"]*f+fit.params["x:f[T.T]"]*x*f
    return 1/(1+np.exp((-1)*z))

plt.scatter(df.x[df.f=="C"],df.y[df.f=="C"],label="C",marker="o",color="pink")
plt.scatter(df.x[df.f=="T"],df.y[df.f=="T"],label="T",marker="x",color="green")

x=np.arange(df.x[df.f=="T"].min(),df.x[df.f=="T"].max(),0.1)

fit=smf.glm("y+I(N-y)~x+f",data=df,family=sm.families.Binomial()).fit()
C_1=logistic_1(x,fit,0)*8
T_1=logistic_1(x,fit,1)*8

fit=smf.glm("y+I(N-y)~x*f",data=df,family=sm.families.Binomial()).fit()
C_2=logistic_2(x,fit,0)*8
T_2=logistic_2(x,fit,1)*8

plt.plot(x,C_1,label="C(x+f)",color="pink")
plt.plot(x,T_1,label="T(x+f)",color="green")
plt.plot(x,C_2,label="C(x*f)",color="red",linestyle='--')
plt.plot(x,T_2,label="T(x*f)",color="blue",linestyle='--')
plt.legend(loc="lower right")
plt.ylabel("# of seed")
plt.xlabel("size")
plt.show()

#相互作用項を追加しても、x+fモデルと変わらない
#むしろ、相互作用項を追加することで、モデルの予測力が悪化した
#AIC(x+f),y+I(N-y)~x+f: 272.21112928522336
#AIC(x*f),y+I(N-y)~x*f: 273.6105967259739
#→むやみに相互作用項を入れない

# %%
#割算値いらずのセットオフ項を用いた手法
#観測値同士の割算による統計モデリングは2つの不具合を生む
#1)情報が失われる
#2)誤差が入った数値の割算値は、確率分布が複雑になる

#y=A*exp(B1+B2*x)　→　exp(B1+B2*x+logA)
#線形予測子z=B1+B2*x+logAを考える

fit_offset=smf.glm("y~x",offset=np.log(df["A"]), data=df, family=sm.families.Poisson())

# %%

#正規分布とその尤度、正規分布の一般化線形モデル（GLM）
import scipy.stats as st
x=np.arange(-5,5,0.1)

p1=[st.norm.pdf(i,0,1) for i in x]
p2=[st.norm.pdf(i,0,3) for i in x]
p3=[st.norm.pdf(i,2,1) for i in x]

plt.plot(x,p1,label="μ=0,σ=1",marker="o",linestyle='--')
plt.plot(x,p2,label="μ=0,σ=3",marker="x",linestyle='--')
plt.plot(x,p3,label="μ=2,σ=1",marker="^",linestyle='--')
plt.legend(loc="upper left")
plt.ylabel("p(y|μ,σ)")
plt.xlabel("y")
plt.show()

# %%
#累積分布関数cdfを用いて、p(1.2<=y<=1.8|0,1)を計算
p=st.norm.cdf(1.8,0,1)-st.norm.cdf(1.2,0,1)
print(p)

# %%
#ガンマ分布のGLM
x=np.arange(0,5,0.1)

#scaleをつける　
p1=[st.gamma.pdf(i,1,scale=1/(1**2)) for i in x]
p2=[st.gamma.pdf(i,5,scale=5/(5**2)) for i in x]
p3=[st.gamma.pdf(i,0.1,scale=0.1/(0.1**2)) for i in x]

plt.plot(x,p1,label="r=1,s=1",marker="o",linestyle='--')
plt.plot(x,p2,label="r=5,s=5",marker="x",linestyle='--')
plt.plot(x,p3,label="r=0.1,s=0.1",marker="^",linestyle='--')
plt.legend(loc="upper right")
plt.ylabel("p(y|t,s)")
plt.xlabel("y")
plt.show()

# %%
