# %%
from matplotlib import colors
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
