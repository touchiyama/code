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

plt.plot(x,p1,label="q=0.1",marker="o",linestyle="--")
plt.plot(x,p2,label="q=0.3",marker="x",linestyle="--")
plt.plot(x,p3,label="q=0.8",marker="^",linestyle="--")
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
#logy=z
#線形予測子z=B1+B2*x+logAを考える

df_setoff=pd.read_csv("/Users/tomoyauchiyama/statisticModel/kubobook_2012/binomial/data4b.csv")
#print(df_setoff.head())
#応答変数yについて考えると、正の離散値で域値の範囲はない　→　ポアソン分布に従うと仮定

plt.scatter(df_setoff.A,df_setoff.y, label="plant")

fit_offset=smf.glm("y~x",offset=np.log(df_setoff.A), data=df_setoff, family=sm.families.Poisson()).fit()
print(fit_offset.summary())

A=np.arange(df_setoff.A.min(),df_setoff.A.max(),(df_setoff.A.max()-df_setoff.A.min())/25)
#x=np.arange(df_setoff.x.min(),df_setoff.x.max(),(df_setoff.x.max()-df_setoff.x.min())/100)

def offset_model(x,A):
    z=fit_offset.params[0]+fit_offset.params[1]*x+np.log(A)
    return(np.exp(z))

plt.plot(A,offset_model(0.1,A),label="x=0.1",marker="*", color="gray")
plt.plot(A,offset_model(0.3,A),label="x=0.3",marker="x", color="gray")
plt.plot(A,offset_model(0.5,A),label="x=0.5",marker="^", color="gray")
plt.plot(A,offset_model(0.7,A),label="x=0.7",marker="D", color="gray")
plt.plot(A,offset_model(0.9,A),label="x=0.9",marker="+", color="gray")
plt.legend(loc="upper left")
plt.xlabel("size of A")
plt.ylabel("# of plant")
plt.show()

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

#花の質量yについてガンマ分布でモデルを作る
#花の質量y（応答変数）、葉の質量x（説明変数）とすると、
#花の質量は、正の値をとり、連続値をとるのでガンマ分布に従うと仮定できる
#ここで、y=A*x**bに従うと仮定すると
#この時、A＝exp(a)とおくと、y=exp(a)*x**b →　y=exp(a+blogx)と表せる
#これより、logy=a+blogx

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import statsmodels.api as sm
import statsmodels.formula.api as smf

df=pd.read_csv("/Users/tomoyauchiyama/statisticModel/kubobook_2012/gamma/data06.csv")
plt.scatter(df.x,df.y,label="weight of flower")

df["logx"]=np.log(df.x)
#print(df.head())
fit_gamma=smf.glm("y~logx",data=df,family=sm.families.Gamma(link=sm.families.links.log)).fit()
print(fit_gamma.summary())

"""
#predictメソッドを使った方法
df["predict"]=fit_gamma.predict()
plt.plot(df.x,df.predict)
"""

#関数を使った方法
x=np.arange(df.x.min(),df.x.max(),(df.x.max()-df.x.min())/100)
def gamma_model(x):
    return(np.exp(fit_gamma.params[0]+fit_gamma.params[1]*(np.log(x))))

def gamma_model2(x):
    return(np.exp(-1.273+0.549*(np.log(x))))

def gamma_model3(x):
    return(np.exp(-0.808+0.817*(np.log(x))))

plt.plot(x,gamma_model(x),color="black",label="coef")
plt.plot(x,gamma_model2(x),color="gray",linestyle="--",label="2.5%")
plt.plot(x,gamma_model3(x),color="gray",linestyle="--",label="97.5%")
plt.xlabel("weight of leaf")
plt.ylabel("weight of flower")
plt.legend(loc="upper left")
plt.show()

# %%
for i in fit_gamma.params:
    print(i)

# %%

#ガンマ分布を用いた予測の評価
#statsmodelsのfitメソッドから得られた予測式を用いて50%(25-75%)と90%(5-95%)の予測区間を図示する

#説明変数xとそれに対応する予測値のデータフレームを作成
df["predict"]=fit_gamma.predict()
#print(df.predict)
dt_test=pd.DataFrame({"x": df.x, "p": df.predict})

#上記の予測式の基となるモデルのガンマ分布のパラメータ(shape(分布の形)、loc(標準化)、scale(=1/r(係数)))
#scipy.statsのfitメソッドで求める
shape, loc, scale=st.gamma.fit(df.y)

#scipy.statsのintervalメソッドで、shape、loc、scaleに基づいて信頼区間75%と95％に当たる確率変数を取得(任意)
lower, upper = st.gamma.interval(alpha=0.75, a=shape, loc=loc, scale=scale)
print('信頼区間95％の下限：', lower)
print('信頼区間95％の上限：', upper)

#scipy.statsのppf (Percent point function) パーセント点関数で
#ガンマ分布 => 平均=shape/rate、scale=1/rate、分散=shape/rate**2
#ppfで必須のパラメーター：shapeとscare
#ここで、平均=exp(a+b*logx)と表されるので、rate=shape/平均=shape/exp(a+b*logx)と表せる
#これより、scale=1/rate=exp(a+b*logx)と表せる
#以下の方法を調べるまで、scaleの値はfitメソッドで求めた値に固定していた。
#モデルの予測式が応答変数の平均を示していることに気づけるかが重要
#説明変数に応じて、ppfの値を変化させなければ、予測式の信頼区間は求められない。

def pred_range(per):
    pred=[]
    for x in df.x:
        a=shape
        rate=a/np.exp(fit_gamma.params[0]+fit_gamma.params[1]*(np.log(x)))
        scale=1/rate
        lower, upper=st.gamma.ppf([per,1-per], a=a, scale=scale)
        pred.append((lower,upper))

    return pred

plt.scatter(df.x,df.y,label="weight of flower")

x=np.arange(df.x.min(),df.x.max(),(df.x.max()-df.x.min())/100)
def gamma_model(x):
    return(np.exp(fit_gamma.params[0]+fit_gamma.params[1]*(np.log(x))))

def gamma_model2(x):
    return(np.exp(-1.273+0.549*(np.log(x))))

def gamma_model3(x):
    return(np.exp(-0.808+0.817*(np.log(x))))

plt.plot(x,gamma_model(x),color="black",label="coef")
plt.plot(x,gamma_model2(x),color="gray",linestyle="--",label="2.5%")
plt.plot(x,gamma_model3(x),color="gray",linestyle="--",label="97.5%")

y1, y2 = map(list, zip(*pred_range(0.25)))
plt.fill_between(x=df.x, y1=y1, y2=y2, alpha=0.2, color="blue")
y3, y4 = map(list, zip(*pred_range(0.05)))
plt.fill_between(x=df.x, y1=y3, y2=y4, alpha=0.1, color="green")

plt.xlabel("weight of leaf")
plt.ylabel("weight of flower")
plt.legend(loc="upper left")
plt.show()


# %%

def gamma_model(x):
    return(np.exp(fit_gamma.params[0]+fit_gamma.params[1]*(np.log(x))))

xx=np.linspace(st.gamma.ppf(q=0.05, a=shape, loc=loc, scale=scale),st.gamma.ppf(q=0.95, a=shape, loc=loc, scale=scale),100)
yy=[gamma_model(st.gamma.pdf(i, a=shape, loc=loc, scale=scale)) for i in xx]

df50=pd.DataFrame({"xx":xx, "yy":yy})

print(df50)

#plt.plot(xx,yy)


#print(st.gamma.ppf(q=0.25, a=shape, loc=loc, scale=scale))

#print(x)
#print(dt_test.head())


# %%

df=pd.read_csv("/Users/tomoyauchiyama/statisticModel/kubobook_2012/binomial/data4a.csv", header=None)
hash={}
for i in range(1,len(df)):
    if(df.iloc[i,3]=="T"):
        hash[df.iloc[i,2]]=df.iloc[i,[1,3]].tolist()

#for i, k in hash.items():
#    print(i, k)

df_h=pd.DataFrame(hash).transpose()
print(df_h)

