#%%

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# グラフとして描画するデータ
x = np.array([1,2,3,4])
y = np.array([2,3,4,5])

# グラフを描画
plt.plot(x, y)
plt.show()
# %%
x = np.array([1,2,3,4])
y = np.array([2,7,8,5])
plt.plot(x, y)
plt.show()

# %%
import csv
import pandas as pd

""""
with open("./02_distribution.csv", "r") as f:
    reader=csv.reader(f)
    header=next(reader)
    for line in reader:
        print(line)
"""
data=pd.read_csv("./02_distribution.csv")

for i in range(len(data)):
    for j in [1]:
        value=data.iloc[i,j]
        print(value)


# %%

#iloc[]で１列目全ての要素を抜き取る
#describe()はデータの基本統計量(https://obgyn.jp/describe/)を表す
summary=data.iloc[0:,1].describe() 
print(summary)

# %%
import matplotlib.pylab as plt

value=data.iloc[0:,1]
print(value)
#plt.hist(value)
# %%
import matplotlib.pyplot as plt

plt.hist(value, range=(-0.5,9.5),color="red",rwidth=0.9)

# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson


x=np.arange(0,9,1)
prob=[poisson.pmf(i,3.56) for i in x]
plt.scatter(x,prob)
plt.show()

#ヒストグラムと確率分布を組み合わせて表示
#http://pineplanter.moo.jp/non-it-salaryman/2018/03/23/python-2axis-graph/

fig,ax1= plt.subplots()
ax2=ax1.twinx()

ax1.hist(value, range=(-0.5,9.5),color="red",rwidth=0.9,label="# of seed")
ax2.plot(x,prob,linestyle="dashed",label="λ=3.56")

plt.xlabel("data")
ax1.set_ylabel("Freq")
ax2.set_ylabel('Prob')

ax1.legend(loc="upper right",bbox_to_anchor=(0.95, 1))
ax2.legend(loc="upper right",bbox_to_anchor=(0.95, 0.9))

plt.show()

# %%
#ポアソン分布の可視化
import numpy as np
from scipy.stats import poisson

#scipy.statsの特徴
#https://qiita.com/y_itoh/items/c388ff82360906240daf
#pdf(Probability density function) 確率密度関数
#pmf(Probability mass function) 確率質量関数

x=np.arange(1,50,1)
p1=[poisson.pmf(i,3.5) for i in x]
p2=[poisson.pmf(i,7.7) for i in x]
p3=[poisson.pmf(i,15.1) for i in x]

plt.plot(x,p1,linestyle="solid",color="red",label="λ=3.5")
plt.plot(x,p2,linestyle="solid",color="green",label="λ=7.7")
plt.plot(x,p3,linestyle="solid",color="blue",label="λ=15.1")

plt.xlabel("data")
plt.ylabel("prob")
plt.legend()
plt.show()

# %%
#乱数の発生と乱数の可視化
import seaborn as sns
from scipy.stats import norm
from scipy.stats import gamma

#正規分布に従う乱数
norm_rvs=norm.rvs(loc=0, scale=1, size=1000) #loc=mean,scale=scale,size=sample_size
sns.distplot(norm_rvs, bins=range(-5,5))

#ガンマ分布に従う乱数
a=5
gamma_rvs=gamma.rvs(a, scale=1, size=1000) #scale=λ,a=α　
#λ固定でαを大きくするとモードが0から中央に移る
sns.distplot(gamma_rvs)


# %%
#パラメータの最尤推定