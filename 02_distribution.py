# %%
import numpy as np
import csv
import pandas as pd
import matplotlib.pyplot as plt
#import pyper
#r=pyper.R()
#r("load(\"/Users/tomoyauchiyama/statisticModel/kubobook_2012/distribution/data.RData\")")

with open("/Users/tomoyauchiyama/statisticModel/kubobook_2012/distribution/02_distribute.csv","r") as df:
    lines=csv.reader(df)
    header=next(lines)
    for line in lines:
        print(line[1])

# %%
import sys
print(sys.version)
print(sys.path)

# %%
#基本統計量　pandasのモジュールに含まれる
df=pd.read_csv("/Users/tomoyauchiyama/statisticModel/kubobook_2012/distribution/02_distribute.csv")
data=df.iloc[0:,1]
summary=data.describe()
print(summary)

var=data.var()
print(var)

# %%
import matplotlib.pyplot as plt

plt.hist(data, range=(-0.5,9.5),color="red",rwidth=0.9)

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

ax1.hist(data, range=(-0.5,9.5),color="red",rwidth=0.9,label="# of seed")
ax2.plot(x,prob,linestyle="dashed",label="λ=3.56")

plt.xlabel("data")
ax1.set_ylabel("Freq")
ax2.set_ylabel('Prob')

ax1.legend(loc="upper right",bbox_to_anchor=(0.95, 1))
ax2.legend(loc="upper right",bbox_to_anchor=(0.95, 0.9))

plt.show()

# %%
from scipy.stats import poisson
#scipy.statsの使い方
#https://qiita.com/supersaiakujin/items/71540d1ecd60ced65add

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
