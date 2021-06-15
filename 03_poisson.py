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
