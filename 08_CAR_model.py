# %%

#
# instrinstic CAR modelの実装
#

# 解析データの用意
import pyper
import pandas as pd

r = pyper.R(use_pandas='True')
r('load(\"/Users/tomoyauchiyama/statisticModel/kubobook_2012/spatial/Y.RData\")')
df = pd.Series(r.get('Y')).astype(int)

# %%

# 隣接情報と重み情報を作成
import numpy as np

N = len(df)

adj = []
for i in range(N):
    if i == 0:
        adj.append([i+1])
    elif i == N-1:
        adj.append([i])
    else:
        adj.append([i-1, i+1])

# リスト型からnumpy.array型に変換
adj = np.array(adj)
print(adj)

# %%

weight = []
for i in range(len(adj)):
    w = [1] * len(adj[i])
    weight.append(w)

# リスト型からnumpy.array型に変換
weight = np.array(weight)
print(weight)

# %%

# 隣接情報と重み情報を行列化
amat2 = np.zeros((N, N))
wmat2 = np.zeros((N, N))
for i, a in enumerate(adj):
    amat2[i, a] = 1
    wmat2[i, a] = weight[i]

# %%

print(amat2)
print(wmat2)

# %%

