# %%
import os
import numpy as np
import pandas as pd
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# %%
# データ生成 --------------------------------
np.random.seed(seed=1)  # 乱数を固定
N = 200  # データの数
K = 3  # 分布の数
T = np.zeros((N, 3), dtype=np.uint8)
X = np.zeros((N, 2))
X_range0 = [-3, 3]  # X0 の範囲 , 表示用
X_range1 = [-3, 3]  # X1 の範囲 , 表示用
Mu = np.array([[-.5, -.5], [.5, 1.0], [1, -.5]])  # 分布の中心
Sig = np.array([[.7, .7], [.8, .3], [.3, .8]])  # 分布の分散
Pi = np.array([0.4, 0.8, 1])  # 各分布への割合
for n in range(N):
    wk = np.random.rand()
    for k in range(K):
        if wk < Pi[k]:
            T[n, k] = 1
            break
    for k in range(2):
        X[n, k] = np.random.randn() * Sig[T[n, :] == 1, k] + \
                   Mu[T[n, :] == 1, k]


# %%
# -------- 2 分類のデータをテスト・訓練データに分割
TrainingRatio = 0.5
X_n_training = int(N * TrainingRatio)
X_train = X[:X_n_training, :]
X_test = X[X_n_training:, :]
T_train = T[:X_n_training, :]
T_test = T[X_n_training:, :]
# %%
print(X_train)
print(X_test)
print(T_train)
print(T_test)

# %%
print(X_train[:, 0].min())
print(X_train[:, 0].max())
print(X_train[:, 1].min())
print(X_train[:, 1].max())

print('---------------------')

print(X_test[:, 0].min())
print(X_test[:, 0].max())
print(X_test[:, 1].min())
print(X_test[:, 1].max())

# %%
def get_colorpalette(n_colors):
    """legendとRBG値を対応させる関数

    Returns:
        rgb (list): RBG値が格納された配列
    """
    palette = sns.hls_palette(n_colors, l=0.6, s=1)
    #palette = sns.color_palette('hls', n_colors)
    #colorpalette = 'hls'
    rgb = ['rgb({},{},{})'.format(*[x*256 for x in rgb]) for rgb in palette]
    return rgb

# %%
colors = get_colorpalette(3)

train_x0 = pd.Series(X_train[:, 0])
train_x1 = pd.Series(X_train[:, 1])
train_c = pd.Series(np.where(T_train == 1)[1])
train_cc = pd.Series([colors[int(i)] for i in train_c])

test_x0 = pd.Series(X_test[:, 0])
test_x1 = pd.Series(X_test[:, 1])
test_c = pd.Series(np.where(T_test == 1)[1])
test_cc = pd.Series([colors[int(i)] for i in test_c])

data = pd.DataFrame(
    [train_x0, train_x1, train_cc, test_x0, test_x1, test_cc]
).transpose()

data.columns = [
    'train_x0',
    'train_x1',
    'train_c',
    'test_x0',
    'test_x1',
    'test_c'
]

print(data)


# %%
colors = get_colorpalette(3)
titles = ['train_data', 'test_data']

n_row = 1
n_col = 2

fig = make_subplots(
    rows=n_row,
    cols=n_col,
    subplot_titles=titles,
    horizontal_spacing=0.1,  # will change
    #vertical_spacing=0.3     # will change
)



ii = -1
for i in range(1, n_row+1):
    for j in range(1, n_col+1):
        ii += 1

        x = data.columns[ii * 3]
        y = data.columns[(ii * 3) + 1]
        c = data.columns[(ii * 3) + 2]

        fig.add_trace(
            go.Scatter(
                x=data.loc[:, x],
                y=data.loc[:, y],
                mode='markers',
                marker=dict(
                    color=data.loc[:, c],
                    size=5,
                    line=dict(color='black', width=1)
                )
            ),
            row=i,
            col=j
        )
        fig.update_layout(
            plot_bgcolor='white',
            height=500,
            width=1000
        )
        min_x = data.loc[:, x].min() - 1
        max_x = data.loc[:, x].max() + 1
        fig.update_xaxes(
            title='x0',
            showline=True,
            linewidth=1,
            linecolor='black',
            mirror=True,
            ticks='inside',
            range=(min_x, max_x),
            dtick=0.5,
            row=i,
            col=j
        )
        min_y = data.loc[:, y].min() - 1
        max_y = data.loc[:, y].max() + 1
        fig.update_yaxes(
            title='x1',
            showline=True,
            linewidth=1,
            linecolor='black',
            mirror=True,
            ticks='inside',
            range=(min_y, max_y),
            dtick=0.5,
            row=i,
            col=j
        )

fig.for_each_xaxis(
    lambda axis: axis.title.update(
        font=dict(
            color='black',
            size=10
        )
    )
)
fig.for_each_yaxis(
    lambda axis: axis.title.update(
        font=dict(
            color='black',
            size=10
        )
    )
)
fig.update_layout(
    title=dict(
        text='classfication',
        x=0.5,
        xanchor='center'
    ),
    showlegend=False
)

fig.update_annotations(
    font=dict(size=10)
)

fig.show()

# %%
"""
MLPClassfierを用いた入力層から中間層のパラメータwの推定から
入力特徴量のデータに及ぼす影響を可視化
"""

from sklearn.datasets import load_breast_cancer
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

# %%
cancer = load_breast_cancer()

X_train, X_test, y_train, y_test = train_test_split(
    cancer.data,
    cancer.target,
    random_state=0
)

# %%
# 各データのスケールの大きさを確認
print(cancer.data.max(axis=0))
# -> オーダがまちまちなので、
# 値のとる範囲を統一したい

# %%
# 訓練データとテストデータの標準化
scaled = StandardScaler()
scaled.fit(X_train)
X_train_scaled = scaled.transform(X_train)

scaled = StandardScaler()
scaled.fit(X_test)
X_test_scaled = scaled.transform(X_test)

# %%
mlp = MLPClassifier(
    max_iter=1000,
    random_state=0
)
mlp.fit(X_train_scaled, y_train)
print("訓練セットの精度: {:.3f}".format(mlp.score(X_train_scaled, y_train)))
print("テストセットの精度: {:.3f}".format(mlp.score(X_test_scaled, y_test)))

# %%
mlp = MLPClassifier(
    max_iter=1000,
    alpha=1,
    #hidden_layer_sizes=[10, 10],
    #activation='tanh'
    random_state=0
)
mlp.fit(X_train_scaled, y_train)
print("訓練セットの精度: {:.3f}".format(mlp.score(X_train_scaled, y_train)))
print("テストセットの精度: {:.3f}".format(mlp.score(X_test_scaled, y_test)))

# %%
#print(mlp.coefs_[0])
#print(cancer.feature_names)
print(mlp.coefs_[0].shape)
# %%
import plotly.graph_objects as go
from plotly.subplots import make_subplots

y, x = mlp.coefs_[0].shape

fig = go.Figure(
    data=go.Heatmap(
        colorbar=dict(
            title='w_value'
        ),
        z=mlp.coefs_[0],
        x=list(range(0, x+1)),
        y=list(cancer.feature_names),
        colorscale='Viridis'
    )
)

fig.update_xaxes(title='w')
fig.update_yaxes(title='features')

fig.update_layout(
    title='paramete w',
    xaxis_nticks=36,
    width=1000,
    height=700,
)
fig.show()

# %%
# FNN（Forward Neural Network）の実装 ----------------------
def sigmoid(x):
    return 1 / (1 + np.exp(-x))

def FNN(wv, M, K, x):
    N, D = x.shape
    w = wv[: M * (D + 1)]
    w = w.reshape((M, D + 1))
    v = wv[K * (M + 1) :]
    v = v.reshape((K, M + 1))
    b = np.zeros((N, M + 1))
    z = np.zeros((N, M + 1))
    a = np.zeros((N, K))
    y = np.zeros((N, K))

    for i in range(N):
        for j in range(M):
            b[i, j] = np.dot(w[j, :] * np.r_[x[n, :], 1])
            z[i, j] = sigmoid(b[i, j])
        z[i, M] = 1
        wkz = 0
        for k in range(K):
            a[i, k] = np.dot(v[k, :], z[i, :])
            wkz += np.exp(a[i, k])
        for k in range(K):
            y[i, k] = np.exp(a[i, k]) / wkz

    return y, a, z, b

# %%
# 平均交差エントロピー誤差　----------------------------------
# yとtを一次元ベクトルに変換して、
# 内積を用いることでfor文を使わなくて済む

def CE_FNN(wv, M, K, x, t):
    y, _, _, _ = FNN(wv, M, K, x)
    ce = np.dot(t.reshape(-1) * np.exp(y.reshape(-1))) / ((-1) * N)

    return ce

# %%
# 平均交差エントロピー誤差微分　-------------------------------
# (1) 数値微分法
# (2) 勾配法(誤差逆伝搬法)
# dE(n)/dVkj = (yk-tk) * zj
# dE(n)/dWji = (h'(bj) * Σ(yk-tk) *Vkj) * xi

def dCEE_FNN_1(wv, M, K, x, t):
    eps = 0.0001
    dwv = np.zeros_like(wv) # 既知行列の全要素を0に初期化
    for i in range(len(wv)):
        # 指定した位置iについて、逐次近似解(2点の傾き)を出す
        wv_mod = wv.copy() # 位置iの値のみをeps分前後にずらす
        wv_mod[i] = wv[i] + eps
        ce1 = CE_FNN(wv_mod, M, K, x, t)
        wv_mod[i] = wv[i] - eps
        ce2 = CE_FNN(wv_mod, M, K, x, t)
        dwv[i] = (ce1 - ce2) / (2 * eps)

        return dwv

def dCEE_FNN_2(wv, M, K, x, t):
    N, D = x.shape
    w = wv[M * (D + 1):]
    w = w.reshape((M, D + 1))
    v = wv[: K * (M + 1)]
    v = v.reshape((K, M + 1))

    y, a, z, b = FNN(wv, M, K, x)

    dwv = np.zeros_like(wv)
    dw = np.zeros((M, D + 1))
    dv = np.zeros((K, M + 1))

    delta1 = np.zeros(M)
    delta2 = np.zeros(K)


    for i in range(N):
        for k in range(K):
            delta2[k] = y[i, k] - t[i, k]
        for j in range(M):
            delta1[j] = b[i, j] * (1 - b[i, j]) * np.dot(delta2, v[:, j])
        for k in range(K):
            dv += delta2[k] * z[i, :] / N
        for j in range(M):
            dw += delta1[k] * np.r_[x[i, :], 1] / N
    dwv = np.c_[
        dw.reshape(1, M * (D + 1)),
        dv.reshape(1, K * (M + 1))
    ]
    dwv = dwv.reshape(-1)
    return dwv

# %%

