# %%

# tensorflowとkerasのインストール ---------------------------

# pip install tensorflow
# pip install keras

"""
tensorboard                        2.7.0
tensorboard-data-server            0.6.1
tensorboard-plugin-wit             1.8.1
tensorflow                         2.7.0
tensorflow-estimator               2.7.0
tensorflow-io-gcs-filesystem       0.23.1
"""

"""
keras                              2.7.0
"""

# %%
import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import keras.optimizers
from keras.models import Sequential
from keras.layers.core import Dense, Activation

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
# optimizerの設定では、tensorflowのkerasを用いる
from tensorflow import keras
from tensorflow.keras import optimizers

# %%
# 2層FNN(フィードフォワードニューラルネット)による3クラス分類modelの作成 -------------
np.random.seed(1)

model = Sequential()
# 中間層の設定
model.add(
    Dense(
        10,
        input_dim=2,
        activation='relu', # sigmoid, tanh
        kernel_initializer='uniform'
    )
)
# 出力層の設定
model.add(
    Dense(
        3,
        activation='softmax',
        kernel_initializer='uniform'
    )
)
# 学習方法の設定
sgd = optimizers.SGD(
    lr=0.5,
    momentum=0.0,
    decay=0.0,
    nesterov=False
)
# コンパイル
model.compile(
    optimizer=sgd,
    loss='categorical_crossentropy',
    metrics=['accuracy']
)
# 学習
startTime = time.time()
history = model.fit(
    X_train,
    T_train,
    epochs=1000,
    batch_size=100,
    verbose=1,
    validation_data=(X_test, T_test)
)
# %%
# モデル評価
score = model.evaluate(
    X_test,
    T_test,
    verbose=1
)
# %%
print(score)
print(startTime)

# %%
print(history.history['loss'])
# %%
# 学習曲線の可視化 -------------------------------------------
fig = make_subplots()
fig.add_trace(
    go.Scatter(
        x=[i for i in range(len(history.history['loss']))],
        y=history.history['loss'],
        mode='lines',
        name='訓練データ',
        line=dict(
            color='black'
        #    size=5,
        #    line=dict(color='black', width=1)
        )
    )
)
fig.add_trace(
    go.Scatter(
        x=[i for i in range(len(history.history['val_loss']))],
        y=history.history['val_loss'],
        mode='lines',
        name='テストデータ',
        line=dict(
            color='red'
        #    size=5,
        #    line=dict(color='black', width=1)
        )
    )
)
fig.update_layout(
    plot_bgcolor='white'
)
fig.update_xaxes(
    title='学習エポック数',
    showline=True,
    linewidth=1,
    linecolor='black',
    mirror=True,
    ticks='inside',
    range=(0, len(history.history['loss'])),
    dtick=250,
)
fig.update_yaxes(
    title='交差エントロピー誤差',
    showline=True,
    linewidth=1,
    linecolor='black',
    mirror=True,
    ticks='inside',
    range=(0, max(history.history['loss']) + 0.1),
    dtick=0.2,
)

# %%
# 精度の可視化 -------------------------------------------
fig = make_subplots()
fig.add_trace(
    go.Scatter(
        x=[i for i in range(len(history.history['accuracy']))],
        y=history.history['accuracy'],
        mode='lines',
        name='訓練データ',
        line=dict(
            color='black'
        #    size=5,
        #    line=dict(color='black', width=1)
        )
    )
)
fig.add_trace(
    go.Scatter(
        x=[i for i in range(len(history.history['val_accuracy']))],
        y=history.history['val_accuracy'],
        mode='lines',
        name='テストデータ',
        line=dict(
            color='red'
        #    size=5,
        #    line=dict(color='black', width=1)
        )
    )
)
fig.update_layout(
    plot_bgcolor='white'
)
fig.update_xaxes(
    title='学習エポック数',
    showline=True,
    linewidth=1,
    linecolor='black',
    mirror=True,
    ticks='inside',
    range=(0, len(history.history['accuracy'])),
    dtick=250,
)
fig.update_yaxes(
    title='正答率',
    showline=True,
    linewidth=1,
    linecolor='black',
    mirror=True,
    ticks='inside',
    range=(0, 1),
    dtick=0.2,
)


# %%
print(history.history.keys())

# %%
print(y[:, 1].reshape(xn, xn).T)


# %%

fig = make_subplots()


xn = 60  # 等高線表示の解像度
x0 = np.linspace(X_range0[0], X_range0[1], xn)
x1 = np.linspace(X_range1[0], X_range1[1], xn)
xx0, xx1 = np.meshgrid(x0, x1)
y = model.predict(np.array([xx0.ravel(), xx1.ravel()]).T)

K = 3
# 境界線表示 --------------------------
for ic in range(K):
    f = y[:, ic]
    f = f.reshape(xn, xn)
    f = f.T
    if ic == 0:
        fig.add_trace(
            go.Contour(
                x=x0,
                y=x1,
                z=f,
                showscale=True,
                colorscale='Viridis',
                contours_coloring='lines',
                contours=dict(
                    showlabels=True,
                    labelfont=dict(
                        size=6.0,
                        color='black'
                    )
                )
            )
        )
    else:
        fig.add_trace(
            go.Contour(
                x=x0,
                y=x1,
                z=f,
                showscale=False,
                colorscale='Viridis',
                contours_coloring='lines',
                contours=dict(
                    showlabels=True,
                    labelfont=dict(
                        size=6.0,
                        color='black'
                    )
                )
            )
        )
fig.update_layout(
    plot_bgcolor='white'
)
fig.update_xaxes(
    title='x0',
    showline=True,
    linewidth=1,
    linecolor='black',
    mirror=True,
    ticks='inside',
    range=(x0.min(), x0.max()),
    dtick=1.0,
)
fig.update_yaxes(
    title='x1',
    showline=True,
    linewidth=1,
    linecolor='black',
    mirror=True,
    ticks='inside',
    range=(x1.min(), x1.max()),
    dtick=1.0,
)
fig.update_layout(
    title=dict(
        text='classfication',
        x=0.5,
        xanchor='center'
    ),
    showlegend=False
)
fig.show()

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

        xn = 60  # 等高線表示の解像度
        x0 = np.linspace(X_range0[0], X_range0[1], xn)
        x1 = np.linspace(X_range1[0], X_range1[1], xn)
        xx0, xx1 = np.meshgrid(x0, x1)
        y = model.predict(np.array([xx0.ravel(), xx1.ravel()]).T)

        K = 3
        # 境界線表示 --------------------------
        for ic in range(K):
            f = y[:, ic]
            f = f.reshape(xn, xn)
            f = f.T
            if ic == 0:
                fig.add_trace(
                    go.Contour(
                        x=x0,
                        y=x1,
                        z=f,
                        showscale=True,
                        colorscale='Viridis',
                        contours_coloring='lines',
                        contours=dict(
                            showlabels=True,
                            labelfont=dict(
                                size=6.0,
                                color='black'
                            )
                        )
                    ),
                    row=i,
                    col=j
                )
            else:
                fig.add_trace(
                    go.Contour(
                        x=x0,
                        y=x1,
                        z=f,
                        showscale=False,
                        colorscale='Viridis',
                        contours_coloring='lines',
                        contours=dict(
                            showlabels=True,
                            labelfont=dict(
                                size=6.0,
                                color='black'
                            )
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
        min_x = x0.min()
        max_x = x0.max()
        fig.update_xaxes(
            title='x0',
            showline=True,
            linewidth=1,
            linecolor='black',
            mirror=True,
            ticks='inside',
            range=(min_x, max_x),
            dtick=1.0,
            row=i,
            col=j
        )
        min_y = x1.min()
        max_y = x1.max()
        fig.update_yaxes(
            title='x1',
            showline=True,
            linewidth=1,
            linecolor='black',
            mirror=True,
            ticks='inside',
            range=(min_y, max_y),
            dtick=1.0,
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
