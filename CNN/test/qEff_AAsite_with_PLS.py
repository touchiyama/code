"""
・線形回帰モデルの問題点
    - 説明変数に相関がある場合、多重共線性が生じる
    - サンプル数が極端に少ない場合、モデル式（回帰係数）が一意に定まらない
    - データに含まれるノイズから大きな影響を受ける
↓
↓

・PLS(部分的最小二乗法)
    回帰分析の精度を高める目的だけでなく, 次元削減あるいは関連因子の抽出といった用法としても注目を集めている
    - 説明変数に相関がある場合、多重共線性が生じる
        → 無相関の説明変数を用いる
    - サンプル数が極端に少ない場合、モデル式（回帰係数）が一意に定まらない
        → 全てのデータを踏まえて生成された潜在変数を少数用いる
    - データに含まれるノイズから大きな影響を受ける
        → 重要度（回帰係数）の高い潜在変数から順番に用いる


・評価
    - 最初にPLSで最適なパラメータを探索し、そのパラメータでVIPスコアを用いた評価を行う。
      次にVIPスコアに基づき選択された入力値で、再度パラメータチューニングを行いPLSで学習
    - RMSEとr（実測値と予測値との相関係数）で評価

"""

# %%
import numpy as np
from pca import pca
import matplotlib.pyplot as plt

# %%
# 入力データXと出力データyの定義
X1 = np.array(
    [
        [-1.12, -0.51, 0.69],
        [-0.43, -1.12, 1.02],
        [0.37, 1.10, -0.98],
        [1.19, 0.53, -0.73]
    ]
)
X2 = np.array(
    [
        [-1.12, -0.51, 0.70],
        [-0.43, -1.12, 1.01],
        [0.36, 1.10, -0.98],
        [1.20, 0.53, -0.73]
    ]
)
y = np.array(
    [0.40, 1.17, -1.14, -0.42]
)

#%%
print(X1)
print('# ---')
print(y)
print('# ---')
yt = y.reshape(-1, 1)
print(yt)

# 行列の行数と列数の提示
print(X1.shape)
# -> サンプル数N x 入力変数Mの行列を考えるので、
# パラメータ重み行列Wは、M x R(潜在変数の数）となる

# %%
# NIPALSによるPLS1モデル

def nipals_pls1(X, y, R):

    # yのベクトル化(転置) ---
    y = y.reshape(-1, 1)

    # パラメータの保存 ---
    W = np.zeros((X.shape[1], R))
    T = np.zeros((X.shape[0], R))
    P = np.zeros((X.shape[1], R))
    D = np.zeros((R, 1))

    for r in range(R):
        w = (X.T @ y) / np.linalg.norm(X.T @ y)
        t = X @ w
        p = (X.T @ t) / (t.T @ t)
        d = (t.T @ y) / (t.T @ t)
        # デフレーション ---
        X = X - t.reshape(-1, 1) @ p.T
        y = y - t @ d

        W[:, r] = w.T
        T[:, r] = t.T
        P[:, r] = p.T
        D[r] = d.T

    beta = W @ np.linalg.inv(P.T @ W) @ D

    return beta, W, T, P, D

# %%
b1 = nipals_pls1(X1, y, 1)[0]
print(b1)
b2 = nipals_pls1(X2, y, 1)[0]
print(b2)

b1 = nipals_pls1(X1, y, 2)[0]
print(b1)
b2 = nipals_pls1(X2, y, 2)[0]
print(b2)

# %%
"""
def nipals_pls1_o(X, y, R):
    NIPALS アルゴリズムを用いて回帰係数を計算します (PLS1).

    パラメータ
    ----------
    X: 入力データ
    y: 出力データ
    R: 潜在変数の数

    戻り値
    -------
    beta: 回帰係数
    W, P, D: PLS 1 モデルのパラメータ
    T: 潜在変数

    # yをベクトル化します
    y = y.reshape(-1, 1)

    # パラメータを保存する変数を作成します
    W = np.zeros((X.shape [1], R))
    P = np.zeros((X.shape [1], R))
    D = np.zeros((R, 1))
    T = np.zeros((X.shape [0], R))

    # NIPALS を計算します
    for r in range(R):
        # 重みを求めます
        w = (X.T @ y)/ np.linalg.norm(X.T @ y)
        # 潜在変数を計算します
        t = X @ w
        # ローディングベクトルを求めます
        p =(X.T @ t)/(t.T @ t)
        # 回帰係数を求めます
        d =(t.T @ y)/(t.T @ t)
        # デフレーションを行います
        X = X - t.reshape(-1, 1) @ p.T
        y = y - t @ d
        # パラメータを格納します
        W[:,r] = w.T
        P[:,r] = p.T
        D[r] = d.T
        T[:,r] = t.T

    # 回帰係数を計算します
    beta = W @ np.linalg.inv(P.T @ W) @ D
    return beta, W, P, D, T

# %%
b1 = nipals_pls1_o(X1, y, 1)[0]
print(b1)
b2 = nipals_pls1_o(X2, y, 1)[0]
print(b2)

b1 = nipals_pls1_o(X1, y, 2)[0]
print(b1)
b2 = nipals_pls1_o(X2, y, 2)[0]
print(b2)
"""

# %%
# NIPALSによるPLS2モデル
def calc_parameter(X, Y, epsilon):
    u = Y[:, 0]
    while True:
        w = (X.T @ u) / np.linalg.norm(X.T @ u)
        t = X @ w
        c = (Y.T @ t) / np.linalg.norm(Y.T @ t)
        u_new = Y @ c
        if np.linalg.norm(u_new - u) < epsilon: break
        u = u_new

    return w, c, t

def nipals_pls2(X, Y, R, epsilon=0.01):

    # yのベクトル化(転置) ---
    y = y.reshape(-1, 1)

    # パラメータの保存 ---
    W = np.zeros((X.shape[1], R))
    T = np.zeros((X.shape[0], R))
    P = np.zeros((X.shape[1], R))
    Q = np.zeros((X.shape[1], R))

    for r in range(R):
        w, c, t = calc_parameter(X, Y, epsilon)
        p = (X.T @ t) / (t.T @ t)
        q = (y.T @ t) / (t.T @ t)
        # デフレーション ---
        X = X - np.outer(t, p)
        Y = Y - np.outer(t, q)

        W[:, r] = w
        T[:, r] = t
        P[:, r] = p
        Q[:, r] = q

    beta = W @ np.linalg.inv(P.T @ W) @ Q.T

    return beta, W, T, P, Q

# %%
def simpls(X, Y, R):
    """
    SIMPLSアルゴリズムを用いて回帰係数を計算
    sklearnに実装されているアルゴリズム

    パラメータ
    ----------
    X: 入力データ
    Y: 出力データ
    R: 潜在変数の数

    戻り値
    -------
    beta: 回帰係数
    W, P, Q: PLS1モデルのパラメータ
    T, U: 潜在変数
    cont: 寄与率
    """
    # yのベクトル化
    Y = Y.reshape(-1, 1)
    # パラメータを保存する変数を作成
    W = np.zeros((X.shape [1], R))
    P = np.zeros((X.shape [1], R))
    Q = np.zeros((Y.shape [1], R))
    T = np.zeros((X.shape [0], R))
    U = np.zeros((Y.shape [0], R))
    # 寄与率の算出で使う
    ssq = np.zeros([R,2])
    ssqX = np.sum(X**2)
    ssqY = np.sum(Y**2)

    for r in range(R):
        # 特異値分をします
        u, s, v = np.linalg.svd(Y.T @ X)
        # 最大特異値に対応する右特異値ベクトルを求めます
        w = v[0,:].T
        # 潜在変数t を求めます
        t = X @ w
        # ローディングベクトルを求めます
        p =(X.T @ t)/(t.T @ t)
        # 回帰係数qを求める
        q =(Y.T @ t)/(t.T @ t)
        # 潜在変数u を求めます
        u = Y @ q
        # デフレーションを行います
        X = X - np.outer(t, p)
        Y = Y - np.outer(t, q)
        ssq [r,0] = np.sum(X**2)/ ssqX
        ssq [r,1] = np.sum(Y**2)/ ssqY
        # パラメータを保存します
        W[:,r] = w.T
        P[:,r] = p.T
        Q[:,r] = q.T
        T[:,r] = t.T
        U[:,r] = u.T

    # 回帰係数を計算
    beta = W @ np.linalg.inv(P.T @ W)@ Q.T

    # 寄与率を計算
    cont = np.zeros([R,2])
    cont [0,:] = 1 - ssq [0,:]
    for r in range(1,R):
        cont [r,:] = ssq [r-1,:] - ssq[r,:]

    return beta, W, P, Q, T, U, cont


# %%
# VIP(Variable Importance for Projection) ---
# PLSモデル向けの変数選択手法
# 各入力変数が潜在変数を通じて出力の推定にどの程度寄与しているのかを定量化
# 式の意味合いとしてはPLSモデルによる出力予測の全変動の中で、m番目の変数が出力予測の変動に影響を及ぼしている割合を示す
# VIPスコアリングでは、閾値V=1以上のスコアの変数をPLSモデルの入力変数として選択
# （全変数Mでかけているのは、全ての変数についての重要度の平均が1になるように設計しているからで、
# 　閾値V=1とすることが多い）

def pls_vip(W, D, T):
    """
    PLSのVIPを計算

    パラメータ
    ----------
    W, D: PLS のモデルパラメータ(W:入力変数の重み、D:回帰係数)
    T: PLS の潜在変数行列

    戻り値
    -------
    sel_var: 選択された入力変数のインデックス
    vips: 入力変数のVIP
    """

    M, R = W.shape
    weight = np.zeros([R])
    vips = np.zeros([M])

    # 分母の計算
    ssy = np.diag(T.T @ T @ D @ D.T) #y=DTと表せる
    total_ssy = np.sum(ssy)

    # 分子の計算
    for m in range(M):
        for r in range(R):
            weight [r] = np.array([(W[m,r] / np.linalg.norm(W[:,r]))**2])

        vips [m] = np.sqrt(M*(ssy @ weight) / total_ssy)

    sel_var = np.arange(M)
    sel_var = sel_var [vips >=1]
    return sel_var, vips


# %%
"""
# K 分割交差検証 -----------------------------
def kfold_gauss_func(x, t, m, k):
    n = x.shape[0]
    mse_train = np.zeros(k)
    mse_test = np.zeros(k)
    for i in range(0, k):
        x_train = x[np.fmod(range(n), k) != i]   # (A)
        t_train = t[np.fmod(range(n), k) != i]   # (A)
        x_test = x[np.fmod(range(n), k) == i]    # (A)
        t_test = t[np.fmod(range(n), k) == i]    # (A)
        #wm = fit_gauss_func(x_train, t_train, m)
        #mse_train[i] = mse_gauss_func(x_train, t_train, wm)
        #mse_test[i] = mse_gauss_func(x_test, t_test, wm)
    return mse_train, mse_test

# %%
t = np.array(
    [
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        [1, 2, 3, 98, 5, 6, 7, 8, 9, 10],
        [1, 2, 3, 4, 92, 6, 7, 8, 9, 10],
        [1, 2, 3, 4, 5, 6, 19, 8, 9, 10]
    ]
) # np.array型にする

# %%
n = t.shape[0]
k = n

for i in range(0, k):
    t_t = t[np.fmod(range(n), k) != i, :]
    print(f'test: {t_t}')
    t_v = t[np.fmod(range(n), k) == i, :]
    print(f'val: {t_v}')
"""

# %%
def autoscale(X):
    meanX = np.mean(X, axis=0)  #行ごとの平均
    sdX = np.var(X, axis=0, ddof=1)
    scaleX = (X - meanX) / sdX

    return scaleX, meanX, sdX

# %%
def scaling(X, meanX, sdX):
    return (X - meanX) / sdX

# %%
def rescaling(scaleX, meanX, sdX):
    return np.multiply(scaleX, sdX) + meanX

# %%
def MSE(y, y_hat):
    return np.mean((y - y_hat)**2)

y = np.array([1, 2, 3, 8, 5, 6, 7, 8, 9, 10])
y_hat = np.array([1, 2, 3, 4, 5, 9, 7, 8, 9, 10])
print(np.mean((y - y_hat)**2))

# %%
def PRESS(y, y_hat):
    return np.diag((y - y_hat).T @ (y - y_hat))


# %%
# パラメータチューニング（最適なパラメータを探索）---
# K=データ数のとき、one-leave-out法となる
# PRESS(Prediction Resudual Sum of Squares)が最小となる時がPLSモデルの最適な潜在変数の数（最適パラメータ）となる。

def pls_cv(X, Y, maxLV, K=N):
    """
    クロスバリデーションでPLSの最適な潜在変数の数を探索
    データの分割の記述方法を修正し、RMSE算出も追加
    ---------- パラメータ
    X: 入力データ
    Y: 出力データ
    maxLV: 探索する潜在変数の最大値
    K: データ分割数（デフォルト10）
    ------- 戻り値
    optR: 最適な潜在変数の数
    press: PRESS
    """
    if type(X).__module__ != 'numpy': X = np.array(X)
    if type(Y).__module__ != 'numpy': Y = np.array(Y)

    n, m = X.shape
    n, l = Y.shape
    R = np.arange(1, maxLV + 1)
    # クロスバリデーションの結果を保存するKxR行列
    res_press = np.matrix(np.zeros((K, len(R))))
    # rmseの結果を保存するKxR行列
    res_mse = np.matrix(np.zeros((K, len(R))))

    # 配列をシャッフルします（重要）
    Z = np.hstack((X, Y))
    rng = np.random.default_rng()
    rng.shuffle(Z, axis = 0)
    X = Z[:,0:m]
    Y = Z[:,m:m+l]

    # 各潜在変数に対してクロスバリデーションにてPRESSを計算
    for i, r in enumerate(R):
        for k in range(K):
            # クロスバリデーション用にデータをK分割し、学習用と検証用データに分ける
            X_train = X[np.fmod(range(n), k) != i, :]
            X_val = X[np.fmod(range(n), k) == i, :]
            Y_train = Y[np.fmod(range(n), k) != i, :]
            Y_val = Y[np.fmod(range(n), k) == i, :]
            # 各データを標準化
            #X_train, meanX, stdX = autoscale(X_train) # 本解析では不要
            Y_train, meanY, stdY = autoscale(Y_train)
            #X_val = scaling(X_val, meanX, stdX)

            # 学習用データを用いて回帰係数を計算　(モデル式の挿入)
            beta, _, _, _, _, _, _ = simpls(X_train, Y_train, r)

            # 計算した回帰係数からX_valの予測値Y_hatを計算し元のスケールに戻す
            Y_hat = X_val @ beta
            J = Y_hat.shape[0]
            for j in range(J):
                Y_hat[j,:] = rescaling(Y_hat[j,:], meanY, stdY)

            # 潜在変数R番目で、k分割目のPRESSとRMSEを計算し保存
            res_press[k, i] = PRESS(Y_val, Y_hat)
            res_mse[k, i] = MSE(Y_val, Y_hat)

    # 各潜在変数について、PRESSとRMSEを求める
    press = np.sum(res_press, axis=0)
    rmse = np.sqrt(np.mean(res_mse, axis=0))
    # PRESS が最小となったとき潜在変数を探索
    optR = R[np.argmin(press)]
    return optR, press, rmse

# %%
"""
def pls_cv(X, Y, maxLV , K=10):
    クロスバリデーションでPLS の最適な潜在変数の数を探索します.

    パラメータ
    ----------
    X: 入力データ
    Y: 出力データ
    maxLV: 探索する潜在変数の最大値
    K: データ分割数（デフォルト10）

    戻り値
    -------
    optR: 最適な潜在変数の数
    press: PRESS

    n, m = X.shape
    n, l = Y.shape
    R = np.arange(1, maxLV + 1)
    all_index = [i for i in range(n)]
    # 分割されたデータのサンプル数
    validation_size = n // K
    # クロスバリデーションの結果を保存する変数
    result = np.matrix(np.zeros((K, len(R))))

    # 配列をシャッフルします
    Z = np.hstack((X, Y))
    rng = np.random.default_rng()
    rng.shuffle(Z, axis = 0)
    X = Z[:,0:m]
    Y = Z[:,m:m+l]

    # 各潜在変数に対してクロスバリデーションにてPRESSを計算
    for i, r in enumerate(R):
        for k in range(K):
            # クロスバリデーション用にデータをK 分割し
            # 学習用データと検証用データを選択します
            if k != K - 1:
                val_index = all_index [k * validation_size : (k+1)* validation_size - 1]
            else :
                val_index = all_index [k * validation_size :]
            train_index = [i for i in all_index if not i in val_index]

            X_train = X[train_index,:]
            X_val = X[val_index,:]
            Y_train = Y[train_index,:]
            Y_val = Y[val_index,:]

            # 各データを標準化します
            X_train, meanX, stdX = autoscale(X_train)
            Y_train, meanY, stdY = autoscale(Y_train)
            X_val = scaling(X_val, meanX, stdX)

            # 学習用データを用いて回帰係数を計算します
            beta, _, _, _, _, _, _ = simpls(X_train, Y_train, r)

            # 計算した回帰係数からX_val の予測値Y_hat を計算し
            # 元のスケールに戻します
            Y_hat = X_val @ beta
            J = Y_hat.shape[0]
            for j in range(J):
                Y_hat[j,:] = rescaling(Y_hat[j,:], meanY, stdY)

            # PRESS を計算し保存します
            press_val = PRESS(Y_val, Y_hat)
            result [k, i] = press_val
            press = np.sum(result, axis = 0)

    # PRESS が最小となったとき潜在変数を探索します
    optR = R[np.argmin(press)]
    return optR, press
"""

# %%
# 回帰モデルの評価 ---
# 回帰モデルの性能評価の指標には、平均的な予測誤差の指標である、
# 根平均二乗誤差RMSE(Root Mean Squered Error)が用いられることが多い

def pred_eval(y, y_hat, plotflg = False):
    """
    RMSEと相関係数を計算

    パラメータ
    ----------
    y: 真値
    y_hat: 予測値
    pltflg: プロットのOn/ Off(デフォルト: False)

    戻り値
    -------
    rmse: RMSE
    r: 相関係数
    """
    rmse = np.linalg.norm(y - y_hat)/np.sqrt(len(y))
    r = np.corrcoef(y, y_hat)[0,1]

    if plotflg:
        # 散布図をプロット
        fig, ax = plt.subplots()
        plt.xlabel('Reference')
        plt.ylabel('Prediction')
        plt.scatter(y, y_hat, color='black')

        # プロットの範囲を取得
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()

        #xymax = max(xmax, ymax)
        #xymin = min(xmin, ymin)
        #ax.set_xlim([xymin, xymax])
        #ax.set_ylim([xymin, xymax])

        plt.plot([xmin,xmax], [ymin,ymax], color = 'gray', linestyle = 'dashed')

        r_text = f'r={r:.2f}'
        rmse_text = f'rmse={ rmse :.2f}'

        posx = (xmax - xmin)*0.01 + xmin
        posy_1 = (ymax - ymin)*0.95 + ymin
        posy_2 = (ymax - ymin)*0.90 + ymin
        ax.text(posx, posy_1, r_text)
        ax.text(posx, posy_2, rmse_text)

        plt.show()

    return rmse, r

# %%
# 数理モデルの構築 ---
import pandas as pd
import scipy.stats as st

# %%
# データの前処理 ---
file = '/Users/tomoyauchiyama/code/CNN/test/PG4796_Ascan_log2_norm.xlsx'
df_norm = pd.read_excel(file, sheet_name='log2_norm')

# %%
X = []
for seq in df_norm['AA_seq']:
    if ('-' in seq) or ('F' in seq) or ('Y' in seq):
        x_list = '-'
    else:
        x = seq.replace('x', '0').replace('A', '1')
        x_list = [int(i) for i in x]
    X.append(x_list)

df_tmp = pd.DataFrame([X]).T
df_tmp.columns = ['X']
df_X = pd.concat([df_norm['Ascan'], df_tmp], axis=1)

# %%
norm_list_1 = []
for name in df_norm.columns:
    if '101_Norm' in name:
        norm_list_1.append(name)

norm_list_2 = []
for name in df_norm.columns:
    if '102_Norm' in name:
        norm_list_2.append(name)

# %%
df_X1 = pd.concat([df_X, df_norm[norm_list_1]], axis=1)
df_X2 = pd.concat([df_X, df_norm[norm_list_2]], axis=1)
df_X1_drop = df_X1[df_X1['X'] != '-'].reset_index(drop=True)
df_X2_drop = df_X2[df_X2['X'] != '-'].reset_index(drop=True)


# %%
"""Lasso回帰モデルの実装
・positionごとの配列の寄与を調べるために、解釈しやすいLasso回帰モデルがよく用いられている。
    - L1正規化により、完全に0となる回帰係数があるため、使われる特徴量が明確
・ハイパーパラメーター(alpha)の値が大きくなると、入力変数が0に近づく。
・alpha = 10^-1 〜 10^-10の範囲で、パラメータチューニングを行う。
・非ゼロの特徴量（回帰係数）の数が最小となるか？ - RMSEとr（実測値と予測値との相関係数）で評価
"""
"""
from sklearn.model_selection import train_test_split
from sklearn.linear_model import Lasso
from sklearn.model_selection import GridSearchCV
import sklearn.metrics as metrics

# %%
# 訓練データとテストデータに分ける ---
# 訓練データを用いて訓練（学習）して予測モデルを構築し、
# その予測モデル構築に用いなかったテストデータをどのくらい正しく予測できるかで性能評価
# パラメータのoverfittingを防ぎたい

X = np.array(df_X1_drop['X'].values.tolist())
y = st.zscore(df_X1_drop['DNA#Lung#101_Norm'].values, ddof=1, axis=0)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=0)


# %%
# グリッドリサーチ（パラメータチューニング）の実行 ---
# scoreはどのように算出されているのか？

N = X_train.shape[0]
grid_param = {'alpha': [0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0, 10]}
grid_search = GridSearchCV(Lasso(), grid_param, cv=10, n_jobs=-1)
grid_search.fit(X_train, y_train)

print('Accuracy Score (train) : ', grid_search.score(X_train, y_train))
print('Accuracy Score (test) : ', grid_search.score(X_test, y_test))

print('Best Parameters : ', grid_search.best_params_)
print('Best CV Score : ', grid_search.best_score_) # trainデータに対するCVの平均CV精度
print('Best Estimator : ', grid_search.best_estimator_)

lasso_best = grid_search.best_estimator_

# %%
# グリッドリサーチに基づく再学習 ---

lasso = lasso_best.fit(X_train, y_train)
y_test_lasso = X_test @ lasso.coef_
print('Accuracy Score (train) : ', lasso.score(X_train, y_train))
print('Accuracy Score (test) : ', lasso.score(X_test, y_test))
print('# of feat used', np.sum(lasso.coef_ != 0))

# %%
# 回帰モデルの評価 ---

#y_test_lasso = grid_search.predict(X_test) # グリッドリサーチに基づく再学習
#print('# of feat used', np.sum(grid_search.coef_ != 0))

print('R2=', metrics.r2_score(y_test, y_test_lasso))
print('RMSE=', np.sqrt(metrics.mean_squared_error(y_test, y_test_lasso)))
#print('MAE=', metrics.mean_absolute_error(y_test, y_test_lasso))

pred_eval(y_test, y_test_lasso, plotflg = True)


# %%
# positionごとの配列の寄与度（回帰係数の大きさ）を可視化 ---
fig, ax = plt.subplots()

plt.bar([1, 2, 3, 4, 5, 6, 7], lasso.coef_, width=0.02, color='black', align='center')
plt.scatter([1, 2, 3, 4, 5, 6, 7], lasso.coef_, color='black')

plt.axhline(y=0, color='gray', linestyle='--')
plt.xlabel('Position')
plt.ylabel('Contributions(Coefficient magntude)')
plt.show()
"""

# %%
"""
# Ridge回帰モデルの実装
from sklearn.model_selection import train_test_split
from sklearn.linear_model import Ridge
from sklearn.model_selection import GridSearchCV
import sklearn.metrics as metrics
import scipy.stats as st

# %%
# 訓練データとテストデータに分ける ---
# 訓練データを用いて訓練（学習）して予測モデルを構築し、
# その予測モデル構築に用いなかったテストデータをどのくらい正しく予測できるかで性能評価
# パラメータのoverfittingを防ぎたい

X = np.array(df_X1_drop['X'].values.tolist())
y = st.zscore(df_X1_drop['DNA#Cortex#101_Norm'].values, ddof=1, axis=0)
#X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=0)

# %%
# グリッドリサーチ（パラメータチューニング）の実行 ---
# scoreはどのように算出されているのか？

N = X_train.shape[0]
grid_param = {'alpha': [0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0, 10, 100, 1000]}
grid_search = GridSearchCV(Ridge(), grid_param, cv=10, n_jobs=-1)
grid_search.fit(X_train, y_train)

print('Accuracy Score (train) : ', grid_search.score(X_train, y_train))
print('Accuracy Score (test) : ', grid_search.score(X_test, y_test))

print('Best Parameters : ', grid_search.best_params_)
print('Best CV Score : ', grid_search.best_score_) # trainデータに対するCVの平均CV精度
print('Best Estimator : ', grid_search.best_estimator_)

ridge_best = grid_search.best_estimator_

# %%
# グリッドリサーチに基づく再学習 ---
ridge = ridge_best.fit(X_train, y_train)
y_test_ridge = X_test @ ridge.coef_
print('Accuracy Score (train) : ', ridge.score(X_train, y_train))
print('Accuracy Score (test) : ', ridge.score(X_test, y_test))
print('# of feat used', np.sum(ridge.coef_ != 0))

# %%
# 回帰モデルの評価 ---

print('R2=', metrics.r2_score(y_test, y_test_ridge ))
print('RMSE=', np.sqrt(metrics.mean_squared_error(y_test, y_test_ridge)))
#print('MAE=', metrics.mean_absolute_error(y_test, y_test_lasso))

pred_eval(y_test, y_test_ridge, plotflg = True)


# %%
# positionごとの配列の寄与度（回帰係数の大きさ）を可視化 ---
fig, ax = plt.subplots()

plt.bar([1, 2, 3, 4, 5, 6, 7], ridge.coef_, width=0.02, color='black', align='center')
plt.scatter([1, 2, 3, 4, 5, 6, 7], ridge.coef_, color='black')

plt.axhline(y=0, color='gray', linestyle='--')
plt.xlabel('Position')
plt.ylabel('Contributions(Coefficient magntude)')
plt.show()
"""

# %%
"""scikit-learnによるPLS回帰モデルの実装
・scikit-learnでは、特異値分解に基づくSIMPLSアルゴリズムが実装されている
・VIPスコアで変数選択する前に、PLS回帰モデルで学習
"""
from sklearn.model_selection import train_test_split
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import GridSearchCV
import sklearn.metrics as metrics

# %%
# 訓練データとテストデータに分ける ---
# 訓練データを用いて訓練（学習）して予測モデルを構築し、
# その予測モデル構築に用いなかったテストデータをどのくらい正しく予測できるかで性能評価
# パラメータのoverfittingを防ぎたい

X = np.array(df_X1_drop['X'].values.tolist())
y = st.zscore(df_X1_drop['DNA#Cortex#101_Norm'].values, ddof=1, axis=0)
#X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=0)

# %%
# グリッドリサーチ（パラメータチューニング）の実行 ---
N = X_train.shape[0]
grid_param = {
    'n_components': [2, 3, 4, 5, 6, 7],
    'scale':[True, False]
}
grid_search = GridSearchCV(PLSRegression(), grid_param, cv=10, n_jobs=-1)
grid_search.fit(X_train, y_train)

print('Accuracy Score (train) : ', grid_search.score(X_train, y_train))
print('Accuracy Score (test) : ', grid_search.score(X_test, y_test))

print('Best Parameters : ', grid_search.best_params_)
print('Best CV Score : ', grid_search.best_score_) # trainデータに対するCVの平均CV精度
print('Best Estimator : ', grid_search.best_estimator_)

pls_best = grid_search.best_estimator_

# %%
# 回帰モデルの評価 ---
pls = pls_best.fit(X_train, y_train)
pls_coef = np.hstack(pls.coef_)
y_test_pls = X_test @ pls_coef

print('R2=', metrics.r2_score(y_test, y_test_pls))
print('RMSE=', np.sqrt(metrics.mean_squared_error(y_test, y_test_pls)))
pred_eval(y_test, y_test_pls, plotflg = True)

# %%
# positionごとの配列の寄与度（回帰係数の大きさ）を可視化 ---
fig, ax = plt.subplots()

plt.bar([1, 2, 3, 4, 5, 6, 7], pls_coef, width=0.02, color='black', align='center')
plt.scatter([1, 2, 3, 4, 5, 6, 7], pls_coef, color='black')

plt.axhline(y=0, color='gray', linestyle='--')
plt.xlabel('Position')
plt.ylabel('Contributions(Coefficient magntude)')
plt.show()


# %%
"""
"""
"""scikit-learnのPLSRegression()メソッドからVIPスコアを算出
パラメータ ----------
    W, D: PLSのモデルパラメータ(W:入力変数の重み、D:回帰係数)
    T: PLSの潜在変数行列

PLS = PLSRegression()とすると、
重みベクトルW -> 右特異的ベクトルに対応（PLS.y_weights_）
回帰係数D -> PLS.coef_
PLSの潜在変数行列 -> X_train @ PLS.y_weights_
"""
"""
W = pls.y_weights_
D = pls.coef_
#T = X_train @ W

#sel_var, vips = pls_vip(W, D, T) # vipスコアの算出

# %%
beta, W, P, Q, T, U, cont = simpls(X_train, y_train, R=2)
sel_var, vips = pls_vip(W, Q.T, T)

# %%
# VIPスコアに基づく回帰係数の評価を可視化 (positionごとの配列の寄与度) ---

fig, ax = plt.subplots()

plt.bar([1, 2, 3, 4, 5, 6, 7], vips, width=0.02, color='black', align='center')
plt.scatter([1, 2, 3, 4, 5, 6, 7], vips, color='black')

plt.axhline(y=1, color='gray', linestyle='--')
plt.xlabel('Position')
plt.ylabel('Contributions(Coefficient magntude)')
plt.show()

# %%
# VIPスコアに基づきデータセットの再構築 --
X_train, X_test, y_train, y_test = train_test_split(X[:, sel_var], y, test_size=0.3, random_state=0)

# %%
# VIPスコアに基づき、グリッドリサーチ（パラメータチューニング）の再実行 ---

grid_param = {
    'n_components': [2, 3, 4, 5, 6],
    'scale':[True, False]
}
grid_search = GridSearchCV(PLSRegression(), grid_param, cv=10, n_jobs=-1)
grid_search.fit(X_train, y_train)

print('Accuracy Score (train) : ', grid_search.score(X_train, y_train))
print('Accuracy Score (test) : ', grid_search.score(X_test, y_test))

print('Best Parameters : ', grid_search.best_params_)
print('Best CV Score : ', grid_search.best_score_) # trainデータに対するCVの平均CV精度
print('Best Estimator : ', grid_search.best_estimator_)

pls_best = grid_search.best_estimator_

# %%
# 回帰モデルの評価 ---

pls = pls_best.fit(X_train, y_train)
pls_coef = np.hstack(pls.coef_)
y_test_pls = X_test @ pls_coef

print('R2=', metrics.r2_score(y_test, y_test_pls))
print('RMSE=', np.sqrt(metrics.mean_squared_error(y_test, y_test_pls)))
#print('MAE=', metrics.mean_absolute_error(y_test, y_test_lasso))
pred_eval(y_test, y_test_pls, plotflg = True)

# %%
# positionごとの配列の寄与度（回帰係数の大きさ）を可視化 ---
fig, ax = plt.subplots()

plt.bar(sel_var+1, np.hstack(pls.coef_), width=0.02, color='black', align='center')
plt.scatter(sel_var+1, np.hstack(pls.coef_), color='black')

#plt.xlim([1, 2, 3, 4, 5, 6, 7])
plt.axhline(y=0, color='gray', linestyle='--')
plt.xlabel('Position')
plt.ylabel('Contributions(Coefficient magntude)')
plt.show()
"""



# %%
import re
import pandas as pd
import matplotlib.pyplot as plt

# %%
file = '/Users/tomoyauchiyama/code/CNN/test/PG4796_Ascan_Summary.xlsx'
df_dna = pd.read_excel(file, sheet_name='DNA')
# %%
col_list = []
for name in df_dna.columns:
    if 'Mean' in name:
        col_list.append(name)
# %%
df_dna_mean = df_dna.loc[:, ['Ascan', 'AA_seq'] + col_list]

# %%
for name in col_list:
    if not 'DNA#Brain_Mean' in name:
        fig, ax = plt.subplots()

        x_label = 'DNA#Brain_Mean'
        y_label = name
        x = np.log2(df_dna_mean[x_label] + 1)
        y = np.log2(df_dna_mean[y_label] + 1)
        xymax = max(np.max(x), np.max(y)) + 1 #xとyの最大値を定義
        ax.set_xlim([0, xymax])
        ax.set_ylim([0, xymax])
        ax.scatter(x, y, color='black', s=10.0)
        #ax.scatter(top100_xy[x_label], top100_xy[y_label], color='red', s=2.0, label='top100')
        ax.plot([0, xymax], [0, xymax], color='gray', linestyle='--', label='FC=1')

        x_name = re.compile(r'DNA#(.*)_Mean').search(x_label).group(1)
        y_name = re.compile(r'DNA#(.*)_Mean').search(y_label).group(1)
        #ax.axvline(5, color='gray', linestyle='--', label=f'{x_name}_G') # qPCRの値
        #ax.axhline(5, color='gray', linestyle='--', label=f'{y_name}_G') # qPCRの値

        x_name = re.compile(r'DNA#(.*)').search(x_label).group(1)
        y_name = re.compile(r'DNA#(.*)').search(y_label).group(1)
        ax.set_xlabel(f'log2(({x_name})+1)')
        ax.set_ylabel(f'log2(({y_name})+1)')
        ax.legend(
            bbox_to_anchor=(1.3, 1.0),
            loc='upper right',
            #borderaxespad=0,
            fontsize=8
        )

        plt.show()

# %%
Gbeta = {
    'Brain': 276.69,
    'Cortex': 209.1825,
    'Heart': 36.7075,
    'Kidney': 285.9275,
    'Liver': 142.7575,
    'Lung': 254.83499999999998,
    'Muscle': 113.065,
    'Pancreas': 1.1175000000000002,
    'Spinal_cord': 302.3875,
    'Spleen': 797.98
}

# %%
file = '/Users/tomoyauchiyama/code/CNN/test/PG4796_Ascan_Summary.xlsx'
df_dna = pd.read_excel(file, sheet_name='DNA')

# %%
col_list_1 = []
for name in df_dna.columns:
    if '101_Norm' in name:
        col_list_1.append(name)

col_list_2 = []
for name in df_dna.columns:
    if '102_Norm' in name:
        col_list_2.append(name)

df_dna_norm_1 = np.log2(df_dna.loc[:, col_list_1] + 1)
df_dna_norm_2 = np.log2(df_dna.loc[:, col_list_2] + 1)

# %%
df_dna_norm_1 = pd.concat([df_dna.loc[:, ['Ascan', 'AA_seq']], df_dna_norm_1], axis=1)
df_dna_norm_2 = pd.concat([df_dna.loc[:, ['Ascan', 'AA_seq']], df_dna_norm_2], axis=1)

# %%
dna_12 = [
    'Ascan', 'AA_seq', 'DNA#Brain', 'DNA#Cortex', 'DNA#Heart', 'DNA#Kidney', 'DNA#Liver',
    'DNA#Lung', 'DNA#Muscle', 'DNA#Spinal_cord', 'DNA#Spleen'
]
df_dna_norm_1.columns = dna_12
df_dna_norm_2.columns = dna_12

# %%
df_Tab_Norm = pd.DataFrame()
for i in range(len(df_dna_norm_1)):
    rep1_Vab = df_dna_norm_1.iloc[i, 2:]
    rep2_Vab = df_dna_norm_2.iloc[i, 2:]
    ave = (rep1_Vab + rep2_Vab) / 2
    total = ((((rep1_Vab - ave) **2) + ((rep2_Vab - ave) **2))) / 2
    std = np.sqrt(total.values.tolist())
    info = df_dna_norm_1.iloc[i, :2].tolist()
    col = rep1_Vab.tolist() + rep2_Vab.tolist() + ave.tolist() + std.tolist()
    df_Tab_Norm = df_Tab_Norm.append(pd.DataFrame(info + col).T, ignore_index=True)

df_Tab_Norm = df_Tab_Norm.replace([np.inf, -np.inf], 0).fillna(0)

# %%
dna_12_mean = [
    'DNA#Brain_Tab_Mean', 'DNA#Cortex_Tab_Mean', 'DNA#Heart_Tab_Mean', 'DNA#Kidney_Tab_Mean', 'DNA#Liver_Tab_Mean',
    'DNA#Lung_Tab_Mean', 'DNA#Muscle_Tab_Mean', 'DNA#Spinal_cord_Tab_Mean', 'DNA#Spleen_Tab_Mean'
]
dna_12_sd = [
    'DNA#Brain_Tab_Sd', 'DNA#Cortex_Tab_Sd', 'DNA#Heart_Tab_Sd', 'DNA#Kidney_Tab_Sd', 'DNA#Liver_Tab_Sd',
    'DNA#Lung_Tab_Sd', 'DNA#Muscle_Tab_Sd', 'DNA#Spinal_cord_Tab_Sd', 'DNA#Spleen_Tab_Sd'
]
df_Tab_Norm.columns =  ['Ascan', 'AA_seq'] + col_list_1 + col_list_2 + dna_12_mean + dna_12_sd


# %%
FC_1 = []
for i, name in enumerate(col_list_1):
    rename = name.replace('_Norm', '_FC')
    FC_1.append(rename)
    FC = df_Tab_Norm[name] / df_Tab_Norm['DNA#Brain#101_Norm']
    df_FC = pd.DataFrame(FC)
    df_FC.columns = [rename]
    df_tmp = pd.concat([df_Tab_Norm['Ascan'], df_FC], axis=1)
    if i == 0:
        df_dna_FC_1 = df_tmp
    else:
        df_dna_FC_1 = pd.merge(df_dna_FC_1, df_tmp, how='left', on='Ascan')


FC_2 = []
for i, name in enumerate(col_list_2):
    rename = name.replace('_Norm', '_FC')
    FC_2.append(rename)
    FC = df_Tab_Norm[name] / df_Tab_Norm['DNA#Brain#102_Norm']
    df_FC = pd.DataFrame(FC)
    df_FC.columns = [rename]
    df_tmp = pd.concat([df_Tab_Norm['Ascan'], df_FC], axis=1)
    if i == 0:
        df_dna_FC_2 = df_tmp
    else:
        df_dna_FC_2 = pd.merge(df_dna_FC_2, df_tmp, how='left', on='Ascan')

# %%
df_dna_FC_1 = df_dna_FC_1.replace([np.inf, -np.inf], 0).fillna(0)
df_dna_FC_2 = df_dna_FC_2.replace([np.inf, -np.inf], 0).fillna(0)


# %%
dna_12 = [
    'Ascan', 'DNA#Brain', 'DNA#Cortex', 'DNA#Heart', 'DNA#Kidney', 'DNA#Liver',
    'DNA#Lung', 'DNA#Muscle', 'DNA#Spinal_cord', 'DNA#Spleen'
]
df_dna_FC_1.columns = dna_12
df_dna_FC_2.columns = dna_12


# %%
df_Tab_FC = pd.DataFrame()
for i in range(len(df_dna_FC_1)):
    rep1_Vab = df_dna_FC_1.iloc[i, 1:]
    rep2_Vab = df_dna_FC_2.iloc[i, 1:]
    ave = (rep1_Vab + rep2_Vab) / 2
    total = ((((rep1_Vab - ave) **2) + ((rep2_Vab - ave) **2))) / 2
    std = np.sqrt(total.values.tolist())
    info = df_dna_FC_1.iloc[i, :1].tolist()
    col = rep1_Vab.tolist() + rep2_Vab.tolist() + ave.tolist() + std.tolist()
    df_Tab_FC = df_Tab_FC.append(pd.DataFrame(info + col).T, ignore_index=True)

df_Tab_FC = df_Tab_FC.replace([np.inf, -np.inf], 0).fillna(0)

# %%
dna_12_mean = [
    'DNA#Brain_Tab_Mean', 'DNA#Cortex_Tab_Mean', 'DNA#Heart_Tab_Mean', 'DNA#Kidney_Tab_Mean', 'DNA#Liver_Tab_Mean',
    'DNA#Lung_Tab_Mean', 'DNA#Muscle_Tab_Mean', 'DNA#Spinal_cord_Tab_Mean', 'DNA#Spleen_Tab_Mean'
]
dna_12_sd = [
    'DNA#Brain_Tab_Sd', 'DNA#Cortex_Tab_Sd', 'DNA#Heart_Tab_Sd', 'DNA#Kidney_Tab_Sd', 'DNA#Liver_Tab_Sd',
    'DNA#Lung_Tab_Sd', 'DNA#Muscle_Tab_Sd', 'DNA#Spinal_cord_Tab_Sd', 'DNA#Spleen_Tab_Sd'
]
df_Tab_FC.columns =  ['Ascan'] + FC_1 + FC_2 + dna_12_mean + dna_12_sd
df_Tab_FC['AA_seq'] = df_Tab_Norm['AA_seq']

# %%
df_Tab_Norm.head()

 # %%
Tab_list = []
upset_data = {}
for name in dna_12_mean:
    if not 'DNA#Brain_Tab_Mean' in name:
        x_label = 'DNA#Brain_Tab_Mean'
        y_label = name
        x_name = 'Brain'
        y_name = re.compile(r'DNA#(.*)_Tab_Mean').search(name).group(1)
        x_cutoff = np.log2(Gbeta[x_name] * 1.0)
        y_cutoff = np.log2(Gbeta[y_name] * 1.0)
        fc_cutoff = 2 * (y_cutoff / x_cutoff)

        df_sel = df_Tab_Norm[(df_Tab_Norm[y_label] >= y_cutoff)]
        df_fc = df_sel[y_label] / df_sel[x_label]
        df_merge = pd.concat([df_sel['Ascan'], df_fc], axis=1)
        df_sel_fc = df_merge[(df_merge[0] >= fc_cutoff)]

        upset_data[y_name] = df_sel_fc['Ascan'].tolist()
        for ascan in df_sel_fc['Ascan']:
            Tab_list.append(ascan)

Tab_ascan = list(set(Tab_list))

# %%
FC_idx = sorted([df_Tab_FC[df_Tab_FC['Ascan'] == name].index[0] for name in Tab_ascan])

# %%
df_Tab_FC.head()

# %%
# 組織で順位付け
re_dna_12_mean = [
    'DNA#Brain_Tab_Mean', 'DNA#Cortex_Tab_Mean', 'DNA#Spinal_cord_Tab_Mean',
    'DNA#Heart_Tab_Mean', 'DNA#Kidney_Tab_Mean', 'DNA#Liver_Tab_Mean',
    'DNA#Lung_Tab_Mean', 'DNA#Muscle_Tab_Mean', 'DNA#Spleen_Tab_Mean'
]
re_dna_12_sd = [
    'DNA#Brain_Tab_Sd', 'DNA#Cortex_Tab_Sd', 'DNA#Spinal_cord_Tab_Sd',
    'DNA#Heart_Tab_Sd', 'DNA#Kidney_Tab_Sd', 'DNA#Liver_Tab_Sd',
    'DNA#Lung_Tab_Sd', 'DNA#Muscle_Tab_Sd',  'DNA#Spleen_Tab_Sd'
]

for i in FC_idx:
    seqID = df_Tab_FC.loc[i, 'Ascan']
    AA_seq = df_Tab_FC.loc[i, 'AA_seq']
    x = [name.replace('_Tab_Mean', '') for name in re_dna_12_mean]
    #x = [i for i in range(len(re_dna_12_mean))]
    y = df_Tab_FC.loc[i, re_dna_12_mean].tolist()
    z = df_Tab_FC.loc[i, re_dna_12_sd].tolist()

    fig, ax = plt.subplots() # figsize=(8, 3)
    error_bar_set = dict(lw=1, capthick=0.5, capsize=2)
    ax.errorbar(
        x, y, yerr=z, capsize=3, elinewidth=1,
        fmt='o', ecolor='black', mfc='None', mec='black'
    )
    ax.axvspan(x[0], x[2], color='yellow', alpha=0.1)
    plt.ylim(0, )
    plt.xticks(rotation=90)
    plt.axhline(y=1, color='red', linestyle='--')
    plt.title(f'{seqID} ({AA_seq})')
    plt.ylabel('Fold-change(Brain=1)')
    #fname = 'PG4699_' + 'tissue_score_' + 'rank_' + str(i+1) + '.png'
    #plt.savefig('A:/PROJECT/PG4699/200_INFO/08_rep_perfomance/' + fname, bbox_inches='tight', dpi=300)
    plt.show()

# %%
outfile = '/Users/tomoyauchiyama/code/CNN/test/PG4796_Ascan_log2_norm.xlsx'
with pd.ExcelWriter(outfile) as writer:
    df_Tab_Norm.to_excel(writer, sheet_name='log2_norm', index=False)
    df_Tab_FC.to_excel(writer, sheet_name='FoldChange(Brain=1)', index=False)


# %%
from upsetplot import from_contents
from matplotlib import pyplot as plt
from upsetplot import plot
# https://upsetplot.readthedocs.io/en/stable/api.html

# %%
ascan = from_contents(upset_data)

fig = plt.figure(figsize=(12, 5))
plot(ascan, sort_by='cardinality', with_lines=False, show_counts=True, fig=fig, element_size=None)
plt.suptitle('Increased element_size')
plt.show()

# %%
fig, ax = plt.subplots()
# 棒の配置位置、ラベルを用意
labels = [
    '1', '2', '3', '4', '5', '6', '7', '8', '9'
]
x = np.array([i+1 for i in range(len(labels))])

# 各系列のデータを用意
data = {}
for i in range(len(labels)):
    data[str(i+1)] = np.random.rand(7)

# マージンを設定
margin = 0.2  #0 <margin< 1
totoal_width = 1 - margin

# 棒グラフをプロット
xx = np.array([i+1 for i in range(7)])

for i, h in enumerate(labels):
  pos = x - (totoal_width * (1 - ((2*i+1)/len(labels))) / 2)
  plt.bar(pos, data[h], width=totoal_width/(len(data[h])))

# ラベルの設定
plt.xticks(x, labels, rotation=90)
plt.xlabel('Position')
plt.show()

# %%
pos = x - (totoal_width * (1 - ((2*i+1)/len(labels))) / 2)


# %%
plt.bar(2, data[h], width=totoal_width/(len(data[h])))