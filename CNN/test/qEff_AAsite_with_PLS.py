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
# 回帰モデルの性能評価の指標には、平均的な予測誤差の指標である
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
        # 散布図をプロッ
        fig, ax = plt.subplots()
        plt.xlabel('Reference ')
        plt.ylabel('Prediction ')
        plt.scatter(y, y_hat)

        # プロットの範囲を取得
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()

        plt.plot([xmin,xmax], [ymin,ymax], color = "darkgreen", linestyle = "dashed")

        r_text = f'r={r:.2f}'
        rmse_text = f'rmse ={ rmse :.2f}'

        posx = (xmax - xmin)*0.1 + xmin
        posy_1 = (ymax - ymin)*0.8 + ymin
        posy_2 = (ymax - ymin)*0.75 + ymin
        ax.text(posx, posy_1, r_text)
        ax.text(posx, posy_2, rmse_text)

        plt.show()

    return rmse, r
