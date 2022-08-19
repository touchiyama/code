# %%
import numpy as np
#from linear_regression import least_squares
from pca import pca

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
    SIMPLS アルゴリズムを用いて回帰係数を計算
    sklearnに実装されているアルゴリズム

    パラメータ
    ----------
    X: 入力データ
    Y: 出力データ
    R: 潜在変数の数

    戻り値
    -------
    beta: 回帰係数
    W, P, Q: PLS 1 モデルのパラメータ
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
import matplotlib.pyplot as plt
import numpy as np
import scipy.io

# %%
CNGATEST_dict = scipy.io.loadmat('/Users/tomoyauchiyama/SmallData/src/CNGATEST.mat')
# %%
print(CNGATEST_dict)
# %%
# データは辞書型として保存されているので，キーを確認します
print(CNGATEST_dict.keys())

# %%
# 辞書型から配列を取り出します
cn_sd_hl = CNGATEST_dict['cn_sd_hl']
cn_sd_ll_a = CNGATEST_dict['cn_sd_ll_a']
cn_sd_ll_b = CNGATEST_dict['cn_sd_ll_b']

cn_y_hl = CNGATEST_dict['cn_y_hl']
cn_y_ll_a = CNGATEST_dict['cn_y_ll_a']
cn_y_ll_b = CNGATEST_dict['cn_y_ll_b']

# 学習用データと検証用データ
# スペクトルデータ
Xtrain = np.vstack([cn_sd_ll_a, cn_sd_hl])
Xval = cn_sd_ll_b
# スペクトルに対応するセタン価
ytrain = np.vstack([cn_y_ll_a, cn_y_hl])
yval = cn_y_ll_b

# 配列のサイズを確認します
print(Xtrain.shape)
# %%
wave = np.arange(750, 1551, 2)
fig ,ax = plt.subplots(facecolor ="w")

plt.xlabel('Wavelength[nm]')
plt.ylabel('Intensity')
ax.plot(wave, Xtrain[26, :], label ='26th,CN=53.9')
ax.plot(wave, Xtrain[32, :], label ='32th,CN=42.8')
ax.legend()
plt.show()
# %%
wave = np.arange(750, 1551, 2)
fig ,ax = plt.subplots(facecolor ="w")

plt.xlabel('Wavelength[nm]')
plt.ylabel('Intensity')
ax.plot(wave, Xval[26, :], label ='26th,CN=53.9')
ax.plot(wave, Xval[32, :], label ='32th,CN=42.8')
ax.legend()
plt.show()

# %%
wave = np.arange(750, 1551, 2)
fig ,ax = plt.subplots(facecolor ="w")

plt.xlabel('Wavelength[nm]')
plt.ylabel('Intensity')
ax.plot(wave, Xval[1, :], label ='26th,CN=53.9')
ax.plot(wave, Xval[100, :], label ='32th,CN=42.8')
ax.legend()
plt.show()

# %%
X = np.array([[3, 2], [1, 2]])
print(X)
print(X.T)

print('# --')
XX = X.T @ X
print(XX)

XXX = np.linalg.inv(X.T @ X)
print(XXX)

# %%
np.linalg.norm(X.T @ X)

# %%
print(Xval)

# %%
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

# %%
X = np.array(
    [
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        [1, 2, 3, 98, 5, 6, 7, 8, 9, 10],
        [1, 2, 3, 4, 92, 6, 7, 8, 9, 10],
        [1, 2, 3, 4, 5, 6, 19, 8, 9, 10]
    ]
) # np.array型にする

Y = np.array(
    [
        [1, 3, 3, 4, 5, 6, 7, 89, 9, 11],
        [1, 9, 3, 98, 5, 6, 7, 92, 9, 16],
        [1, 6, 3, 4, 92, 6, 7, 65, 9, 17],
        [1, 7, 3, 4, 5, 6, 19, 91, 9, 19]
    ]
) # np.array型にする

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
def pls_cv(X, Y, maxLV, K=10):
    """
    クロスバリデーションでPLSの最適な潜在変数の数を探索
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
            X_train, meanX, stdX = autoscale(X_train)
            Y_train, meanY, stdY = autoscale(Y_train)
            X_val = scaling(X_val, meanX, stdX)

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
    # PRESS が最小となったとき潜在変数を探索します
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