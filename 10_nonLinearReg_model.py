# %%

#
# 線形ガウス基底回帰モデル構築
#

import numpy as np
import pandas as pd
import seaborn as sns
from scipy.optimize import minimize
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# %%
#--------------------------------------
# 関数の定義
#--------------------------------------

def model(w, x):
    """線形ガウス基底回帰モデル

    Args:
        w ([type]): [description]
        x ([type]): [description]

    Returns:
        [type]: [description]
    """
    #return w[0]-w[1]*np.exp(-w[2]*x)
    return w[0]*(1-np.exp(-w[1]*x))

def mse_model(w, x, y):
    """誤差2乗平均

    Args:
        w ([type]): [description]
        x ([type]): [description]
        t ([type]): [description]

    Returns:
        [type]: [description]
    """
    m = model(w, x)
    return np.mean((m-y)**2)

def fit_model(w_init, x, t):
    """パラメータ最適化

    Args:
        w_init ([type]): [description]
        x ([type]): [description]
        t ([type]): [description]

    Returns:
        [type]: [description]
    """
    res = minimize(mse_model, w_init, args=(x, t), method="powell")
    return res.x

# %%
#--------------------------------------------
# 欠損値を含まないデータを縦長のデータフレームに加工
#--------------------------------------------

#======  read row data without N.A  ======#
normal_df = pd.read_csv('/Users/tomoyauchiyama/RStanBook/chap08/input/data-conc-2.txt')

print(normal_df)

# %%
#======  convert wide DataFreme to long DataFrame  ======#
"""remove row including N.A
   substitute int type of number to add Time value
   add TimeID column linked to Time value after sorting DataFrame by 'Time'
   sort DataFrame by 'PersonID' and 'TimeID'
   reset index
"""

normal_l_df = pd.melt(normal_df, id_vars='PersonID', ignore_index=True).dropna()
normal_l_df.columns = ('PersonID', 'Time', 'Y')
normal_l_df['Time'] = normal_l_df.loc[:, 'Time'].str.extract(r'Time(\d+)').astype(int)
normal_l_df_sort = normal_l_df.sort_values('Time')
normal_l_df_sort['TimeID'] = normal_l_df_sort['Time'].factorize()[0] + 1
normal_l_df_sort = normal_l_df_sort.sort_values(['PersonID', 'Time']).reset_index()

print(len(normal_l_df_sort))
print(normal_l_df_sort)

# %%
#======  prepare for data visualization ======#
X = {}
Y = {}
for i in range(len(normal_l_df_sort)):
    ID = int(normal_l_df_sort.loc[i, 'PersonID'])
    x_val = str(normal_l_df_sort.loc[i, 'Time'])
    y_val = str(normal_l_df_sort.loc[i, 'Y'])
    if X.get(ID):
        X[ID] += ',' + x_val
        Y[ID] += ',' + y_val
    else:
        X[ID] = x_val
        Y[ID] = y_val

proc_n_df = pd.DataFrame()
for ID in sorted(X.keys()):
    xx = [int(i) for i in X[ID].split(',')]
    yy = [float(i) for i in Y[ID].split(',')]
    col = pd.Series([ID, xx, yy])
    proc_n_df = proc_n_df.append(col, ignore_index=True)

proc_n_df.columns = ['PersonID', 'Time', 'Y']

print(proc_n_df)


# %%
#--------------------------------------
# main関数
#--------------------------------------

x_title = 'Time(hour)'
y_title = 'Conc(mg/ml)'

fig = make_subplots()

max_x = normal_l_df_sort.loc[:, 'Time'].unique().max() + 3
min_x = -0.2
max_y = normal_df.describe().loc['max', :].max() + 3
min_y = 0


y = np.array(proc_n_df.loc[2, 'Y'])
x = np.array(proc_n_df.loc[2, 'Time'])

w_init=[1, 1]
w = fit_model(w_init, x, y)
mse = mse_model(w, x, y)

# %%
pre_x = np.linspace(0, max_x, 60)
pre_y = model(w, pre_x)

fig.add_trace(go.Scatter(
            x=pre_x,
            y=pre_y,
            mode='lines',
            line=dict(color=None, width=2)
            ))
fig.add_trace(go.Scatter(
            x=x,
            y=y,
            mode='markers',
            marker=dict(
                color=None,
                size=5,
                line=dict(color='black', width=1)
                )
            ))
fig.update_layout(plot_bgcolor='white'
                #height=800,
                #width=900
                )
fig.update_xaxes(title=x_title,
                showline=True,
                linewidth=1,
                linecolor='black',
                mirror=True,
                ticks='inside',
                range=(min_x, max_x),
                dtick=5)
fig.update_yaxes(title=y_title,
                showline=True,
                linewidth=1,
                linecolor='black',
                mirror=True,
                ticks='inside',
                range=(min_y, max_y),
                dtick=5)

# %%
#print('w0={0:.1f}, w1={1:.1f}, w2={2:.1f}'.format(w[0], w[1], w[2]))
print('w0={0:.1f}, w1={1:.1f}'.format(w[0], w[1]))
print('SD={0:.2f} mg/ml'.format(np.sqrt(mse)))

# %%
