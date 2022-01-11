# %%
import h5py
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# %%

"""h5ファイル構造
・HDF5ファイルは，階層的にデータを格納することができ，
　行列やテンソルデータをそれぞれの位置で名前付きで格納できる
・with構文直下で、ファイルの中身を操作できる。
・keyでデータを取得
"""

infile = '/Users/tomoyauchiyama/code/CNN/seq.h5'

with h5py.File('seq.h5', 'r') as hf:
    for key in hf.keys():
        print(key, hf[key].shape, hf[key].dtype)

"""output
target_labels (10,) |S29
test_in (500, 131072, 4) bool
test_out (500, 1024, 10) float16
train_in (5000, 131072, 4) bool
train_out (5000, 1024, 10) float16
valid_in (500, 131072, 4) bool
valid_out (500, 1024, 10) float16

・*_outに、10種類のChIP-seqの結果の1024bp毎のカバレッジ値が
　格納(5000セット)
"""

with h5py.File(infile, 'r') as h5:
    #for key in h5.keys(
    print(h5['test_out'][0, :])

"""1セット目の10種のChIP-seqデータ
[[1.417    0.3848   0.121    ... 1.691    0.       5.47    ]
 [2.098    0.4763   0.811    ... 2.877    0.       2.102   ]
 [1.057    0.       0.6177   ... 1.28     0.       1.48    ]
 ...
 [2.781    4.715    0.       ... 0.03041  1.3      1.687   ]
 [0.013245 7.195    0.       ... 1.354    1.794    0.9497  ]
 [0.2449   5.56     0.       ... 0.4526   1.958    1.149   ]]

・縦方向（行）が塩基の位置、横方向（列）が、各種ChIP-seqデータ
　を表している
"""

# %%
with h5py.File(infile, 'r') as h5:
    #for key in h5.keys(
    print(h5['test_out'][:,:,:].min())

# %%
# 3種のChIP-seqデータ(テストデータ)の可視化 ----------------------
fig = make_subplots()
K = 3
with h5py.File(infile, 'r') as h5:
    for i in range(K):
        cov = h5['test_out'][0, :, i]
        x = [i+1 for i in range(len(cov))]
        y = cov
        fig.add_trace(
            go.Scatter(
                x=x,
                y=y,
                name='ChIP_' + str(i+1),
                mode='lines',
                line_width=0.5,
                fill='tozeroy'
                #mode='markers',
                #marker=dict(
                #    size=5,
                #    line=dict(color='black', width=1)
                #)
            )
        )

    fig.update_layout(
        plot_bgcolor='white'
    )
    min_x = min(x)
    max_x = max(x)
    fig.update_xaxes(
        title='Position',
        showline=True,
        linewidth=1,
        linecolor='black',
        mirror=True,
        ticks='inside',
        range=(min_x, max_x)
        #dtick=1.0
    )
    min_y = h5['test_out'][:,:,:].min()
    max_y = h5['test_out'][:,:,:].max()
    fig.update_yaxes(
        title='Covarage',
        showline=True,
        linewidth=1,
        linecolor='black',
        mirror=True,
        ticks='inside',
        range=(min_y, max_y)
        #dtick=1.0
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
            text='test_data',
            x=0.5,
            xanchor='center'
        ),
        showlegend=True
    )

    fig.update_annotations(
        font=dict(size=10)
    )

    fig.show()

# %%

