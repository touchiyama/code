# %%

import re

file = '/Users/tomoyauchiyama/code/vcftools/A01_CC1_TN_10_pass_indel.log'

with open(file, 'r') as gtf:
    for line in gtf.readlines():
        line = line.rstrip('\n')
        if 'Sites' in line:
            num = re.compile(r'kept (\d+) out').search(line).group(1)
            print(f'{num}')
# %%
import os
import re
import glob

indir = '/Users/tomoyauchiyama/code/vcftools'

uniq = {}
legend = {}
snv = {}
indel = {}

for file in glob.glob(os.path.join(indir, '*log')):
    file_name = file.split('/')
    ID = re.compile(r'(.*)_pass_(.*)\.log').search(file_name[-1]).group(1)
    vcf_type = re.compile(r'(.*)_pass_(.*)\.log').search(file_name[-1]).group(2)
    name = re.compile(r'(.*)_(\d+)').search(ID).group(1)
    if legend.get(name):
        legend[name] += ',' + ID
    else:
        legend[name] = ID

    with open(file, 'r') as rf:
        for line in rf.readlines():
            line = line.rstrip('\n')
            if 'Sites' in line:
                variants = re.compile(r'kept (\d+) out').search(line).group(1)
                if vcf_type ==  'snv':
                    snv[ID] = int(variants)
                    print(f'snv: {variants}')
                else:
                    indel[ID] = int(variants)
                    print(f'indel: {variants}')

# %%
for key, value in legend.items():
    for ID in sorted(set([i for i in value.split(',')])):
        sum = snv[ID] + indel[ID]
        print(f'{ID}\t{snv[ID]}\t{indel[ID]}\t{sum}')


# %%

for key, value in sorted(snv.items(), key=lambda x:x[0]):
    print(key, value)

# %%

for key in sorted(snv.keys()):
    print(key)

# %%

for key in snv.keys():
    print(f'\t{key}', end='')
print('\n')

# %%

print('ID\tSNV\tINDEL\tSNV+INDEL')
for ID in sorted(snv.keys()):
    sum = snv[ID] + indel[ID]
    print(f'{ID}\t{snv[ID]}\t{indel[ID]}\t{sum}')

# %%

import seaborn as sns
import pandas as pd


def get_colorpalette(n):
    n_colors = n
    palette = sns.color_palette('hls', n_colors) #colorpalette = 'hls'
    rgb = ['rgb({},{},{})'.format(*[x*256 for x in rgb]) for rgb in palette]
    return rgb


colors = get_colorpalette(1)
json_name = './test'
samples = ['snv', 'indel', 'snvIndel']

for i in range(len(samples)):
    plot_data = pd.DataFrame()
    for j, (name, value) in enumerate(legend.items()):
        X = []
        Y = []
        for ID in sorted(set([i for i in value.split(',')])):
            x = re.compile(r'(.*)_(\d+)').search(ID).group(2)
            X.append(int(x))
            if samples[i] == 'snv':
                Y.append(int(snv[ID]))
                filename = json_name + '_snv'
            elif samples[i] == 'indel':
                Y.append(int(indel[ID]))
                filename = json_name + '_indel'
            else:
                Y.append((int(snv[ID]) + int(indel[ID])))
                filename = json_name + '_snvIndel'

        col = pd.Series([name, Y, X, 'solid', colors[j], 'none'])
        plot_data = plot_data.append(col, ignore_index=True)

        plot_data.columns = ['name', 'y', 'x', 'line_type', 'color', 'fill']

        outfile = filename + '.json'
        plot_data.to_json(outfile, orient='records', indent=4)
        print(f'Output json: {outfile}')

# %%
import numpy as np
import pandas as pd

df = pd.DataFrame(index=['idx'+str(i) for i in range(10)])
for i in range(2):
    df['col' + str(i)] = np.random.rand(10)

print(df)
# %%

res = df.corr()
print(res.loc['col0', 'col1'])

# %%
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

N = 20
x = np.linspace(0, 1, N)
y = np.random.random(N)

#print(x)
#print(y)

plt.scatter(x,y)

# %%
from sklearn.linear_model import LinearRegression

def lm(x, y):
    return LinearRegression().fit(x.reshape(-1, 1), y)

# %%

plt.scatter(x,y)
plt.plot(x, lm(x, y).predict(x.reshape(-1, 1)), linestyle='dashed', color='red')

# %%
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from sklearn.linear_model import LinearRegression
import numpy as np

def lm(x, y):
    return LinearRegression().fit(x.reshape(-1, 1), y)

legend = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']

fig = make_subplots(rows=3,
                    cols=3,
                    subplot_titles=legend,
                    horizontal_spacing=0.15,
                    vertical_spacing=0.15
                    )

x_title = 'X'
y_title = 'Y'

N = 20
x = np.linspace(0, 1, N)
y = np.random.random(N)

for i in range(1, 4):
    for j in range(1, 4):
        fig.add_trace(go.Scatter(x=x,
                                 y=y,
                                 mode='markers'
                                ),
                      row=j,
                      col=i)
        fig.add_trace(go.Scatter(x=x,
                                 y=lm(x, y).predict(x.reshape(-1, 1)),
                                 line_width=2.0,
                                 line=dict(dash='dot', color='gray')
                                ),
                      row=j,
                      col=i)
        fig.update_layout(plot_bgcolor='white',
                          height=900,
                          width=900)
        fig.update_xaxes(title=x_title,
                         color = 'black',
                         showline=True,
                         linewidth=1,
                         linecolor='black',
                         mirror=True,
                         ticks='inside',
                         row=j,
                         col=i)
        fig.update_yaxes(title=y_title,
                         color = 'black',
                         showline=True,
                         linewidth=1,
                         linecolor='black',
                         mirror=True,
                         ticks='inside',
                         row=j,
                         col=i)
        r2 = lm(x, y).score(x.reshape(-1, 1), y)
        text = '*' + 'R^2=' + str(round(r2, 3))
        fig.add_annotation(x=0.13, y=1.05,
                           text=text,
                           font=dict(size=8),
                           showarrow=False,
                           arrowhead=1,
                           row=j,
                           col=i)

fig.for_each_xaxis(lambda axis: axis.title.update(font=dict(size=10)))
fig.for_each_yaxis(lambda axis: axis.title.update(font=dict(size=10)))

fig.update_layout(title=dict(text='TEST',
                  x=0.5,
                  xanchor='center'),
                  showlegend=False)

fig.show()

# %%