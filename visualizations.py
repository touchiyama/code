import os
import re
import sys
import glob
import pandas as pd
import numpy as np
import seaborn as sns
import scipy.stats as st
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots


def groups_freq_scatter(lotID, df, bgs, treats, outdir, prefix):
    logger.info('')
    logger.info('Start making freq_dist_scatter..')

    # cheak num of samples ---
    if len(bgs) != len(treats):
        logger.error('')
        sys.exit()

    # cheak sample name ---

    # plot ---
    fig, ax = plt.subplots() # figsize=(11.8, 8.6)figsize=(9.8, 6.6)

    for i, (bg, tret) in enumerate(zip(bgs, treats), 1):
        bg_name = bg + '_Freq'
        tret_name = tret + '_Fraq'
        x = df[bg_name]
        y = df[tret_name]
        ax.scatter(x, y, alpha=1.0-(0.1*i), s=8, label=f'{tret}(n={len(x)})')

    title = AAV2_name + '_vs_' + CereAAV_name
    ax.set_title(f'{title}')
    ax.set_xlabel(f'Input_count', fontsize=12)
    ax.set_ylabel(f'Treatment_count', fontsize=12)
    ax.legend(
        bbox_to_anchor=(1.43, 1.0),
        loc='upper right',
        #borderaxespad=0,
        fontsize=10
    )
    #plt.show()
    if prefix:
        pngfile = os.path.join(outdir, lotID + prefix + '_freq_dist_scatter.png')
    else:
        pngfile = os.path.join(outdir,lotID + '_freq_dist_scatter.png')

    plt.savefig(pngfile, bbox_inches='tight')


def heatmap(lotID, Zscore, df, samples, outdir, prefix):
    logger.info('')
    logger.info('Start making heatmap..')

    # cheak num of samples ---
    if len(samples) < 2:
        logger.error('more than two samples')
        sys.exit()

    # cheak sample name ---


    if Zscore:
        zscore = st.zscore(df[[samples]], axis=1)
        df_data = pd.concat([df['ID'], pd.DataFrame(zscore, columns=samples)], axis=1)
        df_data = df_data.dropna()
    else:
        df_data = df[['ID'] + samples]

    # dendrogram ---
    data_array = df_data[[samples]]
    fig = ff.create_dendrogram(data_array, orientation='right')
    for i in range(len(fig['data'])):
        fig['data'][i]['yaxis'] = 'y2'

    y = fig['layout']['yaxis']['ticktext']
    y = list(map(int, y))
    x = data_array.columns
    z = df_data.iloc[y, :].reset_index(drop=True)

    # plot ---
    if Zscore:
        title = 'z-score'
    else:
        title = 'values'

    fig = go.Figure(
        data=go.Heatmap(
            colorbar=dict(
                title=title
            ),
            z=z.iloc[:, 1:],
            x=x,
            y=z['ID'],
            colorscale='Viridis'
        )
    )

    fig.update_layout(
        plot_bgcolor='white',
    )

    fig.update_xaxes(title='Samples(normalized)', tickfont=dict(size=10))
    fig.update_yaxes(title='Sequences',  tickfont=dict(size=5))

    fig.update_layout(
        title='heatmap',
        xaxis_nticks=36
    )
    #fig.show()
    if prefix:
        htmlfile = os.path.join(outdir, lotID + prefix + '_heatmap.html')
    else:
        htmlfile = os.path.join(outdir,lotID + '_heatmap.html')

    fig.write_html(htmlfile)
    logger.info('Finish making heatmap..')


def sorted_plot(lotID, df, y_names, y_title, outdir, prefix):
    logger.info('')
    logger.info('Start making sort_plot..')

    # cheak sample name ---

    main_title = y_name
    x = df.index + 1
    min_x = x.min()
    max_x = x.max()

    y = df[y_names]
    min_y = y.min() - 0.5
    max_y = y.max() + 0.5
    text = df['ID']

    fig = make_subplots()
    for y_name in y_names:
        fig.add_trace(
            go.Scatter(
                x=x,
                y=df[y_name].sort_values(ascending=False),
                mode='markers',
                marker=dict(
                    color='gray',
                    size=6.0,
                    line=dict(
                        color='black',
                        width=0.7
                    )
                ),
                text=text
            )
        )

    fig.update_layout(
        plot_bgcolor='white'
        #height=800,
        #width=900
    )
    fig.update_xaxes(
        title = dict(text='Sequences', font=dict(size=40)),
        showline=True,
        linewidth=1,
        linecolor='black',
        mirror=True,
        ticks='inside',
        tickfont=dict(size=15),
        range=(min_x, max_x)
    )
    fig.update_yaxes(
        title=dict(text=y_title, font=dict(size=40)),
        showline=True,
        linewidth=1,
        linecolor='black',
        mirror=True,
        ticks='inside',
        tickfont=dict(size=15),
        range=(min_y, max_y)
    )

    fig.for_each_xaxis(
        lambda axis: axis.title.update(
            font=dict(
                color='black',
                size=15
            )
        )
    )
    fig.for_each_yaxis(
        lambda axis: axis.title.update(
            font=dict(
                color='black',
                size=15
            )
        )
    )
    #fig.update_annotations(
    #    font=dict(size=30),
    #)
    fig.update_layout(
        title=dict(
            text=main_title,
            x=0.5,
            xanchor='center'
        ),
        showlegend=False
    )

    if prefix:
        htmlfile = os.path.join(outdir, lotID + prefix + '_plot.html')
    else:
        htmlfile = os.path.join(outdir,lotID + '_heatmap.html')

    fig.write_html(htmlfile)
    logger.info('Finish making sort_plot..')


def barcode_reps_scatter(df_bc, tissues):
    """_summary_

    Args:
        tissues (_type_): _description_

    tissues = ['Brain', 'Spinal_cord', 'Heart', 'Liver', 'Lung', 'Muscle']
    """
    import seaborn as sns
    from sklearn.linear_model import LinearRegression

    df = pd.concat([df_bc['BC_ID'], np.log2(df_bc.iloc[:, 2:] + 1)], axis=1)

    sns.set_style("white")

    for name in tissues:
        fig, ax = plt.subplots(figsize=(10, 7))

        x_label = f'cDNA#{name}#101_Norm'
        y_label = f'cDNA#{name}#102_Norm'

        sns.scatterplot(data=df, x=x_label, y=y_label, hue='Ascan', palette='gist_rainbow')

        r = df.loc[:, [x_label, y_label]].corr(method='pearson').iloc[0, 1]
        x = df[x_label].to_numpy()
        y = df[y_label].to_numpy()
        lm = LinearRegression().fit(x.reshape(-1, 1), y)
        ax.plot(x, lm.predict(x.reshape(-1, 1)), color='gray', linestyle='--', label='r='+str(round(r, 2)), alpha=0.4)

        xymax = max(np.max(x), np.max(y)) + 1 #xとyの最大値を定義
        ax.set_xlim([0, xymax])
        ax.set_ylim([0, xymax])

        ax.set_title(f'{name}')
        ax.set_xlabel(f'log2((rep1_Norm)+1)', fontsize=12)
        ax.set_ylabel(f'log2((rep2_Norm)+1)', fontsize=12)
        ax.legend(
            bbox_to_anchor=(1.2, 1.0),
            loc='upper right',
            #borderaxespad=0,
            fontsize=9
        )
        fname = 'PG4796_' + name + '_cDNA_rep12_BC_scatter.png'
        plt.savefig('A:/PROJECT/PG4796/900_NOUHIN/report_v03_230331/scatter/' + fname, bbox_inches='tight', dpi=300)
        plt.show()


def scatter_with_adjust_text():
    from adjustText import adjust_text

    file = 'A:/PROJECT/PG4796/900_NOUHIN/report_v03_230331/PG4796_cDNA_Ascan_Summary.xlsx'
    df_dna = pd.read_excel(file, sheet_name='cDNA')

    col_list = []
    for name in df_dna.columns:
        if 'Mean' in name:
            col_list.append(name)

    file = 'A:/PROJECT/PG4796/200_INFO/00_data/PG4796_normalizationGbeta.xlsx'
    df_Gbeta = pd.read_excel(file, sheet_name='Gbeta')
    df_Gbeta_mean = df_Gbeta.loc[:, ['Tissue', 'Mean']].groupby('Tissue').agg('mean').reset_index()

    Gbeta = {}
    for i in range(len(df_Gbeta_mean)):
        tissue = df_Gbeta_mean.loc[i, 'Tissue']
        Gbeta[tissue] = df_Gbeta_mean.loc[i, 'Mean']

    df_dna_mean = df_dna.loc[:, ['Ascan', 'AA_seq'] + col_list]

    for name in col_list:
    if not 'cDNA#Brain_Mean' in name:
        fig, ax = plt.subplots(figsize=(11.8, 8.6))
        x_label = 'cDNA#Brain_Mean'
        y_label = name
        x = np.log2(df_dna_mean[x_label]+1)
        y = np.log2(df_dna_mean[y_label]+1)
        xymax = max(np.max(x), np.max(y)) + 1 #xとyの最大値を定義
        ax.set_xlim([0, xymax])
        ax.set_ylim([0, xymax])
        ax.scatter(x, y, color='black', s=15.0)
        #ax.scatter(top100_xy[x_label], top100_xy[y_label], color='red', s=2.0, label='top100')

        x_name = re.compile(r'cDNA#(.*)_Mean').search(x_label).group(1)
        y_name = re.compile(r'cDNA#(.*)_Mean').search(y_label).group(1)
        ax.axvline(np.log2(Gbeta[x_name] * 1), color='gray', linestyle='--', label=f'{x_name}_Gbeta') # qPCRの値
        ax.axhline(np.log2(Gbeta[y_name] * 1), color='red', linestyle='--', label=f'{y_name}_Gbeta') # qPCRの値

        A = np.log2(Gbeta[y_name] * 1) / np.log2(Gbeta[x_name] * 1)
        binwidth = 0.05 #binの値幅を定義。
        bins = np.arange(0, xymax, binwidth) #binの個数を定義
        ax.plot(bins, A*bins, color='gray', linestyle='--', label=f'qPCR_FC(Brain=1)={round(A, 2)}')
        ax.plot(bins, (2*A)*bins, color='red', linestyle='--', label=f'qPCR_FC(Brain=1)={round((2*A), 2)}')

        #ax.plot([0, xymax], [0, xymax], color='gray', linestyle='--', label='FC=1')

        x_name = re.compile(r'cDNA#(.*)').search(x_label).group(1)
        y_name = re.compile(r'cDNA#(.*)').search(y_label).group(1)
        ax.set_xlabel(f'log2(({x_name})+1)', fontsize=12)
        ax.set_ylabel(f'log2(({y_name})+1)', fontsize=12)
        ax.legend(
            bbox_to_anchor=(1.3, 1.0),
            loc='upper right',
            #borderaxespad=0,
            fontsize=10
        )

        annotations = df_dna_mean['AA_seq'].tolist()
        texts = [plt.text(x[i], y[i], label, ha='center', va='center') for i, label in enumerate(annotations)]
        adjust_text(texts, arrowprops=dict(arrowstyle='->'))
        plt.tick_params(labelsize=12)
        #for i, label in enumerate(annotations):
        #    arrow_dict = dict(arrowstyle='-')
        #    ax.annotate(label, (x[i], y[i]), xytext=(x[i]+0.1, y[i]+0.2), size=10, arrowprops=arrow_dict)
        fname = 'PG4796_cDNA_scatter_' + y_name + '.png'
        plt.savefig('A:/PROJECT/PG4796/900_NOUHIN/report_v03_230331/Secondary_Analysis/' + fname, bbox_inches='tight')
        plt.show()

def paired_samples_dist(lotID, df, outdir, prefix):

    # set up ---
    title_size = 20
    xlab_size = 15
    ylab_size = 20

    np.random.seed(3)
    data = np.random.normal(size=(30, ))

    #Create the dataframe in a wide format with 'Before' and 'After ' as columns
    df = pd.DataFrame({'Before': data, 'After': data+1})

    #Merge dataframe from wide to long for sns.pointplot
    df_long = pd.melt(df, value_vars=['Before','After'])

    #Create a dataframe do display 4 conditions
    df_2 = pd.DataFrame({'Before': data, 'After': data+1, 'Before1': data, 'After1': data-1})

    #Merge dataframe from wide to long for sns.pointplot
    df_long_2 = pd.melt(df_2, value_vars=['Before', 'After', 'Before1', 'After1'])


    np.random.seed(3)
    df_jitter = pd.DataFrame(np.random.normal(loc=0, scale=0.05, size=df.values.shape), columns=df.columns)
    #Update the dataframe with adding a number based on the length on the columns. Otherwise all datapoints would be at the same x-axis location.
    df_jitter += np.arange(len(df.columns))

    # Create empty figure and plot the individual datapoints
    fig, ax = plt.subplots(figsize=(15, 9))

    for col in df:
        ax.plot(
            df_jitter[col],
            df[col],
            'o',
            alpha=0.8,
            zorder=2,
            ms=10,
            mew=1.5
        )

    for idx in df.index:
        ax.plot(
            df_jitter.loc[idx,['Before','After']],
            df.loc[idx,['Before','After']],
            color = 'gray',
            linewidth = 2,
            linestyle = '-',
            alpha = 0.2
        )

    for value in df_long_2:
        sns.violinplot(x='variable', y='value', data=df_long, hue = 'variable', split = True, inner = 'quartile', cut=1, dodge = True)
        sns.boxplot(x='variable', y='value', data=df_long, hue = 'variable', dodge = True, width = 0.2, fliersize = 2)

        #Additonal settings
        ax.set_xticks(range(len(df.columns)))
        ax.set_xticklabels((['Before', 'After']), size=xlab_size)
        ax.set_xlim(-1, len(df.columns))
        ax.set_ylabel('Value', size=ylab_size)
        ax.set_title('individual datapoints with lines, jitter, statistics, box- and violin', size = title_size)
        sns.despine()
        ax.legend_.remove()
        plt.setp(ax.collections, alpha=.1)

    #plt.show()
    if prefix:
        pngfile = os.path.join(outdir, lotID + prefix + '_paired_samples_dist.png')
    else:
        pngfile = os.path.join(outdir, lotID + '_paired_samples_dist.png')

    plt.savefig(pngfile, bbox_inches='tight')
