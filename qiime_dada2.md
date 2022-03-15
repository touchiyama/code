***
## qiime2を通じたdada2の利用 (トリミング長の検討)
***

### 参考文献：<br>
>・[qiime2デモ解析](https://rpubs.com/nishikosh/qiime2_demo) <br>
>・[微生物群集解析ツールQIIME2を使う.Part2.(QC, FeatureTable生成編)](https://qiita.com/kohei-108/items/547b2fbdf28fb04c28ea) <br>
>・[Qiime2 を用いた 16S rRNA 菌叢解析](https://qiita.com/keisuke-ota/items/6399b2f2f7459cd9e418) <br>

### ・qiime2インストール <br>
> ```bash
> # qiime2をダウンロード
> $ wget https://data.qiime2.org/distro/core/qiime2-2021.4-py38-osx-conda.yml
>
> # 仮想環境の作成(yamlファイルの読み込み)
> $ conda env create -n qiime2-2021.4 --file qiime2-2021.4-py38-osx-conda.yml
>
> ```
### ・qiime2の起動 <br>
> ```bash
> # 仮想環境の実行
> $ conda activate qiime2-2021.4
>
> # qiime2のhelp(インストールされているか確認)
> $ qiime --help
> $ qiime --version
> $ qiime tools import --help
>
> # tab補完を有効にする
> $ source tab-qiime
>
> # 仮想環境を抜ける
> $ conda deactivate
> ```

### ・fastqファイルのインポート <br>
> #### manifiest.csvの作成
> ```bash
> # ヘッダ行の書き出し
> echo sample-id","absolute-filepath","direction > manifest.csv
>
> # manifest.csvにサンプルID, fastqファイルパス,ストランドを追記する
> fqs=(*_R1_001.fastq) # fastqファイル名の配列
> wd=`pwd` # 現在いるディレクトリの絶対パス
>
> for r1 in ${fqs[@]}
> do
>   r2=${r1/_R1/_R2} # R2ファイル名
>   fqid=(${r1//_/ }) # サンプルID(fastqファイル名を_で分割して配列に入れた)
>   cpfq_r1=${wd}/${r1}  # R1の絶対パス
>   cpfq_r2=${wd}/${r2}  # R2の絶対パス
>   echo -e ${fqid}","${cpfq_r1}","forward >> manifest.csv # manifest.csvに追記
>   echo -e ${fqid}","${cpfq_r2}","reverse >> manifest.csv
> done
>
> ```
> #### fastqファイルのインポート
> ```bash
> $ qiime tools import \
> --type SampleData[PairedEndSequencesWithQuality] \
> --input-path manifest.csv \
> --output-path sequence.qza \
> --input-format PairedEndFastqManifestPhred33
> ```
> #### output <br>
> ```bash
> # Imported manifest.csv as PairedEndFastqManifestPhred33 to sequence.qza
> ```

### ・qzaファイルの可視化
> 可視化用のファイルqzvの作成 <br>
> ```bash
> $ qiime demux summarize \
> --i-data sequence.qza \
> --o-visualization sequence.qzv
> ```
> qzvファイルを以下のリンク先で表示。<br>
> [https://view.qiime2.org/](https://view.qiime2.org/)

### ・feature tableの作成
> ```bash
> $ qiime feature-table summarize \
> --i-data sequence.qza \
> --o-visualization sequence.qzv \
> --m-sample-metadata-file feature_table.tsv
> ```

### リードのダウンサンプリング
> インポートしたfastqファイルから5%をランダムにサンプリングする <br>
> ```bash
> $ qiime demux subsample-paired \
> --i-sequences sequence.qza \
> --p-fraction 0.05 \
> --o-subsampled-sequences subsample.qza
> ```
>
> 可視化用のファイルqzvの作成 <br>
> ```bash
> $ qiime demux summarize \
> --i-data subsample.qza \
> --o-visualization subsample.qzv
> ```

### ・DADA2処理 (トリミング長の検討)
> ・クオリティが20以下の部分をカットする <br>
> ・--p-trim-left-f x：各シーケンスのforward配列の初めから塩基長xまでを削除する．<br>
> ・--p-trim-left-r y：各シーケンスのReverse配列の初めから塩基長yまでを削除する．<br>
> ・--p-trunc-len-f m：各シーケンスのforward配列から採用する最大配列長mを決定する． mを超える塩基長は除去する．<br>
> ・--p-trunc-len-r n：各シーケンスのreverse配列から採用する最大配列長nを決定する． nを超える塩基長は除去する．<br>
>
> ・PCR産物が約500bpなので、最低でも20bp以上のオーバーラップ領域を加味して、F+Rで520を超える長さで実行する．<br>
> ・f > R かつ F + R = 520 を満たすような条件でトリミングする <br>
> ・non-chimericの割合が最大となる時のトリミング条件を探る <br>
>
> ```bash
> $ for i in {260, 270, 280}; do \
> qiime dada2 denoise-paired \
> --verbose \
> --p-n-threads 4 \
> --p-trim-left-f 0 \
> --p-trim-left-r 0 \
> --p-trunc-len-f $i \
> --p-trunc-len-r 520-$i \
> --i-demultiplexed-seqs subsample.qza \
> --o-table table-dada2-${i}.qza \
> --o-representative-sequences rep-seqs-dada2-${i}.qza \
> --o-denoising-stats stats-dada2-${i}.qza \
> ; done
> ```
>
> このステップで代表配列ファイルであるrep-seqs-dada2-XXX.qza(塩基配列データ)と、サンプル-ASV配列の行列ファイルであるtable-dada2-XXX.qza(リード数データ)が出力される。<br>
>
> ・rep-seqs.qza <br>
> データセット中に存在する塩基配列をまとめたファイル。識別用のIDも割り当てられている。<br>
>
> ・table.qza <br>
> データセット中に存在する塩基配列のリード数をまとめたファイル。各塩基配列はIDによって区別されており、塩基配列の中身は記録されていない。<br>
>
> #### feature tableの作成
> ```bash
> $ qiime feature-table summarize \
> --i-data sequence.qza \
> --o-visualization sequence.qzv \
> --m-sample-metadata-file feature_table.tsv
> ```
> #### データのエクスポート
> ```bash
> $ qiime tools export \
> --input-path table.qza \
> --output-path ./output
> ```
> ```bash
> $ qiime tools export \
> --input-path rep-seqs.qza \
> --output-path ./output
> ```

***
## Trimmomaticを介したdada2の利用 (トリミング長の検討)
***

・250bpがトリミング長として最適であるかを確かめるために、トリミング長変化に伴うnon-chimericリードの割合の変化を調べる。<br>

・主にやること：Trimmomaticでトリムしたfastqファイルをinputとして、DADA2にかける <br>


(1) minlenを210,220,230,240,250,260,270,300と設定し, QV20とQV30で分ける <br>
(2) Trimmomatic終了後、以下のコマンドでダウンサンプリングを行う。<br>
```bash
$ seqkit sample -p 0.05 -s 100000 bud_L001.fastq.gz -o XXXXX_210_sample_trimed.fastq.gz
```
(3) DADA2でfastqファイルのインポート（QV20とQV30で分ける）<br>
・サンプル毎にmanifest.csvの作成 <br>
  > sample-idを「PG4600_01_a_qv20_210」のようにする <br>

・fastqファイルのインポート <br>
```bash
$ qiime tools import \
--type SampleData[PairedEndSequencesWithQuality] \
--input-path PG4600_01_a_manifest.csv \
--output-path PG4600_a_trimed_seq.qza \
--input-format PairedEndFastqManifestPhred33
```

(4) DADA2処理（QV20とQV30で分ける）<br>
・FとRの長さを、Trimmomaticのminlenに合わせる <br>

(5) stats summary tableを作り、可視化 <br>


