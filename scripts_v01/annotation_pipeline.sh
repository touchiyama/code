#!/bin/bash -eux
# Date: 2021/03/19
# Author: Shogo Satoyama
# 実行書式:
# ./annotation_pipeline.sh \
#  PCXXXX \ # Lot ID
#  PCXXXX_01_a,PCXXXX_02_a \ # sanple_ids
#  2 \ # 円グラフの列数
#  2>&1 | tee annotation_pipeline.log
#======================================

#--------------------------------------
# parameter
#--------------------------------------

lot_id=$1
# PCXXXX
samples=$2
# PCXXXX_01_a,PCXXXX_02_a, ...
# スペースなしでカンマ区切りでサンプルIDを入力
sample_ids=(${samples//,/ }) # カンマで区切って配列化
ncol=$3
work_dir=`pwd`

#--------------------------------------
# directory
#--------------------------------------

echo [`date '+%Y-%m-%d %H:%M:%S'`] making directories ...

mkdir -p ${work_dir}/annotation
mkdir -p ${work_dir}/sample_metrics
mkdir -p ${work_dir}/count_insertions
mkdir -p ${work_dir}/annotation_summary
mkdir -p ${work_dir}/plot
mkdir -p ${work_dir}/report
mkdir -p ${work_dir}/report/html

echo [`date '+%Y-%m-%d %H:%M:%S'`] Done.

#--------------------------------------
# Annotation
#--------------------------------------

echo [`date '+%Y-%m-%d %H:%M:%S'`] annotating ...

cd ${work_dir}/annotation

# anootation
for sample_id in ${sample_ids[@]};
do
    Rscript --slave --vanilla \
    /NGSWORK/NGS/DOCUMENT/SOP/GA/SOP-46-Clonality/program/script/annotation/annotatr.R \
    ${work_dir}/${sample_id}/${sample_id}.Results_Homo_Annotated_Sorted.bed \
    ${sample_id} \
    /NGSWORK/NGS/DOCUMENT/SOP/GA/SOP-19_CancerAnalysis/cosmic/grch38/cosmic/v84/cancer_gene_census.csv \
    2>&1 | tee -a ${sample_id}_annotatr.log &
done
wait # 並列実行

# tsvをhtmlに変換
for sample_id in ${sample_ids[@]};
do
    Rscript --slave --vanilla \
    /NGSWORK/NGS/DOCUMENT/SOP/GA/SOP-46-Clonality/program/script/annotation/table2html_top100.R \
    ${work_dir}/annotation/${sample_id}_reads_top100.tsv ${work_dir}/annotation/${sample_id}_reads_top100.html 2>&1 | tee -a ${sample_id}_top100.log &
done
wait # 並列実行

# tsvをxlsxに変換（アノテーションが1行にまとめて記載）
for sample_id in ${sample_ids[@]};
do
    python3 \
    /NGSWORK/NGS/DOCUMENT/SOP/GA/SOP-46-Clonality/program/script/annotation/excelconverter.py \
    -i /${work_dir}/annotation/${sample_id}_annotation_1row.tsv \
    -o /${work_dir}/annotation/${sample_id}_annotation.xlsx | tee -a ${sample_id}_convert_excel_1row.log &
done
wait # 並列実行

# anootation(multi hit reads)
for sample_id in ${sample_ids[@]};
do
    Rscript --slave --vanilla \
    /NGSWORK/NGS/DOCUMENT/SOP/GA/SOP-46-Clonality/program/script/annotation/annotatr.R \
    ${work_dir}/${sample_id}/${sample_id}.repeats.Results_Homo_Annotated_Sorted.bed \
    ${sample_id}_multiplehit \
    /NGSWORK/NGS/DOCUMENT/SOP/GA/SOP-19_CancerAnalysis/cosmic/grch38/cosmic/v84/cancer_gene_census.csv \
    2>&1 | tee -a ${sample_id}_repeats_annotatr.log &
done
wait # 並列実行

# tsvをhtmlに変換(multi hit reads)
for sample_id in ${sample_ids[@]};
do
    Rscript --slave --vanilla \
    /NGSWORK/NGS/DOCUMENT/SOP/GA/SOP-46-Clonality/program/script/annotation/table2html_top100.R \
    ${work_dir}/annotation/${sample_id}_multiplehit_reads_top100.tsv ${work_dir}/annotation/${sample_id}_multiplehit_reads_top100.html 2>&1 | tee -a ${sample_id}_repeats_top100.log &
done
wait # 並列実行

# tsvをxlsxに変換（アノテーションが1行にまとめて記載）(multi hit reads)
for sample_id in ${sample_ids[@]};
do
    python3 \
    /NGSWORK/NGS/DOCUMENT/SOP/GA/SOP-46-Clonality/program/script/annotation/excelconverter.py \
    -i /${work_dir}/annotation/${sample_id}_multiplehit_annotation_1row.tsv \
    -o /${work_dir}/annotation/${sample_id}_multiplehit_annotation.xlsx | tee -a ${sample_id}_convert_excel_multiplehit_1row.log &
done
wait # 並列実行


cd ${work_dir}

echo [`date '+%Y-%m-%d %H:%M:%S'`] Done.

#--------------------------------------
# Sample Metrics
#--------------------------------------

echo [`date '+%Y-%m-%d %H:%M:%S'`] getting sample metrics ...

cd ${work_dir}/sample_metrics

/NGSWORK/NGS/GridEngine/LANG/Python/Python-3.7.7/bin/python3 \
/NGSWORK/NGS/DOCUMENT/SOP/GA/SOP-46-Clonality/program/script/annotation/get_sample_metrics.py \
${lot_id} \
${work_dir} \
${samples} \
2>&1 | tee -a get_sample_metrics.log

Rscript --slave --vanilla \
/NGSWORK/NGS/DOCUMENT/SOP/GA/SOP-46-Clonality/program/script/annotation/table2html_analysis_summary_v2.R \
${lot_id}_map_summary.tsv \
${lot_id}_map_summary.html \
2>&1 | tee -a make_map_summary.log

cd ${work_dir}

echo [`date '+%Y-%m-%d %H:%M:%S'`] Done.

#--------------------------------------
# count insertions
#--------------------------------------

echo [`date '+%Y-%m-%d %H:%M:%S'`] Counting insertions per chromosome ...

cd ${work_dir}/count_insertions

/NGSWORK/NGS/GridEngine/LANG/Python/Python-3.7.7/bin/python3 \
/NGSWORK/NGS/DOCUMENT/SOP/GA/SOP-46-Clonality/program/script/annotation/count_chr_insertions.py \
${lot_id} \
${work_dir} \
${samples} \
2>&1 | tee -a count_chr_insertions.log

Rscript --slave --vanilla \
/NGSWORK/NGS/DOCUMENT/SOP/GA/SOP-46-Clonality/program/script/annotation/table2html.R \
${lot_id}_count_summary.tsv \
${lot_id}_count_summary.html

Rscript --slave --vanilla \
/NGSWORK/NGS/DOCUMENT/SOP/GA/SOP-46-Clonality/program/script/annotation/table2html.R \
${lot_id}_multiplehit_count_summary.tsv \
${lot_id}_multiplehit_count_summary.html

cd ${work_dir}

echo [`date '+%Y-%m-%d %H:%M:%S'`] Done.

#--------------------------------------
# Annotation Summary
#--------------------------------------

echo [`date '+%Y-%m-%d %H:%M:%S'`] Making annotation summary ...

cd ${work_dir}/annotation_summary

/NGSWORK/NGS/GridEngine/LANG/Python/Python-3.7.7/bin/python3 \
/NGSWORK/NGS/DOCUMENT/SOP/GA/SOP-46-Clonality/program/script/annotation/make_annotation_summary.py \
${lot_id} \
${work_dir}/annotation \
${samples} \
2>&1 | tee -a make_annotation_summary.log 

cd ${work_dir}

echo [`date '+%Y-%m-%d %H:%M:%S'`] Done.

#--------------------------------------
# plot
#--------------------------------------

echo [`date '+%Y-%m-%d %H:%M:%S'`] Ploting ...
cd ${work_dir}/plot

# # barchart
# /NGSWORK/NGS/GridEngine/LANG/Python/Python-3.7.7/bin/python3 \
# /NGSWORK/NGS/DOCUMENT/SOP/GA/SOP-46-Clonality/program/script/annotation/barchart.py \
# ${lot_id} \
# ${work_dir}/annotation_summary/${lot_id}_annotation_summary.xlsx \
# ${lot_id}_annotation_summary.html \
# 2>&1 | tee -a barchart.log 

# /NGSWORK/NGS/GridEngine/LANG/Python/Python-3.7.7/bin/python3 \
# /NGSWORK/NGS/DOCUMENT/SOP/GA/SOP-46-Clonality/program/script/annotation/barchart.py \
# ${lot_id} \
# ${work_dir}/annotation_summary/${lot_id}_multiplehit_annotation_summary.xlsx \
# ${lot_id}_multiplehit_annotation_summary.html \
# 2>&1 | tee -a barchart.log 

# piechart
/NGSWORK/NGS/GridEngine/LANG/Python/Python-3.7.7/bin/python3 \
/NGSWORK/NGS/DOCUMENT/SOP/GA/SOP-46-Clonality/program/script/annotation/piechart_read_annots.py \
piechart.html \
${samples} \
${work_dir}/annotation \
${ncol} \
${lot_id} \
2>&1 | tee -a piechart_read_annots.log 

/NGSWORK/NGS/GridEngine/LANG/Python/Python-3.7.7/bin/python3 \
/NGSWORK/NGS/DOCUMENT/SOP/GA/SOP-46-Clonality/program/script/annotation/multiplehit_piechart_read_annots.py \
multiplehit_piechart.html \
${samples} \
${work_dir}/annotation \
${ncol} \
${lot_id} \
2>&1 | tee -a piechart_read_annots.log 

cd ${work_dir}

# echo [`date '+%Y-%m-%d %H:%M:%S'`] Done.

# #--------------------------------------
# # summary report
# #--------------------------------------
# echo [`date '+%Y-%m-%d %H:%M:%S'`] Making ${lot_id} summary report ...
# cd ${work_dir}/report

# /NGSWORK/NGS/GridEngine/LANG/Python/Python-3.7.7/bin/python3 \
# /NGSWORK/NGS/DOCUMENT/SOP/GA/SOP-46-Clonality/program/script/annotation/make_report.py \
# /NGSWORK/NGS/DOCUMENT/SOP/GA/SOP-46-Clonality/program/template/summary_report \
# ${work_dir}/report \
# ${lot_id}_summary_report.html \
# ${samples} \
# ${lot_id} \
# 2>&1 | tee -a make_report.log 

# # copy html files
# for sample_id in ${sample_ids[@]};
# do
#     cp -f ${work_dir}/annotation/${sample_id}_reads_top100.html ${work_dir}/report/html &
# done
# wait # 並列実行

# cp -f ${work_dir}/sample_metrics/${lot_id}_map_summary.html ${work_dir}/report/html
# cp -f ${work_dir}/count_insertions/${lot_id}_count_summary.html ${work_dir}/report/html
# cp -f ${work_dir}/plot/${lot_id}_annotation_summary.html ${work_dir}/report/html
# cp -f ${work_dir}/plot/${lot_id}_hg38_cpg_islands_pie_chart.html ${work_dir}/report/html
# cp -f ${work_dir}/plot/${lot_id}_hg38_genes_exons_pie_chart.html ${work_dir}/report/html
# cp -f ${work_dir}/plot/${lot_id}_hg38_genes_intronexonboundaries_pie_chart.html ${work_dir}/report/html
# cp -f ${work_dir}/plot/${lot_id}_hg38_genes_promoters_pie_chart.html ${work_dir}/report/html

# cd ${work_dir}

echo [`date '+%Y-%m-%d %H:%M:%S'`] Done.

#
# END
#