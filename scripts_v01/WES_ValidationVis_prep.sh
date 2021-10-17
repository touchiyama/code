#!/bin/bash -eu

# run_performance.sh
# Created on 2021-10-13
# Last modified on 2021-10-11
# by Tomoya Uchiyama
# v1.0
# © 2021/10/15 Takara Bio Inc.

day=`date '+%Y_%m_%d_%H_%M_%S'`

#--------------------------------------
# set up LIBRARY_PATH (have to modify)
#--------------------------------------
export PATH=/NGSWORK/NGS/DOCUMENT/SOP/GA/SOP-42_Validation/program/wes_validation/environment/miniconda/envs/py38/bin/:$PATH
script_dir=/NGSWORK/PROJECT
#base_dir=$(cd $(dirname ${BASH_SOURCE:-$0}); pwd)

#--------------------------------------
# get parameter
#--------------------------------------
function usage {
cat <<EOM

usage: $0 -l LotID -s sample1,sample2,.. -t VALUE1,VALUE2,.. -n VALUE1,VALUE2,..

optional arguments:
    -h                        show this help message and exit
    -l [LotID]                input lotID
    -s [sample1,sample2,..]   input sample name based on the experiment
    -t [VALUE1,VALUE2,..]     input tumor_data based on the experiment
    -n [VALUE1,VALUE2,..]     input normal_data based on the experiment

EOM
exit 1
}

#--------------------------------------
# show and set up parameter
#--------------------------------------
if [ $# = 0 ]; then
  usage
  exit 1
fi

while getopts ":l:s:t:n:h" optKey; do
  case "$optKey" in
    l)
      echo "-l = ${OPTARG}"
      lotID=${OPTARG}
      ;;
    s)
      echo "-s = ${OPTARG}"
      sample=${OPTARG}
      ;;
    t)
      echo "-t = ${OPTARG}"
      tumor_data=${OPTARG}
      ;;
    n)
      echo "-n = ${OPTARG}"
      normal_data=${OPTARG}
      ;;
    h)
      usage
      ;;
    *)
      echo " "
      echo "invalid option or syntax: ${OPTARG}"
      usage
      ;;
  esac
done

work_dir=/NGSWORK/PROJECT/${lotID}/200_GENOME_DNA-ABC
out_dir=/NGSWORK/PROJECT/${lotID}/240_GENOME_WES_Vis
prefix=$(basename ${0%.*})
queue=-all.q

legend=(${sample//,/ })
x1=(${tumor_data//,/ })
x2=(${normal_data//,/ })

#=================== setUp file_name  ===================#
name_arr=()
for name in ${legend[@]}
do
  if [[ ${name} =~ ^([A-Z][0-9]{2})_.* ]]
  then
    name_arr+=("${BASH_REMATCH[1]}")
  fi
done
uniq_name=$( printf "%s\n" "${name_arr[@]}" | sort -u )
file_info=`echo ${uniq_name} | tr -d ' '`
file_name=${lotID}${file_info}

#=============  SNV+INDEL-related parameter =============#
snv_dir=${workdir}/260_Somatic_Call/04_snpeff
passVariants_dir=${snv_dir}/pass_variants # have to change
snv_summary=${lotID}_snvIndel_summary.xlsx
markerVars_table=${workdir}/240_GENOME_WES_Vis/${lotID}_markerVars_table.xlsx
bam_dir=${work_dir}/250_MUTATION/DRAGEN

#===============   CNV-related parameter  ===============#
cnv_dir=${work_dir}/270_CNV/02_make_table
cnvCorr_dir=${work_dir}/270_CNV/03_cnv_summary
OncoScan_file=${work_dir}/240_GENOME_WES_Vis/${lotID}_CNV.xlsx # have to change

#--------------------------------------
# directory
#--------------------------------------
log_dir=${out_dir}/log
json_dir=${out_dir}/json
xlsx_dir=${out_dir}/xlsx

mkdir -p ${log_dir}
mkdir -p ${json_dir}
mkdir -p ${xlsx_dir}

#--------------------------------------
# set up log file
#--------------------------------------
log_file=${log_dir}/${prefix}_${day}.log

#--------------------------------------
# run process
#--------------------------------------
echo [`date '+%Y-%m-%d %H:%M:%S'`]  run process ..

#==============  preparation for SNV+INDEL ==============#
python3 ${script_dir}/run_vcftools_indelSnv.py -i ${snv_dir} -o ${passVariants_dir}
sh ${script_dir}/run_vcftools_indelSnv_`date '+%Y%m%d'`.sh
python3 ${script_dir}/make_table.py -l ${lotID} -i ${pass_variants} -o ${xlsx_dir} -n ${snv_summary}

#============= make SNV+INDEL with F1 table =============#
python3 ${script_dir}/03_find_mut_with_F1.py -t ${markerVars_table} -i ${snv_dir} -l ${legend[@]} -b ${bam_dir} -o ${xlsx_dir} -n ${lotID}

#==============  preparation for CNV  ===================#
python3 ${script_dir}/07_calc_cnvCorr.py -i ${cnv_dir} -r ${OncoScan_file} -o ${xlsx_dir} -n ${lotID}

#============= make CNV linearity json file =============#
python3 ${script_dir}/06_plot_cnv_linearity.py -i ${cnv_dir} -x1 10 15 20 30 40 -o ${json_dir} -n ${lotID}

#--------------------------------------
# move log file
#--------------------------------------
mv ${out_dir}/*.log ${log_dir}
mv ${json_dir}/*.log ${log_dir}
mv ${xlsx_dir}/*.log ${log_dir}

echo [`date '+%Y-%m-%d %H:%M:%S'`] done.

#
# End
#

#-rwxr-x--- 1 uchiyamat dgc 11K 10月 12 21:22 WES_valdation_Vis/08_plot_cnvCorr.py

業務目標

TSO500、WESバリデーション
製薬企業向けの新規サービスの情報解析を1件

社内のプログラム開発規約に則り、
作業効率化に寄与するプログラム開発を第4Qまでに2件)以上担当する