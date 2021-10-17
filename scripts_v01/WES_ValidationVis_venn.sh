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

snv_dir=${workdir}/260_Somatic_Call/04_snpeff/pass_variants
cnv_dir=${work_dir}/270_CNV/02_make_table

#--------------------------------------
# directory
#--------------------------------------
log_dir=${out_dir}/log
png_dir=${out_dir}/png

mkdir -p ${log_dir}
mkdir -p ${png_dir}

#--------------------------------------
# set up log file
#--------------------------------------
log_file=${log_dir}/${prefix}_${day}.log

#--------------------------------------
# run process
#--------------------------------------
echo [`date '+%Y-%m-%d %H:%M:%S'`]  run process ..

#===================  make venn   ===================#
qsub -V -q ${queue} -o $log_file -e $log_file -j y -N ${prefix} -pe smp 32 \
  python3 ${script_dir}/04_venn_plot.py -i ${snv_dir} -l ${legend[@]} -y 1 -x1 ${x1[@]} -x2 ${x2[@]} -o ${png_dir} -n ${file_name}

qsub -V -q ${queue} -o $log_file -e $log_file -j y -N ${prefix} -pe smp 32 \
  python3 ${script_dir}/04_venn_plot.py -i ${cnv_dir} -l ${legend[@]} -y 2 -x1 ${x1[@]} -x2 ${x2[@]} -o ${png_dir} -n ${file_name}

#--------------------------------------
# move log file
#--------------------------------------
mv ${out_dir}/*.log ${log_dir}
mv ${png_dir}/*.log ${log_dir}

echo [`date '+%Y-%m-%d %H:%M:%S'`] done.

#
# End
#