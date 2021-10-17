#!/bin/bash -eu

# dragen_somatic_pipeline.sh
# Created on 2021-09-13
# Last modified on 2021-10-11
# by Shogo Satoyama
# v2.0
# © 2021/09/13 Takara Bio Inc.

export PATH=/NGSWORK/NGS/DOCUMENT/SOP/GA/SOP-42_Validation/program/wes_validation/environment/miniconda/envs/py38/bin/:$PATH

base_dir=$(cd $(dirname ${BASH_SOURCE:-$0}); pwd)
source $base_dir/set_env.sh
day=`date '+%Y_%m_%d_%H_%M_%S'`

prefix=${2%.*}
queue=${3:-all.q}
if [ -z "${4+UNDEF}" ]; then
    log_dir=$(dirname $2)
    log_file=$log_dir/${prefix}_${1,,}_${day}.log
else
    log_file=$4
fi

qsub -V -q ${queue} -o $log_file -e $log_file -j y -N ${prefix}_$1 -pe smp 32 \
    python3 $base_dir/dragen_somatic_pipeline_v1.py $1 --yaml $2 --local-scheduler

qsub -V -q ${queue} -o command_list_${day}.log -e command_list_${day}.log -j y -hold_jid ${prefix}_$1 -N ${prefix}_$1_Extract_Commands \
    grep -e '"'"^INFO: qsub"'"' -e '"'"^INFO: cp"'"' $log_file
#
# End
#