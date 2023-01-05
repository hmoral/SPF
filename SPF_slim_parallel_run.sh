#!/bin/bash
SEEDparallel=`echo $RANDOM`
DATE=$(date '+%d_%m_%Y_%H_%M_%S')

Nin=$1
DIRout=$2
outpref=$3
outPref_in=`echo ${DIRout}/${outpref}`
adapP=$4
adapU=$5

CORES=$6
REPS=$7

mkdir -p ${DIRout}

U=0.4
opt1=0.2
optDif=1
collapseN=28

for replicate in $(seq $REPS);do
seed=`echo $RANDOM`
echo ./launch_SPF_slim.sh ${Nin} ${opt1} ${optDif} ${adapP} ${seed} ${collapseN} ${outPref_in} ${U} ${adapU} >> launch_jobs_SPF_slim_${DATE}_${SEEDparallel}.txt
done

mkdir -p tmp_parallel
parallel -j ${CORES} --tmpdir tmp_parallel < launch_jobs_SPF_slim_${DATE}_${SEEDparallel}.txt

