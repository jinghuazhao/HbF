#!/usr/bin/bash

#SBATCH --job-name=deCODE
#SBATCH --account CARDIO-SL0-CPU
#SBATCH --partition cardio
#SBATCH --qos=cardio
#SBATCH --array=1-4907
#SBATCH --mem=10800
#SBATCH --time=1-00:00:00
#SBATCH --output=/rds/user/jhz22/hpc-work/work/_deCODE_%A_%a.out
#SBATCH --error=/rds/user/jhz22/hpc-work/work/_deCODE_%A_%a.err
#SBATCH --export ALL

export TMPDIR=/rds/user/jhz22/hpc-work/work

module load gcc/6

export HbF=${HOME}/COVID-19/HbF
export pgwas=~/rds/results/public/proteomics
export M=0
export p_gwas=1
export deCODE=${pgwas}/deCODE

function make_list()
{
  ls ${deCODE}/*gz | xargs -l basename -s .txt.gz > ${HbF}/work/deCODE.lst
}

export id=$(awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]' ${HbF}/work/deCODE.lst)

function extract()
{
  (
  sed '1d' ${HbF}/work/hbf_hits.txt | cut -f1,3,4,12,13 | grep -v NA | sed 's/, /;/g' | tr '\t' ' ' | \
  parallel -C' ' -j1 --env deCODE --env id '
    export chr={1}
    export pos={2}
    export rsid={3}
    export snpid={4}
    export gene={5}
    export region=$(echo -e "chr${chr}:${pos}-${pos}")
    tabix ${deCODE}/${id}.txt.gz ${region} | \
    awk -v rsid=${rsid} -v snpid=${snpid} -v gene=${gene} -v id=${id} -v p=${p_gwas} -v OFS="\t" "rsid==\$4&&\$8<=p{print rsid,snpid,gene,id,\$0}"
  '
  ) > ${HbF}/deCODE/deCODE-${id}
}

extract
