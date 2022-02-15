#!/usr/bin/bash

export HbF=${HOME}/COVID-19/HbF
export pgwas=~/rds/results/public/proteomics
export M=1e6

function annotate()
{
Rscript -e '
library(pQTLtools)
head(hg19)
#subset(hg19,grepl("^HBG|^HBA|^HBB",SYMBOL))
#        chr   start     end  width strand      ENSEMBL SYMBOL      UNIPROT
#43431 chr11 5246696 5248301   1606      - ENSG0000....    HBB D9YZU5, ....
#43434 chr11 5269502 5271087   1586      - ENSG0000....   HBG1 D9YZU8, ....
#43435 chr11 5269502 5276395   6894      - ENSG0000....   HBG2 D9YZU9, ....
#43436 chr11 5274421 5276395   1975      - ENSG0000....   HBG2 D9YZU9, ....
#43437 chr11 5275522 5667011 391490      - ENSG0000....   HBG2 D9YZU9, ....
#56898 chr16  222846  223709    864      + ENSG0000....   HBA2 D1MGQ2, ....
#56899 chr16  226679  227520    842      + ENSG0000....   HBA1 D1MGQ2, ....
#subset(hg19Tables,grepl("^HBG|^HBA|^HBB",hgncSym))
#      X.chrom chromStart chromEnd strand    acc uniprotName
#11496   chr11    5246830  5248251      - P68871   HBB_HUMAN
#11498   chr11    5269591  5275958      - P69891  HBG1_HUMAN
#11499   chr11    5274509  5275958      - P69892  HBG2_HUMAN
#subset(SomaScanV4.1,grepl("^HBG|^HBA|^HBB",GeneID))
#        #     SeqID    Human.Target.or.Analyte    UniProt.ID   GeneID    Type
#2848 2848   4915-64                 Hemoglobin P69905|P68871 HBA1|HBB Protein
#2849 2849 17137-160    Hemoglobin subunit beta        P68871      HBB Protein
#2852 2852   19774-8 Hemoglobin subunit gamma-2        P69892     HBG2 Protein
'

grep -e HBG ${HbF}/work/hbf_hits.txt
grep -n HBG ${HbF}/work/*tsv | cut -f1 | xargs -l basename | sed 's/\./\t/' | cut -f1 | uniq
grep -e 4915.64 -e 17137.160 -e 19774.8 ${HbF}/work/deCODE.tsv
#https://www.ebi.ac.uk/gwas/studies/GCST003122
}

export deCODE=${pgwas}/deCODE
