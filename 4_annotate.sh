#!/usr/bin/bash

export HbF=${HOME}/COVID-19/HbF
export pgwas=~/rds/results/public/proteomics
export M=1e6

function cis()
(
 for cis in $(awk '$1==11||$1==16 {print $4}' ${HbF}/work/hbf_hits.txt)
  do
    for study in AGES ARIC deCODE Fenland scallop-cvd1 INTERVAL LBC1936 GTEx eQTL
    do
       export f=${HbF}/work/${study}.tsv
       cat <(echo ${cis} ${study}) <(head -1 ${f} | cut -f1-3 --complement) <(cut -f1-3 --complement ${f} | grep -w ${cis})
    done
  done
) > ${HbF}/work/cis.tsv


function all()
(
 for rsid in $(awk 'NR>1{print $4}' ${HbF}/work/hbf_hits.txt)
  do
    for study in AGES ARIC deCODE Fenland scallop-cvd1 INTERVAL LBC1936 GTEx eQTL
    do
       export f=${HbF}/work/${study}.tsv
       case ${study} in
       AGES)
         awk -v study=${study} -v rsid=${rsid} -v OFS='\t' '$1 ~ rsid {print $1,study,$6}' ${f} | \
         sort -k2,2 | uniq
         ;;
       ARIC)
         awk -v study=${study} -v rsid=${rsid} -v OFS='\t' '$1 ~ rsid {print $1,study,$6}' ${f} | \
         sort -k2,2 | uniq
         ;;
       deCODE)
         awk -v study=${study} -v rsid=${rsid} -v OFS='\t' '$1 ~ rsid {print $1,study,$4}' ${f} | \
         sort -k2,2 | uniq
         ;;
       Fenland)
         awk -v study=${study} -v rsid=${rsid} -v OFS='\t' '$1 ~ rsid {print $1,study,$14}' ${f} | \
         sort -k2,2 | uniq
         ;;
       scallop-cvd1)
         awk -v study=${study} -v rsid=${rsid} -v OFS='\t' '$1 ~ rsid {print $1,study,$5}' ${f} | \
         sort -k2,2 | uniq
         ;;
       INTERVAL)
         awk -v study=${study} -v rsid=${rsid} -v OFS='\t' '$1 ~ rsid {print $1,study,$4}' ${f} | \
         sort -k2,2 | uniq
         ;;
       LBC1936)
         awk -v study=${study} -v rsid=${rsid} -v OFS='\t' '$1 ~ rsid {print $1,study,$5}' ${f} | \
         sort -k2,2 | uniq
         ;;
       GTEx)
         awk -v study=${study} -v rsid=${rsid} -v OFS='\t' '$1 ~ rsid {print $1,study,$4}' ${f} | \
         sort -k2,2 | uniq
         ;;
       eQTL)
         awk -v study=${study} -v rsid=${rsid} -v OFS='\t' '$1 ~ rsid {print $1,study,$4}' ${f} | \
         sort -k2,2 | uniq
         ;;
       *)
         ;;
       esac
    done
  done
) > ${HbF}/work/all.tsv

all

function annotate()
{
Rscript -e '
library(pQTLtools)
head(hg19)
> subset(hg19,grepl("^HBG|^HBB|^HBA|^HBD|^HBQ|^HBZ",SYMBOL))
        chr   start     end  width strand      ENSEMBL SYMBOL      UNIPROT
43431 chr11 5246696 5248301   1606      - ENSG0000....    HBB D9YZU5, ....
43432 chr11 5254059 5255858   1800      - ENSG0000....    HBD A0N071, ....
43434 chr11 5269502 5271087   1586      - ENSG0000....   HBG1 D9YZU8, ....
43435 chr11 5269502 5276395   6894      - ENSG0000....   HBG2 D9YZU9, ....
43436 chr11 5274421 5276395   1975      - ENSG0000....   HBG2 D9YZU9, ....
43437 chr11 5275522 5667011 391490      - ENSG0000....   HBG2 D9YZU9, ....
56896 chr16  202854  204504   1651      + ENSG0000....    HBZ       P02008
56898 chr16  222846  223709    864      + ENSG0000....   HBA2 D1MGQ2, ....
56899 chr16  226679  227520    842      + ENSG0000....   HBA1 D1MGQ2, ....
56900 chr16  230333  231178    846      + ENSG0000....   HBQ1 A0A1K0GU....
#subset(hg19Tables,grepl("^HBG|^HBA|^HBB|^HBD|^HBQ|^HBZ",hgncSym))
#      X.chrom chromStart chromEnd strand    acc uniprotName
#11496   chr11    5246830  5248251      - P68871   HBB_HUMAN
#11497   chr11    5254196  5255663      - P02042   HBD_HUMAN
#11498   chr11    5269591  5275958      - P69891  HBG1_HUMAN
#11499   chr11    5274509  5275958      - P69892  HBG2_HUMAN
#15242   chr16     202908   204396      + P02008  HBAZ_HUMAN
#15246   chr16     230485   231104      + P09105  HBAT_HUMAN
#subset(SomaScanV4.1,grepl("Hemoglobin", Human.Target.or.Analyte))
#        #     SeqID    Human.Target.or.Analyte    UniProt.ID   GeneID    Type
#2848 2848   4915-64                 Hemoglobin P69905|P68871 HBA1|HBB Protein
#2849 2849 17137-160    Hemoglobin subunit beta        P68871      HBB Protein
#2850 2850   6992-67   Hemoglobin subunit delta        P02042      HBD Protein
#2851 2851  7136-107 Hemoglobin subunit epsilon        P02100     HBE1 Protein
#2852 2852   19774-8 Hemoglobin subunit gamma-2        P69892     HBG2 Protein
#2853 2853  18198-51 Hemoglobin subunit theta-1        P09105     HBQ1 Protein
#2854 2854   7965-25 Hemoglobin subunit theta-1        P09105     HBQ1 Protein
#2855 2855    6919-3    Hemoglobin subunit zeta        P02008      HBZ Protein
'

grep -e HBG ${HbF}/work/hbf_hits.txt
grep -n HBG ${HbF}/work/*tsv | cut -f1 | xargs -l basename | sed 's/\./\t/' | cut -f1 | uniq
grep -e 4915.64 -e 17137.160 -e 6992.67 -e 7136.107 -e 19774.8 -e 18198.51 -e 7965.25 -e 6919.3 ${HbF}/work/deCODE.tsv
grep -e 4915.64 -e 17137.160 -e 6992.67 -e 7136.107 -e 19774.8 -e 18198.51 -e 7965.25 -e 6919.3 ${HbF}/work/deCODE.tsv | cut -f4 | sort | uniq | grep HB
#17137_160_HBB_Beta_globin
#18198_51_HBQ1_HBAT
#6919_3_HBZ_HBAZ
#6992_67_HBD_HBD

export deCODE=${pgwas}/deCODE
gunzip -c $deCODE/17137_160_HBB_Beta_globin.txt.gz | awk '/chr11/ && $2 >= 5246696 && $2 <= 5248301'
gunzip -c $deCODE/18198_51_HBQ1_HBAT.txt.gz | awk '/chr16/ && $2 >= 230485 && $2 <= 231104'
gunzip -c $deCODE/6919_3_HBZ_HBAZ.txt.gz | awk '/chr16/ && $2 >= 202908 && $2 <= 204396'
gunzip -c $deCODE/6992_67_HBD_HBD.txt.gz | awk '/chr11/ && $2 >= 5254196 && $2 <= 5255663'

#https://www.ebi.ac.uk/gwas/studies/GCST003122
cut -f22 ${HbF}/work/GCST003122.tsv | grep -f - ${HbF}/work/hbf_hits.txt
#6:135097778 rs7776054
grep rs ${HbF}/work/Sardinia3.txt | cut -d' ' -f1 | grep -f - -w ~/INF/work/INTERVAL.rsid
#chr2:60710571_A_G rs13019832
#chr2:60720951_A_G rs4671393
#chr6:135356216_C_G rs11754265
#chr6:135419018_C_T rs9399137
#chr11:5231565_C_T rs12793110
#chr11:5242698_C_G rs11036338
#chr11:5250168_A_G rs7936823
#chr11:5251849_G_T rs7944544
#chr11:5255582_A_C rs35152987
#chr11:5277236_C_T rs2855122
#chr11:5290370_C_G rs67385638
#chr16:149539_A_G rs570013781
#chr16:216593_C_T rs141494605
#chr16:342218_C_T rs148706947
#Trait Loci variant p-value Effect (StdErr) snpid
#HbA1
#α-globin gene cluster
#rs570013781 1.50e-18 -0.2025 (0.023) chr16:149539_A_G
#chr16:391593 3.28e-12 -0.4028 (0.058) NA
#HbA2
#β-globin gene cluster
#rs35152987 8.09e-97 -2.2450 (0.106) chr11:5255582_A_C
#rs7944544 7.90e-42 -1.3020 (0.095) chr11:5251849_G_T
#rs12793110 3.53e-18 -0.1800 (0.021) chr11:5231565_C_T
#rs11036338 5.20e-14 0.1256 (0.017) chr11:5242698_C_G
#rs7936823 5.00e-13 0.1117 (0.015) chr11:5250168_A_G
#α-globin gene cluster
#rs141494605 3.70e-43 -0.3540 (0.026) chr16:216593_C_T
#chr16:391593 1.38e-15 -0.5035 (0.063) NA
#rs148706947 1.04e-08 0.2892 (0.051) chr16:342218_C_T
#HbF
#BCL11A
#rs4671393 3.00e-76 0.4644 (0.025) chr2:60720951_A_G
#rs13019832 9.12e-33 -0.2024 (0.017) chr2:60710571_A_G
#HBS1L-MYB
#rs9399137 7.46e-94 0.5195 (0.025) chr6:135419018_C_T
#rs11754265 5.04e-12 -0.1421 (0.021) chr6:135356216_C_G
#β-globin gene cluster
#rs67385638 4.67e-35 0.3103 (0.025) chr11:5290370_C_G
#rs2855122 2.57e-11 -0.1458 (0.022) chr11:5277236_C_T

grep rs ${HbF}/work/Sardinia7.txt | sed 's/locus[1-9]//;s/(%)//;s/(g\/dl)//;s/(cond.)//;s/^[ ]*//g;s/^[1-9][ ]*//' | awk '{$1=$1};1'
}
