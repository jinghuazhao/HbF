#!/usr/bin/bash

export HbF=${HOME}/COVID-19/HbF
export pgwas=~/rds/results/public/proteomics
export M=1e6

#ARIC, https://sites.cscc.unc.edu/aric/, build 38
export ARIC=${pgwas}/ARIC
ls ${ARIC}/EA/*gz | parallel -C' ' -j15 'tabix -S1 -s1 -b2 -e2 -f {}'
cat <(awk -v OFS="\t" '{print "rsid","snpid","Gene","Somamer","Uniprot","Symbol",$0}' ${ARIC}/glm.linear.hdr) \
    <(
       sed '1d' ${HbF}/work/hbf_hits.txt | cut -f1,3,4,12,13 | grep -v -w NA | sed 's/, /;/g' | \
       while read -r chr pos rsid snpid gene
       do
          export chr=${chr}
          export pos=${pos}
          export rsid=${rsid}
          export snpid=${snpid}
          export gene=${gene}
          export region=$(awk -vchr=${chr} -vpos=${pos} -vM=${M} 'BEGIN{print chr":"pos-M"-"pos+M}')
          sed '1d' ${ARIC}/seqid.txt | cut -f1-3 | tr '\t' ' ' | \
          parallel -C' ' -j15 --env ARIC --env rsid --env snpid --env gene --env region '
             tabix ${ARIC}/EA/{1}.PHENO1.glm.linear.gz ${region} | \
             awk -v rsid=${rsid} -v snpid=${snpid} -v gene=${gene} -v somamer={1} -v uniprot={2} -v symbol={3} -v OFS="\t" "
                 \$13<=1e-5{print rsid,snpid,gene,somamer,uniprot,symbol,\$0}
             "
          '
        done
      ) > ${HbF}/work/ARIC.tsv

#AGES
export AGES=${pgwas}/AGES
cd ${AGES}
(
cat <<END
GCST90086001-GCST90087000
GCST90087001-GCST90088000
GCST90088001-GCST90089000
GCST90089001-GCST90090000
GCST90090001-GCST90091000
END
) | \
parallel -C' ' -j15 '
  lftp -c mirror https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/{}
'
for section in GCST90086001-GCST90087000 GCST90087001-GCST90088000 GCST90088001-GCST90089000 GCST90089001-GCST90090000 GCST90090001-GCST90091000
do
   mv ${section}/GC* ..
done
ls GC*/*tsv | parallel -C' ' -j15 'bgzip {}'
ls GC*/*gz | parallel -C' ' -j15 'tabix -S1 -s3 -b4 -e4 -f {}'
mdir ${AGES}/excluded
cut -f2 ${AGES}/work/AGES.txt | grep -f - -v -w <(ls $AGES/) | parallel -C' ' --env AGES 'mv ${AGES}/{} ${AGES}/excluded'
mv excluded/loc* excluded/Olink-* excluded/README.* excluded/AGES.hdr .
ls Olink-INF | parallel -C' ' 'rm -rf Olink-INF/{}; ln -fs {} Olink-INF/{}'
cd -

R --no-save <<END
  HbF <- Sys.getenv("HbF")
  stables <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-021-27850-z/MediaObjects/41467_2021_27850_MOESM18_ESM.xlsx"
  st16 <- openxlsx::read.xlsx(stables, sheet=1, colNames=TRUE, skipEmptyRows=TRUE, startRow=3)
  write.table(st16,file=file.path(HbF,"work","AGES.txt"),col.names=FALSE,quote=FALSE,row.names=FALSE,sep="\t")
END

cat <(awk -v OFS="\t" '{print "rsid","snpid","Gene","Somamer","GCST","Symbol",$0}' ${AGES}/AGES.hdr) \
    <(
       sed '1d' ${HbF}/work/hbf_hits.txt | cut -f1,3,4,12,13 | grep -v -w NA | sed 's/, /;/g' | \
       while read -r chr pos rsid snpid gene
       do
          export chr=${chr}
          export pos=${pos}
          export rsid=${rsid}
          export snpid=${snpid}
          export gene=${gene}
          export region=$(awk -vchr=${chr} -vpos=${pos} -vM=${M} 'BEGIN{print chr":"pos-M"-"pos+M}')
          awk '{print $1,$2,$7}' ${HbF}/work/AGES.txt | \
          parallel -C' ' -j15 --env AGES --env region '
             tabix ${AGES}/{2}/{2}_buildGRCh37.tsv.gz ${region} | \
             awk -v rsid=${rsid} -v snpid=${snpid} -v gene=${gene} -v somamer={1} -v gcst={2} -v symbol={3} -v OFS="\t" "
                 \$2<=1e-5{print rsid,snpid,gene,somamer,gcst,symbol,\$0}
             "
          '
       done
     ) > ${HbF}/work/AGES.tsv

#deCODE
export deCODE=${pgwas}/deCODE
cat <(awk -v OFS="\t" '{print "rsid","snpid","Gene","Somamer","GCST","Symbol",$0}' ${AGES}/AGES.hdr) \
    <(
       sed '1d' ${HbF}/work/hbf_hits.txt | cut -f1,3,4,12,13 | grep -v -w NA | sed 's/, /;/g' | \
       while read -r chr pos rsid snpid gene
       do
          export chr=${chr}
          export pos=${pos}
          export rsid=${rsid}
          export snpid=${snpid}
          export gene=${gene}
          ls ${deCODE}/*gz | sed 's/.txt.gz//' | \
          parallel -C' ' -j15 --env deCODE --env chr --env pos --env M '
             gunzip -c ${deCODE}/{}.txt.gz | \
             awk -v rsid=${rsid} -v snpid=${snpid} -v gene=${gene} -v id={} -v chr=chr${chr} -v pos=${pos} -v M=${M} -v OFS="\t" "
                 \$8<=1e-5&&\$1==chr&&\$2>=pos-M&&\$2<pos+M{print rsid,snpid,gene,id,\$0}"
          '
       done
     ) > ${HbF}/work/deCODE.tsv

#Fenland
export Fenland=${pgwas}/Fenland
cat <(awk -v OFS="\t" '{print "rsid","snpid","Gene",$0}' ${Fenland}/all.grch37.tabix.hdr) \
    <(
       sed '1d' ${HbF}/work/hbf_hits.txt | cut -f1,3,4,12,13 | grep -v -w NA | sed 's/, /;/g' | \
       while read -r chr pos rsid snpid gene
       do
          export chr=${chr}
          export pos=${pos}
          export rsid=${rsid}
          export snpid=${snpid}
          export gene=${gene}
          export region=$(awk -vchr=${chr} -vpos=${pos} -vM=${M} 'BEGIN{print chr":"pos-M"-"pos+M}')
          ls ${deCODE}/*gz | sed 's/.txt/gz/' | \
          parallel -C' ' -j15 --env Fenland --env region '
             tabix ${Fenland}/all.grch37.tabix.gz ${region} | \
             awk -v rsid=${rsid} -v snpid=${snpid} -v gene=${gene} -v OFS="\t" "
                 \$14<=1e-5{print rsid,snpid,gene,\$0}"
          '
       done
     ) > ${HbF}/work/Fenland.tsv

#SCALLOP-CVD1
export cvd1=${pgwas}/scallop-cvd1
Rscript -e '
  options(width=200)
  HbF <- Sys.getenv("HbF")
  stables <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs42255-020-00287-2/MediaObjects/42255_2020_287_MOESM3_ESM.xlsx"
  st1 <- openxlsx::read.xlsx(stables, sheet=1, colNames=TRUE, skipEmptyRows=TRUE, cols=c(1:23), rows=c(3:95))
  write.table(st1[c("UniProt.ID","Gene","Short_annotation")],file=file.path(HbF,"work","cvd1.txt"),col.names=FALSE,quote=FALSE,row.names=FALSE)
'
cat <(awk -v OFS="\t" '{print "rsid","snpid","Gene","Somamer","Symbol","Prot",$0}' ${AGES}/AGES.hdr) \
    <(
       sed '1d' ${HbF}/work/hbf_hits.txt | cut -f1,3,4,12,13 | grep -v -w NA | sed 's/, /;/g' | \
       while read -r chr pos rsid snpid gene
       do
          export chr=${chr}
          export pos=${pos}
          export rsid=${rsid}
          export snpid=${snpid}
          export gene=${gene}
          cat ${HbF}/work/cvd1.txt | \
          parallel -C' ' -j15 --env deCODE --env chr --env pos --env M '
             gunzip -c ${cvd1}/{3}.gz | \
             awk -v rsid=${rsid} -v snpid=${snpid} -v gene=${gene} -v somamer={1} -v symbol={2} -v prot={3} \
                 -v chr=${chr} -v pos=${pos} -v M=${M} -v OFS="\t" "
                 {split($1,a,":");if(\$8<=1e-5&&a[1]==chr&&a[2]>=pos-M&&a[2]<pos+M){print rsid,snpid,gene,somamer,symbol,prot,\$0}"
          '
       done
     ) > ${HbF}/work/scallop-cvd1.tsv


#LBC1936
export LBC1936=${pgwas}/LBC1936
Rscript -e '
  options(width=200)
  suppressMessages(library(dplyr))
  library(pQTLtools)
  HbF <- Sys.getenv("HbF")
  neu <- subset(Olink_qPCR,Panel=="Neurology")
  target <- neu[["Target"]]
  s <- stringr::str_locate(target,"[(]")[,1]+1
  e <- stringr::str_locate(target,"[)]")[,1]-1
  target.short <- substr(target,s,e)
  target.short <- gsub(" |-","_",target.short)
  overlap <- merge(neu,read.delim(file.path(HbF,"work","hbf_GWAS_top_snps_long.txt"))[c("snpid","gene","acc","rs.ID")],by="gene")
  overlap
'

#GTEx, eQTL Catalog
export GTEx=~/rds/public_databases/GTEx/csv
export eQTLCatalogue=~/rds/public_databases/eQTLCatalogue

gunzip -c $GTEx/Whole_Blood.tsv.gz | head -1 > ${HbF}/work/eQTL.tsv
Rscript -e '
  options(width=200)
  GTEx <- Sys.getenv("GTEx")
  eQTLCatalogue <- Sys.getenv("eQTLCatalogue")
  HbF <- Sys.getenv("HbF")
  suppressMessages(library(dplyr))
  GTEx_files <- setdiff(dir(GTEx),c("coloc","www"))
  GTEx_tissues <- unique(gsub(".tsv.gz|.tbi","",GTEx_files))
  eQTLCatalogue_files <- dir(eQTLCatalogue)
  eQTLCatalogue_tissues <- unique(gsub(".all.tsv.gz|.tbi","",eQTLCatalogue_files))
  all_files <- c(file.path(GTEx,paste0(GTEx_tissues,".tsv.gz")),file.path(eQTLCatalogue,paste0(eQTLCatalogue_tissues,".all.tsv.gz")))
  all_names <- c(GTEx_tissues,eQTLCatalogue_tissues)
  all_indices <- 1:length(all_names)
  write.table(cbind(all_indices,all_names,all_files),row.names=FALSE,col.names=FALSE,quote=FALSE)
' | \
parallel -C' ' --env HbF '
  echo {1} {2} {3}
  sed "1d" ${HbF}/work/hbf_GWAS_top_snps_long.txt | cut -f1,4,5,10 | \
  awk "{gsub(/chr|_[A-Z]*/,\"\",\$1);split(\$1,a,\":\");\$1=\$1\"-\"a[2]};1" | \
  while read -r region rsid gene ensGene
  do
     tabix {3} ${region} | awk -vensGene=${ensGene} "index(\$1,ensGene)||index(\$4,ensGene)"
  done
' >> ${HbF}/work/eQTL.tsv
Rscript -e '
  options(width=200)
  HbF <- Sys.getenv("HbF")
  suppressMessages(library(dplyr))
  eQTL <- read.delim(file.path(HbF,"work","eQTL.tsv"))
  results <- filter(eQTL,!is.na(pvalue) & pvalue<=0.05) %>%
             select(variant,pvalue,molecular_trait_object_id,beta,se,chromosome,position,type,rsid)
  write.table(results,row.names=FALSE,quote=FALSE,sep="\t")
'
awk '$3<1e-5' ${HbF}/work/eQTL.tsv | grep -v gz

#INTERVAL

export INTERVAL=~/rds/results/public/proteomics/somalogic/sun_2018/raw_results/meta

awk -vOFS="\t" 'NR>1{gsub(/chr|:[0-9]*|_[A-Z]*/,"",$1);print}' ${HbF}/work/hbf_GWAS_top_snps_long.txt | cut -f1,5 | tr '\t' ' ' | \
parallel -C' ' --env INTERVAL 'ls ${INTERVAL}/{2}.*/{2}.*chrom_{1}* | grep -v -e tbi -e info' > ${HbF}/work/INTERVAL.lst

cat <(gunzip -c ${INTERVAL}/BACH2.12756.3.3/BACH2.12756.3.3_chrom_6_meta_1.tbl.gz | head -1) \
    <(
     )
