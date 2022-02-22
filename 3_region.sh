#!/usr/bin/bash

export HbF=${HOME}/COVID-19/HbF
export pgwas=~/rds/results/public/proteomics
export M=0
export p_gwas=1

#1. AGES
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
rm -rf ${AGES}/excluded
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
       sed '1d' ${HbF}/work/hbf_hits.txt | cut -f1,2,4,11,13 | grep -v -w NA | sed 's/, /;/g' | \
       while read -r chr pos rsid snpid gene
       do
          export chr=${chr}
          export pos=${pos}
          export rsid=${rsid}
          export snpid=${snpid}
          export gene=${gene}
          export region=$(awk -vchr=${chr} -vpos=${pos} -vM=${M} 'BEGIN{print chr":"pos-M"-"pos+M}')
          awk '{print $1,$2,$7}' ${HbF}/work/AGES.txt | \
          parallel -C' ' -j15 --env AGES --env region --env p_gwas '
             tabix ${AGES}/{2}/{2}_buildGRCh37.tsv.gz ${region} | \
             awk -v rsid=${rsid} -v snpid=${snpid} -v gene=${gene} -v somamer={1} -v gcst={2} -v symbol={3} -v p=${p_gwas} -v OFS="\t" "
                 \$2<=p{print rsid,snpid,gene,somamer,gcst,symbol,\$0}
             "
          '
       done
     ) > ${HbF}/work/AGES.tsv

#2. ARIC, https://sites.cscc.unc.edu/aric/, GRCh38
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
          parallel -C' ' -j15 --env ARIC --env rsid --env snpid --env gene --env region --env p_gwas '
             tabix ${ARIC}/EA/{1}.PHENO1.glm.linear.gz ${region} | \
             awk -v rsid=${rsid} -v snpid=${snpid} -v gene=${gene} -v somamer={1} -v uniprot={2} -v symbol={3} -v p=${p_gwas} -v OFS="\t" "
                 \$13<=p{print rsid,snpid,gene,somamer,uniprot,symbol,\$0}
             "
          '
        done
      ) > ${HbF}/work/ARIC.tsv

#3. deCODE, https://download.decode.is, GRCh38
export deCODE=${pgwas}/deCODE
cat <(cut -f7 --complement ${deCODE}/doc/deCODE.hdr | awk -v OFS="\t" '{print "rsid","snpid","Gene","id",$0,"EAF"}') \
    <(
       while read chr pos rsid snpid gene < <(sed '1d' ${HbF}/work/hbf_hits.txt | cut -f1,3,4,11,13 | sed 's/, /;/g')
       do
          export chr=${chr}
          export pos=${pos}
          export rsid=${rsid}
          export snpid=${snpid}
          export gene=${gene}
          export region=$(awk -vchr=${chr} -vpos=${pos} -vM=${M} 'BEGIN{print "chr"chr":"pos-M"-"pos+M}')
          ls ${deCODE}/*gz | xargs -l basename -s .txt.gz | \
          parallel -C' ' --env deCODE --env p_gwas --env region '
             tabix ${deCODE}/{}.txt.gz ${region} | \
             awk -v rsid=${rsid} -v snpid=${snpid} -v gene=${gene} -v id={} -v p=${p_gwas} -v OFS="\t" "\$8<=p{print rsid,snpid,gene,id,\$0}" | \
             sort -k7,7 | \
             join -17 - <(tabix $deCODE/doc/bgzip/assocvariants.annotated.txt.gz ${region} | sort -k3,3 | cut -f3,7) | \
             cut -f1 --complement | \
             sort -k1,1 -k2,2
          '
       done
     ) > ${HbF}/work/deCODE.tsv

#4. Fenland
export Fenland=${pgwas}/Fenland
cat <(awk -v OFS="\t" '{print "rsid","snpid","Gene",$0}' ${Fenland}/all.grch37.tabix.hdr) \
    <(
       sed '1d' ${HbF}/work/hbf_hits.txt | cut -f1,2,4,11,13 | sed 's/, /;/g' | \
       while read -r chr pos rsid snpid gene
       do
          export chr=${chr}
          export pos=${pos}
          export rsid=${rsid}
          export snpid=${snpid}
          export gene=${gene}
          export region=$(awk -vchr=${chr} -vpos=${pos} -vM=${M} 'BEGIN{print chr":"pos-M"-"pos+M}')
          tabix ${Fenland}/all.grch37.tabix.gz ${region} | \
          awk -v rsid=${rsid} -v snpid=${snpid} -v gene=${gene} -v p=${p_gwas} -v OFS="\t" "\$14<=p{print rsid,snpid,gene,\$0}"
       done
     ) > ${HbF}/work/Fenland.tsv

#SCALLOP-CVD1
export scallop_cvd1=${pgwas}/scallop-cvd1
Rscript -e '
  options(width=200)
  HbF <- Sys.getenv("HbF")
  stables <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs42255-020-00287-2/MediaObjects/42255_2020_287_MOESM3_ESM.xlsx"
  st1 <- openxlsx::read.xlsx(stables, sheet=1, colNames=TRUE, skipEmptyRows=TRUE, cols=c(1:23), rows=c(3:95))
  write.table(st1[c("UniProt.ID","Gene","Short_annotation")],file=file.path(HbF,"work","scallop-cvd1.txt"),col.names=FALSE,quote=FALSE,row.names=FALSE)
'
cat <(awk -v OFS="\t" '{print "rsid","snpid","Gene","UniProt","Symbol","Prot",$0}' ${scallop_cvd1}/cvd1.hdr) \
    <(
       sed '1d' ${HbF}/work/hbf_hits.txt | cut -f1,2,4,11,13 | sed 's/, /;/g' | \
       while read -r chr pos rsid snpid gene
       do
          export chr=${chr}
          export pos=${pos}
          export rsid=${rsid}
          export snpid=${snpid}
          export gene=${gene}
          grep -v -w -e BNP -e IL4 ${HbF}/work/scallop-cvd1.txt | \
          parallel -C' ' -j15 --env deCODE --env chr --env pos --env M --env p_gwas '
             gunzip -c ${scallop_cvd1}/{3}.txt.gz | \
             awk -v rsid=${rsid} -v snpid=${snpid} -v gene=${gene} -v somamer={1} -v symbol={2} -v prot={3} \
                 -v chr=${chr} -v pos=${pos} -v M=${M} -v p=${p_gwas} -v OFS="\t" "
                 {split(\$1,a,\":\");if(\$8<=p&&a[1]==chr&&a[2]>=pos-M&&a[2]<=pos+M){print rsid,snpid,gene,somamer,symbol,prot,\$0}}"
          '
       done
     ) > ${HbF}/work/scallop-cvd1.tsv

#5. LBC1936
export LBC1936=${pgwas}/LBC1936
Rscript -e '
  options(width=200)
  suppressMessages(library(dplyr))
  library(pQTLtools)
  HbF <- Sys.getenv("HbF")
  neu <- subset(Olink_qPCR,Panel=="Neurology") %>%
         mutate(target=Target,s=stringr::str_locate(target,"[(]")[,1]+1,e=stringr::str_locate(target,"[)]")[,1]-1,
         target.short=substr(target,s,e),target.short=gsub(" |-","_",target.short)) %>%
         select(UniProt,gene,target.short)
  write.table(neu,file=file.path(HbF,"work","LBC1936.txt"),col.names=FALSE,row.names=FALSE,quote=FALSE)
'
cat <(awk -v OFS="\t" '{print "rsid","snpid","Gene","UniProt","Symbol","Prot","Ln",$0}' ${LBC1936}/LBC1936.hdr | sed 's/\"//g;s/,/\t/g') \
    <(
       sed '1d' ${HbF}/work/hbf_hits.txt | cut -f1,2,4,11,13 | sed 's/, /;/g' | \
       while read -r chr pos rsid snpid gene
       do
          export chr=${chr}
          export pos=${pos}
          export rsid=${rsid}
          export snpid=${snpid}
          export gene=${gene}
          sed 's/GCP5/GPC5/;s/N_Cdase/N_CDase/;s/, /;/;s/IL_12B,_IL_12A/IL12/' ${HbF}/work/LBC1936.txt | grep -v BDNF | \
          parallel -C' ' -j15 --env LBC1936 --env chr --env pos --env M --env p_gwas '
             gunzip -c ${LBC1936}/{3}.txt.gz | sed "s/\"//g;s/,/\t/g" | \
             awk -v rsid=${rsid} -v snpid=${snpid} -v gene=${gene} -v uniprot={1} -v symbol={2} -v prot={3} \
                 -v chr=${chr} -v pos=${pos} -v M=${M} -v p=${p_gwas} -v OFS="\t" "
                 \$9<=p&&\$3==chr&&\$4>=pos-M&&\$4<=pos+M{print rsid,snpid,gene,uniprot,symbol,prot,\$0}" | \
             grep -v -w -e CPM -e rs138154139
          '
       done
     ) > ${HbF}/work/LBC1936.tsv

#6. INTERVAL
export INTERVAL=~/rds/results/public/proteomics/somalogic/sun_2018/raw_results/meta
cat <(gunzip -c ${INTERVAL}/BACH2.12756.3.3/BACH2.12756.3.3_chrom_6_meta_1.tbl.gz | head -1 |
      awk -v OFS="\t" '{print "rsid","snpid","Gene","id",$0}') \
    <(
       sed '1d' ${HbF}/work/hbf_hits.txt | cut -f1,2,4,11,13 | sed 's/, /;/g' | \
       while read -r chr pos rsid snpid gene
       do
          export chr=${chr}
          export pos=${pos}
          export rsid=${rsid}
          export snpid=${snpid}
          export gene=${gene}
          export region=$(awk -vchr=${chr} -vpos=${pos} -vM=${M} 'BEGIN{print chr":"pos-M"-"pos+M}')
          ls ${INTERVAL} | xargs -l basename | sed 's/_meta_1.tbl.gz//' | \
          parallel -C' ' -j15 --env INTERVAL --env chr --env pos --env M --env p_gwas '
             tabix ${INTERVAL}/{}/{}_chrom_${chr}_meta_1.tbl.gz ${region} | \
             awk -v rsid=${rsid} -v snpid=${snpid} -v gene=${gene} -v id={} -v p=${p_gwas} -v OFS="\t" "10^\$8<=p{print rsid,snpid,gene,id,\$0}"
          '
       done
     ) > ${HbF}/work/INTERVAL.tsv

#7. GTEx, eQTL Catalog, GRCh38
export GTEx=~/rds/public_databases/GTEx/csv
export eQTLCatalogue=~/rds/public_databases/eQTLCatalogue
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
  write.table(cbind(all_indices,all_names,all_files),file=file.path(HbF,"work","eQTL.txt"),row.names=FALSE,col.names=FALSE,quote=FALSE)
'
cat <(gunzip -c ${GTEx}/Adipose_Subcutaneous.tsv.gz | head -1 | awk -v OFS="\t" '{print "rsid","snpid","Gene","id",$0}') \
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
          awk 'NR<50' ${HbF}/work/eQTL.txt | \
          parallel -C' ' -j15 --env rsid --env snpid --env gene --env region --env p_gwas '
             tabix {3} ${region} | \
             awk -v rsid=${rsid} -v snpid=${snpid} -v gene=${gene} -v id={2} -v p=${p_gwas} -v OFS="\t" "\$3<=p{print rsid,snpid,gene,id,\$0}"
          '
        done
      ) > ${HbF}/work/GTEx.tsv
cat <(gunzip -c ${eQTLCatalogue}/Alasoo_2018_ge_macrophage_IFNg.all.tsv.gz | head -1 | awk -v OFS="\t" '{print "rsid","snpid","Gene","id",$0}') \
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
          awk 'NR>49' ${HbF}/work/eQTL.txt | \
          parallel -C' ' -j15 --env rsid --env snpid --env gene --env region --env p_gwas '
             tabix {3} ${region} | \
             awk -v rsid=${rsid} -v snpid=${snpid} -v gene=${gene} -v id={2} -v p=${p_gwas} -v OFS="\t" "\$9<=p{print rsid,snpid,gene,id,\$0}"
          '
        done
      ) > ${HbF}/work/eQTL.tsv

#8. eQTLGen
export eQTLGen=~/rds/public_databases/eQTLGen/tabix
cat <(gunzip -c ${eQTLGen}/cis_full.txt.gz | head -1 |
      awk -v OFS="\t" '{print "rsid","snpid","Gene","id",$0}') \
    <(
       sed '1d' ${HbF}/work/hbf_hits.txt | cut -f1,2,4,11,13 | sed 's/, /;/g' | \
       while read -r chr pos rsid snpid gene
       do
          export chr=${chr}
          export pos=${pos}
          export rsid=${rsid}
          export snpid=${snpid}
          export gene=${gene}
          export region=$(awk -vchr=${chr} -vpos=${pos} -vM=${M} 'BEGIN{print chr":"pos-M"-"pos+M}')
          parallel -C' ' -j15 --env INTERVAL --env chr --env pos --env M --env p_gwas '
             tabix ${eQTLGen}/{}.txt.gz ${region} | \
             awk -v rsid=${rsid} -v snpid=${snpid} -v gene=${gene} -v id={} -v p=${p_gwas} -v OFS="\t" "\$1<=p{print rsid,snpid,gene,id,\$0}"
          ' ::: cis_full
       done
     ) > ${HbF}/work/eQTLGen-cis_full.tsv
cat <(gunzip -c ${eQTLGen}/trans.txt.gz | head -1 |
      awk -v OFS="\t" '{print "rsid","snpid","Gene","id",$0}') \
    <(
       sed '1d' ${HbF}/work/hbf_hits.txt | cut -f1,2,4,11,13 | sed 's/, /;/g' | \
       while read -r chr pos rsid snpid gene
       do
          export chr=${chr}
          export pos=${pos}
          export rsid=${rsid}
          export snpid=${snpid}
          export gene=${gene}
          export region=$(awk -vchr=${chr} -vpos=${pos} -vM=${M} 'BEGIN{print chr":"pos-M"-"pos+M}')
          parallel -C' ' -j15 --env INTERVAL --env chr --env pos --env M --env p_gwas '
             tabix ${eQTLGen}/{}.txt.gz ${region} | \
             awk -v rsid=${rsid} -v snpid=${snpid} -v gene=${gene} -v id={} -v p=${p_gwas} -v OFS="\t" "\$1<=p{print rsid,snpid,gene,id,\$0}"
          ' ::: trans
       done
     ) > ${HbF}/work/eQTLGen-trans.tsv

Rscript -e '
  options(width=200)
  HbF <- Sys.getenv("HbF")
  suppressMessages(library(dplyr))
  GTEx <- read.delim(file.path(HbF,"work","GTEx.tsv")) %>%
          filter(!is.na(pvalue)) %>%
          select(variant,pvalue,molecular_trait_object_id,beta,se,chromosome,position,type,rsid)
  eQTL <- read.delim(file.path(HbF,"work","eQTL.tsv")) %>%
          filter(!is.na(pvalue)) %>%
          select(variant,pvalue,molecular_trait_object_id,beta,se,chromosome,position,type,rsid)
  results <- bind_rows(GTEx,eQTL)
  write.table(results,row.names=FALSE,quote=FALSE,sep="\t")
'
