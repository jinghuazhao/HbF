#!/usr/bin/bash

export HbF=${HOME}/COVID-19/HbF
export pgwas=~/rds/results/public/proteomics

#ARIC, https://sites.cscc.unc.edu/aric/

cat <(awk -v OFS="\t" '{print "Gene",$0}' ${pgwas}/ARIC/glm.linear.hdr) \
    <(
       Rscript -e '
       options(width=200)
       suppressMessages(library(dplyr))
       HbF <- Sys.getenv("HbF")
       pgwas <- Sys.getenv("pgwas")
       snps <- read.delim(file.path(HbF,"work","hbf_GWAS_top_snps_long.txt")) %>%
               left_join(read.delim(file.path(pgwas,"ARIC","seqid.txt")),by=c("acc"="uniprot_id")) %>%
               select(seqid_in_sample,acc,entrezgenesymbol,snpid,snpid38,rs.ID,gene) %>%
               subset(!is.na(seqid_in_sample))
       write.table(select(snps,seqid_in_sample,rs.ID,gene),col.names=FALSE,row.names=FALSE,quote=FALSE)
       ' | \
       parallel -C' ' --env pgwas '
         zgrep -w {2} ${pgwas}/ARIC/EA/{1}.PHENO1.glm.linear.gz | \
         awk -v gene={3} -v OFS="\t" "{print gene,\$0}"
       '
     )

#Gene    #CHROM  POS     ID      REF     ALT     A1      A1_FREQ TEST    OBS_CT  BETA    SE      T_STAT  P       ERRCODE
#BACH2   6       90236760        rs4707609       T       C       C       0.351865        ADD     7213    0.0144076       0.0174087       0.827607        0.40792     .
#MYB     6       134938596       rs7772031       G       A       G       0.36961 ADD     7213    0.0103618       0.0172359       0.601175        0.547742   .
#MYB     6       135097778       rs7776054       A       G       G       0.269721        ADD     7213    0.0292828       0.0187838       1.55894 0.119054   .
#MYB     6       135129617       rs9376095       T       C       C       0.23118 ADD     7213    0.0152785       0.0197551       0.773394        0.439315   .

#AGES

R --no-save <<END
  HbF <- Sys.getenv("HbF")
  stables <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-021-27850-z/MediaObjects/41467_2021_27850_MOESM18_ESM.xlsx"
  st16 <- openxlsx::read.xlsx(stables, sheet=1, colNames=TRUE, skipEmptyRows=TRUE, startRow=3)
  write.table(st16,file=file.path(HbF,"work","AGES.txt"),col.names=FALSE,quote=FALSE,row.names=FALSE,sep="\t")
END

cd ${HbF}/work
sed '1d' hbf_GWAS_top_snps_long.txt | cut -f5 | sort | uniq | grep -f - -w AGES.txt | cut -f2
lftp -c mirror https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90086001-GCST90087000/GCST90086395
lftp -c mirror https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90086001-GCST90087000/GCST90086825
lfpt -c mirror https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90087001-GCST90088000/GCST90087207
bgzip GC*/*tsv
cat ${pgwas}/AGES/AGES.hdr
(
cat <<END
Accession rs.ID gene
GCST90086395 rs10128556 HBG1
GCST90086395 rs398114988 HBG1
GCST90087207 rs4707609 BACH2
GCST90086825 rs7772031 MYB
GCST90086825 rs7776054 MYB
GCST90086825 rs9376095 MYB
END
) | \
sed '1d' | \
parallel -C' ' '
  echo {1}-{2}-{3}
  zgrep -w {2} {1}/{1}_buildGRCh37.tsv.gz
'
cd -
#variant_id      p_value chromosome      base_pair_location      effect_allele   other_allele    effect_allele_frequency beta    standard_error
#GCST90086395-rs10128556-HBG1
#rs10128556      0.919811849322134       11      5263683 T       C       0.3493  -0.00174019     0.0172851
#GCST90086395-rs398114988-HBG1
#GCST90086825-rs7772031-MYB
#rs7772031       0.713143617980774       6       135259734       G       A       0.3357  0.00722091      0.0196404
#GCST90086825-rs7776054-MYB
#rs7776054       0.212966846431685       6       135418916       G       A       0.279   0.0257332       0.0206594
#GCST90086825-rs9376095-MYB
#rs9376095       0.536243449393508       6       135450755       C       T       0.2279  -0.0136041      0.0219938
#GCST90087207-rs4707609-BACH2
#rs4707609       0.0717480331291094      6       90946479        C       T       0.352   0.0340157       0.0188864

#deCODE

export deCODE=${pgwas}/deCODE
ls ${deCODE} | tr '_' '\t' | cut -f3 | sort | uniq | grep -f <(sed '1d' hbf_GWAS_top_snps_long.txt | cut -f5 | sort | uniq) -w -
ls ${deCODE} | grep -e MYB_MYB -e BACH2_BACH2 | grep gz
zgrep -w rs4707609 ${deCODE}/11618_83_MYB_MYB.txt.gz
zgrep -w -e rs7772031 -e rs7776054 -e rs9376095 ${deCODE}/12756_3_BACH2_BACH2.txt.gz
# Chrom   Pos     Name    rsids   effectAllele    otherAllele     Beta    Pval    minus_log10_pval        SE      N       ImpMAF
#chr6    90236760        chr6:90236760:C:T       rs4707609       C       T       0.0172  0.037333        1.42791 0.008261        35335   0.34794
#chr6    134938596       chr6:134938596:A:G      rs7772031       A       G       -0.0087 0.296918        0.52736 0.008341        35363   0.32928
#chr6    135097778       chr6:135097778:G:A      rs7776054       G       A       0.0241  0.006332        2.19846 0.008828        35363   0.27214
#chr6    135129617       chr6:135129617:C:T      rs9376095       C       T       0.0020  0.833272        0.07921 0.009501        35363   0.22820

#Fenland

export Fenland=${pgwas}/Fenland
cat $Fenland/all.grch37.tabix.hdr > ${HbF}/work/Fenland.tsv
Rscript -e '
  options(width=200)
  HbF <- Sys.getenv("HbF")
  suppressMessages(library(dplyr))
  ids <- right_join(pQTLtools::SomaScanV4.1,read.delim(file.path(HbF,"work","hbf_GWAS_top_snps_long.txt")),by=c("GeneID"="gene","UniProt.ID"="acc")) %>%
         mutate(SeqID=gsub("-","_",SeqID),
                chrpos=gsub("chr|_[A-Z]*","",snpid),
                region=paste0(chrpos,"-",unlist(lapply(strsplit(chrpos,":"),"[",2)))) %>%
         select(region,SeqID,UniProt.ID,GeneID,snpid,rs.ID,hgncSym,ensGene) %>%
         filter(!is.na(SeqID))
  write.table(ids,col.names=FALSE,row.names=FALSE,quote=FALSE)
' | \
parallel -C' ' '
  echo {1}-{2}-{3}-{4}-{5}-{6}
  tabix ${Fenland}/all.grch37.tabix.gz {1} | \
  grep -w {2}
' >> ${HbF}/work/Fenland.tsv

#chr	pos	rsid	MarkerName	Allele1	Allele2	Freq1	FreqSE	MinFreq	MaxFreq	Somamer	Effect	StdErr	Pvalue	Direction	HetISq	HetChiSq	HetDf	HetPVal	TotalSampleSize
#2:60718043-60718043-24486_1-Q9H165-BCL11A-chr2:60718043_G_T-rs1427407
#6:1034474346-1034474346-20544_103-Q13002-GRIK2-chr6:1034474346_C_T-rs190118557
#2:60725451-60725451-24486_1-Q9H165-BCL11A-chr2:60725451_C_G-rs7606173
#2:53993037-53993037-21156_5-Q8WUX2-CHAC2-chr2:53993037_T_TTCTTTGGTGCTCCATTGATAAGAGCACCC-rs148978228
#6:90946479-90946479-12756_3-Q9BYV9-BACH2-chr6:90946479_C_T-rs4707609
#11:5263683-5263683-19774_8-P69892-HBG2-chr11:5263683_C_T-rs10128556
#11:5293232-5293232-19774_8-P69892-HBG2-chr11:5293232_T_TA-rs398114988
#6:135418916-135418916-11618_83-P10242-MYB-chr6:135418916_A_G-rs7776054
#6:135450755-135450755-11618_83-P10242-MYB-chr6:135450755_C_T-rs9376095
#6	135450755	rs9376095	chr6:135450755_C_T	t	c	0.7700	0.0035	0.7681	0.7805	11618_83	0.0375	0.0162	0.02076	+++	0.0	0.058	2	0.9715	10708
#6:135259734-135259734-11618_83-P10242-MYB-chr6:135259734_A_G-rs7772031


#SCALLOP-CVD1

Rscript -e '
  options(width=200)
  HbF <- Sys.getenv("HbF")
  stables <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs42255-020-00287-2/MediaObjects/42255_2020_287_MOESM3_ESM.xlsx"
  st1 <- openxlsx::read.xlsx(stables, sheet=1, colNames=TRUE, skipEmptyRows=TRUE, cols=c(1:23), rows=c(3:95))
  write.table(st1[c("UniProt.ID","Gene","Short_annotation","Long.name")],file=file.path(HbF,"work","cvd1.txt"),
              col.names=FALSE,quote=FALSE,row.names=FALSE,sep="\t")
  suppressMessages(library(dplyr))
  ids <- merge(st1[1:4],read.delim(file.path(HbF,"work","hbf_GWAS_top_snps_long.txt"))[c("snpid","gene","acc","rs.ID")],
               by.x=c("Gene"),by.y=c("gene")) %>%
         select(Short_annotation,UniProt.ID,Gene,snpid,rs.ID)
  write.table(ids,col.names=FALSE,row.names=FALSE,quote=FALSE)
  st1[[2]] %in% read.delim(file.path(HbF,"work","hbf_GWAS_top_snps_long.txt"))[["gene"]]
'

#LBC1936

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

#variant pvalue  molecular_trait_object_id       beta    se      chromosome      position        type    rsid
#Brain_Hypothalamus
#chr6_90946479_T_TA      0.025515        ENSG00000112182 -0.167097       0.0739524       6       90946479        INDEL   rs76787307
#Brain_Substantia_nigra
#chr6_90946479_T_TA      0.0266003       ENSG00000112182 0.187447        0.0831331       6       90946479        INDEL   rs76787307
#Ovary
#chr6_90946479_T_TA      0.034979        ENSG00000112182 -0.1029 0.0482772       6       90946479        INDEL   rs76787307
#Testis
#chr6_90946479_T_TA      0.045364        ENSG00000112182 -0.0962325      0.0478606       6       90946479        INDEL   rs76787307

#INTERVAL

export INTERVAL=~/rds/results/public/proteomics/somalogic/sun_2018/raw_results/meta

awk -vOFS="\t" 'NR>1{gsub(/chr|:[0-9]*|_[A-Z]*/,"",$1);print}' ${HbF}/work/hbf_GWAS_top_snps_long.txt | cut -f1,5 | tr '\t' ' ' | \
parallel -C' ' --env INTERVAL 'ls ${INTERVAL}/{2}.*/{2}.*chrom_{1}* | grep -v -e tbi -e info' > ${HbF}/work/INTERVAL.lst

cat <(gunzip -c ${INTERVAL}/BACH2.12756.3.3/BACH2.12756.3.3_chrom_6_meta_1.tbl.gz | head -1) \
    <(
      tabix ${INTERVAL}/BACH2.12756.3.3/BACH2.12756.3.3_chrom_6_meta_1.tbl.gz 6:90236760-90236760
      tabix ${INTERVAL}/HBG1.10665.30.3/HBG1.10665.30.3_chrom_11_meta_1.tbl.gz 11:5242453-5242453
      tabix ${INTERVAL}/HBG1.10665.30.3/HBG1.10665.30.3_chrom_11_meta_1.tbl.gz 11:5272002-5272002
     )

grep '#' 2_lookup.sh | sed 's/#//'
