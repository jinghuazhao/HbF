#!/usr/bin/bash

export HbF=${HOME}/COVID-19/HbF
export pgwas=~/rds/results/public/proteomics

# ARIC, https://sites.cscc.unc.edu/aric/

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

# AGES

export stables=https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-021-27850-z/MediaObjects/41467_2021_27850_MOESM18_ESM.xlsx

R --no-save <<END
  INF <- Sys.getenv("INF")
  HbF <- Sys.getenv("HbF")
  pgwas <- Sys.getenv("pgwas")
  stables <- Sys.getenv("stables")
  st16 <- openxlsx::read.xlsx(stables, sheet=1, colNames=TRUE, skipEmptyRows=TRUE, startRow=3)
  write.table(st16,file=file.path(HbF,"work","AGES.txt"),col.names=FALSE,quote=FALSE,row.names=FALSE,sep="\t")
  suppressMessages(library(dplyr))
  url1 <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90086001-GCST90087000"
  url2 <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90087001-GCST90088000"
  url3 <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90088001-GCST90089000"
  url4 <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90089001-GCST90090000"
  url5 <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90090001-GCST90091000"
  urls <- mutate(st16, tag=gsub("_",".",Study.tag),
                 trait=gsub("Serum levels of protein ","",Reported.trait),
                 Accession=Study.Accession,
                 SOMAMER_ID=paste0(trait,".",tag)) %>%
          left_join(pQTLtools::SomaLogic160410) %>%
          select(Accession,Summary.statistics.file,SOMAMER_ID,UniProt,Target) %>%
          merge(read.delim(file.path(HbF,"work","hbf_GWAS_top_snps_long.txt")),by.x="UniProt",by.y="acc") %>%
          filter(!is.na(UniProt)) %>%
          distinct() %>%
          mutate(url=case_when(
                                 Accession >= "GCST90086001" & Accession <= "GCST90087000" ~ url1,
                                 Accession >= "GCST90087001" & Accession <= "GCST90088000" ~ url2,
                                 Accession >= "GCST90088001" & Accession <= "GCST90089000" ~ url3,
                                 Accession >= "GCST90089001" & Accession <= "GCST90090000" ~ url4,
                                 Accession >= "GCST90090001" & Accession <= "GCST90091000" ~ url5,
                                 TRUE ~ as.character(Accession)
                              ),
                 src=file.path(url,Accession),cmd=paste("lftp -c mirror",src))
  write.table(select(urls,Accession,rs.ID,gene),file=file.path(HbF,"work","AGES.sh"),quote=FALSE,row.names=FALSE)
  write.table(select(urls,Accession,SOMAMER_ID,UniProt,Target) %>% left_join(select(pQTLtools::inf1,-target),by=c('UniProt'='uniprot')),
              file=file.path(HbF,"work","links.txt"),quote=FALSE,row.names=FALSE,sep="\t")
  write.table(select(urls,cmd),file=file.path(HbF,"work","lftp.sh"),quote=FALSE,col.names=FALSE,row.names=FALSE)
END

cd ${HbF}/work
bash lftp.sh
sed '1d' hbf_GWAS_top_snps_long.txt | cut -f5 | sort | uniq | grep -f - -w AGES.txt | cut -f2
lftp -c mirror https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90086001-GCST90087000/GCST90086395
lftp -c mirror https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90086001-GCST90087000/GCST90086825
lfpt -c mirror https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90087001-GCST90088000/GCST90087207
bgzip GC*/*tsv
cat ${pgwas}/AGES/AGES.hdr
sed '1d' ${HbF}/work/AGES.sh | \
parallel -C' ' '
  zgrep -w {2} {1}/{1}_buildGRCh37.tsv.gz
'
zgrep -w -e rs7772031 -e rs7776054 -e rs9376095 ${HbF}/work/GCST90086825/GCST90086825_buildGRCh37.tsv.gz
cd -
#variant_id      p_value chromosome      base_pair_location      effect_allele   other_allele    effect_allele_frequency beta    standard_error
#rs4707609       0.0717480331291094      6       90946479        C       T       0.352   0.0340157       0.0188864
#rs10128556      0.919811849322134       11      5263683 T       C       0.3493  -0.00174019     0.0172851
#rs7772031       0.713143617980774       6       135259734       G       A       0.3357  0.00722091      0.0196404
#rs7776054       0.212966846431685       6       135418916       G       A       0.279   0.0257332       0.0206594
#rs9376095       0.536243449393508       6       135450755       C       T       0.2279  -0.0136041      0.0219938

# deCODE

# Fenland

# GTEx, eQTL Catalog

# SCALLOP-CVD1

# ~/COVID-19/HbF/work/hbf_GWAS_top_snps_long.txt
snpid   Locus   snpid38 rs.ID   gene    p       acc     geneSynonyms    hgncSym ensGene prot
chr2:53993037_T_TTCTTTGGTGCTCCATTGATAAGAGCACCC  1       chr2:53765900_T_TTCTTTGGTGCTCCATTGATAAGAGCACCC  rs148978228     CHAC2   3.93e-15        Q8WUX2     >
chr2:60718043_G_T       2       chr2:60490908_G_T       rs1427407       BCL11A  3.01e-161       Q9H165  CTIP1; EVI9; KIAA1809; ZNF856   BCL11A  ENSG0000011>

# ${pgwas}/ARIC/seqid.txt
seqid_in_sample uniprot_id      entrezgenesymbol        chromosome_name transcription_start_site
SeqId_10000_28  P43320  CRYBB2  22      25212564

# ${pgwas}/AGES/AGES.hdr
# variant_id	p_value	chromosome	base_pair_location	effect_allele	other_allele	effect_allele_frequency	beta	standard_error

# deCODE
