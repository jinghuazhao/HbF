options(width=2000)
require(openxlsx)
suppressMessages(require(dplyr))
suppressMessages(require(gap))
require(tidyr)
hbf_dir <- "/rds/project/jmmh2/rds-jmmh2-projects/olink_proteomics/scallop/HbF"
hbf_xlsx <- "hbf_GWAS_top_snps.xlsx"
if(!dir.exists(file.path(hbf_dir,"work"))) dir.create(file.path(hbf_dir,"work"))
hbf_hits <- read.xlsx(file.path(hbf_dir,hbf_xlsx),sheet=1,colNames=TRUE,startRow=1,skipEmptyRows=TRUE) %>%
            rename(gene="Candidate.gene(s)",b=X9,p.value="p-value") %>%
            mutate(p=pvalue(b/SE))
# only as a reference, chr2:60718043_G_T rs1427407 indicates that it is in GRCh37/hg19 format
hbf_hits
hbf_hits_long <- separate_rows(hbf_hits, gene, sep=",", convert = TRUE) %>%
                 mutate(gene=gsub(" ","",gene),snpid=chr_pos_a1_a2(CHR,BP,REF,ALT)) %>%
                 left_join(data.frame(pQTLtools::hg19Tables),by=c('gene'='geneName')) %>%
                 select(snpid,names(hbf_hits),acc,X.chrom,chromStart,chromEnd,uniprotName,geneSynonyms,hgncSym,ensGene) %>%
                 mutate(prot=gsub("_HUMAN","",uniprotName)) %>%
                 select(-p.value,-uniprotName)
# a more useable form
data.frame(hbf_hits_long)
write.table(hbf_hits,file=file.path(hbf_dir,"work","hbf_GWAS_top_snps.txt"),row.names=FALSE,quote=FALSE,sep="\t")
write.table(hbf_hits_long,file=file.path(hbf_dir,"work","hbf_GWAS_top_snps_long.txt"),row.names=FALSE,quote=FALSE,sep="\t")
