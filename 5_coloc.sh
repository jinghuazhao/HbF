#!/usr/bin/bash

#SBATCH --job-name=_coloc
#SBATCH --account CARDIO-SL0-CPU
#SBATCH --partition cardio
#SBATCH --qos=cardio
#SBATCH --array=1-4
#SBATCH --mem=28800
#SBATCH --time=5-00:00:00
#SBATCH --error=/rds/user/jhz22/hpc-work/work/_coloc_%A_%a.err
#SBATCH --output=/rds/user/jhz22/hpc-work/work/_coloc_%A_%a.out
#SBATCH --export ALL

if [ ! -d ${HbF}/coloc ]; then mkdir ${HbF}/coloc; fi
if [ ! -d ${HbF}/eQTLCatalogue ]; then mkdir ${HbF}/eQTLCatalogue; fi
read chr pos pos38 rsid < \
                        <(sed '1d' ${HbF}/work/hbf_hits.txt | cut -f1-4 | awk 'NR==ENVIRON["SLURM_ARRAY_TASK_ID"]')
if [ ! -f ${HbF}/$1/${prot}-${rsid}.pdf ] || [ ! -f ${HbF}/$1/${prot}-${rsid}.rds ]; then
   cd ${HbF}/$1
   Rscript -e '
      liftRegion <- function(x,chain,flanking=1e6)
      {
        require(GenomicRanges)
        gr <- with(x,GenomicRanges::GRanges(seqnames=chr,IRanges::IRanges(start,end))+flanking)
        seqlevelsStyle(gr) <- "UCSC"
        gr38 <- rtracklayer::liftOver(gr, chain)
        chr <- gsub("chr","",colnames(table(seqnames(gr38))))
        start <- min(unlist(start(gr38)))
        end <- max(unlist(end(gr38)))
        invisible(list(chr=chr,start=start,end=end,region=paste0(chr,":",start,"-",end)))
      }
      sumstats <- function(prot,chr,region37)
      {
        cat("GWAS sumstats\n")
        vcf <- file.path(HbF,"METAL/gwas2vcf",paste0(prot,".vcf.gz"))
        gwas_stats <- gwasvcf::query_gwas(vcf, chrompos = region37)
        gwas_stats <- gwasvcf::vcf_to_granges(gwas_stats) %>%
                      keepSeqlevels(chr) %>%
                      renameSeqlevels(paste0("chr",chr))
        gwas_stats_hg38 <- rtracklayer::liftOver(gwas_stats, chain) %>%
        unlist() %>%
        renameSeqlevels(chr) %>%
        dplyr::as_tibble() %>%
        dplyr::transmute(chromosome = seqnames,
                         position = start, REF, ALT, AF, ES, SE, LP, SS) %>%
        dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>%
        dplyr::mutate(MAF = pmin(AF, 1-AF)) %>%
        dplyr::group_by(id) %>%
        dplyr::mutate(row_count = n()) %>%
        dplyr::ungroup() %>%
        dplyr::filter(row_count == 1)
      }
      gtex <- function(gwas_stats_hg38,ensGene,region38)
      {
        cat("c. GTEx_v8 imported eQTL datasets\n")
        f <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","tabix_ftp_paths_gtex.tsv")
        imported_tabix_paths <- within(read.delim(f, stringsAsFactors = FALSE) %>% dplyr::as_tibble(),
              {ftp_path <- gsub("ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv/GTEx_V8/ge",
                                paste0(HOME,"/rds/public_databases/GTEx/csv"),ftp_path)})
        rnaseq_df <- dplyr::filter(imported_tabix_paths, quant_method == "ge") %>%
                     dplyr::mutate(qtl_id = paste(study, qtl_group, sep = "_"))
        ftp_path_list <- setNames(as.list(rnaseq_df$ftp_path), rnaseq_df$qtl_id)
        hdr <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","column_names.GTEx")
        column_names <- names(read.delim(hdr))
        safe_import <- purrr::safely(import_eQTLCatalogue)
        summary_list <- purrr::map(ftp_path_list,
                                   ~safe_import(., region38, selected_gene_id = ensGene, column_names))
        result_list <- purrr::map(summary_list, ~.$result)
        result_list <- result_list[!unlist(purrr::map(result_list, is.null))]
        result_filtered <- purrr::map(result_list[lapply(result_list,nrow)!=0],
                                      ~dplyr::filter(., !is.na(se)))
        purrr::map_df(result_filtered, ~run_coloc(., gwas_stats_hg38), .id = "qtl_id")
      }
      ge <- function(gwas_stats_hg38,ensGene,region38)
      {
        cat("d. eQTL datasets\n")
        f <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","tabix_ftp_paths_ge.tsv")
        imported_tabix_paths <- within(read.delim(f, stringsAsFactors = FALSE) %>% dplyr::as_tibble(),
              {ftp_path <- gsub("ftp://ftp.ebi.ac.uk/pub/databases/spot/eQTL/csv/Alasoo_2018/ge",
                                paste0(HOME,"/rds/public_databases/eQTLCatalogue"),ftp_path)})
        ftp_path_list <- setNames(as.list(imported_tabix_paths$ftp_path), imported_tabix_paths$unique_id)
        hdr <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","column_names.Alasoo")
        column_names <- names(read.delim(hdr))
        safe_import <- purrr::safely(import_eQTLCatalogue)
        summary_list <- purrr::map(ftp_path_list,
                                   ~safe_import(., region38, selected_gene_id = ensGene, column_names))
        result_list <- purrr::map(summary_list, ~.$result)
        result_list <- result_list[!unlist(purrr::map(result_list, is.null))]
        result_filtered <- purrr::map(result_list[lapply(result_list,nrow)!=0],
                                      ~dplyr::filter(., !is.na(se)))
        purrr::map_df(result_filtered, ~run_coloc(., gwas_stats_hg38), .id = "unique_id")
      }
      gtex_coloc <- function(prot,chr,ensGene,chain,region37,region38,out)
      {
        gwas_stats_hg38 <- sumstats(prot,chr,region37)
        df_gtex <- gtex(gwas_stats_hg38,ensGene,region38)
        if (exists("df_gtex"))
        {
          saveRDS(df_gtex,file=paste0(out,".rds"))
        # dplyr::arrange(df_gtex, -PP.H4.abf)
        # p <- ggplot(df_gtex, aes(x = PP.H4.abf)) + geom_histogram()
        }
        s <- ggplot(gwas_stats_hg38, aes(x = position, y = LP)) + geom_point()
        ggsave(plot = s, filename = paste0(out, "-assoc.pdf"), path = "", device = "pdf",
               height = 15, width = 15, units = "cm", dpi = 300)
      # ggsave(plot = p, filename = paste0(out, "-hist.pdf"), path = "", device = "pdf",
      #        height = 15, width = 15, units = "cm", dpi = 300)
      }
      ge_coloc <- function(prot,chr,ensGene,chain,region37,region38,out)
      {
        gwas_stats_hg38 <- sumstats(prot,chr,region37)
        df_ge <- ge(gwas_stats_hg38,ensGene,region38)
        if (exists("df_ge"))
        {
           saveRDS(df_ge,file=paste0(out,".rds"))
        #  dplyr::arrange(df_ge, -PP.H4.abf)
        #  p <- ggplot(df_ge, aes(x = PP.H4.abf)) + geom_histogram()
        }
        s <- ggplot(gwas_stats_hg38, aes(x = position, y = LP)) + geom_point()
        ggsave(plot = s, filename = paste0(out, "-assoc.pdf"), path = "", device = "pdf",
               height = 15, width = 15, units = "cm", dpi = 300)
      # ggsave(plot = p, filename = paste0(out, "-hist.pdf"), path = "", device = "pdf",
      #        height = 15, width = 15, units = "cm", dpi = 300)
      }
      single_run <- function(r, batch="GTEx")
      {
        sentinel <- sentinels[r,]
        isnpid <- within(gap::inv_chr_pos_a1_a2(sentinel[["SNP"]]),
        {
          chr <- gsub("chr","",chr)
          pos <- as.integer(pos)
          start <- pos-M
          if (start<0) start <- 0
          end <- pos+M
        })
        chr <- with(isnpid,chr)
        region37 <- with(isnpid, paste0(chr,":",start,"-",end))
        ensRegion37 <- with(subset(inf1,prot==sentinel[["prot"]]),
                            {
                              start <- start-M
                              if (start<0) start <- 0
                              end <- end+M
                              paste0(chr,":",start,"-",end)
                            })
        region38 <- with(liftRegion(isnpid,chain),region)
        ensGene <- subset(inf1,prot==sentinel[["prot"]])[["ensembl_gene_id"]]
        ensRegion38 <- with(liftRegion(subset(inf1,prot==sentinel[["prot"]]),chain),region)
        cat(chr,region37,region38,ensGene,ensRegion37,ensRegion38,"\n")
        if (batch=="GTEx")
        {
          f <- file.path(HbF,"coloc",with(sentinel,paste0(prot,"-",SNP)))
          gtex_coloc(sentinel[["prot"]],chr,ensGene,chain,region37,region38,f)
        # gtex_coloc(sentinel[["prot"]],chr,ensGene,chain,ensRegion37,ensRegion38,f)
        } else {
          f <- file.path(HbF,"eQTLCatalogue",with(sentinel,paste0(prot,"-",SNP)))
          ge_coloc(sentinel[["prot"]],chr,ensGene,chain,region37,region38,f)
        # ge_coloc(sentinel[["prot"]],chr,ensGene,chain,ensRegion37,ensRegion38,f)
        }
      }
      collect <- function(coloc_dir="eQTLCatalogue")
      # to collect results when all single runs are done
      {
        df_coloc <- data.frame()
        for(r in 1:nrow(sentinels))
        {
          prot <- sentinels[["prot"]][r]
          snpid <- sentinels[["SNP"]][r]
          rsid <- prot_rsid[["SNP"]][r]
          f <- file.path(HbF,coloc_dir,paste0(prot,"-",snpid,".rds"))
          if (!file.exists(f)) next
          cat(prot,"-",rsid,"\n")
          rds <- readRDS(f)
          if (nrow(rds)==0) next
          df_coloc <- rbind(df_coloc,data.frame(prot=prot,rsid=rsid,snpid=snpid,rds))
        }
        df <- dplyr::rename(df_coloc,H0=PP.H0.abf,H1=PP.H1.abf,H2=PP.H2.abf,H3=PP.H3.abf,H4=PP.H4.abf)
        if (coloc_dir=="coloc") {
           df_coloc <- within(df,{qtl_id <- gsub("GTEx_V8_","",qtl_id)})
           write.table(subset(df,H3+H4>=0.9&H4/H3>=3),file=file.path(HbF,coloc_dir,"GTEx.tsv"),
                       quote=FALSE,row.names=FALSE,sep="\t")
           write.table(df,file=file.path(HbF,coloc_dir,"GTEx-all.tsv"),
                       quote=FALSE,row.names=FALSE,sep="\t")
        } else {
          write.table(subset(df,H3+H4>=0.9&H4/H3>=3),file=file.path(HbF,coloc_dir,"eQTLCatalogue.tsv"),
                      quote=FALSE,row.names=FALSE,sep="\t")
          write.table(df,file=file.path(HbF,coloc_dir,"eQTLCatalogue-all.tsv"),
                      quote=FALSE,row.names=FALSE,sep="\t")
        }
      }
      loop_slowly <- function() for (r in 1:nrow(sentinels)) single_run(r)
    # Environmental variables
      options(width=200)
      HOME <- Sys.getenv("HOME")
      HPC_WORK <- Sys.getenv("HPC_WORK")
      HbF <- Sys.getenv("HbF")
      M <- 5e5
      pkgs <- c("dplyr", "gap", "ggplot2", "readr", "coloc", "GenomicRanges","pQTLtools","seqminer")
      invisible(suppressMessages(lapply(pkgs, require, character.only = TRUE)))
      f <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","hg19ToHg38.over.chain")
      chain <- rtracklayer::import.chain(f)
      gwasvcf::set_bcftools(file.path(HPC_WORK,"bin","bcftools"))
      f <- file.path(path.package("pQTLtools"),"eQTL-Catalogue","tabix_ftp_paths.tsv")
      tabix_paths <- read.delim(f, stringsAsFactors = FALSE) %>% dplyr::as_tibble()
      sentinels <- read.delim(file.path(HbF,"work","hbf_hits.txt")) %>%
                   rename(chr=CHR,pos=BP,rsid=rs.ID)
      r <- as.integer(Sys.getenv("r"))
      single_run(r)
      single_run(r,batch="eQTLCatalogue")
  ' 2>&1 | \
  tee ${prot}-${rsid}.log
# ls *tbi | xargs -I {} bash -c "rm {}"
  cd -
fi

function coloc_gene()
{
Rscript -e '
    suppressMessages(library(dplyr))
    library(openxlsx)
    HbF <- Sys.getenv("HbF")
    metal <- read.delim(file.path(HbF,"work","HbF1.METAL"))
    pqtls <- left_join(metal, pQTLtools::inf1)
    xlsx <- file.path(HbF,"NG","SCALLOP-HbF-ST.xlsx")
    plist <- function(sheet)
             read.xlsx(xlsx,sheet=sheet,startRow=2) %>%
             filter(flag=="x") %>%
             select(Protein) %>%
             distinct()
    gtex <- plist("ST4-GTEx coloc")
    eqtlgen <- plist("EQTLGen_coloc")
    eqtlcat <- plist("EQTL-Catalogue_coloc")
    pall <- unique(bind_rows(gtex,eqtlgen,eqtlcat)) %>%
            rename(target.short=Protein) %>%
            left_join(pQTLtools::inf1) %>%
            select(target.short,uniprot,prot,gene,ensembl_gene_id)
    genes <- pull(pall,gene)
    write.table(genes,row.names=FALSE,col.names=FALSE,quote=FALSE)
'
}
