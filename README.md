# WGS-DDEC
#Load samples 
vcf_dir <- "/Users/sydney/Desktop/WGS_DDEC/"
out_dir <- "O:/Sydney/WGS_DDEC/"    
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Case IDs (must appear in the VCF filename OR match the sample ID in the VCF)
case_ids <- c(
  "GTB_120","GTB_3338","GTB_308","GTB_350","GTB_565","GTB_731","GTB_881",
  "GTB_1491","GTB_1508","GTB_1618","GTB_1761","GTB_1993","GTB_6052","GTB_6350",
  "GTB_6855","GTB_6931","GTB_6955","GTB_9025","GTB_9254","GTB_9459","GTB_9924",
  "GTB_10105","GTB_11316","GTB_11588","GTB_11915"
)

# Target genes, normalize p53 to TP53
raw_genes <- c(
  "POLE","p53","ESR1","PMS2","MSH5","MSH2","MLH1","PBRM1","ARID1A","ARID1B",
  "SMARCA4","SMARCA2","SMARCB1","SMARCD1","SMARCE1","SMARCB1"
)

#Genome build for samples (load the hg38 no alt)
use_hg38 <- TRUE 

# Install/load packages quietly
quietPkg <- function(pkg, bioc=FALSE){
  if (!requireNamespace(pkg, quietly=TRUE)) {
    if (bioc) BiocManager::install(pkg, ask=FALSE, update=FALSE)
    else install.packages(pkg, repos="https://cloud.r-project.org", quiet=TRUE)
  }
  suppressPackageStartupMessages(library(pkg, character.only=TRUE))
}

if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org", quiet=TRUE)

# Core
quietPkg("VariantAnnotation", bioc=TRUE)
quietPkg("GenomicRanges", bioc=TRUE)
quietPkg("GenomicFeatures", bioc=TRUE)
quietPkg("IRanges", bioc=TRUE)
quietPkg("S4Vectors", bioc=TRUE)
quietPkg("rtracklayer", bioc=TRUE)
quietPkg("AnnotationDbi", bioc=TRUE)
install.packages("readr")
quietPkg("readr")


# load hg38 ref genome
install.packages("seqinr")
library(seqinr)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings", force = TRUE)
library(Biostrings)
hg38 <- readDNAStringSet("/Users/sydney/Desktop/ref_genome/hg38_no_alt.fa")

# Tidying and plotting packages
quietPkg("dplyr")
quietPkg("tidyr")
quietPkg("stringr")
quietPkg("purrr")
quietPkg("readr")
quietPkg("ggplot2")
quietPkg("ComplexHeatmap", bioc=TRUE)
quietPkg("circlize", bioc=TRUE)

#normalize gene symbols for consistency
normalize_gene_symbols <- function(genes) {
  symbols <- toupper(genes)
  symbols[symbols %in% c("P53","TP-53")] <- "TP53"
  unique(symbols)
}

# Map symbols 
gene_ranges_from_symbols <- function(symbols, txdb) {
  entrez <- AnnotationDbi::select(org.Hs.eg.db, keys=symbols,
                                  keytype="SYMBOL", columns="ENTREZID") %>%
    dplyr::distinct(SYMBOL, ENTREZID) %>%
    dplyr::filter(!is.na(ENTREZID))
  if (nrow(entrez) == 0) stop("No Entrez IDs found for supplied gene symbols.")
  gr_genes <- genes(txdb)
  gr_sub <- gr_genes[gr_genes$gene_id %in% entrez$ENTREZID]
  id2sym <- entrez$SYMBOL
  names(id2sym) <- entrez$ENTREZID
  gr_sub$SYMBOL <- unname(id2sym[gr_sub$gene_id])
  grl <- split(gr_sub, gr_sub$SYMBOL)
  gr_reduced <- GenomicRanges::reduce(grl)  
  GRanges::GRangesList(gr_reduced) %>%
    unlist() %>%
    (\(gr){ gr$SYMBOL <- rep(names(gr_reduced), lengths(gr_reduced)); gr })()
}

# Pull gene symbols from VEPs
extract_genes_from_info <- function(info, field_names=c("CSQ","ANN")) {
  genes <- character(length(info))
  for (field in field_names) {
    if (field %in% names(info)) {
      raw <- info[[field]]
      v <- as.character(raw)
      first_token <- sub(",.*$", "", v)   
      parts <- strsplit(first_token, "\\|")
      guess <- vapply(parts, function(x){
        if (length(x) >= 4) x[[4]] else NA_character_
      }, character(1))
      genes <- toupper(guess)
      break
    }
  }
  genes
}

variant_key <- function(gr) {
  paste0(as.character(Seqinfo(seqnames(gr))), ":", start(gr), " ",
         mcols(gr)$REF, ">", unlist(mcols(gr)$ALT))
}

# Load one VCF to return tibble of SNVs annotated with target genes
read_case_vcf <- function(vcf_path, target_gr, target_symbols, case_id_guess=NULL) {
  vcf <- VariantAnnotation::readVcf(vcf_path)
  rr <- rowRanges(vcf)
  ref <- as.character(VariantAnnotation::ref(vcf))
  alt_list <- VariantAnnotation::alt(vcf)
  is_snv <- nchar(ref) == 1 & elementNROWS(alt_list) >= 1 & vapply(alt_list, function(a) all(nchar(as.character(a))==1), TRUE)
  rr <- rr[is_snv]
  if (length(rr) == 0) return(dplyr::tibble())
  
  # PASS filter if present
  flt <- VariantAnnotation::fixed(vcf)$FILTER
  if (!is.null(flt)) rr <- rr[which(flt[is_snv] %in% c("PASS",".", NA))]
  if (length(rr) == 0) return(dplyr::tibble())
  
  # Grab gene from INFO if possible
  info <- VariantAnnotation::info(vcf)[is_snv,,drop=FALSE]
  genes_from_info <- extract_genes_from_info(info)
  have_info <- !is.na(genes_from_info) & genes_from_info != ""
  genes_from_info <- genes_from_info[have_info]
  rr_info <- rr[have_info]
  
  tbl_info <- dplyr::tibble(
    gene = genes_from_info,
    chr = as.character(GenomeInfoDb::seqnames(rr_info)),
    pos = GenomicRanges::start(rr_info),
    REF = as.character(VariantAnnotation::ref(vcf)[is_snv][have_info]),
    ALT = vapply(VariantAnnotation::alt(vcf)[is_snv][have_info], function(a) as.character(a)[1], character(1)),
    key = paste0(chr, ":", pos, " ", REF, ">", ALT)
  ) |>
    dplyr::filter(gene %in% target_symbols)
  
  # overlap with genomic ranges if encountering problems
  rr_ov <- rr[!names(rr) %in% rr_info@ranges@NAMES]
  ov <- GenomicRanges::findOverlaps(rr_ov, target_gr, ignore.strand=TRUE)
  tbl_ov <- dplyr::tibble(
    gene = target_gr$SYMBOL[S4Vectors::subjectHits(ov)],
    chr  = as.character(GenomeInfoDb::seqnames(rr_ov))[S4Vectors::queryHits(ov)],
    pos  = GenomicRanges::start(rr_ov)[S4Vectors::queryHits(ov)],
    REF  = as.character(VariantAnnotation::ref(vcf)[is_snv][!have_info])[S4Vectors::queryHits(ov)],
    ALT  = vapply(VariantAnnotation::alt(vcf)[is_snv][!have_info][S4Vectors::queryHits(ov)], function(a) as.character(a)[1], character(1)),
    key = paste0(chr, ":", pos, " ", REF, ">", ALT)
  )
  
  # Case/sample name
  sample_id <- tryCatch({
    s <- colnames(VariantAnnotation::geno(vcf)$GT)
    if (!is.null(colnames(vcf))) colnames(vcf)[1] else if (length(s)) s[1] else NA_character_
  }, error=function(e) NA_character_)
  if (is.na(sample_id) || is.null(sample_id) || sample_id=="") {
    sample_id <- if (is.null(case_id_guess)) basename(vcf_path) else case_id_guess
  }
  
  dplyr::bind_rows(tbl_info, tbl_ov) %>%
    dplyr::mutate(case_id = sample_id) %>%
    dplyr::distinct(case_id, gene, key, .keep_all = TRUE)
}

#Done Parsing, now to the output

message("Normalizing gene list...")
genes <- normalize_gene_symbols(raw_genes)

message("Preparing gene ranges...")
target_gr <- gene_ranges_from_symbols(genes, txdb)
target_symbols <- unique(target_gr$SYMBOL)

message("Locating VCFs...")
vcf_files <- list.files(vcf_dir, pattern = "\\.vcf)?$", full.names = TRUE, ignore.case = TRUE)

# Map case IDs to files 
match_vcfs_for_case <- function(case_id){
  hits <- vcf_files[stringr::str_detect(basename(vcf_files), fixed(case_id))]
  if (length(hits)) return(hits)
  character(0)
}

message("Parsing VCFs...")
results <- purrr::map_dfr(case_ids, function(cid){
  files <- match_vcfs_for_case(cid)
  if (length(files) == 0) {
  }
  if (length(files) == 0) {
    message(sprintf("No VCF found for case %s (by filename).", cid))
    return(dplyr::tibble())
  }
  purrr::map_dfr(files, ~read_case_vcf(.x, target_gr, target_symbols, case_id_guess = cid)) |>
    dplyr::mutate(case_id = cid)  
})

# Write list
mut_list_path <- file.path(out_dir, "mutations_in_target_genes.csv")
readr::write_csv(results, mut_list_path)

# Presence/absence matrix (genes x cases)
mat <- results |>
  dplyr::distinct(case_id, gene) %>%
  tidyr::complete(case_id = case_ids, gene = genes, fill = list(present = TRUE)) %>%
  dplyr::mutate(present = TRUE) %>%
  tidyr::pivot_wider(names_from = case_id, values_from = present, values_fill = FALSE) %>%
  as.data.frame()

rownames(mat) <- mat$gene
mat$gene <- NULL
mat_bin <- as.matrix(mat) * 1 

# Save presence matrix
presence_path <- file.path(out_dir, "mutation_presence_matrix.csv")
readr::write_csv(
  tibble::rownames_to_column(as.data.frame(mat_bin), var="gene"),
  presence_path
)

# Make heatmap
if (nrow(mat_bin) > 0 && ncol(mat_bin) > 0) {
  message("Drawing heatmap...")
  col_fun <- circlize::colorRamp2(c(0,1), c("#FFFFFF","#1F77B4"))
  ht <- ComplexHeatmap::Heatmap(
    mat_bin,
    name = "SNV",
    col = col_fun,
    rect_gp = grid::gpar(col = "grey90", lwd = 0.5),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_title = "Genes",
    column_title = "Cases",
    heatmap_legend_param = list(at=c(0,1), labels=c("Absent","Present"))
  )
  heatmap_file <- file.path(out_dir, "SNV_heatmap.png")
  grDevices::png(heatmap_file, width = 1600, height = 1000, res = 150)
  grid::grid.newpage(); draw(ht)
  grDevices::dev.off()
  message(sprintf("Heatmap saved to %s", heatmap_file))
} else {
  message("No data for heatmap.")
}

message("Done.")
message(sprintf("Tidy SNV list: %s", mut_list_path))
message(sprintf("Presence matrix: %s", presence_path))
