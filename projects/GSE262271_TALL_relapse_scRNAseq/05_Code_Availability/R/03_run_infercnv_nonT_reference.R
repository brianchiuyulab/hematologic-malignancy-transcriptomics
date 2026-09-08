options(stringsAsFactors = FALSE)

user_lib <- Sys.getenv("R_LIBS_USER", unset = "C:/Users/User/Documents/R/win-library/4.4")
if (dir.exists(user_lib)) .libPaths(unique(c(normalizePath(user_lib, winslash = "/"), .libPaths())))

configure_jags <- function() {
  candidates <- unique(c(
    Sys.getenv("JAGS_HOME", unset = NA_character_),
    "C:/Users/User/Documents/JAGS/JAGS-4.3.2_extracted",
    "C:/Program Files/JAGS/JAGS-4.3.2",
    "C:/Program Files (x86)/JAGS/JAGS-4.3.2"
  ))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates) & dir.exists(candidates)]
  if (!length(candidates)) return(invisible(FALSE))
  jags_home <- normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE)
  Sys.setenv(JAGS_HOME = jags_home)
  jags_paths <- file.path(jags_home, c("x64/bin", "bin", "x64/modules", "modules"))
  jags_paths <- jags_paths[dir.exists(jags_paths)]
  if (length(jags_paths)) {
    current_path <- strsplit(Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1]]
    Sys.setenv(PATH = paste(unique(c(jags_paths, current_path)), collapse = .Platform$path.sep))
  }
  invisible(TRUE)
}
configure_jags()

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = "") {
  hit <- match(flag, args)
  if (!is.na(hit) && length(args) >= hit + 1) args[[hit + 1]] else default
}

default_input <- file.path(
  "C:/Users/User/Desktop/Public dataset analysis/GEO dataset/2026/single cell",
  "BTZ 相關",
  "GSE262271_TALL_relapse_scRNAseq",
  "infercnv_inputs"
)
input_dir <- normalizePath(get_arg("--input_dir", default_input), winslash = "/", mustWork = TRUE)

if (!requireNamespace("infercnv", quietly = TRUE)) {
  stop("Package infercnv is not installed. Install infercnv, then rerun this script.", call. = FALSE)
}

counts <- readRDS(file.path(input_dir, "infercnv_raw_counts_matrix.rds"))
annot_file <- file.path(input_dir, "infercnv_cell_annotations.txt")
gene_order_file <- file.path(input_dir, "infercnv_gene_order_from_orgHsEgDb.txt")
annotations <- read.delim(annot_file, header = FALSE, stringsAsFactors = FALSE)
colnames(annotations) <- c("cell", "group")

annotations <- annotations[annotations$group != "ref_Cytotoxic_T_or_NK", , drop = FALSE]
counts <- counts[, annotations$cell, drop = FALSE]
non_t_annot_file <- file.path(input_dir, "infercnv_cell_annotations_nonT_reference.txt")
utils::write.table(annotations, non_t_annot_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

ref_groups <- unique(annotations$group[grepl("^ref_", annotations$group)])
out_dir <- file.path(input_dir, "infercnv_output_nonT_reference")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

infercnv_obj <- infercnv::CreateInfercnvObject(
  raw_counts_matrix = counts,
  annotations_file = non_t_annot_file,
  delim = "\t",
  gene_order_file = gene_order_file,
  ref_group_names = ref_groups
)
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = out_dir,
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = FALSE,
  num_threads = max(1, parallel::detectCores() - 1)
)
saveRDS(infercnv_obj, file.path(out_dir, "infercnv_obj.rds"))
