# Deterministic integrity checks for the published CSV tables, no model fitting.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2)
src <- args[1]
handoff <- args[2]
rd <- function(n) read.csv(file.path(src, paste0(n, '.csv')), check.names = FALSE)
table <- function(n) read.csv(file.path(handoff, '03_Tables', paste0(n, '.csv')), check.names = FALSE)
for (n in c('global25', 'within15', 'within25', 'reference95', 'bulk', 'within15_symbol_sensitivity')) {
  x <- rd(paste0('KEGG_', n))
  stopifnot(all(is.finite(x$pval)), max(abs(x$padj - p.adjust(x$pval, 'BH'))) < 1e-12)
  stopifnot(!anyDuplicated(x$pathway), all(x$tested_sets == nrow(x)))
}
old <- rd('KEGG_historical_rank_high_precision')
new <- rd('KEGG_global25')
cmp <- merge(old, new, by = 'pathway')
stopifnot(nrow(cmp) == 353L,
          max(abs(cmp$NES.x - cmp$NES.y)) < 1e-12,
          max(abs(cmp$pval.x - cmp$pval.y)) < 1e-12)

cx <- table('Table06_cross_dataset_pathway_membership')
shared <- table('Table07_same_direction_intersection')
for (a in c('within15', 'global25')) {
  d <- cx[cx$analysis == a, ]
  expected <- with(d, tested_in_both & sc_FDR < 0.1 & bulk_FDR < 0.1 & sign(sc_NES) == sign(bulk_NES))
  expected[is.na(expected)] <- FALSE
  stopifnot(!anyDuplicated(d$key), sum(d$tested_in_both) == 275L,
            identical(expected, d$shared_FDR010),
            setequal(d$key[expected], shared$key[shared$analysis == a]))
  stopifnot(sum(expected & d$sc_NES < 0, na.rm = TRUE) == 0L)
  stopifnot(sum(expected & d$sc_NES > 0, na.rm = TRUE) == if (a == 'within15') 54L else 23L)
}
within <- rd('KEGG_within15')
bulk <- rd('KEGG_bulk')
stopifnot(sum(within$padj < 0.1 & within$NES > 0) == 83L,
          sum(within$padj < 0.1 & within$NES < 0) == 5L,
          sum(bulk$padj < 0.1 & bulk$NES > 0) == 128L,
          sum(bulk$padj < 0.1 & bulk$NES < 0) == 8L)
t <- within[within$pathway_id %in% c('hsa03050','hsa04137'), ]
stopifnot(nrow(t) == 2L, all(t$NES > 0), all(t$padj < 0.05))
counts <- rd('selection_cell_counts')
stopifnot(sum(counts$cells[counts$analysis == 'within15']) == 1926L,
          min(counts$cells[counts$analysis == 'within15']) == 107L)
cat('PASS: full-collection BH, historical rank reproduction, pathway intersections, directions and cell counts.\n')
