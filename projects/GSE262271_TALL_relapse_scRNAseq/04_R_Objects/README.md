# R object availability

Large Seurat and inferCNV objects are intentionally excluded from GitHub. The reviewed final figures can be regenerated from the compact frozen tables where possible; analyses that require cell-level embeddings or inferCNV burden need the corresponding local objects.

See `05_Code_Availability/README.md` for the reproduction boundary and required R packages.

The local project retains the original GSE262271_seurat_qc_harmony_annotated.rds unchanged.
The current GSE262271_reviewed_lineage_annotation.rds adds reviewed_cell_type,
reviewed_infercnv_score, potential_malignant_global25 and
potential_malignant_within15_exploratory. Its active identities use broad reviewed
lineage labels, while original annotation columns remain available for provenance.
This annotation review is not independent malignant-cell validation.

reviewed/reviewed_cache.rds contains the audited GSEA results,
normalized expression and selections. reviewed/reviewed_metadata_and_selections.rds
contains the smaller annotation and cell-selection record. Neither is uploaded.
