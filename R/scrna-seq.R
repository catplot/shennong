#' Run Seurat Workflow
#'
#' Executes a Seurat workflow for single-cell RNA-seq data analysis. Supports
#' both standard and SCTransform workflows, followed by PCA reduction and
#' post-processing.
#'
#' @inheritParams sn_run_seurat_sctransform
#' @inheritParams sn_run_seurat_postprocessing
#' @inheritParams sn_standardize_genes
#' @param workflow The workflow to use. Can be either "standard" or
#'   "sctransform".
#' @param batch Attribute for batch.
#' @param integration_method Integration method. Can be either "harmony", "mnn",
#'   "cca", "rpca", or "scvi".
#' @param cell_cycle_scoring Whether to score cell cycle phases.
#' @param ... Arguments passed to other methods.
#' @return An updated Seurat object after performing the specified workflow and
#'   post-processing.
#' @importFrom rlang arg_match
#' @importFrom Seurat RunPCA CellCycleScoring NormalizeData RPCAIntegration
#'   HarmonyIntegration CCAIntegration IntegrateLayers
#' @export
#' @concept scrna-seq
sn_run_seurat <- function(object,
                          assay = NULL,
                          workflow = "sctransform",
                          nfeatures = 3000,
                          vars_to_regress = NULL,
                          vst_flavor = "v2",
                          reduction = "pca",
                          batch = NULL,
                          integration_method = "harmony",
                          dims = 1:30,
                          resolution = 0.8,
                          seed = 717,
                          verbose = TRUE,
                          run_tsne = FALSE,
                          cell_cycle_scoring = TRUE,
                          ...) {
  arg_match(workflow, values = c("standard", "sctransform"))
  arg_match(integration_method,
    values = c(
      "harmony", "mnn", "cca", "rpca", "scvi"
    )
  )

  assay <- assay %||% DefaultAssay(object)

  if (cell_cycle_scoring) {
    object <- NormalizeData(
      object = object,
      verbose = verbose
    )
    object <- CellCycleScoring(
      object = object,
      s.features = Seurat::cc.genes.updated.2019$s.genes,
      g2m.features = Seurat::cc.genes.updated.2019$g2m.genes,
      set.ident = FALSE
    )
  }

  if (!is.null(x = batch)) {
    object[[assay]] <- split(
      x = object[[assay]],
      f = object@meta.data[, batch]
    )

    if (integration_method %in% c("mnn", "scvi")) {
      sn_check_package("SeuratWrappers")
    }
    selected_integration_method <- switch(
      EXPR = integration_method,
      harmony = HarmonyIntegration,
      mnn = SeuratWrappers::FastMNNIntegration,
      cca = CCAIntegration,
      rpca = RPCAIntegration,
      scvi = SeuratWrappers::scVIIntegration
    )
  }

  if (workflow == "standard") {
    object <- sn_run_seurat_standard(
      object = object,
      nfeatures = nfeatures,
      vars_to_regress = vars_to_regress,
      verbose = verbose
    )
  } else {
    object <- sn_run_seurat_sctransform(
      object = object,
      nfeatures = nfeatures,
      vars_to_regress = vars_to_regress,
      vst_flavor = vst_flavor,
      seed = seed,
      verbose = verbose
    )
    assay <- "SCT"
  }

  # Perform linear dimensional reduction
  object <- RunPCA(
    object = object,
    seed.use = seed,
    verbose = verbose
  )

  if (!is.null(x = batch)) {
    orig_reducation <- reduction
    reduction <- paste0("integrated.", integration_method)
    object <- IntegrateLayers(
      object = object,
      method = selected_integration_method,
      orig.reduction = orig_reducation,
      assay = assay,
      new.reduction = reduction,
      verbose = verbose,
      ...
    )
  }

  # Perform Seurat post processing
  object <- sn_run_seurat_postprocessing(
    object = object,
    reduction = reduction,
    dims = dims,
    resolution = resolution,
    seed = seed,
    verbose = verbose,
    run_tsne = run_tsne
  )

  return(object)
}

#' Standard Seurat Workflow
#'
#' Performs the standard workflow for Seurat single-cell RNA-seq data analysis.
#' This includes normalization, identification of variable features, and
#' scaling.
#'
#' @param object An object of class Seurat, representing single-cell RNA-seq
#'   data.
#' @param nfeatures Use this many features as variable features after ranking by
#'   residual variance; default is 3000. Only applied if residual.features is
#'   not set.
#' @param vars_to_regress Variables to regress out in a second non-regularized
#'   linear regression. For example, percent.mito. Default is NULL.
#' @param verbose Whether to print messages and progress bars.
#'
#' @return Seurat object after applying the standard workflow.
#' @importFrom Seurat NormalizeData FindVariableFeatures ScaleData
#' @export
#' @concept scrna-seq
sn_run_seurat_standard <- function(object,
                                   nfeatures = 2000,
                                   vars_to_regress = NULL,
                                   verbose = TRUE) {
  object <- NormalizeData(object = object, verbose = verbose)
  object <- FindVariableFeatures(
    object = object,
    nfeatures = nfeatures,
    verbose = verbose
  )
  object <- ScaleData(
    object = object,
    vars.to.regress = vars_to_regress,
    verbose = verbose
  )
  return(object)
}

#' SCTransform Workflow for Seurat
#'
#' Applies SCTransform workflow for single-cell RNA-seq data analysis. This
#' includes normalization and variance stabilization. Requires `glmGamPoi`
#' package for improved speed.
#'
#' @inheritParams sn_run_seurat_standard
#' @param vst_flavor When set to 'v2' sets method = glmGamPoi_offset,
#'   n_cells=2000, and exclude_poisson = TRUE which causes the model to learn
#'   theta and intercept only besides excluding poisson genes from learning and
#'   regularization.
#' @param seed Set a random seed. By default, sets the seed to 1448145. Setting
#'   NULL will not set a seed.
#'
#' @return Returns a Seurat object with a new assay (named SCT by default) with
#'   counts being (corrected) counts, data being log1p(counts), scale.data being
#'   pearson residuals; sctransform::vst intermediate results are saved in misc
#'   slot of the new assay.
#' @importFrom Seurat SCTransform
#' @export
#' @concept scrna-seq
sn_run_seurat_sctransform <- function(object,
                                      nfeatures = 3000,
                                      vars_to_regress = NULL,
                                      vst_flavor = "v2",
                                      seed = 717,
                                      verbose = TRUE) {
  # Check if glmGamPoi package is installed
  sn_check_package("glmGamPoi", error = FALSE)

  object <- SCTransform(
    object = object,
    variable.features.n = nfeatures,
    vars.to.regress = vars_to_regress,
    vst.flavor = vst_flavor,
    seed.use = seed,
    verbose = verbose
  )
  return(object)
}

#' Post-Processing for Seurat Workflow
#'
#' Executes post-processing steps for Seurat workflow including neighbor
#' finding, clustering, and UMAP/t-SNE reduction.
#'
#' @inheritParams sn_run_seurat_sctransform
#' @param reduction Reduction to use as input for building the (S)NN.
#' @param dims Dimensions of reduction to use as input.
#' @param resolution Value of the resolution parameter, use a value above
#'   (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#' @param run_tsne Whether to run t-SNE reduction.
#'
#' @return Seurat object after post-processing.
#' @importFrom Seurat FindNeighbors FindClusters RunUMAP RunTSNE
#' @export
#' @concept scrna-seq
sn_run_seurat_postprocessing <- function(
    object,
    reduction = "pca",
    dims = 1:30,
    resolution = 0.8,
    seed = 717,
    verbose = TRUE,
    run_tsne = FALSE) {
  object <- FindNeighbors(
    object = object,
    reduction = reduction,
    dims = dims,
    verbose = verbose
  )
  object <- FindClusters(
    object = object,
    resolution = resolution,
    random.seed = seed,
    verbose = verbose
  )
  object <- RunUMAP(
    object = object,
    dims = dims,
    reduction = reduction,
    seed.use = seed,
    verbose = verbose
  )
  if (run_tsne) {
    object <- RunTSNE(
      object = object,
      dims = dims,
      reduction = reduction,
      seed.use = seed
    )
  }
  return(object)
}

#' Standardize gene symbols
#'
#' Converts gene names of a Seurat single-cell object to a dictionary of
#' standard symbols. This function is useful prior to integration of datasets
#' from different studies, where gene names may be inconsistent.
#'
#' @param object Seurat object
#' @param assay Name of assay to fetch layer data from or assign layer data to.
#' @param layer Name of layer to fetch or set.
#' @param ensembl_gene_table A data frame of gene name mappings. This should
#'   have the format of Ensembl BioMart tables with fields "Gene name", "Gene
#'   Synonym" and "Gene stable ID" (and optionally "NCBI gene (formerly
#'   Entrezgene) ID").
#'
#' @return Returns a Seurat object or a vector of gene names with standard gene
#'   names. Genes not found in the standard list are removed. Synonyms are
#'   accepted when the conversion is not ambiguous.
#' @importFrom data.table fread
#' @importFrom SeuratObject DefaultAssay CreateSeuratObject GetAssayData
#' @importFrom utils head
#' @export
#' @concept scrna-seq
sn_standardize_genes <- function(
    object,
    assay = NULL,
    layer = NULL,
    ensembl_gene_table = NULL) {
  # If file is given
  if (!is.null(ensembl_gene_table)) {
    ensembl_gene_table <- fread(ensembl_gene_table)
  } else {
    stop("Please provide a gene name mapping table!")
  }

  if (inherits(x = object, what = "Seurat")) {
    assay <- assay %||% DefaultAssay(object = object)
    layer <- layer %||% "counts"

    # Translate Ensembl IDs if necessary
    genes_in <- rownames(
      GetAssayData(
        object = object,
        assay = assay,
        layer = layer
      )
    )
  } else if (is.character(object)) {
    genes_in <- object
  } else {
    stop("Please provide a Seurat object or a vector of gene names!")
  }

  ngenes <- length(x = genes_in)

  ens_count <- length(intersect(
    x = genes_in,
    y = ensembl_gene_table[["Gene stable ID"]]
  ))
  gname_count <- length(intersect(
    x = genes_in,
    y = ensembl_gene_table[["Gene name"]]
  ))

  ncbi_count <- 0
  if ("NCBI gene (formerly Entrezgene) ID" %in% colnames(ensembl_gene_table)) {
    ncbi_count <- length(
      intersect(
        x = genes_in,
        y = ensembl_gene_table[["NCBI gene (formerly Entrezgene) ID"]]
      )
    )
  }

  max <- max(ens_count, gname_count, ncbi_count)
  if (max < length(genes_in) / 2) {
    warning(
      "Over 50% of genes in input object not found in reference gene table"
    )
  }

  gname_format <- FALSE
  if (max == gname_count) {
    gname_format <- TRUE
  }

  if (max == ens_count) {
    # Input object has Ensembl IDs
    to <- "Gene name"
    from <- "Gene stable ID"
    genes_tr <- ensembl_gene_table[[to]][match(
      x = genes_in, table = ensembl_gene_table[[from]]
    )]
    names(x = genes_tr) <- genes_in
    genes_tr <- genes_tr[!is.na(genes_tr) & genes_tr != ""]
  } else if (max == ncbi_count) {
    to <- "Gene name"
    from <- "NCBI gene (formerly Entrezgene) ID"
    genes_tr <- ensembl_gene_table[[to]][match(
      x = genes_in, table = ensembl_gene_table[[from]]
    )]
    names(x = genes_tr) <- genes_in
    genes_tr <- genes_tr[!is.na(genes_tr) & genes_tr != ""]
  } else {
    genes_tr <- genes_in
    names(genes_tr) <- genes_in
  }

  # 1. First match dictionary
  gene_ref_dict <- ensembl_gene_table[["Gene name"]]
  names(x = gene_ref_dict) <- ensembl_gene_table[["Gene Synonym"]]
  gene_ref_dict <- gene_ref_dict[!is.null(x = names(x = gene_ref_dict))]

  message(paste("Number of genes in input object:", ngenes))
  # Keep genes with standard Gene name
  genes_allow_list1 <- genes_tr[
    !is.na(x = genes_tr) &
      genes_tr != "" &
      genes_tr %in% ensembl_gene_table[["Gene name"]]
  ]
  l <- length(x = genes_allow_list1)

  message(
    sprintf(
      "Number of genes with standard symbols: %i (%.2f%%)", l, l / ngenes * 100
    )
  )

  if (l < ngenes && gname_format) {
    message(paste("Examples of non-standard Gene names:"))
    ns <- head(x = genes_tr[!genes_tr %in% ensembl_gene_table[["Gene name"]]])
    message(paste(unname(ns), collapse = ","))
  }

  # 2. Search among synonyms
  # Keep genes with accepted Gene name synonym
  genes_allow_list2 <- genes_tr[
    !genes_tr %in% ensembl_gene_table[["Gene name"]] &
      genes_tr %in% ensembl_gene_table[["Gene Synonym"]]
  ]
  # Translate Gene Synonym to standard Gene name
  genes_allow_list2_gn <- gene_ref_dict[genes_allow_list2]

  message(
    paste(
      "Additional number of genes with accepted Gene name synonym: ",
      length(genes_allow_list2_gn)
    )
  )

  # Names of genes_allow_list contain IDs in matrix
  # elements contain the new names
  genes_allow_list <- c(genes_allow_list1, genes_allow_list2_gn)

  # 3. Check for duplicates
  is_dup <- duplicated(genes_allow_list)
  genes_allow_list <- genes_allow_list[!is_dup]
  message(
    sprintf(
      "Number of duplicated Gene name: %i (%.2f%%)",
      sum(is_dup), sum(is_dup) / ngenes * 100
    )
  )

  l <- length(genes_allow_list)
  message(sprintf("Final number of genes: %i (%.2f%%)", l, l / ngenes * 100))

  if (inherits(x = object, what = "Seurat")) {
    # 4. Subset matrix for allowed genes, and translate names
    matrix <- list()
    matrix[[layer]] <- GetAssayData(
      object = object,
      assay = assay,
      layer = layer
    )
    rows_select <- rownames(x = matrix[[layer]])[
      rownames(x = matrix[[layer]]) %in% names(x = genes_allow_list)
    ]
    matrix[[layer]] <- matrix[[layer]][rows_select, ]
    rownames(matrix[[layer]]) <- unname(obj = genes_allow_list[rows_select])

    counts <- matrix[["counts"]]
    metadata <- object@meta.data
    object <- CreateSeuratObject(
      counts = counts, assay = assay, meta.data = metadata
    )
  } else {
    rows_select <- genes_in[genes_in %in% names(x = genes_allow_list)]
    object <- unname(obj = genes_allow_list[rows_select])
  }

  return(object)
}
