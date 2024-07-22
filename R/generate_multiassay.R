################################################################################
#' Generate a MultiAssayExperiment object from multiple single-omics data
#'
#' Generate a MultiAssayExperiment object from single-omics data matrix
#' and metadata. Currently supports transcriptomics and proteomics data only.
#'
#' @param rawdata_rna A data-frame containing raw RNA count where rows and
#' columns correspond to genes and samples respectively
#' @param rawdata_protein A data-frame containing raw protein abundance where
#' rows and columns correspond to genes and samples respectively
#' @param individual_to_sample Logical whether individual ID and sample names in
#' raw data matrices are the same. If they are different, `mpa_rna` and
#'  `map_protein`
#' dataframes should be provided. Default to `FALSE`.
#' @param map_rna A data-frame of two columns named `primary` and `colname`
#' where primary should contain unique ID name with a link to individual
#' metadata and colname is the column names of the `rawdata_rna` data-frame
#' @param map_protein A data-frame of two columns named `primary` and `colname`
#' where primary should contain unique ID name with a link to individual
#' metadata and colname is the column names of the `rawdata_protein` data-frame
#' @param metadata_rna A data-frame containing `rna` assay specific metadata
#' where rownames are same as the colnames of `rawdata_rna` data.frame
#' @param metadata_protein A data-frame containing `protein` assay specific
#' metadata where rownames are same as the colnames of `rawdata_protein`
#' data.frame
#' @param individual_metadata A data-frame containing sample level metadata
#' @param map_by_column The common column name to link `metadata_rna`
#' and `metadata_protein` to the `individual_metadata`
#' @param rna_qc_data Logical whether to add rna QC data.
#' @param rna_qc_data_matrix A data.frame containing RNA qc level data where
#' rownames are same as the colnames of `rawdata_rna` data.frame.
#' @param organism The organism type for transcripts and proteins.
#' Possible values are `human` and `mouse`. Default to `human`.
#'
#' @return a MultiAssayExperiment object
#'
#' @family Pre-processing
#'
#' @importFrom MultiAssayExperiment MultiAssayExperiment listToMap colData
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom magrittr set_names
#' @importFrom cli cli_alert_danger style_bold cli_alert_success
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @importFrom EnsDb.Hsapiens.v86 EnsDb.Hsapiens.v86
#' @export
generate_multiassay <- function(assay_list, # names = rna_raw, protein_raw, IDENTIFIER [generic_raw]
                                individual_metadata,
                                map_list = NULL, # order & names must match the assay list
                                sample_metadata_list = NULL,
                                map_by_column = NULL, ## make optional as above; defaults to rownames
                                qc_data_list = NULL, # names & order need to match those in assay_list
                                organism = "human") {


  # check for >= 2 assays provided
  list_length <- length(assay_list)
  if (list_length < 2) {
    stop(cli::cli_alert_danger(
      paste("At least 2 assays should be provided in",
            cli::style_bold("assay_list"),
            sep = " "
      )
    ))
  }

  ## get assay types
  assay_names <- names(assay_list)
  assay_types <- gsub(".*\\[(.*)\\].*", "\\1", assay_names)

  ## check assay types match allowed
  if (!all(assay_types %in% c("rna_raw", "protein_raw", "generic_raw"))) {
    stop(cli::cli_alert_danger(
      "Assays should be of class rna_raw, protein_raw, or generic_raw"
      ))
  }

  ## generate maps if not provided
  if (is.null(map_list)) {
    map_list <- lapply(assay_list, function(x) {
      data.frame(
        primary = colnames(x),
        colname = colnames(x),
        stringsAsFactors = FALSE
      ) -> df
      rownames(df) <- colnames(x)
      return(df)
    })
  }

  ## check assay & map order match
  if (!identical(names(map_list), assay_names)) {
    stop(cli::cli_alert_danger(
      paste("Names of", cli::style_bold("map_list"),
            "are not identical to", cli::style_bold("assay_list"),
            sep = " "
      )
    ))
  }

  ## generate sample_metadata_list if not provided
  if (is.null(sample_metadata_list)) {
    sample_metadata_list <-
      lapply(1:list_length, function(x) individual_metadata)
    names(sample_metadata_list) <- assay_names
  }

  ## check assay & metadata order match

  if (!identical(names(sample_metadata_list), assay_names)) {
    stop(cli::cli_alert_danger(
      paste("Names of", cli::style_bold("sample_metadata_list"),
            "are not identical to", cli::style_bold("assay_list"),
            sep = " "
      )
    ))
  }

  if (!is.null(qc_data_list) &
      !identical(names(qc_data_list), assay_names)) {
    stop(cli::cli_alert_danger(
      paste("Names of", cli::style_bold("qc_data_list"),
            "are not identical to rna_raw assays in ",
            cli::style_bold("assay_list"),
            sep = " "
      )
    ))
  }

  ## function to check agreement of sample names and return unmatched
  name_match <- function(list_1, list_2, f1 = "colnames", f2 = "rownames") {
    l1 <- eval(as.name(list_1))
    l2 <- eval(as.name(list_2))
    sapply(1:length(l1),
           function(i) {
             !identical(
               eval(as.name(f1))(l1[[i]]),
               eval(as.name(f2))(l2[[i]])
             )
           }) -> sample_check

    if (any(sample_check)) {
      stop(cli::cli_alert_danger(
        paste0("The ", f1, " in ", cli::style_bold(list_1), " elements: \"",
              paste(names(l1)[sample_check], collapse = ", "),
              "\" should match the ", f2, " of ", cli::style_bold(list_2))
      ))
    }
  }

  ## check assay_list samples
  name_match("assay_list", "sample_metadata_list")

  ## check map_list samples
  name_match("map_list", "sample_metadata_list", f1 = 'rownames')

  ## check qc_data_list samples
  if (!is.null(qc_data_list)) {
    name_match("qc_data_list", "sample_metadata_list",  f1 = 'rownames')
  }

  ## check map_list columns and return unmatched
  sapply(map_list,
         function(x) {
           !identical(colnames(x), c("primary", "colname"))
         }) -> map_check

  if (any(map_check)) {
    stop(cli::cli_alert_danger(
      paste(cli::style_bold("primary"), "and", cli::style_bold("colname"),
            "columns are not found in",
            paste(names(map_list)[map_check], collapse = ", "))
    ))
  }

  ## load reference
  if (organism == "human") {
    EnsDb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  }
  if (organism == "mouse") {
    EnsDb <- org.Mm.eg.db::org.Mm.eg.db
  }

  ## generate summarized experiment object function
  generate_se_list <- function(data_l, metadata_l, assay_types) {
    lapply(1:length(data_l), function(i) {
      name <- names(metadata_l)[[i]]
      type <- assay_types[i]
      data <- data_l[[i]]
      metadata <- metadata_l[[i]]

      data <- data[!duplicated(rownames(data)), ]

      se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(name = data),
        colData = S4Vectors::DataFrame(metadata),
        metadata = metadata
      )

      if (type %in% c("rna_raw", "protein_raw")) {
        cli::cli_alert_success(
          paste("Retrieval of feature annotation:", cli::style_bold(name))
        )
        id_type <- get_ID_type(rownames(data))
        SummarizedExperiment::rowData(se) <-
          magrittr::set_names(
            S4Vectors::DataFrame(gene_id = as.character(rownames(data))),
            id_type
          )
        annot <- rownames(se)
        keytypes <- c(
          "gene_name" = "SYMBOL", "ensembl_gene_id" = "GENEID",
          "uniprot_id" = "UNIPROTID", "entrez_gene_id" = "ENTREZID",
          "gene_biotype" = "GENEBIOTYPE"
        )
        keytype = keytypes[id_type]
        columns <- c("SYMBOL", "GENEID", "ENTREZID", "GENEBIOTYPE")
        AnnotationDbi::select(
          EnsDb,
          keys = annot,
          keytype = keytype,
          columns = columns
        ) -> annot_df
        names(annot_df) <- names(keytypes)[match(names(annot_df), keytypes)]
        annot_df <- annot_df[match(annot, annot_df[, names(keytype)]), ]
        annot_df[, names(keytype)] <- annot
        rowData(se) <- annot_df
      } else {
        cli::cli_alert_info(
          paste(
            "Skipping feature annotation for generic data class:",
            cli::style_bold(name)
          )
        )
      }

      return(se)
    }) -> se_list
    names(se_list) <- names(data_l)
    return(se_list)
  }

  se_list <- generate_se_list(assay_list, sample_metadata_list, assay_types)

  ## add qc data to metadata
  if (!is.null(qc_data_list)) {
    lapply(1:list_length, function(i){
      se <- se_list[[i]]
      qc <- qc_data_list[[i]]
      se@metadata$rna_qc_data <- qc[colnames(se), ]
    }) -> se_list
  }

  ## alter rownames of individual metadata to provided column
  if (!is.null(map_by_column)) {
    rownames(individual_metadata) <- individual_metadata[, map_by_column]
  }

  ## generate multiassay
  suppressMessages({
    dfmap <- MultiAssayExperiment::listToMap(map_list)
    multiassay <- MultiAssayExperiment::MultiAssayExperiment(
      se_list,
      colData = S4Vectors::DataFrame(individual_metadata),
      sampleMap = S4Vectors::DataFrame(dfmap)
    )
  })

  cli::cli_alert_success("MultiAssayExperiment object generated!")

  return(multiassay)
}
