#' Predicts transcription factors driving differential gene expression based on a list of
#' differentially expressed genes
#'
#' @param gene_ids vector of Entrez Gene IDs of differentially expressed genes
#' @param species three-letter species ID (currently, function works only with human genes ('hsa'))
#'
#' @return data frame with transcription factors ranked according to their probability scores
#'
#' @examples
#' tfpred(c(5714, 135114, 26958, 5899, 494514, 26502, 388325, 79007, 4688, 1428, 259307))
#'
#' @import glmnet
#'
#' @export
tfpred <- function(gene_ids, species = 'hsa'){

  #Checks the input data
  gene_ids <- unlist(gene_ids)
  gene_ids <- unique(gene_ids)
  if (!is.null(dim(gene_ids))) stop('The list of gene IDs should be flat')
  if (length(gene_ids) < 10) stop('The list of genes is too short')
  gene_ids <- as.integer(gene_ids)
  if (any(is.na(gene_ids))) stop("'gene_ids' should be a list of numeric Entrez IDs")
  if (species!='hsa') stop("Only 'hsa' option is currently supported.")

  # Used data
  # pathways_hsa - table of pathways attributed to human genes
  # go_hsa - table of GO terms attributed to human genes
  # tf_names - list of TFs
  # lr_model - logistic regression model

  if (nrow(pathways_hsa[pathways_hsa$gene_id %in% gene_ids,]) < 10) {
    stop ("Too few or none of provided gene IDs were found in the database.")
  }

  # Data frame to store results
  full_data_df <- data.frame(TF = character(),
                             TF_id = character(),
                             p_share = double(),
                             p_all = double(),
                             p_score = double(),
                             go_share = double(),
                             go_all = double(),
                             go_score = double())

  max_pw_count = max(table(pathways_hsa$pathway))
  max_go_count = max(table(go_hsa$GO_term))

  for (i in seq_len(nrow(tf_names))){
    tf_id <- tf_names[[i, 'index']]
    tf_name <- tf_names[[i, 'gene']]

    # Get unique pathways and terms related to the TF and their counts
    tf_ptws <- pathways_hsa[pathways_hsa$gene_id==tf_id,]$pathway
    tf_terms <- go_hsa[go_hsa$gene_id==tf_id,]$GO_term
    tf_pw_count <- table(tf_ptws)
    tf_term_count <- table(tf_terms)
    p_all <- length(tf_pw_count)
    go_all <- length(tf_term_count)

    # Get pathways and terms related to differentially expressed genes
    case_ptws <- pathways_hsa[pathways_hsa$gene_id %in% gene_ids,]$pathway
    case_terms <- go_hsa[go_hsa$gene_id %in% gene_ids,]$GO_term

    # Calculate relative pathway/GO term importance
    pw_df <- as.data.frame(table(case_ptws))
    pw_df['importance'] <- sapply(pw_df$Freq, function(f) dbeta(f/max_pw_count, 2, 2) )
    go_df <- as.data.frame(table(case_terms))
    go_df['importance'] <- sapply(go_df$Freq, function(f) dbeta(f/max_go_count, 2, 2) )

    # Keep only common pathways and terms
    pw_df <- pw_df[pw_df$case_ptws %in% tf_ptws,]
    go_df <- go_df[go_df$case_terms %in% tf_terms,]

    # Calculate relative pathway/term overlap and scores based on 'importance'
    p_common <- go_common <- p_score <- go_score <- p_share <- go_share <- 0.0
    if (p_all > 0 & nrow(pw_df) > 0) {
      p_common <- length(pw_df$case_ptws[pw_df$case_ptws %in% names(tf_pw_count)])
      p_score <- sum(sapply(seq_len(nrow(pw_df)), function(f)  100 * pw_df$importance[f] * pw_df$Freq[f] / length(gene_ids)))
      p_share <- p_common/p_all
    }

    if (go_all > 0 & nrow(go_df) > 0) {
      go_common <- length(go_df$case_terms[go_df$case_terms %in% names(tf_term_count)])
      go_score <- sum(sapply(seq_len(nrow(go_df)), function(f) 100 * go_df$importance[f] * go_df$Freq[f] / length(gene_ids)))
      go_share <- go_common/go_all
    }

    full_data_df <- rbind(full_data_df,
                          data.frame(TF = as.character(tf_name),
                                     TF_id = as.character(tf_id),
                                     p_share = as.double(p_share),
                                     p_all = as.double(p_all),
                                     p_score = as.double(p_score),
                                     go_share = as.double(go_share),
                                     go_all = as.double(go_all),
                                     go_score = as.double(go_score)),
                          stringsAsFactors=FALSE)
  }

  # Check whether there are enough pathways/terms associated with each TF
  full_data_df <- full_data_df[full_data_df$p_all>2,]
  full_data_df <- full_data_df[full_data_df$go_all>2,]

  # Predictions of Logistic Regression model
  lr_pred <- predict(lr_model,
                     as.matrix(full_data_df[, c('p_share',
                                                'p_score',
                                                'go_share',
                                                'go_score')]),
                     type="response")
  full_data_df$LR_prob <- as.vector(lr_pred)

  # Preparing results
  full_data_df <- full_data_df[,c('TF', 'TF_id', 'LR_prob')]
  full_data_df <- full_data_df[!full_data_df$TF_id %in% gene_ids,]
  full_data_df <- full_data_df[order(full_data_df$LR_prob, decreasing = TRUE),]
  rownames(full_data_df) <- 1:nrow(full_data_df)

  full_data_df
}
