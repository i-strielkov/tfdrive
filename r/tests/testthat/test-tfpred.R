context('tfpred')

test_that('data intact', {
  load('sysdata.rda')
  expect_true(all(c("gene_id","pathway") %in% names(pathways_hsa)))
  expect_equal(nrow(pathways_hsa), 52415)
  expect_equal(names(go_hsa), c('gene_id', 'GO_term'))
  expect_equal(nrow(go_hsa), 152579)
  expect_true(all(c('gene', 'index') %in% names(tf_names)))
  expect_equal(nrow(tf_names), 174)
})

test_that("tfpred scores correct", {
  load('sysdata.rda')
  gene_ids <- c(5714, 135114,  26958,   5899, 494514,  26502, 388325,  79007,
                4688,   1428, 259307,   2537,   1117,  55096,    953,    823,
                51031,   6892,    705,  81622,   9466, 374928,  10768,   6774,
                1178,   4494,   2769, 266971,  55651,   5168,   6813,    475,
                22906,   2793,   2766,   5269,  55179,  27339, 140771,   1901,
                2882,  79064,  10184,   6844,  51701,  57827,   1687, 122509,
                4967,  23522,  85437,   9852,   9459,    581, 729230,  57407,
                10985,   7094,   7414,   7791,   2648,   9470)

  # Data frame to store results
  full_data_df <- data.frame(TF = character(),
                             TF_id = character(),
                             p_share = double(),
                             p_all = double(),
                             p_score = double(),
                             go_share = double(),
                             go_all = double(),
                             go_score = double())

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
    pw_df['importance'] <- sapply(pw_df$Freq, function(f) dbeta(f/max(pw_df$Freq), 3, 5) )
    go_df <- as.data.frame(table(case_terms))
    go_df['importance'] <- sapply(go_df$Freq, function(f) dbeta(f/max(go_df$Freq), 3, 5) )

    # Keep only common pathways and terms
    pw_df <- pw_df[pw_df$case_ptws %in% tf_ptws,]
    go_df <- go_df[go_df$case_terms %in% tf_terms,]

    # Calculate relative pathway/term overlap and scores based on 'importance'
    p_common <- go_common <- p_score <- go_score <- p_share <- go_share <- 0.0
    if (p_all > 0 & nrow(pw_df) > 0) {
      p_common <- length(pw_df$case_ptws[pw_df$case_ptws %in% names(tf_pw_count)])
      p_score <- sum(sapply(seq_len(nrow(pw_df)), function(f) pw_df$importance[f] * pw_df$Freq[f]))
      p_share <- p_common/p_all
    }

    if (go_all > 0 & nrow(go_df) > 0) {
      go_common <- length(go_df$case_terms[go_df$case_terms %in% names(tf_term_count)])
      go_score <- sum(sapply(seq_len(nrow(go_df)), function(f) go_df$importance[f] * go_df$Freq[f]))
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

  expect_equal(as.character(full_data_df$TF[1]), 'STAT4')
  expect_equal(as.character(full_data_df$TF_id[1]), '6775')
  expect_equal(round(full_data_df$LR_prob[1], 6), 0.482965)
})
