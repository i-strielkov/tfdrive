context('tfpred')

test_that('data intact', {
  load('sysdata.rda')
  expect_equal(names(pathways_hsa), c('gene_id', 'pathway'))
  expect_equal(nrow(pathways_hsa), 50801)
  expect_equal(names(go_hsa), c('gene_id', 'GO_term'))
  expect_equal(nrow(go_hsa), 139684)
  expect_equal(names(tf_names), c('gene', 'gene_id'))
  expect_equal(nrow(tf_names), 174)
  expect_equal(names(tf_families), c('ids', 'family'))
  expect_equal(nrow(tf_families), 81)
  expect_equal(ncol(transfac_df), 202)
  expect_equal(nrow(transfac_df), 13216)
  expect_equal(ncol(chea_df), 200)
  expect_equal(nrow(chea_df), 21585)
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
  tf_id <- tf_names[[1, 'gene_id']]
  tf_name <- tf_names[[1, 'gene']]

  # Data frame to store results
  full_data_df <- data.frame(TF = character(),
                             TF_id = character(),
                             p_share = double(),
                             p_all = double(),
                             p_score = double(),
                             go_share = double(),
                             go_all = double(),
                             go_score = double())

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

  # If GO terms values are missing, insert those from KEGG pathway values
  for (i in seq_len(nrow(full_data_df))){
    if (full_data_df$go_all[i] == 0){
      full_data_df$go_share[i] <- full_data_df$p_share[i]
      full_data_df$go_all[i] <- full_data_df$p_all[i]
      full_data_df$go_score[i] <- full_data_df$p_score[i]
    }
  }

  # Calculate relative interaction frequencies for each family for each group
  tf_group_scores_df <- data.frame(family = character(),
                                   transfac = double(),
                                   chea = double())
  gene_ids_str <- as.character(gene_ids)

  for (f in seq_len(nrow(tf_families))){
    family <- tf_families$family[f]
    tf_ids <- tf_families$ids[[f]]

    # Transfac
    t_value <- 0.0
    t_tf_ids <- names(transfac_df)[names(transfac_df) %in% tf_ids]
    t_g_ids <-  gene_ids_str[gene_ids_str %in% transfac_df$gene_id]
    if (length(t_tf_ids) > 0 & length(t_g_ids) > 0) {
      t_value <- (sum(transfac_df[transfac_df$gene_id %in% t_g_ids,][t_tf_ids]) / length(t_g_ids)) /
        (sum(transfac_df[,t_tf_ids]) / nrow(transfac_df))
    }

    # ChEA
    c_value <- 0.0
    c_tf_ids <- names(chea_df)[names(chea_df) %in% tf_ids]
    c_g_ids <-  gene_ids_str[gene_ids_str %in% chea_df$gene_id]
    if (length(c_tf_ids) > 0 & length(c_g_ids) > 0) {
      c_value <- (sum(chea_df[chea_df$gene_id %in% c_g_ids,][c_tf_ids]) / length(c_g_ids)) /
        (sum(chea_df[,c_tf_ids]) / nrow(chea_df))
    }
    tf_group_scores_df <- rbind(tf_group_scores_df,
                                list(family = family,
                                     transfac = as.double(t_value),
                                     chea = as.double(c_value)),
                                stringsAsFactors=FALSE)
  }
  # Create TF to TF family association table
  tf_family_pairs <- data.frame(family = character(), id = character())
  for (i in seq_len(nrow(tf_families))){
    family <- tf_families$family[i]
    ids <- tf_families$ids[[i]]
    for (g in ids){
      tf_family_pairs <- rbind(tf_family_pairs,
                               list(family = family,
                                    id = g),
                               stringsAsFactors=FALSE)
    }
  }

  # Assign values to full_data_df
  full_data_df$transfac <- 0.0
  full_data_df$chea <- 0.0
  for (i in seq_len(nrow(full_data_df))){
    tf_id <- as.character(full_data_df$TF_id[i])
    if (tf_id %in% tf_family_pairs$id){
      family <- tf_family_pairs[tf_family_pairs$id==tf_id,]$family
      full_data_df$transfac[i] <- tf_group_scores_df[tf_group_scores_df$family==family,]$transfac
      full_data_df$chea[i] <- tf_group_scores_df[tf_group_scores_df$family==family,]$chea
    } else {
      # If TF is not in the TF families, calculate individual values
      if (tf_id %in% names(transfac_df) & length(t_g_ids) > 0) {
        t_value <- (sum(transfac_df[transfac_df$gene_id %in% t_g_ids,][tf_id]) / length(t_g_ids)) /
          (sum(transfac_df[,tf_id]) / nrow(transfac_df))
        full_data_df$transfac[i] <- t_value
      }
      if (tf_id %in% names(chea_df) & length(c_g_ids) > 0) {
        c_value <- (sum(chea_df[chea_df$gene_id %in% c_g_ids,][tf_id]) / length(c_g_ids)) /
          (sum(chea_df[,tf_id]) / nrow(chea_df))
        full_data_df$chea[i] <- c_value
      }
    }
  }

  # Fill zeroes if only one ('chea' or 'transfac') value is available
  full_data_df$transfac <- sapply(seq_len(nrow(full_data_df)), function(i) ifelse(full_data_df$transfac[i]==0,
                                                                                  full_data_df$chea[i],
                                                                                  full_data_df$transfac[i]))
  full_data_df$chea <- sapply(seq_len(nrow(full_data_df)), function(i) ifelse(full_data_df$chea[i]==0,
                                                                              full_data_df$transfac[i],
                                                                              full_data_df$chea[i]))
  # Predictions of Logistic Regression model
  lr_pred <- predict(lr_model,
                     as.matrix(full_data_df[, c('p_share','p_all','p_score',
                                                'go_share','go_all','go_score',
                                                'transfac','chea')]),
                     type="response")
  full_data_df$LR_prob <- as.vector(lr_pred)


  expect_equal(as.character(full_data_df$TF[1]), 'AHR')
  expect_equal(as.character(full_data_df$TF_id[1]), '196')
  expect_equal(round(full_data_df$p_share[1], 6), 0.25)
  expect_equal(full_data_df$p_all[1], 4)
  expect_equal(round(full_data_df$p_score[1], 6), 4.609053)
  expect_equal(round(full_data_df$go_share[1], 6), 0.5)
  expect_equal(full_data_df$go_all[1], 20)
  expect_equal(round(full_data_df$go_score[1], 6), 24.375322)
  expect_equal(round(full_data_df$transfac[1], 6), 0.409608)
  expect_equal(round(full_data_df$chea[1], 6), 1.491789)
  expect_equal(round(full_data_df$LR_prob[1], 6), 0.384725)
})
