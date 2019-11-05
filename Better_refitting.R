library(tidyverse)
library(magrittr)
library(MutationalPatterns)
library(pracma)
library(vegan)

#Slightly modifief from MutationalPatterns to plot contribution of normalized (colSum 1) signatures
plot_contribution_nonmf = function (contribution, index = c(), coord_flip = FALSE, 
          mode = "relative", palette = c()) 
{
  if (!(mode == "relative" | mode == "absolute")) 
    stop("mode parameter should be either 'relative' or 'absolute'")
  if (length(index > 0)) {
    contribution = contribution[, index]
  }
  Sample = NULL
  Contribution = NULL
  Signature = NULL
  if (mode == "relative") {
    m_contribution = contribution %>% 
      as.data.frame() %>% 
      rownames_to_column("Signature") %>% 
      tidyr::gather(key = "Sample", value = "Contribution", -Signature)
    plot = ggplot(m_contribution, aes(x = factor(Sample), y = Contribution, fill = factor(Signature), order = Sample)) +
      geom_bar(position = "fill", stat = "identity", colour = "black") + 
      labs(x = "", y = "Relative contribution") + 
      theme_classic() + 
      theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 90)) + 
      theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
  }
  else {
    m_contribution =    m_contribution = contribution %>% 
      as.data.frame() %>% 
      rownames_to_column("Signature") %>% 
      tidyr::gather(key = "Sample", value = "Contribution", -Signature)
      plot = ggplot(m_contribution, aes(x = factor(Sample), y = Contribution, fill = factor(Signature), order = Sample)) + 
        geom_bar(stat = "identity", colour = "black") + 
        labs(x = "", y = "Absolute contribution \n (no. mutations)") + 
        theme_classic() + 
        theme(panel.grid.minor.x = element_blank(), panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 90)) + 
        theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank())
  }
  if (length(palette) > 0) 
    plot = plot + scale_fill_manual(name = "Signature", values = palette)
  else plot = plot + scale_fill_discrete(name = "Signature")
  if (coord_flip) 
    plot = plot + coord_flip() + xlim(rev(levels(factor(m_contribution$Sample))))
  else plot = plot + xlim(levels(factor(m_contribution$Sample)))
  return(plot)
}

#Function to plot the cosine similarities between samples or between signatures themselves.
plot_inner_cosheat = function(mat){
  cos_sim_matrix = cos_sim_matrix(mat, mat)
  
  clust = as.dist(1-cos_sim_matrix) %>% hclust()
  #clust = dist(t(mat)) %>% hclust()
  order = colnames(mat)[clust$order]
  
  
  dhc = as.dendrogram(clust)
  ddata = dendro_data(dhc, type = "rectangle")
  sig_dendrogram = ggplot(segment(ddata)) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
    theme_dendro()
  
  sample_order = order
  Cosine.sim = NULL
  Signature = NULL
  Sample = NULL
  x = NULL
  y = NULL
  xend = NULL
  yend = NULL
  cos_sim_matrix.m = cos_sim_matrix %>% as.data.frame() %>% rownames_to_column("Sample") %>% gather(-Sample, key = "Sample2", value = "Cosine.sim")
  cos_sim_matrix.m$Sample2 = factor(cos_sim_matrix.m$Sample2, levels = order)
  cos_sim_matrix.m$Sample = factor(cos_sim_matrix.m$Sample, levels = sample_order)
  heatmap = ggplot(cos_sim_matrix.m, aes(x = Sample2, y = Sample, 
                                         fill = Cosine.sim, order = Sample)) + geom_tile(color = "white") + 
    scale_fill_distiller(palette = "YlGnBu", direction = 1, 
                         name = "Cosine \nsimilarity", limits = c(0, 1.000001)) + 
    theme_bw() + theme(axis.text.x = element_text(angle = 90, 
                                                  hjust = 1)) + labs(x = NULL, y = NULL)
  
  heatmap = plot_grid(sig_dendrogram, heatmap, align = "v", rel_heights = c(0.2, 1), axis = "lr", nrow = 2, scale = c(1,1))
  return(heatmap)
}


#Faster mut_matrix function
C_TRIPLETS = c(
  "ACA", "ACC", "ACG", "ACT",
  "CCA", "CCC", "CCG", "CCT",
  "GCA", "GCC", "GCG", "GCT",
  "TCA", "TCC", "TCG", "TCT")

T_TRIPLETS = c(
  "ATA", "ATC", "ATG", "ATT",
  "CTA", "CTC", "CTG", "CTT",
  "GTA", "GTC", "GTG", "GTT",
  "TTA", "TTC", "TTG", "TTT")

CONTEXTS_96 = c(rep(C_TRIPLETS, 3), rep(T_TRIPLETS, 3))
SUBSTITUTIONS = c('C>A','C>G','C>T','T>A','T>C','T>G')
SUBSTITUTIONS_96 = rep(SUBSTITUTIONS, each=16)
TRIPLETS_96 = paste(substr(CONTEXTS_96,1,1), "[", SUBSTITUTIONS_96, "]", substr(CONTEXTS_96,3,3), sep = "")
mut_96_occurrences2 = function(type_context, gr_sizes){
  cats = tibble("categories" = factor(TRIPLETS_96, levels = TRIPLETS_96))
  full_context = paste0(substr(type_context$context, 1, 1), "[", type_context$types, "]", substr(type_context$context, 3, 3)) %>% factor(levels = TRIPLETS_96)
  
  if (is.null(names(gr_sizes))){
    n = length(gr_sizes)
    names(gr_sizes) = seq(n)
  }
  
  sample_vector = rep(names(gr_sizes), gr_sizes) %>% factor(levels = names(gr_sizes))
  
  #Count the mutations per type and per sample
  counts = tibble("categories" = full_context, "sample" = sample_vector) %>% dplyr::filter(!is.na(categories)) %>% dplyr::group_by(categories, sample) %>% dplyr::summarise(count = n())
  counts = left_join(cats, counts, by = "categories")
  
  #Transform the data into a mutation matrix
  counts = spread(counts, key = sample, value = count, fill = 0)
  unnecesary_cols = which(colnames(counts) == "<NA>")
  mut_mat = as.matrix(counts[,-c(1, unnecesary_cols)])
  rownames(mut_mat) = counts$categories
  return(mut_mat)
  
  # cats = tibble("categories" = TRIPLETS_96)
  # full_context = paste0(substr(type_context$context, 1, 1), "[", type_context$types, "]", substr(type_context$context, 3, 3))
  #counts = table(full_context) %>% enframe(name = "categories", value = "count") %>% as.data.frame()
  #counts = left_join(cats, counts, by = "categories")
  #counts[is.na(counts$count), "count"] = 0
  #vector = counts$count
  #names(vector) = counts$categories
  #return(vector)
}

mut_matrix2 = function (grl, ref_genome) {
  if (class(grl)[[1]] == "CompressedGRangesList"){
    gr_sizes = elementNROWS(grl)
    gr = unlist(grl)
  } else if (class(grl)[[1]] == "GRanges"){
    gr = grl
    gr_sizes = length(gr)
    names(gr_sizes) = "My_sample"
  } else{
    stop("This function requires either a GRanges or a CompressedGRangesList object as its first argument.")
  }
  type_context = type_context(gr, ref_genome)
  mut_mat = mut_96_occurrences2(type_context, gr_sizes)
  return(mut_mat)
}


plot_96_sigs = function(signatures){
  sig_figs = lapply(1:ncol(signatures), function(x){
    plot_96_profile(signatures[,x, drop = F], ymax = 0.4)
  })
  return(sig_figs)
}


#Slightly modified from mutational patterns, so that non 96 feature signatures can be used.
fit_to_signatures = function (mut_matrix, signatures)
{
  n_feats = nrow(mut_matrix)
  n_samples = dim(mut_matrix)[2]
  n_signatures = dim(signatures)[2]
  lsq_contribution = matrix(NA, nrow = n_signatures, ncol = n_samples)
  lsq_reconstructed = matrix(NA, nrow = n_feats, ncol = n_samples)
  for (i in 1:ncol(mut_matrix)) {
    y = mut_matrix[, i]
    lsq = lsqnonneg(signatures, y)

    #ratio = sum(y) / sum(lsq$x)
    #lsq$x = lsq$x * ratio

    lsq_contribution[, i] = lsq$x
    lsq_reconstructed[, i] = signatures %*% as.matrix(lsq$x)
  }
  sample_names = colnames(mut_matrix)
  signature_names = colnames(signatures)
  mut_type_names = rownames(signatures)
  colnames(lsq_contribution) = sample_names
  rownames(lsq_contribution) = signature_names
  colnames(lsq_reconstructed) = sample_names
  rownames(lsq_reconstructed) = mut_type_names
  res = list(lsq_contribution, lsq_reconstructed)
  names(res) = c("contribution", "reconstructed")
  return(res)
}

#Function to merge similar signatures based on their cosine similarity
merge_signatures = function(signatures, sig_sim_lim = 0.8){
  my_signatures = signatures
  nr_lowercase = colnames(my_signatures) %>% grepl("_", .) %>% sum()
  if (nr_lowercase > 0){
    stop("Please remove all lowercases from your signature names.")
  }
  
  n_sigs = ncol(signatures)
  for (i in 2:n_sigs){
    sig_sig = cos_sim_matrix(my_signatures, my_signatures)
    diag(sig_sig) = 0
    
    max = max(sig_sig)
    if (max <= sig_sim_lim){#Stop merging signatures when the maximum similarity is below the cutoff.
      break
    }
    max_index = order(sig_sig, decreasing = T)[1]
    max_loc = arrayInd(max_index, dim(sig_sig), useNames = T)
    sigs_left = my_signatures[,-max_loc, drop = F]
    sigs_to_combi = my_signatures[,max_loc, drop = F]
    weights = colnames(sigs_to_combi) %>% str_count("_") + 1 #Signatures that have already been merged and thus exist of multiple signatures are weighted accordingly.
    combi_sig = sigs_to_combi %*% diag(weights)
    combi_sig = combi_sig %>% rowSums() %>% matrix()
    combi_sig = combi_sig / sum(weights)
    colnames(combi_sig) = paste(colnames(sigs_to_combi), collapse = "_")
    
    my_signatures = cbind(sigs_left, combi_sig)
    merged_sig_names = paste0(colnames(sigs_to_combi), collapse = " ")
    print(paste0("Combined the following two signatures: ", merged_sig_names))
  }
  return(my_signatures)
}

fit_to_signatures_permuted = function(mut_mat, signatures, n_boots = 1000, max_delta = 0.05, method = c("selection", "selection_build", "mutpatterns_default", "mutpatterns_strict")){
  method = match.arg(method)
  
  min_nr_muts = colSums(mut_mat) %>% min()
  if (min_nr_muts <= 10){
    warning("At least one of your samples has less than 10 mutations. This will negatively impact the signature refitting. Please consider removing or pooling this sample.")
  }
  
  
  #Perform permutations. Both margins are kept fixed.
  perm = permatfull(mut_mat, times = n_boots)
  
  sig_names_tb = tibble("sigs" = colnames(signatures))
  contri_list = vector("list", n_boots)
  for (i in 1:n_boots){
    
    #Get a single permutation of the matrix
    mut_mat_perm = perm$perm[[i]]
    
    
    
    if (method == "selection"){
      refit_out = fit_to_signatures_selection(mut_mat_perm, signatures, max_delta = max_delta)
      contri = refit_out$fit_res$contribution
    }
    else if (method == "selection_build"){
      fit_res = fit_to_signatures_regularized(mut_mat_resampled, signatures, regularized = T, min_delta = max_delta)
      contri = fit_res$contribution
    }
    else if (method == "mutpatterns_default"){
      fit_res = fit_to_signatures(mut_mat_perm, signatures)
      contri = fit_res$contribution
    }
    else if (method == "mutpatterns_strict"){
      fit_res = fit_to_signatures(mut_mat_perm, signatures)
      index = rowSums(fit_res$contribution >= 10) != 0 #Check whether a signature has at least 10 mutations in a single sample
      contri = fit_res$contribution[index,]
    }
    colnames(contri) = paste(colnames(contri), i, sep = "_")
    contri_tb = contri %>% as.data.frame() %>% rownames_to_column("sigs")
    contri_tb = left_join(sig_names_tb, contri_tb, by = "sigs") %>% dplyr::select(-sigs)
    contri_list[[i]] = contri_tb
    
    if (i%%10 == 0){
      print(paste0("Performed ", i, " of ", n_boots, " iterations"))
    }
  }
  contri_boots = do.call(cbind, contri_list) %>% t()
  colnames(contri_boots) = sig_names_tb$sigs
  return(contri_boots)
}

#Function to resample a mutation matrix
resample_mut_mat = function(mut_mat){
  mut_mat_resampled = apply(mut_mat, 2, function(x){
    total_muts = sum(x)
    sample_weights = x / total_muts
    feature_rows = sample(1:length(x), total_muts, replace = T, prob = sample_weights)
    row_counts = table(feature_rows)
    index = as.numeric(names(row_counts))
    x[index] = as.vector(row_counts)
    x[-index] = 0
    return(x)
  })
  return(mut_mat_resampled)
}


#Function to plot how often each signature was selected.
plot_fraction_contri = function(contri_boots){
  samples = str_replace(rownames(contri_boots), "_\\d*", "")
  contri_boots_l = split(as.data.frame(contri_boots), samples, drop = F)
  nr_na_l = purrr::map(contri_boots_l, function(x) {
    purrr::map_int(x, ~sum(is.na(.)))
  })
  nr_na = do.call(rbind, nr_na_l)
  nr_boots = nrow(contri_boots) / length(unique(samples))
  ratio_sel = 1 - nr_na / nr_boots
  ratio_sel = ratio_sel[,colSums(ratio_sel) != 0, drop = F]
  ratio_sel_df = ratio_sel %>% 
    as.data.frame() %>% 
    rownames_to_column("sample") %>% 
    tidyr::gather(key = "sig", value = "fraction_selected", -sample) %>% 
    dplyr::mutate(sig = factor(sig, levels = unique(sig)))

  ratio_sel_df_l = split(ratio_sel_df, ceiling(group_indices(ratio_sel_df,sample)/5))
  figs = purrr::map(ratio_sel_df_l, function(x) {ggplot(data = x, aes(x = sig, y = fraction_selected)) +
    geom_bar(stat = "identity") +
    facet_grid(sample ~ .) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Signature", y = paste0("Selected (n = ", nr_boots, ")")) +
    theme(axis.text.x = element_text(angle = 90, size = 18), text = element_text(size=24), strip.text.y = element_text(size = 14))
    })
  return(figs)
}

#Function to get signatures that were selected during the bootstrapping at least x percent of the time.
get_common_signatures = function(contri_boots, signatures, frac_sel_lim = 0.5){
  nr_na = apply(contri_boots, 2, function(x) sum(is.na(x)))
  nrow_contri = nrow(contri_boots)
  ratio_sel = 1 - nr_na / nrow_contri
  sel_sigs_index = ratio_sel > frac_sel_lim
  sel_signatures = signatures[, sel_sigs_index, drop = F]
  return(sel_signatures)
}

#Helper function to replace NAs and remove unselected signatures from a bootstrapped contributions table.
clean_contri_boots = function(contri_boots){
  contri_boots[is.na(contri_boots)] = 0
  contri_boots = contri_boots[,colSums(contri_boots) != 0, drop = F]
  return(contri_boots)
}

#Function to plot bootstrapped signature contribution.
plot_boots_contri = function(contri_boots, mode = c("absolute", "relative")){
  mode = match.arg(mode)
  contri_boots = clean_contri_boots(contri_boots)
  ylab_text = "Mean nr contributed mutations"
  jitter_height = 0.2
  if (mode == "relative"){
    contri_boots = contri_boots / rowSums(contri_boots)
    ylab_text = "Relative mutation contribution"
    jitter_height = 0.02
  }
  contri_tb = contri_boots %>% as.data.frame() %>% rownames_to_column("exp") %>% gather(key = "sig", value = "contri", -exp) %>% mutate(sample = gsub("_[^_]+$", "", exp), sig = factor(sig, levels = unique(sig)))
  contri_tb_l = split(contri_tb, ceiling(group_indices(contri_tb,sample)/5))
  contri_boots_fig1 = purrr::map(contri_tb_l, function(x) {ggplot(data = x, aes(x = sig, y = contri, color = sig)) +
    geom_jitter(stat = "identity", height = jitter_height, size = 0.3) +
    facet_grid(sample ~ .) +
    scale_color_discrete(guide = F) +
    labs(x = "Signature", y = ylab_text) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, size = 18), text = element_text(size=24), strip.text.y = element_text(size = 14))
    })

  contri_tb2 = contri_tb %>% 
    group_by(sample, sig) %>% 
    dplyr::summarise(mean = mean(contri), lower = quantile(contri, 0.025), upper = quantile(contri, 0.975)) %>% 
    ungroup()
  contri_tb2_l = split(contri_tb2, ceiling(group_indices(contri_tb2,sample)/5))
  contri_boots_fig2 = purrr::map(contri_tb2_l, function(x) {ggplot(data = x, aes(x = sig, y = mean, fill = sig)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
    facet_grid(sample ~ .) +
    scale_fill_discrete(guide = F) +
    labs(x = "Signature", y = ylab_text) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, size = 18), text = element_text(size=24), strip.text.y = element_text(size = 14))
    })

  contri_boots_fig3 = purrr::map(contri_tb_l, function(x) {ggplot(data = x, aes(x = sig, y = contri, fill = sig)) +
    geom_violin(aes(fill = sig), scale = "width") +
    facet_grid(sample ~ .) +
    labs(x = "Signature", y = ylab_text) +
    scale_fill_discrete(guide = F) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, size = 18), text = element_text(size=24), strip.text.y = element_text(size = 14))
    })
  
  
  return(list(contri_boots_fig1, contri_boots_fig2, contri_boots_fig3))
}


plot_perm_contri = function(contri_boots, diff_sigs, mode = c("absolute", "relative")){
  mode = match.arg(mode)
  
  #Get cleaned up and relative data
  contri_boots = clean_contri_boots(contri_boots)
  
  #Add total mutations per sample to diff_sigs
  diff_sigs %<>% 
    dplyr::group_by(sample) %>% 
    dplyr::mutate(total_muts = sum(main_contri)) %>% 
    ungroup()
  
  
  #Set settings specific for relative or absolute mode.  
  if (mode == "relative"){
    ylab_text = "Relative mutation contribution"
    diff_sigs %<>% dplyr::mutate(plot_contri = rel_contri, mean_contri = mean_contri / total_muts, lower = lower / total_muts, upper = upper / total_muts)
    contri_boots = contri_boots / rowSums(contri_boots)
    ann = diff_sigs %>% dplyr::mutate(max_plot = ifelse(max > main_contri, max / total_muts, main_contri / total_muts))
  } else{
    ylab_text = "Mutation contribution"
    diff_sigs %<>% dplyr::mutate(plot_contri = main_contri)
    ann = diff_sigs %>% dplyr::mutate(max_plot = ifelse(max > main_contri, max, main_contri))
  }
  max_y = 1.05 * max(ann$max_plot)
  
  #Finish pvalue annotation
  ann %<>%
    dplyr::mutate(my_text = paste0("pval: ", pval)) %>% 
    dplyr::group_by(sig) %>% 
    dplyr::mutate(contri = max_y, mean_contri = contri) %>% 
    ungroup() %>% 
    dplyr::select(sig, sample, my_text, contri, mean_contri)
  
  #Calculate number of boots
  n_samples = rownames(contri_boots) %>% gsub("_[^_]+$", "", .) %>% unique() %>% length()
  nrow_contri = nrow(contri_boots)
  n_boots = nrow_contri / n_samples
  
  contri_tb = contri_boots %>% as.data.frame() %>% rownames_to_column("exp") %>% gather(key = "sig", value = "contri", -exp) %>% mutate(sample = gsub("_[^_]+$", "", exp), sig = factor(sig, levels = unique(sig)))
  contri_perm_fig1 = ggplot(data = contri_tb, aes(x = sample, y = contri, color = sig)) +
    geom_jitter(stat = "identity", height = 0, size = 0.3) +
    geom_point(data = diff_sigs, aes(x = sample, y = plot_contri, fill = "Measured"), color = "black", shape = 21) +
    geom_text(data = ann, label = ann$my_text, color = "black") +
    facet_grid(sig ~ ., scales = "free_y") +
    scale_fill_manual(values = "black") +
    guides(color = guide_legend(override.aes = list(size = 5)), fill = guide_legend(override.aes = list(size = 5))) +
    labs(x = "Sample", y = ylab_text, color = paste0("H0 (permutations n = " , n_boots, ")"), fill = "") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, size = 18), text = element_text(size=24)) +
    coord_cartesian(ylim = c(0, max_y))
  
  
  #contri_tb2 = contri_tb %>% group_by(sample, sig) %>% dplyr::summarise(mean = mean(contri), lower = quantile(contri, 0.05), upper = quantile(contri, 0.95))
  contri_perm_fig2 = ggplot(data = diff_sigs, aes(x = sample, y = mean_contri, fill = sig)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
    geom_point(aes(y = plot_contri, color = "Measured")) +
    geom_text(data = ann, label = ann$my_text, color = "black") +
    facet_grid(sig ~ ., scales = "free_y") +
    scale_color_manual(values = "black") +
    labs(x = "Sample", y = ylab_text, fill = paste0("H0 (permutations n = " , n_boots, ")"), color = "") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, size = 18), text = element_text(size=24)) +
    coord_cartesian(ylim = c(0, max_y))
  
  
  return(list(contri_perm_fig1, contri_perm_fig2))
}


#Function to calculate difference between samples in bootstrapped sig estimates.
calc_diff_sigs_boots = function(contri_boots, sample1, sample2){
  contri_boots = clean_contri_boots(contri_boots)
  contri_boots_rel = contri_boots / rowSums(contri_boots)
  
  contri_tb = contri_boots_rel %>% as.data.frame() %>% rownames_to_column("exp") %>% gather(key = "sig", value = "contri", -exp) %>% mutate(sample = gsub("_[^_]+$", "", exp), sig = factor(sig, levels = unique(sig)))
  group1 = contri_tb %>% dplyr::filter(sample == sample1)
  group2 = contri_tb %>% dplyr::filter(sample == sample2)
  diff = group1$contri - group2$contri
  diff_tb = tibble("sig" = group1$sig, diff) %>% dplyr::group_by(sig) %>% dplyr::summarise(mean_diff = mean(diff), lower = quantile(diff, 0.025), upper = quantile(diff, 0.975), n = n(), nr_offdirection = ifelse(mean_diff > 0, sum(diff <= 0), sum(diff >= 0)), chance = nr_offdirection / n)
  return(diff_tb)
}

#Function to calculate difference between samples in permutated sig estimates.
calc_diff_sigs_perm = function(contri_boots, mut_mat, signatures, max_delta = 0.05, method = c("selection", "mutpatterns_default", "mutpatterns_strict")){
  method = match.arg(method)
  
  contri_boots = clean_contri_boots(contri_boots)
  contri_tb = contri_boots %>% as.data.frame() %>% rownames_to_column("exp") %>% gather(key = "sig", value = "contri", -exp) %>% mutate(sample = gsub("_[^_]+$", "", exp), sig = factor(sig, levels = unique(sig)))
  
  if (method == "selection"){
    refit_out = fit_to_signatures_selection(mut_mat, signatures, max_delta = max_delta)
    contri = refit_out$fit_res$contribution
  } else if (method == "mutpatterns_default"){
    fit_res = fit_to_signatures(mut_mat, signatures)
    contri = fit_res$contribution
  } else if (method == "mutpatterns_strict"){
    fit_res = fit_to_signatures(mut_mat, signatures)
    index = rowSums(fit_res$contribution >= 10) != 0 #Check whether a signature has at least 10 mutations in a single sample
    contri = fit_res$contribution[index,]
  }
  #refit_out = fit_to_signatures_selection(mut_mat, signatures, max_delta = max_delta)
  
  main_contri = contri %>% as.data.frame() %>% rownames_to_column("sig")
  all_id_sigs = tibble("sig" = colnames(contri_boots))
  main_contri = left_join(all_id_sigs, main_contri, by = "sig") %>% gather(key = "sample", value = "main_contri", -sig) %>% mutate(sig = factor(sig, levels = unique(sig)))
  main_contri[is.na(main_contri)] = 0
  
  contri_tb = left_join(contri_tb, main_contri, by = c("sig", "sample"))
  diff_sigs = contri_tb %>% dplyr::group_by(sig, sample) %>% dplyr::summarise(main_contri = main_contri[1], mean_contri = mean(contri), lower = quantile(contri, 0.025), upper = quantile(contri, 0.975), min = min(contri), max = max(contri), n = n(), nr_offdirection =  ifelse(main_contri > mean_contri, sum(main_contri <= contri), sum(main_contri >= contri)), pval = 2*nr_offdirection / n)
  diff_sigs %<>% dplyr::ungroup() %>% mutate(fdr = p.adjust(pval, method = "fdr")) %>% dplyr::group_by(sample) %>% dplyr::mutate(rel_contri = main_contri / sum(main_contri)) %>% ungroup()
  return(diff_sigs)
}

#Helper function to plot the correlation of signature contribution between the signatures.
plot_sig_contri_cor_sample = function(clean_contri_boots, sample){
  sig_cor = cor(clean_contri_boots)
  sig_cor_tb = sig_cor %>% as.data.frame() %>% rownames_to_column("sig_row") %>% gather(key = "sig_column", value = "cor", -sig_row)
  sig_cor_tb %<>% mutate(sig_row = factor(sig_row, levels = unique(sig_row)), sig_column = factor(sig_column, levels = unique(sig_column)))
  sig_cor_fig = ggplot(data = sig_cor_tb, aes(x = sig_column, y = sig_row, fill = cor), order = NULL) +
    geom_tile(color = "white") +
    scale_fill_distiller(palette = "RdYlBu", direction = -1, name = "Correlation") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size=24), panel.border = element_rect(colour = "black", fill=NA, size=1)) +
    labs(x = NULL, y = NULL, title = sample)
  return(sig_cor_fig)
}


#Function to plot the correlation of signature contribution between the signatures. This is done for the entire contribution table and per sample seperately.
plot_sig_contri_cor = function(contri_boots){
  contri_boots = clean_contri_boots(contri_boots)
  
  samples = rownames(contri_boots) %>% gsub("_[^_]+$", "", .)
  unique_samples = unique(samples)
  n_samples = length(unique_samples)
  figs = vector("list", n_samples +1)
  for (i in 1:n_samples){
    sample = unique_samples[i]
    index = sample == samples
    contri_boots_sample = contri_boots[index,, drop = F]
    fig = plot_sig_contri_cor_sample(contri_boots_sample, sample)
    figs[[i]] = fig
  }
  figs[[n_samples+1]] = plot_sig_contri_cor_sample(contri_boots, "all")
  return(figs)
}

#Function to plot the number of signatures that was selected during the different bootstrapping rounds
plot_nr_signatures_selected = function(contri_boots){
  samples = contri_boots %>% rownames() %>% gsub("_[^_]+$", "", .) %>% unique()
  n_samples = length(samples)
  selected_sigs_tb_l = vector("list", n_samples)
  for (i in seq(1, n_samples)){
    index = seq(i, nrow(contri_boots), n_samples)
    contri_boots_sample = contri_boots[index,]
    nr_removed_sigs = is.na(contri_boots_sample) %>% rowSums()
    nr_selected_sigs = ncol(contri_boots_sample) - nr_removed_sigs
    nr_selected_sigs_tb = tibble("nr_selected_sigs" = nr_selected_sigs, "sample" = samples[i])
    selected_sigs_tb_l[[i]] = nr_selected_sigs_tb
  }
  nr_selected_sigs_tb = do.call(rbind, selected_sigs_tb_l)
  nr_selected_sigs_tb_l = split(nr_selected_sigs_tb, ceiling(group_indices(nr_selected_sigs_tb,sample)/5))
  
  fig = purrr::map(nr_selected_sigs_tb_l, function(x) {ggplot(data = x, aes(x = nr_selected_sigs)) +
    geom_bar(fill = "darkred") +
    facet_grid(sample ~ .) +
    labs(x = "Number of selected signatures", y= "Frequency") +
    theme(text = element_text(size=24), strip.text.y = element_text(size = 14))
  })
  return(fig)
}


#Function to get the cosine similarity between a mutation matrix and its reconstructed fit.
get_cos_sim_ori_vs_rec = function(mut_mat, fit_res){
  cos_sim_all = cos_sim_matrix(mut_mat, fit_res$reconstructed)
  cos_sim = diag(cos_sim_all)
  mean_cos_sim = mean(cos_sim)
}


#Fucntion to show which signatures were removed in the selection process and how this affected cosine similarity
plot_sim_decay = function(sims, removed_sigs, max_delta){
  sims = sims[!isEmpty(sims)] %>% unlist()
  removed_sigs = removed_sigs[!isEmpty(removed_sigs)] %>% unlist()
  tb = tibble("Cosine_similarity" = sims, "Removed_signatures" = factor(removed_sigs, levels = removed_sigs))
  
  #Determine if the final removed signature exceeded the cutoff.
  sims_l = length(sims)
  col = rep("low_delta", sims_l)
  final_delta = sims[sims_l-1] - sims[sims_l]
  if (final_delta > max_delta){
    col[sims_l] = "high_delta"
  }
  
  fig = ggplot(data = tb, aes(x = Removed_signatures, y = Cosine_similarity, fill = col)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(limits = c("low_delta", "high_delta"), values = c("grey", "red"), guide = F) +
    labs(x = "Removed signatures", y = paste0("Cosine similarity (max delta: ", max_delta, ")")) +
    theme(axis.text.x = element_text(angle = 90, size = 18), text = element_text(size=24))
  return(fig)
}

#Function to plot the cosine similarities between a mutation matrix and the mean signature contributions based on bootstrapping.
plot_cosine_bootstrapped = function(mut_mat, contri_boots, signatures){
  contri_boots[is.na(contri_boots)] = 0
  contri_boots %<>% as.data.frame() %>% rownames_to_column("exp") %>% mutate(sample = gsub("_[^_]+$", "", exp))
  mean_contri = contri_boots %>% dplyr::group_by(sample) %>% dplyr::select(-exp) %>% dplyr::summarise_all("mean")
  mean_contri_m = mean_contri %>% dplyr::select(-sample) %>% as.matrix() %>% t()
  colnames(mean_contri_m) = mean_contri$sample
  reconstructed = signatures %*% mean_contri_m
  reconstructed = reconstructed[,colnames(mut_mat), drop = F]
  cos_sim_all = cos_sim_matrix(mut_mat, reconstructed)
  cos_sim = diag(cos_sim_all) %>% enframe(name = "sample", value = "cos_sim")
  
  orivsrec_fig = ggplot(cos_sim, aes(y = cos_sim, x = sample)) +
    geom_bar(stat = "identity", fill = "skyblue3") +
    coord_cartesian(ylim = c(0.6, 1), expand = F) +
    labs(y = "Cosine similarity\n original VS reconstructed", x = "") +
    theme_classic() +
    theme(text = element_text(size = 20), axis.text.x = element_text(angle = 90)) +
    geom_hline(aes(yintercept = 0.95))
  return(orivsrec_fig)
}

#Taken from stack overflow user whuber https://stats.stackexchange.com/questions/321542/how-can-i-draw-a-value-randomly-from-a-kernel-density-estimate?rq=1
rdens <- function(n, density=z, data=x, kernel="gaussian") {
  width <- z$bw                              # Kernel width
  rkernel <- function(n) rnorm(n, sd=width)  # Kernel sampler
  sim = sample(x, n, replace=TRUE) + rkernel(n)    # Here's the entire algorithm
  sim = as.vector(sim)
}

simulate_new_signatures = function(signatures, n){
  new_sigs = apply(signatures, 1, function(x){
    # Compute a kernel density estimate.
    # It returns a kernel width in $bw as well as $x and $y vectors for plotting.
    z <- density(x, bw="SJ", kernel="gaussian")
    y = rdens(n, z, x)
    return(y)
  })
  new_sigs %<>% t() %>% prop.table(2)
  colnames(new_sigs) = paste0("SimulatedSig", 1:ncol(new_sigs))
  return(new_sigs)
}

#Taken from stack overflow user whuber https://stats.stackexchange.com/questions/321542/how-can-i-draw-a-value-randomly-from-a-kernel-density-estimate?rq=1
plot_sig_simulation = function(signatures_row){
  x = signatures_row
  z <- density(x, bw="SJ", kernel="gaussian")
  y = rdens(n, z, x)
  
  # Plot the sample.
  h.density <- hist(y, breaks=60, plot=FALSE)
  
  # Plot the KDE for comparison.
  h.sample <- hist(x, breaks=h.density$breaks, plot=FALSE)
  
  # Display the plots side by side.
  histograms <- list(Sample=h.sample, Density=h.density)
  y.max <- max(h.density$density) * 1.25
  par(mfrow=c(1,2))
  for (s in names(histograms)) {
    h <- histograms[[s]]
    plot(h, freq=FALSE, ylim=c(0, y.max), col="#f0f0f0", border="Gray",
         main=paste("Histogram of", s))
    lines(z$x, z$y, col="Red", lwd=2)# KDE of data
    
  }
  par(mfrow=c(1,1))
}

####Refit signatures with signature selection
fit_to_signatures_selection = function(mut_mat, signatures, max_delta = 0.05, nr_simulated_sigs = 0, signatures_simulation_basis = F){
  
  #if (dim(mut_mat)[1] != 96) 
  #    stop(paste("Mutation count matrix input should have", 
  #               "dimensions 96 X n samples"))
  #if (dim(signatures)[1] != 96) 
  #    stop("Signatures input should have dimensions 96 X n signatures")
  
  #Add simulated signatures to be used for refitting.
  if (nr_simulated_sigs != 0){
    if (signatures_simulation_basis == F){
      signatures_simulation_basis = signatures
    }
    sim_sigs = simulate_new_signatures(signatures_simulation_basis, n = nr_simulated_sigs)
    signatures = cbind(signatures, sim_sigs)
  }
  
  #Remove signatures with zero contribution
  fit_res = fit_to_signatures(mut_mat, signatures)
  sig_pres = rowSums(fit_res$contribution) != 0
  my_signatures_total = signatures[,sig_pres, drop = F]
  nsigs = ncol(my_signatures_total)
  
  #perform signature selection per sample
  all_results = vector("list", ncol(mut_mat))
  for (i in seq(1, ncol(mut_mat))){
    my_signatures = my_signatures_total
    mut_mat_sample = mut_mat[, i, drop = F]
    #Fit again
    fit_res = fit_to_signatures(mut_mat_sample, my_signatures)
    sim = get_cos_sim_ori_vs_rec(mut_mat_sample, fit_res)
    
    #Keep track of the cosine similarity and which signatures are removed.
    sims = vector("list", nsigs)
    sims[[1]] = sim
    removed_sigs = vector("list", nsigs)
    removed_sigs[[1]] = "None"
    
    #Sequentially remove the signature with the lowest contribution
    for (j in seq(2, nsigs)){
      
      #Remove signature with the weakest relative contribution
      contri_order = fit_res$contribution %>% prop.table(2) %>% rowSums() %>% order()
      weakest_sig_index = contri_order[1]
      weakest_sig = colnames(my_signatures)[weakest_sig_index]
      removed_sigs[[j]] = weakest_sig
      signatures_sel = my_signatures[,-weakest_sig_index, drop = F]
      
      
      #Fit with new signature selection
      fit_res = fit_to_signatures(mut_mat_sample, signatures_sel)
      sim_new = get_cos_sim_ori_vs_rec(mut_mat_sample, fit_res)
      
      if (is.nan(sim_new) == TRUE){
        sim_new = 0
        warning("New similarity between the original and the reconstructed spectra after the removal of a signature was NaN. It has been converted into a 0. This happened with the following fit_res:")
        print(fit_res)
      }
      sims[[j]] = sim_new
      
      #Check if the loss in cosine similarity between the original vs reconstructed after removing the signature is below the cutoff.
      delta = sim - sim_new
      if(delta <= max_delta){
        my_signatures = signatures_sel
        sim = sim_new
      }
      else{
        break
      }
  }
  
  #Combine data from the different samples
  sim_decay_fig = plot_sim_decay(sims, removed_sigs, max_delta)
  fit_res = fit_to_signatures(mut_mat_sample, my_signatures) #Perform final fit on selected signatures
  results = list("sim_decay_fig" = sim_decay_fig, "fit_res" = fit_res)
  all_results[[i]] = results
  }
  decay_figs = map(all_results, "sim_decay_fig")
  fit_res = map(all_results, "fit_res")
  contribution = map(fit_res, "contribution") %>% map(function(x) rownames_to_column(as.data.frame(x)))
  contribution = purrr::reduce(contribution, full_join, by = "rowname")
  #contribution[is.na(contribution)] = 0
  rownames(contribution) = contribution$rowname
  contribution = contribution %>% dplyr::select(-rowname) %>% as.matrix()
  
  reconstructed = map(fit_res, "reconstructed")
  reconstructed = do.call(cbind, reconstructed)
  
  fit_res = list("contribution" = contribution, "reconstructed" = reconstructed)
  results = list("sim_decay_fig" = decay_figs, "fit_res" = fit_res)
  
  return(results)
}


#Function to perform bootstrapping simulations on the signature refitting. The mutation matrix is resampled after which refitting is performed.
fit_to_signatures_bootstrapped = function(mut_mat, signatures, n_boots = 1000, max_delta = 0.05, method = c("selection", "selection_build", "mutpatterns_default", "mutpatterns_strict"), nr_simulated_sigs = 0, signatures_simulation_basis = F){
  method = match.arg(method)
  
  min_nr_muts = colSums(mut_mat) %>% min()
  if (min_nr_muts <= 10){
    warning("At least one of your samples has less than 10 mutations. This will negatively impact the signature refitting. Please consider removing or pooling this sample.")
  }
  
  # if (dim(mut_mat)[1] != 96) 
  #     stop(paste("Mutation count matrix input should have", 
  #                "dimensions 96 X n samples"))
  # if (dim(signatures)[1] != 96) 
  #     stop("Signatures input should have dimensions 96 X n signatures")
  
  n_samples = ncol(mut_mat)
  sig_names_tb = tibble("sigs" = colnames(signatures))
  if (nr_simulated_sigs != 0){
    sig_names_tb = rbind(sig_names_tb, tibble("sigs" = "SimulatedSigs"))
  }
  contri_list = vector("list", n_boots)
  for (i in 1:n_boots){
    mut_mat_resampled = resample_mut_mat(mut_mat)
    
    
    if (method == "selection"){
      refit_out = fit_to_signatures_selection(mut_mat_resampled, signatures, max_delta = max_delta, nr_simulated_sigs = nr_simulated_sigs, signatures_simulation_basis = signatures_simulation_basis)
      contri = refit_out$fit_res$contribution
    }
    else if (method == "selection_build"){
      fit_res = fit_to_signatures_regularized(mut_mat_resampled, signatures, regularized = T, min_delta = max_delta)
      contri = fit_res$contribution
    }
    else if (method == "mutpatterns_default"){
      fit_res = fit_to_signatures(mut_mat_resampled, signatures)
      contri = fit_res$contribution
    }
    else if (method == "mutpatterns_strict"){
      fit_res = fit_to_signatures(mut_mat_resampled, signatures)
      index = rowSums(fit_res$contribution >= 10) != 0 #Check whether a signature has at least 10 mutations in a single sample
      contri = fit_res$contribution[index,]
    }
    colnames(contri) = paste(colnames(contri), i, sep = "_")
    contri_tb = contri %>% as.data.frame() %>% rownames_to_column("sigs")
    simsigs = grepl("SimulatedSig", contri_tb$sigs)
    if (sum(simsigs) > 0){
      sum_simsigs = contri_tb %>% dplyr::filter(simsigs) %>% dplyr::select(-sigs) %>% dplyr::summarise_all(list(~sum)) %>% dplyr::mutate(sigs = "SimulatedSigs")
      contri_tb_other = contri_tb[!simsigs,]
      contri_tb = rbind(contri_tb_other, sum_simsigs)
    }
    contri_tb = left_join(sig_names_tb, contri_tb, by = "sigs") %>% dplyr::select(-sigs)
    contri_list[[i]] = contri_tb
    
    if (i%%10 == 0){
      print(paste0("Performed ", i, " of ", n_boots, " iterations"))
    }
  }
  contri_boots = do.call(cbind, contri_list) %>% t()
  colnames(contri_boots) = sig_names_tb$sigs
  return(contri_boots)
}


#Function to pool samples from a mutation matrix
pool_mut_mat = function(mut_mat, grouping){
  grouping = factor(grouping)
  mut_mat_group = mut_mat %>% t(.) %>% as_tibble() %>% dplyr::mutate(factor = grouping) %>% group_by(factor) %>% summarise_all(sum) %>% dplyr::select(-factor) %>% t(.)
  colnames(mut_mat_group) = levels(grouping)
  return(mut_mat_group)
}

#Functio to create and original vs reconstructed cosine similarity plot.
create_orivsrec_fig = function(mut_mat, signatures){
  cos_sim_ori_rec_full = cos_sim_matrix(mut_mat, signatures)
  cos_sim_ori_rec = as.data.frame(diag(cos_sim_ori_rec_full))
  colnames(cos_sim_ori_rec) = "cos_sim"
  cos_sim_ori_rec$sample = row.names(cos_sim_ori_rec_full)
  orivsrec_fig = ggplot(cos_sim_ori_rec, aes(y = cos_sim, x = factor(sample, levels = sample))) +
    geom_bar(stat = "identity", fill = "skyblue3") +
    coord_cartesian(ylim = c(0.6, 1), expand = F) +
    labs(y = "Cosine similarity\n original VS reconstructed", x = "") +
    theme_classic() +
    theme(text = element_text(size = 12), axis.text.x = element_text(angle = 90)) +
    geom_hline(aes(yintercept = 0.95))
  return(orivsrec_fig)
}

#Functions to remove indels from data
remove_indels_gr = function(gr, sample_name = "sample"){
  indel_f = width(gr$REF) != 1 | width(unlist(gr$ALT)) != 1
  nr_indel = sum(indel_f)
  gr = gr[!indel_f]
  if (nr_indel >= 1){
    print(paste0("Removed ", nr_indel, " indels from sample: ", sample_name))
  }
  return(gr)
}


remove_indels = function(grl){
  if (class(grl)[[1]] == "CompressedGRangesList"){
    gr_list = mapply(remove_indels_gr, grl, names(grl))
    grl = GRangesList(gr_list)
    return(grl)
  } else if (class(grl)[[1]] == "GRanges"){
    gr = remove_indels_gr(grl)
    return(gr)
  }
}

