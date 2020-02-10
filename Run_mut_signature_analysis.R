#Script to perform mutational signature analysis.
#Please copy this script and then change the paths, to the data you wish to analyze
#You will need to change some settings based on your data / research question.
#Settings that may need to be changed are indicated in the comments with: CHANGE
#v0.12.0

# -._    _.--'"`'--._    _.--'"`'--._    _.--'"`'--._    
# '-:`.'|`|"':-.  '-:`.'|`|"':-.  '-:`.'|`|"':-.  '.` :   
#   '.  '.  | |  | |'.  '.  | |  | |'.  '.  | |  | |'.  
#   : '.  '.| |  | |  '.  '.| |  | |  '.  '.| |  | |   
# '   '.  `.:_ | :_.' '.  `.:_ | :_.' '.  `.:_ | :_.'  
#          `-..,..-'       `-..,..-'       `-..,..-'      

# Load libraries ----------------------------------------------------------
library(tidyverse)
library(VariantAnnotation)
library(MutationalPatterns)
library(GenomicRanges)
library(magrittr)
library(pracma)
library(vegan)
library(ggdendro)
library(cowplot)
source("~/surfdrive/Shared/Boxtel_General/Scripts/pmc_vanboxtel/Freek_misc/Run_mut_signature_analysis/Better_refitting.R")


# Reading in data ---------------------------------------------------------
#Output directory. Please choose a directory where you want to store the results of your analysis.
dir = "~/surfdrive//Shared/Projects/Freek/"

#Creates a separate output directory which will be placed in the directory you choose.
dir_out = file.path(dir, "mut_signature_analysis")
if (!dir.exists(dir_out)){
    dir.create(dir_out)
}
setwd(dir_out)

#Choose your genome version here
library(BSgenome.Hsapiens.UCSC.hg19)
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"

#The basic name of your signatures. Change this to ID for indels and DBS for double base substitutions.
base_sig_name = "SBS"

#Read in your vcf files. Please provide your own vcf files and sample names here. If you already loaded your vcfs you can skip this step.
vcfs = list.files("~/surfdrive/Shared/Boxtel_General/Data/Mutation_data/SNVs/Mirjam_HSCT/", full.names = T)
vcfs = vcfs[grepl("_CALLABLE.vcf", vcfs)]
samples = gsub(".*/(.*)_.*_Q50.*", "\\1", vcfs)
grl = read_vcfs_as_granges(vcfs, samples, genome = ref_genome)
grl = remove_indels(grl) #Remove any indels you might have in your data.

#Read in the existing signatures.
cosmic_sig_fname = "~/surfdrive/Shared/Boxtel_General/Data/Sigs/SBS/sigProfiler_SBS_working_signatures_incl_hspc.txt"
signatures = read.table(cosmic_sig_fname, sep = "\t", header = T)
signatures = as.matrix(signatures[,-c(1,2)])

signatures_serena = read.table("~/surfdrive/Shared/Boxtel_General/Data/Sigs/SBS/Mutagen53_sub_signature.txt", header = T, sep = "\t")
signatures_serena = as.matrix(signatures_serena[,-c(1,2)])

#Read in the existing mutation matrix for nmf
mut_matrix_existing_fname = "~/surfdrive/Shared/Boxtel_General/Data/MutationMatrix/NMF_matrix.txt"
mut_mat_existing = read.table(mut_matrix_existing_fname, sep = "\t", header = T)
write.table(mut_mat_existing, "mut_matrix_existing.txt", sep = "\t", quote = F, col.names = T, row.names = T) #The existing data matrix is also written out, so you can easily look back which one was used.


# Determine mutation contexts ---------------------------------------------
#Create mutation matrix of your data (mut_matrix2 is a faster version of the MutationalPatterns mut_matrix)
mut_mat = mut_matrix2(grl, ref_genome)
write.table(mut_mat, "mut_matrix.txt", sep = "\t", quote = F, col.names = T, row.names = T)
grouping = F

#OPTIONAL: Pool mutations that are in the same group
grouping = c(rep("group1", 11), rep("group2", 11))
mut_mat = pool_mut_mat(mut_mat, grouping)
write.table(mut_mat, "pooled_mut_matrix.txt", sep = "\t", quote = F, col.names = T, row.names = T)

# NMF ---------------------------------------------------------------------
#Add existing data for NMF.
mut_mat_nmf = cbind(mut_mat, mut_mat_existing) + 0.0001

#Create plots showing how many signatures you should generate with nmf. This takes a long time.
estimate = nmf(mut_mat_nmf, rank=2:7, method="brunet", nrun=10, seed=123456)
estimate_fig = plot(estimate)
print(estimate_fig)

pdf("nmf_estimate.pdf", useDingbats = F)
estimate_fig
dev.off()

#The rank is the number of de novo signatures that nmf will generate. Set it based on the estimate plots.
#If there is an 'elbow' in these plots, this is what your rank should be. It can be usefull to try out multiple values and see what works best.
#The mut_mat_existing data will generally result in the inclusion of signatures: SBS1, SBS5, SBS18 and HSPC.
#Using rank 5 will include a 5th often novel signature.
###CHANGE this based on the plots. DO NOT just use the default value###
rank_nmf = 5
nmf_res = extract_signatures(mut_mat_nmf, rank = rank_nmf, nrun = 10)



# Combine NMF signatures with existing ones  ---------------------------------
#If any of the signatures generated with NMF is already known, this signature will receive the known name
rownames(nmf_res$contribution) = 1:nrow(nmf_res$contribution)
colnames(nmf_res$signatures) = 1:ncol(nmf_res$signatures)
sim_matrix = cos_sim_matrix(signatures, nmf_res$signatures)

#Cutof at which NMF and existing signature are considered identical. Can be CHANGED if needed
cutof_same_sig = 0.85

j = 0
for (i in 1:ncol(sim_matrix)){
    sig_sim = sim_matrix[,i]
    cossim = max(sig_sim)
    if (cossim > cutof_same_sig){
        row = which.max(sig_sim)
        sig = names(sig_sim)[row]
        rownames(nmf_res$contribution)[i] = sig
        colnames(nmf_res$signatures)[i] = sig
    } else{
        j = j + 1
        rownames(nmf_res$contribution)[i] = paste0(base_sig_name, LETTERS[j])
        colnames(nmf_res$signatures)[i] = paste0(base_sig_name, LETTERS[j])
    }
}
dupli_sigs = colnames(nmf_res$signatures) %>% 
    duplicated() %>% 
    any()
if (dupli_sigs){
    warning("You have multiple NMF signatures that are linked to the same existing signature.\n
            Please use a lower rank in the NMF or increase the cutoff at which a NMF and \n
            existing signature are considered identical.")
}
# Plot figures of the NMF results -----------------------------------------
#Similarity of nmf signatures with existing signatures
sig_cosim_heat_fig = plot_sim_with_signatures(nmf_res$signatures, signatures)

#Contribtution of signatures according to nmf
pc1 = plot_contribution(nmf_res$contribution, nmf_res$signature, mode = "relative") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
pc2 = plot_contribution(nmf_res$contribution, nmf_res$signature, mode = "absolute") +
    theme(axis.text.x = element_text(angle = 90))

contri_fig = plot_grid(pc1, pc2, align = "v", nrow = 2)

contri_heat_fig = plot_contribution_heatmap(nmf_res$contribution, cluster_samples = F) +
    theme_classic() +
    theme(text = element_text(size = 6))

#96 profiles of nmf signatures
nmf_sigs_fig = plot_96_profile(nmf_res$signatures, condensed = T)

#Similarity between reconstruction and original profiles.
orivsrec_fig = create_orivsrec_fig(mut_mat_nmf, nmf_res$reconstructed)

pdf("nmf_results.pdf", useDingbats = F)
print(nmf_sigs_fig)
print(sig_cosim_heat_fig)
print(plot_grid(pc1, pc2, align = "v", nrow = 2))
print(contri_heat_fig)
print(orivsrec_fig)
dev.off()

# Plot the similarity of the mut_mat with existing signatures -----------------------------
sig_cosim_heat_fig = plot_sim_with_signatures(mut_mat, signatures)
serena_sig_cosim_heat_fig = plot_sim_with_signatures(mut_mat, signatures_serena)
inner_cosheat_fig = plot_inner_cosheat(mut_mat)

pdf("similarity_signatures.pdf", useDingbats = F, width = 10)
print(sig_cosim_heat_fig)
print(serena_sig_cosim_heat_fig)
print(inner_cosheat_fig)
dev.off()

# Signature refitting -----------------------------------------------------
#Get signatures that will be used for refitting. (If a signature similar to SBS1 is found in the NMF, then SBS1 is used not the similar signature)
nmf_signatures = nmf_res$signatures
new_signatures = nmf_signatures[,!colnames(nmf_signatures) %in% colnames(signatures), drop = F]
new_signatures = prop.table(new_signatures, 2)

if (ncol(new_signatures) > 0){
    new_sigs_fig = plot_96_sigs(new_signatures)
    pdf("new_signatures.pdf", useDingbats = F)
    print(new_sigs_fig)
    dev.off()
    
    write.table(new_signatures, "new_signatures.txt", sep = "\t", quote = F, col.names = T, row.names = F)
}

old_signatures = signatures[,colnames(signatures) %in% colnames(nmf_signatures), drop = F]
signatures = cbind(old_signatures, new_signatures)

#Increasing the max_delta will increase the strictness of the signature refitting, resulting in less signatures. (and vice versa).
#CHANGE this if you want to change the strictness of the signature refitting
max_delta = 0.05

#Fit to the signatures.
refit_out = fit_to_signatures_selection(mut_mat, signatures, max_delta = max_delta)
sim_decay_figs = refit_out$sim_decay_fig #Show how the cosine similarity is decreased as signatures are removed
contri = refit_out$fit_res$contribution
contri[is.na(contri)] = 0
refit_heatmap_fig = plot_contribution_heatmap(contri, cluster_samples = F) #Show a contribution heatmap of the refitting.
pc1 = plot_contribution_nonmf(contri, mode = "relative") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
pc2 = plot_contribution_nonmf(contri, mode = "absolute")
contri_fig = plot_grid(pc1, pc2, align = "v", nrow = 2)

orivsrec_fig = create_orivsrec_fig(mut_mat, refit_out$fit_res$reconstructed) #Plot the cosime similarity of the original vs the reconstructed profile.

pdf("signature_refitting.pdf", useDingbats = F)
print(contri_fig)
print(refit_heatmap_fig)
print(orivsrec_fig)
print(sim_decay_figs)
dev.off()

#Fit to the signatures with bootstrapping
contri_boots = fit_to_signatures_bootstrapped(mut_mat, signatures, method = "selection", max_delta = max_delta)
nr_selected_fig = plot_nr_signatures_selected(contri_boots) #Plot how many signatures were selected during the different bootstrapping iterations.
fraction_fig = plot_fraction_contri(contri_boots) #Plot how often each signature was selected.
contri_boots_fig = plot_boots_contri(contri_boots) #Plot the contribution of bootstrapped selected signatures
contri_boots_fig2 = plot_boots_contri(contri_boots, mode = "relative") #Plot the contribution of bootstrapped selected signatures
ori_vs_rec_fig = plot_cosine_bootstrapped(mut_mat, contri_boots, signatures) #Plot the cosime similarity of the original vs the reconstructed profile. The reconstructed profile is based on the mean bootstrapped signature contributions
sig_cor_figs = plot_sig_contri_cor(contri_boots) #Plot the signature contribution correlations between the signatures.

pdf("bootstrapped_signature_refitting.pdf", useDingbats = F)
print(nr_selected_fig)
print(fraction_fig)
print(contri_boots_fig)
print(contri_boots_fig2)
print(ori_vs_rec_fig)
print(sig_cor_figs)
dev.off()

#Generate a text file containing the settings you have used. This way you can always lookup how you did your analysis.
settings_tb = tibble("Settings" = c("ref_genome", "base_sig_name", "sig_file", "mut_matrix_existing", "vcfs", "samples", "grouping", "rank_nmf", "max_delta", "cutof_same_sig"), 
       "Used" = c(ref_genome, base_sig_name, cosmic_sig_fname, mut_matrix_existing_fname, str_c(vcfs, collapse = ";"), str_c(samples, collapse = ";"), str_c(grouping, collapse = ";"), rank_nmf, max_delta, cutof_same_sig))
write_tsv(settings_tb, "settings.txt")
