# Packages
library(tximport)
library(tximeta)
library(GenomicFeatures)
library(data.table)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(org.Hs.eg.db)
library(biomaRt)
library(plyranges)
library(readr)
library(ggtext)
library(formattable)
library(flextable)
library(officer)
library(readxl)
library(openxlsx)


source("/work/users/s/e/seyoun/colabs/LCSP/seq_pipelines/scripts/rna_workflow/utils/diff_rna.r")

# Data
base_dir <- "/work/users/s/e/seyoun/colabs/LCSP/seq_pipelines"
rna_dir <- file.path(base_dir,"rna_output")
differential_dir <- file.path(rna_dir,"00_differential_expression")
diff_data_dir <- file.path(differential_dir,"data")
diff_plot_dir <- file.path(differential_dir,"plots")

if(!dir.exists(differential_dir)){dir.create(differential_dir)}
if(!dir.exists(diff_data_dir)){dir.create(diff_data_dir)}
if(!dir.exists(diff_plot_dir)){dir.create(diff_plot_dir)}

# Normoxia ----------------------------------------------------------------------------------------------------
sampleInfo <- fread(file.path(base_dir,"samplesheet.txt"))
DonorInfo <- fread(file.path(base_dir,"donorinfo.txt"))
donorinfo_sub <- DonorInfo |> dplyr::select(Donor, Sex, Age,Date_of_collection)

sampleInfo_sub <- sampleInfo |> mutate(names = paste0(Proj, "_", Donor, "_", Condition, "_", Treatment)) |>
                    dplyr::select(names,Donor, Condition,Treatment,RIN,RNA_Isolation_Date) |> 
                    left_join(donorinfo_sub, by="Donor")

sampleNMs <- sampleInfo_sub$names

sample_norm  <- sampleInfo_sub |> filter(Condition == "norm")
sampleNMs_norm <- sample_norm$names


# Add quant paths and names
coldata_norm <- sample_norm
coldata_norm$files <- file.path(rna_dir,"quant",sampleNMs_norm,"quant.sf")
file.exists(coldata_norm$files)

## Import data with tximeta & summarize to gene
se_norm <- tximeta(coldata_norm)
gse_norm <- summarizeToGene(se_norm)

# Gene level counts
colData(gse_norm)[] <- lapply(colData(gse_norm), factor)
colData(gse_norm)$Treatment <- factor(colData(gse_norm)$Treatment, levels= c("PBS","DOXO"))
dds_norm <- DESeqDataSet(gse_norm, design = ~ Donor + Treatment)

## Filter out lowely expressed genes 10 for at least 40% of total samples
keep <- rowSums(counts(dds_norm) >= 10) >= ceiling(nrow(colData(gse_norm))* 0.4)
dds_norm_sub <- dds_norm[keep,]

# Fit model
dds_norm_fit <- DESeq(dds_norm_sub)

# Get results
resRNA_norm <- results(dds_norm_fit)
#resultsNames(dds_norm_fit)
res_Shrink_norm_gr <- lfcShrink(dds_norm_fit, coef='Treatment_DOXO_vs_PBS', type="apeglm", format =  "GRanges") |>
  plyranges::names_to_column("gene_id")

save(res_Shrink_norm_gr,dds_norm_fit,resRNA_norm, sample_norm, file = file.path(diff_data_dir,"deseq_results_norm_tximeta.Rdata"))

load(file = file.path(diff_data_dir,"deseq_results_norm_tximeta.Rdata")) # res_Shrink_norm_gr,dds_norm_fit,resRNA_norm, sample_norm

# PCA plot for normoxia
res_Shrink_norm_gr <- res_Shrink_norm_gr[!is.na(res_Shrink_norm_gr$padj),]

norm_dds_norm <- normTransform(dds_norm_fit)
rlog_dds_norm <- rlog(dds_norm_fit)
vst_dds_norm <- vst(dds_norm_fit)

save(norm_dds_norm, rlog_dds_norm, vst_dds_norm, file = file.path(diff_data_dir,"deseq_transformed_norm_tximeta.Rdata"))

load(file = file.path(diff_data_dir,"deseq_transformed_norm_tximeta.Rdata")) # norm_dds_norm, rlog_dds_norm, vst_dds_norm

# plotPCA(vst_dds_norm, c("Treatment","Sex"),ntop=18973)

pca_vst_re <- prcomp(t(assay(vst_dds_norm)))
pca_df <- as.data.frame(pca_vst_re$x)
variance_explained <- summary(pca_vst_re)$importance[2, ]
pca_df$Donor <- colData(vst_dds_norm)$Donor
pca_df$Treatment <-colData(vst_dds_norm)$Treatment
peak_count <- nrow(assay(vst_dds_norm))

normPCAplot <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Donor, color = Treatment)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("PBS" = "#0C7BDC", "DOXO" = "#FFC20A")) +
  labs(
    title = paste0("Total peaks (n=", peak_count, ") with all samples"),
    x = paste0("PC1 Variance: ", round(variance_explained[1] * 100, 2), "%"),
    y = paste0("PC2 Variance: ", round(variance_explained[2] * 100, 2), "%")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill='transparent'))

save(normPCAplot, file= file.path(diff_plot_dir,"PCA_normoxia_deseq.rda"))

# hypoxia ----------------------------------------------------------------------------------------------------

sample_hypo  <- sampleInfo_sub |> filter(Condition == "hypo")
sampleNMs_hypo <- sample_hypo$names


# Add quant paths and names
coldata_hypo <- sample_hypo
coldata_hypo$files <- file.path(rna_dir,"quant",sampleNMs_hypo,"quant.sf")
file.exists(coldata_hypo$files)

## Import data with tximeta & summarize to gene
se_hypo <- tximeta(coldata_hypo)
gse_hypo <- summarizeToGene(se_hypo)

# Gene level counts
colData(gse_hypo)[] <- lapply(colData(gse_hypo), factor)
colData(gse_hypo)$Treatment <- factor(colData(gse_hypo)$Treatment, levels= c("PBS","DOXO"))
dds_hypo <- DESeqDataSet(gse_hypo, design = ~ Donor + Treatment)
## Filter out lowely expressed genes 10 for at least 40% of total samples
keep <- rowSums(counts(dds_hypo) >= 10) >= ceiling(nrow(colData(gse_hypo))* 0.4)
dds_hypo_sub <- dds_hypo[keep,]

# Fit model
dds_hypo_fit <- DESeq(dds_hypo_sub)

# Get results
resRNA_hypo <- results(dds_hypo_fit)
#resultsNames(dds_hypo_fit)
res_Shrink_hypo_gr <- lfcShrink(dds_hypo_fit, coef='Treatment_DOXO_vs_PBS', type="apeglm", format =  "GRanges") |>
  plyranges::names_to_column("gene_id")

save(res_Shrink_hypo_gr,dds_hypo_fit,resRNA_hypo, sample_hypo, file = file.path(diff_data_dir,"deseq_results_hypo_tximeta.Rdata"))

# PCA plot for normoxia
res_Shrink_hypo_gr <- res_Shrink_hypo_gr[!is.na(res_Shrink_hypo_gr$padj),]

norm_dds_hypo <- normTransform(dds_hypo_fit)
rlog_dds_hypo <- rlog(dds_hypo_fit)
vst_dds_hypo <- vst(dds_hypo_fit)
save(norm_dds_hypo, rlog_dds_hypo, vst_dds_hypo, file = file.path(diff_data_dir,"deseq_transformed_hypo_tximeta.Rdata"))
load(file = file.path(diff_data_dir,"deseq_transformed_hypo_tximeta.Rdata")) #norm_dds_hypo, rlog_dds_hypo, vst_dds_hypo
# plotPCA(vst_dds_hypo, c("Treatment","Sex"),ntop=19232)

pca_vst_re <- prcomp(t(assay(vst_dds_hypo)))
pca_df <- as.data.frame(pca_vst_re$x)
variance_explained <- summary(pca_vst_re)$importance[2, ]
pca_df$Donor <- colData(vst_dds_hypo)$Donor
pca_df$Treatment <-colData(vst_dds_hypo)$Treatment
peak_count <- nrow(assay(vst_dds_hypo))

hypoPCAplot <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Donor, color = Treatment)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("PBS" = "#0C7BDC", "DOXO" = "#FFC20A")) +
  labs(
    title = paste0("Total peaks (n=", peak_count, ") with all samples"),
    x = paste0("PC1 Variance: ", round(variance_explained[1] * 100, 2), "%"),
    y = paste0("PC2 Variance: ", round(variance_explained[2] * 100, 2), "%")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill='transparent'))

save(hypoPCAplot, file= file.path(diff_plot_dir,"PCA_hypoxia_deseq.rda"))


#--------------------------------------------------------------------
# Differential gene expression 
# Normoxia
load(file = file.path(diff_data_dir,"deseq_results_norm_tximeta.Rdata")) # res_Shrink_norm_gr,dds_norm_fit,resRNA_norm, sample_norm

diff_normoxia_sig_p05log2fc1 <- res_Shrink_norm_gr |> plyranges::filter(abs(log2FoldChange) > 1 & padj < 0.05) # 1281
diff_normoxia_sig_p01log2fc1 <- res_Shrink_norm_gr |>  plyranges::filter(abs(log2FoldChange) > 1 & padj < 0.01) # 1153
diff_normoxia_sig_p05log2fc15 <- res_Shrink_norm_gr |> plyranges::filter(abs(log2FoldChange) > 1.5 & padj < 0.05) # 632


diff_normoxia_all  <- diff_normoxia_sig_p05log2fc1
#finiding up sig and down sig
up_norm <- diff_normoxia_all[diff_normoxia_all$log2FoldChange > 0 & diff_normoxia_all$padj < 0.05,]
down_norm <- diff_normoxia_all[diff_normoxia_all$log2FoldChange < 0 & diff_normoxia_all$padj < 0.05,]
static_norm <- res_Shrink_norm_gr[!(res_Shrink_norm_gr$gene_id %in% c(up_norm$gene_id, down_norm$gene_id))]


# Making all the background genes and adding the class
res_Shrink_norm_df <- as.data.frame(res_Shrink_norm_gr)
res_Shrink_norm_df$class <- "static"
res_Shrink_norm_df$class[res_Shrink_norm_df$gene_id %in% up_norm$gene_id] <- "gained"
res_Shrink_norm_df$class[res_Shrink_norm_df$gene_id %in% down_norm$gene_id] <- "lost"

res_Shrink_norm_df_hgnc <- res_Shrink_norm_df |> add_symbol(gene_annot_clean)
res_Shrink_norm_df_hgnc |> filter(grepl("ENSG", SYMBOL, ignore.case = TRUE)) |> dim() # 2912 doens't have symbol

# # Adding the HGNC symbol
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# res_Shrink_norm_df_cleanVer <- res_Shrink_norm_df |> dplyr::mutate(gene_id = gsub("\\..*$", "", gene_id))

# annots <- AnnotationDbi::select(org.Hs.eg.db, res_Shrink_norm_df_cleanVer$gene_id,
#                                 columns=c("SYMBOL"), keytype="ENSEMBL")
# annots <- dplyr::rename(annots, gene_id = ENSEMBL)
# annots_unique <- annots %>% distinct(gene_id, .keep_all = TRUE)
# res_Shrink_norm_df_hgnc <- left_join(res_Shrink_norm_df_cleanVer, annots_unique, by = "gene_id") # 2865 of them couldn't find the HGNC
save(res_Shrink_norm_df_hgnc, file= file.path(diff_data_dir,"shrink_normoxia_hgnc_tximeta.Rdata"))
write.table(res_Shrink_norm_df_hgnc, file= file.path(diff_data_dir,"deseq_normoxia_results_tximeta.tsv"),
            sep="\t",quote=F,col.names=TRUE, row.names = FALSE)

#  Save to run GO and pathway analysis
# up genes
res_Shrink_norm_df_hgnc |> filter(class == "gained") |> dplyr::mutate(gene_id = gsub("\\..*$", "", gene_id)) |>
    dplyr::select(gene_id) |> 
    distinct() |> 
    write.table(file.path(diff_data_dir,"diff_normoxia_upregulated_padj05log2fc1_genes.txt"),
    sep="\t",quote=F,col.names=FALSE, row.names = FALSE)

# down genes
res_Shrink_norm_df_hgnc |> filter(class == "lost")  |> dplyr::mutate(gene_id = gsub("\\..*$", "", gene_id)) |>
    dplyr::select(gene_id) |> 
    distinct() |> 
    write.table(file.path(diff_data_dir,"diff_normoxia_downregulated_padj05log2fc1_genes.txt"),
    sep="\t",quote=F,col.names=FALSE, row.names = FALSE)

#background genes
res_Shrink_norm_df_hgnc  |> dplyr::mutate(gene_id = gsub("\\..*$", "", gene_id)) |>
    dplyr::select(gene_id) |> 
    distinct() |> 
    write.table(file.path(diff_data_dir,"diff_normoxia_background_padj05log2fc1_genes.txt"),
    sep="\t",quote=F,col.names=FALSE, row.names = FALSE)

homer_rna <- file.path(base_dir,"scripts/rna_workflow/01_DEG/run_homer.sh")

up_norm_cmd <- paste(
  homer_rna,
  file.path(diff_data_dir, "diff_normoxia_upregulated_padj05log2fc1_genes.txt"),
  file.path(diff_data_dir, "diff_normoxia_background_padj05log2fc1_genes.txt"),
  file.path(diff_data_dir, "deseq_homer", "normoxia_up_LFC1_padj05")
)

down_norm_cmd <- paste(
  homer_rna,
  file.path(diff_data_dir, "diff_normoxia_downregulated_padj05log2fc1_genes.txt"),
  file.path(diff_data_dir, "diff_normoxia_background_padj05log2fc1_genes.txt"),
  file.path(diff_data_dir, "deseq_homer", "normoxia_down_LFC1_padj05")
)

#system(up_norm_cmd)
#system(down_norm_cmd)


# Hypoxia ---------
load(file = file.path(diff_data_dir,"deseq_results_hypo_tximeta.Rdata")) #res_Shrink_hypo_gr,dds_hypo_fit,resRNA_hypo, sample_hypo

diff_hypoxia_sig_p05log2fc1 <- res_Shrink_hypo_gr |> plyranges::filter(abs(log2FoldChange) > 1 & padj < 0.05) # 593
diff_hypoxia_sig_p01log2fc1 <- res_Shrink_hypo_gr |> plyranges::filter(abs(log2FoldChange) > 1 & padj < 0.01) # 525
diff_hypoxia_sig_p05log2fc15 <- res_Shrink_hypo_gr |> plyranges::filter(abs(log2FoldChange) > 1.5 & padj < 0.05) # 299


diff_hypoxia_all  <- diff_hypoxia_sig_p05log2fc1
#finiding up sig and down sig
up_hypo <- diff_hypoxia_all[diff_hypoxia_all$log2FoldChange > 0 & diff_hypoxia_all$padj < 0.05,] # 316
down_hypo <- diff_hypoxia_all[diff_hypoxia_all$log2FoldChange < 0 & diff_hypoxia_all$padj < 0.05,] # 279
static_hypo <- res_Shrink_hypo_gr[!(res_Shrink_hypo_gr$gene_id %in% c(up_hypo$gene_id, down_hypo$gene_id))] # 18637


# Making all the background genes and adding the class
res_Shrink_hypo_df <- as.data.frame(res_Shrink_hypo_gr)
res_Shrink_hypo_df$class <- "static"
res_Shrink_hypo_df$class[res_Shrink_hypo_df$gene_id %in% up_hypo$gene_id] <- "gained"
res_Shrink_hypo_df$class[res_Shrink_hypo_df$gene_id %in% down_hypo$gene_id] <- "lost"


res_Shrink_hypo_df_hgnc <- res_Shrink_hypo_df |> add_symbol(gene_annot_clean)
res_Shrink_hypo_df_hgnc |> filter(grepl("ENSG", SYMBOL, ignore.case = TRUE)) |> dim() # 3078 doens't have symbol


# Adding the HGNC symbol
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# res_Shrink_hypo_df_cleanVer <- res_Shrink_hypo_df |> dplyr::mutate(gene_id = gsub("\\..*$", "", gene_id))
# annots <- AnnotationDbi::select(org.Hs.eg.db, res_Shrink_hypo_df_cleanVer$gene_id,
#                                 columns=c("SYMBOL"), keytype="ENSEMBL")
# annots <- dplyr::rename(annots, gene_id = ENSEMBL)
# annots_unique <- annots %>% distinct(gene_id, .keep_all = TRUE)
# res_Shrink_hypo_df_hgnc <- left_join(res_Shrink_hypo_df_cleanVer, annots_unique, by = "gene_id") # 3027 of them couldn't find the HGNC
save(res_Shrink_hypo_df_hgnc, file= file.path(diff_data_dir,"shrink_hypoxia_hgnc_tximeta.Rdata"))

write.table(res_Shrink_hypo_df_hgnc, file= file.path(diff_data_dir,"deseq_hypoxia_results_tximeta.tsv"),
            sep="\t",quote=F,col.names=TRUE, row.names = FALSE)

#  Save to run GO and pathway analysis
# up genes
res_Shrink_hypo_df_hgnc |> filter(class == "gained") |> dplyr::mutate(gene_id = gsub("\\..*$", "", gene_id)) |>
    dplyr::select(gene_id) |> 
    distinct() |> 
    write.table(file.path(diff_data_dir,"diff_hypoxia_upregulated_padj05log2fc1_genes.txt"),
    sep="\t",quote=F,col.names=FALSE, row.names = FALSE)

# down genes
res_Shrink_hypo_df_hgnc |> filter(class == "lost") |> dplyr::mutate(gene_id = gsub("\\..*$", "", gene_id)) |>
    dplyr::select(gene_id) |> 
    distinct() |> 
    write.table(file.path(diff_data_dir,"diff_hypoxia_downregulated_padj05log2fc1_genes.txt"),
    sep="\t",quote=F,col.names=FALSE, row.names = FALSE)

#background genes
res_Shrink_hypo_df_hgnc |> dplyr::mutate(gene_id = gsub("\\..*$", "", gene_id)) |>
    dplyr::select(gene_id) |>  
    distinct() |> 
    write.table(file.path(diff_data_dir,"diff_hypoxia_background_padj05log2fc1_genes.txt"),
    sep="\t",quote=F,col.names=FALSE, row.names = FALSE)

homer_rna <- file.path(base_dir,"scripts/rna_workflow/01_DEG/run_homer.sh")

up_hypo_cmd <- paste(
  homer_rna,
  file.path(diff_data_dir, "diff_hypoxia_upregulated_padj05log2fc1_genes.txt"),
  file.path(diff_data_dir, "diff_hypoxia_background_padj05log2fc1_genes.txt"),
  file.path(diff_data_dir, "deseq_homer", "hypoxia_up_LFC1_padj05")
)

down_hypo_cmd <- paste(
  homer_rna,
  file.path(diff_data_dir, "diff_hypoxia_downregulated_padj05log2fc1_genes.txt"),
  file.path(diff_data_dir, "diff_hypoxia_background_padj05log2fc1_genes.txt"),
  file.path(diff_data_dir, "deseq_homer", "hypoxia_down_LFC1_padj05")
)

# system(up_hypo_cmd)
# system(down_hypo_cmd)


#------------------------------------------------------------------------------------------------------------------------------------
# Interaction between condition:treatment 
#------------------------------------------------------------------------------------------------------------------------------------
sampleInfo <- fread(file.path(base_dir,"samplesheet.txt"))
DonorInfo <- fread(file.path(base_dir,"donorinfo.txt"))
donorinfo_sub <- DonorInfo |> dplyr::select(Donor, Sex, Age,Date_of_collection)

sampleInfo_sub <- sampleInfo |> mutate(names = paste0(Proj, "_", Donor, "_", Condition, "_", Treatment)) |>
                    dplyr::select(names,Donor, Condition,Treatment,RIN,RNA_Isolation_Date) |> 
                    left_join(donorinfo_sub, by="Donor")

sampleNMs <- sampleInfo_sub$names

# Add quant paths and names
coldata <- sampleInfo_sub
coldata$files <- file.path(rna_dir,"quant",sampleNMs,"quant.sf")
file.exists(coldata$files)

## Import data with tximeta & summarize to gene
se <- tximeta(coldata)
gse <- summarizeToGene(se)

# Gene level counts
colData(gse)[] <- lapply(colData(gse), factor)
colData(gse)$Condition <- factor(colData(gse)$Condition, levels= c("norm","hypo"))
colData(gse)$Treatment <- factor(colData(gse)$Treatment, levels= c("PBS","DOXO"))
dds <- DESeqDataSet(gse, design = ~ Donor + Condition + Treatment + Condition:Treatment)

## Filter out lowely expressed genes 10 for at least 40% of total samples
keep <- rowSums(counts(dds) >= 10) >= ceiling(nrow(colData(gse))* 0.4)
dds_sub <- dds[keep,]

# Fit model
dds_sub_fit <- DESeq(dds_sub)

# Get results
resRNA_inter <- results(dds_sub_fit)

#resultsNames(dds_sub_fit)
res_Shrink_gr <- lfcShrink(dds_sub_fit, coef='Conditionhypo.TreatmentDOXO', type="apeglm", format =  "GRanges") |>
  plyranges::names_to_column("gene_id")

save(res_Shrink_gr,dds_sub_fit,resRNA_inter, sampleInfo_sub, file = file.path(diff_data_dir,"deseq_results_interaction_tximeta.Rdata"))
load(file = file.path(diff_data_dir,"deseq_results_interaction_tximeta.Rdata")) # res_Shrink_gr,dds_sub_fit,resRNA_inter, sampleInfo_sub

# PCA plot for normoxia
res_Shrink_inter_gr <- res_Shrink_gr[!is.na(res_Shrink_gr$padj),]

norm_dds_inter <- normTransform(dds_sub_fit)
rlog_dds_inter <- rlog(dds_sub_fit)
vst_dds_inter <- vst(dds_sub_fit)

save(norm_dds_inter, rlog_dds_inter, vst_dds_inter, file = file.path(diff_data_dir,"deseq_transformed_interaction_tximeta.Rdata"))
load(file = file.path(diff_data_dir,"deseq_transformed_interaction_tximeta.Rdata")) #norm_dds_inter, rlog_dds_inter, vst_dds_inter

# plotPCA(vst_dds_inter, c("Donor"),ntop=18868) # Donor 589 might be slightly outlier?

pca_vst_re <- prcomp(t(assay(vst_dds_inter)))
pca_df <- as.data.frame(pca_vst_re$x)
variance_explained <- summary(pca_vst_re)$importance[2, ]
pca_df$Donor <- colData(vst_dds_inter)$Donor
pca_df$Treatment <- colData(vst_dds_inter)$Treatment
pca_df$Condition <- colData(vst_dds_inter)$Condition
gene_count <- nrow(assay(vst_dds_inter))

interactionPCAplot <- ggplot(pca_df, aes(x=PC1, y=PC2, label=Donor, color=Treatment,shape = Condition)) +
  geom_point(size=2) +
  scale_color_manual(values = c("PBS" = "#0C7BDC", "DOXO" = "#FFC20A")) +
  scale_shape_manual(values = c("norm" = 16, "hypo" = 17)) +
  labs(
    title = paste0("Total genes (n=", gene_count, ") with all samples"),
    x = paste0("PC1 Variance: ", round(variance_explained[1]*100, 2), "%"),
    y = paste0("PC2 Variance: ", round(variance_explained[2]*100, 2), "%")
  ) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        plot.background = element_rect(fill='transparent'))

save(interactionPCAplot, file = file.path(diff_plot_dir, "PCA_interaction_deseq_tximeta.rda"))

# Interaction finding the difference ---------
load(file = file.path(diff_data_dir,"deseq_results_interaction_tximeta.Rdata")) #res_Shrink_gr,dds_sub_fit,resRNA_inter, sampleInfo_sub
load(file = file.path(diff_data_dir,"deseq_transformed_interaction_tximeta.Rdata")) #norm_dds_inter, rlog_dds_inter, vst_dds_inter

diff_inter_sig_p05log2fc1 <- res_Shrink_gr |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.05) # 294
diff_inter_sig_p01log2fc1 <- res_Shrink_gr |> dplyr::filter(abs(log2FoldChange) > 1 & padj < 0.01) # 199
diff_inter_sig_p05log2fc15 <- res_Shrink_gr |> dplyr::filter(abs(log2FoldChange) > 1.5 & padj < 0.05) # 133


diff_inter_all  <- diff_inter_sig_p05log2fc1
#finiding up sig and down sig
up_inter <- diff_inter_all[diff_inter_all$log2FoldChange > 0 & diff_inter_all$padj < 0.05,] # 151 
down_inter <- diff_inter_all[diff_inter_all$log2FoldChange < 0 & diff_inter_all$padj < 0.05,] # 143
static_inter <- res_Shrink_gr[!(res_Shrink_gr$gene_id %in% c(up_inter$gene_id, down_inter$gene_id))] # 18574

# Making all the background genes and adding the class
res_Shrink_inter_df <- as.data.frame(res_Shrink_gr)
res_Shrink_inter_df$class <- "static"
res_Shrink_inter_df$class[res_Shrink_inter_df$gene_id %in% up_inter$gene_id] <- "gained"
res_Shrink_inter_df$class[res_Shrink_inter_df$gene_id %in% down_inter$gene_id] <- "lost"

res_Shrink_inter_df_hgnc <- res_Shrink_inter_df |> add_symbol(gene_annot_clean)
res_Shrink_inter_df_hgnc |> filter(grepl("ENSG", SYMBOL, ignore.case = TRUE)) |> dim() # 2866 of them couldn't find the HGNC

# Adding the HGNC symbol
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# res_Shrink_inter_df_cleanVer <- res_Shrink_inter_df |> dplyr::mutate(gene_id = gsub("\\..*$", "", gene_id))
# annots <- AnnotationDbi::select(org.Hs.eg.db, res_Shrink_inter_df_cleanVer$gene_id,
#                                 columns=c("SYMBOL"), keytype="ENSEMBL")
# annots <- dplyr::rename(annots, gene_id = ENSEMBL)
# annots_unique <- annots %>% distinct(gene_id, .keep_all = TRUE)
# res_Shrink_inter_df_hgnc <- left_join(res_Shrink_inter_df_cleanVer, annots_unique, by = "gene_id") # 2816 of them couldn't find the HGNC

save(res_Shrink_inter_df_hgnc, file= file.path(diff_data_dir,"shrink_interaction_hgnc_tximeta.Rdata"))
 write.table(res_Shrink_inter_df_hgnc, file= file.path(diff_data_dir,"deseq_interaction_results_tximeta.tsv"),
            sep="\t",quote=F,col.names=TRUE, row.names = FALSE)

#  Save to run GO and pathway analysis
# up genes
res_Shrink_inter_df_hgnc |> filter(class == "gained") |> dplyr::mutate(gene_id = gsub("\\..*$", "", gene_id)) |>
    dplyr::select(gene_id) |> 
    distinct() |> 
    write.table(file.path(diff_data_dir,"diff_interaction_upregulated_padj05log2fc1_genes.txt"),
    sep="\t",quote=F,col.names=FALSE, row.names = FALSE)

# down genes
res_Shrink_inter_df_hgnc |> filter(class == "lost") |> dplyr::mutate(gene_id = gsub("\\..*$", "", gene_id)) |>
    dplyr::select(gene_id) |> 
    distinct() |> 
    write.table(file.path(diff_data_dir,"diff_interaction_downregulated_padj05log2fc1_genes.txt"),
    sep="\t",quote=F,col.names=FALSE, row.names = FALSE)

#background genes
res_Shrink_inter_df_hgnc  |> dplyr::mutate(gene_id = gsub("\\..*$", "", gene_id)) |>
    dplyr::select(gene_id) |> 
    distinct() |> 
    write.table(file.path(diff_data_dir,"diff_interaction_background_padj05log2fc1_genes.txt"),
    sep="\t",quote=F,col.names=FALSE, row.names = FALSE)

homer_rna <- file.path(base_dir,"scripts/rna_workflow/01_DEG/run_homer.sh")

up_inter_cmd <- paste(
  homer_rna,
  file.path(diff_data_dir, "diff_interaction_upregulated_padj05log2fc1_genes.txt"),
  file.path(diff_data_dir, "diff_interaction_background_padj05log2fc1_genes.txt"),
  file.path(diff_data_dir, "deseq_homer", "interaction_up_LFC1_padj05")
)

down_inter_cmd <- paste(
  homer_rna,
  file.path(diff_data_dir, "diff_interaction_downregulated_padj05log2fc1_genes.txt"),
  file.path(diff_data_dir, "diff_interaction_background_padj05log2fc1_genes.txt"),
  file.path(diff_data_dir, "deseq_homer", "interaction_down_LFC1_padj05")
)

 system(up_inter_cmd)
 system(down_inter_cmd)


#----------------------------------------------------------------------------------------------------------
# GO and Pathway output
# normoxia ------------------------------------------------------
upsig_GO_data_norm <- read_delim(file.path(diff_data_dir, "deseq_homer", "normoxia_up_LFC1_padj05","biological_process.txt"), delim = "\t") |>
    mutate(pval = exp(1)^logP) |> dplyr::filter(pval < 0.01)

upsig_norm_go <- reduceGO(upsig_GO_data_norm,category = "Upregulated")


## Format and write to table
upgo_norm_table <- write_GO_table(upsig_norm_go)
write_csv(upgo_norm_table, file = file.path(diff_data_dir,"deseq_homer","tables","upsigRNA_GO_normoxia.csv"))


downsig_GO_data_norm <- read_delim(file.path(diff_data_dir, "deseq_homer", "normoxia_down_LFC1_padj05","biological_process.txt"), delim = "\t") |>
    mutate(pval = exp(1)^logP) |> dplyr::filter(pval < 0.01)
downsig_norm_go <- reduceGO(downsig_GO_data_norm, category = "Downregulated")

## Format and write to table
downgo_norm_table <- write_GO_table(downsig_norm_go)
write_csv(downgo_norm_table, file = file.path(diff_data_dir,"deseq_homer","tables","downsigRNA_GO_normoxia.csv"))

# Select 5 each for plotting
upsig_go_norm_plotting <- upsig_norm_go |>
  filter(Term == parentTerm) |>
  filter(parentTerm %in% c("multicellular organismal process",
                           "cell communication",
                           "response to stimulus",
                           "system development",
                           "defense response")) |>
  arrange(`-log10pval`)

downsig_go_norm_plotting <- downsig_norm_go |>
  filter(Term == parentTerm) |>
  filter(parentTerm %in% c("nuclear chromosome segregation",
                           "mitotic nuclear division",
                           "signaling",
                           "regulation of cell population proliferation",
                           "extracellular structure organization")) |>
  arrange(`-log10pval`) 

# Combine into one
go_normoxia_plotting <- bind_rows(upsig_go_norm_plotting, downsig_go_norm_plotting)
go_normoxia_plotting$parentTerm <- factor(go_normoxia_plotting$parentTerm, levels = go_normoxia_plotting$parentTerm)
go_normoxia_plotting$category <- factor(go_normoxia_plotting$category, levels = c("Upregulated", "Downregulated"))

get_fill_scale <- function(condition = c("normoxia", "hypoxia", "interaction")) {
  condition <- match.arg(condition)
  color_map <- switch(condition,
    normoxia = c(
      Downregulated = "#bebfe3",
      Upregulated = "#f2bcd5"
    ),
    hypoxia = c(
      Downregulated = "#8AA3C5",
      Upregulated = "#febfa6"
    ),
    interaction = c(
      Downregulated = "#b5d0be",
      Upregulated = "#f5d7a3"
    ))
}

# Plot all in barplot
normoxia_GO_barplots <- gg_GO_barplot(go_normoxia_plotting, get_fill_scale("normoxia"), title = "GO Terms")

ggsave(filename = file.path(diff_plot_dir,"normoxia_GO_barplots.pdf"),
       plot = normoxia_GO_barplots, width = 5, height = 8, units = "in")
save(normoxia_GO_barplots, file = file.path(diff_plot_dir,"normoxia_GO_barplots.rda"))

#kEGG pathway-------------------------------------
# Read in from Homer
norm_upsig_kegg_data <- read_delim(file.path(diff_data_dir,"deseq_homer","normoxia_up_LFC1_padj05","kegg.txt")) |>
  mutate(pval = exp(1)^logP) |>
  filter(pval < 0.01) |>
  distinct(Term, .keep_all = TRUE) |>
  mutate(`-log10pval` = -log10(pval)) |>
  mutate(category = "Upregulated")

## Format and write to table
norm_upsig_kegg_table <- norm_upsig_kegg_data |>
  dplyr::select(-logP, -pval, -`Entrez Gene IDs`, -category) |>
  relocate(`-log10pval`, .after = Enrichment) |>
  arrange(desc(`-log10pval`))

write_csv(norm_upsig_kegg_table, file = file.path(diff_data_dir,"deseq_homer","tables","KEGG_upsig_normoxia.csv"))


norm_downsig_kegg_data <- read_delim(file.path(diff_data_dir,"deseq_homer","normoxia_down_LFC1_padj05","kegg.txt")) |>
  mutate(pval = exp(1)^logP) |>
  filter(pval < 0.01) |>
  distinct(Term, .keep_all = TRUE) |>
  mutate(`-log10pval` = -log10(pval)) |>
  mutate(category = "Downregulated")

## Format and write to table
norm_downsig_kegg_table <- norm_downsig_kegg_data |>
  dplyr::select(-logP, -pval, -`Entrez Gene IDs`, -category) |>
  relocate(`-log10pval`, .after = Enrichment) |>
  arrange(desc(`-log10pval`))

write_csv(norm_downsig_kegg_table, file = file.path(diff_data_dir,"deseq_homer","tables","KEGG_downsig_normoxia.csv"))


# Plot top 5 significant for each category
norm_upsig_kegg_plotting <- norm_upsig_kegg_data |>
  filter(Term %in% c("Systemic lupus erythematosus",
                     "Cytokine-cytokine receptor interaction",
                     "Neuroactive ligand-receptor interaction",
                     "ErbB signaling pathway",
                     "Wnt signaling pathway")) |>
  mutate(category = "Upregulated") |> arrange(`-log10pval`)

norm_downsig_kegg_plotting <- norm_downsig_kegg_data |>
  filter(Term %in% c("Protein digestion and absorption",
  "PI3K-Akt signaling pathway", 
  "Staphylococcus aureus infection",
  "Cell cycle", 
  "ECM-receptor interaction")) |>
  mutate(category = "Downregulated") |> arrange(`-log10pval`)  
  


# Combine into one
kegg_norm_plotting <- bind_rows(norm_upsig_kegg_plotting, norm_downsig_kegg_plotting)
kegg_norm_plotting$Term <- factor(kegg_norm_plotting$Term, levels = kegg_norm_plotting$Term)
kegg_norm_plotting$category <- factor(kegg_norm_plotting$category, levels = c("Upregulated", "Downregulated"))

# Plot all in barplot
norm_kegg_barplots <- kegg_barplot(kegg_norm_plotting, get_fill_scale("normoxia"), title = "KEGG Pathways")

ggsave(filename = file.path(diff_plot_dir,"normoxia_KEGG_barplots.pdf"),
       plot = norm_kegg_barplots, width = 5, height = 8, units = "in")
save(norm_kegg_barplots, file = file.path(diff_plot_dir,"normoxia_KEGG_barplots.rda"))

# Hypoxia ------------------------------------------------------

upsig_GO_data_hyp <- read_delim(file.path(diff_data_dir, "deseq_homer", "hypoxia_up_LFC1_padj05","biological_process.txt"), delim = "\t") |>
    mutate(pval = exp(1)^logP) |> dplyr::filter(pval < 0.01)

upsig_hyp_go <- reduceGO(upsig_GO_data_hyp,category = "Upregulated")

## Format and write to table
upgo_hyp_table <- write_GO_table(upsig_hyp_go)
write_csv(upgo_hyp_table, file = file.path(diff_data_dir,"deseq_homer","tables","upsigRNA_GO_hypoxia.csv"))


downsig_GO_data_hyp <- read_delim(file.path(diff_data_dir, "deseq_homer", "hypoxia_down_LFC1_padj05","biological_process.txt"), delim = "\t") |>
    mutate(pval = exp(1)^logP) |> dplyr::filter(pval < 0.01)
downsig_hyp_go <- reduceGO(downsig_GO_data_hyp, category = "Downregulated")

## Format and write to table
downgo_hyp_table <- write_GO_table(downsig_hyp_go)
write_csv(downgo_hyp_table, file = file.path(diff_data_dir,"deseq_homer","tables","downsigRNA_GO_hypoxia.csv"))

# Select 5 each for plotting
upsig_go_hyp_plotting <- upsig_hyp_go |>
  filter(Term == parentTerm) |>
  filter(parentTerm %in% c("multicellular organismal process",
                           "response to stimulus",
                           "positive regulation of developmental process",
                           "anatomical structure development",
                           "cell communication")) |>
  arrange(`-log10pval`)

downsig_go_hyp_plotting <- downsig_hyp_go |>
  filter(Term == parentTerm) |> 
  filter(parentTerm %in% c("mitotic cell cycle process",
                           "cell division",
                           "organelle fission",
                           "microtubule cytoskeleton organization involved in mitosis",
                           "metaphase chromosome alignment")) |>
  arrange(`-log10pval`)

# Combine into one
go_hypoxia_plotting <- bind_rows(upsig_go_hyp_plotting, downsig_go_hyp_plotting)
go_hypoxia_plotting$parentTerm <- factor(go_hypoxia_plotting$parentTerm, levels = go_hypoxia_plotting$parentTerm)
go_hypoxia_plotting$category <- factor(go_hypoxia_plotting$category, levels = c("Upregulated", "Downregulated"))

# Plot all in barplot
hypoxia_GO_barplots <- gg_GO_barplot(go_hypoxia_plotting, get_fill_scale("hypoxia"), title = "GO Terms")

ggsave(filename = file.path(diff_plot_dir,"hypoxia_GO_barplots.pdf"),
       plot = hypoxia_GO_barplots, width = 5, height = 8, units = "in")
save(hypoxia_GO_barplots, file = file.path(diff_plot_dir,"hypoxia_GO_barplots.rda"))

#kEGG pathway-------------------------------------
# Read in from Homer
hyp_upsig_kegg_data <- read_delim(file.path(diff_data_dir,"deseq_homer","hypoxia_up_LFC1_padj05","kegg.txt")) |>
  mutate(pval = exp(1)^logP) |>
  filter(pval < 0.01) |>
  distinct(Term, .keep_all = TRUE) |>
  mutate(`-log10pval` = -log10(pval)) |>
  mutate(category = "Upregulated")

## Format and write to table
hyp_upsig_kegg_table <- hyp_upsig_kegg_data |>
  dplyr::select(-logP, -pval, -`Entrez Gene IDs`, -category) |>
  relocate(`-log10pval`, .after = Enrichment) |>
  arrange(desc(`-log10pval`))

write_csv(hyp_upsig_kegg_table, file = file.path(diff_data_dir,"deseq_homer","tables","KEGG_upsig_hypoxia.csv"))


hyp_downsig_kegg_data <- read_delim(file.path(diff_data_dir,"deseq_homer","hypoxia_down_LFC1_padj05","kegg.txt")) |>
  mutate(pval = exp(1)^logP) |>
  filter(pval < 0.01) |>
  distinct(Term, .keep_all = TRUE) |>
  mutate(`-log10pval` = -log10(pval)) |>
  mutate(category = "Downregulated")

## Format and write to table
hyp_downsig_kegg_table <- hyp_downsig_kegg_data |>
  dplyr::select(-logP, -pval, -`Entrez Gene IDs`, -category) |>
  relocate(`-log10pval`, .after = Enrichment) |>
  arrange(desc(`-log10pval`))

write_csv(hyp_downsig_kegg_table, file = file.path(diff_data_dir,"deseq_homer","tables","KEGG_downsig_hypoxia.csv"))


# Plot top 5 significant for each category
hyp_upsig_kegg_plotting <- hyp_upsig_kegg_data |>
  filter(Term %in% c("Systemic lupus erythematosus", "Alcoholism", 
  "Viral carcinogenesis", "p53 signaling pathway", "ABC transporters")) |>
  mutate(category = "Upregulated") |> arrange(`-log10pval`)

hyp_downsig_kegg_plotting <- hyp_downsig_kegg_data |>
  filter(Term %in% c("Cell cycle",
  "Oocyte meiosis","Cell cycle - G2/M transition","Progesterone-mediated oocyte maturation",
  "p53 signaling pathway")) |>
  mutate(category = "Downregulated") |> arrange(`-log10pval`)  
  
# Combine into one
kegg_hypoxia_plotting <- bind_rows(hyp_upsig_kegg_plotting, hyp_downsig_kegg_plotting)
kegg_hypoxia_plotting <- kegg_hypoxia_plotting |> mutate(Term_unique = paste0(Term, " (", category, ")"))

# Order by category first, then by -log10pval within each category
kegg_hypoxia_plotting <- kegg_hypoxia_plotting |>
  arrange(category, `-log10pval`) |>
  mutate(Term_unique = factor(Term_unique, levels = Term_unique))

kegg_hypoxia_plotting$category <- factor(kegg_hypoxia_plotting$category, 
                                          levels = c("Upregulated", "Downregulated"))


# Plot all in barplot
hypoxia_kegg_barplots <- kegg_barplot(kegg_hypoxia_plotting, get_fill_scale("hypoxia"), title = "KEGG Pathways")
hypoxia_kegg_barplots

ggsave(filename = file.path(diff_plot_dir,"hypoxia_KEGG_barplots.pdf"),
       plot = hypoxia_kegg_barplots, width = 5, height = 8, units = "in")
save(hypoxia_kegg_barplots, file = file.path(diff_plot_dir,"hypoxia_KEGG_barplots.rda"))

#-----------------------------------------------------------------------------------------
# Interaction------------------------------------------------------------------------------
upsig_GO_data_inter <- read_delim(file.path(diff_data_dir, "deseq_homer", "interaction_up_LFC1_padj05","biological_process.txt"), delim = "\t") |>
    mutate(pval = exp(1)^logP) |> dplyr::filter(pval < 0.01)

upsig_inter_go <- reduceGO(upsig_GO_data_inter,category = "Upregulated")

## Format and write to table
upgo_inter_table <- write_GO_table(upsig_inter_go)
write_csv(upgo_inter_table, file = file.path(diff_data_dir,"deseq_homer","tables","upsigRNA_GO_interaction.csv"))

downsig_GO_data_inter <- read_delim(file.path(diff_data_dir, "deseq_homer", "interaction_down_LFC1_padj05","biological_process.txt"), delim = "\t") |>
    mutate(pval = exp(1)^logP) |> dplyr::filter(pval < 0.01)
downsig_inter_go <- reduceGO(downsig_GO_data_inter, category = "Downregulated")

## Format and write to table
downgo_inter_table <- write_GO_table(downsig_inter_go)
write_csv(downgo_inter_table, file = file.path(diff_data_dir,"deseq_homer","tables","downsigRNA_GO_interaction.csv"))

# Select 5 each for plotting
upsig_go_inter_plotting <- upsig_inter_go |>
  filter(Term == parentTerm) |>
  filter(parentTerm %in% c("multicellular organismal process",
  "ossification","animal organ development","response to estradiol",
  "regulation of cell population proliferation"  )) |>
  arrange(`-log10pval`)

downsig_go_inter_plotting <- downsig_inter_go |>
  filter(Term == parentTerm) |> 
  filter(parentTerm %in% c("multicellular organismal process",
  "regulation of BMP signaling pathway","signaling","negative regulation of viral process","system development")) |>
  arrange(`-log10pval`)

# Combine into one
go_interaction_plotting <- bind_rows(upsig_go_inter_plotting, downsig_go_inter_plotting)
go_interaction_plotting <- go_interaction_plotting |> mutate(Term_unique = paste0(Term, " (", category, ")"))

# Order by category first, then by -log10pval within each category
go_interaction_plotting <- go_interaction_plotting |>
  arrange(category, `-log10pval`) |>
  mutate(Term_unique = factor(Term_unique, levels = Term_unique))

go_interaction_plotting$category <- factor(go_interaction_plotting$category, levels = c("Upregulated", "Downregulated"))

# Plot all in barplot
interaction_GO_barplots <- gg_GO_barplot(go_interaction_plotting, get_fill_scale("interaction"), title = "GO Terms")
 
ggsave(filename = file.path(diff_plot_dir,"interaction_GO_barplots.pdf"),
       plot = interaction_GO_barplots, width = 5, height = 8, units = "in")
save(interaction_GO_barplots, file = file.path(diff_plot_dir,"interaction_GO_barplots.rda"))

#kEGG pathway-------------------------------------
# Read in from Homer
inter_upsig_kegg_data <- read_delim(file.path(diff_data_dir,"deseq_homer","interaction_up_LFC1_padj05","kegg.txt")) |>
  mutate(pval = exp(1)^logP) |>
  filter(pval < 0.01) |>
  distinct(Term, .keep_all = TRUE) |>
  mutate(`-log10pval` = -log10(pval)) |>
  mutate(category = "Upregulated")

## Format and write to table
inter_upsig_kegg_table <- inter_upsig_kegg_data |>
  dplyr::select(-logP, -pval, -`Entrez Gene IDs`, -category) |>
  relocate(`-log10pval`, .after = Enrichment) |>
  arrange(desc(`-log10pval`))

write_csv(inter_upsig_kegg_table, file = file.path(diff_data_dir,"deseq_homer","tables","KEGG_upsig_interaction.csv"))


inter_downsig_kegg_data <- read_delim(file.path(diff_data_dir,"deseq_homer","interaction_down_LFC1_padj05","kegg.txt")) |>
  mutate(pval = exp(1)^logP) |>
  filter(pval < 0.01) |>
  distinct(Term, .keep_all = TRUE) |>
  mutate(`-log10pval` = -log10(pval)) |>
  mutate(category = "Downregulated")

## Format and write to table
inter_downsig_kegg_table <- inter_downsig_kegg_data |>
  dplyr::select(-logP, -pval, -`Entrez Gene IDs`, -category) |>
  relocate(`-log10pval`, .after = Enrichment) |>
  arrange(desc(`-log10pval`))

write_csv(inter_downsig_kegg_table, file = file.path(diff_data_dir,"deseq_homer","tables","KEGG_downsig_interaction.csv"))


# Plot top 5 significant for each category
inter_upsig_kegg_plotting <- inter_upsig_kegg_data |>
  filter(Term %in% c("Cytokine-cytokine receptor interaction","PI3K-Akt signaling pathway","Melanoma",
  "Staphylococcus aureus infection","Neuroactive ligand-receptor interaction" )) |>
  mutate(category = "Upregulated") |> arrange(`-log10pval`)

inter_downsig_kegg_plotting <- inter_downsig_kegg_data |>
  filter(Term %in% c("ECM-receptor interaction","Axon guidance","Ras signaling pathway",
  "Cell adhesion molecules (CAMs)","Synthesis and degradation of ketone bodies")) |>
  mutate(category = "Downregulated") |> arrange(`-log10pval`)  
  
# Combine into one
kegg_interaction_plotting <- bind_rows(inter_upsig_kegg_plotting, inter_downsig_kegg_plotting)
kegg_interaction_plotting <- kegg_interaction_plotting |> mutate(Term_unique = paste0(Term, " (", category, ")"))

# Order by category first, then by -log10pval within each category
kegg_interaction_plotting <- kegg_interaction_plotting |>
  arrange(category, `-log10pval`) |>
  mutate(Term_unique = factor(Term_unique, levels = Term_unique))

kegg_interaction_plotting$category <- factor(kegg_interaction_plotting$category, 
                                          levels = c("Upregulated", "Downregulated"))


# Plot all in barplot
interaction_kegg_barplots <- kegg_barplot(kegg_interaction_plotting, get_fill_scale("interaction"), title = "KEGG Pathways")
interaction_kegg_barplots

ggsave(filename = file.path(diff_plot_dir,"interaction_KEGG_barplots.pdf"),
       plot = interaction_kegg_barplots, width = 5, height = 8, units = "in")
save(interaction_kegg_barplots, file = file.path(diff_plot_dir,"interaction_KEGG_barplots.rda"))


#=========================================================================================
# save GO and pathwat as EXCEL file
table_dir <- "/work/users/s/e/seyoun/colabs/LCSP/seq_pipelines/rna_output/00_differential_expression/data/deseq_homer/tables"

# GO Terms
go_files <- list(
  "GO_Norm_Up" = "upsigRNA_GO_normoxia.csv",
  "GO_Norm_Down" = "downsigRNA_GO_normoxia.csv",
  "GO_Hypo_Up" = "upsigRNA_GO_hypoxia.csv",
  "GO_Hypo_Down" = "downsigRNA_GO_hypoxia.csv",
  "GO_Inter_Up" = "upsigRNA_GO_interaction.csv",
  "GO_Inter_Down" = "downsigRNA_GO_interaction.csv"
)

# KEGG Pathways
kegg_files <- list(
  "KEGG_Norm_Up" = "KEGG_upsig_normoxia.csv",
  "KEGG_Norm_Down" = "KEGG_downsig_normoxia.csv",
  "KEGG_Hypo_Up" = "KEGG_upsig_hypoxia.csv",
  "KEGG_Hypo_Down" = "KEGG_downsig_hypoxia.csv",
  "KEGG_Inter_Up" = "KEGG_upsig_interaction.csv",
  "KEGG_Inter_Down" = "KEGG_downsig_interaction.csv"
)
# combine all the files
all_files <- c(go_files, kegg_files)

# list and read
excel_data_list <- lapply(all_files, function(f) {
  full_path <- file.path(table_dir, f)
  if (file.exists(full_path)) {
    return(read_csv(full_path))
  } else {
    return(data.frame(Message = "File not found"))
  }
})

# save as excel
output_path <- file.path(table_dir, "GO_Pathway_Analysis_Results.xlsx")
write.xlsx(excel_data_list, file = output_path)


# DESEQ Reqsults table
table_dir <- "/work/users/s/e/seyoun/colabs/LCSP/seq_pipelines/rna_output/00_differential_expression/data/"
deseq_files <- list(
  "Normoxia_treatment" = "deseq_normoxia_results_tximeta.tsv",
  "Hypoxia_treatment" = "deseq_hypoxia_results_tximeta.tsv",
  "Interaction_condition_treatment" = "deseq_interaction_results_tximeta.tsv"
)

# list and read
lapply(names(deseq_files), function(name) {
  # read files
  input_path <- file.path(table_dir, deseq_files[[name]])
  if (file.exists(input_path)) {
    df <- readr::read_tsv(input_path)
    
    # save as different name
    output_filename <- paste0("DESeq2_Result_", name, ".xlsx")
    write.xlsx(df, file = file.path(table_dir, output_filename))
  }
})

# venndiagram

library(eulerr)
library(ggVennDiagram)

load(file= file.path(diff_data_dir,"shrink_normoxia_hgnc_tximeta.Rdata")) # res_Shrink_norm_df_hgnc
load(file= file.path(diff_data_dir,"shrink_hypoxia_hgnc_tximeta.Rdata")) # res_Shrink_hypo_df_hgnc
load(file= file.path(diff_data_dir,"shrink_interaction_hgnc_tximeta.Rdata")) # res_Shrink_inter_df_hgnc

diff_norm <- res_Shrink_norm_df_hgnc |> filter(class != "static")  # 1281
diff_hyp <- res_Shrink_hypo_df_hgnc |> filter(class != "static") # 593
diff_inter <- res_Shrink_inter_df_hgnc |> filter(class != "static") # 294

diff_norm_ensg <- diff_norm$gene_id |> unique()
diff_hyp_ensg <- diff_hyp$gene_id |> unique()
diff_inter_ensg <- diff_inter$gene_id |> unique()
onlyInteraction  <- setdiff(diff_inter_ensg, c(diff_norm_ensg, diff_hyp_ensg)) |> as.matrix()
norm_sub <- res_Shrink_norm_df_hgnc |> filter(gene_id  %in% onlyInteraction) 


venn_list <- list(
  "Normoxia" = diff_norm_ensg,
  "Hypoxia" = diff_hyp_ensg,
  "Interaction" = diff_inter_ensg
)

# Fit euler diagram (area-proportional)
venn_fit <- euler(venn_list)

# Plot with custom styling
deg_venn_grob <- plot(
  venn_fit,
  fills = list(
    fill = c("#A7D3D4", "#df91a3", "#f5d7a3"),
    alpha = 0.4
  ),
  edges = list(col = NA),
  quantities = list(cex = 0.8, col = "black", font = 2),
  labels = list(
    col = "black",
    font = 2,
    cex = 0.9
  )
)
deg_venn_grob

# Convert to ggplot for saving/combining
deg_venn_plot <- ggplot() +
  annotation_custom(
    grob = grid::grobTree(deg_venn_grob)
  ) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(0, 0, 0, 0)
  )

save(deg_venn_plot, file = file.path(diff_plot_dir, "deg_venn_plot.rda"))



#--------------Finding genes 

norm_hyp <- intersect(diff_norm_ensg, diff_hyp_ensg)
norm_inter <- intersect(diff_norm_ensg, diff_inter_ensg)
hyp_inter <- intersect(diff_hyp_ensg, diff_inter_ensg)
all_three <- Reduce(intersect, list(diff_norm_ensg, diff_hyp_ensg, diff_inter_ensg))

# Get only interaction genes
onlyInteraction <- setdiff(diff_inter_ensg, union(diff_norm_ensg, diff_hyp_ensg))

# Function to convert ENSG to SYMBOL (or keep ENSG if no SYMBOL)
ensg_to_symbol <- function(ensg_ids, df_with_symbols) {
  symbols <- sapply(ensg_ids, function(id) {
    symbol <- df_with_symbols$SYMBOL[df_with_symbols$gene_id == id][1]
    if (is.na(symbol) || symbol == "") {
      return(id)  # Return ENSG if no SYMBOL
    } else {
      return(symbol)
    }
  })
  return(symbols)
}

all_three_symbols <- ensg_to_symbol(all_three, res_Shrink_inter_df_hgnc)
onlyInteraction_symbols <- ensg_to_symbol(onlyInteraction, res_Shrink_inter_df_hgnc)

# Separate ENSG and non-ENSG for all_three
ensg_all_three <- all_three_symbols[grep("^ENSG", all_three_symbols)]
non_ensg_all_three <- all_three_symbols[!grepl("^ENSG", all_three_symbols)]

# Sort and combine
sorted_all_three <- c(sort(non_ensg_all_three), ensg_all_three)

# Separate ENSG and non-ENSG for onlyInteraction
ensg_only_inter <- onlyInteraction_symbols[grep("^ENSG", onlyInteraction_symbols)]
non_ensg_only_inter <- onlyInteraction_symbols[!grepl("^ENSG", onlyInteraction_symbols)]

# Sort and combine
sorted_only_inter <- c(sort(non_ensg_only_inter), ensg_only_inter)
#---------------------------------------------------------------------------------------------------------------------

# Combine all the graphs
pdf(file=file.path(diff_plot_dir, "deseq_result_figure.pdf"),width=10.5,height=10.2,bg="transparent")

pageCreate(width = 10.5, height = 10.2, showGuides = FALSE)

plotText("a", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

load(file.path(diff_plot_dir,"PCA_interaction_deseq_tximeta.rda"))
plotGG(plot=interactionPCAplot, x=0.35, y=0.5, height=3, width=4.5)


# Venndiagram
plotText("b", x = 0.1, y = 3.75, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
load(file.path(diff_plot_dir,"deg_venn_plot.rda"))
plotGG(plot=deg_venn_plot, x=0.4, y=4, height=2.5, width=3)

x_pos <- 1
y_start <- 7.5
genes_per_column <- 10
y_spacing <- 0.15
x_spacing <- 0.65

for (i in 1:length(sorted_all_three)) {
  col_num <- ceiling(i / genes_per_column)
  row_num <- (i - 1) %% genes_per_column + 1
  x_current <- x_pos + (col_num - 1) * x_spacing
  y_current <- y_start + (row_num - 1) * y_spacing
  
  plotText(
    label = sorted_all_three[i],
    x = x_current,
    y = y_current,
    fontsize = 7,
    fontfamily = "Helvetica",
    just = "left",
    fontface = "italic",
    fontcolor = "#ca562c"
  )
}

x_pos <- 3.75
y_start <- 4
genes_per_column <- 33
y_spacing <- 0.15
x_spacing <- 0.65

for (i in 1:length(sorted_only_inter)) {
  col_num <- ceiling(i / genes_per_column)
  row_num <- (i - 1) %% genes_per_column + 1
  x_current <- x_pos + (col_num - 1) * x_spacing
  y_current <- y_start + (row_num - 1) * y_spacing
  
  plotText(
    label = sorted_only_inter[i],
    x = x_current,
    y = y_current,
    fontsize = 7,
    fontfamily = "Helvetica",
    just = "left",
    fontface = "italic",
    fontcolor = "#f2bc40"
  )
}

# GO and KEGG

plotText("c", x = 5.2, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

load(file.path(diff_plot_dir,"normoxia_GO_barplots.rda"))
load(file.path(diff_plot_dir,"normoxia_KEGG_barplots.rda"))

x= 5.5
y=0.5
w = 2.4
h = 3

plotText("Dox/PBS (Normoxia)", x = (x+x+w+w)/2+0.25, y = 0.3, just = c("center", "top"), fontfamily = "Helvetica",
         fontsize = 10, fontface = "bold")
plotGG(plot=normoxia_GO_barplots, x=x, y=y, height=h, width=w)
plotGG(plot = norm_kegg_barplots, x= x+ w , y=y, height=h, width=w)


load(file.path(diff_plot_dir,"hypoxia_GO_barplots.rda"))
load(file.path(diff_plot_dir,"hypoxia_KEGG_barplots.rda"))

x= 5.5
y= y + h +0.3 
plotText("Dox/PBS (Hypoxia)", x = (x+4*(w))/2+0.25, y = y-0.2, just = c("center", "top"), fontfamily = "Helvetica",
         fontsize = 10, fontface = "bold")
plotGG(plot=hypoxia_GO_barplots, x=x, y=y, height=h, width=w)
plotGG(plot = hypoxia_kegg_barplots, x= x+ w , y=y, height=h, width=w)

load(file.path(diff_plot_dir,"interaction_GO_barplots.rda"))
load(file.path(diff_plot_dir,"interaction_KEGG_barplots.rda"))

x= 5.5
y= y + h +0.3
plotText("Condition:Treatment interaction", x = (x+4*(w))/2+0.25, y = y-0.2, just = c("center", "top"), fontfamily = "Helvetica",
         fontsize = 10, fontface = "bold")
plotGG(plot=interaction_GO_barplots, x=x, y=y, height=h, width=w)
plotGG(plot = interaction_kegg_barplots, x= x+ w , y=y, height=h, width=w)

dev.off()
