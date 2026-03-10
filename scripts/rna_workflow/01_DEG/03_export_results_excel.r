# Creating excel file for the 
# 1) counts data normoxia, hypoxia
# 2) normalized expression in vst_dds_norm
# 3) pathway analysis results for the reactome and wikipath
# 4) GO for the cc and MF 

library(DESeq2)
library(openxlsx)
library(dplyr)
library(readr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GenomicFeatures)

library(rtracklayer)
library(dplyr)
library(data.table)

gtf_path <- "/proj/phanstiel_lab/Reference/human/hg38/annotations/gencode.v49.primary_assembly.annotation.gtf"
gtf <- rtracklayer::import(gtf_path)
save(gtf, file= file.path("/users/s/e/seyoun/Ref/hg38/gencode/v49", "gencode_v49_gtf.Rdata"))

genes <- gtf[gtf$type == "gene"]

gene_annot <- as.data.frame(genes) |>
  dplyr::select(gene_id, gene_name) |>
  distinct(gene_id, .keep_all = TRUE) |>
  dplyr::rename(gene_id = gene_id, SYMBOL = gene_name)

dup_counts <- gene_annot |>
  dplyr::count(SYMBOL, sort = TRUE) |>
  dplyr::filter(n > 1) # 483 genes have same gene symbol they are pseudogene etc other reasons

gene_annot_clean <- gene_annot |>
  add_count(SYMBOL, name = "symbol_count") |>
  mutate(SYMBOL_UNIQUE = ifelse(symbol_count > 1, gsub("\\..*$", "", gene_id), SYMBOL)) |>
  select(-c("symbol_count", "SYMBOL")) |> 
  dplyr::rename(SYMBOL = SYMBOL_UNIQUE)

save(gene_annot_clean, file= file.path("/users/s/e/seyoun/Ref/hg38/gencode/v49", "gencode_v49_gene_annot.Rdata"))


# Paths
base_dir <- "/work/users/s/e/seyoun/colabs/LCSP/seq_pipelines"
rna_dir <- file.path(base_dir, "rna_output")
diff_data_dir <- file.path(rna_dir, "00_differential_expression", "data")
homer_dir <- file.path(diff_data_dir, "deseq_homer")

# Load data
load(file.path(diff_data_dir, "deseq_results_norm_tximeta.Rdata"))
load(file.path(diff_data_dir, "deseq_results_hypo_tximeta.Rdata"))
load(file.path(diff_data_dir, "deseq_results_interaction_tximeta.Rdata"))
load(file.path(diff_data_dir, "deseq_transformed_norm_tximeta.Rdata"))
load(file.path(diff_data_dir, "deseq_transformed_hypo_tximeta.Rdata"))
load(file.path(diff_data_dir, "deseq_transformed_interaction_tximeta.Rdata"))

#=============================================================================
# Function to add HGNC symbol
#=============================================================================

add_hgnc <- function(df) {
  df <- df |>  tibble::rownames_to_column("gene_id") 
  df_hgnc <- df |> add_symbol(gene_annot_clean)
  return(df_hgnc)
}

# 1) Raw Counts with HGNC
counts_norm <- as.data.frame(counts(dds_norm_fit)) |>  add_hgnc()
counts_hypo <- as.data.frame(counts(dds_hypo_fit)) |>  add_hgnc()
counts_inter <- as.data.frame(counts(dds_sub_fit)) |>  add_hgnc()

counts_list <- list(
  "Normoxia_counts" = counts_norm,
  "Hypoxia_counts" = counts_hypo,
  "Interaction_samples_counts" = counts_inter
)

#write.xlsx(counts_list, file = file.path(diff_data_dir, "DESeq2_raw_counts.xlsx"))

for (name in names(counts_list)) {
  save_path <- file.path(diff_data_dir, paste0(name, ".xlsx"))
  openxlsx::write.xlsx(counts_list[[name]], file = save_path, rowNames = FALSE)
}

# 2) VST Normalized Expression with HGNC

vst_norm <- as.data.frame(assay(vst_dds_norm)) |> add_hgnc()
vst_hypo <- as.data.frame(assay(vst_dds_hypo)) |> add_hgnc()
vst_inter <- as.data.frame(assay(vst_dds_inter)) |> add_hgnc()

vst_list <- list(
  "Normoxia_VST" = vst_norm,
  "Hypoxia_VST" = vst_hypo,
  "Interaction_samples_VST" = vst_inter
)

# write.xlsx(vst_list, file = file.path(diff_data_dir, "DESeq2_VST_normalized.xlsx"))

for (name in names(vst_list)) {
  save_path <- file.path(diff_data_dir, paste0(name, ".xlsx"))
  openxlsx::write.xlsx(vst_list[[name]], file = save_path, rowNames = FALSE)
}

#----------------------------------------------------------------------------------------
# 3) GO & Pathways
# Function to read and process HOMER output

read_homer_go <- function(file_path, category, ont = "BP") {
  if (!file.exists(file_path)) {return(data.frame(Message = "File not found"))}
  
  df <- read_delim(file_path, delim = "\t", show_col_types = FALSE) |>
    mutate(pval = exp(1)^logP) |>
    dplyr::filter(pval < 0.01)
  
  if (nrow(df) == 0) {return(data.frame(Message = "No significant terms"))}

  reduced <- reduceGO(df, category = category, ont = ont)
  table_out <- write_GO_table(reduced)
  return(table_out)
}


# Pathway
read_homer_pathway <- function(file_path, category) {
  if (!file.exists(file_path)) {return(data.frame(Message = "File not found"))}
  
  df <- read_delim(file_path, delim = "\t", show_col_types = FALSE) |>
    mutate(pval = exp(1)^logP) |>
    filter(pval < 0.01) |>
    distinct(Term, .keep_all = TRUE) |>
    mutate(`-log10pval` = -log10(pval)) |>
    mutate(category = category)
  
  if (nrow(df) == 0) {return(data.frame(Message = "No significant terms"))}
  
  table_out <- df |>
    dplyr::select(-logP, -pval, -`Entrez Gene IDs`, -category) |>
    relocate(`-log10pval`, .after = Enrichment) |>
    arrange(desc(`-log10pval`))
  
  return(table_out)
}

# Conditions and types

conditions <- c("normoxia", "hypoxia", "interaction")
directions <- list(up = "up_LFC1_padj05",down = "down_LFC1_padj05")

go_types <- list(
  BP = list(file = "biological_process.txt", ont = "BP"),
  CC = list(file = "cellular_component.txt", ont = "CC"),
  MF = list(file = "molecular_function.txt", ont = "MF")
)

pathway_types <- list(
  KEGG = "kegg.txt",
  Reactome = "reactome.txt",
  WikiPathways = "wikipathways.txt"
)

# Both GO and pathway creating excels

for (cond in conditions) {
   print(cond)
  go_list <- list()
  
  for (dir_name in names(directions)) {
    print(dir_name)
    dir_suffix <- directions[[dir_name]]
    folder_path <- file.path(homer_dir, paste0(cond, "_", dir_suffix))
    category_label <- ifelse(dir_name == "up", "Upregulated", "Downregulated")
    
    for (go_name in names(go_types)) {
      print(go_name)
      file_path <- file.path(folder_path, go_types[[go_name]]$file)
      ont <- go_types[[go_name]]$ont
      sheet_name <- paste0(dir_name, "_", go_name)
      go_list[[sheet_name]] <- read_homer_go(file_path, category_label, ont)
    }
  }
  
  go_output <- file.path(homer_dir, "tables", paste0("GO_", cond, ".xlsx"))
  write.xlsx(go_list, file = go_output)
  message("Saved: ", go_output)
  
  pathway_list <- list()
  
  for (dir_name in names(directions)) {
    dir_suffix <- directions[[dir_name]]
    folder_path <- file.path(homer_dir, paste0(cond, "_", dir_suffix))
    category_label <- ifelse(dir_name == "up", "Upregulated", "Downregulated")
    
    for (pw_name in names(pathway_types)) {
      file_path <- file.path(folder_path, pathway_types[[pw_name]])
      sheet_name <- paste0(dir_name, "_", pw_name)
      pathway_list[[sheet_name]] <- read_homer_pathway(file_path, category_label)
    }
  }
  
  pathway_output <- file.path(homer_dir, "tables", paste0("Pathway_", cond, ".xlsx"))
  write.xlsx(pathway_list, file = pathway_output)
  message("Saved: ", pathway_output)
}