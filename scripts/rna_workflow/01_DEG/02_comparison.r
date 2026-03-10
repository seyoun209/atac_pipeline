library(eulerr)
library(ggVennDiagram)
library(ggplot2)
library(DESeq2)
library(ggrepel)
library(dplyr)
library(EnhancedVolcano)

load(file= file.path(diff_data_dir,"shrink_normoxia_hgnc_tximeta.Rdata")) # res_Shrink_norm_df_hgnc
load(file= file.path(diff_data_dir,"shrink_hypoxia_hgnc_tximeta.Rdata")) # res_Shrink_hypo_df_hgnc
load(file= file.path(diff_data_dir,"shrink_interaction_hgnc_tximeta.Rdata")) # res_Shrink_inter_df_hgnc

# ENSG00000111215 this ensg should be PRH4 but annotation did wrong as PRH1. 

# all_res <- list(
#   res_Shrink_norm_df_hgnc = res_Shrink_norm_df_hgnc,
#   res_Shrink_hypo_df_hgnc = res_Shrink_hypo_df_hgnc,
#   res_Shrink_inter_df_hgnc = res_Shrink_inter_df_hgnc)

# for (i in names(all_res)) {
#   # Locate the row index
#   target_idx <- which(all_res[[i]]$gene_id == "ENSG00000111215")
  
#   if (length(target_idx) > 0) {
#     all_res[[i]]$SYMBOL[target_idx] <- "PRH4"
#   }
# }

# list2env(all_res, envir = .GlobalEnv)

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



# Volcano plot for the differential normoxia -------------------------------------------
## Normoxia  (Doxorubicin vs pbs)
  # EnhancedVolcano(res_Shrink_norm_df_hgnc,
  #   lab = res_Shrink_norm_df_hgnc$SYMBOL,
  #   x = 'log2FoldChange',
  #   y = 'padj')



make_volcano <- function(df, up_color, down_color, title = NULL) {
  df <- df |> filter(!is.na(padj)) |> mutate(
      label_name = ifelse(is.na(SYMBOL), gene_id, SYMBOL),
      neg_log10_padj = -log10(padj),
      color_group = case_when(
        padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
        padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
        TRUE ~ "NS"
      )
    )
  
  sig_up <- df |> filter(color_group == "Upregulated") |>  arrange(desc(log2FoldChange)) |> head(20)
  sig_down <- df |>  filter(color_group == "Downregulated") |>  arrange(log2FoldChange) |> head(20)
  select_genes <- c(sig_up$label_name, sig_down$label_name)
  
  df <- df |>  mutate(label = ifelse(label_name %in% select_genes, label_name, NA))
  y_axis_p <- df$neg_log10_padj |> max() |> round(digits = -1)

  p <- ggplot(df, aes(x = log2FoldChange, y = neg_log10_padj, color = color_group)) +
    geom_point(size = 1, alpha = 0.7) +
    geom_text_repel(aes(label = label),
      size = 2,fontface = "italic",
      color = "black",max.overlaps = Inf,na.rm = TRUE) +
    scale_color_manual(
      values = c("Upregulated" = up_color, "Downregulated" = down_color, "NS" = "grey70")) +
    # scale_x_continuous(breaks = seq(-10, 10, 2), expand = c(0, 0)) +
    scale_x_continuous(limits = c(-10, 10), breaks = seq(-10, 10, 2)) + 
    scale_y_continuous(limits = c(0, y_axis_p), breaks = seq(0, y_axis_p, length.out = 5)) +
    # coord_cartesian(xlim = c(-10, 10), ylim = c(0, 14), clip = "off") +
    labs(x = expression(Log[2]~Fold~Change), y = expression(-Log[10]~adjusted~p-value)) +
    theme(
      axis.line.y = element_line(linewidth = 0.25),
      axis.line.x = element_line(linewidth = 0.25),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_line(color = "black", linewidth = 0.25),
      axis.ticks.length.y = unit(-0.1, "cm"),
      axis.title.x = element_text(size = 8, family = "Helvetica", margin = margin(t = 5)),
      axis.title.y = element_text(size = 8, family = "Helvetica", margin = margin(r = 5)),
      text = element_text(family = "Helvetica"),
      axis.text = element_text(color = "black", size = 7),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.grid = element_blank(),
      legend.position = c(0.9, 0.89),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "transparent", color = NA),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, "cm")
    )
  
  return(p)
}

p_normoxia <- make_volcano(res_Shrink_norm_df_hgnc, 
                            up_color = "#f2bcd5", 
                            down_color = "#bebfe3")

p_hypoxia <- make_volcano(res_Shrink_hypo_df_hgnc, 
                           up_color = "#febfa6", 
                           down_color = "#8AA3C5")

# For this one change the make_volcano for length.out = 5 instead by =10
p_interaction <- make_volcano(res_Shrink_inter_df_hgnc, 
                               up_color = "#f5d7a3", 
                               down_color = "#b5d0be")


save(p_normoxia, p_hypoxia, p_interaction, file = file.path(diff_plot_dir, "volcanoplots.rda"))

pdf(file=file.path(diff_plot_dir, "deseq_volcano_figure.pdf"),width=7.5,height=13.5,bg="transparent")

pageCreate(width = 7.5, height = 13.5, showGuides = FALSE)

plotText("a", x = 0.1, y = 0.1, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

load(file.path(diff_plot_dir, "volcanoplots.rda"))
plotGG(plot=p_normoxia, x=0.35, y=0.5, height=4, width=7)

plotText("b", x = 0.1, y = 4.5, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")
plotGG(plot=p_hypoxia, x=0.35, y=4.9, height=4, width=7)

plotText("c", x = 0.1, y = 8.9, just = c("left", "top"), fontfamily = "Helvetica",
         fontsize = 14, fontface = "bold")

plotGG(plot=p_interaction, x=0.35, y=9.3, height=4, width=7)

dev.off()