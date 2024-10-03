#########VOLCANO PLOT####################
getwd()
setwd("/home/rohan/Documents/Sakshi/plots")
# Load necessary libraries
library(ggplot2)
library(ggrepel)

# res object has my DESeq2 results
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

#  thresholds for significance
log2fc_threshold <- 1
pvalue_threshold <- 0.05

# Create a new column to classify genes
res_df$Expression <- "Non-significant"
res_df$Expression[res_df$padj < pvalue_threshold & res_df$log2FoldChange > log2fc_threshold] <- "Upregulated"
res_df$Expression[res_df$padj < pvalue_threshold & res_df$log2FoldChange < -log2fc_threshold] <- "Downregulated"

# Selected top 5 upregulated genes
top_upregulated <- res_df[res_df$Expression == "Upregulated", ]
top_upregulated <- top_upregulated[order(top_upregulated$log2FoldChange, decreasing = TRUE), ]
top_upregulated <- head(top_upregulated, 5)

# Selected top 7 downregulated genes
top_downregulated <- res_df[res_df$Expression == "Downregulated", ]
top_downregulated <- top_downregulated[order(top_downregulated$log2FoldChange), ]
top_downregulated <- head(top_downregulated, 5)

# Combine these into a single data frame for labeling
label_data <- rbind(top_upregulated, top_downregulated)


volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Expression)) +
   geom_point(alpha = 0.5, size = 3) +
   scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "#00FFFF", "Non-significant" = "#A9A9A9")) +
   labs(title = "",
        x = "Log2FC",
        y = "-Log10 Padj") +
   theme_minimal() +
   theme(legend.position = "top") +
   geom_label_repel(data = label_data, aes(label = gene), size = 3, color = "black", 
                    box.padding = 0.5, point.padding = 7.5,
                    max.overlaps = Inf) +  
   geom_vline(xintercept = log2fc_threshold, linetype = "dashed", color = "#DCDCDC") +
   geom_vline(xintercept = -log2fc_threshold, linetype = "dashed", color = "#DCDCDC") +
   geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "#DCDCDC")


# Print the plot
print(volcano_plot)

# Save the plot as PDF
ggsave("volcano_plot.pdf", plot = volcano_plot, width = 10, height = 8, dpi = 300, device = cairo_pdf)

# Save the plot as PNG
ggsave("volcano_plot.png", plot = volcano_plot, width = 10, height = 8, dpi = 300)

# Save the plot as JPEG
ggsave("volcano_plot.jpeg", plot = volcano_plot, width = 10, height = 8, dpi = 300, device = "jpeg")
