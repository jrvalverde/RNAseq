file_names <- c("raw_PFU:_wt_x_pfu000.1_annotated.tab", 
				"raw_PFU:_wt_x_pfu001.0_annotated.tab", 
                "raw_PFU:_wt_x_pfu010.0_annotated.tab", 
                "raw_PFU:_wt_x_pfu100.0_annotated.tab")

file_names <- c("signif_sorted_PFU:_wt_x_pfu000.1.tab", 
				"signif_sorted_PFU:_wt_x_pfu001.0.tab", 
                "signif_sorted_PFU:_wt_x_pfu010.0.tab", 
                "signif_sorted_PFU:_wt_x_pfu100.0.tab")

data_list <- list()

for (file in file_names) {
  file_data <- read.table(file, header=TRUE)
  data <-  data.frame(gene = rownames(file_data), l2fc=file_data$log2FoldChange)

  name <- sub(".*_wt", "wt", file)
  name <- sub(".tab", "_l2fc", name)
  colnames(data) <- c("gene", name)
  data_list[[file]] <- data
}

combined_data <- data_list[[1]] 
for (i in 2:length(data_list)) {
  combined_data <- merge(combined_data, data_list[[i]], by="gene", all = TRUE)
  print(head(combined_data))
}

rownames(combined_data) <- combined_data$rows
head(combined_data)


#all_log2FC <- write.table (combined_data)
#write.table(combined_data, file = "all_log2FC.tab", sep = "\t", quote = FALSE, row.names = TRUE)

write.table(combined_data, file = "all_log2FC.txt", 
		    sep = "\t", 
            quote = FALSE, 
            row.names=F, col.names=T)


####
combined_data <- do.call(cbind, data_list)

combined_data <- merge (data_list, by=row.names, all=TRUE)

colnames(combined_data) <- c("pfu000.1",
                             "pfu001.0",
                             "pfu010.0",
                             "pfu100.0")

colnames(combined_data, col_name)
