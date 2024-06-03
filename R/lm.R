df_transpose <- t(combined_data)
rownames(df_transpose) <- NULL
# dim 5, 17839
colnames(df_transpose) <- df_transpose[1, ]
df_transpose <- df_transpose [-1, ]
df_numerico <- as.data.frame(apply(df_transpose, 2, as.numeric))
# dim 4, 17839

nueva_columna <- c(0.1, 1.0, 10.0, 100.0)
df <- cbind(nueva_columna, df)
nombre_columna <- "PFU"
colnames(df)[1] <- nombre_columna
# dim 4, 17840

modelo <- lm(PFU ~ ., data = df)
# Error: protect(): protection stack overflow

# Con los archivos signif y cambiando los NA por 0...
df_sin_na <- replace(df, is.na(df_numerico), 0)
write.table(df, 
	        file = "data/coturnix/rnaseq-test/both_ends/DESeq2/signif/inf_vs_genes.txt", 
	        sep = "\t", 
            quote = FALSE, 
            row.names=F, col.names=T)
