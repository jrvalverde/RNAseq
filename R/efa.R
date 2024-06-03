df_transpose <- t(combined_data)

write.table(df_transpose, 
	    file = "data/coturnix/rnaseq-test/both_ends/DESeq2/all_log2FC_transpose.txt", 
	    sep = "\t", 
            quote = FALSE, 
            row.names=T, col.names=F)

# Eliminar los rownames
rownames(df_transpose) <- NULL
df_transpose <- df_transpose [-1, ] #eliminar la fila 1

# Convertir a valores numéricos
df_numerico <- as.data.frame(apply(df_transpose, 2, as.numeric))

## TRATAMIENTO DE VALORES NULOS
# A) Reemplazar valores NA por 0 en todo el dataframe
df_sin_na <- replace(df, is.na(df_numerico), 0)

# B) Eliminar las columnas que contienen valores NA
df_sin_na <- df_transpose[, colSums(is.na(df_transpose)) == 0]

mydata <- df_sin_na

## EVALUAR LA MATRIZ DE CORRELACIÓN
write.csv(cor(mydata)>0.8, file="Suspect_Correlations.csv") # correlación de Pearson
write.csv(cor(mydata), file="Correlation_Values.csv")  # probar con Spearman ?

## DETERMINAR EL NÚMERO DE FACTORES
ev <- eigen(cor(mydata)) # obtener autovalores
ev$values

scree(mydata, pc=FALSE)
fa.parallel(mydata, fa="fa")

## EXTRAER FACTORES
# If you wish to begin with the assumption of correlated (non-independent) factors, 
# you have two oblique rotation option: promax - oblimin
# If you wish to begin with the assumption of uncorrelated (independent) factors, 
# you have a few orthogonal rotation options: varimax - quartimax - equamax

Nfacs <- 4  # cambiar el número de factores
fit <- factanal(mydata, Nfacs, rotation="promax") # seleccionar método rotación
print(fit, digits=2, cutoff=0.3, sort=TRUE)

#plots
load <- fit$loadings[,1:2]
plot(load,type="n") # set up plot
text(load,labels=names(mydata),cex=.7)

library(psych)
loads <- fit$loadings
fa.diagram(loads)






