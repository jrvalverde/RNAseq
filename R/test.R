for i in df.columns
	if 0 not in set(df[i])
	print(df[[i]])




tc.i.rrf[,colSums(tc.i.rrf) != 0]

tc.i.rrf[rowSums(tc.i.rrf), != 0]

for i in 

sum(dataframe$column_name)

vec = is.zero(df[,1]) 

for (i in names(tc.i.rrf))

	tc.i.rrf$i[tc.i.rrf$i != 0, ]
	
	
results_ca <- data.frame(Gen=c(0,0))
for (i in length(lista_df_ordenados)){
	results_ca <- merge(x=results_ca, y=df[[i]], by='Gen', all=T)
}
	
if (exists(quote(tc.i.lasso)) && ! is.null(tc.i.lasso)) {
    	# convert list output to a table that we can annotate
		lasso.table <- data.frame(var=c('(Intercept)'))
		for (i in names(tc.i.lasso)) { 
    		if (length(names(tc.i.lasso[[i]])) == 0) names(tc.i.lasso[[i]]) <- '(Intercept)'

    		df <- data.frame(var=names(tc.i.lasso[[i]]), coef=tc.i.lasso[[i]])
    		lasso.table <- merge(x=lasso.table, y=df, by='var', all=T)
			colnames(lasso.table)[ dim(lasso.table)[2] ]  <- i

		}
		lasso.table[ is.na(lasso.table) ] <- 0
		# annotate
		elasso.table <- enrich_genes(lasso.table, bio.ann,
									 'var', 'entrezgene_accession')
        # save
		write.table(elasso.table, 'impgenes/imp_lasso_tc.txt', sep='\t', row.names=F)
    }


# for dalex
# select non-negligible
# and remove first value (total)
nn <- subset(tc.i.dalex, mean_dropout_loss > median(tc.i.dalex$mean_dropout_loss))[-1,]

# select most significant among the non-negligible
quantile(nn, 0.8)
ms <- subset(tc.i.dalex, mean_dropout_loss > quantile(nn$mean_dropout_loss, 0.8))[-1,]

