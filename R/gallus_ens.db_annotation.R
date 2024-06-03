geneSymbols <- mapIds(org.db, keys=cmp.df.a$geneid, column="GENENAME",keytype="GENENAME",multiVals="first")

# get ENSMBLID
res <- mapIds(ens.db, keys=d$geneid, column=c("GENEID"),keytype="GENENAME")
sum(!is.na(res))
[1] 11470

# is in ann?
ser <- cmp.df.a$geneid %in% bio.ann$entrezgene_accession
sum(ser)
[1] 12098

sum(!is.na(res) | ser)
[1] 13517

# get same gene name
res <- mapIds(org.db, keys=cmp.df.a$geneid, column="GENENAME",keytype="GENENAME",multiVals="first")
sum(!is.na(res))
[1] 23936

# all columns
res <- AnnotationDbi::select(ens.db, keys=d$geneid, column=columns(ens.db),keytype="GENENAME")
write.table(res, 'rnaseq-gallus/both_ends/ens.db.ann.tab', col.names=T, sep='\t')
dim(res)
[1] 3825517      37

res.1 <- AnnotationDbi::select(ens.db, keys=d$geneid, column=columns(ens.db),keytype="GENENAME", multiVals='first')
dim(res.1)



ens.ann <- AnnotationDbi::select(ens.db, 
                  column='GENENAME', keytype= 'GENENAME', keys=rownames(countData), 
                  columns= c('SEQNAME', 'SYMBOL', 'DESCRIPTION', # no longer available
                             'GENENAME', 'GENEID', 'ENTREZID', # empty
                             'TXID', 'TXBIOTYPE', # these make the call fail
                             'PROTEINID', 'UNIPROTID' # no longer available
                            ))

ann <- AnnotationDbi::select(ens.db, 
              keytype= 'GENENAME', keys=rownames(countData), 
              columns= columns(ens.db))

ens.ann.1 <- ens.ann[ ! duplicated(ens.ann$GENEID), ]
ann.1 <- ann[ ! duplicated(ann$GENEID), ]
dim(ens.ann)
[1] 21362    10

dim(ann)
[1] 3825517      37

dim(ens.ann.1)
[1] 11485    10

dim(ann.1)
[1] 11485    37


nc <- as.data.frame(normalized_counts)
nc$gene.name <- rownames(nc)
sum(nc$gene.name %in% ann.1$GENENAME)
anc <- nc
anc <- cbind(anc, ann.1[match(anc$gene.name, ann.1$GENENAME), ])
sum(!is.na(anc$SYMBOL))
for (i in colnames(anc)) cat(i, sum(!is.na(anc[ , i])), '\n')
# annotates 11457
write.table(anc, paste(folder, "normalized_counts.ens.annot.tab", sep='/'), sep='\t', col.names=T)

anc <- cbind(anc, bio.ann[match(anc$gene.name, bio.ann$entrezgene_accession), ])
# increases annotation to 12098
write.table(anc, paste(folder, "normalized_counts.ens.biomart.annot.tab", sep='/'), sep='\t', col.names=T)



oann <- oannotationDbi::select(org.db, 
              keytype= 'GENENAME', keys=rownames(countData), 
              columns= columns(org.db))

oann.1 <- oann[ ! duplicated(oann$GENENAME), ]

dim(oann)
[1] 1780109      23

dim(oann.1)
[1] 23936    23

nc <- as.data.frame(normalized_counts)
dim(nc)
dim(nc)
[1] 23936    18

nc$gene.name <- rownames(nc)
sum(nc$gene.name %in% oann.1$GENENAME)

anc <- nc
anc <- cbind(anc, oann.1[match(anc$gene.name, oann.1$GENENAME), ])
sum(!is.na(anc$SYMBOL))
for (i in colnames(anc)) cat(i, sum(!is.na(anc[ , i])), '\n')
# only 33 get annotation
