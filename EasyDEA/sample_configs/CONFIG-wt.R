ALIGN 					<- TRUE
PAIRED                  <- FALSE
BOTH 					<- FALSE
USE.ONLINE.ANNOTATION 	<- TRUE
USE.EDGER 				<- FALSE
USE.DESEQ2 				<- TRUE
reference 				<- "./ref/Gallus_gallus_genome.fa"
annotation 				<- "./ref/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.112.gtf"
release 				<- "GRCg7b"
target.organism 		<- 'Gallus gallus'
ens.version 			<- '112'
mart.name 				<- 'ggallus_gene_ensembl'
org.package 			<- "org.Gg.eg.db"
ncbi.taxid              <- '9031'       # Human 9606
kegg.organism			<- "gga"
n.genes 				<- 1000         # 'all' => q<0.05
fastq.dir 				<- 'fastq-qc-named'  # needed if ALIGN=TRUE
alignment.dir 			<- 'aligned-named-cnb'
feature.count.dir		<- "feature_counts"
my.name 				<- 'J. R. Valverde'
my.email 				<- '<jrvalverde@cnb.csic.es>'
my.user 				<- 'sci'
my.password 			<- "password"
metadata 				<- './SampleInfo.tab'
cpm.threshold 			<- 0.5
significance.threshold 	<- 1
#design.column 			<- c("sample", "status", "type", "persistent", "recovered", "wt". "resistant"
design.column 			<- "wt"     # should not contain '_' !!!
rnaseq.out 				<- paste('rnaseq', "by", design.column, 'cnb', sep='.')
#rnaseq.out              <- "rnaseq.out"
config.file 			<- NULL
INTERACTIVE 			<- FALSE
VERBOSE 				<- TRUE
#print(alignment.dir)
#print(rnaseq.out)
