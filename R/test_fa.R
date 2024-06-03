# NOTE: requires starting R with
#	R --max-ppsize=500000
# to be able to make linear models.
library(caret)
library(randomForest)
library(xgboost)
library(relaimpo)
library(earth)
library(Boruta)
library(rFerns)
library(randomForest)
library(DALEX)
library(vita)
library(corrplot)
library(factoextra)
library(ctv)
library(psych)
library(psychTools)
library(nFactors)
# the next ones are no longer maintained, do not use
# library(devtools)
#    install_github("tomasgreif/woe")
#    install_github("tomasgreif/riv")
# library(woe)
# library(riv)

# whether we want to test the functions in this file
test_imp.R <- TRUE

#
# Load the data
#

# Load log2FC of significantly altered genes:
#df <- read.table("data/coturnix/rnaseq-test/both_ends/DESeq2/signif/all_log2FC_transpose.txt",header=T)
df <- read.table("all_log2FC.txt",header=T)
fd <- read.table("all_log2FC_transpose.txt",header=T)
# set name of first column
names(fd)[1] <- 'wt.vs.PFU'
fd[,1] <- c(0.1, 1.0, 10.0, 100.0)
# define working data sets
dat <- fd[ , colSums(is.na(fd))==0]			# set all NAs to zero
# dat with zero row (wt.vs.wt, inf.PFU=0, log2FC=0)
datz <- as.data.frame(rbind(rep(0,dim(dat)[2]),dat))		


# load normalized counts of all genes
nc <- read.table(file="normalized_counts.tab", header=T)
# prepend a row with infection levels
nc <- rbind(c(rep(0, 3), rep(0.1, 3), rep(1.0,3), rep(10.0,3), rep(100.0,3)), nc)
rownames(nc)[1] <- 'inf.PFU'
# transpose
tnc <- t(nc)
tnc <- as.data.frame(tnc)
# select only counts of sigfigicantly altered genes:
sig.tc <- cbind(inf.PFU=tnc[ ,1], tnc[ , colnames(tnc) %in% colnames(dat)] )
tc <- sig.tc




dat.cor <- cor(dat, use="pairwise.complete.obs")
corrplot(dat.cor, order='hclust', tl.col="black")	
# the correlation plot should show if there is structure

dat.fa.paral <- fa.parallel(x=dat.cor, fm="minres", fa="fa")
print(dat.fa.paral)
plot(dat.fa.paral)
# this should show a recommended number of factors
dat.scree <- nScree(x=dat.cor, model="factors")
plot(dat.scree)

dat.fa.1 <- fa(dat)		# nfactors defaults to 1 (see help for meaning)
dat.fa.6 <- fa(dat, nfactors=6, fm='uls')	# does not converge
dat.fa.6 <- fa(dat, nfactors=6, fm='minres')# does not converge
dat.fa.6 <- fa(dat, nfactors=6, fm='minres',
               rotate='varimax', scores='tenBerge')  # results are suspect
dat.fa.6 <- fa(dat, nfactors=6, fm='fa',
               rotate='varimax', scores='tenBerge')  # results are suspect
dat.fa.6 <- fa(dat, nfactors=6, fm='minres',
               rotate='quartimax', scores='tenBerge')# suspect
dat.fa.6 <- fa(dat, nfactors=6, fm='minres',
               rotate='bentlerT', scores='tenBerge') # singular
dat.fa.6 <- fa(dat, nfactors=6, fm='wls',
               rotate='bentlerT', scores='tenBerge') # singular

dat.fa <- fa(r=dat.cor, nfactors=6)
# this may need several attempts to get it to converge
# (it may fail with lack of convergence after 1000 iterations)
dat.fa <- fa(r=dat.cor, nfactors=6, fm="minres")
# faster, does not get confidence intervals:
dat.fac <- fac(r=dat.cor, nfactors=6, 
            fm=c"minres", rotate='varimax', scores='tenBerge')
fa.plot(dat.fa)
fa.plot(dat.fac)
fa.diagram(dat.fa)
fa.diagram(dat.fac)


tc.cor <- cor(tc, use="pairwise.complete.obs")
corrplot(tc.cor, order='hclust', tl.col="black")	
# the correlation plot should show if there is structure

tc.fa.paral <- fa.parallel(x=tc.cor, fm="minres", fa="fa")
print(tc.fa.paral)
plot(tc.fa.paral)
vss(tc)
# this should show a recommended number of factors
tc.scree <- nScree(x=tc.cor, model="factors")
plot(tc.scree)

tc.fa <- fa(r=tc.cor, nfactors=6)
# this may need several attempts to get it to converge
# (it may fail with lack of convergence after 1000 iterations)
tc.fa <- fa(r=tc.cor, nfactors=6, fm="minres")
tc.fa <- fa(r=tc.cor, nfactors=6, 
            fm="minres", rotate='varimax', scores='tenBerge')
# faster, does not get confidence intervals:
tc.fac <- fac(r=tc.cor, nfactors=6, 
            fm="minres", rotate='varimax', scores='tenBerge')

fa.plot(tc.fa)
fa.plot(tc.fac)

fa.graph(tc.fa)
fa.graph(tc.fac)

fa.diagram(tc.fa)
fa.diagram(tc.fac)

# use a hierarchical clustering solution to find omega
iclust(tc)
omega(tc)
bassAkward(tc)

# find alpha as an estimate of reliability
alpha(tc)


#
# Bayesian EFA
#
library(BayesFM)

befa(datz, Kmax=6, burnin=10, iter=10)
datz.befa <- post.column.switch(datz.befa)
datz.befa <- post.sign.switch(datz.befa)
plot(datz.befa, what='hppm')
summary(datz.befa)

befa( inf.PFU ~ . , tnc, Kmax=6 )
