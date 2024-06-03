

#df <- read.table("data/coturnix/rnaseq-test/both_ends/DESeq2/signif/all_log2FC_transpose.txt",header=T)
df <- read.table("all_log2FC.txt",header=T)
fd <- read.table("all_log2FC_transpose.txt",header=T)

head(names(fd))
head(fd[1])
names(fd)[1] <- 'wt.vs.PFU'
fd[,1] <- c(0.1, 1.0, 10.0, 100.0)

dat <- fd[ , colSums(is.na(fd))==0]
datz <- as.data.frame(rbind(rep(0,1182),dat))		# dat with zero row

# try random forest
# Random forest can be very effective to find a set of predictors that best
# explains the variance in the response variable.
library(caret)
library(randomForest)
library(varImp)
# do RF (too few dependent variable levels)
regressor <- randomForest(wt.vs.PFU ~ . , data=dat, importance=TRUE)  # fit the random forest with default parameter
### alternative: make a GLM
### regressor <- glm(wt.vs.PFU ~ ., data=dat)


caret::varImp(regressor) 		# get variable importance, based on mean decrease in accuracy
caret::varImp(regressor, conditional=TRUE) # conditional=True, adjusts for correlations between predictors
varImp::varimpAUC(regressor)    # more robust towards class imbalance.

# try xgboost
library(caret)
library(xgboost)
regressor=train(wt.vs.PFU ~ ., data=dat, 
                method = "xgbTree", 
                trControl = trainControl("cv", number = 10),
                scale=T)
varImp(regressor)

# try relimp
library(relaimpo)
regressor <- lm(wt.vs.PFU ~ . , data=dat) # fit lm() model
relImportance <- calc.relimp(regressor, type = "lmg", rela = TRUE) # calculate relative importance scaled to 100
sort(relImportance$lmg, decreasing=TRUE) # relative importance

# try MARS (earth package)
# The earth package implements variable importance based on Generalized cross
# validation (GCV), number of subset models the variable occurs (nsubsets) and
# residual sum of squares (RSS).
library(earth)
regressor <- earth(wt.vs.PFU ~ . , data=dat) # build model
ev <- evimp (regressor) # estimate variable importance
plot (ev)



# 5. Step-wise Regression Method
# 
# If you have large number of predictors , split the Data in chunks of 10
# predictors with each chunk holding the responseVar.

base.mod <- lm(wt.vs.PFU ~ 1 , data=dat) # base intercept only model
all.mod <- lm(wt.vs.PFU ~ . , data=dat) # full model with all predictors
stepMod <- step(base.mod, 
		scope = list(lower = base.mod, upper = all.mod), 
                direction = "both", 
                trace = 1, 
                steps = 1000) # perform step-wise algorithmshortlistedVars <- names(unlist(stepMod[[1]])) # get the shortlisted variable.shortlistedVars <- shortlistedVars[!shortlistedVars %in% "(Intercept)"] # remove intercept

# The output might include levels within categorical variables, since
# 'stepwise' is a linear regression based technique.
# 
# If you have a large number of predictor variables, the above code may need to
# be placed in a loop that will run stepwise on sequential chunks of
# predictors. The shortlisted variables can be accumulated for further analysis
# towards the end of each iteration. This can be very effective method, if you
# want to
# 
# * Be highly selective about discarding valuable predictor variables. 
# 
# * Build multiple models on the response variable.


# 6. Boruta Method
# 
# The 'Boruta' method can be used to decide if a variable is
# important or not.


library(Boruta)
# Decide if a variable is important or not using Boruta
boruta_output <- Boruta(wt.vs.PFU ~ . , data=dat, doTrace=2) # perform Boruta search
boruta_signif <- names(boruta_output$finalDecision[boruta_output$finalDecision %in% c("Confirmed", "Tentative")]) 
# collect Confirmed and Tentative variables
# for faster calculation(classification only)

library(rFerns)
boruta.train <- Boruta(factor(wt.vs.PFU)~., data=dat, 
                       doTrace = 2, 
                       getImp=getImpFerns, 
                       holdHistory = F)
print(boruta.train)
 
boruta_signif <- names(boruta.train$finalDecision[boruta.train$finalDecision %in% c("Confirmed", "Tentative")]) 
# collect Confirmed and Tentative variables
print(boruta_signif)

#
#
getSelectedAttributes(boruta_signif, withTentative = F)
boruta.df <- attStats(boruta_signif)
print(boruta.df)


# 7. Information value and Weight of evidence Method
# 		Use for binary outcomes (try e.g. with 0 vs 100)
library(devtools)
   install_github("tomasgreif/woe")
library(woe)
library(riv)

iv_df <- iv.mult(dat, y="wt.vs.PFU", summary=TRUE, verbose=TRUE)
iv <- iv.mult(dat, y="wt.vs.PFU", summary=FALSE, verbose=TRUE)
iv_dfiv.plot.summary(iv_df) # Plot information value summary
# Calculate weight of evidence variables
data_iv <- iv.replace.woe(dat, iv, verbose=TRUE) # add woe variables to original data frame.

# The newly created woe variables can alternatively be in place of the original
# factor variables.

# 
# 8. Learning Vector Quantization (LVQ) Method

library(caret)
control <- trainControl(method="repeatedcv", number=10, repeats=3)
# train the model
regressor<- train(wt.vs.PFU ~ ., data=dat, 
                  method="lvq", 
                  preProcess="scale", 
                  trControl=control)
# estimate variable importance
importance <- varImp(regressor, scale=FALSE)

#
# 9. Recursive Feature Elimination RFE Method

library(caret)
# define the control using a random forest selection function
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# run the RFE algorithm
results <- rfe(dat[,-1], dat[,1], sizes=c(1:8), rfeControl=control)
# summarize the results
# list the chosen features
predictors(results)
# plot the results
plot(results, type=c("g", "o"))

#
#
# 10. DALEX Method

library(randomForest)
library(DALEX)
regressor <- randomForest(wt.vs.PFU ~ . , data=dat, importance=TRUE) 
# fit the random forest with default parameter
### alternate: make a biomial GLM
### regressor <- glm(wt.vs.PFU ~ ., data=dat, family="binomial")
# Variable importance with DALEX
explained_rf <- explain(regressor, data=dat, y=dat$wt.vs.PFU)
# Get the variable importances
varimps <- variable_dropout(explained_rf, type='raw')
print(varimps)
plot(varimps)
feature_importance(explained_rf, type='raw')


#
# 11. VITA

library(vita)
regressor <- randomForest(wt.vs.PFU ~ . , data=dat, importance=TRUE) 
# fit the random forest with default parameter
pimp.varImp.reg <- PIMP(dat, dat$wt.vs.PFU, regressor,S =10, parallel=TRUE)
print(pimp.varImp.reg)
print(pimp.varImp.reg$VarImp)

pimp.varImp.reg$VarImp
sort(pimp.varImp.reg$VarImp,decreasing = T)

#
# 12. Genetic Algorithm

library(caret)
# Define control function
ga_ctrl <- gafsControl(functions = rfGA,  method = "cv", repeats = 3)
# another option is `caretGA`
# Genetic Algorithm feature selection
ga_obj <- gafs(x=dat[, -1],
               y=dat[, 1],
               iters = 3,  # normally much higher (100+)        
               gafsControl = ga_ctrl)

print(ga_obj)
# Optimal variables
print(ga_obj$optVariables)

#
# 13. Simulated Annealing

library(caret)
# Define control function
sa_ctrl <- safsControl(functions = rfSA,
                       method = "repeatedcv",            
                       repeats = 3,            
                       improve = 5)
# improve = n iterations without improvement before a reset
# Simulated Annealing Feature Selection
set.seed(20240419)
sa_obj <- safs(x=dat[, -1],
               y=dat[, 1],
               safsControl = sa_ctrl)
print(sa_obj)
# Optimal variables
print(sa_obj$optVariables)

#
# 14. Correlation Method

library(caret)
# calculate correlation matrix
correlationMatrix <- cor(dat[,-1])
# summarize the correlation matrix
print(correlationMatrix)
# find attributes that are highly correlated (ideally >0.75)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.5)
# print indexes of highly correlated attribute
print(highlyCorrelated)




# EFA

dat=fd
dat.zero <- dat  
dat.zero[is.na(dat.zero)] <- 0
dat.std <- data.frame(scale(dat.zero , center=TRUE, scale=TRUE))
# describe(dat.std)

# correlation plot from 'corrplot'
{
    library('corrplot')
    corrplot(cor(dat, use="complete.obs"), 
             order = "original", 		# also try 'hclust'
             tl.col='black', tl.cex=.75) 
}


{
    # EFA with no rotation 
    # (also try rotation = varimax, promax, oblimin)
    # (also try scores="regression", "Barlett")
    efa.result <- factanal( ~ . , data=dat, factors=10, rotation="none", na.action=na.exclude) #note the formula specification allows NA

    # look at sums of squared loadings (eigenvalues, variance explained by factor)
    # factor 1
    loadings_fac1 = efa.result$loadings[,1] 
    eigenv_fac1 <- sum(loadings_fac1 ^ 2)

    # uniqueness (1 - communality)
    #    where communality = sum(ss of all factor loadings of a variable)
    efa.result$uniquenesses
    # variable 1
    # Calculate communality
    loadings_V1 <- efa.result$loadings[1,]  #loadings for first variable (1st row of Lambda)
    communality_V1 <- sum(loadings_V1 ^ 2)  #SS of factor loadings
    uniqueness_V1 <- 1-communality_V1

    # plot factor loadings for first two factors
    load <- efa.result$loadings[,1:2]
    plot(load, type="n") # set up plot 
    text(load,labels=names(dat),cex=.7) # add variable names

}

{
    library(psych)

    cor.dat <- cor(dat[, -1], use="pairwise.complete.obs")
    cor.dat.zero<- cor(dat[, -1], use="pairwise.complete.obs")
    cor.dat.std <- cor(dat[, -1], use="pairwise.complete.obs")

   # get parallel analysis scree plots
    fa.parallel(x=cor.dat, fm="minres", fa="fa")
    
    library(nfactors)
    nScree(x=cor.dat, model="factors")
    plot(nScree(x=cor.dat,model="factors"))

    
    efa.res <- fa(r=cor.dat, nfactors=10, rotate="oblimin", fm="pa")
efa.res <- fa(r=cor.dat.zero, nfactors=10, rotate="oblimin", fm="pa")
    efa.res <- fa(r=cor.dat.std, nfactors=10, rotate="oblimin", fm="pa")

}

{
    library("devtools")
    devtools::install_github('hfgolino/EGA')

    library(EGAnet)
    ega.res <- EGA(data=dat.std, model='glasso')	# fails (FORTRAN error)
    ega.res <- EGA(data=dat.std)
    
}
