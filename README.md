## Integrating multiple Omics from GBM data

The following scripts illustrate how to implement the models presented in *Bernal Rubio et al., G3, 2018*, the scripts are also provided at: [https://github.com/anainesvs/BERNAL_RUBIO_etal_G3_2018](https://github.com/anainesvs/BERNAL_RUBIO_etal_G3_2018), please refer to that webpage for updates.

**Contact**: avazquez@msu.edu
**Contact**: jlbernalr@gmail.com

#### (1) Installing BGLR

The code below illustrates how to install and load the necessary package from CRAN using `install.packages()`.

```R
   install.packages(pkg='BGLR')    # install BGLR
   library(BGLR); 
```  

#### (2)  Loading data

   * `XF`: a matrix for clinical covariates (age, is.white, ... name them in order, is.white is a 0/1 variable, where 1 indicates whether subject is white caucassian). 
   * `X.ge`: a matrix for gene expression. 
   * `X.meth`: a matrix for methylation values at various sites.
   * `X.cnv`: a matrix of mean CNV intensity per gen.
   * `y`: a matrix with 4 columns. The first column 'time' is the time to last follow up or event, columns 2 and 3 are variables `a` and `b` required by BGLR to analyze censored data, see description on this analysis on BGLR can be found on [https://cran.r-project.org/web/packages/BGLR/vignettes/BGLR-extdoc.pdf](https://cran.r-project.org/web/packages/BGLR/vignettes/BGLR-extdoc.pdf). Finally, the column 4 is an indicator variable (0/1) with 0 indicating whether the patient is alive at the 5th month after diagnosis.
   * `batch`: A vector containing the batch by subject.
   * `folds`: A vector containing fold numbers (1 and 2) for cross validations. The two folds are balanced so they have the same proportions of dead/alive individuals. 
  
After downloading, files can be loaded into R with:

```R 
  load('OMIC_DATA.rda')
  #Checking files present and dimensions
  ls()
  #[1] "X.ge"   "X.meth" "X.cnv"  "XF"     "y"   
  
```

Clinical covariates, as well as omic datasets should be edited by removing outliers (for instance, observations outside of mean +- 3SD) and ordering samples id equally across datasets. Also, incidence matrix should be conformable in dimensions for subsequent analysis.

##### (3) Pre-adjust gene expression matrix by batch

THIS WILL PROBABLY BE DONE BEFORE, THEN RM THE SECTION
ALSO, MENTION BATCH EF CORRECTION IN (2) DATA

For demostration, in this section we show how to pre-correct the incidence matrix corresponding to gene expression by batch effects using mixed models. The input used is the matrix "X.ge" centered, and the batch code for each sample. Using a linear mixed model, the pre-correction regresses each column of the Xge matrix into the corresponding batch effect, assuming them as random effects.

```R
  Xge<- scale(X.ge, scale=TRUE, center=FALSE) 

  #Load library lme4 to fit a linear mixed model
   library(lme4)
   adj_batch<-list()
   for(i in 1:ncol(X.ge)){
      adj_batch[[i]]<-residuals(lmer(as.vector(X.ge[,i])~1| batch_id))
   }

   X.ge_batch<-t(do.call(rbind,adj_batch))
   rownames(X.ge_batch)<-rownames(X.ge)
   colnames(X.ge_batch)<-colnames(X.ge)
```

#### (4) Computing similarity matrices

 Some of the models fitted in the study use similarity matrices of the form G=XX' 
 computed from omics. The following code illustrates how to compute this matrix for 
 SNP genotypes. A similar code could be use to compute a G-matrix for methylation 
 or other omics.
 
```R 
   #Computing a similarity matrix for SNP genotypes
   X.meth_sc<- scale(X.meth, scale=TRUE, center=TRUE) #centering and scaling
   Gmeth<-tcrossprod(X.meth)                          #computing crossproductcts
   Gmeth<-Gmeth/mean(diag(Gmeth))                     #scales to an average diagonal value of 1.
   
   #The same procedure was followed for other omics.
```
**NOTE**: for larger data sets it may be more convinient to use the `getG()` function available in [BGData](https://github.com/quantgen/BGData) R-package. This function allows computing G without loading all the data in RAM and offers methods for multi-core computing. 

#### (4)  Exploring demographic effects in the principal components (e.g., Methylation-derived PC vs age)

The following code shows how to search for effects of demographics such as age on the omic derived PCs, for instance from Methylation.

```R
   #Gmeth has been obtained as described in previous section
   EVD<-eigen(Gmeth,symmetric=TRUE)
   pc_meth<-EVD$vectors
   rownames(pc_meth)<-colnames(pc_meth)<-rownames(Gmeth)

```

Principal component coefficients can be plotted to visualize stratification in the dataset by gender as follows:


```R
   library(ggplot2)
   #Plot of Methylation-derived PCs by age, for instance PC1 vs PC2
   tmp<-as.data.frame(pc_meth)
   
   ***fix here below, eg: color certain older age
   
   tmp$gender<-ifelse(XF$gender==0,"female","male")
   ggplot(temp, aes(temp[,1], temp[,2], col=gender)) + 
   geom_point(aes(shape=gender), size=2) +   # draw points
   labs(title="SNP Clustering", 
       subtitle="PC1 and PC2 by gender")
```

A similar approach can be implemented for additional omics (gene expression and methylation).

#### (5) Variance explained by omics

The following code illustrates how to use BGLR to fit a fixed effects model. The incidence matrix "XF" as well as the object "y" described in section (3) are used. There is no column for intercept in XF because BGLR adds the intercept automatically. Predictors are given to BGLR in the form a two-level list. The argument `save_at` can be used to provide a path and a pre-fix to be added to the files saved by BGLR. For further details see [Perez-Rodriguez and de los Campos, Genetics, 2014](http://www.genetics.org/content/genetics/198/2/483.full.pdf). The code also shows how to retrieve estimates of effects and of success probabilities. In the examples below we fit the model using the default number of iterations (1,500) and burn-in (500). In practice, longer chains are needed and therefore, the user can increase the number of iterations or the burn-in using the arguments `nIter` and `burnIn` of `BGLR`.

```R
   ### Inputs: XF, G.meth and y
   ETA.cov_meth<-list(cov = list(~XF[,2:3?   FIX VECTOR OF COLUMNS WITH TMZ AND AGE  ],model='FIXED'),
                      meth = list(X=PC.meth,model="RKHS"))
   # Fitting the model
   fm.meth<-BGLR(y=y[,4],a=y[,5],b=y[,6],ETA=ETA.cov_meth,saveAt="Example_")
   # Retrieving posterior variances
   var.meth <- fm.meth$ETA$meth$varU
   var.error <- fm.meth$varE
   
```

```R
   ### Inputs: XF, G.meth, G.ge and y
   Aca con 2 omics !!!!!!!!!!!
   ETA.meth<-list(list(~XF[,2:3?   FIX VECTOR OF COLUMNS WITH TMZ AND AGE  ],model='FIXED'),
                  list(       do this   ))
   # Fitting the model
   fm<-BGLR(y=y[,4],a=y[,5],b=y[,6],ETA=ETA_cov,saveAt="Example_")
   # Retrieving estimates
   fm$coefficients
```


#### (8)  Evaluating prediction accuracy.

The following code illustrates how to use BGLR to calculate prediction accuracy of omic models. Note that in Bernal-Rubio et al 2018, AUC was calculated by estimating the proportion of deaths at different time points, equivalent to every month after diagnosis (e.g. month 1 after diagnosis, month 2 after diagnosis, ...). Below is an example for month 5 after the diagnosis (vital status available in y[, 4]).

The following illustrates how to select a validation set using the model `covariates+methylation` as example.

```R
#Installing and loading library pROC to compute Area Under the ROC Curve.
   install.packages(pkg='pROC')    # install pROC
   library(pROC);
   n <- length(y$ycen)
   status.5months <- y[,4]
# Setting training and testing sets
   trn <- (1:n)[folds == 1]
   tst <- -trn
   yNA <- y$status
   yNA[tst] <- NA
# Fit the model only in the training set
   fm.meth<-BGLR(y=y[,4],a=y[,5],b=y[,6],ETA=ETA.cov_meth,saveAt="Example_")
 # Find probability of survival for the testing set
   pred <-fm.meth$probs[tst,2]
# Estimate AUC
   AUC <- auc(status.5months[tst],pred)
#For the first individual, area under the standard normal curve (CDF) 
#of estimated y from full model:
   pnorm(fm.meth$yHat[1])
```

