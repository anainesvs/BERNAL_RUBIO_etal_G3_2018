## Integrating multiple Omics for Prediction of BC Survival

The following scripts illustrate how to fit some of the models presented in *Bernal Rubio et al., G3, 2018*, the scripts are also provided at: [https://github.com/anainesvs/BERNAL_RUBIO_etal_G3_2018](https://github.com/anainesvs/glioblastoma), please refer to that webpage for updates.

**Contact**: avazquez@msu.edu
**Contact**: jlbernalr@gmail.com

#### (1) Installing BGLR

The code below illustrates how to install and load the necessary package from CRAN using `install.packages()`.
```R
   install.packages(pkg='BGLR')    # install BGLR
   library(BGLR); 
```  

#### (2) Loading data
**Data**: The code assumes that the user has saved in the file `OMIC_DATA.rda` 
the objects that contain the phenotypic information, clinical covariates, and omic data. 
The code assumes that the file `OMIC_DATA.rda` contain the following objects:

   * `XF`: an incidence matrix for clinical covariates.
   * `Xge`: an incidence matrix for gene expression. 
   * `Xcnv`: an incidence matrix with CNV mean intensity summarized per gene. 
   * `Xmt`: an incidence matrix for methylation values at various sites.
   * `Xsnp`: an incidence matrix for SNP genotypes.
   * `y`: a vector with the response, in this case a 0/1 where 0 denotes alive.
   
The code below assumes that all the predictors were edited by removing outliers 
and predictors that did not vary in the sample, transformed if needed, and 
missing values were imputed. Also, incidence matrix should be conformable in dimensions for
subsequent analysis.

#### (3) Computing similarity matrices

 Some of the models fitted in the study use similarity matrices of the form G=XX' 
 computed from omics. The following code illustrates how to compute this matrix for 
 gene expression. A similar code could be use to compute a G-matrix for methylation 
 or other omics (see (6)).
 
 ```R 
  load('OMIC_DATA.rda')
  #Computing a similarity matrix for gene-expression data
   Xge<- scale(Xge, scale=true, center=TRUE) #centering and scaling
   Gge<-tcrossprod(Xge)                      #computing crossproductcts
   Gge<-Gge/mean(diag(Gge))                  #scales to an average diagonal value of 1.
```
**NOTE**: for larger data sets it may be more convinient to use the `getG()` function available in [BGData](https://github.com/quantgen/BGData) R-package. This function allows computing G without loading all the data in RAM and offers methods for multi-core computing. 

##### (4) Pre-adjust gene expression matrix by batch

In this section we show how to pre-correct the incidence matrix corresponding to gene expression by batch effects. The input used is the matrix "Xge" already centered, and the batch code for each sample. Using a linear mixed model, the pre-correction regresses each column of the Xge matrix into the corresponding batch effect, assuming them as random effects.

```R
#Load library lme4 to fit a linear mixed model
   library(lme4)
   adj_batch<-list()
   for(i in 1:ncol(Xge)){
      adj_batch[[i]]<-residuals(lmer(as.vector(Xge[,i])~1|batcheffect))
}

Xge_batch<-t(do.call(rbind,adj_batch))
rownames(Xge_batch)<-rownames(Xge)
colnames(Xge_batch)<-colnames(Xge)

```
The output "Xge_batch" can be used to compute the similarity matrix for gene expression according to section (3).

#### (5)  Exploring demographic effects in the principal comoponents (e.g., methylation-derived PC vs age at diagnosis)

The following code shows how search for effects of demographics sucha as age on the omic derived PCs, for instance Methylation.

```R
   #Gmt has been obtained as described in section (3)
   EVD<-eigen(Gmt,symmetric=TRUE)
   pc_mt<-EVD$vectors
   rownames(pc_mt)<-colnames(pc_mt)<-rownames(Gmt)

#First 10 PC derived from methylation by age at diagnosis
   npc<-10
   #Retrieve coefficient estimates, p-values and Rsquare 
   nresults<-3          
   OUT=matrix(nrow=10,ncol=3)
   colnames(OUT)=c('estimate','p-value','R2')
   rownames(OUT)=paste0('PC_',1:10)
   for(i in 1:10){
      fm=summary(lm(pc_mt[,i]~XF$age_diagnosis))
      OUT[i,1]=fm$coef[2,1]
      OUT[i,2]=fm$coef[2,4]
      OUT[i,3]=fm$r.squared
   }
   OUT
```

#### (6)  Fitting a binary regression for (the "fixed effects" of) Clinical Coavariates using BGLR (COV)

The following code illustrates how to use BGLR to fit a fixed effects model. The incidence matrix "XF" described in section (2) is used. There is no column for intercept in XF because BGLR adds the intercept automatically. The response variable `y` is assumed to be coded with two lables (e.g., 0/1), the argument `response_type` is used to indicate to BGLR that the response is ordinal (the binary case is a special case with only two levels). Predictors are given to BGLR in the form a two-level list. The argument `save_at` can be used to provide a path and a pre-fix to be added to the files saved by BGLR. For further details see [Pérez-Rodriguez and de los Campos, Genetics, 2014](http://www.genetics.org/content/genetics/198/2/483.full.pdf). The code also shows how to retrieve estimates of effects and of success probabilities. In the examples below we fit the model using the default number of iterations (1,500) and burn-in (500). In practice, longer chains are needed and therefore, the user can increase the number of iterations or the burn-in using the arguments `nIter` and `burnIn` of `BGLR`.

```R
### Inputs
# Centering and scaling the incidence matrix for fixed effects.
   XF<- scale(XF, scale=FALSE, center=TRUE) 
   ETA.COV<-list( COV=list(X=XF, model='FIXED') )
# Fitting the model
   fm=BGLR(y=y, ETA=ETA.COV, saveAt='cov_', response_type='ordinal')
# Retrieving estimates
   fm$ETA$COV$b      # posterior means of fixed effects
   fm$ETA$COV$SD.b   # posteriro SD of fixed effects
   head(fm$probs)    # estimated probabilities for the 0/1 outcomes.
```

#### (7)  Fitting a binary model for fixed effects and one omic (for instance, whole genome gene expression) using BGLR.

The following code illustrates how to use BGLR to fit a mixed effects model that accomodates both clinical covariates and whole-genome-gene expression. 

```R
# Setting the linear predictor to fit the model COV+GE
  ETA.COV.GE<-list( COV=list(X=XF, model='FIXED'), GE=list(K=Gge, model='RKHS'))
# Fitting the model
  fm.COV.GE<- BGLR(y=y, ETA=ETA.COV.GE, response_type='ordinal',saveAt='cov_ge_')
#  Retrieving predictors
  fm.COV.GE$mu            # intercept
  fm.COV.GE$ETA$COV$b     # effects of covariates
  fm$COV.GE$ETA$GE$varU   # variance associated to GE SD.varU gives posterior SD
  fm.COV.GE$ETA$GE$u      # random effects associated to gene expression
  plot(scan('cov_ge_ETA_GE_varU.dat'),type='o',col=4) # trace plot of variance of GE.
```
**NOTE**: to fit a similar model for covariates and a different omic, for instance COV+METH, the user just needs to change the inputs in the defintiion of the linear predictor by providing Gmt instead of Gge.

#### (8)  Fitting a binary model for fixed effects covariates and more than one omic (COV+GE+METH)

The following code shows how to extend the model `COV+GE` for inclusion of an additional omic dataset (methylation).

```R
# Computing a similarity matrix for methylation data
   Xmt<- scale(Xmt, scale=TRUE, center=TRUE)  #centering and scaling
   Gmt<-tcrossprod(Xmt)                       #computing crossproductcts
   Gmt<-Gmt/mean(diag(Gmt))                   #scales to an average diagonal value of 1.

ETA.COV.GE.MT<-list( COV=list(X=XF, model='FIXED'),
                     GE=list(K=Gge, model='RKHS'),
                     METH=list(K=Gmt, model='RKHS'))
# Fitting models 
fm.COV.GE.MT<- BGLR(y=y, ETA=ETA.COV.GE.MT, 
                 response_type='ordinal',saveAt='cov_ge_mt_')
```

#### (9) Validation

The following illustrates how to select a validation set using the model `COV` as example.

```R
#Installing and loading library pROC to compute Area Under the ROC Curve.
   install.packages(pkg='pROC')    # install pROC
   library(pROC);
   n <- length(y)
# Randomly select a 20% of the data to be the testing set 
   tst<- runif(n) <0.2
   yNA = y; yNA[tst] <-NA
# Fit the model only in the training set
   fm.COVtr<- BGLR(y=yNA, ETA=ETA.COV, response_type='ordinal')
# Find probability of survival for the testing set
   pred <-fm.COVtr$probs[tst,2]
# Estimate AUC
   AUC_train<-auc(y[!tst],fm.COVtr$yHat[!tst])
   AUC_test<-auc(y[tst], pred)
#For the first individual, area under the standard normal curve (CDF) 
#of estimated y from full model:
   pnorm(fm.COVtr$yHat[1])
```

