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

#### (2) Simulating survival data, demographics and omic information

The following script allows to simulate

```R
  #Simulating a SNP matrix.
  p <- runif(1)
  X.snp <- matrix(sample(0:2, prob = c(p^2, 2*p*(1-p), (1-p)^2), replace = T,size=100*1000),100,1000)
  X.snp[1:4,1:3]

  #Simulating GE matrix (w/counts).
  X.ge <- matrix(rpois(n = 100*300, lambda = 50),100,300)
  X.ge[1:4,1:3]

  #Simulating DM matrix (beta values).
  X.meth <- matrix(rbeta(100*2000, shape1 = 0.5, shape2 = 0.5),100,2000)
  X.meth[1:4,1:3]

  #Simulating survival data.
  library(simsurv)
  covs <- as.data.frame(cbind(X.ge[,sample(1:300,4)],
                            X.meth[,sample(1:2000,10)],
                            X.snp[,sample(1:300,5)]))
  s.fm <- simsurv(lambdas = 0.08, gammas = 5,x = covs, maxt = 2)

  #Simulating gender.
  X.cov <- model.matrix(~sample(factor(c("female","male")), prob=c(p,1-p),replace=T,size=100))
  gender<-X.cov[,2]
  s.fm$gender<-ifelse(gender==0,"female","male")

  #Simulating batch codes
  s.fm$batch_code<-sample(50,100,replace=TRUE)

  #Save objects editing names of rows and columns.
  XF<-s.fm[,c(1,4:5)]
  y<-s.fm[,c(1,2,3)]
  
  y$ycen <- y$eventtime
  n <- nrow(y)
  y$b <- rep(Inf,n)
  y$a <- rep(-Inf,n)
  isCensored <- y$status==0
  y$a[isCensored] <- y$ycen[isCensored]
  y$ycen[isCensored] <- NA
  
  head(y)

  rownames(X.snp)<-rownames(X.ge)<-rownames(X.meth)<-as.character(s.fm$id)
  colnames(X.meth)<-c(paste0("CPG_",seq(1:ncol(X.meth))))
  colnames(X.ge)<-c(paste0("ge_",seq(1:ncol(X.ge))))
  colnames(X.snp)<-c(paste0("SNP_",seq(1:ncol(X.snp))))

  save(X.ge,X.meth,X.snp,y,XF,file="OMIC_DATA.rda")
``` 

#### (3) Loading data
**Data**: The code assumes that the user has downloaded the file `OMIC_DATA.rda`, or executed the previous simulation script which contains:

   * `XF`: a matrix for clinical covariates (age).
   * `X.ge`: a matrix for gene expression. 
   * `X.meth`: a matrix for methylation values at various sites.
   * `X.snp`: a matrix of Single Nucleotide Polymorphism (SNP) genotypes.
   * `y`: a vector with the response variable (time to death 0/1 where 0 denotes alive) and variables `a` and `b` required by BGLR to analyze censored data. More description on this analysis on BGLR can be found on [https://cran.r-project.org/web/packages/BGLR/vignettes/BGLR-extdoc.pdf](https://cran.r-project.org/web/packages/BGLR/vignettes/BGLR-extdoc.pdf)
  
After downloading, files can be loaded into R with:

```R 
  load('OMIC_DATA.rda')
  #Checking files present and dimensions
  ls()
  #[1] "X.ge"   "X.meth" "X.snp"  "XF"     "y"   
  
```

Clinical covariates, as well as omic datasets should be edited by removing outliers (for instance, observations outside of mean $\pm$ 3SD) and ordering samples id equally across datasets. Also, incidence matrix should be conformable in dimensions for subsequent analysis.

#### (4) Description of SNP data

Within the omic datasets used in this code, we describe briefly the structure of SNP file. SNP genotypes correspond to a matrix of sample ids in rows and SNP names in columns. The values included in the matrix (discrete values of 0,1,2 for observed genotypes or continuous values between 0 and 2 for imputed genotypes) represent the number of copies of the reference allele observed or imputed for an individual in a particular SNP. For instance:

```R 
   X.snp[1:5,1:5]
   #  SNP_1 SNP_2 SNP_3 SNP_4 SNP_5
   #1     1     0     0     0     2
   #2     1     0     1     0     0 
   #3     0     0     0     2     0
   #4     1     0     0     1     1
   #5     0     2     0     0     0
```

#### (5) Computing similarity matrices

 Some of the models fitted in the study use similarity matrices of the form G=XX' 
 computed from omics. The following code illustrates how to compute this matrix for 
 SNP genotypes. A similar code could be use to compute a G-matrix for methylation 
 or other omics (see (6)).
 
```R 
   #Computing a similarity matrix for SNP genotypes
   Xsnp_sc<- scale(X.snp, scale=TRUE, center=TRUE) #centering and scaling
   Gsnp<-tcrossprod(Xsnp_sc)                       #computing crossproductcts
   Gsnp<-Gsnp/mean(diag(Gsnp))                     #scales to an average diagonal value of 1.
   
   #Computing a similarity matrix for GE genotypes
   Xge_sc<- scale(X.ge, scale=TRUE, center=TRUE) #centering and scaling
   Gge<-tcrossprod(Xge_sc)                       #computing crossproductcts
   Gge<-Gge/mean(diag(Gge))                      #scales to an average diagonal value of 1.
```
**NOTE**: for larger data sets it may be more convinient to use the `getG()` function available in [BGData](https://github.com/quantgen/BGData) R-package. This function allows computing G without loading all the data in RAM and offers methods for multi-core computing. 

#### (6)  Exploring demographic effects in the principal components (e.g., SNP-derived PC vs gender)

The following code shows how to search for effects of demographics such as gender on the omic derived PCs, for instance from SNP genotypes

```R
   #Gsnp has been obtained as described in previous section
   EVD<-eigen(Gsnp,symmetric=TRUE)
   pc_snp<-EVD$vectors
   rownames(pc_snp)<-colnames(pc_snp)<-rownames(Gsnp)

```

Principal component coefficients can be plotted to visualize stratification in the dataset by gender as follows:


```R
   library(ggplot2)
   #Plot of SNP-derived PCs by gender, for instance PC1 vs PC2
   temp<-as.data.frame(pc_snp)
   temp$gender<-ifelse(XF$gender==0,"female","male")
   ggplot(temp, aes(temp[,1], temp[,2], col=gender)) + 
   geom_point(aes(shape=gender), size=2) +   # draw points
   labs(title="SNP Clustering", 
       subtitle="PC1 and PC2 by gender")
```

A similar approach can be implemented for additional omics (gene expression and methylation).

#### (7) Variance explained by omics

For comparison of variance explained between omics, a linear regression of omic-derived PC on demographics can be done. For instance, the following script allows to obtain the estimates, p-values and variance explained (given by R-squared) by the first 10 principal components derived from SNP genotypes. 

```R

#First 10 PC derived from SNP genotypes by gender
   npc<-10
   #Retrieve coefficient estimates, p-values and Rsquare 
   nresults<-3          
   OUT=matrix(nrow=10,ncol=3)
   colnames(OUT)=c('estimate','p-value','R2')
   rownames(OUT)=paste0('PC_',1:10)
   for(i in 1:10){
      fm=summary(lm(pc_snp[,i]~XF$gender))
      OUT[i,1]=fm$coef[2,1]
      OUT[i,2]=fm$coef[2,4]
      OUT[i,3]=fm$r.squared
   }
   OUT
```
The PCs estimates, p-values and R-squared will be printed (OUT object). The variance explained can be visualized as follows: 


```R
   temp<-as.data.frame(OUT)
   temp$PC<-as.character(rownames(temp))
   ggplot(temp, aes(x=PC, y=R2)) + 
   geom_bar(stat="identity", width=.5, fill="tomato3") + 
   labs(title="Variance explained by PC", 
       subtitle="First 10 PC SNP-derived") 
```
       
##### (8) Pre-adjust gene expression matrix by batch

In this section we show how to pre-correct the incidence matrix corresponding to gene expression by batch effects. The input used is the matrix "X.ge" centered, and the batch code for each sample. Using a linear mixed model, the pre-correction regresses each column of the Xge matrix into the corresponding batch effect, assuming them as random effects.

```R
  Xge<- scale(X.ge, scale=TRUE, center=FALSE) 

  #Load library lme4 to fit a linear mixed model
   library(lme4)
   adj_batch<-list()
   for(i in 1:ncol(Xge)){
      adj_batch[[i]]<-residuals(lmer(as.vector(Xge[,i])~1|XF$batch_code))
   }

   Xge_batch<-t(do.call(rbind,adj_batch))
   rownames(Xge_batch)<-rownames(Xge)
   colnames(Xge_batch)<-colnames(Xge)

```
The output "Xge_batch" can be used to compute the similarity matrix for gene expression according to section (5).


#### (9)  Fitting a binary regression for (the "fixed effects" of) Clinical Covariates or demographics using BGLR.

The following code illustrates how to use BGLR to fit a fixed effects model. The incidence matrix "XF" as well as the object "y" described in section (3) are used. There is no column for intercept in XF because BGLR adds the intercept automatically. Predictors are given to BGLR in the form a two-level list. The argument `save_at` can be used to provide a path and a pre-fix to be added to the files saved by BGLR. For further details see [PÃ©rez-Rodriguez and de los Campos, Genetics, 2014](http://www.genetics.org/content/genetics/198/2/483.full.pdf). The code also shows how to retrieve estimates of effects and of success probabilities. In the examples below we fit the model using the default number of iterations (1,500) and burn-in (500). In practice, longer chains are needed and therefore, the user can increase the number of iterations or the burn-in using the arguments `nIter` and `burnIn` of `BGLR`.

```R
   ### Inputs: XF and y
   ETA_cov<-list(list(~XF[,2],model='FIXED'))
   # Fitting the model
   fm<-BGLR(y=y[,4],a=y[,5],b=y[,6],ETA=ETA_cov,saveAt="Example_")
   # Retrieving estimates
   fm$coefficients
```

#### (10)  Fitting a binary model for fixed effects and one omic (for instance, whole genome gene expression) using BGLR.

The following code illustrates how to use BGLR to fit a mixed effects model that accomodates both demographics (or clinical if available) covariates and whole-genome-gene expression. The input related to gene expression corresponds to the similarity matrix based on this omic "Gge", which can be obtained considering step (5).

```R
   # Setting the linear predictor to fit the model gender+Gene expression
   ETA_cov_ge<-list(list(X=model.matrix(~XF[,c('gender')])[,-1],model='FIXED'),
	              list(K=Gge,model='RKHS',saveEffects=TRUE))
   # Fitting the model
   fm<-BGLR(y=y[,4],a=y[,5],b=y[,6],ETA=ETA_cov_ge,saveAt="Example_")
   # Retrieving estimates
   fm$mu                   # intercept
   fm$ETA[[1]]$b           # effects of covariates
   fm$ETA[[2]]$varU        # variance associated to GE SD.varU gives posterior SD
   fm$ETA[[2]]$u           # random effects associated to gene expression
   plot(scan('Example_ETA_2_varU.dat'),type='o',col=4) # trace plot of variance of GE.
```

**NOTE**: to fit a similar model for covariates and a different omic, for instance gender+Methylation, the user just needs to change the inputs in the definition of the linear predictor by providing the similarity matrix for methylation instead of Gge.

#### (11)  Fitting a binary model for fixed effects covariates and more than one omic (COV+GE+METH)

The following code shows how to extend the previous model `COV+GE` for inclusion of an additional omic dataset (methylation).

```R
# Computing a similarity matrix for methylation data following step (5)
   Xmt<- scale(X.meth, scale=TRUE, center=TRUE)  #centering and scaling
   Gmt<-tcrossprod(Xmt)                          #computing crossproducts
   Gmt<-Gmt/mean(diag(Gmt))                      #scales to an average diagonal value of 1.
   

   # Setting the linear predictor to fit the model gender+Gene expression + methylation
   ETA_cov_ge_mt <-list(list(X=model.matrix(~XF[,c('gender')])[,-1],model='FIXED'),
	              list(K=Gge,model='RKHS',saveEffects=TRUE),
	              list(K=Gmt,model='RKHS',saveEffects=TRUE))
   # Fitting the model
   fm<-BGLR(y=y[,4],a=y[,5],b=y[,6],ETA=ETA_cov_ge_mt,saveAt="Example_")

```

#### (12) Validation

The following illustrates how to select a validation set using the model `Gender+Gene expression+ methylation` as example.

```R
#Installing and loading library pROC to compute Area Under the ROC Curve.
   install.packages(pkg='pROC')    # install pROC
   library(pROC);
   n <- length(y$ycen)
# Randomly select a 20% of the data to be the testing set 
   sample(which(y$status==1),)
   tst<- (1:n) %in% c(sample(which(y$status==0),size=ceiling(0.2*sum(y$status==0))),
           sample(which(y$status==1),size=ceiling(0.2*sum(y$status==1))))
   yNA = y$status; yNA[tst] <-NA
# Setting the linear predictor to fit the model gender+Gene+Methylation expression
   ETA_cov_ge_mt <-list(list(X=model.matrix(~XF[,c('gender')])[,-1],model='FIXED'),
	              list(K=Gge,model='RKHS',saveEffects=TRUE),
	              list(K=Gmt,model='RKHS',saveEffects=TRUE))
# Fit the model only in the training set
   fm.tr<- BGLR(y=yNA, ETA=ETA_cov_ge_mt, response_type='ordinal')
# Find probability of survival for the testing set
   pred <-fm.tr$probs[tst,2]
# Estimate AUC
   AUC <- auc(y$status[tst],pred)
#For the first individual, area under the standard normal curve (CDF) 
#of estimated y from full model:
   pnorm(fm.tr$yHat[1])
```

