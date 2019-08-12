library(dplyr)
library(lfmm)
library(RSpectra)
require(gdata)
library(vegan)
library(tidyverse)
library(LEA)
library(Hmisc)
library(robustHD)
#devtools::install_github("bcm-uga/LEA")
#################################################################
#                                                               #
#                  making pops                                  #
#                                                               #
#################################################################
setwd("~/Documents/UWS/GEA meta/repository/gea_code_test")
#load in file

############################################################
#
#   3 input files. Need to be exactly the right format. 
#
############################################################

#Meta data
        #emicr.meta = ./rawdata/JordanEmicroMetadataBySample.csv
        #marri.meta = ./rawdata/MarriMetadataBySample.csv
names.df.2 <- marri.meta   #read in file with sample column 1, pop column 3, lat longs columns 3 and 4, env data all other columns
n.df <- names.df.2[,c(1,2)]                                     #this effectively creates a file with just sample names and pops

#SNP data
        #emicr = ./rawdata/JordanEmicroSNP_012.tsv
        #marri = ./rawdata/marri.tsv
full.data.frame <- marri  #snp data; samples column 1; snps all remaining columns; we use the 012 format
full.data.frame[1:20,1:20]                                          #check the data frame, make sure it looks ok.

###############################################################################
###### input order for the choose.pop function
#1 pop dataframe (Labeled "Sample" and "Population")
#2 snp dataframe (sample labeled "V1")
#3 pop.size    #number of samples in a population
#4 numb   ## name of species; works best if keep to 5 characters, particularly in the baypass function (i.e. marri or emacu or yebox)
#5 maf.x  ## minor allele frequencies 
#6 prop   ## proportion of intact data (i.e. 0.9 = 10% missing data)
#7 env    ##  whole environmental data frame use the names.df.2 file
#8 bioclimvar   ## which bioclim variables to choose from the 19

#identify the parameters for all pops
minor.allele <- c(0.1)
prop.md <- c(0.8, 0.9)
ind.pop <- c(10, 6)
bio.var <- c("BIO1","BIO5", "BIO12", "BIO14")   #choose which world clim variables you want, needs to be all caps like this.
name.files <- "marri"                           #beginning of the name for all files. must be 5 characters.
#lfmm number of clusters
kclust <- 3  #need to run a pca as described in the manual to identify number of ancestral clusters for LFMM.

#for loop that uses the above parameters to clean the data and create many different data sets
#samples are randomly chosen
#pops, snps, pop.size, name (5 characters), maf.x, prop, env file, environments to choose - see line 50.
for (i in minor.allele) {
  for(j in prop.md) {
    for(k in ind.pop) {
      choose.pop(n.df, full.data.frame, k, name.files, i, j, names.df.2, bio.var) 
      #third is the label of the written file, and probably the species
    }
  }
}


#################################################################
#                                                               #
#                  lfmm                                         #
#                                                               #
#################################################################
dir.create('./1_lfmm2_data/')
pth <- "./1_lfmm2_data/" #specify the path holding the datafiles
#read in file names
filenames.lfmm.format <- list.files(pattern="*.txt")
for (i in 1:length(filenames.lfmm.format)) {
  genotype.file <- read.table(paste0(filenames.lfmm.format[i]),header = F, row.names=1)
  lfmm.file <- genotype.file[,-1]
  write.table(lfmm.file,paste0(pth,filenames.lfmm.format[i],".lfmm"), quote = FALSE, col.names=FALSE, row.names = FALSE)
}

filenames.lfmm <- list.files(path="./1_lfmm2_data/", pattern="*.lfmm")
prs.list <- list()  #create a list for the polygenic scores to be written into
df.out <- data.frame(matrix(vector(), 0, length(bio.var),  #this is the final table with all of the polygenic scores (with 4 environments)
                            dimnames=list(c(), c(unlist(bio.var)))), #
                     stringsAsFactors=F) 
dir.create('./1_lfmm2_data/calpval')
#run for loop
for (i in 1:length(filenames.lfmm)) {
  ##### load lfmm file #####
  x<-read.table(paste0(pth,filenames.lfmm[i]), header = FALSE)
  rows <- nrow(x)
  E <- read.table(paste0(name.files,"_",rows,".lfmm.env"), header = FALSE)  ### env files have to be named exactly right.
  #### impute #####
  
  for (j in 1:ncol(E)) {   #number of environments
    x <- impute(x, fun=median)
     #### Run LFMM2 #####
    mod.lfmm <- lfmm_ridge(Y = x, 
                           X = E[,j],
                           K = kclust) #chose K based on PCA
    ## performs association testing using the fitted model:
    pv <- lfmm_test(Y = x, 
                    X = E[,j], 
                    lfmm = mod.lfmm, 
                    calibrate = "gif")
    ## record the calibrated pvalues
    pvalues <- pv$calibrated.pvalue 
    write.table(pvalues, paste0(pth,"calpval/", filenames.lfmm[i],"_env_",j,".calpval"), quote = FALSE)
    #polygenic risk scores
    x.pred <- predict_lfmm(x,E[,1], mod.lfmm, fdr.level = 0.05, newdata = NULL)
    
    #perform lm on the predicted values vs env (i.e. can snps predict environment)
    x.lm <- lm(x.pred$pred ~ scale(E[,j], scale = FALSE))
    #calculate r2 for the polygenic risk score
    prs <- summary(x.lm)$adj.r.squared
    prs.list[[j]] <- prs #put all vectors in the list
  }
  df <- t(do.call("rbind",prs.list)) #throw the list of PRS into the dataframe
  #extract polygenic scores to dataframe
  df.out[i,] <- df
  write.table(df.out, paste0("./1_lfmm2_data/polygenicScores.pgs"), quote = FALSE)
}
paste("Finished LFMM")
#################################################################
#                                                               #
#                  baypass                                      #
#                                                               #
#################################################################
#(INPUTFILE,LONGLAT,BAYPASSFILE)
#baypass must be active in this folder (through $PATH or actual executable) 
#Need to create this directory: ./baypass_dir
dir.create('./2_baypass/')
filenames.bp.format <- list.files(pattern="*.txt")
for (g in 1:length(filenames.bp.format)) {
  genotype.file <- read.table(paste0(filenames.bp.format[g]),header = F, row.names=1)
  merged.gf <- merge(genotype.file, n.df, by.x = "row.names", by.y = "Sample")
  baypass.pops <- merged.gf %>%
    select(Row.names, Population) %>% 
    rename(Sample = Row.names)
  new.string <- str_remove(paste0(filenames.bp.format[g]), ".txt")
  baypass_format(paste0(filenames.bp.format[g]), baypass.pops, paste0("./2_baypass/",new.string,".baypass"))
}
paste("Finished writing baypass format")

# (POPNUM,DATAFILE,NTHREADS,NVAL,BURN,NPILOT,PILOTLENGTH,ENV))
filenames.baypass <- list.files("./2_baypass/", pattern="*.baypass")
for (f in 1:length(filenames.baypass)) {
  x<-read.table(paste0("./2_baypass/",filenames.baypass[f]), header = FALSE)
  n.pop <- ncol(x)/2
  n.ind<- substr(filenames.baypass[f], 7, 9)
  baypass(n.pop, paste0("./2_baypass/",filenames.baypass[f]), 1, 1000, 2500, 20, 1000, paste0(name.files,"_",n.ind,".baypass.env")) #variaables real low for testing
}

paste("Finished baypass")
#################################################################
#                                                               #
#                  RDA                                          #
#                                                               #
#################################################################
#parameters
n.permutations <-  99

n.cpus <- 4
## Setup working data and directory
dir.create('./3_rda_out')

## Data input/formatting and sanity checking


#Now make an allele frequency matrix.
filenames.rda.format <- list.files(pattern="*.txt")
for (i in 1:length(filenames.rda.format)) {
  rda.file <- read.table(paste0(filenames.rda.format[i]),header = F, row.names=1)
  names <- rownames(rda.file)
  merge.env <- merge(as.data.frame(names), names.df.2, by.x="names", by.y = "Sample" )
  af = apply(rda.file, 2, function(snp) {
    split(snp, merge.env$Population) %>%
      sapply(function(genos) {
        n = sum(!is.na(genos)) * 2
        s = sum(genos, na.rm=T)
        return (s / n)
      })
  })
  rda.string <- str_remove(paste0(filenames.rda.format[i]), ".txt")
  write.table(af, paste0("./3_rda_out/",rda.string,".xrda"), quote = FALSE, sep = "\t", col.names = FALSE)

}

dim(af)
#RDA requires complete dataframes (no NAs).  Check for NAs

n.na = sum(is.na(af))
if (n.na > 0) warning("RDA requires no missing data")
#Now, let's just look at the genetic pca to check our allele freq matrix looks OK.


#pc = prcomp(af, retx = T)
#pc$x %>%
#  as.data.frame() %>%
#  rownames_to_column("ID") %>%
#  ggplot(aes(PC1, PC2)) +
#  geom_text(aes(label=ID)) +
#  theme_bw()

files.rda <- list.files("./3_rda_out/", pattern="*.xrda")

for (h in 1:length(files.rda)) {
## RDA proper ###################################################################

# Run RDA (If you want to analyze factors, have to write out full equation
# - see Forester - otherwise can just use ~ . shorthand)
# Constrained axes will equal # vars in model
# Proportion Constrained = proportion of variance explained by EV predictors
# Most SNPs will be neutral, so low var explained is expected
# RDA - Borcard et al. (2011)
  
#extract the # of individuals to read in correct env data file. (although they should all be the same)
n.ind.rda<- substr(files.rda[h], 7, 9)

#Get environmental data
env.pop <- read.table(paste0(name.files,"_",n.ind.rda,".rda.env"),header = T, row.names = 1)
freq <- as.matrix(read.table(paste0("./3_rda_out/",files.rda[h]),header = F, row.names=1))
freq[is.na(freq)] <- 0 # this is just to run a small dataset.
euc.rda = rda(freq ~ ., data=env.pop, scale=T)
# Calculate adjusted R2 (R2 should be adjusted based on # of predictors)
RsquareAdj(euc.rda)
# Summarize variance explained by each canonical axis (eigenvalues)
summary(eigenvals(euc.rda, model="constrained"))

#Assess model significance
# Assess significance of full model w/ CCA F-statistics (Legendre, Oksanen, & terBraak 2010)
signif.full = anova.cca(euc.rda, parallel=n.cpus,
                        permutations=n.permutations)
signif.full


#And now we test the significance of each axis

# Assess significance of individual canonical RDA axes to determine which axes to use for investigating candidate loci
signif.axis <- anova.cca(euc.rda, by="axis", parallel=n.cpus,
                         permutations = n.permutations)
signif.axis

#Check Variance Inflation Factors (VIFs) among co-variates. We want VIF<10 or more conservatively <5 (Zuur et al. 2010).
vif.cca(euc.rda)


# Identify extreme outlier samples
snp_loadings =  scores(euc.rda, choices=c(1:3), display="species")
#use the outlier function from Forrester et al to identify the significant snps
cand1 <- outliers(snp_loadings[,1],3)
cand2 <- outliers(snp_loadings[,2],3)
cand3 <- outliers(snp_loadings[,3],3)

ncand <- length(cand1) + length(cand2) + length(cand3)

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

#combine all significant snps into one data frame, with loadings
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

#create input dataframe for each climate variable used
foo <- matrix(nrow=(ncand), ncol=ncol(env.pop))  # 8 columns for 8 predictors
colnames(foo) <- c(unlist(bio.var))

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- freq[,nam]
  foo[i,] <- apply(env.pop,2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  

ncol1 <- (ncol(cand)+1)
ncol2 <- (ncol(cand)+2)
#Change the numbers based on how many climate variables you use.
for (j in 1:length(cand$snp)) {
  bar <- cand[j,]
  cand[j,ncol1] <- names(which.max(abs(bar[4:(length(bio.var)+3)]))) # gives the variable
  cand[j,ncol2] <- max(abs(bar[4:(length(bio.var)+3)]))              # gives the correlation
}

#Change the numbers based on how many climate variables you use.
colnames(cand)[ncol1] <- "predictor"
colnames(cand)[ncol2] <- "correlation"

table(cand$predictor) 
write.table(cand, paste0("./3_rda_out/",files.rda[h],".out"), quote=FALSE)
}

#################################
#
#write.csv(snp_loadings, file.path("./rda_out", "snp_loadings.csv"))
#sample_loadings =  scores(euc.rda)$sites
#write.csv(sample_loadings, file.path("./rda_out", "sample_loadings.csv"))
#
#
#save(list = ls(), file=file.path("./rda_out", "data.Rdat"))
#
#################################
