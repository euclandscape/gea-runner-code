#####################################################
#   function
#####################################################


###############################
### choosing populations from a large data set
##############################
choose.pop <- function(pops, snps, pop.size, numb, maf.x, prop, env, bioclimvar){
  size.x <- pop.size
  count.pops<-pops %>% group_by(Population) %>% count(Population)
  keep<-count.pops %>% filter(n > size.x)
  indata <- subset(pops, Population %in% keep$Population)
  sub.10 <- indata %>%
    group_by(Population) %>%
    sample_n(size = size.x)
  outdata <- subset(pops, !Population %in% keep$Population)
  final <- rbind(as.data.frame(sub.10), as.data.frame(outdata))
  final.merge.snp <- merge(final, snps, by.x= "Sample", by.y = "V1")
  final.merge.snp <- final.merge.snp[,-2]
  #write.csv(final.merge.snp, paste0(numb, pop.size, ".csv"))
  out <- final.merge.snp[,-1]
  rownames(out) <- final.merge.snp[,1]
  
  #clean snps
  geno <- t(out)
  # look at counts across samples
  samplecalls<-apply(geno,1, function(x) sum(!is.na(x)));
  samp.thres <- prop*ncol(geno)
  
  g.snp <- geno[samplecalls > samp.thres,]
  ###### minor allele frequency ###########
  maf<-rowSums(g.snp,na.rm=T)/((ncol(g.snp)-rowSums(is.na(g.snp)))*2)
  #table(maf>mafThresh)
  maf<-g.snp[maf>maf.x,]
  #major allele freq
  jaf<-rowSums(maf,na.rm=T)/((ncol(g.snp)-rowSums(is.na(maf)))*2)
  jafThresh<-1-maf.x
  #table(jaf<jafThresh)
  snps<-maf[jaf<jafThresh,]
  ind <- ncol(snps)
  ####### create env files ##########
  r.names <- row.names(t(snps))
  df.r.names <- as.data.frame(r.names)
  merged.meta <- merge(df.r.names, env[,c(1,5:ncol(env))], by.x= "r.names", by.y = "Sample")
  st.env.meta <- standardize(merged.meta[,-1], centerFun = mean, scaleFun = sd)
  st.env.meta <- st.env.meta %>%
    select(bioclimvar)
  env.pop.rda = names.df.2 %>%
    select(-Sample) %>%
    unique() %>%
    remove_rownames %>%
    column_to_rownames(var = "Population")  %>%
    select(bioclimvar)  %>%
    scale(center=T, scale=T) %>%
    as.data.frame()
  ####### write table ###########
  write.table(t(snps), paste0(numb,"_",ind,"_",k,"_",j,"_",i,".txt"), quote = FALSE, col.names = FALSE)
  write.table(st.env.meta, paste0(numb,"_",ind,".lfmm.env"), col.names = F, row.names = F, quote = F) #lfmm env file
  write.table(t(env.pop.rda), paste0(numb,"_",ind,".baypass.env"), col.names = F, row.names = F, quote = F)
  write.table(env.pop.rda, paste0(numb,"_",ind,".rda.env"), quote = F)
}


###############################
### Baypass format and  outlier function
##############################
#Input file needs to be one row/line per individual, first column is id, then one column per loci, seperated by tab
#Output has pops in alphabetical order, need to make sure env file for baypass matches this

baypass_format<-function(INPUTFILE,LONGLAT,BAYPASSFILE){
  gen<-read.table(INPUTFILE, header = F, row.names=1)
  lonlat = LONGLAT
  allele1 = apply(gen, 2, function(snp) {
    split(snp, lonlat$Population) %>%
      sapply(function(genos) {
        n = sum(!is.na(genos)) * 2
        s = sum(genos, na.rm=T)
        return (s)
      })
  })
  allelecount = apply(gen, 2, function(snp) {
    split(snp, lonlat$Population) %>%
      sapply(function(genos) {
        n = sum(!is.na(genos)) * 2
        s = sum(genos, na.rm=T)
        return (n)
      })
  })
  allelefreq = apply(gen, 2, function(snp) {
    split(snp, lonlat$Population) %>%
      sapply(function(genos) {
        n = sum(!is.na(genos)) * 2
        s = sum(genos, na.rm=T)
        return (s / n)
      })
  })
  allele2<-allelecount-allele1
  whole<-interleave(allele1,allele2)
  baypass<-t(whole)
  write.table(baypass,BAYPASSFILE, col.names = F, row.names = F, sep="\t")
}


baypass<-function(POPNUM,DATAFILE,NTHREADS,NVAL,BURN,NPILOT,PILOTLENGTH,ENV){
  #Run baypass on core model to get covariance matrix, repeat 4 times
  SUB <- substr(filenames.baypass[f], 1, 9) #extract name species (5 letters) _ number of individuals (3 characters))
  name <- paste0(SUB,"_",f)
  for (i in 1:4) {
    cmd<-paste("./baypass -npop ",POPNUM," -gfile ",DATAFILE," -outprefix ./2_baypass/",name,"_",i," -seed $RANDOM -nthreads ",NTHREADS," -nval ",NVAL," -burnin ",BURN," -npilot ",NPILOT," -pilotlength ",PILOTLENGTH, sep="")
    system(cmd)
  }
  
  ###########################################################
  ###Calculate means from initial baypass runs on core model and produce   covariance matrix for initital core run with all SNPs
  myList<- vector(mode = "list", length = 4)
  
  for (i in 1:4) {
    myList[[i]]<-as.matrix(read.table(paste("./2_baypass/",name,"_",i,"_mat_omega.out",sep="")))
  }
  #mean omega matrix
  bfmean <- sapply(1:ncol(myList[[1]]), function(j) {apply(do.call(cbind,lapply(myList,`[`,,j)), 1, mean)})
  
  write.table(bfmean, file=paste0("./2_baypass/",name,"_mean_mat_omega.out"), row.names=FALSE, col.names=FALSE, sep="\t")
  
  #######################################################################################################################
  ###Run aux model of baypass for assocaitions using created covairance matrix 
  MATRIXFILE=paste0("./2_baypass/",name,"_mean_mat_omega.out")
  CORROUT=paste0("./2_baypass/",name)
  
  #aux model all snps
  for (i in 1:2) {
    cmd<-paste("./baypass -auxmodel -npop ",POPNUM," -gfile ./2_baypass/",DATAFILE," -efile ",ENV,"-outprefix ",CORROUT,"_",i," -omegafile ",MATRIXFILE," -seed $RANDOM -nthreads ",NTHREADS," -nval ",NVAL," -burnin ",BURN," -npilot ",NPILOT," -pilotlength ",PILOTLENGTH, sep="")
    system(cmd)
  }
  
  #######################################################################################################################
  ###Calculate means for aux runs and output a csv file with only those above BF20
  
  myList<- vector(mode = "list", length = 2)

  for (i in 1:2) {
    myList[[i]]<-as.matrix(read.table(paste("./2_baypass/",name,"_",i,"_summary_betai.out",sep=""), header=T))
  }
  
  piXtXmean <- sapply(1:ncol(myList[[1]]), function(j) {apply(do.call(cbind,lapply(myList,`[`,,j)), 1, mean)})
  colnames(piXtXmean) <-colnames(myList[[1]])
  
  write.table(piXtXmean, file=paste0("./2_baypass/mean_summary_betai_",name,".out"), row.names=FALSE, col.names=T, sep="\t", quote=F)
  
  newtable<-as.data.frame(piXtXmean)
  
  BF3 <- newtable[newtable$`BF.dB.` > 3,]
  write.csv(BF3, file=paste0("./2_baypass/BF3_",name,".csv"), row.names=F)
  
}


###############################
### RDA outlier function
##############################
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}
