#Note: Make sure the genotype and phenotype file has same number of samples
# Filter Genotype file using Plink
## select.txt = list of subjectID/Sample to be selected
# plink --make-bed --bfile chr --fam chr.fam --out chr_select --keep select.txt
#####################################################################################
#Use this script for running multinomial GWAS:


library(snpStats) # for reading plink files
library(nnet) # for fitting the models

gen <- read.plink('chr_select')
phen <- read.table('all_pheno_mult.txt',header=T)

outfile = "nnet_output.txt"

for (i in 1:length(gen$genotypes[1,])){
    print(i)
   
    #fit the alternative model:
    md1 <- multinom(phen$ENDOTYPE ~ as.numeric(gen$genotypes[,i]))
   
    #fit the null model
    md0 <- multinom(phen$ENDOTYPE ~ 1)
   
    #do the LRT
    test <- anova(md1,md0)
   
    #write the output
    cat(paste(gen$map[i,],collapse="\t"),file=outfile,append=T)
    cat('\t',file=outfile,append=T)
    cat(paste(exp(coef(md1)[,2]),collapse="\t"),file=outfile,append=T)
    cat('\t',file=outfile,append=T)
    cat(paste(-0.5*test$"Resid. Dev",collapse="\t"),file=outfile,append=T)
    cat('\t',file=outfile,append=T)
    cat(test$P[2],file=outfile,append=T)
    cat("\n",file=outfile,append=T)
}

