
setwd('/Users/sbharathkumar/Documents/Models/CT_model')    #the directory with all the target files
library(ggplot2)
library(bigsnpr)

#Loading the bim file of target (target files are in PLINK format - .bed,.fam,.bim files)
bim <- read.table("target.bim")                                   #from reference of UK biobank, A1 is alternate allele,and A2 is reference allele
colnames(bim) <- c("chr", "snp", "GD", "pos", "a1", "a0")        #a1 is alternate allele, a0 is reference allele

# Read in the GWAS data and matching the SNPs
pain <- list()     #creating an empty list, which can load all the base files as data frames at each position(index)
base_snp <- list()

#the loop below loads base file, adds column names as required, perform SNP match with target data
#then, extract the required fields of base data from the obtained snp matched data 
#changing the order of columns to known order earlier
#modifying few names of columns to previous existing names(changing alleles back to A1 and A2 from a1 and a0 respectively)

for (i in 1:17){
  
  pain[[i]] <- read.table(paste0("/Users/sbharathkumar/Downloads/Pleio GWASs/base_pleio",i,".txt"),header = T,stringsAsFactors = F,sep="\t")   #reading all the base files 
  colnames(pain[[i]]) <- c("chr", "pos","snp","a1", "a0","se","beta","n","p")
  info_snp<- snp_match(pain[[i]],bim)  
  base_snp_matched <- info_snp[,(1:9)]
  base_snp[[i]] <- base_snp_matched[, c("chr","pos","snp.ss","a1","a0","se","beta","n","p")]    #order of columns as required                      
  colnames(base_snp[[i]]) <- c("chr","pos","snp","A1","A2","se","beta","n","p")                 ##renaming as required
  
}

names(bim)[names(bim) == "a1"] <- "A1"
names(bim)[names(bim) == "a0"] <- "A2"

#the loop below will create file for all the new transformed base data obtained from snp_match

for  (i in 1:17){
  write.table(base_snp[[i]],paste0("base",i,".transformed"),quote =FALSE,row.names = FALSE,col.names = TRUE,sep = '\t')  #creating a new base file obtained after matching SNPs
}

######################## execute this code in terminal###########################
#generating multiple threshold ranges to use them later for obtaining different profiles with different thresholds
echo "0.001 0 0.001" > range_list 
echo "0.05 0 0.05" >> range_list
echo "0.1 0 0.1" >> range_list
echo "0.2 0 0.2" >> range_list
echo "0.3 0 0.3" >> range_list
echo "0.4 0 0.4" >> range_list
echo "0.5 0 0.5" >> range_list

#loop performing clumping and threshold method
#the first plink command will form the clumps of SNPs in the target data
#it will generate target+base_number.clumped files consisting the data of clumps formed

for i in {1..17}
do
plink \
--bfile target \
--clump-p1 1 \
--clump-r2 0.1 \
--clump-kb 250 \
--clump base${i}.transformed \
--clump-snp-field snp \
--clump-field p \
--out base${i}

#extracting Valid SNP from the generated clumped file
awk 'NR!=1{print $3}' base${i}.clumped >  base${i}.valid.snp

#file with SNPs and P values from base data
awk '{print $3,$9}' base${i}.transformed > SNP.pvalue.base${i}

#generating PRS score 
plink \
--bfile target \
--score base${i}.transformed 3 4 7 header \
--q-score-range range_list SNP.pvalue.base${i} \
--extract base${i}.valid.snp \
--out base${i}

done

#the prs score is generated for different threshold values present in the range_list for each of the base taken.
#So, for each base, there will be 7 files generated in format - target+base_number.threshold.profile

#population stratification
plink \  
--bfile target \
--indep-pairwise 200 50 0.25 \
--out target

#obtaining principal components, taken 3 because of very small dataset
#checked by changing different pca values (6,5,4,3,2) and chosen the best value as 3, which is giving better results compared to others

plink \
--bfile target \
--extract target.prune.in \
--pca 3 \ 
--out target

#################################### in R ############################################

#modeling the best fit PRS from the data

p.threshold <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)

# Reading the phenotype file 
data <- read.table("target.fam", header = FALSE)
phenotype <- data[, c(1,2,5,6)]                           # Extracting 5th and 6th columns and assign column names (which are the sex and pain)
colnames(phenotype) <- c("FID", "IID","sex", "pain")

# Reading the PCs
pcs <- read.table("target.eigenvec", header=F)
colnames(pcs) <- c("FID", "IID", paste0("PC",1:3))        # (1:3 because there are 3 PCs)
pheno <- merge(phenotype, pcs, by=c("FID","IID"))         # Now merging all the phenotype data

######### Null model, model using only the covariates(sex here), principal components(PCs) ################
# We then calculate the null model using a linear regression 
null.model <- lm(pain~., data=pheno[,!colnames(pheno)%in%c("FID","IID")])
null.r2 <- summary(null.model)$r.squared             # R2 of the null model

######### Linear model, model using only the covariates(sex here), principal components(PCs), and PRS Score################

#this loop will iterate through all the profiles generated for each base.
#each profile is for a p value threshold for different bases.
#model is built for each profile, the results are noted for all the profiles of a base in a file base_number.target_type.prs.result 

prs.result<-list()

for (i in 1:17){
  
  prs.result<- NULL
  
  for (j in p.threshold){
    
    prs <- read.table(paste0("base",i,".",j,".profile"), header=T)
    pheno.prs <- merge(pheno, prs[,c("FID","IID", "SCORE")], by=c("FID", "IID"))
    model <- lm(opioid~., data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")])
    model.r2 <- summary(model)$r.squared              # model R2 is obtained
    prs.r2 <- model.r2-null.r2                        # R2 of PRS is simply calculated as the model R2 minus the null R2 
    prs.coef <- summary(model)$coeff["SCORE",]        # We can also obtain the coeffcient and p-value of association of PRS as follow
    prs.beta <- as.numeric(prs.coef[1])
    prs.se <- as.numeric(prs.coef[2])
    prs.p <- as.numeric(prs.coef[4])
    
    # We can then store the results in the prs.result table
    prs.result <- rbind(prs.result, data.frame(Threshold=j, R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se))
    
  }
  write.table(prs.result,paste0('base.',i,'.prs.result'),sep="\t",row.names = FALSE)
}

#the loop below is used to generate the plots for each base data chosen.

for (i in 1:17){
  prs.result <- read.table(paste0('base.',i,'.prs.result'),header =TRUE)
  prs.result$print.p <- round(prs.result$P, digits = 3)
  prs.result$print.p[!is.na(prs.result$print.p) & prs.result$print.p == 0] <- format(prs.result$P[!is.na(prs.result$print.p) & prs.result$print.p == 0], digits = 2)
  prs.result$print.p <- sub("e", "*x*10^", prs.result$print.p)
  
  ggplot(data = prs.result, aes(x = factor(Threshold), y = R2)) +
    # Specifying we want to print p-value on top of the bars
    geom_text(
      aes(label = paste('P:',print.p)),
      vjust = -1.5,
      hjust = 0,
      angle = 45,
      cex = 4,
      parse = T
    )  +
    # Specifying the range of the plot, *1.25 to provide enough space for the p-values
    scale_y_continuous(limits = c(0, max(prs.result$R2) * 1.25)) +
    # Specifying the axis labels
    xlab(expression(italic(P) - value ~ threshold ~ (italic(P)[T]))) +
    ylab(expression(paste("PRS model fit:  ", R ^ 2))) +
    # Drawing a bar plot
    geom_bar(aes(fill = -log10(P)), stat = "identity") +
    # Specifying the colors
    scale_fill_gradient2(
      low = "dodgerblue",
      high = "firebrick",
      mid = "dodgerblue",
      midpoint = 1e-4,
      name = bquote(atop(-log[10] ~ model, italic(P) - value),)
    ) +
    # Some beautification of the plot
    theme_classic() + theme(
      axis.title = element_text(face = "bold", size = 18),
      axis.text = element_text(size = 14),
      legend.title = element_text(face = "bold", size =
                                    18),
      legend.text = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust =
                                   1)
    )
  
  ggsave(paste0("base.",i,".bar.png"), height = 7, width = 7)
  
}
