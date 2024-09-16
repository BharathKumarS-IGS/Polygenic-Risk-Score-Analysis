setwd("/Users/sbharathkumar/Documents/Models/LDpred2_model")

library(bigsnpr)
library(bigreadr)
library(data.table)
(NCORES <- nb_cores()) 

map_ldref <- readRDS("map_hm3_plus.rds")     #I have downloaded the file from the site: "https://figshare.com/ndownloader/files/37802721"

snp_readBed("target.bed")                   # Read from bed/bim/fam, it generates .bk and .rds files.
obj.bigSNP <- snp_attach("target.rds")
map_test <- setNames(obj.bigSNP$map[-3], c("chr", "rsid","pos", "a1", "a0"))
G <- obj.bigSNP$genotypes                   #genotypes we use later in 'big_prodVec' to calculate the prs score

#######################################building the model################################################
# Read in the phenotype file 
data <- read.table("target.fam", header = FALSE)
phenotype <- data[, c(1,2,5,6)]                           # Extracting 5th and 6th columns and assign column names (which are the sex and pain)
colnames(phenotype) <- c("FID", "IID","sex", "pain")

# Read in the PCs
pcs <- read.table("target.eigenvec", header=F)
colnames(pcs) <- c("FID", "IID", paste0("PC",1:3))        # (1:3 because there are 3 PCs)

pheno <- merge(phenotype, pcs, by=c("FID","IID"))         # Now merging all the phenotype data

######### Null model, model using only the covariates(sex here), principal components(PCs) ################
# We then calculate the null model using a linear regression 
null.model <- lm(pain~., data=pheno[,!colnames(pheno)%in%c("FID","IID")])
null.r2 <- summary(null.model)$r.squared             # R2 of the null model
###########################################################################################################
result<-NULL     #to store the R2 of PRS for 17 traits

#loop to go over through all the traits 
for (i in 1:17){
  #loading summary statistics
  sumstats <- fread2(paste0("/Users/sbharathkumar/Downloads/Pleio GWASs/base_pleio",i,".txt"),col.names = c("chr", "pos","snp","a1", "a0","beta_se","beta","n_eff","p"))
  
  #performing SNP matching between the variants of summary statistics and the LD reference 
  info_snp <- snp_match(sumstats, map_ldref)
  (info_snp <- tidyr::drop_na(tibble::as_tibble(info_snp)))    #dropping any NA's in the data
  
  in_test <- vctrs::vec_in(info_snp[, c("chr", "pos")], map_test[, c("chr", "pos")])   #filtering the SNPs(matching) with target data
  df_beta <- info_snp[in_test, ]
  
  tmp <- tempfile(tmpdir = "tmp-data")
  
  #loop to get the LD correlation matrix 
  for (chr in 1:22) {
    
    cat(chr, ".. ", sep = "")
    ## indices in 'df_beta'
    ind.chr <- which(df_beta$chr == chr)
    ## indices in 'map_ldref'
    ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
    ## indices in 'corr_chr'
    ind.chr3 <- match(ind.chr2, which(map_ldref$chr == chr))
    ind.chr3 <- ind.chr3[!is.na(ind.chr3)]  # Remove NAs
    if (length(ind.chr3) > 0) {  # Ensure non-empty indices
      corr_chr <- readRDS(paste0("/Users/sbharathkumar/Downloads/ldref_hm3_plus/LD_with_blocks_chr", chr, ".rds"))[ind.chr3, ind.chr3]
      
      if (chr == 1) {
        corr <- as_SFBM(corr_chr, tmp, compact = TRUE)
      } else {
        corr$add_columns(corr_chr, nrow(corr))
      }
    }
  }
  
  #ld_score regression to obtain the h2 estimate for the trait
  (ldsc <- with(df_beta, snp_ldsc(ld, ld_size = nrow(map_ldref),
                                  chi2 = (beta / beta_se)^2,
                                  sample_size = n_eff,
                                  ncores = NCORES)))
  
  h2_est <- ldsc[["h2"]]   
  
  print(paste0("The heritability estimate of trait ",i," is: ",h2_est))
  
  #performing the LDpred2 auto
  multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                                 vec_p_init = seq_log(1e-4, 0.2, length.out = 30),
                                 allow_jump_sign = FALSE, shrink_corr = 0.95,
                                 ncores = NCORES)
  
  
  (range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est))))     #gives how much the correlation estimates vary (their range) in each item of the list 'multi_auto'
  (keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE))))        #contains the indices of the elements in range that are greater than 95% of the 95th percentile value. This effectively selects the top values within the range that are above this threshold.
  
  beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))  #beta_auto : vector containing the average effect size estimates across the chosen length of list(30). 
  
  #calculating the prs score 
  pred_auto <- big_prodVec(G, beta_auto,ind.row = rows_along(G),ind.col = 1:length(beta_auto), ncores = NCORES)
  
  reg.dat     <- pheno
  reg.dat$PRS <- pred_auto      #adding prs column to the data
  
  #sum(reg.dat$PRS>0)       #to see how many individuals have positive prs score, or negative if <
  
  auto.model <- lm(opioid~., data = reg.dat[,!colnames(reg.dat)%in%c("FID","IID")])   #linear model along with PRS
  result[i] <- data.table(auto = summary(auto.model)$r.squared - null.r2,null = null.r2)  #R2 of PRS
  
  print(paste0("result of trait ",i," is ",result[i]))   
  
  file.remove(paste0(tmp, ".sbk"))      #removing the temporary file created
}

write.table(result,"result_ldpred2.txt") #to save the results


###############################for visualizing the R2 for all traits###############################
library(ggplot2)
library(dplyr)
library(data.table)
library(readxl)

# Loading the data from the excel containing the trait name for each trait number
pain_data <- read_excel("/Users/sbharathkumar/Downloads/Pleio GWASs/codes_correspondence.xlsx", col_names = FALSE, sheet = 1)

# adding r2 score to the data frame
pain_data$r2 <- unlist(result)

# Renaming columns
colnames(pain_data) <- c("trait_number", "pain_type", "r2_score")

# Finding the indices of the top 3 bars with the highest R2 scores, to make the interpretation easier
top_3_indices <- order(pain_data$r2_score, decreasing = TRUE)[1:3]

# Creating the bar plot
ggplot(pain_data, aes(x = pain_type, y = r2_score)) +
  geom_bar(stat = 'identity', aes(fill = factor(ifelse(1:nrow(pain_data) %in% top_3_indices, "top", "other")))) +
  scale_fill_manual(values = c("top" = "darkblue", "other" = "skyblue"), guide = FALSE) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(size = 0.5, color = "gray"),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "Traits", y = "R2 score", title = "Bar Plot of Traits and their R2")
