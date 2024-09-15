# Polygenic Risk Score (PRS) Analysis for Pain-Related Traits

This repository contains R scripts and workflows used for **Polygenic Risk Score (PRS)** analysis on pain-related traits for a cohort of 150 individuals who underwent knee surgery. The project explores three methods: **Clumping + Thresholding (C+T)**, **LDpred2**, and **PRSCS** to compute and compare the variance explained by PRS in different pain traits.

## Project Overview

### 1. Literature Review and Familiarization
In the first stage, we reviewed literature on **Genome-Wide Association Studies (GWAS)**, **quality control (QC)** steps, and **Polygenic Risk Scores (PRS)**. To understand the methodologies, we followed the [PRS tutorial](https://choishingwan.github.io/PRS-Tutorial).

### 2. Data Description
- **Target Data**: A cohort of 150 individuals who underwent knee surgery, with the target phenotype being **pain** and the covariate as **sex**.
- **Base Data**: Summary statistics for 17 different pain-related traits (e.g., chest pain, headache, hip pain) and medications (e.g., ibuprofen, paracetamol).

### 3. Methodology

#### Clumping + Thresholding (C+T)
This method selects genetic variants based on their statistical significance (P-value). The variants are “clumped” to ensure they are independent, and thresholds are applied to filter out variants with low significance. We tested different P-value thresholds to identify the best PRS.

#### LDpred2
A Bayesian method that adjusts for the correlation (linkage disequilibrium) between genetic variants. It computes a more accurate PRS by incorporating information from all variants, even those with low effect sizes.

#### PRSCS
Another Bayesian approach that improves on LDpred2 by allowing us to specify certain parameters. It shrinks the effect sizes of many genetic variants to zero while keeping the most important ones, which can enhance the accuracy of PRS.

### 4. Data Analysis
For each method, we computed PRS and incorporated it into a linear regression model to assess the variance explained (R2). We compared the R2 values for:
1. Models without PRS
2. Models with PRS

This allowed us to evaluate the contribution of PRS to the prediction of pain traits.

### 5. Visualization and Results
- For **C+T**, bar plots were generated to visualize the maximum R2 value across different **P-value thresholds** (e.g., 0.001, 0.01, 0.05, 0.1).
- For **LDpred2** and **PRSCS**, we plotted R2 values for each trait.

#### Key Results:
- **C+T**: Top traits—Chest pain, chronic headache, and hip pain.
- **LDpred2**: Top traits—Chronic back pain, hip pain, neck shoulder.
- **PRSCS**: Top traits—Hip pain, neck shoulder, chronic headache.

### 6. Future Work
We aim to:
- Explore PRS performance using other base summary statistics (e.g., depression or opioid use).
- Develop a unified model to generate a **single PRS value** across multiple traits.

## Repository Structure
- **`scripts/`**: R scripts used for PRS computation, data analysis, and visualization.
    - `PRS_CT_method.R`: Clumping + Thresholding method implementation.
    - `PRS_LDpred2_method.R`: LDpred2 method implementation.
    - `PRS_PRSCS_method.R`: PRSCS method implementation.
- **`results/`**: Output files with R2 values and visualizations.
  
## Learnings
Throughout this project, I gained valuable insights and skills, including:
- Understanding **Genome-Wide Association Studies (GWAS)** and **quality control (QC)** steps in genomic research.
- Gaining proficiency in **Polygenic Risk Scores (PRS)** and working with **PLINK** software.
- Developing programming and visualization skills in **R**.
- Learning about DNA data, particularly how SNP arrays are structured and utilized.

## Challenges
The project presented several exciting challenges, including:
- **Data Quality Control**: Working with real data involved navigating complex quality control steps. Understanding the data was crucial for implementing accurate methods.
- **Data Understanding**: Gaining a thorough understanding of the data was essential for effectively applying different PRS methods and ensuring reliable results.



This is my Summer Internship project at Institute for Genome Sciences(IGS) with Evelina Mocci, Research Associate at IGS.
