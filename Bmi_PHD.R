# Perform two-sample MR for BMI and PHD in European populations
# Perform analysis: Overall, males, females

setwd("C:/Users/valirien/OneDrive/1 DPhil Population Health/MR course Cambridge/BMI and PHD")


# Import libraries 
require(tidyverse)
library(data.table)
library(TwoSampleMR)
require(LDlinkR)
library(MendelianRandomization)


# Load outcome data -- 
phd_gwas <- fread('415_PheCode.v1.0.fastGWA.tsv')
# write.csv(phd_gwas, 'phd_gwas.csv')

# Import PHD data -- 
phd_outcome <- TwoSampleMR::read_outcome_data(
  filename='phd_gwas.csv',
  sep=",",
  snp_col="variant_id",
  beta_col="beta",
  se_col="standard_error",
  effect_allele_col="effect_allele",
  other_allele_col='other_allele',
  eaf_col="effect_allele_frequency",
  pval_col="p_value")

# Import BMI data --- 
bmi_gwas <- read_table("BMI_Giant_Combined.txt") 

# Filter significant variants ---
bmi_unclumped <- mutate(
  bmi_gwas,
  SNP=map_chr(str_split(SNP, pattern=":"), ~.x[1])) %>% 
  filter(P<5e-8)

dim(bmi_unclumped) # Dimension before joining: 42620    11

# Find SNPs present in the outcome data 
# The harmonised function will help you with this
bmi_unclumped_matched <- inner_join(
    bmi_unclumped,
    select(phd_gwas, chromosome, base_pair_location, variant_id), 
    by=c("CHR"="chromosome", "POS"="base_pair_location"))

# Save file (these are unclumped) - this saves to the working directory 
write.csv(bmi_unclumped_matched, "unclumped_matched.csv")

#Import unclumped exposure BMI data into TwoSampleMR package 
bmi_unclumped_ivs <- TwoSampleMR::read_exposure_data(
      filename="unclumped_matched.csv",
      sep=",",
      snp_col="SNP",
      beta_col="BETA",
      se_col="SE",
      samplesize_col="N",
      effect_allele_col="Tested_Allele",
      other_allele_col='Other_Allele',
      eaf_col="Freq_Tested_Allele",
      pval_col="P")

# Clump at r2 = 0.001
# Do not use loops, it is likely to fail
# To retrieve the data 
bmi_iv <- clump_data(bmi_unclumped_ivs, clump_r2=0.001, pop="EUR")

# Save clumped data -- 
write.csv(bmi_iv, "clumped_data.csv")


# Check instruments associated with confounders --- 
## Confounders considered 
## Find BMI SNPs that are also present in outcome data  
conf_data <- extract_outcome_data(
  outcomes=c(
     "ukb-b-18103",     # Pulse rate 
     # "ukb-b-20175",     # SBP
     # "ebi-a-GCST90038633", # Diabetes
     "ieu-a-1283",       # Alcohol intake 
     "ebi-a-GCST90093322" # Accelerometer derived PA
     ), 
  snps=bmi_iv$SNP
  )

## Harmonise data
conf_harm <- harmonise_data(
  exposure_dat=bmi_iv, 
  outcome_dat=conf_data)

n_snps <- conf_harm %>% 
  group_by(outcome) %>% 
  group_split(.) %>% 
  map_dbl(~nrow(.x)) # Get the number of snp and use that to perform bonferroni correction

# Get SNPs that are associated with confounders at GWAS significance 
res <- conf_harm %>% 
  group_by(outcome) %>% 
  group_split(.) %>% 
  map(~mr_singlesnp(.x) %>% filter(p<5e-8))

## Get SNPs associated with confounders 
## By default, these SNPs are also associated with BMI
## There considered pleiotropic variants 
snp_pleio <- map(
  .x=res, 
  ~filter(.x, !grepl(SNP, pattern="All")) %>% pull(SNP)
  ) %>% 
  unlist(.) %>% 
  unique(.)

## Filter pleitropic SNP 
bmi_iv <- filter(bmi_iv, !SNP %in% snp_pleio)

### Save dataset without pleitropic SNPs
write.csv(bmi_iv, "Instrument_pleio_snp_removed.csv")

# Harmonise data BMI and PHD 
bmi_phd_harmonised <- harmonise_data(bmi_iv, phd_outcome)

# F statistics
bmi_phd_harmonised <- mutate(
  bmi_phd_harmonised, 
  fstatistic=beta.exposure^2/se.exposure^2)

median(bmi_phd_harmonised$fstatistic)

# MR analysis 
(MR_bmi_phd <- mr(
  bmi_phd_harmonised, 
  method_list = c("mr_ivw",
                  "mr_weighted_median", 
                  "mr_egger_regression"))
)

mr_scatter_plot(MR_bmi_phd, bmi_phd_harmonised)


MRObject <- with(
    bmi_phd_harmonised, 
    MendelianRandomization::mr_input(bx=beta.exposure, 
                                     bxse=se.exposure, 
                                     by=beta.outcome, 
                                     byse=se.outcome,
                                     snps=SNP))

## Fixed and random effect IVW 
map(.x=c("fixed", "random"), ~mr_ivw(MRObject, model= .x))


Per_snp <- TwoSampleMR::mr_singlesnp(bmi_phd_harmonised)
mr_forest_plot(Per_snp)
mr_plot(MRObject)

# single snp plot
mr_forest(MRObject, ordered=TRUE)

# leave one out plot
mr_loo(MRObject)

# funnel plot
mr_funnel(MRObject)

## Heterogeneity 
(MR_heterogeneity <- mr_heterogeneity(bmi_phd_harmonised))

