
##------------------------------
## Install Packages 
##------------------------------

## CRAN dependencies:

install.packages("Rcpp")
install.packages("remotes")
install.packages("tensorA")

## support packages from Mark Bravington's repo:

install.packages(pkgs = c("atease", "mvbutils", "vecless"), repos = "https://markbravington.github.io/Rmvb-repo")

install.packages(pkgs = c("gbasics"), repos = "https://markbravington.github.io/Rmvb-repo")

## kinference and gbasics: note that you'll need to specify your own filepaths here. The files are attached.

remotes::install_local("C:/Users/samue/Desktop/Honours/Chapter_2_Relatedness_paper/Kinference/kinference-master/kinference")
remotes::install_local("C:/Users/samue/Desktop/Honours/Chapter_2_Relatedness_paper/Kinference/gbasics-master/gbasics")
remotes::install_github("thierrygosselin/radiator")
install.packages("radiator")

## Okay, I think we're good to go...

##--------------------------------
## Load said packages 
##--------------------------------

library(Rcpp)
library(remotes)
library(atease)
library(mvbutils)
library(vecless)
library(kinference)
library(dartRverse)
library(radiator)
library(SNPRelate)
library(mvbutils)
library(gbasics) ## object specification for `snpgeno` objects - shouldn't need to explicitly call it
options(repos = unique( c(
  mvb = 'https://markbravington.r-universe.dev',
  getOption( 'repos')[ 'CRAN']
)))
install.packages( "gbasics")
library(gbasics)

# the general workflow: gl -> snpgds -> snpgeno

##----------------------------------
## Load in data, run conversions, tune final form
##----------------------------------

gl <- get(load("C:/Users/samue/Desktop/daly_geno_raw.Rdata")); gl

## OR 

gl <- get(load("C:/Users/samue/Desktop/daly_geno_clean.Rdata")); gl

dartR.base::gl2vcf(gl, 
                   plink.bin.path = "C:/Users/samue/Desktop/Other_Research/Adspersa_writeup/plink", 
                   outpath = "C:/Users/samue/Desktop/Honours", 
                   outfile = "Pristis_test",
                   v = 5)

pristis <- radiator::genomic_converter(data = "C:/Users/samue/Desktop/Honours/Pristis_test.vcf", 
                                       output = "gds", 
                                       parallel.core = 1,
                                       filename = "pristis_test")



SNPRelate::snpgdsVCF2GDS(vcf.fn = "C:/Users/samue/Desktop/Honours/Pristis_test.vcf", 
                         out.fn = "C:/Users/samue/Desktop/Honours/pristis_test.gds", 
                         method = "biallelic.only",
                         verbose = T)

## Read the snpgds file in using SNPRelate 

snpgdsOpen(filename = "C:/Users/samue/Desktop/Honours/pristis_test.gds", allow.duplicate = TRUE)

if( require(SNPRelate)){
  df <- gbasics::read_snpgds2snpgeno(filename = "C:/Users/samue/Desktop/Honours/pristis_test.gds", 
                                       locusID = "snp.id", 
                                       sampleID = "sample.id", 
                                       locinfoFields = c("snp.rs.id", "snp.position", "snp.chromosome", "snp.allele"))
}

#### It finally worked hooray ####

str(df)

## Estimate the Allele Frequencies 

head(df$locinfo)

df$locinfo$pbonzer <- re_est_ALF(df)$locinfo$pambig

##------------------------------------
## Data cleaning + QAQC
##------------------------------------

## 1. Check 4 and 6 (i.e., hwe check for the 4 allele model and the 6 allele model)

pvals <- check6and4(df, thresh_pchisq_6and4 = c(1e-4, 1e-6), show6 = FALSE)

## Remove the worst loci p < 1x10^-10

df_1 <- df[, pvals$pval4 > 1e-6]

pvals_1 <- check6and4(df_1, thresh_pchisq_6and4 = c(1e-4, 1e-6), show6 = FALSE)

## Check for duplicated samples 

dups <- find_duplicates(df_1, max_diff_loci = 50) # none that are alarming yay!

## Check for fish from different pops/species

ilglk <- ilglk_geno(df_1, list(nclass = 1000, xlim = c(-2, 2)))

abline(v = -2, col = "red") ## do I even need a cutoff?!

## Estimate power to find kin i.e., how well will this set of markers work 

df_2 <- kinference::kin_power(df_1, k = 0.5)

head(df_2$locinfo) ## Check that the info we want is there...looks good. 

## Examine the heterozygosities of our inds: 
## Importnant for nulls + contamination (low + high)

hetz_rich <- hetzminoo_fancy(df_2, "rich", list(xlim = c(-1,0.4)))

hetz_poor <- hetzminoo_fancy(df_2, "poor")

## Quick clean up 

df_3 <- df_2[(hetz_rich > 0), ]

## Identify degraded samples by their lack of heterozygotes and excess of
## double-nulls.

hetz_poor <- hetzminoo_fancy(df_3, "poor")

# sanity check # sanity checkdf_2 

## Revisit our allele frequencies

pvals_2 <- check6and4(df_3, thresh_pchisq_6and4 = c(1e-4, 1e-6), show6 = FALSE)

## It appears we have only removed our good loci? 

lociPvals <- data.frame(pvals = c(pvals_2$pval4, pvals_1$pval4),
                        stage = c(rep("after", length(pvals_2$pval4)),
                                  rep("before", length(pvals_1$pval4)))
)

lm1 <- lm(pvals ~ stage, data = lociPvals)
summary(lm1) 

## The p-values improved by aout 3.3%

###-------------------------------
### Kin Finding 
###-------------------------------

## Kin finding workflow ##

# Whilst unlikely, check data for POPS first

pops <- find_POPs(df_3, limit_pairs = choose(nrow(df_3), 2), keep_thresh = 1, nbins = 100)
hist(pops$wpsex, breaks = 30, xlim = c(0.018, 0.15), ylim = c(0, 1000), xlab = "wpsex")
## There's a bit of a dip around 0.04, but it's not _clear_ separation
abline(v = 0.04, lty = 2, col = "red")
box()

with(pops, plot.default(nABOO ~ wpsex))
abline(v = 0.04, col = "red", lty = 2)


splitPOPs <- split_FSPs_from_POPs(df_3, pops[pops$wpsex < 0.06, ], gerr = 0.01)

str(splitPOPs)

ggplot() + 
  geom_histogram(data = splitPOPs, aes(x = PLOD_FP), bins = 50, fill = "grey", col = "black") + 
  geom_vline(xintercept = splitPOPs@E_FSP, col = "orange") + 
  geom_vline(xintercept = splitPOPs@E_POP, col = "navyblue") + 
  theme_bw()

## find HSPs-or-stronger using find_HSPs
HSPs_1 <- find_HSPs(df_3, keep_thresh = -10, limit_pairs = choose(nrow(df_3),2))
PLOD_loghisto(HSPs_1)

HSP_histo(HSPs_1, fullsib_cut = 150, lb = -10, bin = 10)

## The hist looks weird asf with high freq of HSP? 

POPsOrFSPs <- HSPs_1[HSPs_1$PLOD > 150,]
POPsOrFSPs$ij <- paste(POPsOrFSPs$i, POPsOrFSPs$j) ## a pair ID to match the one in POPs
FSPs <- POPsOrFSPs[!POPsOrFSPs$ij %in% POPs$ij,] ## find our full-sibling pairs
FSPs

nrow(FSPs)

## This should give us all the FSPs in our data? 

HSPs_1$ij <- paste(HSPs_1$i, HSPs_1$j)
HSPs_2 <- HSPs_1[!HSPs_1$ij %in% c(POPs$ij, FSPs$ij),]

## plot the empirical PLOD scores with the expected distribution
## of the PLOD for half-sibling pairs

thresh <- autopick_threshold(df_3, kin = HSPs_2,
                             fitrange_PLOD = c(10, 150),
                             FPtol_pairs = 1, 
                             use4th = TRUE, 
                             selecto = "ML", 
                             plot_bins = 5)


HSPs <- HSPs_2[HSPs_2$PLOD > thresh,]
fpos_rate <- thresh@FPtol_pairs / nrow(HSPs); fpos_rate ## approx 2% of our 'HSPs' are false-positives
fneg_rate <- thresh@info["Pr_FNeg"]; fneg_rate ## we're missing about 0.5% of HSPs to false-negatives
nrow(HSPs) ## 52 inferred half-sibling pairs

## Add intividual metadata to the kin pairs ##

POPs # This is a false positive
FSPs
HSPs

kinf_kin <- rbind(POPs, FSPs, HSPs); kinf_kin
dim(kinf_kin)

kinf_kin$rel <- ifelse(kinf_kin$PLOD > 300, "FSP", ifelse(kinf_kin$PLOD > 0, "HSP", "POP")) # this needs to be adjusted to fit the datasets...

kinf_kin

ind_names <- df$info 
ind_names <- ind_names %>%
  mutate(Our_sample = sub("^pop1_", "", Our_sample)) %>%
  mutate(row_id = row_number())

ind_names

kinf_kin <- kinf_kin %>%
  left_join(ind_names, by = c("i" = "row_id" )) %>%
  rename(ind_i = Our_sample)%>%
  left_join(ind_names, by = c("j" = "row_id" )) %>%
  rename(ind_j = Our_sample)

kinf_kin

ind_metadata <- read.csv("C:/Users/samue/Desktop/10658_prelim_report/Pristis_metadata_SRA.csv")

kinf_kin2 <- kinf_kin %>%
  left_join(ind_metadata, by = c("ind_i" = "id")) %>%
  rename(birth_year_1 = birth_year, billabong_1 = Site_ID, capture_year1 = year_caught) %>%
  left_join(ind_metadata, by = c("ind_j" = "id")) %>%
  rename(birth_year_2 = birth_year, billabong_2 = Site_ID, capture_year2 = year_caught) %>%
  mutate(
    birth_year_diff = abs(birth_year_1 - birth_year_2),  # Absolute age difference
    Within_Billabong = as.integer(billabong_1 == billabong_2), # returns true false as 1,0
    Within_Cohort = as.integer(capture_year1 == capture_year2)) %>% #returns true false as 1,0
  dplyr::select(ind_i, ind_j, rel,billabong_1, billabong_2, birth_year_1, birth_year_2, capture_year1, capture_year2, birth_year_diff, Within_Billabong, Within_Cohort)

kin_pairs <- tibble(kinf_kin2)

print(kin_pairs)

## Miscelaneous shit


PLINK_BIN_PATH <- "C:/Users/samue/Desktop/Other_Research/Adspersa_writeup/plink" # <<<EDIT ME>>>
# <<<EDIT ME>>>
GENLIGHT_RDATA <- "C:/Users/samue/Desktop/daly_geno_raw.Rdata"                    # <<<EDIT ME>>>
VCF_FILE       <- file.path(WORK_DIR, "Pristis_test.vcf")


if (!dir.exists(WORK_DIR)) dir.create(WORK_DIR, recursive = TRUE)

# Expecting an object named `gl` inside the .Rdata
gl <- get(load(GENLIGHT_RDATA))
gl

dartR.base::gl2vcf(
  gl,
  plink.bin.path = PLINK_BIN_PATH,
  outpath        = WORK_DIR,
  outfile        = tools::file_path_sans_ext(basename(VCF_FILE)),
  v              = 5
)

SNPRelate::snpgdsVCF2GDS(
  vcf.fn = VCF_FILE,
  out.fn = GDS_FILE,
  method = "biallelic.only",
  verbose = TRUE
)

snpg <- SNPRelate::snpgdsOpen(filename = GDS_FILE, allow.duplicate = TRUE)
SNPRelate::snpgdsClose(snpg)

