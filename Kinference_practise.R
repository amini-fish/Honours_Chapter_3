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

## Okay, I think we're good to go...

library(Rcpp)
library(remotes)
library(atease)
library(mvbutils)
library(vecless)
library(kinference)

## Packages loaded now for analysis/eda

data(dropbears)

kinPalette()

glimpse(dropbears)
head(dropbears)

head(dropbears@info)
head(dropbears@locinfo)

## a boring but necessary data-prep step:
dropbears <- kin_power(dropbears, k = 0.5)

dropbears

plod_HU <- find_HSPs(dropbears, keep_thresh = -10)
PLOD_loghisto(plod_HU)
  


## Need to work out how to convert gl to snpgeno
install.packages("radiator")

library(dartRverse)

remotes::install_github("thierrygosselin/radiator")

library(radiator)

# the general workflow: gl -> snpgds -> snpgeno

#### Read pilot SNP data into R

gl <- get(load("C:/Users/samue/Desktop/Honours/Chapter_2_Relatedness_paper/analysis/daly_geno_clean.Rdata")); gl

dartR.base::gl2vcf(gl, 
                   plink.bin.path = "C:/Users/samue/Desktop/Other_Research/Adspersa_writeup/plink", 
                   outpath = "C:/Users/samue/Desktop/Honours", 
                   outfile = "Pristis_test",
                   v = 5)

pristis <- radiator::genomic_converter(data = "C:/Users/samue/Desktop/Honours/Pristis_test.vcf", 
                                       output = "gds", 
                                       parallel.core = 1,
                                       filename = "pristis_test")

library(SNPRelate)

SNPRelate::snpgdsVCF2GDS(vcf.fn = "C:/Users/samue/Desktop/Honours/Pristis_test.vcf", 
                         out.fn = "C:/Users/samue/Desktop/Honours/pristis_test.gds", 
                         method = "biallelic.only",
                         verbose = T)

## Read the snpgds file in using SNPRelate 

snpgdsOpen(filename = "C:/Users/samue/Desktop/Honours/pristis_test.gds", allow.duplicate = TRUE)

library(gbasics) ## object specification for `snpgeno` objects - shouldn't need to explicitly call it
library(mvbutils)


options(repos = unique( c(
  mvb = 'https://markbravington.r-universe.dev',
  getOption( 'repos')[ 'CRAN']
)))
install.packages( "gbasics")

library(gbasics)

if( require(SNPRelate)){
  df <- gbasics::read_snpgds2snpgeno(filename = "C:/Users/samue/Desktop/Honours/pristis_test.gds", 
                                       locusID = "snp.id", 
                                       sampleID = "sample.id", 
                                       locinfoFields = c("snp.rs.id", "snp.position", "snp.chromosome", "snp.allele"))
}

## It finally worked hooray ## 

str(df)

## Estimate the Allele Frequencies 

head(df@locinfo)
head(df$locinfo)

df$locinfo$pbonzer <- re_est_ALF(df)$locinfo$pambig

View(df$locinfo$pbonzer)


# Check for all inds for matching genotype where 200 bp is the max allowed sharing? 
dups <- find_duplicates(df, max_diff_loci = 200) # none that are alarming yay

Lfish <- ilglk_geno(df)

Lfish <- ilglk_geno(df, list(nclass = 3, xlim = c(-1, 0)))
Lthreshs <- c(-1, 0)
abline(v = Lthreshs, col = "red")

## 
df_b <- kinference::kin_power(df, k = 0.5)

head(df_b$locinfo)

hetz_rich <- hetzminoo_fancy(df_b, "rich")

plod_HU <- find_HSPs(df_b, keep_thresh = -10)
PLOD_loghisto(plod_HU)

HSP_histo(plod_HU, fullsib_cut = 120, lb = -10, bin = 10)

## Implement the workflow in the vignette 

strongKin <- find_POPs(df_b, limit_pairs = choose(nrow(df_b), 2), keep_thresh = 1000)
hist(strongKin$wpsex, breaks = 50)

hist(strongKin$wpsex[strongKin$wpsex < 0.06], breaks = 50)

UPcut <- 0.055

# note that this step will be handy in the final analysis but will only be of interest here
splitPOPs <- split_FSPs_from_POPs(df_b, strongKin[strongKin$wpsex < UPcut,], gerr = 0.01)

kinPalette()

hist(splitPOPs$PLOD_FP, breaks = 50, main = "", xlab = "PLOD_FP", xlim = c(-200, 200))
abline(v = 164, col = "orange", lwd = 2)
abline(v = -116, col = "navyblue", lwd = 2)
## colours 'FSPcol' and 'POPcol' are defined by kinPalette(), for fullsibs
## and parent-offspring pairs, respectively
legend("topright", legend = c("POP expectation", "FSP expectation"), lty = c(1,1), lwd = c(2,2),
       col = c(POPcol, FSPcol), bg = "white")

## Interesting result worth discussing with PF - why are the PLODs centred around zero with only a few samples conforming to the true FSP thresholds...these must be out true FSPs with everything else being random? 

POPs <- splitPOPs[splitPOPs$PLOD_FP < 0,] ## our first set of definite kin pairs
POPs$ij <- paste(POPs$i, POPs$j) ## a pair ID


HSPs_1 <- find_HSPs(df_b, keep_thresh = -10)
PLOD_loghisto(HSPs_1)
HSP_histo(HSPs_1, fullsib_cut = 120, lb = -10, bin = 10)


POPsOrFSPs <- HSPs_1[HSPs_1$PLOD > 400,]
POPsOrFSPs$ij <- paste(POPsOrFSPs$i, POPsOrFSPs$j) ## a pair ID to match the one in POPs
FSPs <- POPsOrFSPs[!POPsOrFSPs$ij %in% POPs$ij,] ## find our full-sibling pairs

HSPs_1$ij <- paste(HSPs_1$i, HSPs_1$j)
HSPs_2 <- HSPs_1[!HSPs_1$ij %in% c(POPs$ij, FSPs$ij),]

View(HSPs_2)

## plot the empirical PLOD scores with the expected distribution
## of the PLOD for half-sibling pairs


#### Find root function ####
"find.root" <-
  function(
    f, 
    start, step, 
    fdirection = ifelse(f(start + step, ...) > f0, 
                        "increasing", "decreasing"), 
    target = 0, 
    min.x =  -Inf, max.x = Inf, 
    args.to.uniroot = list(), 
    ...
  ){
    ## Copied directly into 'kinference' from package 'handy2' to avoid faffing
    ## No longer so necessary thx2 changes in stats::uniroot, but don't wanna risk breaking anything in 'kinference' with last-minute changes
    ## Prolly should move to 'mvbutils' or something
    
    # 2024: added some sanity checks and to avoid array warnings  
    target <- as.vector( target)
    f0 <- as.vector( f(start, ...))
    stopifnot( 
      length( target) == 1,
      length( f0) == 1
    )
    
    if( f0==target)
      return( start)
    
    step <- abs(step) * ifelse(xor(f0 < target, fdirection == "decreasing"), 1, -1)
    bound <- ifelse(xor.thing <- (step < 0), min.x, max.x)
    repeat {
      new <- start + step
      if(xor(new > bound, xor.thing))
        new <- (start + bound)/2
      f1 <- as.vector( f(new, ...))
      if(xor(f0 < target, f1 < target))
        break
      start <- new
      f0 <- f1
      step <- step * 2
    }
    o <- order(x <- c(start, new))
    fvals <- c(f0, f1) - target
    dotnames <- paste(c(names(f)[1], names(list(...))), collapse = ",")
    ff <- function(...) as.vector( f(...)) - target 
    
    #  ff <- c(f[ - length(f)], list(FFFF = 1, target = 1), parse(text = "FFFF(" %&% dotnames %&% ") - target"))
    #  mode(ff) <- "function"
    
    do.call("uniroot", c(list(ff, x[o]), 
                         # f.lower = fvals[o[1]], f.upper = fvals[o[2]]), 
                         args.to.uniroot, list(...)))$root
  }
      
#### Autopick threshold ####

thresh <- autopick_threshold(df_b, kin = HSPs_2, fitrange_PLOD = c(110, 250),
                             FPtol_pairs = 2, use4th = TRUE, selecto = "ML", plot_bins = 5)


HSPs <- HSPs_2[HSPs_2$PLOD > thresh,]
fpos_rate <- thresh@FPtol_pairs / nrow(HSPs) ## approx 2% of our 'HSPs' are false-positives
fneg_rate <- thresh@info["Pr_FNeg"] ## we're missing about 0.5% of HSPs to false-negatives
nrow(HSPs) ## 7 inferred half-sibling pairs - seems larger than other methods

#### Data cleaning tutorial ####

pvals <- check6and4(df, thresh_pchisq_6and4 = c(1e-4, 1e-6), show6 = FALSE)


df_1 <- df_b[, pvals$pval4 > 1e-10]
pvals_1 <- check6and4(df_1, thresh_pchisq_6and4 = c(1e-4, 1e-6), show6 = FALSE)

hetzrich <- hetzminoo_fancy(df_1, "rich") ## a couple of outlying low samples here
df_2 <- df_1[hetzrich < 0.28,]

hetzpoor <- hetzminoo_fancy(df_2,
                            "poor",
                            hist_pars = list(breaks = 20, col = "steelblue", main = "Heterozygosity Histogram", xlab = "Observed Heterozygosity"),
                            showPlot = TRUE) ## also a couple here

pvals_2 <- check6and4(df_2, thresh_pchisq_6and4 = c(1e-4, 1e-6), show6 = FALSE)

lociPvals <- data.frame(pvals = c(pvals_2$pval4, pvals_1$pval4),
                        stage = c(rep("after", length(pvals_2$pval4)),
                                  rep("before", length(pvals_1$pval4)))
)

lm1 <- lm(pvals ~ stage, data = lociPvals)
summary(lm1)




