## ----, message=FALSE-----------------------------------------------------
library(saps)

## ------------------------------------------------------------------------
followup <- rep(1, 25)

## ------------------------------------------------------------------------
time <- c(25, 27, 24, 21, 26, sample(1:3, 20, TRUE))*365

## ------------------------------------------------------------------------
dat <- matrix(rnorm(25*100), nrow=25, ncol=100)
colnames(dat) <- as.character(1:100)

## ------------------------------------------------------------------------
set1 <- sample(colnames(dat), 5)
set2 <- sample(colnames(dat), 5)

## ------------------------------------------------------------------------
genesets <- rbind(set1, set2)

## ----firstrun, echo=FALSE, message=FALSE, cache=TRUE---------------------
results <- saps(genesets, dat, time, followup, random.samples=100)

## ------------------------------------------------------------------------
saps_table <- results$saps_table
saps_table[1:7]

## ------------------------------------------------------------------------
dat2 <- dat
dat2[1:5, set1] <- dat2[1:5, set1]+10

## ----secondrun, echo=FALSE, message=FALSE, cache=TRUE--------------------
results <- saps(genesets, dat2, time, followup, random.samples=100)
saps_table <- results$saps_table
saps_table[1:7]

## ------------------------------------------------------------------------
saps_table[8:11]

## ----, message=FALSE-----------------------------------------------------
set1 <- results$genesets[["set1"]]

plotKM(set1, time/365, followup, x.label="Overall survival (years)")

## ------------------------------------------------------------------------
plotRandomDensity(set1)

## ------------------------------------------------------------------------
plotEnrichment(set1, results$rankedGenes)

## ----qvalue, message=FALSE, cache=TRUE-----------------------------------
results <- saps(genesets, dat2, time, followup, random.samples=100, 
                compute_qvalue=TRUE, qvalue.samples=10)

saps_table <- results$saps_table

saps_table[, c("saps_qvalue", "saps_qvalue_adj")]

## ------------------------------------------------------------------------
set1 <- results$genesets[["set1"]]

plotSapsScoreDensity(set1)

