#' @export
#' @import survcomp
#' @import survival
#' @title Plot Kaplan-Meier curves for a gene set
#' @description Plots Kaplan-Meier survival curves for a given gene set using the
#'     cluster labels generated during the computation of \code{p_pure} to
#'     stratify patients into two survival groups. The function is a wrapper for
#'     \code{\link[survcomp]{km.coxph.plot}} in the \pkg{survcomp} package.
#' @param geneset A geneset as returned by \code{\link{saps}}.
#' @param survivalTimes A vector of survival times, as used in a call to
#'     \code{\link{saps}}.
#' @param followup A vector of 0 or 1 values, indicating whether the patient was
#' lost to followup (0) or not (1), as used in a call to \code{\link{saps}}.
#' @param title The plot title. Defaults to "Kaplan-Meier curves for geneset
#'     [\code{geneset["name"]}]".
#' @param y.label The y-axis label. Defaults to "Probability of survival".
#' @param x.label The x-axis label. Defaults to "Overall survival".
#' @param p.text Text to display in the lower left hand corner. Defaults to
#'     displaying \code{p_pure} and \code{p_pure_adj}.
#' @param ... Additional arguments to be passed to \code{\link[survcomp]{km.coxph.plot}}
#' @examples
#' # 25 patients, none lost to followup
#' followup <- rep(1, 25)
#'
#' # first 5 patients have good survival (in days)
#' time <- c(25, 27, 24, 21, 26, sample(1:3, 20, TRUE))*365
#'
#' # create data for 100 genes, 25 patients
#' dat <- matrix(rnorm(25*100), nrow=25, ncol=100)
#' colnames(dat) <- as.character(1:100)
#'
#' # create two random genesets of 5 genes each
#' set1 <- sample(colnames(dat), 5)
#' set2 <- sample(colnames(dat), 5)
#'
#' genesets <- rbind(set1, set2)
#'
#' # run saps
#' \dontrun{
#' results <- saps(genesets, dat, time, followup, random.samples=100)
#'
#' set <- results$genesets[["set1"]]
#'
#' # KM plots should not seperate at this point
#' plotKM(set, time/365, followup, x.label="Overall survival (years)")
#'
#' # increase expression levels for set1 for first 5 patients
#' dat[1:5, set1] <- dat[1:5, set1]+10
#'
#' # run saps again
#' results <- saps(genesets, dat, time, followup, random.samples=100)
#'
#' set <- results$genesets[["set1"]]
#'
#' # KM plots should now seperate
#' plotKM(set, time/365, followup, x.label="Overall survival (years)")
#' }
#' @seealso \code{\link{saps}} \code{\link{calculatePPure}}
#'     \code{\link[survcomp]{km.coxph.plot}}
plotKM <- function(geneset, survivalTimes, followup, title=NA, y.label=NA,
                   x.label=NA, p.text=NA, ...) {

  name <- geneset["name"]

  if (is.na(title))
    title <- paste("Kaplan-Meier curves for geneset ", name)

  if (is.na(y.label))
    y.label <- "Probability of survival"

  if (is.na(x.label))
    x.label <- "Overall survival"

  cluster <- geneset$cluster
  p_pure <- geneset$saps_unadjusted["p_pure"]
  p_pure_adj <- geneset$saps_adjusted["p_pure"]

  dd <- data.frame("time"=survivalTimes, "event"=followup, "cluster"=cluster)

  text <- paste("p_pure = ", round(p_pure, digits=6), ", p_pure_adj = ",
                round(p_pure_adj, digits=6))

  survcomp::km.coxph.plot(
    formula.s = survival::Surv(survivalTimes, followup) ~ cluster,
    data.s = dd, main.title=title, y.label=y.label, x.label=x.label,
    o.text=text, ...)

}


#' @export
#' @title Draw density plot of \code{p_pure} values for random gene sets
#' @description This function retrieves the \code{p_pure} values for the
#'     random gene sets generated during the computation of \code{p_random} for
#'     a given gene set. These are drawn as a density plot, with the value of
#'     \code{p_pure} for the gene set indicated. The value of \code{p_random}
#'     for the gene set is displayed as well.
#' @param geneset A geneset as returned by \code{\link{saps}}.
#' @param ... Additional arguments to be passed to \code{\link{plot}}
#' @examples
#' # 25 patients, none lost to followup
#' followup <- rep(1, 25)
#'
#' # first 5 patients have good survival (in days)
#' time <- c(25, 27, 24, 21, 26, sample(1:3, 20, TRUE))*365
#'
#' # create data for 100 genes, 25 patients
#' dat <- matrix(rnorm(25*100), nrow=25, ncol=100)
#' colnames(dat) <- as.character(1:100)
#'
#' # create two random genesets of 5 genes each
#' set1 <- sample(colnames(dat), 5)
#' set2 <- sample(colnames(dat), 5)
#'
#' genesets <- rbind(set1, set2)
#'
#' # run saps
#' \dontrun{
#' results <- saps(genesets, dat, time, followup, random.samples=100)
#'
#' set <- results$genesets[["set1"]]
#'
#' # should not be significant
#' plotRandomDensity(set)
#'
#' # increase expression levels for set1 for first 5 patients
#' dat[1:5, set1] <- dat[1:5, set1]+10
#'
#' # run saps again
#' results <- saps(genesets, dat, time, followup, random.samples=100)
#'
#' set <- results$genesets[["set1"]]
#'
#' # now it should be significant
#' plotRandomDensity(set)
#' }
#' @seealso \code{\link{saps}} \code{\link{calculatePRandom}}
plotRandomDensity <- function(geneset,  ...) {

  name <- geneset["name"]
  p_pure <- geneset$saps_unadjusted["p_pure"]
  p_random <- geneset$saps_unadjusted["p_random"]
  p_random_adj <- round(geneset$saps_adjusted["p_random"], 3)
  random_p_pures <- geneset[["random_p_pures"]]

  d <- density(-log10(random_p_pures))

  title <- paste("Significance of p_pure for ", name, "vs. random gene sets")

  plot(d, main=title, xlab="-log10 p_pure of random gene sets", cex=0.85)

  polygon(d, col=rgb(1, 0, 0, 0.5))

  arrows(x0=-log10(p_pure),x1=-log10(p_pure),y0=.25,y1=0.01,lwd=2)

  text(paste("-log10 p_pure = ", round(-log10(p_pure), digits=3)),
       x=-log10(p_pure), y=0.35, cex=0.8)

  legend <- paste("-log 10 p_pure > ", sum(random_p_pures > p_pure), " of ",
                  length(random_p_pures), " random gene sets   \n (p_random = ",
                  round(p_random, digits=5), ", p_random_adj = ",
                  round(p_random_adj, digits=5), ")   ", sep="")

  mtext(legend, side=3, line=-2.5, adj=1)

}


#' @export
#' @title Draw density plot of \code{saps_score} values for random gene sets
#' @description This function retrieves the \code{saps_score} values for the
#'     random gene sets generated during the computation of \code{saps_qvalue} for
#'     a given gene set. These are drawn as a density plot, with the value of
#'     \code{saps_score} for the gene set indicated.
#' @param geneset A geneset as returned by \code{\link{saps}}.
#' @param ... Additional arguments to be passed to \code{\link{plot}}
#' @examples
#' # 25 patients, none lost to followup
#' followup <- rep(1, 25)
#'
#' # first 5 patients have good survival (in days)
#' time <- c(25, 27, 24, 21, 26, sample(1:3, 20, TRUE))*365
#'
#' # create data for 100 genes, 25 patients
#' dat <- matrix(rnorm(25*100), nrow=25, ncol=100)
#' colnames(dat) <- as.character(1:100)
#'
#' # create two random genesets of 5 genes each
#' set1 <- sample(colnames(dat), 5)
#' set2 <- sample(colnames(dat), 5)
#'
#' genesets <- rbind(set1, set2)
#'
#' # increase expression levels for set1 for first 5 patients
#' dat[1:5, set1] <- dat[1:5, set1]+10
#'
#' # run saps and compute q-values
#' \dontrun{
#' results <- saps(genesets, dat, time, followup, random.samples=100,
#'                compute_qvalue=TRUE, qvalue.samples=10)
#'
#' set <- results$genesets[["set1"]]
#'
#' # qvalue.samples=10 is too small to achieve significance
#' plotSapsScoreDensity(set)
#' }
#' @seealso \code{\link{saps}} \code{\link{calculateQValue}}
plotSapsScoreDensity <- function(geneset,  ...) {

  name <- geneset["name"]
  saps_score <- round(geneset$saps_unadjusted["saps_score"], 3)
  saps_score_abs <- abs(saps_score)
  saps_score_adj <- round(geneset$saps_adjusted["saps_score"], 3)
  random_saps_scores <- abs(geneset[["random_saps_scores"]])

  d <- density(-log10(random_saps_scores))

  title <- paste("Significance of saps score for ", name, "vs. random gene sets")

  plot(d, main=title, xlab="-log10 saps score of random gene sets", cex=0.85)

  polygon(d, col=rgb(1, 0, 0, 0.5))

  arrows(x0=-log10(saps_score_abs),x1=-log10(saps_score_abs),y0=.25,y1=0.01,lwd=2)

  text(paste("-log10 \n          saps score = ", round(-log10(saps_score_abs),
                        digits=3)), x=-log10(saps_score_abs), y=0.3, cex=0.8)

  legend <- paste("-log 10 saps_score > ", sum(random_saps_scores > saps_score_abs),
                  " of ", length(random_saps_scores),
                  " random gene sets   \n (saps_score = ",
                  saps_score, ", saps_score_adj = ", saps_score_adj, ")   ",
                  sep="")

  mtext(legend, side=3, line=-2.5, adj=1)

}


#' @export
#' @title Plot concordance indices for a geneset
#' @description This function draws concordance indices for a given geneset
#'     relative to the concordance indices for all the genes in the dataset
#'     (i.e., the degree of enrichment for the geneset).
#' @param geneset A geneset as returned by \code{\link{saps}}.
#' @param rankedGenes A vector of concordance index z-scores. Usually this
#'     will be the \code{rankedGenes} element returned by \code{\link{saps}}.
#' @examples
#' # 25 patients, none lost to followup
#' followup <- rep(1, 25)
#'
#' # first 5 patients have good survival (in days)
#' time <- c(25, 27, 24, 21, 26, sample(1:3, 20, TRUE))*365
#'
#' # create data for 100 genes, 25 patients
#' dat <- matrix(rnorm(25*100), nrow=25, ncol=100)
#' colnames(dat) <- as.character(1:100)
#'
#' # create two random genesets of 5 genes each
#' set1 <- sample(colnames(dat), 5)
#' set2 <- sample(colnames(dat), 5)
#'
#' genesets <- rbind(set1, set2)
#'
#' # run saps
#' \dontrun{
#' results <- saps(genesets, dat, time, followup, random.samples=100)
#'
#' set <- results$genesets[["set1"]]
#'
#' # p_enrich should not be significant
#' plotEnrichment(set, results$rankedGenes)
#'
#' # increase expression levels for set1 for first 5 patients
#' dat[1:5, set1] <- dat[1:5, set1]+10
#'
#' # run saps again
#' results <- saps(genesets, dat, time, followup, random.samples=100)
#'
#' set <- results$genesets[["set1"]]
#'
#' # now it should be significant
#' plotEnrichment(set, results$rankedGenes)
#' }
#' @seealso \code{\link{saps}} \code{\link{rankConcordance}}
plotEnrichment <- function(geneset, rankedGenes) {

  name <- geneset["name"]
  geneset_genes <- geneset[["genes"]]

  gene_vals <- rankedGenes[geneset_genes]

  dummy_y1 <- rep(1, length(rankedGenes))
  dummy_y2 <- rep(1.1, length(gene_vals))

  title <- paste("Relative ranking of geneset ", name, "vs. all genes")

  plot(rankedGenes, dummy_y1, pch=1, col=rgb(0, 0, 0, 0.2), main=title,
       xlab="Concordance Index", ylab="", yaxt="n")

  points(gene_vals, dummy_y2, pch=8, col="red")

  legend("topright", c("geneset genes", "all genes"), pch=c(8, 1),
         inset=.1 ,col=c("red","black"), bty="n")

  p_enrich <- round(geneset$saps_unadjusted["p_enrich"], 3)
  p_enrich_adj <- round(geneset$saps_adjusted["p_enrich"], 3)

  text <- paste("   p_enrich = ", p_enrich, ", p_enrich_adj = ",
                p_enrich_adj, sep="")

  mtext(text, side=1, line=-2, adj=0)


}
