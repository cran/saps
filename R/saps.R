#' Implements Significance Analysis of Prognostic Signatures (SAPS), a
#' robust method for determining prognostically significant gene sets
#'
#' \code{\link{saps}} will usually be the only function needed.
#'
#'
#' @references Beck AH, Knoblauch NW, Hefti MM, Kaplan J, Schnitt SJ, et al.
#' (2013) Significance Analysis of Prognostic Signatures. PLoS Comput Biol 9(1):
#' e1002875.doi:10.1371/journal.pcbi.1002875
#' @name saps-package
#' @docType package
NULL


#' @export
#' @title Compute SAPS statistics
#' @description This is the main user interface to the \pkg{saps} package, and is
#' usually the only function needed.
#' @details saps provides a robust method for identifying biologically significant
#' gene sets associated with patient survival. Three basic statistics are computed.
#' First, patients are clustered into two survival groups based on differential
#' expression of a candidate gene set. p_pure is calculated as the probability of
#' no survival difference between the two groups.
#'
#' Next, the same procedure is applied to randomly generated gene sets, and p_random
#' is calculated as the proportion achieving a p_pure as significant as the candidate
#' gene set. Finally, a pre-ranked Gene Set Enrichment Analysis (GSEA) is performed
#' by ranking all genes by concordance index, and p_enrich is computed to indicate
#' the degree to which the candidate gene set is enriched for genes with univariate
#' prognostic significance.
#'
#' A saps_score is calculated to summarize the three statistics, and
#' optionally a saps_qvalue is computed to estimate the significance of the
#' saps_score by calculating the saps_score for random gene sets.
#' @param candidateGeneSets A matrix with at least one row, where each row represents
#' a gene set, and the column values are gene identifiers. The row names should contain
#' unique names for the gene sets. The column values may contain \code{NA} values, since
#' in general gene sets will have differing lengths.
#' @param dataSet A matrix, where the column names are gene identifiers
#' (in the same format as the values in \code{candidateGeneSets})
#' and the values are gene expression levels. Each row should contain data for a single
#' patient.
#' @param survivalTimes A vector of survival times. The length must equal the number of
#' rows (i.e. patients) in \code{dataSet}.
#' @param followup A vector of 0 or 1 values, indicating whether the patient was
#' lost to followup (0) or not (1). The length must equal the number of rows
#' (i.e. patients) in \code{dataSet}.
#' @param random.samples An integer that specifies how many random gene sets to sample
#' when computing P_random. Defaults to 1000.
#' @param cpus An integer that specifies the number of cpus/cores to be used when
#' calculating P_enrichment. If greater than 1 (the default), the \pkg{snowfall}
#' package must be installed or an error will occur.
#' @param gsea.perm The number of permutations to be used when calculating
#' p_enrich. This is passed to the \code{\link[piano]{runGSA}} function in the
#' \pkg{piano} package. Defaults to 1000.
#' @param compute_qvalue A boolean indicating whether to include calculation
#' of the saps q_value. Setting this to \code{TRUE} will significantly
#' increase the computational time.
#' @param qvalue.samples An integer that specifies how many random gene sets to
#' sample when computing the saps q_value. Defaults to 1000.
#' @param verbose A boolean indicating whether to display status messages during
#' computation. Defaults to \code{TRUE}.
#' @return The function returns a list with the following elements:
#'
#' \item{rankedGenes}{Vector of concordance index z-scores for the genes in
#'    \code{dataSet}, named by gene identifier.}
#' \item{geneset.count}{The number of gene sets analyzed.}
#' \item{genesets}{A list of genesets (see below).}
#' \item{saps_table}{A dataframe summarizing the adjusted and unadjusted
#'    saps statistics for each geneset analyzed. The dataframe contains
#'    the following columns: \code{size, p_pure, p_random, p_enrich,
#'    direction, saps_score, saps_qvalue, p_pure_adj, p_random_adj,
#'    p_enrich_adj, saps_score_adj, saps_qvalue_adj}. Each row summarizes
#'    a single geneset. Note that the saps statistics are stored
#'    with each individual \code{geneset} as well; this table is
#'    provided simply for convenience.}
#' \code{genesets} is in turn a list with the following elements:
#'
#' \item{name}{The name of the geneset.}
#' \item{size}{The number of genes in the geneset.}
#' \item{genes}{Vector of gene labels for this geneset.}
#' \item{saps_unadjusted}{Vector with elements \code{p_pure}, \code{p_random},
#'     \code{p_enrich}, \code{saps_score}, and \code{saps_qvalue} containing
#'     the respective unadjusted p-values.}
#' \item{saps_adjusted}{Vector with elements \code{p_pure}, \code{p_random},
#'     \code{p_enrich}, \code{saps_score}, and \code{saps_qvalue} containing
#'     the respective p-values adjusted for multiple comparisons.}
#' \item{cluster}{Vector of assigned cluster (1 or 2) for each patient using this
#'     candidate geneset.}
#' \item{random_p_pures}{Vector of p_pure values for each random geneset generated
#'     during the computation of p_random.}
#' \item{random_saps_scores}{Vector of saps_score values for each random geneset
#'     generated during the computation of saps_qvalue.}
#' \item{direction}{Direction (-1 or 1) of the enrichment association for this geneset.}
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
#' # compute saps
#' \dontrun{
#' results <- saps(genesets, dat, time, followup, random.samples=100)
#'
#' # check results
#' saps_table <- results$saps_table
#' saps_table[1:7]
#'
#' # increase expression levels for set1 for first 5 patients
#' dat[1:5, set1] <- dat[1:5, set1]+10
#'
#' # run again, should get significant values for set1
#' results <- saps(genesets, dat, time, followup, random.samples=100)
#'
#' # check results
#' saps_table <- results$saps_table
#' saps_table[1:7]
#' }
#' @seealso \code{\link[survival]{survdiff}} \code{\link[survcomp]{concordance.index}}
#'     \code{\link[piano]{runGSA}}
#' @references Beck AH, Knoblauch NW, Hefti MM, Kaplan J, Schnitt SJ, et al.
#' (2013) Significance Analysis of Prognostic Signatures. PLoS Comput Biol 9(1):
#' e1002875.doi:10.1371/journal.pcbi.1002875
saps <- function(candidateGeneSets, dataSet, survivalTimes,
                 followup, random.samples=1000, cpus=1, gsea.perm=1000,
                 compute_qvalue=FALSE, qvalue.samples=1000, verbose=TRUE) {

  if ((cpus > 1) & (!is.installed("snowfall")))
    stop("'snowfall' package not found (required for multiple CPU support)")

  candidateSetCount <- nrow(candidateGeneSets)

  # prepare results
  results <- list("rankedGenes"=NA, "geneset.count"=candidateSetCount,
                  "genesets"=list(), "saps_table"=NA)

  # get concordance index

  if (verbose)
    message("Calculating concordance index...", appendLF=FALSE)

  ci <- rankConcordance(dataSet, survivalTimes, followup)

  # use z-scores of CI
  rankedGenes <- ci[, -1]

  results$rankedGenes <- rankedGenes

  if (verbose)
    message("done.")


  setNames <- rownames(candidateGeneSets)
  geneNames <- colnames(dataSet)

  # prepare the matrices that will hold the adjusted and
  # unadjusted SAPS statistics
  saps_unadjusted <- matrix(nrow=candidateSetCount, ncol=7,
                            dimnames=list(setNames,
                                          c("size", "p_pure",
                                            "p_random", "p_enrich",
                                            "direction", "saps_score",
                                            "saps_qvalue")))

  saps_adjusted <- saps_unadjusted

  # compute saps statistics for each candidate gene set
  for (i in 1:candidateSetCount) {

    candidateGeneSet <- candidateGeneSets[i,,drop=FALSE]

    setName <- setNames[i]

    # get candidate genes
    candidateGenes <- as.character(candidateGeneSet[!is.na(candidateGeneSet)])
    commonGenes <- intersect(geneNames, candidateGenes)

    candidateSetSize <- length(commonGenes)

    saps_unadjusted[setName, "size"] <- candidateSetSize
    saps_adjusted[setName, "size"] <- candidateSetSize

    # prepare results for this geneset
    saps_vec <- vector(mode="numeric", length=5)
    names(saps_vec) <- c("p_pure", "p_random", "p_enrich", "saps_score",
                           "saps_qvalue")

    set_results <- list("name" = setName, "size" = candidateSetSize,
                        "genes" = commonGenes, "cluster" = NA,
                        "random_p_pures" = NA, "random_saps_scores" = NA,
                        "direction" = NA, "saps_unadjusted" = saps_vec,
                        "saps_adjusted" = saps_vec)


    if(candidateSetSize == 0) {

      warning(c("No gene data found for gene set ", setName, ", cannot compute SAPS."))

    }
    else {

      if (verbose) {
        message(c("Using gene set ", setName, ", size = ", candidateSetSize))
        message(c("gene set #", i, " of ", candidateSetCount))
      }

      scaledData <- scale(dataSet[, commonGenes])

      if (verbose)
        message("Calculating P_pure...", appendLF=FALSE)

      pure <- calculatePPure(scaledData, survivalTimes, followup)

      p_pure <- pure[["p_pure"]]

      saps_unadjusted[setName, "p_pure"] <- p_pure
      set_results["cluster"] <- pure["cluster"]

      if (verbose)
        message("done.")

      if (verbose)
        message("Calculating P_random...", appendLF=FALSE)

      random <- calculatePRandom(dataSet, candidateSetSize, p_pure, survivalTimes,
                                   followup, random.samples)

      p_random <- random[["p_random"]]

      # adjust 0 value
      if (p_random == 0)
        p_random <- 1/(random.samples+1)

      saps_unadjusted[setName, "p_random"] <- p_random
      set_results["random_p_pures"] <- random["p_pures"]

      if (verbose)
        message("done.")

      if (verbose)
        message("Calculating P_enrichment...", appendLF=FALSE)

      gsa_results <- calculatePEnrichment(rankedGenes, candidateGeneSet,
                                          cpus, gsea.perm)

      p_enrich <- gsa_results$P_enrichment
      direction <- gsa_results$direction

      # adjust 0 values
      if (p_enrich == 0)
        p_enrich <- 1/(gsea.perm+1)

      saps_unadjusted[setName, "p_enrich"] <- p_enrich
      saps_unadjusted[setName, "direction"] <- direction
      saps_adjusted[setName, "direction"] <- direction

      if (verbose)
        message("done.")

      # compute q_value if requested
      if (compute_qvalue) {

        saps_score <- -log10(max(p_pure, p_random, p_enrich)) * direction

        if (verbose)
          message("Calculating q_value...", appendLF=FALSE)

        qval <- calculateQValue(dataSet, candidateSetSize, survivalTimes,
                                   followup, saps_score, random.samples,
                                   qvalue.samples, cpus, gsea.perm, rankedGenes)

        q_value <- qval[["q_value"]]
        saps_scores <- qval["random_saps_scores"]

        # adjust 0 values
        if (q_value == 0)
          q_value <- 1/(qvalue.samples+1)

        saps_unadjusted[setName, "saps_qvalue"] <- q_value

        set_results["random_saps_scores"] <- saps_scores

        if (verbose)
          message("done.")

      }

    }

    results$genesets[[setName]] <- set_results

  }

  # adjust p-values (if needed)
  if (candidateSetCount > 1) {

    if (verbose)
      message("Adjusting p-values...", appendLF=FALSE)

    if (compute_qvalue)
      to_adjust <- saps_unadjusted[, c("p_pure", "p_random",
                                      "p_enrich", "saps_qvalue")]
    else
      to_adjust <- saps_unadjusted[, c("p_pure", "p_random", "p_enrich")]

    adjusted <- apply(to_adjust, 2, p.adjust, "BH")

    saps_adjusted[, "p_pure"] <- adjusted[, "p_pure"]
    saps_adjusted[, "p_random"] <- adjusted[, "p_random"]
    saps_adjusted[, "p_enrich"] <- adjusted[, "p_enrich"]

    if (compute_qvalue)
      saps_adjusted[, "saps_qvalue"] <- adjusted[, "saps_qvalue"]

    if (verbose)
      message("done.")

  }

  # calculate saps scores
  if (verbose)
    message("Calculating saps scores...", appendLF=FALSE)

  p_max <- pmax(saps_unadjusted[,"p_pure"],
               saps_unadjusted[,"p_random"],
               saps_unadjusted[,"p_enrich"])

  direction <- saps_unadjusted[,"direction"]

  saps_scores <- -log10(p_max) * direction

  saps_unadjusted[, "saps_score"] <- saps_scores

  p_max <- pmax(saps_adjusted[,"p_pure"],
               saps_adjusted[,"p_random"],
               saps_adjusted[,"p_enrich"])

  saps_scores <- -log10(p_max) * direction

  saps_adjusted[, "saps_score"] <- saps_scores

  if (verbose)
    message("done.")

  # finally, store each saps statistic with the corresponding
  # geneSet in the results list (also save the statistics in the
  # saps_test convenience table)
  if (verbose)
    message("Saving SAPS statistics...", appendLF=FALSE)

  saps_table <- data.frame("size"=saps_unadjusted[, "size"],
                           "p_pure"=saps_unadjusted[,"p_pure"],
                           "p_random"=saps_unadjusted[, "p_random"],
                           "p_enrich"=saps_unadjusted[, "p_enrich"],
                           "direction"=saps_unadjusted[, "direction"],
                           "saps_score"=saps_unadjusted[, "saps_score"],
                           "saps_qvalue"=saps_unadjusted[, "saps_qvalue"],
                           "p_pure_adj"=saps_adjusted[, "p_pure"],
                           "p_random_adj"=saps_adjusted[, "p_random"],
                           "p_enrich_adj"=saps_adjusted[, "p_enrich"],
                           "saps_score_adj"=saps_adjusted[, "saps_score"],
                           "saps_qvalue_adj"=saps_adjusted[, "saps_qvalue"])

  rownames(saps_table) <- rownames(saps_unadjusted)

  results$genesets <- lapply(results$genesets, save_saps,
                            saps_unadjusted, saps_adjusted)


  results$saps_table <- saps_table

  if (verbose)
    message("done.")


  return(results)

}


#' @export
#' @import reshape2
#' @import piano
#' @title Compute P_enrichment
#' @description This function performs a pre-ranked gene set enrichment
#' analysis (GSEA) to evaluate the degree to which a candidate gene set is
#' overrepresented at the top or bottom extremes of a ranked list of
#' concordance indices. This function is normally called by
#' \code{\link{saps}}.
#' @param rankedGenes An \emph{nx1} matrix of concordance indices for \emph{n}
#' genes. Generally this will be the z-score returned by
#' \code{\link{rankConcordance}}. The row names should contain gene identifiers.
#' @param candidateGeneSet A \emph{1xp} matrix of \emph{p} gene identifiers.
#' The row name should contain a name for the gene set.
#' @param cpus This value is passed to the \code{\link[piano]{runGSA}} function
#' in the \pkg{piano} package. For multi-core CPUs, this value should be set to
#' the number of cores (which will significantly improve the computational time).
#' @param gsea.perm The number of permutations to be used in the GSEA. This
#' value is passed to \code{\link[piano]{runGSA}}.
#' @return The function returns a matrix with the following columns:
#'
#'  \item{P_enrichment}{the enrichment score}
#'  \item{direction}{either 1 or -1 depending on the direction of association}
#'
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
#' # tweak data for first 5 patients for set1
#' dat[1:5, set1] <- dat[1:5, set1]+10
#'
#' # rank all genes by concordance index
#' \dontrun{
#' ci <- rankConcordance(dat, time, followup)[,"z"]
#'
#' # set1 should achieve significance
#' p_enrich <- calculatePEnrichment(ci, genesets["set1",,drop=FALSE], cpus=1)
#' p_enrich
#'
#' # set2 should not
#' p_enrich <- calculatePEnrichment(ci, genesets["set2",,drop=FALSE], cpus=1)
#' p_enrich
#' }
#' @seealso \code{\link{saps}} \code{\link[piano]{runGSA}}
#' @references Beck AH, Knoblauch NW, Hefti MM, Kaplan J, Schnitt SJ, et al.
#' (2013) Significance Analysis of Prognostic Signatures. PLoS Comput Biol 9(1):
#' e1002875.doi:10.1371/journal.pcbi.1002875
#'
#' Subramanian A, Tamayo P, Mootha VK, Mukherjee S, Ebert BL, et al. (2005)
#' Gene set enrichment analysis: a knowledge-based approach for interpreting
#' genome-wide expression profiles. Proc Natl Acad Sci USA 102: 15545-15550.
calculatePEnrichment <- function(rankedGenes, candidateGeneSet,
                                 cpus, gsea.perm=1000) {

  # reshape candidate gene set into long form for piano input
  candidateSetLong <- (reshape2::melt(as.matrix(candidateGeneSet),
                                      na.rm=TRUE))[, c("value", "Var1")]

  gsc <- piano::loadGSC(candidateSetLong, type="data.frame")

  # run GSEA using piano package
  gsa <- piano::runGSA(rankedGenes, geneSetStat="gsea", signifMethod="geneSampling",
                adjMethod="fdr", gsc=gsc, gsSizeLim=c(1, 250), nPerm=gsea.perm,
                ncpus=cpus, verbose=FALSE)

  # get p-values for up and down-regulated gene sets
  gsa_results <- cbind(names(gsa$gsc), gsa$pDistinctDirUp, gsa$pDistinctDirDn)

  # build direction vector
  direction <- ifelse(is.na(gsa_results[, 2]), -1, 1)
  direction <- cbind(direction)
  rownames(direction) <- gsa_results[, 1]

  # merge the up and down p-values
  neg_pos <- cbind(gsa_results[, 2])

  neg_na <- is.na(neg_pos[, 1])

  neg_pos[neg_na] <- (gsa_results[, 3])[neg_na]

  gsa_results <- cbind(gsa_results[, 1], neg_pos)

  rownames(gsa_results) <- gsa_results[, 1]

  gsea_p <- cbind("P_enrichment"=as.numeric(gsa_results[, 2]))
  rownames(gsea_p) <- rownames(gsa_results)

  # merge p-values and direction vector
  results <- merge(gsea_p, direction, by="row.names")

  # restore row names
  rownames(results) <- results$Row.names
  results <- results[, -1]

  return(results)

}


#' @export
#' @import survival
#' @title Compute P_pure
#' @description This function stratifies patients into two groups via k-means
#' clustering (k=2) on an \emph{nxp} matrix consisting of \emph{n} patients
#' and \emph{p} genes in the candidate prognostic set. It is normally called
#' by \code{\link{saps}}.
#' @param geneData An \emph{nxp} matrix consisting of \emph{n} patients
#' and \emph{p} genes in the candidate prognostic geneset.
#' @param survivalTimes A vector of survival times. The length must equal
#' the number of rows \emph{n} in \code{geneData}.
#' @param followup A vector of 0 or 1 values, indicating whether the patient was
#' lost to followup (0) or not (1). The length must equal the number of rows
#' (i.e. patients) in \code{geneData}.
#' @return A list with the following elements:
#' \item{p_pure}{A log-rank p-value indicating the probability that the two groups
#'     show no survival difference.}
#' \item{cluster}{Vector of assigned cluster (1 or 2) for each patient using the
#'     supplied candidate prognostic geneset.}
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
#' # create random genesets of 5 genes
#' set1 <- sample(colnames(dat), 5)
#'
#' # get gene data for set1
#' set1_data <- dat[, set1]
#'
#' # shouldn't achieve significance
#' p_pure <- calculatePPure(set1_data, time, followup)
#' p_pure$p_pure
#'
#' # alter expression data for first 5 patients for set1
#' dat[1:5, set1] <- dat[1:5, set1]+10
#'
#' set1_data <- dat[, set1]
#'
#' # now p_pure should be significant
#' p_pure <- calculatePPure(set1_data, time, followup)
#' p_pure$p_pure
#' @seealso \code{\link{saps}}
#' @references Beck AH, Knoblauch NW, Hefti MM, Kaplan J, Schnitt SJ, et al.
#' (2013) Significance Analysis of Prognostic Signatures. PLoS Comput Biol 9(1):
#' e1002875.doi:10.1371/journal.pcbi.1002875
calculatePPure <- function(geneData, survivalTimes, followup) {

  # cluster the patients on the candidate genes
  cluster <- kmeans(geneData, 2)$cluster

  # compute probability of no survival difference
  survtest <- survival::survdiff(Surv(survivalTimes, followup) ~ cluster)

  p_pure <- 1 - pchisq(survtest$chisq, 1)

  return(list("p_pure"=p_pure, "cluster"=cluster))

}


#' @export
#' @title Compute P_random
#' @description This function randomly samples gene sets, and calculates
#' P_pure (via \code{\link{calculatePPure}}) for each one. P_random is the
#' proportion of randomly sampled gene sets achieving a P_pure at least as
#' significant as the provided \code{p_pure}. This function is normally called
#' by \code{\link{saps}}.
#' @param dataSet A matrix, where the column names are gene identifiers
#' and the values are gene expression levels. Each row should contain data for a
#' single patient.
#' @param sampleSize The desired size for the randomly sampled gene sets.
#' @param p_pure The candidate P_pure against which to compare the P_pure
#' values for the randomly generated gene sets.
#' @param survivalTimes A vector of survival times. The length must equal
#' the number of rows in \code{dataSet}.
#' @param followup A vector of 0 or 1 values, indicating whether the patient was
#' lost to followup (0) or not (1). The length must equal the number of rows
#' (i.e. patients) in \code{dataSet}.
#' @param random.samples The number of random gene sets to sample.
#' @return A list with the following elements:
#' \item{p_random}{The proportion of randomly sampled gene sets with a calculated
#'     p_pure at least as significant as the provided \code{p_pure}.}
#' \item{p_pures}{A vector of calculated p_pure values for each randomly
#'     generated geneset.}
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
#' # relatively low threshold
#' p_pure <- 0.05
#' \dontrun{
#' p_random <- calculatePRandom(dat, 5, p_pure, time, followup, random.samples=100)
#' p_random$p_random
#' hist(p_random$p_pures)
#' length(p_random$p_pures[p_random$p_pures <= p_pure])
#'
#' # set a more stringent threshold
#' p_pure <- 0.001
#'
#' p_random <- calculatePRandom(dat, 5, p_pure, time, followup, random.samples=100)
#' length(p_random$p_pures[p_random$p_pures <= p_pure])
#' }
#' @seealso \code{\link{saps}}
#' @references Beck AH, Knoblauch NW, Hefti MM, Kaplan J, Schnitt SJ, et al.
#' (2013) Significance Analysis of Prognostic Signatures. PLoS Comput Biol 9(1):
#' e1002875.doi:10.1371/journal.pcbi.1002875
calculatePRandom <- function(dataSet, sampleSize, p_pure, survivalTimes, followup, random.samples=10000) {

  geneNames <- colnames(dataSet)

  if (sampleSize == 1)
    sampleSize <- 2

  calculateRandomPPures <- function() {

    randomGeneNames <- sample(geneNames, sampleSize)

    randomGeneSet <- scale(dataSet[, randomGeneNames])

    return(calculatePPure(randomGeneSet, survivalTimes, followup)[["p_pure"]])

  }

  # calculate p_pures for randomly generated genesets
  p_pures <- replicate(random.samples, calculateRandomPPures(), simplify=TRUE)

  # p_random is the proportion of the random p_pures at least as significant
  # as the p_pure for the candidate geneset
  p_random <- sum(p_pures <= p_pure)/random.samples

  return (list("p_random"=p_random, "p_pures"=p_pures))

}


#' @export
#' @title Compute saps q-value
#' @description This function computes the saps q-value for a candidate prognostic
#' geneset by computing the saps score for randomly generated genesets and
#' determining the proportion at least as significant as the saps score for
#' the candidate set. This function is normally called by \code{\link{saps}}.
#' @param dataSet A matrix, where the column names are gene identifiers
#' and the values are gene expression levels. Each row should contain data for a
#' single patient.
#' @param sampleSize The desired size for the randomly sampled gene sets.
#' @param survivalTimes A vector of survival times. The length must equal
#' the number of rows in \code{dataSet}.
#' @param followup A vector of 0 or 1 values, indicating whether the patient was
#' lost to followup (0) or not (1). The length must equal the number of rows
#' (i.e. patients) in \code{dataSet}.
#' @param saps_score The saps score for the candidate geneset. The q-value is
#' calculated as the proportion of random saps scores whose absolute value
#' is >= to the provided saps score.
#' @param random.samples The number of random gene sets to sample during
#' the calculation of \code{P_random}.
#' @param qvalue.samples The number of random gene sets to sample for
#' purposes of computing the q-value.
#' @param cpus An integer that specifies the number of cpus/cores to be used when
#' calculating \code{P_enrichment}. If greater than 1, the
#' \pkg{snowfall} package must be installed or an error will occur.
#' @param gsea.perm The number of permutations to be used when calculating
#' \code{P_enrich}. This is passed to the \code{\link[piano]{runGSA}} function
#' in the \pkg{piano} package.
#' @param rankedGenes A vector of ranking scores for each gene in \code{dataSet}.
#' Ordinarily this will be the z-scores obtained by a call to
#' \code{\link{rankConcordance}}.
#' @return The function returns a list with two elements:
#'
#' \item{q_value}{the calculated q-value.}
#' \item{random_saps_scores}{a vector of individual saps scores for each
#' randomly generated geneset.}
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
#' # borderline significant saps score
#' saps_score <- 1.3
#'
#' \dontrun{
#' rankedGenes <- rankConcordance(dat, time, followup)[,"z"]
#'
#' q_value <- calculateQValue(dat, 5, time, followup, saps_score, random.samples=100,
#'      qvalue.samples=10, cpus=1, gsea.perm=1000, rankedGenes)
#'
#' q_value$q_value
#' random_scores <- abs(q_value$random_saps_scores)
#' hist(random_scores)
#' length(random_scores[random_scores > saps_score])
#' }
#' @seealso \code{\link{saps}}
#' @references Beck AH, Knoblauch NW, Hefti MM, Kaplan J, Schnitt SJ, et al.
#' (2013) Significance Analysis of Prognostic Signatures. PLoS Comput Biol 9(1):
#' e1002875.doi:10.1371/journal.pcbi.1002875
calculateQValue <- function(dataSet, sampleSize, survivalTimes, followup,
                            saps_score, random.samples, qvalue.samples,
                            cpus, gsea.perm, rankedGenes) {

  geneNames <- colnames(dataSet)

  if (sampleSize == 1)
    sampleSize <- 2

  calculateRandomSapsScores <- function() {

    # sample random geneset
    randomGeneNames <- sample(geneNames, sampleSize)

    randomGeneSet <- scale(dataSet[, randomGeneNames])

    # put random geneset into correct format for piano input
    random_name <- paste("random", sample(1:1000, 1), sep="")
    gsa_geneset <- matrix(randomGeneNames, nrow=1, dimnames=list(random_name))

    # calculate p_pure, p_random, p_enrich, saps_score
    p_pure <- calculatePPure(randomGeneSet, survivalTimes, followup)[["p_pure"]]

    p_random <- calculatePRandom(dataSet, sampleSize, p_pure, survivalTimes,
                               followup, random.samples)[["p_random"]]

    if (p_random == 0)
      p_random <- 1/(random.samples+1)

    gsa_results <- calculatePEnrichment(rankedGenes, gsa_geneset, cpus, gsea.perm)

    p_enrich <- gsa_results$P_enrichment
    direction <- gsa_results$direction

    if (p_enrich == 0)
      p_enrich <- 1/(gsea.perm+1)

    random_saps_score <- -log10(max(p_pure, p_random, p_enrich)) * direction

    return(random_saps_score)

  }

  # calculate saps scores for randomly generated genesets
  saps_scores <- replicate(qvalue.samples, calculateRandomSapsScores(),
                           simplify=TRUE)

  # q_value is the proportion of the random saps scores at least as significant
  # as the saps score for the candidate geneset
  q_value <- sum(abs(saps_scores) >= abs(saps_score))/qvalue.samples

  return(list("q_value"=q_value, "random_saps_scores"=saps_scores))

}



#' @export
#' @import survcomp
#' @title Compute concordance indices
#' @description Computes concordance indices for a gene expression data set,
#' and returns the concordance index and the z-score.
#' @details This function is a wrapper for
#' \code{\link[survcomp]{concordance.index}} in the \pkg{survcomp}
#' package. It applies the latter over the columns of \code{dataset} to
#' calculate concordance indices and the corresponding z-score for each gene.
#' @param dataset A matrix, where the column names are gene identifiers
#' and the values are gene expression levels. Each row should contain data for
#' a single patient.
#' @param survivalTimes A vector of survival times. The length must equal the number of
#' rows (i.e. patients) in \code{dataset}.
#' @param followup A vector of 0 or 1 values, indicating whether the patient was
#' lost to followup (0) or not (1). The length must equal the number of rows
#' (i.e. patients) in \code{dataSet}.
#' @return The function returns a matrix with two columns:
#'
#' \item{cindex}{concordance index estimate.}
#' \item{z}{z-score of the concordance index estimate.}
#'
#' and as many rows as \code{dataset}. The row names contain the gene identifiers.
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
#' ci <- rankConcordance(dat, time, followup)
#' z <- ci[,"z"]
#' range(z)
#' hist(z)
#' @seealso \code{\link{saps}} \code{\link[survcomp]{concordance.index}}
rankConcordance <- function(dataset, survivalTimes, followup) {

  concordance_f <- function(x) {

    tt <- survcomp::concordance.index(x, surv.time=survivalTimes,
                  surv.event=followup, method="noether", na.rm=TRUE)

    return (c("cindex"=tt$c.index, "z"=(tt$c.index - 0.5)/tt$se))

  }

  return (t(apply(dataset, 2, concordance_f)))

}


is.installed <- function(pkg) {
  return (is.element(pkg, installed.packages()[,1]))
}


save_saps <- function(geneSet, saps_unadj, saps_adj) {

  name <- geneSet[["name"]]

  size <- saps_unadj[name, "size"]
  direction <- saps_unadj[name, "direction"]

  p_pure <- saps_unadj[name, "p_pure"]
  p_random <- saps_unadj[name, "p_random"]
  p_enrich <- saps_unadj[name, "p_enrich"]
  saps_score <- saps_unadj[name, "saps_score"]
  saps_qvalue <- saps_unadj[name, "saps_qvalue"]

  p_pure_adj <- saps_adj[name, "p_pure"]
  p_random_adj <- saps_adj[name, "p_random"]
  p_enrich_adj <- saps_adj[name, "p_enrich"]
  saps_score_adj <- saps_adj[name, "saps_score"]
  saps_qvalue_adj <- saps_adj[name, "saps_qvalue"]

  geneSet["direction"] <- direction

  geneSet$saps_unadjusted["p_pure"] <- p_pure
  geneSet$saps_unadjusted["p_random"] <- p_random
  geneSet$saps_unadjusted["p_enrich"] <- p_enrich
  geneSet$saps_unadjusted["saps_score"] <- saps_score
  geneSet$saps_unadjusted["saps_qvalue"] <- saps_qvalue

  geneSet$saps_adjusted["p_pure"] <- p_pure_adj
  geneSet$saps_adjusted["p_random"] <- p_random_adj
  geneSet$saps_adjusted["p_enrich"] <- p_enrich_adj
  geneSet$saps_adjusted["saps_score"] <- saps_score_adj
  geneSet$saps_adjusted["saps_qvalue"] <- saps_qvalue_adj

  return(geneSet)

}
