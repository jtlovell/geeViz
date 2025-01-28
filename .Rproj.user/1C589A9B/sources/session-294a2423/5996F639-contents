#' @title Classify the genome into major groups
#' @description
#' \code{classify_genome} Neural network classification of large-scale genomic
#' features (telomeres, subtelomeric (aka arms), pericentromeric, centromeres,
#' and assembly gaps) trained on feature density, typically gene and repeats,
#' but could include recombination breakpoints or other information.
#'
#' @param featuresGr GRanges object giving the ranges of the features to use for
#' training. Must have a "class" or similar column with the identity of each
#' feature type. Processed hierarchically so that a single base is only ever
#' assigned to the first feature in the Granges object.
#' @param trainingGr Granges object giving positions to use for training. Must
#' have a metadata column 'grp' with at least two occurances of 'peric',
#' 'subtelo' and 'centromere'. If supplied, `centromeres` is ignored.
#' @param dnaSS DNAStringSet with the genomic sequence
#' @param width integer of length 1, specifying the window size in basepairs
#' @param step integer of length 1, specifying the step size in basepairs
#' @param teloKmers DNAStringSet or character vector of DNA kmers that are used
#' to identify the telomeres.
#' @param stripChrName regular expression used to process chromosome names. By
#' default, removes all text after the first whitespace.
#' @param teloBuffer integer specifying the distance from the telomeres where
#' the training region for chromosome arms should begin
#' @param centroBuffer integer specifying the distance from the centromeres
#' where the training region for pericentromeres should begin
#' @param minGeneDensity4train numeric giving the minimum proportion of genes
#' in the training regions for subtelomeric sequence
#' @param trainWidth integer specifying the width of the training intervals
#' @param minProb2assign numeric specifying the minimum assignment probability
#' to call a region non-ambiguous
#'
#' @details Neural network prediction of regions of interest. Explore the plot
#' and if certain regions appear to be problematic, exclude them and provide
#' another run with `trainingGr` with only the good regions.
#'
#' @return list of two elements: the first 'trained' gives the probability of
#' assignment for the training intervals; the second gives the "predicted" data
#' for each window in the genome.
#'
#' @examples
#' \dontrun{
#' # coming soon.
#' }
#'
#' @import GenomicRanges
#' @import data.table
#' @importFrom neuralnet neuralnet predict.nn
#' @importFrom GENESPACE add_rle
#' @export
classify_genome <- function(featuresGr,
                            trainingGr = NULL,
                            dnaSS,
                            centromeres,
                            maskGr = NULL,
                            width = 1e6,
                            step = 50e3,
                            teloKmers = c("CCCGAAA", "CCCTAAA"),
                            stripChrName = "\\s.*",
                            teloBuffer = 50e3,
                            centroBuffer = 1e6,
                            subteloMinGeneDensity = .1,
                            pericenMaxGeneDensity = 0.25,
                            trainWidth = width * 2,
                            minProb2assign = .75,
                            ignoreCentromeres = TRUE,
                            minRunLength = ceiling((0.5 *width)/step),
                            verbose = TRUE){

  md <- elementMetadata(featuresGr)[,1]
  featuresGr <- makeGRangesFromDataFrame(as.data.frame(featuresGr))
  featuresGr$type <- md

  # ignoreCentromeres <- TRUE
  ##############################################################################
  # 1. Make the training intervals for subtelomeres
  subtelos <- makeGRangesFromDataFrame(rbind(
    data.table(
      grp = "subteloLeft",
      seqnames = names(dnaSS),
      start = teloBuffer,
      end = teloBuffer + trainWidth),
    data.table(
      grp = "subteloRight",
      seqnames = names(dnaSS),
      start = width(dnaSS) - (teloBuffer + trainWidth),
      end = width(dnaSS) -  teloBuffer)), keep.extra.columns = T)

  ##############################################################################
  # 2. Make the training intervals for the pericentromeres
  if(!is.null(trainingGr)){
    trainingRegions <- trainingGr
  }else{
    tmp <- data.table(as.data.frame(centromeres))
    perics <- makeGRangesFromDataFrame(rbind(
      with(tmp, data.table(
        grp = "pericLeft",
        seqnames = seqnames,
        start = start - (centroBuffer + trainWidth),
        end = start - centroBuffer)),
      with(tmp, data.table(
        grp = "pericRight",
        seqnames = seqnames,
        start = end + centroBuffer,
        end = end + (centroBuffer + trainWidth)))), keep.extra.columns = T)

    # -- Join the training regions
    trainingRegions <- sort(c(subtelos, perics, cens))
  }

  ##############################################################################
  # 3. Build windows for the training sets
  # -- determine training window sizes/step (if smaller than the default window)
  if(width > min(width(trainingRegions))){
    widthTrain <- min(width(trainingRegions))
  }else{
    widthTrain <- width
  }

  stepTrain <- floor(widthTrain / 2)

  # -- window the training regions
  trainingWindows <- window_gr(
    x = trainingRegions,
    width = widthTrain,
    step = stepTrain)
  trainingWindows <<- trainingWindows
  featuresGr <<- featuresGr

  ##############################################################################
  # 4. Count the features in the training regions
  trainingCounts <- count_overlapsByGroup(
    regions = trainingWindows,
    features = featuresGr)

  ##############################################################################
  # 5. Count the features in sliding window regions
  # -- make the sliding windows
  grid <- window_gr(
    x = dnaSS,
    width = width,
    step = step)
  grid$region <- 1:length(grid)

  # -- count in each
  testCounts <- count_overlapsByGroup(
    regions = grid,
    features = featuresGr)

  ##############################################################################
  # 6. Build the training data
  # -- specify the levels in alphabetical order for consistency with nnet
  if(ignoreCentromeres){
    levs <- c("subtelo", "peric")
    levs <- levs[order(levs)]
  }else{
    levs <- c("subtelo", "peric", "centromere")
    levs <- levs[order(levs)]
  }

  # -- convert granges to data.table
  if(!is.null(maskGr)){
    fos <- findOverlaps(trainingCounts, maskGr)
    if(length(fos) > 0){
      trainingCounts <- trainingCounts[-unique(queryHits(fos))]
    }
  }

  trainDat <- data.table(as.data.frame(trainingCounts))

  # -- remove levt/right specification
  trainDat[,fac2predict := gsub("Right|Left", "", region)]
  trainDat <- subset(trainDat, fac2predict %in% levs)
  trainDat[,fac2predict := factor(fac2predict, levels = levs)]

  # -- reformat to wide
  trainDat <- dcast(
    trainDat,
    seqnames + start + end + fac2predict ~ grp,
    value.var = "prop")

  ##############################################################################
  # 7. Drop values outside of expected gene density ranges
  # -- if minGeneDensity4train, subset these out
  if(subteloMinGeneDensity > 0){
    cn <- colnames(trainDat)
    wh <- cn[grepl("CDS|GENE|MRNA|INTRON", toupper(cn))]
    trainDat[,sumGene := rowSums(.SD), .SDcols = wh]
    trainDat <- subset(
      trainDat, sumGene >= subteloMinGeneDensity | fac2predict != "subtelo")
  }

  if(pericenMaxGeneDensity < 1){
    cn <- colnames(trainDat)
    wh <- cn[grepl("CDS|GENE|MRNA|INTRON", toupper(cn))]
    trainDat[,sumGene := rowSums(.SD), .SDcols = wh]
    trainDat <- subset(
      trainDat, sumGene <= pericenMaxGeneDensity | fac2predict != "peric")
  }

  ##############################################################################
  # 8. Build the neural network model
  u <- unique(trainingCounts$grp)
  form <- sprintf("fac2predict ~ %s", paste(u, collapse = " + "))
  nnModel <- neuralnet(
    as.formula(form),
    data = trainDat,
    # hidden = c(4,2),
    linear.output = FALSE)

  ##############################################################################
  # 9. Build the test data

  if(!is.null(maskGr)){
    fos <- findOverlaps(testCounts, maskGr)
    if(length(fos) > 0){
      testCounts <- testCounts[-unique(queryHits(fos))]
    }
  }

  testDat <- data.table(as.data.frame(testCounts))

  # -- reformat to wide
  testDat <- dcast(
    testDat,
    seqnames + start + end ~ grp,
    value.var = "prop")

  ##############################################################################
  # 10. Do the prediction on the train data, return correspondence
  conf <- data.table(predict(nnModel, trainDat))
  setnames(conf, levs)
  conf <- data.table(trainDat, conf)
  conf[,bestAssign := levs[apply(.SD, 1, which.max)],
       .SDcols = levs]

  ##############################################################################
  # 11. Do the prediction on the test data
  pred <- data.table(predict(nnModel, testDat))
  setnames(pred, levs)
  predout <- data.table(testDat, pred)

  ##############################################################################
  # 12. Remove anything below the best assignment probability threshold
  predout[,maxProb := apply(.SD, 1, max), .SDcols = levs]
  predout[,bestAssign := levs[apply(.SD, 1, which.max)], .SDcols = levs]
  predout$bestAssign[predout$maxProb < minProb2assign] <- "ambiguous"

  ##############################################################################
  # 13. convert to run-length eqivs
  torl <- data.table(as.data.frame(subset(predout, bestAssign != "ambiguous")))
  setkey(torl, seqnames, start, end)
  torl[,pos := (start + end)/2]
  for(i in 2:minRunLength){
    torl[,rl := add_rle(bestAssign), by = "seqnames"]
    torl <- subset(torl, rl >= i)
  }
  torl[,blkID := add_rle(bestAssign, which = "id"), by = "seqnames"]

  blks <- torl[,list(start = min(start), end = max(end)),
               by = c("seqnames", "blkID", "bestAssign")]

  ##############################################################################
  # 14. for overlapping blocks, choose the midpoint as the break and reformat
  blks[,offLeft := c(start[-1] - end[-.N], NA)/2, by = "seqnames"]
  blks[,offRight := c(NA, end[-.N] - start[-1])/2, by = "seqnames"]
  blks[is.na(blks)] <- 0
  blks[,`:=`(start = ceiling(blks$start + offRight),
             end = floor(blks$end + offLeft))]

  blks <- blks[,c("seqnames", "start", "end", "bestAssign")]

  return(list(trained = conf, predicted = predout, blks = blks))
}

