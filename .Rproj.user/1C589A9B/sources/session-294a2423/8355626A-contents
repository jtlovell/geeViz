#' @title Generic internal functions used by geeviz
#' @description
#' \code{kmerUtils} Convenience functions for geeviz, not meant to be called
#' directly by the user. Little documentation support provided, use at your own
#' risk.
#' @name kmerUtils
#'
#' @param x likely a genomicRanges object, but check the specific function
#' @param gff3File file.path of length 1, specifying the gff3 file to use
#' @param extraCol2keep character vector of extra columns to keep
#' @param reduceIt logical, should overlapping ranges be reduced?
#' \cr
#' If called, \code{kmerUtils} returns its own arguments.
#'
#'

#' @title pull kmer blocks
#' @description
#' \code{pull_kmerBlocks} generic function to pull granges of kmer exact hits to
#' an assembly
#' @rdname kmerUtils
#' @import GenomicRanges
#' @importFrom Biostrings matchPDict PDict
#' @export
pull_kmerBlocks <- function(kmerSS,
                            assemblySS,
                            maxGap2join = 100){

  if(length(unique(width(kmerSS))) > 1){
    minSeqLen <- min(width(kmerSS))
    kmerSS <- kmerize_ss(ss = kmerSS, width = minSeqLen)
    kmerSS <- kmerSS[!duplicated(kmerSS)]
  }
  names(kmerSS) <- 1:length(kmerSS)
  kmerLen <- width(kmerSS)[1]
  pdictKmer <- PDict(kmerSS)

  kmerLocs <- lapply(names(assemblySS), function(i){
    mdict <- makeGRangesFromDataFrame(data.table(
      chr = i,
      as.data.frame(unlist(matchPDict(
        pdict = pdictKmer,
        subject = assemblySS[[i]])))[,1:2]))

    out <- join_gappedRanges(
      x = mdict,
      min.gapwidth = maxGap2join)
    return(out)
  })
  suppressWarnings(out <- BiocGenerics::do.call(c, kmerLocs))
  return(out)
}
