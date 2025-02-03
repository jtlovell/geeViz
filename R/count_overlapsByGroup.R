#' @title count overlapping bases
#'
#' @description
#' \code{count_overlapsByGroup} Hierarchically classify the genome by coordinates in
#' the list of bed files. Sequences in the first bed file are masked in the
#' second and so on.
#'
#' @param regions GRanges object with at least one elementMetadata column. Ranges
#' contain the bounds in which ranges are counted (e.g. windows, blocks, etc.)
#' @param features GRanges object with at least one elementMetadata column. Ranges
#' regions to be summed (e.g. genes, repeats, etc.)
#'
#' @details intersects, gaps and reduces overlapping intervals
#'
#' @return a list of granges matchin the listOfBeds
#'
#' @examples
#' \dontrun{
#' # coming soon.
#' }
#'
#' @import GenomicRanges
#' @importFrom data.table data.table foverlaps CJ setnames setkey :=
#' @importFrom stats complete.cases
#' @export
count_overlapsByGroup <- function(features, regions){

  x <- features
  y <- regions

  if(ncol(elementMetadata(y)) == 0){
    warning("no classification (elementMetadata) for regions, adding one per range\n")
    y$region <- 1:length(y)
  }

  if(ncol(elementMetadata(x)) == 0){
    warning("no classification (elementMetadata) for features, adding one per range\n")
    x$class <- "noclass"
  }

  # -- format features granges, condensing overlapping ranges within groups
  x <- reduce_grByGroup(x)

  # -- enforce hierarchical so that no group overlaps another group
  x <- make_hierarchical(x)

  # -- convert features to data.table
  dtx <- data.table(
    chr = as.character(seqnames(x)),
    start = start(x),
    end = end(x),
    classx = elementMetadata(x)[[1]],
    key = c("chr", "start", "end"))

  # -- convert regions to data.table, adding an 'index' for each unique region
  dty <- data.table(
    chr = as.character(seqnames(y)),
    start = start(y),
    end = end(y),
    classy = elementMetadata(y)[[1]],
    index = sprintf("region%s", 1:length(y)),
    key = c("chr", "start", "end"))

  # -- overlap join the regions and features
  fo <- foverlaps(dtx, dty)
  # print(subset(dty, index == "region184"))
  # print(dtx)
  fo <- subset(fo, complete.cases(fo))

  # -- do the truncation of the ranges
  i.start <- i.end <- start <- end <- NULL
  whs <- with(fo, i.start < start)
  whe <- with(fo, i.end > end)
  fo$i.start[whs] <- fo$start[whs]
  fo$i.end[whe] <- fo$end[whe]

  # -- ensure that they are no duplicates
  # classy <- classx <- NULL
  # fogr <- makeGRangesFromDataFrame(fo)
  # fogr$u <- with(fo, paste(index, classx))

  # -- get counts by each grouping
  i.start <- i.end <- NULL
  foo <- fo[,list(nbp = sum(i.end - i.start)),
            by = c("chr", "start", "end", "classy", "classx", "index")]

  # -- ensure all unique combinations of classes are in there
  cls <- as.character(unique(dtx$classx))
  ucomb <- dty[,list(classx = cls, nbp = 0),
               by = c("chr", "start", "end", "classy", "index")]
  out <- rbind(foo, ucomb)
  out <- subset(out, !duplicated(out[,1:6]))
  setkey(out, chr, start, end, classx, classy)

  # -- add 0's where no values exist
  out$nbp[is.na(out$nbp)] <- 0

  # -- get "missing" counts
  nbp <- end <- start <- width <- tot <- NULL
  out[,tot := sum(nbp), by = c("classy", "chr", "start", "end", "index")]
  out[,width := 1 + (end[1] - start[1]), by = c("classy", "chr", "start", "end", "index")]
  out[,missing := width - tot]
  out[,`:=`(tot = NULL, width = NULL)]

  # -- add missing counts back to the main object
  missing <- classy <- NULL
  outm <- subset(out, !duplicated(index))
  outm[,`:=`(nbp = missing, classx = "missing")]
  outm[,missing := NULL]
  out[,missing := NULL]
  out <- rbind(outm, out)

  # -- calculate proportion membership
  out[,tot := sum(nbp), by = "index"]
  out[,prop := nbp / tot]
  out[,tot := NULL]

  # -- re-order to follow genome organization
  chr <- start <- end <- classy <- classx <- NULL
  setkey(out, chr, start, end, classy, classx)

  # -- rename grouping variables by input
  setnames(
    out,
    c("classy", "classx"),
    c(names(elementMetadata(y))[1],
      names(elementMetadata(x))[1]))

  out[,uniqueRegionID := index]
  out[,index := NULL]

  # -- covert to granges and return
  out <- makeGRangesFromDataFrame(
    out,
    keep.extra.columns = TRUE,
    ignore.strand = TRUE)
  return(out)
}
