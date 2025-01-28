#' @title Generic internal functions used by geeviz
#' @description
#' \code{utils} Convenience functions for geeviz, not meant to be called
#' directly by the user. Little documentation support provided, use at your own
#' risk.
#' @name utils
#'
#' @param x likely a genomicRanges object, but check the specific function
#' @param gff3File file.path of length 1, specifying the gff3 file to use
#' @param extraCol2keep character vector of extra columns to keep
#' @param reduceIt logical, should overlapping ranges be reduced?
#' \cr
#' If called, \code{utils} returns its own arguments.
#'
#'

#' @title kmerize stringset
#' @description
#' \code{kmerize_ss} Split sequences into 1bp offset sliding windows
#' @rdname utils
#' @import GenomicRanges
#' @importFrom Biostrings seqinfo
#' @export
kmerize_ss <- function(ss, width = 10, dropDuplicates = TRUE){
  if(is.null(names(ss)))
    names(ss) <- 1:length(ss)
  gr <- GRanges(seqinfo(ss))
  sw <- unlist(slidingWindows(gr, width = width, step = 1))
  sso <- ss[sw]
  sso <- sso[width(sso) == width]
  if(dropDuplicates)
    sso <- sso[!duplicated(sso)]
  names(sso) <- NULL
  return(sso)
}

#' @title join adjacent ranges
#' @description
#' \code{join_gappedRanges} Reduces ranges, returning a granges with proportion
#' of sequence and total sequence covered.
#' @rdname utils
#' @import GenomicRanges
#' @export
join_gappedRanges <- function(x, min.gapwidth = 0){
  gr <- sort(x)

  gr <- reduce(
    gr,
    drop.empty.ranges = TRUE,
    ignore.strand = TRUE)

  out <- reduce(
    gr,
    min.gapwidth = min.gapwidth,
    drop.empty.ranges = TRUE,
    ignore.strand = TRUE,
    with.revmap = TRUE)

  out1 <- subset(out, sapply(out$revmap, length) == 1)
  out1$matchBases <- width(out1)
  out1$matchProp <- 1

  out2 <- subset(out, sapply(out$revmap, length) > 1)
  if(length(out2) > 0){
    rmap <- as.list(out2$revmap)
    wds <- sapply(rmap, function(x) sum(width(gr[x])))
    out2$matchBases <- wds
    out2$matchProp <- wds/width(out2)

    out1 <- c(out1, out2)
  }

  out1$revmap <- NULL
  return(out1)
}

#' @title Reduce granges by group
#' @description
#' \code{reduce_grByGroup} Join overlapping ranges within a group
#' @rdname utils
#' @import GenomicRanges
#' @export
reduce_grByGroup <- function(x){
  tmp <- x
  splb <- elementMetadata(x)[[1]]
  spl <- split(x, splb)
  x <- unlist(reduce(spl))
  ren <- names(elementMetadata(tmp))[1]
  x$tmp <- names(ranges(x))
  names(attr(x, "elementMetadata")) <- ren
  names(ranges(x)) <- NULL
  return(x)
}

#' @title Make granges hierarchical
#' @description
#' \code{make_hierarchical} Ensure ranges present in the previous group are
#' masked
#' @rdname utils
#' @import GenomicRanges
#' @importFrom methods as
#' @export
make_hierarchical <- function(x){

  if(ncol(elementMetadata(x)) == 0)
    return(x)

  grpDat <- elementMetadata(x)[,1]
  if(length(unique(grpDat)) == 1)
    return(x)

  if(is.factor(grpDat)){
    u <- levels(grpDat)
  }else{
    u <- unique(grpDat)
  }

  x <- makeGRangesFromDataFrame(as.data.frame(x))
  x$grp <- as.character(grpDat)

  spl <- split(x, x$grp)[u]
  spl <- spl[sapply(spl, length) > 0]
  out <- spl[[1]]
  grps <- names(spl)
  for(i in grps[-1]){
    xi <- spl[[i]]
    tmp <-  setdiff(xi, out)
    if(length(tmp) > 0){
      tmp$grp <- i
      out <- c(out, tmp)
    }
  }
  out$grp <- factor(out$grp, levels = grps)
  out <- out[order(out$grp),]
  return(out)
}

#' @title Read a gff3 file into GRanges
#' @description
#' \code{read_gffAsGr} Import a gff3-formatted annotation file into GRanges
#' @rdname utils
#' @import GenomicRanges
#' @importFrom rtracklayer readGFF
#' @export
read_gffAsGr <- function(gff3File,
                         extraCol2keep = NULL,
                         reduceIt = FALSE,
                         ...){
  tmp <- as.data.frame(readGFF(gff3File, ...))
  if(!is.null(extraCol2keep)){
    if(all(extraCol2keep %in% colnames(tmp))){
      tmp <- tmp[,c("seqid", "start", "end", extraCol2keep)]
    }else{
      tmp <- tmp[,c("seqid", "start", "end")]
    }
  }else{
    tmp <- tmp[,c("seqid", "start", "end")]
  }

  gr <- makeGRangesFromDataFrame(
    tmp,
    keep.extra.columns = TRUE,
    ignore.strand = TRUE)
  if(reduceIt){
    if(length(names(attr(gr, "elementMetadata"))) == 1){
      gr <- reduce_grByGroup(gr)
    }else{
      gr <- reduce(gr)
    }
  }
  return(gr)
}

#' @title Calculate the linear positions of chromosomes
#' @description
#' \code{linearize_chrCoords} Produce a data.table with the linear positions for
#' plotting
#' @rdname utils
#' @import GenomicRanges
#' @importFrom Biostrings fasta.seqlengths
#' @importFrom GenomeInfoDb seqlengths
#' @export
linearize_chrCoords <- function(x, gapSize = 0, stripChrName = ""){

  # -- get the input data formats into a vector
  if(is(x, "DNAStringSet")){
    si <- width(x)
    names(si) <- names(x)
  }else{
    if(is(x, "Seqinfo")){
      si <- seqlengths(x)
    }else{
      if(is(x, "character") && length(x) == 1){
        if(file.exists(x)){
          si <- fasta.seqlengths(x)
        }else{
          stop("cannot find file: ", x, "\n")
        }
      }else{
        stop("x must be a DNAStringSet, Seqinfo, or path to a fasta file\n")
      }
    }
  }

  # -- convert to data table
  dt <- data.table(chr = names(si), width = si)

  # -- get the starts and ends
  dt[,gp := gapSize]
  dt[,xStart := c(1, cumsum(width[-.N]) + cumsum(gp[-.N]) + 1)]
  dt[,`:=`(xEnd = xStart + width, gp = NULL, mid = (xStart + width + xStart)/2)]

  # -- rename chrs if necessary
  if(stripChrName != ""){
    dt[,chr := gsub(stripChrName, "", chr)]
  }
  return(dt)
}
