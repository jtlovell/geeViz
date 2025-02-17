% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/window_gr.R
\name{window_gr}
\alias{window_gr}
\title{flexible sliding window constructor}
\usage{
window_gr(x, width, step, minWindowSize = width/2)
}
\arguments{
\item{x}{GRanges, DNAStringSet or SeqInfo object. The start and end
coordinates for the windows are taken from this.}

\item{width}{integer of length 1, specifying the window size in basepairs}

\item{step}{integer of length 1, specifying the step size in basepairs}

\item{minWindowSize}{numeric of length 1, specifying the minimum window size
to be retained. Defaults to 50% of the width.}
}
\value{
GRanges object with the window ranges without considering strand
}
\description{
\code{window_gr} Build a granges sliding window from a seqinfo, stringset or
granges object. Identical output to GenomicRanges::slidingWindows(), except
allowing non-overlapping windows.
}
\details{
Constructs the windows in data.table
}
\examples{
\dontrun{
# coming soon.
}

}
