% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classify_genome.R
\name{classify_genome}
\alias{classify_genome}
\title{Classify the genome into major groups}
\usage{
classify_genome(
  featuresGr,
  trainingGr = NULL,
  dnaSS,
  centromeres,
  maskGr = NULL,
  width = 1e+06,
  step = 50000,
  teloKmers = c("CCCGAAA", "CCCTAAA"),
  stripChrName = "\\\\s.*",
  teloBuffer = 50000,
  centroBuffer = 1e+06,
  subteloMinGeneDensity = 0.1,
  pericenMaxGeneDensity = 0.25,
  trainWidth = width * 2,
  minProb2assign = 0.75,
  ignoreCentromeres = TRUE,
  minRunLength = ceiling((0.5 * width)/step),
  verbose = TRUE
)
}
\arguments{
\item{featuresGr}{GRanges object giving the ranges of the features to use for
training. Must have a "class" or similar column with the identity of each
feature type. Processed hierarchically so that a single base is only ever
assigned to the first feature in the Granges object.}

\item{trainingGr}{Granges object giving positions to use for training. Must
have a metadata column 'grp' with at least two occurances of 'peric',
'subtelo' and 'centromere'. If supplied, `centromeres` is ignored.}

\item{dnaSS}{DNAStringSet with the genomic sequence}

\item{width}{integer of length 1, specifying the window size in basepairs}

\item{step}{integer of length 1, specifying the step size in basepairs}

\item{teloKmers}{DNAStringSet or character vector of DNA kmers that are used
to identify the telomeres.}

\item{stripChrName}{regular expression used to process chromosome names. By
default, removes all text after the first whitespace.}

\item{teloBuffer}{integer specifying the distance from the telomeres where
the training region for chromosome arms should begin}

\item{centroBuffer}{integer specifying the distance from the centromeres
where the training region for pericentromeres should begin}

\item{trainWidth}{integer specifying the width of the training intervals}

\item{minProb2assign}{numeric specifying the minimum assignment probability
to call a region non-ambiguous}

\item{minGeneDensity4train}{numeric giving the minimum proportion of genes
in the training regions for subtelomeric sequence}
}
\value{
list of two elements: the first 'trained' gives the probability of
assignment for the training intervals; the second gives the "predicted" data
for each window in the genome.
}
\description{
\code{classify_genome} Neural network classification of large-scale genomic
features (telomeres, subtelomeric (aka arms), pericentromeric, centromeres,
and assembly gaps) trained on feature density, typically gene and repeats,
but could include recombination breakpoints or other information.
}
\details{
Neural network prediction of regions of interest. Explore the plot
and if certain regions appear to be problematic, exclude them and provide
another run with `trainingGr` with only the good regions.
}
\examples{
\dontrun{
# coming soon.
}

}
