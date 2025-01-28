## This is geeViz - or whatever its going to be called - v.0.0.1

# Do not use this package!!

This release is a placeholder until a full package+documentation can be built

### Classify genome

This function splits genomic sequences into two bins: pericentromere and chromosome arms (aka subtelomeres). Future releases will:

1. Cluster centromeres and extract the position of gaps and telomeres
2. Integrate with DEEPSPACE riparian plots and GENESPACE contig maps. The latter will likely migrate to geeViz in the first full release
3. Make the below script an easy executable
4. XXX add more to do's

Installation requires the normal call to github:

```
devtools::install_github("jtlovell/geeViz")
library(geeViz)
```

To run the commands below, we also need some other packages, which should install along with `geeViz`

```
library(Biostrings)
library(GenomicRanges)
library(data.table)
library(ggplot2)
```

To run `classify_genome`, you need the following:

#### Some parameters

```
# regex to strip off text after and including the first whitespace
chrNamesStrip <- "\\s.*" 

```

#### Assembly fasta file, read in as a DNAStringSet object

```
assemFaFile <- "/path/to/assembly_v1.0.fa.gz"
ss <- Biostrings::readDNAStringSet(assemFaFile)
names(ss) <- gsub(chrNamesStrip, "", names(ss))
```

#### Protein-coding gene annotations file, read in as a Granges object

Optionally can split out introns and exons. Here we do this, then make the Granges hierarchical so than any sequence within CDS is masked in the gene types, converting those to introns/UTRs

```
geneGff3File <- "/path/to/genomeAnnotation.gene.gff3.gz"
geneGr <- read_gffAsGr(
  gff3File = geneGff3File,
  filter = list(type = c("CDS", "gene")),
  reduceIt = TRUE, 
  extraCol2keep = "type")
geneGr$type[geneGr$type == "gene"] <- "intron"
```

#### Repeat gff file, read in as a Granges object

```
# -- 3. read in the repeats
repGr <- read_gffAsGr(
  gff3File = repGff3File,
  extraCol2keep = "class",
  reduceIt = TRUE)
repGr$type <- repGr$class
repGr$class <- NULL
```

#### Combine the granges and organize 

Same as for the genes, but more levels that again need to be hierarchical. Also complicated because we don't want to train the model on very rare elements and the heirarchy needs to have "unknown" repeats at the end. 

```
feats <- c(geneGr, repGr)
featureOrder <- unique(feats$type)

wh <- c(grep("CDS|EXON", toupper(featureOrder)),
        grep("CENK", toupper(featureOrder)),
        grep("CDS|EXON|CENK|UNKNOWN|INTRON", toupper(featureOrder), invert = T),
        grep("INTRON", toupper(featureOrder)),
        grep("UNKNOWN", toupper(featureOrder)))
wh <- wh[!duplicated(wh)]
if(length(wh) < length(featureOrder)){
  al <- 1:length(featureOrder)
  mis <- al[!al %in% wh]
  wh <- c(wh, mis)
}
featureOrder <- featureOrder[wh]
  
feats$type <- factor(feats$type, levels = featureOrder)
feats <- feats[order(feats$type),]
feats$type <- as.character(feats$type)
feats$type <- gsub("[^[:alnum:] ]", "", feats$type)

featsH <- make_hierarchical(feats)
```

#### Finally, read in the centromeres for training:

```
cenRegFile <- "/path/to/entromere_candidate_regions.csv"
cenRegs <- fread(
  cenRegFile,
  select = c(1:4), col.names = c("seqnames", "start", "end", "genome"))
cenRegs[,width := end - start]
setorder(cenRegs, genome, seqnames, -width)
cenRegs <- subset(cenRegs, !duplicated(paste(genome, seqnames)))
cens <- makeGRangesFromDataFrame(subset(cenRegs, genome == "MN106"))
cens$grp <- "centromere"
```

Now we are ready to run classify_genome:

```
classd <- classify_genome(
  featuresGr = featsH,
  trainingGr = NULL,
  dnaSS = ss,
  # maskGr = cens2mask,
  centromeres = cens,
  width = 200e3,
  step = 10e3, 
  teloBuffer = 10e3, 
  centroBuffer = 1e6,
  trainWidth = 3e6,
  minProb2assign = .75, 
  ignoreCentromeres = T)
```

For viz, we slide a window across the features
```
linChrs <- linearize_chrCoords(
  x = ss,
  stripChrName = "//s.*",
  gapSize = 5e6)
  
grid <- window_gr(x = ss, width = 500e3, step = 100e3)
grid$grp <- 1
feats2p <- feats
feats2p$type <- ifelse(grepl("Copia", feats2p$type), "Copia",
                     ifelse(grepl("Gypsy", feats2p$type), "Ty3",
                            ifelse(grepl("CDS|intron", feats2p$type), "gene", "otherRepeat")))
feats2p$type <- factor(feats2p$type, levels = c("gene", "otherRepeat", "Ty3", "Copia"))

grid$grp <- NULL
sw <- count_overlapsByGroup(regions = grid, features = feats2p)
dt <- data.table(as.data.frame(sw))
dt[,mid := (start + end)/2]
xv <- linChrs$xStart; names(xv) <- linChrs$chr
dt[,x := mid + xv[seqnames]]
dt[,class := factor(grp, levels = c("gene", "missing", "otherRepeat", "Ty3", "Copia"))]

blks <- classd$blks
blks[,`:=`(xStart = xv[seqnames] + start, xEnd = xv[seqnames] + end)]
setkey(dt, class, seqnames, mid)
dt[,grp2 := factor(paste(seqnames, class), levels = unique(paste(seqnames, class)))]
```


Then we can plot the output:

```
p <- ggplot()+
  geom_rect(data = linChrs, aes(xmin = xStart, xmax = xEnd, ymin = -0.2, ymax = -0.01))+
  geom_rect(data = subset(blks, bestAssign == "subtelo"),
            aes(xmin = xStart, xmax = xEnd, ymin = -0.18, ymax = -0.03),
            col = NA, fill = "darkorange")+
  geom_rect(data = subset(blks, bestAssign == "peric"),
            aes(xmin = xStart, xmax = xEnd, ymin = -0.18, ymax = -0.03),
            col = NA, fill = "darkblue")+
  scale_y_continuous(limits = c(-.2, 1), expand = c(0,0))+
  scale_x_continuous(expand = c(.001, .001), breaks = linChrs$mid, labels = linChrs$chr)+
  geom_area(data = dt, aes(x = x, y = prop, group = grp2, fill = class))+
  scale_fill_manual(values = c("darkorange3", "grey95", "skyblue", "dodgerblue", "dodgerblue4"))+
  GENESPACE::theme_genespace()
```

