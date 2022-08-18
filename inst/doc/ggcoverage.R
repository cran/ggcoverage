## ----style, echo = FALSE, results = 'asis'------------------------------------
BiocStyle::markdown()

## ----setup, echo=FALSE, warning=FALSE-----------------------------------------
library(knitr)
htmltools::tagList(rmarkdown::html_dependency_font_awesome())

# set dpi
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dpi=60
)

## ----install, eval=FALSE------------------------------------------------------
#  # install via CRAN
#  install.package("ggcoverage")
#  
#  # install via Github
#  # install.package("remotes")   #In case you have not installed it.
#  remotes::install_github("showteeth/ggcoverage")

## ----library, message=FALSE, warning=FALSE------------------------------------
library("rtracklayer")
library("graphics")
library("ggcoverage")

## ----load_metadata------------------------------------------------------------
# load metadata
meta.file <- system.file("extdata", "RNA-seq", "meta_info.csv", package = "ggcoverage")
sample.meta = read.csv(meta.file)
sample.meta

## ----load_track---------------------------------------------------------------
# track folder
track.folder = system.file("extdata", "RNA-seq", package = "ggcoverage")
# load bigwig file
track.df = LoadTrackFile(track.folder = track.folder, format = "bw",
                         meta.info = sample.meta)
# check data
head(track.df)

## ----prepare_mark-------------------------------------------------------------
# create mark region
mark.region=data.frame(start=c(21678900,21732001,21737590),
                       end=c(21679900,21732400,21737650),
                       label=c("M1", "M2", "M3"))
# check data
mark.region

## ----load_gtf-----------------------------------------------------------------
gtf.file = system.file("extdata", "used_hg19.gtf", package = "ggcoverage")
gtf.gr = rtracklayer::import.gff(con = gtf.file, format = 'gtf')

## ----basic_coverage, eval=FALSE-----------------------------------------------
#  basic.coverage = ggcoverage(data = track.df, color = "auto",
#                              mark.region = mark.region, range.position = "out")
#  basic.coverage

## ----basic_coverage_plot, echo=FALSE, fig.height = 6, fig.width = 12, fig.align = "center"----
knitr::include_graphics("../man/figures/README-basic_coverage-1.png")

## ----basic_coverage_2, eval=FALSE---------------------------------------------
#  basic.coverage = ggcoverage(data = track.df, color = "auto",
#                              mark.region = mark.region, range.position = "in")
#  basic.coverage

## ----basic_coverage_2_plot, echo=FALSE, fig.height = 6, fig.width = 12, fig.align = "center"----
knitr::include_graphics("../man/figures/README-basic_coverage_2-1.png")

## ----gene_coverage, eval=FALSE------------------------------------------------
#  basic.coverage +
#    geom_gene(gtf.gr=gtf.gr)

## ----gene_coverage_plot, echo=FALSE, fig.height = 8, fig.width = 12, fig.align = "center"----
knitr::include_graphics("../man/figures/README-gene_coverage-1.png")

## ----transcript_coverage, eval=FALSE------------------------------------------
#  basic.coverage +
#    geom_transcript(gtf.gr=gtf.gr,label.vjust = 1.5)

## ----transcript_coverage_plot, echo=FALSE, fig.height = 12, fig.width = 12, fig.align = "center"----
knitr::include_graphics("../man/figures/README-transcript_coverage-1.png")

## ----ideogram_coverage_1, eval=FALSE------------------------------------------
#  basic.coverage +
#    geom_gene(gtf.gr=gtf.gr) +
#    geom_ideogram(genome = "hg19",plot.space = 0)

## ----ideogram_coverage_1_plot, echo=FALSE, fig.height = 10, fig.width = 12, fig.align = "center"----
knitr::include_graphics("../man/figures/README-ideogram_coverage_1-1.png")

## ----ideogram_coverage_2, eval=FALSE------------------------------------------
#  basic.coverage +
#    geom_transcript(gtf.gr=gtf.gr,label.vjust = 1.5) +
#    geom_ideogram(genome = "hg19",plot.space = 0)

## ----ideogram_coverage_2_plot, echo=FALSE, fig.height = 14, fig.width = 12, fig.align = "center"----
knitr::include_graphics("../man/figures/README-ideogram_coverage_2-1.png")

## ----load_bin_counts----------------------------------------------------------
# track file
track.file = system.file("extdata", "DNA-seq", "CNV_example.txt", package = "ggcoverage")
track.df = read.table(track.file, header = TRUE)
# check data
head(track.df)

## ----basic_coverage_dna, eval=FALSE-------------------------------------------
#  basic.coverage = ggcoverage(data = track.df,color = NULL, mark.region = NULL,
#                              region = 'chr4:61750000-62,700,000', range.position = "out")
#  basic.coverage

## ----basic_coverage_dna_plot, echo=FALSE, fig.height = 6, fig.width = 12, fig.align = "center"----
knitr::include_graphics("../man/figures/README-basic_coverage_dna-1.png")

## ----gc_coverage, eval=FALSE--------------------------------------------------
#  # load genome data
#  library("BSgenome.Hsapiens.UCSC.hg19")
#  # create plot
#  basic.coverage +
#    geom_gc(bs.fa.seq=BSgenome.Hsapiens.UCSC.hg19) +
#    geom_gene(gtf.gr=gtf.gr) +
#    geom_ideogram(genome = "hg19")

## ----gc_coverage_plot, echo=FALSE, fig.height = 10, fig.width = 12, fig.align = "center"----
knitr::include_graphics("../man/figures/README-gc_coverage-1.png")

## ----load_single_nuc----------------------------------------------------------
# prepare sample metadata
sample.meta <- data.frame(
  SampleName = c("tumorA.chr4.selected"),
  Type = c("tumorA"),
  Group = c("tumorA")
)
# load bam file
bam.file = system.file("extdata", "DNA-seq", "tumorA.chr4.selected.bam", package = "ggcoverage")
track.df <- LoadTrackFile(
  track.file = bam.file,
  meta.info = sample.meta,
  single.nuc=TRUE, single.nuc.region="chr4:62474235-62474295"
)
head(track.df)

## ----base_color_scheme, warning=FALSE, fig.height = 2, fig.width = 6, fig.align = "center"----
# color scheme
nuc.color = c("A" = "#ff2b08", "C" = "#009aff", "G" = "#ffb507", "T" = "#00bc0d")
opar <- graphics::par() 
# create plot
graphics::par(mar = c(1, 5, 1, 1))
graphics::image(
  1:length(nuc.color), 1, as.matrix(1:length(nuc.color)),
  col = nuc.color,
  xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n"
)
graphics::text(1:length(nuc.color), 1, names(nuc.color))
graphics::mtext(
  text = "Base", adj = 1, las = 1,
  side = 2
)

# reset par default
graphics::par(opar)

## ----aa_color_scheme, warning=FALSE, fig.height = 9, fig.width = 10, fig.align = "center"----
aa.color = c(
  "D" = "#FF0000", "S" = "#FF2400", "T" = "#E34234", "G" = "#FF8000", "P" = "#F28500",
  "C" = "#FFFF00", "A" = "#FDFF00", "V" = "#E3FF00", "I" = "#C0FF00", "L" = "#89318C",
  "M" = "#00FF00", "F" = "#50C878", "Y" = "#30D5C8", "W" = "#00FFFF", "H" = "#0F2CB3",
  "R" = "#0000FF", "K" = "#4b0082", "N" = "#800080", "Q" = "#FF00FF", "E" = "#8F00FF",
  "*" = "#FFC0CB", " " = "#FFFFFF", " " = "#FFFFFF", " " = "#FFFFFF", " " = "#FFFFFF"
)

graphics::par(mar = c(1, 5, 1, 1))
graphics::image(
  1:5, 1:5, matrix(1:length(aa.color),nrow=5),
  col = rev(aa.color),
  xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n"
)
graphics::text(expand.grid(1:5,1:5), names(rev(aa.color)))
graphics::mtext(
  text = "Amino acids", adj = 1, las = 1,
  side = 2
)

# reset par default
graphics::par(opar)

## ----base_aa_coverage, eval=FALSE---------------------------------------------
#  ggcoverage(data = track.df, color = "grey", range.position = "out", single.nuc=T, rect.color = "white") +
#    geom_base(bam.file = bam.file,
#              bs.fa.seq = BSgenome.Hsapiens.UCSC.hg19) +
#    geom_ideogram(genome = "hg19",plot.space = 0)

## ----base_aa_coverage_plot, echo=FALSE, fig.height = 10, fig.width = 12, fig.align = "center"----
knitr::include_graphics("../man/figures/README-base_aa_coverage-1.png")

## ----load_metadata_chip-------------------------------------------------------
# load metadata
sample.meta = data.frame(SampleName=c('Chr18_MCF7_ER_1','Chr18_MCF7_ER_2','Chr18_MCF7_ER_3','Chr18_MCF7_input'),
                         Type = c("MCF7_ER_1","MCF7_ER_2","MCF7_ER_3","MCF7_input"),
                         Group = c("IP", "IP", "IP", "Input"))
sample.meta

## ----load_track_chip----------------------------------------------------------
# track folder
track.folder = system.file("extdata", "ChIP-seq", package = "ggcoverage")
# load bigwig file
track.df = LoadTrackFile(track.folder = track.folder, format = "bw",
                         meta.info = sample.meta)
# check data
head(track.df)

## ----prepare_mark_chip--------------------------------------------------------
# create mark region
mark.region=data.frame(start=c(76822533),
                       end=c(76823743),
                       label=c("Promoter"))
# check data
mark.region

## ----basic_coverage_chip, eval=FALSE------------------------------------------
#  basic.coverage = ggcoverage(data = track.df, color = "auto", region = "chr18:76822285-76900000",
#                              mark.region=mark.region, show.mark.label = FALSE)
#  basic.coverage

## ----basic_coverage_chip_plot, echo=FALSE, fig.height = 6, fig.width = 12, fig.align = "center"----
knitr::include_graphics("../man/figures/README-basic_coverage_chip-1.png")

## ----peak_coverage, eval=FALSE------------------------------------------------
#  # get consensus peak file
#  peak.file = system.file("extdata", "ChIP-seq", "consensus.peak", package = "ggcoverage")
#  
#  basic.coverage +
#    geom_gene(gtf.gr=gtf.gr) +
#    geom_peak(bed.file = peak.file) +
#    geom_ideogram(genome = "hg19",plot.space = 0)

## ----peak_coverage_plot, echo=FALSE, fig.height = 10, fig.width = 12, fig.align = "center"----
knitr::include_graphics("../man/figures/README-peak_coverage-1.png")

## ----session------------------------------------------------------------------
sessionInfo()

