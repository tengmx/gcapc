#' @title CATplot of Peaks
#'
#' @description
#' Plot the consistancy between two peak lists by their significance.
#'
#' @param x A GRanges of identified peaks from one method or one replicate.
#' At least one meta column should be included to allow for significance
#' ranking of peaks.
#'
#' @param y A GRanges of identified peaks from compared method or
#' anoter replicate. At least one meta column should be included to
#' allow for significance ranking of peaks.
#'
#' @param ranks A non-negative integer vector specifying the ranks to
#' be used for CAT plot.
#'
#' @param exclude A GRanges object specifying regions to be excluded for
#' CAT plot, such as the blacklist regions proposed by ENCODE Consortium.
#'
#' @param seqinfo A vector of chromosome names to limit the CAT plot to
#' selected chromosomes. Chromosome names here must be in the same format
#' as \code{seqnames} in \code{x} and \code{y}. This parameter also helps
#' exclude uncommon chromosomes, e.g. using value
#' \code{paste0('chr',c(1:22,'X','Y'))} for human genome. Default: NULL
#' means no limit to chromosomes.
#'
#' @param esx A non-negative integer specifying which meta column of
#' \code{x} to be used to rank peak significance. Larger values in this
#' column should indicate higher significance.
#'
#' @param esy A non-negative integer specifying which meta column of
#' \code{y} to be used to rank peak significance. Larger values in this
#' column should indicate higher significance.
#'
#' @param add A logical vector which, when TRUE, adds the current plotting
#' line to existing plots. FALSE will generate a new plot.
#'
#' @param ... Other parameters passed to \code{plot} or \code{lines}.
#'
#' @return A CAT plot.
#'
#' @import S4Vectors
#' @import GenomicRanges
#' @importFrom graphics plot
#' @importFrom graphics lines
#'
#' @export
#' @examples
#' bam <- system.file("extdata", "chipseq.bam", package="gcapc")
#' cov <- read5endCoverage(bam)
#' bdw <- bindWidth(cov)
#' gcb1 <- gcEffects(cov, bdw, samp=0.15, plot=FALSE)
#' peaks1 <- gcapcPeaks(cov, gcb1, bdw)
#' gcb2 <- gcEffects(cov, bdw, samp=0.2, plot=FALSE)
#' peaks2 <- gcapcPeaks(cov, gcb2, bdw)
#' peaksCAT(peaks1, peaks2, ranks=seq(100,200,5), ylim=c(0,1))

peaksCAT <- function(x,y,ranks=seq(200,min(length(x),length(y),20000),50),
                     exclude=GRanges(),seqinfo=NULL,
                     esx=1,esy=1,add=FALSE,...){
    if(!is.null(seqinfo)) {
        seqlevels(x,force=TRUE) <- seqinfo
        seqlevels(y,force=TRUE) <- seqinfo
        ranks <- ranks[ranks <= min(length(x),length(y))]
    }
    if(length(exclude)>=1){
        fox <- findOverlaps(x,exclude)
        foy <- findOverlaps(y,exclude)
        x <- x[setdiff(seq_along(x),queryHits(fox))]
        y <- y[setdiff(seq_along(y),queryHits(foy))]
        ranks <- ranks[ranks <= min(length(x),length(y))]
    }
    fo <- findOverlaps(x,y)
    xfo <- queryHits(fo)
    yfo <- subjectHits(fo)
    xo <- order(mcols(x)[,esx],decreasing=TRUE)
    yo <- order(mcols(y)[,esy],decreasing=TRUE)
    xycat <- sapply(ranks,function(t){
        idx <- seq_len(t)
        sum(xfo %in% xo[idx] & yfo %in% yo[idx])/t
    })
    if(add) lines(ranks,xycat,...)
    else plot(ranks,xycat,type='l',...)
}
