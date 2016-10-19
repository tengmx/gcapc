#' @title CATplot of Peaks
#'
#' @description
#' Plot the consistancy between two peak lists.
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
#' bam <- system.file("extdata/chipseq.bam",package="gcapc")
#' cov <- rc5end(bam)
#' bdw <- bdwidth(cov)
#' gcb1 <- gcbias(cov,bdw,samp=0.15,plot=FALSE)
#' peaks1 <- callpeaks(cov,gcb1,bdw)
#' gcb2 <- gcbias(cov,bdw,samp=0.1,plot=FALSE)
#' peaks2 <- callpeaks(cov,gcb2,bdw)
#' peakcat(peaks1,peaks2,ranks=seq(100,200,5),ylim=c(0,1))

peakcat <- function(x,y,ranks=seq(200,20000,50),esx=1,esy=1,add=FALSE,...){
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
