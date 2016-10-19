#' @title Reads Coverage Using 5-end Base
#'
#' @description
#' Reads coverage in single base pair resolution using only 5-prime end
#' of BAM file records. Coverages are reported for forward and reverse
#' strands separately. Options for customized filtering of BAM records
#' are provided.
#'
#' @param bam The path to a BAM file, which is sorted and indexed.
#'
#' @param chroms NULL or a vector of chromosome names that compatible with
#' the provided BAM file. Reads coverage will be generated for these
#' chromosomes. Default (NULL) will use all chromosomes in BAM file.
#'
#' @param mapq A non-negative integer specifying the minimum mapping
#' quality to include. BAM records with mapping qualities less
#' than \code{mapq} are discarded.
#'
#' @param duplicate A logical vector which, when FALSE (Default), returns
#' maximum coverage of 1 for every base pair. Reads that start at the same
#' position but on different strands are not treated as duplicates.
#'
#' @param flag A returned object by \code{Rsamtools::scanBamFlag}.
#' Additional options for BAM records filtering.
#'
#' @return A list of two objects by \code{GenomicRanges::coverage}
#' \item{fwd}{Coverage object for forward strand.}
#' \item{rev}{Coverage object for reverse strand.}
#'
#' @import GenomeInfoDb
#' @import S4Vectors
#' @import IRanges
#' @import GenomicRanges
#' @importFrom BiocGenerics strand
#' @importFrom Rsamtools BamFile
#' @importFrom Rsamtools scanBamFlag
#' @importFrom Rsamtools ScanBamParam
#' @importFrom GenomicAlignments readGAlignments
#'
#' @export
#' @examples
#' bam <- system.file("extdata/chipseq.bam",package="gcapc")
#' cov <- rc5end(bam)

rc5end <- function(bam,chroms=NULL,mapq=30L,duplicate=FALSE,
                   flag=scanBamFlag(isUnmappedQuery=FALSE,
                       isSecondaryAlignment=FALSE,
                       isNotPassingQualityControls=FALSE)){
    baminfo <- seqinfo(BamFile(bam))
    if(!is.null(chroms)) {
        if(sum(!(chroms %in% seqlevels(baminfo)))>0)
            stop("chromosome names not compatible with bam files")
        gal <- readGAlignments(bam,param=ScanBamParam(what = "mapq",
                                         flag = flag, mapqFilter = mapq,
                                         which = GRanges(baminfo)[chroms]))
        gr <- granges(gal)
        seqlevels(gr, force=TRUE) <- chroms
    }else{
        gal <- readGAlignments(bam,param=ScanBamParam(what = "mapq",
                                         flag = flag, mapqFilter = mapq))
        gr <- granges(gal)
    }
    grf <- resize(gr[strand(gr)=="+"],1)
    grr <- resize(gr[strand(gr)=="-"],1)
    covf <- coverage(grf)
    covr <- coverage(grr)
    if(!duplicate) {
        runValue(covf)[runValue(covf)>1] <- 1
        runValue(covr)[runValue(covr)>1] <- 1
    }
    list(fwd=covf,rev=covr)
}
