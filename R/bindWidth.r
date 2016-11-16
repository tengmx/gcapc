#' @title ChIP-seq Binding Width And Peak Window Size Estimation
#'
#' @description
#' ChIP-seq experiments usually use crosslinking strategy to capture
#' sequencing fragments. The fragment location is affected by at least but
#' not limited to two factors, protein real binding and crosslinking
#' operation. This function estimate size of binding part in crosslinked
#' DNA-protein complexes, and denoted that as ChIP-seq binding width.Also,
#' the peak detection window half size is estimated based on binding width.
#'
#' @param cov A list object returned by function \code{read5endCoverage}.
#'
#' @param range A non-nagative integer vector with length 2. This vector
#' set the range within which binding width and peak window size are 
#' estimated. Default c(50,500) represents most ChIP-seq experiments.
#'
#' @param step A non-negative integer to set the resolution of binding
#' width estimation within \code{range}. This value will be tuned if
#' \code{auto} is TRUE. Default 50 is based on default value of 
#' \code{range}.
#' 
#' @param odd A logical vector which, when TRUE, only allows return odd
#' number of binding width, which is preferred by the 
#' effective GC content estimation. Default: TRUE.
#'
#' @return A numeric vector with 2 elements:
#' Estimated binding width and half size of peak detection window.
#'
#' @import S4Vectors
#' @import IRanges
#' @import GenomicRanges
#'
#' @export
#' @examples
#' bam <- system.file("extdata", "chipseq.bam", package="gcapc")
#' cov <- read5endCoverage(bam)
#' bindWidth(cov)

bindWidth <- function(cov,range=c(50L,500L),step=50L,odd=TRUE){
    cat("Starting to estimate bdwidth.\n")
    if(!is.list(cov) || length(cov)!=2)
        stop("bdwidth: cov is not a list of 2 elements\n")
    if(sum(names(cov) %in% c("fwd","rev"))!=2)
        stop("bdwidth: names of cov is not correct\n")
    if(length(cov$fwd)!=length(cov$rev))
        stop("bdwidth: chroms differ between fwd and rev strands\n")
    if(step<1) stop("bdwidth: step must be at least 1")
    step <- round(step)
    ## cross correlation estimation
    readsum <- sapply(cov$fwd,sum) + sapply(cov$rev,sum)
    chromlen <- sapply(cov$fwd,length)
    cycle <- 1
    rangein <- range
    w2 <- 0
    shifts <- seq(range[1],range[2],step)
    repeat{
        cat("...... Cycle",cycle,"for bind width estimation\n")
        regionpos <- regionneg <- IRangesList()
        for(i in seq_along(cov$fwd)){
            regionpos[[i]] <- IRanges(start=rep(1,length(shifts)),
                              end=chromlen[i]-shifts)
            regionneg[[i]] <- IRanges(start=shifts+1,
                              end=rep(chromlen[i],length(shifts)))
        }
        viewpos <- Views(cov$fwd,regionpos)
        viewneg <- Views(cov$rev,regionneg)
        cors <- sapply(seq_along(viewpos),function(i)
                       cor(viewpos[[i]],viewneg[[i]]))
        corsall <- colSums(t(cors) * readsum / sum(readsum))
        #if(cycle==1) plot(shifts,corsall)
        if(step>1){
          range <- shifts[which.max(corsall)]+c(-step,step)
          step <- max(round(step/5),1)
          cycle <- cycle + 1
          shifts <- seq(range[1],range[2],step)
        }else if(w2==0){
          w1 <- shifts[which.max(corsall)]
          cc <- max(corsall)
          tmp <- seq(0,rangein[2],5)
          shifts <- tmp[tmp > w1]
          cycle <- 'final'
          w2 <- 1
        }else{
          ccw <- (cc-corsall[length(corsall)])/3+corsall[length(corsall)]
          w2 <- min(max(shifts[corsall>=ccw]),w1*2,rangein[2])
          break
        }
    }
    if(odd && w1%%2==0) w1 <- w1 + 1
    cat("...... Estimated bind width as",w1,"\n")
    cat("...... Estimated peak window half size as",w2,"\n")
    ### refine of peak window half size by region correlation
    cat("Refining peak window half size by region from two strands\n")
    halfbdw <- floor(w1/2)
    seqs <- sapply(cov$fwd,length)
    seqs <- floor(seqs/w1-2)*w1
    starts <- lapply(seqs, function(i) seq(1+w1*2, i, w1))
    ends <- lapply(seqs, function(i) seq(w1*3, i, w1))
    chrs <- rep(names(seqs), times=sapply(starts, length))
    sampidx <- sort(sample.int(length(chrs),ceiling(length(chrs)*0.05)))
    region <- GRanges(chrs[sampidx], IRanges(start=unlist(starts)[sampidx],
                                             end=unlist(ends)[sampidx]))
    pdwhs <- seq(w2,w1,-5)
    corr <- c()
    for(pdwh in pdwhs){
      flank <- pdwh-w1+halfbdw
      regionsp <- resize(split(region,seqnames(region)),pdwh)
      rcfwd <- unlist(viewSums(Views(cov$fwd,
                                   ranges(shift(regionsp,-flank)))))
      rcrev <- unlist(viewSums(Views(cov$rev,
                                   ranges(shift(regionsp,halfbdw)))))
      corr <- c(corr,cor(rcfwd,rcrev))
      cat('.')
    }
    #plot(pdwhs,corr)
    cat('\n')
    w2 <- pdwhs[which.max(corr)]
    cat("...... Refined peak window half size as",w2,"\n")
    c(w1,w2)
}
