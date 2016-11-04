#' @title ChIP-seq Binding Width Estimation
#'
#' @description
#' ChIP-seq experiments usually use crosslinking strategy to capture
#' sequencing fragments. The fragment location is affected by at least but
#' not limited to two factors, protein real binding and crosslinking
#' operation. This function estimate size of binding part in crosslinked
#' DNA-protein complexes, and denoted that as ChIP-seq binding width.
#'
#' @param cov A list object returned by function \code{read5endCoverage}.
#'
#' @param range A non-nagative integer vector with length 2. This vector
#' set the range within which binding width are estimated.
#' Default c(50,600) represents most ChIP-seq experiments.
#'
#' @param step A non-negative integer to set the resolution of binding
#' width estimation within \code{range}. This value will be tuned if
#' \code{auto} is TRUE. Default 50 is based on default value of \code{range}.
#'
#' @param auto A logical vector which, when TRUE, allows to automatically
#' tune to higher resolution of binding width estimation. The highest
#' resolution of 2bp is allowed. Default: TRUE.
#' 
#' @param odd A logical vector which, when TRUE, only allows return odd
#' number of binding width, which is preferred by the 
#' weighted GC content estimation. Default: TRUE.
#'
#' @return A numeric indicating estimated binding width.
#'
#' @import S4Vectors
#' @import IRanges
#'
#' @export
#' @examples
#' bam <- system.file("extdata", "chipseq.bam", package="gcapc")
#' cov <- read5endCoverage(bam)
#' bindWidth(cov)

bindWidth <- function(cov,range=c(50L,600L),step=50L,auto=TRUE,odd=TRUE){
    cat("Starting to estimate bdwidth.\n")
    if(!is.list(cov) || length(cov)!=2)
        stop("bdwidth: cov is not a list of 2 elements\n")
    if(sum(names(cov) %in% c("fwd","rev"))!=2)
        stop("bdwidth: names of cov is not correct\n")
    if(length(cov$fwd)!=length(cov$rev))
        stop("bdwidth: chroms differ between fwd and rev strands\n")
    readsum <- sapply(cov$fwd,sum) + sapply(cov$rev,sum)
    chromlen <- sapply(cov$fwd,length)
    cycle <- 1
    repeat{
        cat("...... cycle",cycle,"for bind width estimation\n")
        shifts <- seq(range[1],range[2],step)
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
        if(!auto | step==1){
            break
        }
        range <- shifts[which.max(corsall)]+c(-step,step)
        step <- max(round(step/5),1)
        cycle <- cycle + 1
    }
    w <- shifts[which.max(corsall)]
    if(odd && w%%2==0) w <- w + 1
    cat("...... estimated bind width as",w,"\n")
    w
}
