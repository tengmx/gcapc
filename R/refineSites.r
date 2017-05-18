#' @title Adjust ChIP-seq Read Count Table 
#'
#' @description
#' For a given set of sites with the same/comparable width, their
#' read count table from multiple samples are adjusted based on
#' potential GC effects. For each sample separately, GC effects are
#' estimated based on their effective GC content and
#' reads count using generalized linear mixture models. Then, count
#' table is adjusted based on estimated GC effects.
#' It it important that the given sites includes both foreground and
#' background regions, see \code{sites} below.
#'
#' @param counts A count matrix with each row corresponding to each element
#' in \code{sites} and each column corresponding to one sample. Every value
#' in the matrix indicates the read counts for one site in one sample. It is
#' noted that since effective GC content is used in this function, it is
#' important to extend either original reads or original \code{sites} to
#' consider reads that 5' starting in \code{flank} regions, when counting
#' sequencing reads.
#'
#' @param sites A GRanges object with length equivalent to number of rows
#' in \code{counts} matrix. It is preferable that every GRange have the same
#' width; otherwise, the mixture model is modeling different things with
#' wider GRanges certainly have more reads. However, it is OK if only a
#' minority of GRanges have different width, since the model is pretty robust
#' to outliers. Also, it is important that \code{sites} including both
#' foreground and background regions in each sample, otherwise the mixture
#' model will fail to fit two components. Fortunately, if you are inputing
#' a large collection of samples, foreground sites in one sample may play
#' the role as background in other samples. In this case, manually selecting
#' real background is not necessary.
#'
#' @param flank A non-negative integer specifying the flanking width of
#' ChIP-seq binding. This parameter provides the flexibility that reads
#' appear in flankings by decreased probabilities as increased distance
#' from binding region. This paramter helps to define effective GC
#' content calculation.
#'
#' @param outputidx A logical vector with the length equivalent to number
#' of rows in \code{counts}. This provides which subset of adjusted count
#' matrix should be outputed. This would be extremely useful if you have
#' manually collected background sites and want to only export the sites
#' you care about.
#' 
#' @param gcrange A non-nagative numeric vector with length 2. This vector
#' set the range of GC content to filter regions. For human, most regions
#' have GC content between 0.3 and 0.8, which is set as the default. Other
#' regions with GC content beyond this range will be ignored.
#'
#' @param emtrace A logical vector which, when TRUE (default), allows to 
#' print the trace of log likelihood changes in EM iterations.
#'
#' @param plot A logical vector which, when TRUE (default), returns miture
#' fitting plot.
#'
#' @param model A character specifying the distribution model to be used in
#' generalized linear model fitting. The default is negative
#' binomial(\code{nbinom}), while \code{poisson} is also supported currently.
#' More details see \code{gcEffects}.
#' 
#' @param mu0 A non-negative numeric initiating read count signals for
#' background sites. This is treated as the starting value of background mean
#' for poisson/nbinom fitting. 
#'
#' @param mu1 A non-negative numeric initiating read count signals for
#' foreground sites. This is treated as the starting value of foreground mean
#' for poisson/nbinom fitting.
#'
#' @param theta0 A non-negative numeric initiating the shape parameter of
#' negative binomial model for background sites. For more detail, see
#' theta in \code{\link[MASS]{glm.nb}} function.
#'
#' @param theta1 A non-negative numeric initiating the shape parameter of
#' negative binomial model for foreground sites. For more detail, see
#' theta in \code{\link[MASS]{glm.nb}} function.
#' 
#' @param p A non-negative numeric specifying the proportion of foreground
#' sites in all estimated sites. This is treated as a starting value for
#' EM algorithm. 
#'
#' @param converge A non-negative numeric specifying the condition of EM
#' algorithm termination. EM algorithm stops when the ratio of log likelihood
#' increment to whole log likelihood is less or equivalent to 
#' \code{converge}.
#'
#' @param genome A \link[BSgenome]{BSgenome} object containing the sequences
#' of the reference genome that was used to align the reads, or the name of
#' this reference genome specified in a way that is accepted by the
#' \code{\link[BSgenome]{getBSgenome}} function defined in the \pkg{BSgenome}
#' software package. In that case the corresponding BSgenome data package
#' needs to be already installed (see \code{?\link[BSgenome]{getBSgenome}} in
#' the \pkg{BSgenome} package for the details).
#' 
#' @param gctype A character vector specifying choice of method to calculate
#' effective GC content. Default \code{ladder} is based on uniformed fragment
#' distribution. A more smoother method based on tricube assumption is also
#' allowed. However, tricube should be not used if \code{flank} is too large.
#'
#' @return The count matrix after GC adjustment. The matrix values are not
#' integer any more.
#'
#' @import S4Vectors
#' @import IRanges
#' @import GenomicRanges
#' @import Biostrings
#' @import MASS
#' @importFrom BSgenome getBSgenome
#' @importFrom BSgenome getSeq
#' @importFrom matrixStats colMedians
#' @importFrom splines ns
#' @importFrom grDevices rgb
#' @importFrom graphics plot
#' @importFrom graphics axis
#' @importFrom graphics lines
#' @importFrom stats dpois
#' @importFrom stats dnbinom
#' @importFrom stats glm
#' @importFrom stats predict
#'
#' @export


refineSites <- function(counts,sites,flank=250L,
                        outputidx=rep(TRUE,nrow(counts)),
                        gcrange=c(0.3,0.8),emtrace=TRUE,plot=TRUE,
                        model=c('nbinom','poisson'),
                        mu0=1,mu1=50,theta0=mu0,theta1=mu1,
                        p=0.2,converge=1e-4,
                        genome="hg19",gctype=c("ladder","tricube")){
    ### input sanity check
    genome <- getBSgenome(genome)
    model <- match.arg(model)
    gctype <- match.arg(gctype)
    sitew <- median(width(sites))
    if(sum(gcrange<0)>0 || sum(gcrange>1)>0 || sum(is.na(gcrange))>0)
        stop("Parameter 'gcrange' error.\n")
    if(mu0<=0 || mu1<=0 || mu0>=mu1)
        stop("Parameter 'mu0' or 'mu1' error.\n")
    if(model=='nbinom' && (theta1<=0 || theta0<=0))
        stop("Parameter 'theta0' or 'theta1' error in nbinom model.\n")
    if(p<=0 || p>=1)
        stop("'p' must be in (0,1).\n")
    if(converge<=0 || converge>0.1)
        stop("'converge' must be in (0,0.1].\n")
    if(gctype=="ladder"){
        weight <- c(seq_len(flank),rep(flank+1,sitew),rev(seq_len(flank)))
        weight <- weight/sum(weight)
    }else if(gctype=="tricube"){
        w <- flank+floor(sitew/2)
        weight <- (1-abs(seq(-w,w)/w)^3)^3
        weight <- weight/sum(weight)
    }
    ### effective gc content
    cat("Start to estimate GC effects.\n")
    cat("...... Calculating GC content with flanking",flank,"\n")
    nr <- shift(resize(sites,width(sites) + flank*2),-flank)
    seqs <- getSeq(genome,nr)
    gcpos <- startIndex(vmatchPattern("S", seqs, fixed="subject"))
    gc <- round(sapply(gcpos,function(x) sum(weight[x])),3)
    rm(nr,seqs,gcpos)
    ### em algorithms
    cat("...... Estimating GC effects\n")
    fitmu0 <- fitmu1 <- fitz <- matrix(NA,sum(outputidx),ncol(counts))
    for(rep in seq_len(ncol(counts))){
        cat("......... Estimating sample",rep,"\n")
        rc <- counts[,rep]
        idx <- gc>=gcrange[1] & gc<=gcrange[2] & !is.na(rc) & !is.na(gc)
        dat <- data.frame(y=rc[idx],gc=gc[idx])
        theta1E <- theta1
        theta0E <- theta0
        if(model=='poisson'){
            logp1 <- dpois(dat$y, lambda = mu1, log = TRUE)
            logp0 <- dpois(dat$y, lambda = mu0, log = TRUE)
        }else{
            logp1 <- dnbinom(dat$y, size=theta1E, mu=mu1, log = TRUE)
            logp0 <- dnbinom(dat$y, size=theta0E, mu=mu0, log = TRUE)
        }
        z <- 1/(1+exp(logp0-logp1)*(1-p)/p)
        llf <- sum(z*(logp1+log(p)) + (1-z)*(logp0+log(1-p)))
        llgap <- llf
        i <- 0
        while(abs(llgap) > (abs(llf) * converge) && i < 100){
            p <- (2+sum(z))/(2*2+length(z))
            dat1 <- dat[z>=0.5,]
            dat0 <- dat[z<0.5,]
            if(model=='poisson'){
                lmns0 <- glm(y ~ ns(gc, df = 2), data=dat0, family="poisson")
                lmns1 <- glm(y ~ ns(gc, df = 2), data=dat1, family="poisson")
                predY0 <- predict(lmns0,data.frame(gc=dat$gc),type="response")
                predY1 <- predict(lmns1,data.frame(gc=dat$gc),type="response")
                logp1 <- dpois(dat$y, lambda = predY1, log = TRUE)
                logp0 <- dpois(dat$y, lambda = predY0, log = TRUE)
            }else{
                lmns0 <- glm.nb(y ~ ns(gc, df=2),data=dat0,init.theta=theta0E)
                lmns1 <- glm.nb(y ~ ns(gc, df=2),data=dat1,init.theta=theta1E)
                predY0 <- predict(lmns0,data.frame(gc=dat$gc),type="response")
                predY1 <- predict(lmns1,data.frame(gc=dat$gc),type="response")
                theta1E <- lmns1$theta
                theta0E <- lmns0$theta
                logp1 <- dnbinom(dat$y, size=theta1E, mu=predY1, log = TRUE)
                logp0 <- dnbinom(dat$y, size=theta0E, mu=predY0, log = TRUE)
            }
            z <- 1/(1+exp(logp0-logp1)*(1-p)/p)
            if(sum(z>=0.5) < length(gc)*0.0005 | sum(z<0.5) < length(gc)*0.0005)
                break;
            lli <- sum(z*(logp1+log(p)) + (1-z)*(logp0+log(1-p)))
            llgap <- lli - llf
            llf <- lli
            i <- i + 1
            if(emtrace)
                cat("......... Iteration",i,'\tll',llf,'\tincrement',llgap,'\n')
        }
        if(plot){
            idx0 <- sample.int(nrow(dat),min(50000,nrow(dat)))
            plot(dat$gc[idx0],dat$y[idx0]+0.5,col=rgb(0,0,0,alpha=0.1),
                 xlim=gcrange,pch=20,main=paste("rep",rep),log='y',yaxt='n',
                 xlab='Effective GC content',ylab="Read counts")
            idx00 <- sample.int(nrow(dat),min(1000,nrow(dat)))
            idx00 <- idx00[order(dat$gc[idx00])]
            lines(dat$gc[idx00],predY1[idx00]+0.5,col='red',lwd=3)
            lines(dat$gc[idx00],predY0[idx00]+0.5,col='blue',lwd=3)
            axis(side=2, at=c(0,2^(0:10))+0.5, labels=c(0,2^(0:10)))
        }
        ### gc effects
        fitmu0[idx[outputidx],rep] <- predY0[outputidx[idx]]
        fitmu1[idx[outputidx],rep] <- predY1[outputidx[idx]]
        fitz[idx[outputidx],rep] <- z[outputidx[idx]]
        cat("......... Sample",rep,"finished\n")
    }
    gce <- log2(t(t(fitmu1)/colMedians(fitmu1[fitz>=0.5],na.rm=T))) * fitz +
           log2(t(t(fitmu0)/colMedians(fitmu0[fitz>=0.5],na.rm=T))) * (1-fitz)
    counts[outputidx,] / 2^gce
}
