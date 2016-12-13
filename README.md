# gcapc: GC effects aware peak caller

### Introduction
ChIP-seq has been widely utilized as the standard technology to detect 
protein binding regions, where peak calling algorithms were developed 
particularly to serve the analysis. Existing peak callers lack of power 
on ranking peaks' significance due to sequencing technology might undergo
sequence context biases, *e.g.* GC bias. **gcapc** is designed to address 
this deficiency by modeling GC effects into peak calling.

### Installation

**gcapc** is an R-package. It can be installed into R by following code
using the GitHub source here.
```s
source("https://bioconductor.org/biocLite.R")
biocLite(c("BiocGenerics","GenomeInfoDb","S4Vectors","IRanges","Biostrings",
	"BSgenome","GenomicRanges","Rsamtools","GenomicAlignments","splines"))

library(devtools)
install_github("tengmx/gcapc")
```

Alternatively, it can be installed through Bioconductor by following code
(available soon).
```s
source("https://bioconductor.org/biocLite.R")
biocLite("gcapc")
```

After installation, the package can be loaded into R.

```s
library(gcapc)
```

### Using gcapc

Details of using this package, please see the 
[vignette](https://github.com/tengmx/gcapc/blob/master/vignettes/gcapc.Rmd).

#### Note:

1. This package is currently only available in the development verion of
Bioconductor, and will be added into release version in April 2017. Thus,
develop version of Bioconductor is required for installation through
Bioconductor repository. 
2. For the purpose of synchronization between GitHub
and Bioconductor, development version of R is required if installing through
GitHub. Or, manually change required version number of R in the DESCRIPTION
file and install package locally if preferring working on release version of
R.