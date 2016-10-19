# gcapc: GC effects aware peak caller

### Introduction
ChIP-seq has been widely utilized as the standard technology to detect 
protein binding regions, where peak calling algorithms were developed 
particularly to serve the analysis. Existing peak callers lack of power 
on ranking peaks' significance due to sequencing technology might undergo
sequence context biases, *e.g.* GC bias. **gcapc** is designed to address 
this deficiency by modeling GC effects into peak calling.

### Installation

R-package **gcapc** can be installed:
```s
library(devtools)
install_github("tengmx/gcapc")
```
After installation, the package can be loaded into R.

```s
library(gcapc)
```

### Using gcapc

Details of using this package, please see the 
[vignette](https://github.com/tengmx/gcapc/blob/master/vignettes/gcapc.pdf).