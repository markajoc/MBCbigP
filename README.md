## MBCbigP: Model-based Clustering for High-dimensional Data  

[![BuildStatus](https://travis-ci.org/markajoc/condvis.svg?branch=devel)](https://travis-ci.org/markajoc/MBCbigP)
[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)

### Fit a mixture of Gaussians using sequential conditioning

Installation:
```r
devtools::install_github("markajoc/MBCbigP")
```

Example to get started:  
```r
library(MBCbigP)
data(banknote, package = "mclust")
o <- mbcbigp(banknote[, -1], batches = 2)
```
