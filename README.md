# alignProMises: Alignment by the Procrustes von Mises-Fisher model in R
 
 [![DOI](https://zenodo.org/badge/290147620.svg)](https://zenodo.org/badge/latestdoi/290147620)

 
## Installation

You can install the released version of alignProMises with:

``` r
devtools::install_github("angeella/alignProMises")
``` 
 
## Example 
 
```r
data<- array(rnorm(24576*60*3), dim= c(24576,60,3))
#data <- list(data[,,1],data[,,2],data[,,3])

system.time(out <-ProMisesModel(data, maxIt = 20, t = 1/100, k = 1,
                            scaling = FALSE, reflection = FALSE, 
                            subj = FALSE, centered = TRUE))
out$Xest
plot(out$dist, type = 'l')
```

# References

Andreella, A., & Finos, L. (2022). Procrustes analysis for high-dimensional data. psychometrika, 87(4), 1422-1438.

Andreella, A., Finos, L., & Lindquist, M. A. (2023). Enhanced hyperalignment via spatial prior information. Human Brain Mapping, 44(4), 1725-1740.

Gower, J. C., & Dijksterhuis, G. B. (2004). Procrustes problems (Vol. 30). OUP Oxford.

Haxby, J. V., Guntupalli, J. S., Connolly, A. C., Halchenko, Y. O., Conroy, B. R., Gobbini, M. I., ... & Ramadge, P. J. (2011). A common, high-dimensional model of the representational space in human ventral temporal cortex. Neuron, 72(2), 404-416.

Goodall, C. (1991). Procrustes methods in the statistical analysis of shape. Journal of the Royal Statistical Society: Series B (Methodological), 53(2), 285-321.

# Did you find some bugs?

Please write to angela.andreella[\at]unive[\dot]it or insert a reproducible example using [reprex](https://github.com/tidyverse/reprex) on my [issue github page](https://github.com/angeella/alignProMises/issues).
