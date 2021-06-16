# ProMises model
 
## Installation

You can install the released version of ProMises with:

``` r
devtools::install_github("angeella/ProMises")
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

