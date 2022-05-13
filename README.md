
<!-- README.md is generated from README.Rmd. Please edit that file -->

# geohclust

<!-- badges: start -->
<!-- badges: end -->

geohclust offers two functions `?geohclust_poly` and `?geohclust_graph`
that enable the clustering of spatial data such as polygons with a
hclust type approach but taking advantages of contiguity constraints.
The contiguity naturally create a sparsely connected graph that can be
leveraged to speed-up the calculations and deal with more than 30000
polygons in seconds.

## Installation

You can install the development version of geohclust from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("comeetie/geohclust")
```

## Example

This is a basic example, we first prepare some spatial polygons data,
here the results at the municipality level in one french department for
the :

``` r
library(geohclust)
library(dplyr)
library(sf)
data("modesshare")
deplist = c(37)
dep = modesshare |> 
  filter(DEP %in% deplist) 
```

Do the clustering and use the classical function from `?hclust`
(`?plot.hclust` and `?cutree`):

``` r
hc=geohclust_poly(dep)
#> Warning: Some features were not numeric and have been removed from the
#> clustering.
plot(hc)
cutree(hc,k=30) |> head(20)
#>  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 
#>  1  2  2  3  4  1  5  6  6  2  6  7  6  1  2  8  1  9  1  2
```

<img src="man/figures/README-clustering-1.png" width="100%" />

You may also use the `?geocutree` function which build directly a
spatial data.frame with the clustering results:

``` r
plot(geocutree(hc,k=30))
```

<img src="man/figures/README-plot-1.png" width="100%" />
