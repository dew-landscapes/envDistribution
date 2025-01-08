
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `envDistribution`: an R package for working with species distributions

<!-- badges: start -->
<!-- badges: end -->

The goal of `envDistribution` is to help create, source and work with
species distributions (i.e. geographic ranges), where the ultimate goal
is to use distributions to clean species occurrence data for use in
Species Distribution Models (SDMs) or IUCN Red List threatened status
type assessments.

## Installation

`envDistribution` is not on [CRAN](https://CRAN.R-project.org).

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Acanthiza/envDistribution")
```

## Contents of `envDistribution`

The following functions and data sets are provided in `envDistribution`.

| object                                      | class    | description                                                                                                           |
|:--------------------------------------------|:---------|:----------------------------------------------------------------------------------------------------------------------|
| `envDistribution::dists_source()`           | function | Find file path to relevant taxa distribution (geographic range) parquet files.                                        |
| `envDistribution::filter_by_distribution()` | function | Filter species occurrence records by spatial layers representing their distributions at state and or national scales. |
| `envDistribution::fix_spatial_taxonomy()`   | function | Fix taxonomies based spatial layers (i.e. in particular geographic areas).                                            |
| `envDistribution::make_grd()`               | function | Grid around records                                                                                                   |
| `envDistribution::make_mcp()`               | function | Minimum convex polygon around records                                                                                 |
| `envDistribution::reg_cont()`               | function | Regional contribution to national range (AOO and EOO)                                                                 |
