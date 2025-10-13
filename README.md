
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
devtools::install_github("Calamanthus/envDistribution")
```

## Contents of `envDistribution`

The following functions and data sets are provided in `envDistribution`.

| object | class | description |
|:---|:---|:---|
| `envDistribution::bi_to_tri()` | function | Convert binomial (species) level occurrences to trinomial based on distributions. |
| `envDistribution::calc_environ()` | function | Calculate a species environment category based on occurrences. |
| `envDistribution::clusts()` | function | Identify clusters of occurrence points based on distance. |
| `envDistribution::dists_source()` | function | Find file path to relevant taxa distribution (geographic range) parquet files. |
| `envDistribution::filter_by_distribution()` | function | Filter species occurrence records by spatial layers representing their distributions at state and or national scales. |
| `envDistribution::fix_spatial_taxonomy()` | function | Fix taxonomies based on spatial layers (i.e. in particular geographic areas). |
| `envDistribution::make_grd()` | function | Make grid around occurrence records |
| `envDistribution::make_mcp()` | function | Make minimum convex polygon (MCP) around occurrence records |
| `envDistribution::reg_cont()` | function | Regional contribution to range (AOO and EOO) |
| `envDistribution::reg_taxa()` | function | Return relevant taxa to a region of interest based on occurrences and or distributions. |
