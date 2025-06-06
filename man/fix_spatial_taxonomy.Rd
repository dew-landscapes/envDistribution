% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fix_spatial_taxonomy.R
\name{fix_spatial_taxonomy}
\alias{fix_spatial_taxonomy}
\title{Fix taxonomies based on spatial layers (i.e. in particular geographic areas).}
\usage{
fix_spatial_taxonomy(
  bio_df,
  coords = c("long", "lat"),
  crs = 4326,
  taxa_df,
  taxa_ds,
  change_if_intersects = FALSE,
  taxa_col = "taxa",
  taxonomy,
  levels_to_change = c("species", "subspecies")
)
}
\arguments{
\item{bio_df}{Data frame with taxa column, and lat & long columns.}

\item{coords}{Character vector. Spatial coordinate columns in order of c("x","y"), e.g. c("long", "lat").}

\item{crs}{Integer. Coordinate reference system that the coords correspond to (epsg code).}

\item{taxa_df}{Data frame with taxa1 and taxa2 columns indicating the names of the taxa pairs to be fixed.
Single taxa can be input to fix only their records based on to their distribution alone
by entering the single taxa name into the taxa1 column, and leaving taxa2 as NA.
The taxa names need to correspond to the taxa distributions in \code{taxa_ds}.}

\item{taxa_ds}{Data frame. Taxa distribution sources data frame containing 'taxa' and 'file' columns.
File column should reference the file path of the distribution for the associated taxa in the 'taxa' column corresponding to taxa1 and taxa2.}

\item{change_if_intersects}{Logical. If TRUE, change taxonomy if intersects with layer,
or if FALSE, change if doesn't intersect with layer.}

\item{taxa_col}{Character. Taxa column in bio_df and taxa_ds (must be the same).}

\item{taxonomy}{Taxonomy object returned by envClean::make_taxonomy in relation to the taxa in bio_df and taxa_ds.}

\item{levels_to_change}{Character vector. At what level of the taxonomic hierarchy are the taxonomic fixes (either 'species' or 'subspecies' or both).
Default is c("species","subspecies").}
}
\value{
Original data frame with corrected taxa names for taxa input into taxa1 and taxa2.
}
\description{
If >=2 range layers exist (i.e. one for each taxa), the closest will be used to update the taxonomic name.
If only one layer exists, the taxa1 will be changed to the taxa2 according to the change_if_intersects argument.
}
