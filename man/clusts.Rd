% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusts.R
\name{clusts}
\alias{clusts}
\title{Identify clusters of occurrence points based on distance.}
\usage{
clusts(presence, h, pres_x = "long", pres_y = "lat", ...)
}
\arguments{
\item{presence}{Cleaned and filtered data frame of presences.}

\item{h}{Distance used to cluster points, i.e. points within individual clusters
will not be further apart than this distance.}

\item{pres_x, pres_y}{Character. Name of the columns in \code{presence} that have
the x and y coordinates.}
}
\value{
Vector of cluster ID numbers for each row in \code{presence},
indicating which cluster each row belongs.
}
\description{
Can be used as a precursor to make_mcp or make_ahull for creating refined distribution polygons
for species with geographically separated populations/subspecies.
}
