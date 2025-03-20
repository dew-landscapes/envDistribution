#' Identify clusters of occurrence points based on distance.
#'
#' Can be used as a precursor to make_mcp or make_ahull for creating refined distribution polygons
#' for species with geographically separated populations/subspecies.
#'
#' @param presence Cleaned and filtered dataframe of presences
#' @param h Distance used to cluster points, i.e. points within individual clusters
#' will not be further apart than this distance.
#'
#' @return data frame saved as rds to `out_file` path containing environ and land_pc (percent land overlap ) fields for the taxa.
#'
#' @export
#'

clusts <- function(presence, h, ...){

  if(nrow(presence)>1){

    prep <- presence %>%
      dplyr::mutate(rownames = dplyr::row_number()) |>
      dplyr::select(rownames,X,Y)

    # Using hierarchical clustering on a distance matrix to create a dendogram
    dend <- hclust(dist(prep), method="single")

    # Cut tree to obtain clusters within specified distance
    clust <- data.frame(rownames = prep$rownames
                        , clust = cutree(dend, h = h)) # Distance with a 50000 m threshold

  } else { # if only one species record in landscape (can't do clustering above)

    clust <- data.frame(rownames = 1
                        , clust = 1
    )
  }

  clust <- clust %>%
    dplyr::arrange(rownames) %>%
    dplyr::pull(clust)

}
