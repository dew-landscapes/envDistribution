#' Identify clusters of occurrence points based on distance.
#'
#' Can be used as a precursor to make_mcp or make_ahull for creating refined distribution polygons
#' for species with geographically separated populations/subspecies.
#'
#' @param presence Cleaned and filtered data frame of presences.
#' @param h Distance used to cluster points, i.e. points within individual clusters
#' will not be further apart than this distance.
#' @param pres_x,pres_y Character. Name of the columns in `presence` that have
#' the x and y coordinates.
#'
#' @return Vector of cluster ID numbers for each row in `presence`,
#' indicating which cluster each row belongs.
#'
#' @export
#'

clusts <- function(presence
                   , h
                   , pres_x = "long"
                   , pres_y = "lat"
                   , ...
){

  prep <- presence %>%
    dplyr::mutate(rownames = dplyr::row_number()) |>
    dplyr::select(rownames, !!rlang::ensym(pres_x), !!rlang::ensym(pres_y))

  if(nrow(prep)>1){

    # Using hierarchical clustering on a distance matrix to create a dendogram
    dend <- hclust(dist(prep), method = "single")

    # Cut tree to obtain clusters within specified distance
    clust <- tibble::tibble(rownames = prep$rownames
                            , clust = cutree(dend, h = h)) # Distance with a 50000 m threshold

  } else { # if only one unique occurrence, then can't do clustering above and only one cluster.

    clust <- prep |>
      dplyr::mutate(rownames = 1
                    , clust = 1
      )

  }

  clust <- clust %>%
    dplyr::arrange(rownames) %>%
    dplyr::pull(clust)

  return(clust)

}
