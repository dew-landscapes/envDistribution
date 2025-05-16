
#' Return relevant taxa to a region of interest based on occurrences and or distributions.
#'
#' Useful for species level analyses seeking to link to relevant trinomials (subspecies, races, varieties, forms, etc),
#' where there may be other irrelevant trinomials occurring outside the region of interest.
#'
#' @param presence Cleaned and filtered dataframe of presences.
#' @param distrib_file Character. File path/s for the relevant distributions per taxa. If no relevant distribution, then use NA.
#' Currently, only parquet files are accepted.
#' Use sfarrow::st_write_parquet to write sf objects to parquet.
#' @param taxa Character. Name of the taxa. Needed for joining to taxonomy to obtain species for trinomials (see below).
#' @param pres_x,pres_y Character. Name of the columns in `presence` that have
#' the x and y coordinates.
#' @param region_bound sf object. Region boundary polygon.
#' @param buf Distance in metres to buffer the `region_bound`.
#' @param in_crs epsg code for coordinates in `presence`.
#' @param use_crs Integer. Coordinate reference system (epsg code) to use as the standard for all objects in spatial joins or intersections.
#' Using a projected coordinate system removes the sf warning around assuming coordinates are planar if they're not, and
#' can sometimes help with spherical geometry issues. If NULL, the crs of the region_bound will be used as the standard
#' (i.e. the presence coords will be converted to this crs before joining or intersecting).
#' @param taxonomy Taxonomy object returned by envClean::make_taxonomy in relation to the taxa.
#' Must have subspecies level taxonomy.
#' @param remove Logical. Remove out of region taxa?
#' @param use_mcp Logical. Use a minimum convex polygon (MCP) around presences to determine region overlap?
#' @param mcp_file Character. If use_mcp == TRUE, optional file path of existing mcp to be used.
#' If use_mcp == TRUE and mcp_file == NULL, a new mcp will be constructed from `presence`.
#' Currently, only geoparquet files are accepted. Use sfarrow to convert sf objects to geoparquet.
#'
#' @return Data frame with species column, subspecies column, and 'pres_dist' & 'distrib_dist' columns indicating
#' distance in metres from region to closest presence and distribution respectively. If remove == TRUE, rows with
#' out of region subspecies will not appear.
#'
#' @export
#'

reg_taxa <- function(presence
                     , distrib_file = NA
                     , region_bound
                     , taxa
                     , pres_x = "long"
                     , pres_y = "lat"
                     , pres_crs = 4326
                     , use_crs = NULL
                     , buf = 0
                     , taxonomy
                     , remove = FALSE
                     , use_mcp = FALSE
                     , mcp_file = NULL
) {

  region_bound <- region_bound %>%
    {if(!is.null(use_crs)) sf::st_transform(., crs = use_crs) else .}

  # Distance of closest presence to region ----
  pres_dist <- presence %>%
    dplyr::distinct(!!rlang::ensym(pres_x), !!rlang::ensym(pres_y)) %>%
    sf::st_as_sf(coords = c(pres_x, pres_y), crs = pres_crs) %>%
    {if(!is.null(use_crs)) sf::st_transform(., crs = use_crs)
      else sf::st_transform(crs = sf::st_crs(., region_bound))
    } %>%
    nngeo::st_nn(region_bound
                 , .
                 , k = 1
                 , maxdist = Inf
                 , returnDist = TRUE
                 , progress = FALSE
                 , parallel = 1
    ) %>%
    .$dist %>%
    unlist() %>%
    tibble::as_tibble_col(column_name = "pres_dist") %>%
    dplyr::mutate(taxa = taxa)

  # Distance of distribution(s) to region ----
  if(!any(is.na(distrib_file)|is.null(distrib_file))){

    distrib_dist <- distrib_file %>%
      purrr::map(\(x) sfarrow::st_read_parquet(x)) %>%
      purrr::list_rbind() %>%
      sf::st_sf() %>%
      sf::st_make_valid() %>%
      {if(!is.null(use_crs)) sf::st_transform(., crs = use_crs)
        else sf::st_transform(crs = sf::st_crs(., region_bound))
      } %>%
      sf::st_make_valid() %>%
      dplyr::summarise() %>%
      sf::st_make_valid() %>%
      nngeo::st_nn(region_bound
                   , .
                   , k = 1
                   , maxdist = Inf
                   , returnDist = TRUE
                   , progress = FALSE
                   , parallel = 1
      ) %>%
      .$dist %>%
      unlist() %>%
      tibble::as_tibble_col(column_name = "distrib_dist") %>%
      dplyr::mutate(taxa = taxa)

  }

  # MCP overlap with region ----
  if(use_mcp){

    if(!is.null(mcp_file)) {

      mcp <- sfarrow::st_read_parquet(mcp_file)

    } else {

      mcp <- envDistribution::make_mcp(presence = presence
                                       , out_file = tempfile()
                                       , force_new = FALSE
                                       , pres_x = pres_x
                                       , pres_y = pres_y
                                       , in_crs = pres_crs
                                       , out_crs = use_crs
                                       , buf = 0
                                       , clip = NULL
      )

    }

    mcp_dist <- mcp |>
      sf::st_transform(crs = use_crs) %>%
      nngeo::st_nn(region_bound
                   , .
                   , k = 1
                   , maxdist = Inf
                   , returnDist = TRUE
                   , progress = FALSE
                   , parallel = 1
      ) %>%
      .$dist %>%
      unlist() %>%
      tibble::as_tibble_col(column_name = "mcp_dist") %>%
      dplyr::mutate(taxa = taxa)

  }

  # Relevant taxa to region ----
  pres_distrib <- pres_dist %>%
    dplyr::mutate(in_region = pres_dist <= buf) %>%
    {if(!any(is.na(distrib_file)|is.null(distrib_file))) dplyr::full_join(., distrib_dist)
      else dplyr::mutate(., in_region = distrib_dist <= buf|in_region)
    } %>%
    {if(use_mcp) dplyr::full_join(., mcp_dist)
      else dplyr::mutate(., in_region = mcp_dist <= buf|in_region)
    } %>%
    dplyr::select(tidyr::any_of(c("taxa", "pres_dist", "distrib_dist", "in_region"))) %>%
    dplyr::distinct() %>%
    {if(remove) dplyr::filter(., in_region) else .}

  rel_taxa <- pres_distrib |>
    dplyr::left_join(taxonomy$subspecies$lutaxa |>
                       dplyr::select(taxa, returned_rank)
    ) |>
    dplyr::left_join(taxonomy$subspecies$taxonomy |>
                       dplyr::select(taxa, species)
    ) |>
    dplyr::mutate(subspecies = ifelse(returned_rank == "subspecies" & stringr::str_count(taxa, "\\w+") > 2, taxa, NA)) |> # count of words needed to overcome erroneous subspecies e.g. arising from hybrids - Acacia provincialis/retinodes
    dplyr::select(tidyr::any_of(c("species", "subspecies", "pres_dist", "distrib_dist", "in_region"))) %>%
    dplyr::distinct() %>%
    {if(remove) dplyr::select(., -in_region) else .}

  return(rel_taxa)

}
