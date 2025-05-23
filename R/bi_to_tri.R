
#' Convert binomial (species) level occurrences to trinomial based on distributions.
#'
#' @param presences Data frame of presences relevant to the species and subspecies.
#' Must contain columns named 'species', 'subspecies', and x, y coordinate columns corresponding to the names
#' in `pres_x` & `pres_y`.
#' @param distrib_files Data frame of relevant distribution file paths for all the subspecies.
#' If no relevant distribution, then use NA. Currently, only geoparquet files are accepted.
#' Use sfarrow::st_write_parquet to write sf objects to parquet.
#' @param use_mcp Logical. Use a minimum convex polygon (MCP) around presences to determine region overlap?
#' @param mcp_files If use_mcp == TRUE, optional data frame of file paths of existing mcp files to be used.
#' If use_mcp == TRUE and mcp_file == NULL, a new mcp will be constructed from `presences`.
#' Currently, only geoparquet files are accepted. Use sfarrow to convert sf objects to geoparquet.
#' @param pres_x,pres_y Character. Name of the columns in `presences` that have
#' the x and y coordinates.
#' @param in_crs Integer. Coordinate reference system (epsg code) for coordinates in `presences`.
#' @param use_crs Integer. Coordinate reference system (epsg code) to use as the standard for all objects in spatial joins or intersections.
#' Using a projected coordinate system removes the sf warning around assuming coordinates are planar if they're not, and
#' can sometimes help with spherical geometry issues. If NULL, the crs of the distribution will be used as the standard
#' (i.e. the presence coords will be converted to this crs before joining or intersecting).
#' @param buf Integer. Distance in metres to buffer the distribution.
#'
#' @details
#' Will only update binomial occurrence records that fall within a distribution polygon to trinomial if no other
#' trinomial presences or distributions for a species are found intersecting the distribution polygon.
#' This avoids erroneously updating binomials to the wrong trinomial where distributions of subspecies overlap.
#'
#' @return Data frame with updated names in taxa_col for taxa where distributions were found and
#' did not overlap with multiple trinomials in presence or distributions.
#'
#' @export
#'

bi_to_tri <- function(presences
                      , distrib_files = NULL
                      , use_mcp = FALSE
                      , mcp_files = NULL
                      , pres_x = "long"
                      , pres_y = "lat"
                      , pres_crs = 4326
                      , use_crs = NULL
                      , buf = 0
) {

  # presences prep ----

  all_pres <- presences %>%
    dplyr::distinct(subspecies, !!rlang::ensym(pres_x), !!rlang::ensym(pres_y))

  # all_pres |>
  #   dplyr::count(subspecies)

  bi_pres <- all_pres %>%
    dplyr::filter(is.na(subspecies)) |>
    dplyr::select(-subspecies)

  tri_pres <- all_pres %>%
    dplyr::filter(!is.na(subspecies))

  # distributions prep ----

  if(!is.null(distrib_files)) {

    tri_dists <- distrib_files %>%
      dplyr::distinct(subspecies) |>
      dplyr::pull(subspecies) |>
      purrr::map(\(x) {

        dist <- distrib_files |>
          dplyr::filter(subspecies == x) |>
          dplyr::pull(file) |>
          purrr::map(\(f) sfarrow::st_read_parquet(f)) |>
          dplyr::bind_rows() |>
          dplyr::summarise() |>
          sf::st_make_valid() |>
          dplyr::mutate(subspecies = x)

      }
      ) |>
      dplyr::bind_rows() |>
      sf::st_make_valid()

  }

  # mcp prep ----

  if(use_mcp) {

    if(is.null(mcp_files) & nrow(tri_pres)) {

      mcp_prep <- tri_pres %>%
        dplyr::rename(tidyr::any_of(c("long" = "cell_long", "lat" = "cell_lat"))) |> # to overcome odd error of pres_x and pres_y not being accepted in make_mcp below
        dplyr::add_count(subspecies) |>
        dplyr::filter(n >= 3) |>
        dplyr::select(-n) |>
        tidyr::nest(data = -c(subspecies))

      tri_mcps <- purrr::pmap(list(mcp_prep$data
                                   , mcp_prep$subspecies
      )
      , \(x,y) {

        mcp <- make_mcp(presence = x
                        , out_file = tempfile()
                        , force_new = FALSE
                        , pres_x = "long"
                        , pres_y = "lat"
                        , in_crs = pres_crs
                        , out_crs = use_crs
                        , buf = buf
                        , clip = NULL
        ) |>
          dplyr::mutate(subspecies = y)

      }
      ) |>
        dplyr::bind_rows() |>
        sf::st_make_valid()

    } else if(!is.null(mcp_files)) {

      tri_mcps <- mcp_files %>%
        dplyr::distinct(subspecies) |>
        dplyr::pull(subspecies) |>
        purrr::map(\(x) {

          mcp <- mcp_files |>
            dplyr::filter(subspecies == x) |>
            dplyr::pull(file) |>
            purrr::map(\(f) sfarrow::st_read_parquet(f)) |>
            dplyr::bind_rows() |>
            dplyr::summarise() |>
            sf::st_make_valid() |>
            dplyr::mutate(subspecies = x)

        }
        ) |>
        dplyr::bind_rows() |>
        sf::st_make_valid()

    } else {

      use_mcp <- FALSE

    }

  }

  # trinomials per distribution ----

  if(!is.null(distrib_files)) {

    dist_tri <- tri_dists |>
      dplyr::pull(subspecies) |>
      purrr::map(\(x) {

        ## distributions overlap ----
        dist <- tri_dists |>
          dplyr::filter(subspecies == x) |>
          dplyr::select(-subspecies)

        other_dists <- tri_dists |>
          dplyr::filter(subspecies != x)

        dist_overlap <- other_dists |>
          sf::st_join(dist, left = FALSE)

        ## pres overlap ----
        other_pres <- tri_pres |>
          dplyr::filter(subspecies != x) %>%
          sf::st_as_sf(coords = c(pres_x, pres_y), crs = pres_crs) %>%
          {if(!is.null(use_crs)) sf::st_transform(., crs = use_crs)
            else sf::st_transform(crs = sf::st_crs(., dist))
          }

        pres_overlap <- other_pres |>
          sf::st_join(dist %>%
                        {if(!is.null(use_crs)) sf::st_transform(., crs = use_crs) else .}
                      , left = FALSE
          )

        pres_overlap <- other_pres |>
          sf::st_join(dist %>%
                        {if(!is.null(use_crs)) sf::st_transform(., crs = use_crs) else .}
                      , left = FALSE
          )

        res <- tibble::tibble(subspecies = x
                              , pres_overlap = nrow(pres_overlap)
                              , dist_overlap = nrow(dist_overlap)
        ) |>
          dplyr::mutate(overlap = dplyr::if_any(tidyr::contains("overlap"), \(x) x > 0))

        return(res)

      }
      ) |>
      dplyr::bind_rows()

  }

  # trinomials per mcp ----

  if(use_mcp) {

    mcp_tri <- tri_mcps |>
      dplyr::pull(subspecies) |>
      purrr::map(\(x) {

        ## mcp overlap ----
        mcp <- tri_mcps |>
          dplyr::filter(subspecies == x) |>
          dplyr::select(-subspecies)

        other_mcps <- tri_mcps |>
          dplyr::filter(subspecies != x)

        mcp_overlap <- other_mcps |>
          sf::st_join(mcp, left = FALSE)

        ## pres overlap ----
        other_pres <- tri_pres |>
          dplyr::filter(subspecies != x) %>%
          sf::st_as_sf(coords = c(pres_x, pres_y), crs = pres_crs) %>%
          {if(!is.null(use_crs)) sf::st_transform(., crs = use_crs)
            else sf::st_transform(crs = sf::st_crs(., mcp))
          }

        pres_overlap <- other_pres |>
          sf::st_join(mcp %>%
                        {if(!is.null(use_crs)) sf::st_transform(., crs = use_crs) else .}
                      , left = FALSE
          )

        ## distributions overlap ----
        other_dists <- tri_dists |>
          dplyr::filter(subspecies != x)

        dist_overlap <- other_dists %>%
          {if(!is.null(use_crs)) sf::st_transform(., crs = use_crs)
            else sf::st_transform(crs = sf::st_crs(., mcp))
          } |>
          sf::st_join(mcp, left = FALSE)

        res <- tibble::tibble(subspecies = x
                              , pres_overlap = nrow(pres_overlap)
                              , mcp_overlap = nrow(mcp_overlap)
                              , dist_overlap = nrow(dist_overlap)
        ) |>
          dplyr::mutate(overlap = dplyr::if_any(tidyr::contains("overlap"), \(x) x > 0))

        return(res)

      }
      ) |>
      dplyr::bind_rows()

  }

  # update binomials ----

  ## list of trinomials ----
  # trinomials can have +/- dist or mcp, so need to find the right combo that incorporates them all

  tri <- tri_pres %>%
    {if(!is.null(distrib_files)) dplyr::bind_rows(., dist_tri) else .} %>%
    {if(use_mcp) dplyr::bind_rows(., mcp_tri) else .} |>
    dplyr::distinct(subspecies) |>
    dplyr::pull(subspecies)

  ## update per trinomial ----

  res <- tri |>
    purrr::map(\(x) {

      ### choose polygon ----

      distrib <- dist_tri |>
        dplyr::filter(subspecies == x) |>
        dplyr::pull(overlap) |>
        any() |>
        isFALSE()

      mcp <- mcp_tri |>
        dplyr::filter(subspecies == x) %>%
        {if(nrow(.)) dplyr::pull(., overlap) %>%
            any() %>%
            isFALSE()
          else FALSE
        }

      if(distrib & mcp) {

        tri_poly <- tri_dists |>
          dplyr::filter(subspecies == x) %>%
          {if(!is.null(use_crs)) sf::st_transform(., crs = use_crs) else .} %>%
          dplyr::bind_rows(tri_mcps |>
                             dplyr::filter(subspecies == x) %>%
                             {if(!is.null(use_crs)) sf::st_transform(., crs = use_crs)
                               else sf::st_transform(crs = sf::st_crs(., tri_dists))
                             }
          ) |>
          dplyr::summarise() |>
          sf::st_make_valid() |>
          dplyr::mutate(subspecies = x)

      } else if(distrib & !mcp) {

        tri_poly <- tri_dists |>
          dplyr::filter(subspecies == x) %>%
          {if(!is.null(use_crs)) sf::st_transform(., crs = use_crs) else .}

      } else if(!distrib & mcp) {

        tri_poly <- tri_mcps |>
          dplyr::filter(subspecies == x) %>%
          {if(!is.null(use_crs)) sf::st_transform(., crs = use_crs) else .}

      }

      ### update names ----

      new_names <- bi_pres %>%
        sf::st_as_sf(coords = c(pres_x, pres_y), crs = pres_crs, remove = FALSE) %>%
        {if(!is.null(use_crs)) sf::st_transform(., crs = use_crs)
          else sf::st_transform(crs = sf::st_crs(., tri_poly))
        } |>
        sf::st_join(tri_poly, left = FALSE) |>
        sf::st_set_geometry(NULL) |>
        dplyr::mutate(bi_to_tri = TRUE)

    }
    ) |>
    dplyr::bind_rows()

  if(FALSE) {

    library(tmap)
    tmap_mode("view")

    test <- bi_pres %>%
      dplyr::left_join(res) |>
      sf::st_as_sf(coords = c(pres_x, pres_y), crs = pres_crs, remove = FALSE)

    tm_shape(test)+tm_dots("subspecies")+tm_shape(tri_dists)+tm_borders("subspecies")+tm_shape(tri_mcps)+tm_borders("subspecies")

    dplyr::count(test |> sf::st_set_geometry(NULL), subspecies) |>
      dplyr::mutate(pc = n/nrow(test)*100)

  }

  return(res)

}
