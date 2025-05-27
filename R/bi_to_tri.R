
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
#' @param clust_dist Integer. Distance to base clusters on in clust function if `use_clust` = TRUE.
#' @param overrides Data frame containing subspecies and poly columns indicating any subspecies to
#' override the overlap checking and force binomials to be updated within the polygon type
#' specified in 'poly', where 'poly' is a vector of values containing any of 'dist',
#' 'mcp' or 'clust' for distribution, minimum convex polygon or cluster respectively.
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
                      , use_clust = FALSE
                      , pres_x = "long"
                      , pres_y = "lat"
                      , pres_crs = 4326
                      , use_crs = NULL
                      , buf = 0
                      , clust_dist = 50000
                      , overrides = tibble::tribble(~subspecies, ~ poly,
                                                    "Stipiturus malachurus intermedius", c("dist", "mcp")
                      )
) {

  # presences ----

  all_pres <- presences %>%
    dplyr::distinct(subspecies, !!rlang::ensym(pres_x), !!rlang::ensym(pres_y))

  bi_pres <- all_pres %>%
    dplyr::filter(is.na(subspecies)) |>
    dplyr::select(-subspecies)

  tri_pres <- all_pres %>%
    dplyr::filter(!is.na(subspecies))


  # polygons ----

  ## distributions ----

  if(!is.null(distrib_files)) {

    tri_dist <- distrib_files %>%
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

  ## mcp ----

  if(use_mcp) {

    if(is.null(mcp_files) & nrow(tri_pres)) {

      mcp_prep <- tri_pres %>%
        dplyr::rename(tidyr::any_of(c("long" = "cell_long", "lat" = "cell_lat"))) |> # to overcome odd error of pres_x and pres_y not being accepted in make_mcp below
        dplyr::add_count(subspecies) |>
        dplyr::filter(n >= 3) |>
        dplyr::select(-n) |>
        tidyr::nest(data = -c(subspecies))

      tri_mcp <- purrr::pmap(list(mcp_prep$data
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

      tri_mcp <- mcp_files %>%
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

  ## clusters ----

  if(use_clust) {

    clust_prep <- all_pres %>%
      dplyr::select(-subspecies) %>%
      dplyr::rename(tidyr::any_of(c("long" = "cell_long", "lat" = "cell_lat"))) %>% # to overcome odd error of pres_x and pres_y not being accepted in make_clust below
      sf::st_as_sf(coords = c("long", "lat"), crs = pres_crs, remove = FALSE) %>%
      sf::st_transform(crs = use_crs) %>%
      dplyr::mutate(X = sf::st_coordinates(.)[,1]
                    , Y = sf::st_coordinates(.)[,2]
      ) %>%
      sf::st_set_geometry(NULL) %>%
      dplyr::mutate(clust = clusts(presence = .
                                   , h = 50000
                                   , pres_x = "X"
                                   , pres_y = "Y"
      )
      ) |>
      dplyr::add_count(clust) |>
      dplyr::filter(n >= 3) |>
      tidyr::nest(.by = clust, data = c(long, lat))

    tri_clust <- purrr::pmap(list(clust_prep$data
                                  , clust_prep$clust
    )
    , \(x,y) {

      clust_poly <- make_mcp(presence = x
                             , out_file = tempfile()
                             , force_new = FALSE
                             , pres_x = "long"
                             , pres_y = "lat"
                             , in_crs = pres_crs
                             , out_crs = use_crs
                             , buf = buf
                             , clip = NULL
      ) |>
        dplyr::mutate(clust = y)

    }
    ) |>
      dplyr::bind_rows() |>
      sf::st_make_valid() |>
      sf::st_join(tri_pres %>%
                    sf::st_as_sf(coords = c(pres_x, pres_y), crs = pres_crs) %>%
                    sf::st_transform(., crs = use_crs)
                  , left = FALSE
      ) |>
      dplyr::group_by(subspecies) |>
      dplyr::summarise() |>
      sf::st_make_valid()

  }

  # overlaps ----

  overlap_prep <- tri_pres |>
    dplyr::distinct(subspecies) |>
    dplyr::pull(subspecies) |>
    purrr::set_names() |>
    purrr::map(\(x) {

      tri_pres |>
        dplyr::filter(subspecies != x) |>
        dplyr::select(!subspecies)

    }
    ) |>
    purrr::list_rbind(names_to = "subspecies") |>
    tidyr::nest(.by = subspecies, other_pres = tidyr::all_of(c(pres_x, pres_y))) %>%
    {if(exists("tri_dist")) dplyr::left_join(., tri_dist |>
                                               tidyr::nest(.by = subspecies, dist = tidyr::everything())
    ) else .} %>%
    {if(exists("tri_mcp")) dplyr::left_join(., tri_mcp |>
                                              tidyr::nest(.by = subspecies, mcp = tidyr::everything())
    ) else .} %>%
    {if(exists("tri_clust")) dplyr::left_join(., tri_clust |>
                                                tidyr::nest(.by = subspecies, clust = tidyr::everything())
    ) else .} |>
    tidyr::pivot_longer(cols = c(dist, mcp, clust), values_to = "poly") %>%
    dplyr::mutate(poly = dplyr::select(., name, poly) |> tibble::deframe()) |>
    dplyr::group_by(subspecies, other_pres) |>
    dplyr::summarise(poly = list(poly)) |>
    dplyr::ungroup()

  ## overlap ----

  overlaps <- purrr::pmap(list(overlap_prep$subspecies
                               , overlap_prep$other_pres
                               , overlap_prep$poly
  )
  , \(x,y,z) {

    overlap <- purrr::imap(z, \(p, idp) {

      ### poly ----

      if(!is.null(p)) {

        poly <- p |>
          dplyr::filter(subspecies == x) |>
          dplyr::select(-subspecies)

        other_poly <- p |>
          dplyr::filter(subspecies != x)

        poly_overlap <- other_poly |>
          sf::st_join(poly, left = FALSE)

      } else {

        poly_overlap <- tibble::tibble()

      }

      ### pres ----

      if(!is.null(y) & !is.null(p)) {

        pres_overlap <- y %>%
          sf::st_as_sf(coords = c(pres_x, pres_y), crs = pres_crs) %>%
          {if(!is.null(use_crs)) sf::st_transform(., crs = use_crs)
            else sf::st_transform(crs = sf::st_crs(., z[[1]]))
          } %>%
          sf::st_join(poly %>%
                        {if(!is.null(use_crs)) sf::st_transform(., crs = use_crs) else .}
                      , left = FALSE
          )

      } else {

        pres_overlap <- tibble::tibble()

      }

      res <- tibble::tibble(subspecies = x
                            , poly = idp
                            , pres_overlap = nrow(pres_overlap)
                            , poly_overlap = nrow(poly_overlap)
      ) |>
        dplyr::mutate(overlap = dplyr::if_any(tidyr::contains("overlap"), \(o) o > 0))

      return(res)

    }
    ) |>
      dplyr::bind_rows()

    return(overlap)

  }
  ) |>
    dplyr::bind_rows()

  # update binomials ----

  ## prep ----

  update_prep <- overlaps |>
    dplyr::filter(!overlap) |>
    dplyr::group_by(subspecies) |>
    dplyr::summarise(poly = list(poly)) |>
    dplyr::ungroup() %>%
    dplyr::bind_rows(overrides |>
                       dplyr::filter(subspecies %in% overlaps$subspecies)
    ) |>
    dplyr::rename(use_poly = poly) |>
    dplyr::left_join(overlap_prep) |>
    dplyr::select(-other_pres)

  ## polygons ----

  update_poly <- purrr::pmap(list(update_prep$subspecies
                                  , update_prep$use_poly
                                  , update_prep$poly
  )
  , \(x,y,z) {

    run <- all(purrr::map_lgl(unname(z[y]), \(r) !is.null(r)))

    if(run) {

      tri_poly <- z[y] %>%
        purrr::map(\(p) {

          if(!is.null(p)) {

            p |>
              dplyr::filter(subspecies == x) %>%
              {if(!is.null(use_crs)) sf::st_transform(., crs = use_crs)
                else sf::st_transform(., crs = z[[1]])
              }

          } else {

            dplyr::slice(tri_dist, 0) %>%
              {if(!is.null(use_crs)) sf::st_transform(., crs = use_crs)
                else sf::st_transform(., crs = z[[1]])
              }

          }

        }
        ) |>
        dplyr::bind_rows() |>
        dplyr::summarise() |>
        dplyr::mutate(subspecies = x)

    } else {

      tri_poly <- dplyr::slice(tri_dist, 0) %>%
        {if(!is.null(use_crs)) sf::st_transform(., crs = use_crs)
          else sf::st_transform(., crs = z[[1]])
        }

    }

    return(tri_poly)

  }
  ) |>
    dplyr::bind_rows()

  ## update names ----
  # done separately to above to enable easy mapping of polygons used vs updated occurrences
  update_names <- update_poly |>
    tidyr::nest(.by = subspecies, poly = c(subspecies, geometry)) |>
    dplyr::pull(poly) |>
    purrr::map(\(p) {

      new_names <- bi_pres %>%
        sf::st_as_sf(coords = c(pres_x, pres_y), crs = pres_crs, remove = FALSE) %>%
        {if(!is.null(use_crs)) sf::st_transform(., crs = use_crs)
          else sf::st_transform(crs = sf::st_crs(., p))
        } |>
        sf::st_join(p, left = FALSE) |>
        sf::st_set_geometry(NULL) |>
        dplyr::mutate(bi_to_tri = TRUE)

      return(new_names)

    }

    ) |>
    dplyr::bind_rows()

  ## checks ----
  if(FALSE) {

    library(tmap)
    tmap_mode("view")

    test <- bi_pres %>%
      dplyr::left_join(update_names) |>
      sf::st_as_sf(coords = c(pres_x, pres_y), crs = pres_crs, remove = FALSE)

    tm_shape(test)+tm_dots("subspecies")+tm_shape(update_poly)+tm_borders("subspecies")

    test |>
      sf::st_set_geometry(NULL) |>
      dplyr::count(subspecies) |>
      dplyr::mutate(pc = n/nrow(test)*100)

  }

  return(update_names)

}
