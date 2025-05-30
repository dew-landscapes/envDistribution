
#' Convert binomial (species) level occurrences to trinomial based on distributions.
#'
#' @param species Character. Name of the species (binomial) to update occurrences to trinomial for.
#' @param presences Data frame of all binomial and trinomial presences relevant to the species.
#' Must contain a 'subspecies' column and x, y coordinate columns corresponding to the names
#' in `pres_x` & `pres_y`. For binomial presences 'subspecies' should be NA.
#' @param distrib_files Data frame of relevant distribution file paths for all the subspecies
#' (e.g. for distributions sourced from redlist, epbc etc). If no relevant distribution, then use NA.
#' Currently, only geoparquet files are accepted. Use sfarrow::st_write_parquet to write sf objects to parquet.
#' @param use_mcp Logical. Use a minimum convex polygon (MCP) around presences to determine relevant trinomials?
#' @param mcp_files If use_mcp == TRUE, optional data frame of file paths of existing mcp files to be used.
#' If use_mcp == TRUE and mcp_file == NULL, a new mcp will be constructed from `presences`.
#' Currently, only geoparquet files are accepted. Use sfarrow to convert sf objects to geoparquet.
#' @param use_clust Logical. Use clusters around around presences to determine relevant trinomials?
#' If TRUE, `use_crs` must correspond to a projected crs to allow clustering on a positive set of numbers.
#' @param pres_x,pres_y Character. Name of the columns in `presences` that have
#' the x and y coordinates.
#' @param in_crs Integer. Coordinate reference system (epsg code) for coordinates in `presences`.
#' @param use_crs Integer. Coordinate reference system (epsg code) to use as the standard for all objects in spatial joins.
#' Using a projected coordinate system removes the sf warning around assuming coordinates are planar if they're not, and
#' can sometimes help with spherical geometry issues. If NULL, the crs of the distribution will be used as the standard
#' (i.e. the `presences` will be converted to this crs before joining or intersecting).
#' @param buf Integer. Distance in metres to buffer the distribution.
#' @param clust_dist Integer. Distance to base clusters on in clust function if `use_clust` = TRUE.
#' Needs to be in the units of the projected crs supplied in `use_crs`, e.g. metres.
#' @param overrides Character vector containing any trinomials that should override the overlap checking
#' and force binomials within their distribution polygons to be updated to that trinomial even if the
#' binomials overlap distribution polygons of other trinomials.
#' @param skip_bi Character vector of any binomials (species) to be skipped and not have binomials
#' updated to trinomial.
#' @param skip_tri Character vector of any trinomials (subspecies, varieties, forms etc) to be skipped
#' and not have binomials updated to these trinomials.
#'
#' @details
#' Binomial occurrence records won't be updated to trinomial if they fall within multiple trinomial distributions.
#' This avoids updating binomials to the wrong trinomial where distributions of subspecies overlap.
#' Distribution refers generally to all the polygons used to represent distributions,
#' i.e. distribution files (e.g. from redlist or epbc), MCPs, or cluster polygons.
#'
#' @return Data frame equivalent to `presences` with 'subspecies' updated for taxa where
#' distributions were found and did not overlap with multiple trinomials distributions,
#' and an additional 'bi_to_tri' column indicating which records were updated.
#'
#' @export
#'

bi_to_tri <- function(species
                      , presences
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
                      , overrides = c("Stipiturus malachurus intermedius")
                      , skip_bi = NULL
                      , skip_tri = NULL
) {

  # presences ----

  all_pres <- presences %>%
    dplyr::distinct(subspecies, !!rlang::ensym(pres_x), !!rlang::ensym(pres_y))

  bi_pres <- all_pres %>%
    dplyr::filter(is.na(subspecies)) |>
    dplyr::select(-subspecies)

  tri_pres <- all_pres %>%
    dplyr::filter(!is.na(subspecies))

  if(nrow(bi_pres) & any(nrow(tri_pres), !is.null(distrib_files)) & !species %in% skip_bi) {

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
            dplyr::mutate(subspecies = x
                          , poly = "dist"
            ) |>
            sf::st_buffer(buf)

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

        if(nrow(mcp_prep)) {

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
              dplyr::mutate(subspecies = y
                            , poly = "mcp"
              ) |>
              sf::st_buffer(buf)

          }
          ) |>
            dplyr::bind_rows() |>
            sf::st_make_valid()

        }

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
              dplyr::mutate(subspecies = x
                            , poly = "mcp"
              ) |>
              sf::st_buffer(buf)

          }
          ) |>
          dplyr::bind_rows() |>
          sf::st_make_valid()

      }

    }

    ## clusters ----

    if(all(use_clust, nrow(tri_pres), nrow(bi_pres))) {

      # Note: error in `hclust()`: ! size cannot be NA nor exceed 65536
      # may need to have a catch for this that tries a larger h value if it errors or just skip

      clust_prep <- tryCatch(
        {
          all_pres %>%
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

        },
        error = function(e) {

          tibble::tibble()

        }
      )

      if(nrow(clust_prep)) {

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
            dplyr::mutate(clust = y) |>
            sf::st_buffer(buf)

        }
        ) |>
          dplyr::bind_rows() |>
          sf::st_make_valid() |>
          sf::st_join(tri_pres %>%
                        sf::st_as_sf(coords = c(pres_x, pres_y), crs = pres_crs) %>%
                        sf::st_transform(., crs = use_crs)
                      , left = FALSE
          ) |>
          dplyr::group_by(clust) |>
          dplyr::summarise(subspecies = stringr::str_flatten_comma(sort(unique(subspecies)))) |>
          sf::st_make_valid() |>
          dplyr::filter(!grepl(",", subspecies)) |>
          dplyr::group_by(subspecies) |>
          dplyr::summarise() |>
          sf::st_make_valid() |>
          dplyr::mutate(poly = "clust")

      }

    }

    ## all ----

    polys <- mget(ls(pattern = "^tri_dist|^tri_mcp|^tri_clust")) |>
      purrr::map(\(x) {
        if(!is.null(use_crs)) { sf::st_transform(x, crs = use_crs)
        } else sf::st_transform(x, crs = sf::st_crs(polys))
      }
      ) |>
      dplyr::bind_rows() %>%
      {if(nrow(.)) dplyr::group_by(., subspecies) |>
          dplyr::summarise(poly = stringr::str_flatten_comma(sort(unique(poly)))) |>
          sf::st_make_valid() else .}

    # run if polys exist
    if(nrow(polys)) {

      # update binomials ----
      ## new names ----
      new_names <- bi_pres %>%
        sf::st_as_sf(coords = c(pres_x, pres_y), crs = pres_crs, remove = FALSE) %>%
        {if(!is.null(use_crs)) sf::st_transform(., crs = use_crs)
          else sf::st_transform(crs = sf::st_crs(., polys))
        } |>
        sf::st_join(polys, left = FALSE) |>
        sf::st_set_geometry(NULL) |>
        dplyr::group_by(!!rlang::ensym(pres_x), !!rlang::ensym(pres_y)) |>
        dplyr::summarise(subspecies = stringr::str_flatten_comma(sort(unique(subspecies)))) |>
        dplyr::ungroup() |>
        dplyr::mutate(subspecies = ifelse(grepl(paste0(overrides$subspecies, collapse = "|")
                                                , subspecies) & grepl(",", subspecies)
                                          , stringr::str_extract(subspecies
                                                                 , paste0(overrides$subspecies, collapse = "|")
                                          )
                                          , subspecies
        )
        ) |>
        dplyr::filter(!grepl(",", subspecies)) |>
        dplyr::mutate(subspecies = ifelse(subspecies %in% skip_tri
                                          , NA
                                          , subspecies
        )
        , bi_to_tri = TRUE
        ) |>
        dplyr::full_join(bi_pres, by = c(pres_x, pres_y)) |>
        dplyr::bind_rows(tri_pres)

      ## checks ----

      if(FALSE) {

        library(tmap)
        tmap_mode("view")

        test <- bi_pres %>%
          dplyr::left_join(new_names) |>
          sf::st_as_sf(coords = c(pres_x, pres_y), crs = pres_crs, remove = FALSE)

        tm_shape(test)+tm_dots("subspecies")+tm_shape(polys)+tm_borders("subspecies")

        test |>
          sf::st_set_geometry(NULL) |>
          dplyr::count(subspecies) |>
          dplyr::mutate(pc = n/nrow(test)*100)

      }

      ## summary ----

      bi_updated <- new_names |>
        dplyr::filter(bi_to_tri)

      pc <- nrow(bi_updated)/nrow(bi_pres)*100

      message(species, ": ", nrow(bi_updated), " of ", nrow(bi_pres)
              ," binomial occurrences updated to trinomial (", round(pc, 1), "%)")

    } else {

      new_names <- all_pres

      message(species, ": ", "No binomials updated to trinomial due to no trinomial distributions and unable to construct mcp or clust polygons due to < 3 occurrences")

    }

  } else {

    new_names <- all_pres


    if(species %in% skip_bi) {

      message(species, ": ", "No binomials updated to trinomial as skip specified")

    } else if(!nrow(bi_pres)) {

      message(species, ": ", "No binomials updated to trinomial due to no binomial occurrences")

    } else {

      message(species, ": ", "No binomials updated to trinomial due to no trinomial distributions or occurrences")

    }

  }

  res <- new_names |>
    dplyr::mutate(species = species) |>
    dplyr::relocate(species)

  rm(list=setdiff(ls(), "res"))

  return(res)

}
