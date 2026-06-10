
#' Convert binomial (species) level occurrences to trinomial based on distributions.
#'
#' @param species Character. Name of the species (binomial) to update occurrences to trinomial for.
#' @param presences Data frame of all binomial and trinomial presences relevant to the species.
#' Must contain a 'subspecies' column and x, y coordinate columns corresponding to the names
#' in `pres_x` & `pres_y`. For binomial presences 'subspecies' should be NA.
#' @param distrib_files Data frame of relevant distribution file paths for all the subspecies
#' (e.g. for distributions sourced from redlist, epbc etc). This must contain a 'subspecies' field, and
#' state' and 'national' fields for state and national distribution files. If no relevant distribution, then use NA.
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
#' @param buf_override Character vector. Vector of character strings representing distribution sources reflected in file paths
#' for sources that should not be buffered. Used to override the buffer for particular distribution sources where buffering may not be relevant.
#' For example, EPBC distributions already include a 'may occur' area, and if these are stored in a directory called 'epbc_dists',
#' then can use this string to override the buffer for taxa with these distributions.
#' @param state_poly sf object. Optional polygon representing the state boundary for erasing the state
#' section of a national distribution and replacing it with the state distribution.
#' For use where state distributions are of greater resolution or accuracy.
#' @param clust_dist Integer. Distance to base clusters on in clust function if `use_clust` = TRUE.
#' Needs to be in the units of the projected crs supplied in `use_crs`, e.g. metres.
#' @param overlap_thres Integer. Percent threshold corresponding to the area of overlap between distributions
#' at which to not update binomial occurrences to trinomial, i.e. with a threshold of 50%, if a trinomial
#' distribution overlaps any others by more than 50% then no binomials will be updated to that trinomial.
#' @param overrides Character vector containing any trinomials that should override the overlap checking
#' and force binomials within their distribution polygons to be updated to that trinomial even if the
#' binomials overlap distribution polygons of other trinomials.
#' @param skip_bi Character vector of any binomials (species) to be skipped and not have binomials
#' updated to trinomial.
#' @param skip_tri Character vector of any trinomials (subspecies, varieties, forms etc) to be skipped
#' and not have binomials updated to these trinomials.
#'
#' @details
#' Binomial occurrence records won't be updated to trinomial if they fall within multiple
#' trinomial distributions and the overlap is beyond the `overlap_thres`.
#' This avoids updating binomials to the wrong trinomial where distributions of subspecies overlap.
#' Distribution refers generally to all the polygons used to represent distributions,
#' i.e. distribution files (e.g. from redlist or epbc), MCPs, or cluster polygons.
#'
#' @return Data frame equivalent to `presences` with 'subspecies' updated for taxa where
#' distributions were found and did not overlap beyond the `overlap_thres` with multiple
#' trinomials distributions, and an additional 'bi_to_tri' column indicating which records
#' were updated.
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
                      , buf_override = NULL
                      , state_poly = NULL
                      , clust_dist = 50000
                      , overlap_thres = 50
                      , overrides = NULL
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

          dist_prep <- distrib_files |>
            dplyr::filter(subspecies == x)

          if(!is.na(dist_prep$state)) {

            dist <- dist_prep |>
              dplyr::pull(state) |>
              purrr::map(\(f) sfarrow::st_read_parquet(f)) |>
              dplyr::bind_rows() |>
              dplyr::summarise() |>
              sf::st_make_valid() |>
              sf::st_make_valid() %>%
              {if(!is.null(use_crs)) sf::st_transform(.,crs = use_crs) else .} %>%
              {if(!is.null(state_poly))
                sf::st_intersection(., state_poly |>
                                      dplyr::select(geometry) |>
                                      sf::st_make_valid() %>%
                                      {if(!is.null(use_crs)) sf::st_transform(., crs = use_crs) else .}
                ) |>
                  sf::st_make_valid() |>
                  nngeo::st_remove_holes() |>
                  sf::st_make_valid() else .} |>
              dplyr::mutate(dist_type = "state")

          }

          if(!is.na(dist_prep$national)) {

            dist <- dist_prep |>
              dplyr::pull(national) %>%
              .[[1]] |>
              purrr::map(\(f) sfarrow::st_read_parquet(f)) |>
              dplyr::bind_rows() |>
              dplyr::summarise() |>
              sf::st_make_valid() %>%
              {if(!is.null(use_crs)) sf::st_transform(.,crs = use_crs) else .} |>
              sf::st_make_valid() %>%
              {if(!is.null(state_poly) & !is.na(dist_prep$state))
                sf::st_difference(., state_poly |>
                                    dplyr::select(geometry) |>
                                    sf::st_make_valid() %>%
                                    {if(!is.null(use_crs)) sf::st_transform(., crs = use_crs) else .}
                ) |>
                  dplyr::bind_rows(dist) |>
                  dplyr::summarise() |>
                  sf::st_make_valid() |>
                  # remove any polygons that don't intersect the state dist (removes any polygons created by buffering slivers that can't be removed by other means, e.g. Malurus cyaneus leggei)
                  sf::st_cast("MULTIPOLYGON") |>
                  sf::st_cast("POLYGON") |>
                  sf::st_join(dist, left = FALSE) |>
                  dplyr::summarise() |>
                  sf::st_make_valid() |>
                  nngeo::st_remove_holes() |>
                  sf::st_make_valid()
                else .} |>
              dplyr::mutate(dist_type = "national")
          }

          dist <- dist |>
            #sf::st_buffer(buf) |>
            dplyr::mutate(subspecies = x
                          , poly = "dist"
            )

        }
        ) |>
        dplyr::bind_rows() |>
        sf::st_make_valid()

      ### buffer non-overlapping parts of distribution polygons ----

      if(all(nrow(tri_dist) > 1, buf > 0)) {

        tri_dist <- tri_dist |>
          sf::st_drop_geometry() |>
          dplyr::distinct(subspecies) |>
          dplyr::pull(subspecies) |>
          purrr::map(\(s) {

            dist_prep <- distrib_files |>
              dplyr::filter(subspecies == s)

            if(any(grepl(paste(buf_override, collapse = "|"), dist_prep$national))
               |any(grepl(paste(buf_override, collapse = "|"), dist_prep$state))
            ) buf <- 0 else buf

            this_ssp_poly <- tri_dist |>
              dplyr::filter(subspecies == s)

            other_ssp_polys <- tri_dist |>
              dplyr::filter(subspecies != s) |>
              dplyr::summarise() |>
              sf::st_make_valid()

            if(FALSE) {

              library(tmap)
              tmap_mode("view")
              tm_shape(other_ssp_polys)+tm_polygons(fill = "blue")+tm_shape(this_ssp_poly)+tm_polygons(fill = "red")

            }

            within <- any(sf::st_within(this_ssp_poly
                                        , other_ssp_polys |>
                                          sf::st_buffer(1) # small buffer to encapsulate single regions completely within multi-regions that otherwise doesn't trigger st_within, e.g. KI & MLR herbarium regions within whole state.
                                        , sparse = FALSE
            )
            )

            if(!within) {

              this_ssp_buf <- this_ssp_poly %>%
                # try sf::st_difference if rmapshaper::ms_erase errors with 'Not compatible with STRSXP: [type=list].'
                # can't use sf::st_difference for all, as doesn't return the correct result for some ssp
                {tryCatch(expr = rmapshaper::ms_erase(., other_ssp_polys, remove_slivers = TRUE)
                          , error = function(e) sf::st_difference(., other_ssp_polys))} |> # need an extra erase here as otherwise can buffer some edges that overlap other dists
                sf::st_make_valid() |>
                sf::st_buffer(buf) |>
                sf::st_make_valid() %>%
                {tryCatch(expr = rmapshaper::ms_erase(., other_ssp_polys, remove_slivers = TRUE)
                          , error = function(e) sf::st_difference(., other_ssp_polys))} |>
                sf::st_make_valid() |>
                dplyr::bind_rows(this_ssp_poly) |>
                dplyr::group_by(subspecies, poly, dist_type) |>
                dplyr::summarise() |>
                dplyr::ungroup() |>
                sf::st_make_valid()

            } else {

              this_ssp_buf <- this_ssp_poly

            }

          }
          ) |>
          dplyr::bind_rows()

      } else {

        tri_dist <- tri_dist |>
          sf::st_buffer(buf)

      }

    }

    ## mcp ----

    if(use_mcp) {

      if(is.null(mcp_files) & nrow(tri_pres)) {

        mcp_prep <- tri_pres %>%
          dplyr::rename(tidyr::any_of(c("long" = "cell_long", "lat" = "cell_lat"))) |> # to overcome odd error of pres_x and pres_y not being accepted in make_mcp below
          dplyr::add_count(subspecies) |>
          dplyr::filter(n >= 3) |>
          dplyr::select(-n) %>%
          # don't use mcp for subspecies with national distributions
          {if(exists("tri_dist")) dplyr::anti_join(., tri_dist |>
                                                     sf::st_drop_geometry() |>
                                                     dplyr::filter(dist_type == "national")
          ) else .} |>
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
              ) #|>
            #sf::st_buffer(buf)

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

    # don't use clust where all ssp (in presences & distributions) have national distributions
    if(use_clust) {

      use_clust <- tri_pres |>
        dplyr::distinct(subspecies) |>
        dplyr::bind_rows(tri_dist |>
                           sf::st_drop_geometry() |>
                           dplyr::distinct(subspecies)
        ) |>
        dplyr::anti_join(tri_dist |>
                           sf::st_drop_geometry() |>
                           dplyr::filter(dist_type == "national") |>
                           dplyr::distinct(subspecies)
        ) %>%
        nrow(.) > 0

    }

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
                                         , h = clust_dist
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
            dplyr::mutate(clust = y) #|>
          #sf::st_buffer(buf)

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

    ## overlaps ----

    # run if polys exist
    if(nrow(polys)) {

      overlaps <- polys |>
        dplyr::mutate(pc_overlap =
                        purrr::map_dbl(subspecies
                                       , \(x) {

                                         poly <- polys |>
                                           dplyr::filter(subspecies == x)

                                         if(nrow(poly) > 1) {

                                           other_poly <- polys |>
                                             dplyr::filter(subspecies != x) |>
                                             dplyr::summarise() |>
                                             sf::st_make_valid()

                                           if(nrow(other_poly)) {

                                             poly_area <- poly %>%
                                               dplyr::mutate(area = as.numeric(sf::st_area(.))) |>
                                               sf::st_set_geometry(NULL) |>
                                               dplyr::pull(area)

                                             overlap_area <- poly |>
                                               sf::st_intersection(other_poly) |>
                                               sf::st_make_valid() %>%
                                               dplyr::mutate(area = as.numeric(sf::st_area(.))) |>
                                               sf::st_set_geometry(NULL) |>
                                               dplyr::pull(area) %>%
                                               {if(length(.)) . else 0}

                                             overlap_area/poly_area*100

                                           } else 0

                                         } else 0

                                       }
                        )
        ) |>
        sf::st_set_geometry(NULL)

      # update binomials ----
      ## new names ----
      new_names <- bi_pres %>%
        sf::st_as_sf(coords = c(pres_x, pres_y), crs = pres_crs, remove = FALSE) %>%
        {if(!is.null(use_crs)) sf::st_transform(., crs = use_crs)
          else sf::st_transform(crs = sf::st_crs(., polys))
        } |>
        sf::st_join(polys, left = FALSE) |>
        sf::st_set_geometry(NULL) %>%
        {if(nrow(.)) dplyr::group_by(., !!rlang::ensym(pres_x), !!rlang::ensym(pres_y)) |>
            dplyr::summarise(subspecies = stringr::str_flatten_comma(sort(unique(subspecies)))) |>
            dplyr::ungroup() |>
            dplyr::left_join(overlaps) |>
            dplyr::mutate(subspecies = dplyr::case_when(grepl(paste0(overrides, collapse = "|")
                                                              , subspecies) & grepl(",", subspecies)
                                                        ~ stringr::str_extract(subspecies
                                                                               , paste0(overrides, collapse = "|")
                                                        )
                                                        , pc_overlap > overlap_thres ~ NA
                                                        , subspecies %in% skip_tri ~ NA
                                                        , .default = subspecies
            )
            ) |>
            dplyr::filter(!grepl(",", subspecies)
                          , !is.na(subspecies)
            ) |>
            dplyr::mutate(bi_to_tri = TRUE)
          else dplyr::mutate(., bi_to_tri = NA)
        } |>
        dplyr::select(tidyr::any_of(c("subspecies", "cell_long", "cell_lat", "bi_to_tri")))

      res <- new_names |>
        dplyr::full_join(bi_pres, by = c(pres_x, pres_y)) |>
        dplyr::bind_rows(tri_pres)

      ## checks ----

      if(FALSE) {

        library(tmap)
        tmap_mode("view")

        test <- bi_pres %>%
          dplyr::left_join(new_names) |>
          sf::st_as_sf(coords = c(pres_x, pres_y), crs = pres_crs, remove = FALSE)

        tm_shape(test)+tm_dots("subspecies")+tm_shape(polys)+tm_borders("subspecies")+tm_scale_bar()

        test |>
          sf::st_set_geometry(NULL) |>
          dplyr::count(subspecies) |>
          dplyr::mutate(pc = n/nrow(test)*100)

      }

      ## summary ----

      bi_updated <- res |>
        dplyr::filter(bi_to_tri)

      pc <- nrow(bi_updated)/nrow(bi_pres)*100

      message(species, ": ", nrow(bi_updated), " of ", nrow(bi_pres)
              ," binomial occurrences updated to trinomial (", round(pc, 1), "%)")

    } else {

      res <- all_pres

      message(species, ": ", "No binomials updated to trinomial due to no trinomial distributions and unable to construct mcp or clust polygons due to < 3 occurrences")

    }

  } else {

    res <- all_pres


    if(species %in% skip_bi) {

      message(species, ": ", "No binomials updated to trinomial as skip specified")

    } else if(!nrow(bi_pres)) {

      message(species, ": ", "No binomials updated to trinomial due to no binomial occurrences")

    } else {

      message(species, ": ", "No binomials updated to trinomial due to no trinomial distributions or occurrences")

    }

  }

  res <- res |>
    dplyr::mutate(species = species) |>
    dplyr::relocate(species)

  rm(list=setdiff(ls(), "res"))

  return(res)

}
