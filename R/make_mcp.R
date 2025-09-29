
#' Make minimum convex polygon (MCP) around occurrence records
#'
#' Usually used as predict boundary or an Extent of Occurrence (EOO) estimate for a taxa.
#'
#' @param presence Cleaned and filtered dataframe of presences
#' @param out_file Character. Path for the mcp to be saved. Will be saved by
#' `sfarrow::st_write_parquet()`. Currently will not work very well with any
#' full stop in the path. Other file types are changed to .parquet
#' @param force_new Logical. If `out_file` exists, recreate it?
#' @param pres_x,pres_y Character. Name of the columns in `presence` that have
#' the x and y coordinates
#' @param in_crs epsg code for coordinates in `presence`
#' @param out_crs epsg code for coordinates in output mcp. Usually the same as
#' predictors
#' @param buf Distance in metres to buffer the mcp
#' @param clip sf to clip the mcp
#' @param dens_int Numeric. `terra::densify()` `interval` argument ("positive
#' number, specifying the desired minimum distance between nodes. The unit is
#' meter for lonlat data, and in the linear unit of the crs for planar data").
#' Set to `NULL` to not densify.
#'
#' @return sf. `out_file` saved.
#'
#' @export
#'

make_mcp <- function(presence
                     , out_file
                     , force_new = FALSE
                     , pres_x = "long"
                     , pres_y = "lat"
                     , in_crs = 4326
                     , out_crs = in_crs
                     , buf = 0
                     , clip = NULL
                     , dens_int = NULL
                     ) {

  run <- if(file.exists(out_file)) force_new else TRUE

  if(run) {

    fs::dir_create(dirname(out_file))

    suppressMessages({

      sf::sf_use_s2(FALSE)

      res <- presence %>%
        dplyr::distinct(!!rlang::ensym(pres_y), !!rlang::ensym(pres_x)) %>%
        sf::st_as_sf(coords = c(pres_x, pres_y)
                     , crs = in_crs
        ) %>%
        sf::st_union() |>
        sf::st_convex_hull() %>%
        {if(!is.null(dens_int)) (.) |> terra::vect() |> terra::densify(dens_int) |> sf::st_as_sf() else (.)} %>%
        {if(buf != 0) sf::st_buffer(., buf) else .} %>%
        sf::st_make_valid()

      if(!is.null(clip)) {

        if(!is.null(dens_int)) {

          clip <- clip |>
            terra::vect() |>
            terra::densify(dens_int) |>
            sf::st_as_sf()

        }

        if(all(sf::st_intersects(res, sf::st_transform(clip, crs = in_crs), sparse = FALSE))) {

          res <- res %>%
            sf::st_intersection(clip %>%
                                  sf::st_transform(crs = in_crs) %>%
                                  sf::st_make_valid() %>%
                                  dplyr::rename(geometry = attr(., "sf_column")) |>
                                  sf::st_set_geometry("geometry") %>%
                                  sf::st_make_valid()
            )

        }

      }

    })

    res <- res %>%
      sf::st_transform(crs = out_crs) %>%
      sf::st_make_valid()

    if(isTRUE(nrow(res)>0)) suppressWarnings(sfarrow::st_write_parquet(res, out_file)) else res <- NA

  } else {

    res <- sfarrow::st_read_parquet(out_file)

  }

  return(res)

}
