#' Add spatial attribute/s
#'
#' Adds attributes from a spatial layer to a data frame,
#' e.g. add indigenous or vagrant status from a distribution layer to a data frame of taxa occurrences.
#'
#' @param df Data frame containing coordinate columns.
#' @param lyr sf object with polygons containing the attributes to be added to `df`.
#' @param att_cols Character vector of attribute columns from `lyr` to add to `df`.
#' @param renames Named character vector of attribute columns to rename for the output,
#' e.g. c("spatial_ind" = "ind", "spatial_ind" = "isIndigenous").
#' @param df_x,df_y Character. Name of the columns in `df` that have the x and y coordinates.
#' @param df_crs Integer. Coordinate reference system that the coordinates in `df` correspond to (epsg code).
#' @param use_crs Integer. Coordinate reference system (epsg code) to use as the standard for all objects in spatial joins.
#' Using a projected coordinate system removes the sf warning around assuming coordinates are planar if they're not, and
#' can sometimes help with spherical geometry issues. If NULL, the crs of the distribution will be used as the standard
#' (i.e. the `df` coords and `lyr` will be converted to this crs before joining).
#' @param maxdist_km Integer. Maximum distance in kilometres to attribute occurrence records outside `lyr`.
#' Occurrence records that don't fall within the polygons in `lyr` but within `maxdist_km` will be attributed by
#' the attributes in the nearest polygon in `lyr`.
#'
#' @return Data frame equivalent to `df` with specified attributes from `lyr`.
#'
#'
#' @export
#'

add_spt_att <- function(df,
                        lyr,
                        att_cols,
                        renames = NULL,
                        df_x = "long",
                        df_y = "lat",
                        df_crs = 4326,
                        use_crs = NULL,
                        maxdist_km = NULL

){

  xy_att <- df %>%
    dplyr::distinct(!!rlang::ensym(df_x), !!rlang::ensym(df_y)) |>
    sf::st_as_sf(coords = c(df_x, df_y)
                 , crs = df_crs
                 , remove = FALSE
    ) |>
    sf::st_transform(crs = use_crs) |>
    sf::st_join(lyr |>
                  sf::st_transform(crs = use_crs) |>
                  sf::st_make_valid()
    )

  out_of_lyr <- xy_att |>
    dplyr::filter(dplyr::if_any(tidyr::any_of(att_cols), \(x) is.na(x))) |>
    dplyr::select(tidyr::all_of(c(df_x, df_y)))

  if(!is.null(maxdist_km) & maxdist_km > 0 & nrow(out_of_lyr) > 0) {

    library(sf)

    out_of_lyr <- out_of_lyr %>%
      sf::st_join(lyr |>
                    dplyr::summarise() |>
                    sf::st_buffer(maxdist_km*1000)
                  , left = FALSE
      ) |>
      sf::st_join(lyr, join = st_nearest_feature) |>
      sf::st_drop_geometry()

  }

  xy_att <- xy_att |>
    sf::st_drop_geometry() |>
    {if(!is.null(maxdist_km) & maxdist_km > 0 & nrow(out_of_lyr) > 0) dplyr::anti_join(., out_of_lyr
                                                                                       , by = c(df_x, df_y)
    ) |>
        dplyr::bind_rows(out_of_lyr) else .} |>
    dplyr::distinct(dplyr::across(tidyr::any_of(c(df_x, df_y, att_cols)))) %>%
    {if(!is.null(renames)) dplyr::rename(., tidyr::any_of(renames)) else .}

  res <- df |>
    dplyr::left_join(xy_att)

  return(res)

}
