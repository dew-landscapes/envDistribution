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
#' @param buffer_km Integer. Distance in kilometres used to buffer the `lyr` and attribute records in that radius.
#' Use 0 for no buffer.
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
                        buffer_km = 0

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
                  sf::st_buffer(buffer_km) |>
                  sf::st_make_valid()
    ) |>
    sf::st_drop_geometry() |>
    dplyr::distinct(across(!!rlang::ensym(df_x), !!rlang::ensym(df_y)
                           , tidyr::any_of(att_cols)
    )
    ) %>%
    {if(!is.null(renames)) dplyr::rename(., tidyr::any_of(renames)) else .}

  res <- df |>
    dplyr::left_join(xy_att)

  return(res)

}
