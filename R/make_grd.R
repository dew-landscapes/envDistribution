
#' Make grid around occurrence records
#'
#' Usually used as Area of Occupancy (AOO) estimate for a taxa.
#'
#' @param presence Cleaned and filtered dataframe of presences.
#' @param out_file Character. Path for the mcp to be saved. Will be saved by
#' `sfarrow::st_write_parquet()`. Currently will not work very well with any
#' full stop in the path. Other file types are changed to .parquet.
#' @param force_new Logical. If `out_file` exists, recreate it?
#' @param pres_x,pres_y Character. Name of the columns in `presence` that have
#' the x and y coordinates.
#' @param pres_taxa Character. Name of the column in `presence` that
#' contains the taxa names. Needed for ConR (see below).
#' @param in_crs epsg code for coordinates in `presence`.
#' @param out_crs epsg code for coordinates in output grid.
#' @param use_ConR Use the ConR package for creating the grid.
#' This uses ConR::AOO.computing to find the grid with the lowest number of cells,
#' but can be very slow for common taxa at national level.
#' Otherwise, sf::st_make_grid is used. Default is FALSE.
#' @param cell_size Grid cell size in km,
#' e.g. the standard IUCN redlist cell size is 2 km (i.e. a 2 x 2 km square),
#' which is the default.
#' @param num_rast_pos Number of raster start positions when ConR = TRUE.
#' See nbe.rep.rast.AOO parameter in ConR::AOO.computing for details.
#'
#' @return sf. `out_file` saved.
#'
#' @export
#'

make_grd <- function(presence
                     , out_file
                     , force_new = FALSE
                     , pres_x = "long"
                     , pres_y = "lat"
                     , pres_taxa = "taxa"
                     , in_crs = 4326
                     , out_crs = in_crs
                     , use_ConR = FALSE
                     , cell_size = 2
                     , num_rast_pos = 0
) {

  run <- if(file.exists(out_file)) force_new else TRUE

  out_file <- gsub(tools::file_ext(out_file), "", out_file)
  out_file <- gsub("\\.$", "", out_file)
  out_file <- paste0(out_file, ".parquet")

  if(run) {

    fs::dir_create(dirname(out_file))

    res <- presence %>%
      dplyr::distinct(!!rlang::ensym(pres_y), !!rlang::ensym(pres_x), !!rlang::ensym(pres_taxa)) %>%
      {if(use_ConR) ConR::AOO.computing(XY = .
                                        , cell_size_AOO = cell_size
                                        , show_progress = FALSE
                                        , export_shp = TRUE
                                        , proj_type = out_crs
                                        , nbe.rep.rast.AOO = num_rast_pos
      ) %>%
          .$AOO_poly
        else sf::st_as_sf(.,coords = c(pres_x, pres_y)
                       , crs = in_crs
          ) %>%
          sf::st_make_grid()
      } %>%
      sf::st_transform(crs = out_crs) %>%
      dplyr::select(tidyr::any_of(contains(c("geometry","shape"))))

    sfarrow::st_write_parquet(res, out_file)

  } else {

    res <- sfarrow::st_read_parquet(out_file)

  }

  return(res)

}
