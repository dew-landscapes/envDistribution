#' Calculate a species environment category based on occurrences.
#'
#' !Experimental! In future this can also be based on distributions.
#'
#' @param presence Cleaned and filtered dataframe of presences
#' @param land_poly sf polygon object indicating the land in the area of interest.
#' @param thresh Integer. Threshold for percent records to determine if a species is terrestrial, marine or both.
#' @param out_file Character. File path for the results files to be saved.
#' @param force_new Logical. If `out_file` exists, recreate it?
#' @param pres_x,pres_y Character. Name of the columns in `presence` that have
#' the x and y coordinates.
#' @param in_crs epsg code for coordinates in `presence`.
#' @param out_crs epsg code for coordinates in output mcp.
#' @param buf Distance in metres to buffer the `land_poly`.
#'
#' @return data frame saved as rds to `out_file` path containing environ and land_pc (percent land overlap ) fields for the taxa.
#'
#' @export
#'

calc_environ <- function(presence
                         , land_poly
                         , thresh = 95
                         , out_file
                         , force_new = FALSE
                         , pres_x = "long"
                         , pres_y = "lat"
                         , in_crs = 4326
                         , out_crs = in_crs
                         , buf = 0
)
{

  run <- if(file.exists(out_file)) force_new else TRUE

  out_file <- gsub(paste0(tools::file_ext(out_file),"$"), "", out_file)
  out_file <- gsub("\\.$", "", out_file)
  out_file <- paste0(out_file, ".rds")

  if(run) {

    fs::dir_create(dirname(out_file))

    pres <- presence %>%
      dplyr::distinct(!!rlang::ensym(pres_y), !!rlang::ensym(pres_x)) %>%
      sf::st_as_sf(coords = c(pres_x, pres_y)
                   , crs = in_crs
      ) %>%
      sf::st_transform(crs = out_crs)

    land <- land_poly %>%
      sf::st_transform(crs = out_crs) |>
      dplyr::mutate(environ = "land") |>
      sf::st_buffer(buf)

    res <- pres |>
      sf::st_join(land) |>
      sf::st_set_geometry(NULL) |>
      dplyr::mutate(environ = ifelse(is.na(environ), "ocean", environ)
                    , total = length(environ)
                    ) |>
      dplyr::count(environ, total) %>%
      {if(!"land" %in% .$environ) dplyr::bind_rows(., tibble::tibble(environ = "land", total = unique(.$total))) else .} %>%
      {if(!"ocean" %in% .$environ) dplyr::bind_rows(., tibble::tibble(environ = "ocean", total = unique(.$total))) else .} |>
      dplyr::group_by(environ, total) |>
      dplyr::summarise(n = sum(n, na.rm = TRUE)) |>
      dplyr::ungroup() |>
      dplyr::mutate(percent = n/total*100) |>
      tidyr::pivot_wider(names_from = environ, values_from = c(n, percent), values_fill = 0) |>
      dplyr::mutate(environ = dplyr::case_when(percent_land >= thresh ~ "terrestrial",
                                               percent_ocean >= thresh ~ "marine",
                                               percent_land < thresh & percent_ocean < thresh ~ "both",
                                               .default = NA
                                               )
                    )

    rio::export(res, out_file, format = "rds")

  } else {

    res <- rio::import(out_file, format = "rds")

  }

  return(res)
}
