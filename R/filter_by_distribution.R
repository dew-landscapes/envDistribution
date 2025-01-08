#' Filter species occurrence records by spatial layers representing their distributions at state and or national scales.
#'
#' e.g. IUCN Red List or EPBC (Environmental Protection and Biodiversity Conservation act) distributions.
#'
#' @param presence Dataframe of presences for a taxa (i.e. containing x,y coordinate columns).
#' @param state Character. File path/s for the relevant state distributions per taxa.
#' If no relevant state distribution, then use NA and filtering will only be conducted using the national distribution. Currently, only parquet files are accepted.
#' Use sfarrow::st_write_parquet to write sf objects to parquet.
#' @param national Character. File path/s for the relevant national distributions per taxa.
#' If no relevant national distribution, then use NA and filtering will only be conducted using the state distribution,
#' but will maintain all interstate records (if any exist). As with `state`, only parquet files are currently accepted.
#' @param out_file Character. Path for a parquet file to save a pre-filtered dataframe per taxa with records flagged as
#' dist_in = 0 or 1. If previous outputs exist, these outputs will be used to determine what records are new
#' and need to be assessed/filtered, which can save processing time in datasets with many taxa
#' and or complex distributions.
#' @param pres_x,pres_y Character. Name of the columns in `presence` that have the x and y coordinates.
#' @param pres_crs Integer. Coordinate reference system that the presence coordinates correspond to (epsg code).
#' @param use_crs Integer. Coordinate reference system (epsg code) to use as the standard for all objects in spatial joins or intersections.
#' Using a projected coordinate system removes the sf warning around assuming coordinates are planar if they're not, and
#' can sometimes help with spherical geometry issues. If NULL, the crs of the distribution will be used as the standard
#' (i.e. the presence coords and the data_ext_poly will be converted to this crs before joining or intersecting).
#' @param buffer_km Integer. Distance in kilometres used to buffer the distribution and keep records in that radius.
#' Use 0 for no buffer.
#' @param buffer_override Character vector. Vector of character strings representing distribution sources reflected in file paths
#' for sources that should not be buffered. Used to override the buffer for particular distribution sources where buffering may not be relevant.
#' For example, EPBC distributions already include a 'may occur' area, and if these are stored in a directory called 'epbc_dists',
#' then can use this string to override the buffer for taxa with these distributions.
#' @param remove Logical. Remove out of distribution records? If FALSE, out of distribution records will be flagged (see below).
#' @param state_poly sf object. Polygon representing the state boundary for state and interstate record delineation.
#' @param data_ext_poly sf object. Polygon representing the data extent to use for broad level filtering of taxa
#' with distributions that do not intersect the area of interest.
#' @param force_new Logical. Force new taxa outputs? If TRUE, any previous taxa outputs will be ignored and
#' new filtering/flagging of records conducted.
#'
#' @return Presence dataframe with out of range records either removed (if remove == TRUE),
#' or flagged with '1' in the dist_out field (if remove == FALSE). In addition, regardless of `remove`,
#' a pre-filtered dataframe will be written to the `out_file` path that contains records flagged as in or out (1 or 0)
#' in the 'dists_in' field.
#'
#'
#' @export
#'

filter_by_distribution <- function(presence,
                                   state,
                                   national,
                                   out_file,
                                   pres_x = "long",
                                   pres_y = "lat",
                                   pres_crs = 4326,
                                   use_crs = NULL,
                                   buffer_km = 10,
                                   buffer_override = NULL,
                                   remove = TRUE,
                                   state_poly,
                                   data_ext_poly,
                                   force_new = FALSE

){

  # Setup ----

  ## Unique xy data ----
  xy <- presence %>%
    dplyr::distinct(!!rlang::ensym(pres_x), !!rlang::ensym(pres_y))

  # New records function ----

  new_rec <- function(df, out_file){

    old_file <- fs::path(out_file)

    if(file.exists(old_file) & !force_new){

      old <- rio::import(old_file)

      new <- df %>%
        dplyr::anti_join(old, by = c(pres_x, pres_y)) %>%
        dplyr::distinct(!!rlang::ensym(pres_x), !!rlang::ensym(pres_y))

    } else {

      new <- df

    }

    return(new)

  }

  # Distributions join function ----

  dists_join <- function(df, file, interstate, ...){

    if(interstate) sf::sf_use_s2(FALSE) # turn off spherical geometry to overcome loop is not valid errors in international distributions. Even though this was done in creation of the distribution in envVec, it is assumedly needed again here due to the additional summarise by taxa.

    ### distribution ----
    dists <- file %>%
      .[[1]] %>%
      purrr::map(\(x) sfarrow::st_read_parquet(x)) %>% # use map in case of multiple files returned per distribution source due to taxonomy lumping (e.g. Tringa brevipes and Heteroscelus brevipes both have distributions in epbc but are now the same taxa (Tringa bevipes) - gbif is correct to combine them!)
      purrr::list_rbind() %>%
      sf::st_sf() %>%
      sf::st_make_valid() %>%
      {if(!is.null(use_crs)) sf::st_transform(.,crs=use_crs) else .} %>%
      sf::st_make_valid() %>%
      dplyr::summarise() %>%
      sf::st_make_valid() %>%
      dplyr::mutate(dist_in = 1) %>%
      sf::st_set_agr("constant") # set attributes of sf objects to be used in joins and intersections to constant to avoid warnings (from https://github.com/r-spatial/sf/issues/406)

    ### overlap with data extent ----

    data_ext_poly <- data_ext_poly %>%
      {if(!is.null(use_crs)) sf::st_transform(.,crs=use_crs) else sf::st_transform(.,crs=sf::st_crs(dists))} %>%
      sf::st_set_agr("constant")

    dists_de <- dists %>%
      sf::st_join(data_ext_poly
                  , left = FALSE
      ) %>%
      sf::st_set_agr("constant")

    if(nrow(dists_de)>0 & interstate) {

      # extra step to remove any international polygons to improve join processing time

      dists <- dists_de %>%
        sf::st_cast("MULTIPOLYGON") %>%
        sf::st_cast("POLYGON") %>%
        sf::st_join(data_ext_poly, left = FALSE) %>%
        sf::st_make_valid() %>%
        sf::st_set_agr("constant") %>%
        sf::st_intersection(data_ext_poly) %>%
        sf::st_make_valid() %>%
        dplyr::summarise() %>%
        sf::st_make_valid() %>%
        dplyr::mutate(dist_in = 1) %>%
        sf::st_set_agr("constant")

    }

    ### taxa xy ----
    xy_taxa <- df %>%
      sf::st_as_sf(coords = c(pres_x, pres_y), crs = pres_crs, remove=FALSE) %>%
      sf::st_make_valid() %>%
      {if(!is.null(use_crs)) sf::st_transform(.,crs=use_crs) else sf::st_transform(.,crs=sf::st_crs(dists))}

    ### join taxa xy to distribution ----

    if(any(grepl(paste(buffer_override, collapse = "|"), file))) dist_buffer <- 0 else dist_buffer <- buffer_km*1000

    # if distribution overlaps data_ext_poly then join to taxa records, otherwise give taxa records NA to enable filtering later
    xy_dists <- xy_taxa %>%
      {if(nrow(dists_de)>0) sf::st_join(., dists, join = st_is_within_distance, dist = dist_buffer)
        else dplyr::mutate(., dist_in = NA)
      }

    # Update xy_dists with dist_found=1 and convert to tbl df
    xy_dists <- xy_dists %>%
      sf::st_set_geometry(NULL) %>%
      dplyr::mutate(dist_found = 1)

    #rm(list = setdiff(ls(), "xy_dists"))

    return(xy_dists)

  }

  # All distributions join function ----

  all_dists_join <- function(df, out_file, state, national, ...){

    # New records? ----

    if(!force_new) df <- new_rec(df, out_file)

    # Taxon distribution scenarios ----

    if(!nrow(df)>0){

      ### No new records ----

      xy_dists <- rio::import(out_file)

    } else if(all(is.na(state) & is.na(national))){ # use all() to cope with list cols

      ### No distributions ----

      xy_dists <- df %>%
        dplyr::mutate(dist_found = 0,
                      dist_in = NA
        )

    } else if(all(!is.na(national)) & all(is.na(state))){

      ### National distribution only ----

      xy_dists <- dists_join(df
                             , file = national
                             , interstate = TRUE
                             )


    } else if(all(!is.na(state))){

      ### State and or national distributions ----

      state_poly <- state_poly %>%
        {if(!is.null(use_crs)) sf::st_transform(., crs = use_crs) else .} %>%
        dplyr::mutate(state_in = 1)

      state_in_out <- df %>%
        sf::st_as_sf(coords = c(pres_x, pres_y), crs = pres_crs, remove=FALSE) %>%
        sf::st_make_valid() %>%
        {if(!is.null(use_crs)) sf::st_transform(., crs = use_crs) else sf::st_transform(., crs = sf::st_crs(state_poly))} %>%
        sf::st_join(state_poly) %>%
        sf::st_set_geometry(NULL) %>%
        dplyr::mutate(state_in = ifelse(is.na(state_in), 0, state_in)) %>%
        tidyr::nest(df = -c(state_in)) %>%
        dplyr::mutate(file = ifelse(state_in == 1, state, national)
                      , interstate = ifelse(state_in == 1, FALSE, TRUE)
                      , rec = purrr::map_int(df, nrow)
                      )

      state_rec <- 1 %in% state_in_out$state_in # taxa has state data?

      interstate_rec <- 0 %in% state_in_out$state_in # taxa has interstate data?

      if(all(all(is.na(national)), interstate_rec, state_rec)) {

        ## Interstate data but no national distribution ----
        # with state data & state distribution

        xy_dists1 <- state_in_out %>%
          dplyr::filter(state_in == 1) %>%
          purrr::pmap(dists_join) %>%
          purrr::list_rbind()

        xy_dists2 <- state_in_out %>%
          dplyr::filter(state_in == 0) %>%
          .$df %>%
          .[[1]] %>%
          dplyr::mutate(dist_found = 0,
                        dist_in = NA
          )

        xy_dists <- dplyr::bind_rows(xy_dists1, xy_dists2)

      } else {

        ## State and national distributions + equivalent data ----
        # covers state data & distribution only, and both state & national data & distributions

        xy_dists <- state_in_out %>%
          purrr::pmap(dists_join) %>%
          purrr::list_rbind()

      }

    }

    # Join back to previous df ----

    if(all(file.exists(out_file), nrow(df)>0, !force_new)){

      xy_dists <- rio::import(out_file) %>%
        dplyr::anti_join(xy_dists, by = c(pres_x, pres_y)) %>%
        dplyr::bind_rows(xy_dists) %>%
        dplyr::distinct()

    }

    # Write parquet ----
    if(nrow(df)>0) {

      fs::dir_create(basename(out_file))

      arrow::write_parquet(xy_dists, out_file)

    }

    # Clean up to avoid memory issues ----
    rm(list=setdiff(ls(), "xy_dists"))
    gc()

    # Output df ----
    return(xy_dists)

  }


  # Run all_dists_join function per taxa ----
  xy_dists <- all_dists_join(df = xy
                             , out_file = out_file
                             , state = state
                             , national = national
                             )

  # Rejoin to original df & filter or flag outliers ----
  df_filt <- xy_dists %>%
    dplyr::mutate(keep=dplyr::case_when(dist_found==0|is.na(dist_found) ~ 1,
                                        dist_found==1 & is.na(dist_in) ~ 0,
                                        .default = 1
    )
    ) %>%
    {if(remove) dplyr::filter(.,keep==1) else dplyr::mutate(.,dist_out=ifelse(keep==1,0,1))} %>%
    dplyr::select(-tidyr::any_of(c("dist_in","keep")))

  # Clean up to avoid memory issues ----
  rm(list=setdiff(ls(), "df_filt"))
  gc()

  return(df_filt)

}
