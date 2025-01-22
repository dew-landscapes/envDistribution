#' Fix taxonomies based on spatial layers (i.e. in particular geographic areas).
#'
#' If >=2 range layers exist (i.e. one for each taxa), the closest will be used to update the taxonomic name.
#' If only one layer exists, the taxa1 will be changed to the taxa2 according to the change_if_intersects argument.
#'
#' @param bio_df Data frame with taxa column, and lat & long columns.
#' @param coords Character vector. Spatial coordinate columns in order of c("x","y"), e.g. c("long", "lat").
#' @param crs Integer. Coordinate reference system that the coords correspond to (epsg code).
#' @param taxa_df Data frame with taxa1 and taxa2 columns indicating the names of the taxa pairs to be fixed.
#' Single taxa can be input to fix only their records based on to their distribution alone
#' by entering the single taxa name into the taxa1 column, and leaving taxa2 as NA.
#' The taxa names need to correspond to the taxa distributions in `taxa_ds`.
#' @param taxa_ds Data frame. Taxa distribution sources data frame containing 'taxa' and 'file' columns.
#' File column should reference the file path of the distribution for the associated taxa in the 'taxa' column corresponding to taxa1 and taxa2.
#' @param change_if_intersects Logical. If TRUE, change taxonomy if intersects with layer,
#'  or if FALSE, change if doesn't intersect with layer.
#' @param taxa_col Character. Taxa column in bio_df and taxa_ds (must be the same).
#' @param taxonomy Taxonomy object returned by envClean::make_taxonomy in relation to the taxa in bio_df and taxa_ds.
#' @param levels_to_change Character vector. At what level of the taxonomic hierarchy are the taxonomic fixes (either 'species' or 'subspecies' or both).
#' Default is c("species","subspecies").
#'
#' @return Original data frame with corrected taxa names for taxa input into taxa1 and taxa2.
#'
#' @export
#'

fix_spatial_taxonomy <- function(bio_df,
                                 coords = c("long","lat"),
                                 crs = 4326,
                                 taxa_df,
                                 taxa_ds,
                                 change_if_intersects = FALSE,
                                 taxa_col = "taxa",
                                 taxonomy,
                                 levels_to_change = c("species","subspecies")
)
{

  # Find taxonomic fixes (correct names according to distribution)
  taxa_fixes <- purrr::pmap_dfr(taxa_df, function(taxa1, taxa2, ...){

    # Extract relevant distributions ----
    # for each taxa pair (or single)

    lyr <- c(taxa1,taxa2) %>%
      purrr::map(function(x){

        dists <- taxa_ds %>%
          dplyr::filter(!!rlang::ensym(taxa_col) == x) %>%
          dplyr::pull(file) %>%
          purrr::map(~ sfarrow::st_read_parquet(.x) %>%
                       sf::st_transform(crs = crs) %>%
                       dplyr::summarise() %>%
                       sf::st_make_valid() %>%
                       {if(attr(., "sf_column") != "geometry") dplyr::rename(., geometry = attr(., "sf_column")) %>%
                           sf::st_set_geometry("geometry")
                         else .
                       }
          ) %>% # use map in case of multiple files returned per distribution source
          purrr::list_rbind() %>%
          sf::st_sf() %>%
          dplyr::summarise() %>%
          sf::st_make_valid() %>%
          dplyr::mutate(near_name = x
                        , overlap = 1
          )

        return(dists)

      }
      ) %>%
      purrr::list_rbind() %>%
      sf::st_sf()

    # relevant taxa data ----

    taxa_xy <- bio_df %>%
      {if(nrow(lyr)>1) dplyr::filter(.,grepl(paste(c(taxa1,taxa2),collapse = "|"),!!rlang::ensym(taxa_col))) else .} %>%
      {if(nrow(lyr)==1) dplyr::filter(.,grepl(taxa1,!!rlang::ensym(taxa_col))) else .} %>%
      dplyr::distinct(dplyr::pick(tidyr::all_of(c(taxa_col, coords)))) %>%
      {if("subspecies" %in% levels_to_change) dplyr::left_join(., taxonomy$subspecies$lutaxa, by = taxa_col)
        else dplyr::left_join(., taxonomy$species$lutaxa, by = taxa_col)
      } %>%
      dplyr::filter(returned_rank %in% levels_to_change) %>%
      dplyr::distinct(dplyr::pick(tidyr::all_of(c(taxa_col,coords)))) %>%
      sf::st_as_sf(coords = coords, remove = FALSE, crs = crs)

    # find records in any overlapping areas between two different distributions ----

    if(nrow(lyr)>1){

      lyr1 <- lyr %>%
        dplyr::slice(1)

      lyr2 <- lyr %>%
        dplyr::slice(2)

      if(any(sf::st_intersects(lyr1,lyr2, sparse = FALSE))) {

        lyr_overlap <- sf::st_intersection(lyr1, lyr2) %>%
          sf::st_make_valid() %>%
          dplyr::summarise() %>%
          sf::st_make_valid()

        rec_in_overlap <- taxa_xy %>%
          sf::st_join(lyr_overlap, left = FALSE) %>%
          sf::st_set_geometry(NULL) %>%
          dplyr::mutate(near_name = !!rlang::ensym(taxa_col)
                        , overlap = 1
                        )

      } else {

        rec_in_overlap <- data.frame()

      }

    } else {

      rec_in_overlap <- data.frame()

    }

    # taxa record and distribution overlap (proximity) ----

    taxa_overlap <- taxa_xy %>%
      {if(nrow(rec_in_overlap)>0) dplyr::anti_join(., rec_in_overlap) else .} %>% # remove any records in areas where distributions overlap, so won't have their name changed (as cannot decide which taxa is correct in an overlapping area)
      {if(nrow(lyr)>1) sf::st_join(.,lyr, join = st_nearest_feature) else sf::st_join(.,lyr)} %>%
      sf::st_set_geometry(NULL) %>%
      {if(nrow(rec_in_overlap)>0) dplyr::bind_rows(., rec_in_overlap) else .} # re-instate any records in areas where distributions overlap to preserve original taxa names

    # correct names ----

    if(nrow(lyr)>1){

      taxa_fix <- taxa_overlap %>%
        dplyr::mutate(.,correct_name = stringr::str_replace(!!rlang::ensym(taxa_col)
                                                            , pattern = paste(c(taxa1,taxa2),collapse = "|")
                                                            , replacement = near_name
        )
        ) %>%
        dplyr::select(tidyr::any_of(c(coords, taxa_col, "correct_name")))

    } else {

      taxa_fix <- taxa_overlap %>%
        {if(change_if_intersects) dplyr::mutate(.,correct_name = dplyr::if_else(overlap==1,
                                                                                stringr::str_replace(!!rlang::ensym(taxa_col)
                                                                                                     , pattern = paste(c(taxa1,
                                                                                                                         taxa2),
                                                                                                                       collapse = "|")
                                                                                                     , replacement = near_name
                                                                                ),
                                                                                !!rlang::ensym(taxa_col)
        )
        ) else .} %>%
        {if(!change_if_intersects) dplyr::mutate(.,correct_name = dplyr::if_else(is.na(overlap),
                                                                                 stringr::str_replace(!!rlang::ensym(taxa_col)
                                                                                                      , pattern = paste(c(taxa1,
                                                                                                                          taxa2),
                                                                                                                        collapse = "|")
                                                                                                      , replacement = near_name
                                                                                 ),
                                                                                 !!rlang::ensym(taxa_col)
        )
        ) else .} %>%
        dplyr::select(tidyr::any_of(c(coords, taxa_col, "correct_name")))

    }

    return(taxa_fix)
  }
  )

  # Join back to original data and update names ----

  tax_fix_all <- bio_df %>%
    dplyr::left_join(taxa_fixes) %>%
    dplyr::mutate(!!rlang::ensym(taxa_col) := dplyr::if_else(!is.na(correct_name),
                                                             correct_name,
                                                             !!rlang::ensym(taxa_col)
    )
    ) %>%
    dplyr::select(-correct_name)

  return(tax_fix_all)

}
