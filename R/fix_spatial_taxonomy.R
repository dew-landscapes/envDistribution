#' Fix taxonomies based on spatial layers (i.e. in particular geographic areas).
#'
#' If >=2 range layers exist (i.e. one for each taxa), the closest will be used to update the taxonomic name.
#' If only one layer exists, the taxa1 will be changed to the taxa2 according to the change_if_intersects argument.
#'
#' @param bio_df Data frame with taxa column, and lat & long columns.
#' @param coords Character vector. Spatial coordinate columns in order of c("x","y"), e.g. c("long", "lat").
#' @param crs Integer. Coordinate reference system that the coords correspond to (epsg code).
#' @param taxa1 Character. First taxa to change (paired with taxa 2).
#' @param taxa2 Character. Second taxa to change (paired with taxa 1).
#' @param taxa_ds Data frame. Taxa distribution sources data frame containing 'taxa' and 'file' columns.
#' File column should reference the file path of the distribution for the associated taxa in the 'taxa' column corresponding to taxa1 and taxa2.
#' @param change_if_intersects Logical. If TRUE, change taxonomy if intersects with layer,
#'  or if FALSE, change if doesn't intersect with layer.
#' @param taxa_col Character. Taxa column in bio_df.
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
                                 taxa1,
                                 taxa2,
                                 taxa_ds,
                                 change_if_intersects = FALSE,
                                 taxa_col = "taxa",
                                 taxonomy,
                                 levels_to_change = c("species","subspecies")
)
{

  # Find taxonomic fixes (correct names according to distribution)
  taxa_fixes <- purrr::pmap_dfr(list(taxa1,taxa2), function(taxa1,taxa2){

    # Extract relevant distributions for each taxa pair (or single)
    lyr <- c(taxa1,taxa2) %>%
      purrr::map(function(x){

        dists <- taxa_ds %>%
          dplyr::filter(original_name==x) %>%
          dplyr::pull(file) %>%
          purrr::map(~ sfarrow::st_read_parquet(.x)) %>% # use map in case of multiple files returned per distribution source due to gbif taxonomy lumping (e.g. Tringa brevipes and Heteroscelus brevipes both have distributions in epbc but are now the same taxa (Tringa bevipes) - gbif is correct to combine them!)
          purrr::list_rbind() %>%
          sf::st_sf() %>%
          dplyr::summarise() %>%
          sf::st_make_valid() %>%
          sf::st_transform(crs = crs) %>%
          dplyr::mutate(near_name = x
                        , overlap = 1
                        )

        return(dists)

      }
      ) %>%
      purrr::list_rbind() %>%
      sf::st_sf()

    taxa_df <- bio_df %>%
      {if(nrow(lyr)>1) dplyr::filter(.,grepl(paste(c(taxa1,taxa2),collapse = "|"),!!rlang::ensym(taxa_col))) else .} %>%
      {if(nrow(lyr)==1) dplyr::filter(.,grepl(taxa1,!!rlang::ensym(taxa_col))) else .} %>%
      dplyr::select(-tidyr::any_of("returned_rank")) %>%
      {if("subspecies" %in% levels_to_change) dplyr::left_join(.,taxonomy$subspecies$lutaxa,by=setNames("original_name",taxa_col))
        else dplyr::left_join(.,taxonomy$species$lutaxa,by=setNames("original_name",taxa_col))
      } %>%
      dplyr::filter(returned_rank %in% levels_to_change)

    taxa_overlap <- taxa_df %>%
      dplyr::distinct(dplyr::pick(tidyr::all_of(c(taxa_col,coords)))) %>%
      sf::st_as_sf(coords = coords, remove = FALSE, crs = crs) %>%
      {if(nrow(lyr)>1) sf::st_join(.,lyr, join = st_nearest_feature) else sf::st_join(.,lyr)} %>%
      sf::st_set_geometry(NULL)

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

  # Join back to original data and update names
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
