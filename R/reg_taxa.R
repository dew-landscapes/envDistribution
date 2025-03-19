
#' Return relevant taxa to a region of interest based on occurrences and or distributions.
#'
#' Useful for species level analyses seeking to link to relevant trinomials (subspecies, races, varieties, forms, etc),
#' where there may be other irrelevant trinomials occurring outside the region of interest.
#'
#' @param bio_df Data frame with taxa column, and lat & long columns.
#' @param out_dir Character. Directory path for the results files to be saved (summary_stats.rds and mcp_noreg.parquet).
#' mcp_noreg will be saved by `sfarrow::st_write_parquet()`. Currently will not work very well with any
#' full stop in the path. Other file types are changed to .parquet.
#' @param region_bound sf object. Region boundary polygon.
#' @param force_new Logical. If any of the results files in the `out_dir` exists, recreate it?
#' @param pres_x,pres_y Character. Name of the columns in `presence` that have
#' the x and y coordinates.
#' @param in_crs epsg code for coordinates in `presence`.
#' @param out_crs epsg code for coordinates in output grid.
#' @param buf Distance in metres to buffer the `region_bound`.
#' @param taxa_col Character. Taxa column in bio_df and taxa_ds (must be the same).
#' @param taxonomy Taxonomy object returned by envClean::make_taxonomy in relation to the taxa in bio_df and taxa_ds.
#'
#' @return summary_stats.rds with aoo and eoo stats,
#' and mcp_noreg.parquet with the mcp excluding regional presences.
#'
#' @export
#'

reg_taxa <- function(bio_df
                     , out_dir
                     , region_bound
                     , force_new = FALSE
                     , pres_x = "long"
                     , pres_y = "lat"
                     , in_crs = 4326
                     , out_crs = in_crs
                     , buf = 0
                     , taxonomy = targets::tar_read(taxonomy
                                                    , store = fs::path("H:/dev/out/envCleaned"
                                                                       , "sa_br_dissolve______0__P50Y__sa_br_dissolve"
                                                                       , "90__90__P1Y__subspecies"
                                                                       , "clean"
                                                    ))
                     , listed_df = arrow::read_parquet(fs::path("H:","data","taxonomy","all_status_galah.parquet"))
                     , taxa_ds = dists_source(distrib_dir = data_dir,
                                               sources = c("epbc","expert"),
                                               source_rank = FALSE,
                                               datatype = "vector"
                     )
) {

  reg_taxa_file <- fs::path(out_dir, "reg_taxa.rds")

  run <- if(file.exists(reg_taxa_file)) force_new else TRUE

  if(run) {

    # All data ----

    taxa_data <- bio_df %>%
      dplyr::distinct(!!rlang::ensym(pres_y), !!rlang::ensym(pres_x), !!rlang::ensym(taxa_col))

    # Taxa bins ----
    taxa_bins <- taxa_data %>%
      dplyr::left_join(taxonomy$subspecies$lutaxa) %>%
      dplyr::left_join(taxonomy$species$lutaxa %>%
                         dplyr::rename(species=taxa) %>%
                         dplyr::distinct(original_name,species)
      ) %>%
      dplyr::distinct() %>%
      dplyr::filter(returned_rank <= settings$analysis_rank)

    # Taxa of interest ----

    aoi_taxa <-

    listed_taxa <- listed_df %>%
      dplyr::rename_with(~ gsub("rating","status",.x),contains("rating")) %>%
      dplyr::select(-tidyr::any_of("rank"))

    taxa_int <- taxa_bins %>%
      dplyr::inner_join(aoi_taxa %>%
                          dplyr::left_join(taxa$species$lutaxa) %>%
                          dplyr::select(original_name
                                        ,species=taxa
                          )
      ) %>%
      dplyr::inner_join(listed_taxa) %>%
      dplyr::distinct(species,taxa,original_name,long,lat)

    # Distance of closest record to study area ----
    taxa_rec_dist <- taxa_fbd %>%
      tidyr::nest(.by = c(species,taxa,original_name), data = c(long,lat)) %>%
      dplyr::mutate(dist=furrr::future_map_dbl(data, ~ nngeo::st_nn(settings$aoi
                                                                    , .x %>%
                                                                      sf::st_as_sf(coords=c("long","lat"),crs=settings$latlon_epsg) %>%
                                                                      sf::st_transform(crs = sf::st_crs(settings$aoi))
                                                                    , k = 1
                                                                    , maxdist = Inf
                                                                    , returnDist = TRUE
                                                                    , progress = FALSE
                                                                    , parallel = 1
      ) %>%
        .$dist %>%
        unlist()
      , .options = furrr::furrr_options(globals = "settings"
                                        , seed = TRUE
      )
      )
      ) %>%
      dplyr::distinct(species,taxa,original_name,dist)

   # Distance of distribution to study area ----

    # Available distribution details for all taxa
    taxa_ds <- taxa_ds %>%
      dplyr::left_join(taxa$subspecies$lutaxa) %>%
      dplyr::filter(returned_rank <= settings$analysis_rank) %>%
      dplyr::distinct(taxa,ds,file)

    # Find distance to relevant distribution per taxa
    taxa_distrib_dist <- taxa_fbd %>%
      dplyr::distinct(species, taxa, original_name) %>%
      dplyr::inner_join(taxa_ds) %>%
      dplyr::mutate(dist = furrr::future_map_dbl(file
                                                 , \(x) nngeo::st_nn(settings$aoi
                                                                     , x %>%
                                                                       arrow::open_dataset() %>%
                                                                       sfarrow::read_sf_dataset() %>%
                                                                       sf::st_transform(crs = sf::st_crs(settings$aoi)) %>%
                                                                       sf::st_make_valid() %>%
                                                                       dplyr::summarise() %>%
                                                                       sf::st_make_valid()
                                                                     , k = 1
                                                                     , maxdist = Inf
                                                                     , returnDist = TRUE
                                                                     , progress = FALSE
                                                                     , parallel = 1
                                                 ) %>%
                                                   .$dist %>%
                                                   unlist()
                                                 , .options = furrr::furrr_options(globals = "settings"
                                                                                   , seed = TRUE
                                                 )
      )
      ) %>%
      dplyr::distinct(species,taxa,original_name,dist)


    listed_timer <- timer("distance to distribution"
                          , time_df = listed_timer
    )

    # Relevant listed taxa to aoi ----
    # Relevant taxa data to aoi (for use in multiple locations below)
    taxa_rel <- taxa_rec_dist %>%
      dplyr::bind_rows(taxa_distrib_dist) %>%
      dplyr::filter(dist <= settings$use_buffer) %>%
      dplyr::distinct(species,taxa,original_name)

    listed_timer <- timer("aoi relevant listed taxa"
                          , notes = paste0(length(unique(taxa_rel$taxa)), " aoi listed original names")
                          , time_df = listed_timer
    )

    # Add listed binomials ----
    # i.e. taxa listed at species level
    # To keep listed binomials with MCP overlap in aoi_taxa
    # Really need to calculate MCPs at ssp level to avoid this
    # If so, aoi_listed would be performing all the same functions of aoi_taxa at spp & ssp levels
    bi_listed <- aoi_taxa %>%
      dplyr::inner_join(listed_taxa) %>%
      dplyr::filter(rated_rank=="sp") %>%
      dplyr::select(tidyr::any_of(contains(c("taxa","original_name","status","listed"))))

    listed_timer <- timer("bi_listed"
                          , notes = paste0(length(unique(bi_listed$taxa)), " aoi listed binomials")
                          , time_df = listed_timer
    )

    # All relevant listed taxa with listing details ----
    aoi_listed <- taxa_rel %>%
      dplyr::left_join(listed_taxa) %>%
      dplyr::select(-c(taxa,rated_rank)) %>%
      dplyr::rename(taxa=species) %>%
      dplyr::bind_rows(bi_listed) %>%
      dplyr::arrange(taxa) %>%
      dplyr::distinct() %>%
      dplyr::filter(!taxa %in% settings$delete_taxa$aoi_listed) # temporary fix for error in future listed database

    saveRDS(aoi_listed, reg_taxa_file)

  } else {

    ss <- readRDS(reg_taxa_file)

  }

  return(ss)

}
