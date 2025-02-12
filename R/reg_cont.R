
#' Regional contribution to national range (AOO and EOO)
#'
#' @param taxa Character. Taxa name for use in ConR and outputs.
#' @param presence Cleaned and filtered dataframe of presences.
#' @param out_dir Character. Directory path for the results files to be saved (summary_stats.rds and mcp_noreg.parquet).
#' mcp_noreg will be saved by `sfarrow::st_write_parquet()`. Currently will not work very well with any
#' full stop in the path. Other file types are changed to .parquet.
#' @param mcp_file Character. Path for the national mcp. If the file exists it will be read in,
#'  otherwise it will be created using make_mcp.
#' @param grd_file Character. Path for the national grid. If the file exists it will be read in,
#'  otherwise it will be created using make_grd.
#' @param region_bound sf object. Region boundary polygon.
#' @param force_new Logical. If any of the results files in the `out_dir` exists, recreate it?
#' @param pres_x,pres_y Character. Name of the columns in `presence` that have
#' the x and y coordinates.
#' @param in_crs epsg code for coordinates in `presence`.
#' @param out_crs epsg code for coordinates in output grid.
#' @param nearest_pop Logical. Calculate distance to nearest extra-regional population?
#' @param use_ConR,cell_size,num_rast_pos Passed to make_grd if national grid doesn't exist or force_new = TRUE.
#' @param clip Passed to make_mcp if national mcp doesn't exist or force_new = TRUE.
#'
#' @return summary_stats.rds with aoo and eoo stats,
#' and mcp_noreg.parquet with the mcp excluding regional presences.
#' @export
#'

reg_cont <- function(taxa
                     , presence
                     , out_dir
                     , mcp_file = NULL
                     , grd_file = NULL
                     , region_bound
                     , force_new = FALSE
                     , pres_x = "long"
                     , pres_y = "lat"
                     , in_crs = 4326
                     , out_crs = in_crs
                     , nearest_pop = FALSE
                     , use_ConR = TRUE
                     , cell_size = 2
                     , num_rast_pos = 0
                     , clip = NULL
) {

  sum_stats_file <- fs::path(out_dir,"summary_stats.rds")

  run <- if(file.exists(sum_stats_file)) force_new else TRUE

  if(run) {

    # Base dataset ----
    # Dataset for generating both grids (aoo) and mcps (eoo)
    df <- presence %>%
      dplyr::distinct(!!rlang::ensym(pres_y), !!rlang::ensym(pres_x)) %>%
      sf::st_as_sf(coords = c(pres_x, pres_y)
                   , crs = in_crs
                   , remove = FALSE
      ) %>%
      sf::st_join(region_bound %>%
                    sf::st_transform(crs = in_crs) %>%
                    dplyr::mutate(region="inreg")
      ) %>%
      dplyr::mutate(region=ifelse(!is.na(region),"inreg","noreg")) %>% # for generating in and out of region stats
      sf::st_make_valid()

    # output dir ----
    fs::dir_create(out_dir)

    # AOO ----

    ## load/create national grid ----
    if(file.exists(grd_file) & !force_new) {

      grd <- sfarrow::st_read_parquet(grd_file) %>%
        sf::st_transform(crs=out_crs)

    } else {

      grd <- make_grd(presence = df %>%
                        dplyr::mutate(taxa = taxa)
                      , out_file = grd_file
                      , force_new = FALSE
                      , in_crs = in_crs
                      , out_crs = out_crs
                      , use_ConR = use_ConR
                      , cell_size = cell_size
                      , num_rast_pos = num_rast_pos
      )

    }


    ## aoo summary stats ----
    aoo <- grd %>%
      sf::st_transform(crs = out_crs) %>%
      dplyr::mutate(cell=dplyr::row_number()) %>%
      sf::st_join(region_bound %>%
                    sf::st_transform(crs = out_crs) %>%
                    dplyr::mutate(region="inreg")
      ) %>%
      dplyr::mutate(region=ifelse(!is.na(region),"inreg","noreg")) %>% # for generating in and out of region stats
      sf::st_set_geometry(NULL) %>%
      dplyr::group_by(region) %>%
      dplyr::summarise(cells=dplyr::n_distinct(cell),
                       AOO=cells*(cell_size^2) # instead of rounding error produced by st_area
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(-cells) %>%
      tidyr::pivot_wider(names_from = "region",names_prefix = "AOO_",values_from = "AOO") %>%
      {if(!"AOO_noreg" %in% names(.)) dplyr::mutate(.,AOO_noreg=0) else .} %>%
      {if(!"AOO_inreg" %in% names(.)) dplyr::mutate(.,AOO_inreg=0) else .} %>%
      dplyr::mutate(AOO_tot=AOO_inreg+AOO_noreg,
                    AOO_regpc=AOO_inreg/AOO_tot*100
      ) %>%
      dplyr::relocate(AOO_tot,AOO_noreg,AOO_inreg,AOO_regpc)


    # EOO ----

    ## load/create national mcp ----
    if(file.exists(mcp_file) & !force_new) {

      mcp <- sfarrow::st_read_parquet(mcp_file) %>%
        sf::st_transform(crs=out_crs) %>%
        dplyr::mutate(EOO=as.numeric(sf::st_area(.)/1e+6), # area in square km
                      EOO=round(EOO,2),
                      type="tot",
                      EOO_label=paste0(as.character(EOO)," km","\U00B2")
        )

    } else if(nrow(df)>=3 & (!file.exists(mcp_file)|force_new)) {

      mcp <- make_mcp(presence = df
                      , out_file = mcp_file
                      , force_new = FALSE
                      , in_crs = in_crs
                      , out_crs = out_crs
                      , buf = 0
                      , clip = clip
      ) %>%
        dplyr::mutate(EOO=as.numeric(sf::st_area(.)/1e+6), # area in square km
                      EOO=round(EOO,2),
                      type="tot",
                      EOO_label=paste0(as.character(EOO)," km","\U00B2")
        )

    } else {

      mcp <- tibble::tibble(EOO=NA,
                            type="tot",
                            EOO_label=NA
      )

    }

    ## national mcp without region ----

    mcp_noreg_file <- fs::path(out_dir,"mcp_noreg.parquet")

    if(file.exists(mcp_noreg_file) & !force_new) {

      mcp_no_region <- sfarrow::st_read_parquet(mcp_noreg_file) %>%
        sf::st_transform(crs=out_crs) %>%
        dplyr::mutate(EOO=as.numeric(sf::st_area(.)/1e+6),
                      EOO=round(EOO,2),
                      type="noreg"
        )

    } else if(nrow(dplyr::filter(df,region=="noreg"))>=3) {

      mcp_no_region <- df %>%
        sf::st_transform(crs=out_crs) %>%
        dplyr::filter(region=="noreg") %>%
        make_mcp(out_file = mcp_noreg_file
                 , force_new = force_new
                 , in_crs = in_crs
                 , out_crs = out_crs
                 , buf = 0
                 , clip = clip
        )

      if(isTRUE(nrow(mcp_no_region)>0)) { # catch where make_mcp clip removes whole polygon and returns NA

        mcp_no_region <- mcp_no_region %>%
          dplyr::mutate(EOO=as.numeric(sf::st_area(.)/1e+6),
                        EOO=round(EOO,2),
                        type="noreg"
          )

      } else {

        mcp_no_region <- tibble::tibble(EOO=NA,
                                        type="noreg"
        )

      }


    } else {

      mcp_no_region <- tibble::tibble(EOO=NA,
                                      type="noreg"
      )

    }

    ## eoo summary stats ----
    eoo <- mcp %>%
      dplyr::select(-EOO_label) %>%
      dplyr::bind_rows(mcp_no_region) %>%
      {if("geometry" %in% names(.)) sf::st_set_geometry(.,NULL) else .} %>%
      tidyr::pivot_wider(names_from = "type",names_prefix = "EOO_",values_from = "EOO") %>%
      dplyr::mutate(EOO_inreg=EOO_tot-EOO_noreg, # region contribution to national EOO
                    EOO_regpc=EOO_inreg/EOO_tot*100
      )

    # Distance to nearest extra-regional population ----
    if(nearest_pop) {

      if(nrow(dplyr::filter(df,region=="noreg"))>0 & nrow(dplyr::filter(df,region=="inreg"))>0) {

        near_pop <- df %>%
          dplyr::filter(region=="inreg") %>%
          dplyr::summarise() %>%
          sf::st_make_valid() %>%
          nngeo::st_nn(df %>% dplyr::filter(region=="noreg"),
                       sparse = TRUE,
                       k = 1,
                       maxdist = Inf,
                       returnDist = TRUE,
                       progress = FALSE,
                       parallel = 1
          ) %>%
          .$dist %>%
          unlist()

      } else {

        near_pop=NA

      }

    }

    # Summary stats ----
    ss <- aoo %>%
      dplyr::bind_cols(eoo) %>%
      dplyr::mutate(reg_rec=nrow(dplyr::filter(df,region=="inreg"))
                    , taxa = taxa
                    ) %>%
      {if(nearest_pop) dplyr::mutate(.,dist_near_pop=near_pop/1000) else .} %>% # convert distance in metres from st_nn to km
      dplyr::relocate(taxa)

    saveRDS(ss, sum_stats_file)

  } else {

    ss <- readRDS(sum_stats_file)

  }

  return(ss)

}
