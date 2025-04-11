#' Find file path to relevant taxa distribution (geographic range) parquet files.
#'
#' For use inside filter_by_distribution and fix_spatial_taxonomy functions,
#' and for any other situation needing to return relevant distribution file paths per taxa (e.g. for aoi-taxa and dist in envPIA).
#'
#' @param distrib_dir Directory path. Directory containing distribution layers.
#' @param sources Character vector. Distribution sources to use in the distrib_dir, i.e. EPBC, Expert, Birds SA, Redlist, or other.
#' @param source_scales Character vector. Geographic scale of the distribution sources, e.g. state or national. Needs to be the same length as sources,
#' and the position needs to match the position of the source it is referring to. For example, with sources = c("epbc","expert"),
#' the corresponding scale_sources would be c("national","state") where epbc is national and expert is state.
#' Used where different distributions are required at different scales, e.g. state and national distributions both required in filter_by_distribution.
#' @param standardise_taxonomy Logical. Does the taxonomy need to be standardised using a taxonomy object from envClean::make_taxonomy?
#' @param target_ranks Character vector. Target taxonomic ranks for taxonomy standardisation. Either 'species' or 'subspecies', or both.
#' @param taxonomy Taxonomy object returned by envClean::make_taxonomy to use for taxonomy standardisation.
#' @param rm_ssp_mismatches Logical. Remove subspecies mismatches, i.e. where original name was subspecies
#' but was matched to species level by envClean::make_taxonomy.
#' Setting to TRUE is useful where the result is input into filter_by_distribution,
#' as it avoids potentially erroneous filtering of binomial records by trinomial distributions.
#' The downside of setting to TRUE is some relevant distributions may be lost where the match of ssp to spp is correct,
#' e.g. a subspecies is now a new species. Ideally, fixes should be provided for these cases in envClean::make_taxonomy
#' to create the taxonomy object and this parameter set to FALSE.
#' @param source_rank Logical. Rank the distribution sources and filter top ranked? Ranking will be in the order given by sources vector.
#' If !is.null(source_scales) & source_rank == TRUE, then source ranking will done in source_scale groups,
#' i.e. the highest ranked source per scale will be retained.
#' Also, if standardise_taxonomy == TRUE, then ranking will be according to the standardised taxonomy.
#' @param datatype Character. Either 'vector' (default) or 'raster'.
#'
#' @return Data frame with relevant distribution sources (ds) and file path (file) per taxa (original_name).
#' If !is.null(source_scales), there will be a file column per scale (e.g. national, state), instead of the single source (ds) column.
#' If standardise_taxonomy == TRUE, the output taxa column will be 'taxa' and not 'original_name'.
#'
#' @export
#'

dists_source <- function(distrib_dir = fs::path("H:","data"),
                         sources = c("other","epbc","expert","bsa","redlist"),
                         source_scales = NULL,
                         standardise_taxonomy = FALSE,
                         target_ranks = c("species", "subspecies"),
                         taxonomy = NULL,
                         rm_ssp_mismatches = FALSE,
                         source_rank = FALSE,
                         datatype = "vector"
){


  # Find all sources for all original_name ----
  taxa_ds <- sources %>%
    purrr::set_names() %>%
    purrr::map( \(x)

                 if(datatype=="vector"){

                   ds_path <- fs::path(distrib_dir, datatype, "distribution") %>%
                     fs::dir_info(type = "directory", regexp = x) |>
                     dplyr::filter(basename(path) == x) |> # remove any old variants or mismatches
                     dplyr::pull(path)

                   if(any(fs::dir_exists(ds_path))){ # ensures no erroneous entries where source dir/files do not exist (e.g. for sources other than epbc and expert in raster)

                     ds <- ds_path %>%
                       fs::dir_ls(type = "file", regexp = "\\.parquet$", recurse = TRUE) %>%
                       tibble::as_tibble() %>%
                       dplyr::rename(file = value) %>%
                       dplyr::mutate(original_name = gsub(paste0("^.*", x, "/(.*)/part-.*"), "\\1", file)
                                     , original_name = gsub("_", " ", original_name)
                                     , original_name = stringr::str_squish(original_name)
                                     , ds = x
                       ) %>%
                       dplyr::relocate(original_name, ds, file)

                   } else {

                     ds <- tibble::tibble()

                   }

                 } else {

                   ds_path <- fs::path(distrib_dir,datatype,"distribution","sa_ibrasub_xn____0__90",x)

                   if(fs::dir_exists(ds_path)){ # ensures no erroneous entries where source dir/files do not exist (e.g. for sources other than epbc and expert in raster)

                     ds <- ds_path %>%
                       fs::dir_ls(type = "file", regexp="\\.tif$", recurse = TRUE) %>%
                       tibble::as_tibble() %>%
                       dplyr::rename(file = value) %>%
                       dplyr::mutate(original_name = gsub(paste0("^.*/",x,"/(.*)\\.tif"), "\\1", file),
                                     original_name = stringr::str_squish(original_name),
                                     ds = x
                       ) %>%
                       dplyr::relocate(original_name, ds, file)

                   } else {

                     ds <- tibble::tibble()

                   }

                 }

    ) %>%
    purrr::list_rbind()

  # Add source scale ----
  if(!is.null(source_scales)){

    src_scale_lu <- tibble::tibble(ds = sources
                                   , scale = source_scales
    )

    taxa_ds <- taxa_ds %>%
      dplyr::left_join(src_scale_lu, by = "ds") %>%
      dplyr::arrange(original_name)

  }

  # Standardise taxonomy ----
  # need to do before source rank to obtain rankings per standardised taxa (as standardising groups some taxa)
  if(standardise_taxonomy){

    if(!is.null(taxonomy)){

      taxa_ds <- taxa_ds %>%
        {if("subspecies" %in% target_ranks) dplyr::left_join(.,taxonomy$subspecies$lutaxa, by="original_name")
          else dplyr::left_join(.,taxonomy$species$lutaxa, by="original_name")
        } %>%
        dplyr::filter(returned_rank %in% target_ranks) %>%
        {if(rm_ssp_mismatches) dplyr::filter(.,!(original_is_tri & returned_rank == "species")) else .} %>%
        dplyr::select(tidyr::any_of(c("taxa","ds","file","scale")))

    } else {

      warning("standardise_taxonomy == TRUE but no taxonomy object provided. Taxonomy will not be standardised.")

    }

  }

  # taxa_col for following steps
  if(standardise_taxonomy) taxa_col <- "taxa" else taxa_col <- "original_name"

  # Find most authoritative/detailed distribution per original_name ----
  if(source_rank){

    taxa_ds <- taxa_ds %>%
      dplyr::mutate(ds = factor(ds, levels = sources, ordered = TRUE)) %>%
      {if(!is.null(source_scales)) dplyr::group_by(.,dplyr::pick(tidyr::all_of(c(taxa_col, "scale")))) else dplyr::group_by(.,dplyr::pick(tidyr::all_of(taxa_col)))} %>%
      dplyr::filter(ds == min(ds)) %>%
      dplyr::ungroup() %>%
      dplyr::select(tidyr::any_of(c(taxa_col, "ds", "file","scale"))) %>%
      dplyr::arrange(taxa_col)

  }

  # One file column per scale ----
  if(!is.null(source_scales)){

    taxa_ds <- taxa_ds %>%
      tidyr::pivot_wider(id_cols = taxa_col, names_from = scale, values_from = file, values_fill = NULL, values_fn = list) %>%
      dplyr::mutate(across(tidyr::all_of(unique(source_scales)), ~replace(., lengths(.) == 0, NA)))

  }

  return(taxa_ds)

}
