
#' Error check for loci in input data
#'
#' @param loci_pr User provided loci info.
#' @param loci_mix Loci info from mixture.
#' @param loci_base Loci info from baseline.
#'
#' @noRd
#'
check_loci <- function(loci_pr, loci_mix, loci_base) {

  if (!is.null(loci_pr)) { # check loci if provided

    if (!setequal(loci_mix, loci_pr)) {
      if (length(setdiff(loci_pr, loci_mix) > 0)) {
        return(c("These provided loci names are not in mixture: ",
                 paste0(setdiff(loci_pr, loci_mix), ", ")))
      } else {
        return(c("Unidentified loci in mixture (not used): ",
                 paste0(setdiff(loci_mix, loci_pr), ", ")))
      }
    }

    if (!setequal(loci_base, loci_pr)) {
      if (length(setdiff(loci_pr, loci_base) > 0)) {
        return(c("These provided loci names are not in baseline: ",
                 paste0(setdiff(loci_pr, loci_base), ", ")))
      } else {
        return(c("Unidentified loci in baseline (not used): ",
                 paste0(setdiff(loci_base, loci_pr), ", ")))
      }
    }

  } else if (!setequal(loci_base, loci_mix)) {

    if (length(setdiff(loci_base, loci_mix)) > 0) {
      return(c("These loci are only in baseline: ",
               paste0(setdiff(loci_base, loci_mix), ", ")))
    } else {
      return(c("These loci are only in mixture: ",
               paste0(setdiff(loci_mix, loci_base), ", ")))
    }

  }

  return("good")

}


#' Preparing MAGMA input data
#'
#' @param wd Directory where you set up the *data* folder.
#' @param age_classes Hard code class categories for group ages.
#' @param loci_names Optional. String containing loci names.
#' @param save_data Logical (with default = `TRUE`). Option to save the data in the *data* folder.
#'
#' @return A list objects as the input data for `msgsi_mdl()`
#'
#' @importFrom magrittr %>%
#' @importFrom doRNG %dorng%
#' @importFrom foreach %dopar%
#'
#' @examples
#' \dontrun{
#' wd <- getwd() # path to data folder
#' magma_data <- magmatize_data(wd = wd, save_data = FALSE)
#' }
#'
#' @export
magmatize_data <-
  function(wd, age_classes = "all", fishery = NULL, loci_names = NULL, save_data = TRUE) {

    start_time <- Sys.time()

    load(paste0(wd, "/data/mixture.RData"))
    load(paste0(wd, "/data/baseline.RData"))

    message("Compiling input data, may take a minute or two...")

    #### Allele frequency function #### ----

    allefreq <- function(gble_in, loci, alleles, collect_by = SillySource) {

      n_alleles = alleles %>%
        dplyr::group_by(locus) %>%
        dplyr::summarise(n_allele = max(as.numeric(altyp)), .groups = "drop")

      scores_cols = sapply(loci, function(locus) {
        c(locus, paste0(locus, ".1"))
      }) %>%
        as.vector()

      gble_in %>%
        dplyr::select(c({{ collect_by }}, dplyr::all_of(scores_cols))) %>%
        tidyr::pivot_longer(
          cols = -{{ collect_by }},
          names_to = "locus",
          values_to = "allele"
        ) %>%
        dplyr::mutate(
          locus = stringr::str_replace(string = locus, pattern = "\\.1$", replacement = "")
        ) %>%
        dplyr::left_join(alleles,
                         by = c("locus" = "locus", "allele" = "call"),
                         keep = FALSE) %>%
        dplyr::group_by({{ collect_by }}, locus) %>%
        dplyr::count(altyp, .drop = FALSE) %>%
        dplyr::filter(!is.na(altyp)) %>%
        dplyr::left_join(n_alleles,
                         by = c("locus" = "locus"),
                         keep = FALSE) %>%
        dplyr::filter(as.numeric(altyp) <= n_allele) %>%
        dplyr::select(-n_allele) %>%
        tidyr::unite("altyp", c(locus, altyp)) %>%
        tidyr::pivot_wider(names_from = altyp, values_from = n) %>%
        dplyr::ungroup()

    }


    #### Age Data #### ----

    metadat0 <-
      utils::read.table(
        file = paste0(wd, "/data/metadata.txt"),
        sep = "\t",
        header = TRUE,
        stringsAsFactors = FALSE,
        row.names = 1
      )

    euro_age <- stats::setNames(metadat0$AGE_EUROPEAN, rownames(metadat0))

    euro_age[!is.na(euro_age)] <-
      sapply(euro_age[!is.na(euro_age)], function(chr) {
        paste(c(rep(0, 2 - nchar(chr)), chr), collapse = "")
      }) # add 0's to 0x class

    fw_age_range <-
      range(as.integer(substr(
        euro_age[!is.na(euro_age)], start = 1, stop = 1
      )))
    sw_age_range <-
      range(as.integer(substr(
        euro_age[!is.na(euro_age)], start = 2, stop = 2
      )))
    fw_ages <- seq.int(fw_age_range[1], fw_age_range[2])
    sw_ages <- seq.int(sw_age_range[1], sw_age_range[2])

    euro_ages <-
      apply(expand.grid(fw_ages, sw_ages), 1 , paste, collapse = "")

    extra_ages <-
      age_classes[which(!age_classes %in% c(euro_ages, "0X", "other", "all"))]

    euro_ages <-
      sort(c(extra_ages, euro_ages))

    C <- length(euro_ages)

    if ("all" %in% age_classes) {
      age_classes <- euro_ages
    } else if (any(!stats::na.omit(euro_age) %in% age_classes) & !"other" %in% age_classes) {
      stop("Unspecified age class(es) found in metadata. Include an 'other' age class to catch unspecified age classes.")
    } else if (!all(c("0X", "other") %in% age_classes) & length(age_classes) > length(euro_ages)) {
      euro_ages <- age_classes
      C <- length(euro_ages)
    }

    A <- length(age_classes)

    age_class <- stats::setNames(rep(A, C), euro_ages)

    if ('0X' %in% age_classes) {
      age_class[substring(euro_ages, 1, 1) == "0"] <- 1L
      if (!any(substring(euro_ages, 1, 1) == "0")) warning("You have 0X age class, but no 0 age detected in metadata.")
    } # for instances when we do not want to hard code age 0 composite class.

    age_class[euro_ages %in% age_classes] <-
      match(euro_ages[euro_ages %in% age_classes], age_classes)

    #### Baseline Data #### ----

    groups0 <-
      utils::read.table(
        file = paste0(wd, "/data/groups", fishery, ".txt"),
        header = TRUE,
        row.names = 1,
        sep = "\t",
        stringsAsFactors = FALSE
      )

    group_names <-
      utils::read.table(
        file = paste0(wd, "/data/group_names", fishery, ".txt"),
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE
      )

    hatcheries <- sort(unique(metadat0$SOURCE[metadat0$SOURCE != "WILD"]))

    # identify hatcheries show up in metadat but missing from groups
    miss_hatch_grp <- hatcheries[!hatcheries %in% rownames(groups0)]

    G <- apply(group_names, 2, function(g) {
      ifelse("Other" %in% g, which(g == "Other"), 0)
    }) # group number for "Other" group

    if (sum(G) == 0 & length(miss_hatch_grp) > 0) {
      stop(paste("unidentified hatchery groups:", miss_hatch_grp, sep = " "))
    }

    groups <-
      rbind(groups0,
            matrix(
              rep(G, length(miss_hatch_grp)),
              nrow = length(miss_hatch_grp),
              ncol = length(G),
              byrow = TRUE,
              dimnames = list(miss_hatch_grp, names(groups0))
            ) ) # append missing hatcheries

    wildpops <-
      rownames(groups)[grepl("([0-9]+).*$", rownames(groups)) & sapply(rownames(groups), function(gn) nchar(gn) > 4)] # look for pops with numbers (sample year) and > 4 characters, not sure this works all the time?

    # identify hatcheries in groups that didn't show up in metadat (sample)
    miss_hatch_meta <- rownames(groups)[!rownames(groups) %in% c(wildpops, hatcheries)]

    hatcheries <- c(hatcheries, miss_hatch_meta)

    all_pops <- c(wildpops, hatcheries)

    groups <- groups[all_pops, , drop = FALSE]
    # make sure groups and baseline are in the same row order
    # force groups as an array, otherwise groups is converted to a vector
    # under single district scenarios, crashing downstream analyses.

    base_data <-
      sapply(wildpops, function(pop) {
        get(paste0(pop, ".gcl"), pos = -1)
      }, simplify = FALSE) %>%
      dplyr::bind_rows()

    # loci check

    mix_sillys <-
      unique(sapply(stringr::str_split(rownames(metadat0), "_"), "[", 1))

    mixture_data <-
      sapply(mix_sillys, function(silly) {
        get(paste0(silly, ".gcl"), pos = -1)
      }, simplify = FALSE) %>%
      dplyr::bind_rows() # all gcl objects into a list, for parallel process

    loci_base <-
      dplyr::tibble(locus = names(base_data)) %>%
      dplyr::filter(grepl("\\.1$", locus)) %>%
      dplyr::mutate(locus = substr(locus, 1, nchar(locus) - 2)) %>%
      dplyr::pull(locus)

    loci_mix <-
      dplyr::tibble(locus = names(mixture_data)) %>%
      dplyr::filter(grepl("\\.1$", locus)) %>%
      dplyr::mutate(locus = substr(locus, 1, nchar(locus) - 2)) %>%
      dplyr::pull(locus)

    error_message <- check_loci(loci_names, loci_mix, loci_base)

    if (grepl("Unidentified", error_message[1])) {
      warning(error_message)
    } else if (!"good" %in% error_message) {
      stop(error_message)
    }

    if (is.null(loci_names)) {
      loci <- loci_mix
    } else {
      loci <- loci_names
    }

    alleles_tib <- lapply(loci, function(loc) {
      dplyr::tibble(locus = loc,
                    call = base_data %>%
                      dplyr::select(dplyr::all_of(loc), paste0(loc, ".1")) %>%
                      unlist() %>% unique() %>% .[!is.na(.)],
                    altyp = seq.int(dplyr::n_distinct(call)) %>% factor())
    }) %>% dplyr::bind_rows()

    nalleles <- alleles_tib %>%
      dplyr::group_by(locus) %>%
      dplyr::summarise(n_allele = max(as.numeric(altyp)), .groups = "drop") %>%
      tibble::deframe()

    y0 <- allefreq(base_data, loci, alleles_tib, collect_by = SILLY_CODE)

    K <- length(wildpops)

    H <- length(hatcheries)

    age <- as.integer(factor(euro_age[rownames(metadat0)], levels = euro_ages))

    traits <- c(loci, "age")
    nstates <-
      c(nalleles[loci], stats::setNames(C, "age"))
    trait_fac <-
      factor(rep(traits, nstates), levels = traits)
    states_of_traits <-
      paste(trait_fac, unlist(lapply(nstates, seq)), sep = "_")
    alleles <- states_of_traits[-seq(sum(nalleles) + 1, sum(nstates))]

    y <-
      matrix(
        0L,
        nrow = K + H,
        ncol = sum(nstates),
        dimnames = list(all_pops, states_of_traits)
      )

    y[wildpops, alleles] <-
      as.matrix(
        y0[match(wildpops, y0$SILLY_CODE),
           match(colnames(y)[seq(sum(nalleles))], colnames(y0))]
      )


    #### Mixture #### ----

    if (length(mix_sillys) > 1) {

      ncores <- min(length(mix_sillys), 4)
      cl <- parallel::makePSOCKcluster(ncores)
      doParallel::registerDoParallel(cl, cores = ncores)

      x0 <- foreach::foreach(mix = mix_sillys, .packages = c("magrittr", "dplyr")) %dopar% {
        do_gcl <- dplyr::filter(mixture_data, SILLY_CODE == mix)
        allefreq(do_gcl, loci, alleles_tib)

      } %>% dplyr::bind_rows()

      parallel::stopCluster(cl)

    } else {

      x0 <- allefreq(mixture_data, loci, alleles_tib)

    }

    x <-
      matrix(
        0L,
        nrow = nrow(metadat0),
        ncol = sum(nstates),
        dimnames = list(rownames(metadat0), states_of_traits)
        # ncol = sum(nalleles),
        # dimnames = list(rownames(metadat0), states_of_traits[seq(sum(nalleles))])
      ) # x rows are in the same order as metadat
    x[x0$SillySource, seq(sum(nalleles))] <-
      as.matrix(x0[, match(colnames(x)[seq(sum(nalleles))], colnames(x0))])
    x[!is.na(age), seq(sum(nalleles) + 1, sum(nstates))] <-
      t(sapply(seq(sum(!is.na(age))),
               function(mm) {
                 ag = age[!is.na(age)][mm]
                 ags = rep(0L, C)
                 ags[+ag] = 1L
                 ags
               }))
    # x <- x %>%
    #   dplyr::bind_cols(
    #     {t(sapply(age, function(ag) {
    #       if (is.na(ag)) rep(0, C)
    #       else c(rep(0, ag - 1), 1, rep(0, C - ag))
    #     })) %>%
    #         dplyr::as_tibble(.name_repair = "minimal") %>%
    #         stats::setNames(paste0("age_", seq(C)))}
    #     )


    #### metadata stuff #### ----

    i <- as.integer(factor(metadat0$SOURCE, levels = hatcheries)) + K

    districts <-
      stats::setNames(as.character(sort(unique(
        metadat0$DISTRICT
      ))), sort(unique(metadat0$DISTRICT)))
    D <- length(districts)
    district <-
      as.integer(factor(metadat0$DISTRICT, levels = districts))

    subdistricts <-
      lapply(districts, function(d) {
        stats::setNames(sort(unique(
          subset(metadat0, DISTRICT == d)$SUBDISTRICT
        )), sort(unique(
          subset(metadat0, DISTRICT == d)$SUBDISTRICT
        )))
      })
    S = sapply(subdistricts, length)
    subdistrict <-
      as.integer(unlist(lapply(districts, function(d) {
        as.integer(
          factor(metadat0$SUBDISTRICT[metadat0$DISTRICT == d],
                 levels = subdistricts[[d]])
        )
      })))

    stat_weeks <- sort(unique(metadat0$STAT_WEEK))
    W <- length(stat_weeks)
    stat_week <-
      as.integer(factor(metadat0$STAT_WEEK, levels = stat_weeks))

    metadat <-
      data.frame(
        district = district,
        subdist = subdistrict,
        week = stat_week,
        iden = i,
        age = age,
        stringsAsFactors = FALSE,
        row.names = rownames(metadat0)
      )

    harvest <-
      utils::read.table( paste0(wd, "/data/harvest", fishery, ".txt"),
                  header = TRUE,
                  sep = "\t",
                  stringsAsFactors = FALSE
      )
    harvest$DISTRICT <-
      as.integer(factor(harvest$DISTRICT, levels = districts))
    harvest$SUBDISTRICT <-
      sapply(seq(nrow(harvest)), function(rw) {
        as.integer(
          factor(harvest$SUBDISTRICT[rw],
                 levels = subdistricts[[harvest$DISTRICT[rw]]])
        )
      })
    harvest$STAT_WEEK <-
      as.integer(factor(harvest$STAT_WEEK, levels = stat_weeks))
    harvest$HARVEST <- as.numeric(harvest$HARVEST)


    #### save data and format output #### ----

    dat_out <- list(
      x = x,
      y = y,
      metadat = metadat,
      harvest = harvest,
      nstates = nstates,
      nalleles = nalleles,
      C = C,
      groups = groups,
      group_names = group_names,
      age_class = age_class,
      age_classes = age_classes,
      wildpops = wildpops,
      hatcheries = hatcheries,
      districts = districts,
      subdistricts = subdistricts,
      stat_weeks = stat_weeks
    )

    # if (save_data) {
    #   magma_data_names <-
    #     c("y", # baseline of allele freq
    #       "x", # mixture data
    #       "metadat", # metadata, where, when, age, otolith assignment (i)
    #       "harvest", # harvest
    #       "groups", # groupvec with all hatcheries
    #       "group_names", # all reporting groups
    #       "K", # wild pops
    #       "H", # hatchery pops
    #       "A", # age classes
    #       "C", # total possible age classes
    #       # "Yr", # year
    #       "D", # District
    #       "S", # Subdistrict
    #       "W", # weeks
    #       "age_class", # length = C; all possible
    #       "age_classes", # length = A; age reporting groups
    #       # "years2run", # year
    #       "districts", # list of district names
    #       "subdistricts", # subdistricts for each district
    #       "stat_weeks", # which stat weeks
    #       "loci", # list of loci
    #       "nalleles", # number of possible alleles for each loci
    #       "nstates", # number of possible alleles for each loci + number of possible ages
    #       "hatcheries", # Hatcheries
    #       "wildpops", # Wild populations
    #       "magma_data_names",
    #       "dat_out")
    #
    #   save(list = magma_data_names,
    #        file = paste0(wd, "/data/magma_data", fishery, ".RData"))
    # }

    if(isTRUE(save_data)) {
      saveRDS(dat_out, file = paste0(wd, "/data/magma_data", fishery, ".Rds"))
    }

    if(length(miss_hatch_grp) == 0) {
      message("No missing hatcheries")
    } else {
      message(
        paste("Missing hatcheries found and added to groups:", miss_hatch_grp, sep = " ")
      )
    } # added catch to make aware of any missing hatcheries

    print(Sys.time() - start_time)

    return(dat_out)

  }


#' Prepare malia data
#'
#' Take a magma data object and save the items as text files in a designated directory.
#'
#' @param magma_data Magma data object made using `magmatize_data()` function.
#' @param path A designated directory where you set up for malia data set.
#'
#' @return Items in magma data object saved as text files
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' prep_malia_data(magma_data, "D:/bobby_adfg/projects/magma/malia/data/bb2022")
#' }
#'
#' @export
prep_malia_data <- function(magma_data, path) {

  magma_data$metadat %>%
    tibble::rownames_to_column(var = "name") %>%
    tibble::as_tibble() %>%
    tidyr::replace_na(list(iden=0, age=0)) %>%
    as.matrix %>%
    trimws %>%
    utils::write.table(., paste0(path, "/metadat.txt"), row.names = FALSE)

  readr::write_lines(magma_data$C, paste0(path, "/c.txt"))

  tibble::enframe(magma_data$nstates) %>%
    dplyr::rename(n_state = value) %>%
    as.matrix %>%
    trimws %>%
    utils::write.table(., paste0(path, "/nstates.txt"), row.names = FALSE)

  tibble::enframe(magma_data$nalleles) %>%
    dplyr::rename(n_allele = value) %>%
    as.matrix %>%
    trimws %>%
    utils::write.table(., paste0(path, "/nalleles.txt"), row.names = FALSE)

  readr::write_lines(magma_data$wildpops, paste0(path, "/wildpops.txt"))

  readr::write_lines(magma_data$hatcheries, paste0(path, "/hatcheries.txt"))

  magma_data$groups %>%
    tibble::rownames_to_column(var = "grp_name") %>%
    as.matrix %>%
    trimws %>%
    utils::write.table(., paste0(path, "/groups.txt"), row.names = FALSE)

  tibble::enframe(magma_data$age_class) %>%
    dplyr::rename(class = value) %>%
    as.matrix %>%
    trimws %>%
    utils::write.table(., paste0(path, "/age_class.txt"), row.names = FALSE)

  utils::write.table(magma_data$x, paste0(path, "/x.txt"), row.names = FALSE)
  utils::write.table(magma_data$y, paste0(path, "/y.txt"), row.names = FALSE)

}


utils::globalVariables(c("SillySource", "locus", "n", "altyp", "n_allele",
                         "SILLY_CODE", "mix", ".",#"fishery",
                         "collection", "DISTRICT"))















