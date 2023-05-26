
#' TBR output format wrapper step 2
#'
#' @param out1 Output from step 1
#' @param summ_level Summarize at district or subdistrict level
#' @param type Identify "pop" or "age" to summarize only populations or age class.
#'
#' @return Model output in process and metadata as a list object.
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' # format data
#' wd <- "D:/bobby_adfg/projects/magma/test_TBR" # path to data folder
#' magma_data <- magmatize_data(wd = wd, save_data = FALSE)
#'
#' # model run
#' magma_out <- magmatize_mdl(magma_data,
#'   nreps = 50, nburn = 25, thin = 1, nchains = 3)
#'
#' # summary steps 1 and 2
#' tbr1 <- magmatize_summ_tbr1(which_dist = 3,
#'   ma_out = magma_out,
#'   ma_dat = magma_data)
#'
#' tbr2 <- magmatize_summ_tbr2(tbr1, summ_level = "district", type = "pop")
#'
magmatize_summ_tbr2 <- function(out1, summ_level, type) {

  nreps <- out1$sub_dat$nreps
  nburn <- out1$sub_dat$nburn
  thin <- out1$sub_dat$thin
  nchains <- out1$sub_dat$nchains
  keep_burn <- out1$sub_dat$keep_burn
  summ_level -> out1$sub_dat$summ_level
  type -> out1$sub_dat$type

  C <- out1$sub_dat$C # number of age classes
  D <- dplyr::n_distinct(out1$sub_dat$metadat$district) # number of districts
  S <- out1$sub_dat$metadat %>%
    dplyr::group_by(district) %>%
    dplyr::summarise(S = dplyr::n_distinct(subdist), .groups = "drop") %>%
    dplyr::pull(S) # number of subdistricts
  W <- max(dplyr::n_distinct(out1$sub_dat$metadat$week), length(out1$sub_dat$stat_weeks)) # number of stats weeks

  if (type == "pop" & summ_level == "district") {
    out_a <- format_district_pop_a(out1$holder, out1$sub_dat, nreps, nburn, thin, nchains, keep_burn, C, D, S, W)
  } else if (type == "pop" & summ_level == "subdistrict") {
    out_a <- format_subdistrict_pop_a(out1$holder, out1$sub_dat, nreps, nburn, thin, nchains, keep_burn, C, D, S, W)
  } else if (type == "age" & summ_level == "district") {
    out_a <- format_district_age_a(out1$holder, out1$sub_dat, nreps, nburn, thin, nchains, keep_burn, C, D, S, W)
  } else if (type == "age" & summ_level == "subdistrict") {
    out_a <- format_subdistrict_age_a(out1$holder, out1$sub_dat, nreps, nburn, thin, nchains, keep_burn, C, D, S, W)
  }

  out <- list(holder = out_a,
              sub_dat = out1$sub_dat)

  return(out)

}

# pop ----
# Summarize population part of MAGMA output (part a)

#' Format district populations step a
#'
#' @param outraw A subset of MAGMA output.
#' @param dat_in A subset of MAGMA data.
#' @param nreps The same as *nreps* in MAGMA model run.
#' @param nburn The same as *nburn* in MAGMA model run.
#' @param thin The same as *thin* in MAGMA model run.
#' @param nchains The same as *nchains* in MAGMA model run.
#' @param keep_burn The same as *keep_burn* in MAGMA model run.
#' @param C Number of age classes.
#' @param D Number of districts.
#' @param S Number of subdistricts.
#' @param W Number of weeks.
#'
#' @importFrom magrittr %>%
#'
#' @noRd
#'
format_district_pop_a <- function(outraw, dat_in, nreps, nburn, thin, nchains, keep_burn, C, D, S, W) {

  ### data input ### ----
  metadat <- dat_in$metadat # age and iden info
  harvest <- dat_in$harvest %>% # for calculating p0
    dplyr::mutate(n = apply(.[, 2:4], 1, function(aa) {
      dplyr::filter(metadat,
                    district == aa[1],
                    subdist == aa[2],
                    week == aa[3]) %>% nrow()
    }) ) %>%
    dplyr::mutate(HARVEST = HARVEST * (n != 0)) %>%
    dplyr::select(-n)

  groups <- dat_in$groups # vector id for reporting groups (aka groupvec)
  group_names <- dat_in$group_names # reporting groups

  if (isFALSE(keep_burn)|keep_burn == "false") keep_burn = FALSE

  ### prepare output ### ----
  # organization of out_list:
  # [[chain]][[yr]][[dist]][[sub]][[week]][age, pop, itr]

  ap_prop <- ap_prop2 <- ap_prop2b <- ap_prop2c <- ap_combo_subdis <- ap_combo_dis <- list()
  ap_prop_all_1 <- ap_prop_all_2 <- list()

  for (d_idx in 1:D) {
    ap_prop[[d_idx]] <- list()
    ap_prop2b[[d_idx]] <- list()
    ap_prop2c[[d_idx]] <- list()
    ap_prop_all_1[[d_idx]] <- ap_prop_all_2[[d_idx]] <- list()

    idx <- 0
    for (w_idx in 1:W) {
      ap_prop[[d_idx]][[w_idx]] <- list()
      ap_prop2b[[d_idx]][[w_idx]] <- list()
      ap_prop2c[[d_idx]][[w_idx]] <- list()

      for (s_idx in 1:S[d_idx]) {

        # go in each chain and grab age-pop output according to specified sampling period
        # new! ap_prop: [[d]][[s]][[w]][[chain]][age*itr, c(pop, agevec)]
        # pivot_longer %>% [pop*age*itr, c(itr, age, grp, value)]

        # ap_combo_dis or ap_combo_subdis
        # combine age classes and pop groups
        # weighted by harvest
        # this is done for each week in each district
        # new! [[chain]][[d]][[w]][group*itr, c(itr, grpvec, ages)]

        idx <- idx + 1
        ap_prop[[d_idx]][[w_idx]][[s_idx]] <-
          apply(outraw, MARGIN = 6,
                function(ol) ol[ , , w_idx, s_idx, d_idx] %>%
                  as.data.frame() %>%
                  tidyr::pivot_longer(cols = 1:(ncol(.)-2)) %>%
                  dplyr::mutate(
                    grpvec = rep(group_names[groups[, d_idx], d_idx],
                                 C* (nreps- nburn* isFALSE(keep_burn))/ thin)
                  ) %>%
                  dplyr::select(-name) %>%
                  dplyr::rename(itr = 1, agevec = 2))

        ap_prop2b[[d_idx]][[w_idx]][[s_idx]] <-
          lapply(ap_prop[[d_idx]][[w_idx]][[s_idx]],
                 function(ap) {
                   ap %>%
                     dplyr::mutate(value = value* {harvest %>%
                         dplyr::group_by(DISTRICT, STAT_WEEK) %>%
                         dplyr::mutate(prop_harv = prop.table(HARVEST)) %>%
                         tidyr::replace_na(list(prop_harv = 0)) %>%
                         dplyr::filter(DISTRICT == d_idx,
                                       SUBDISTRICT == s_idx,
                                       STAT_WEEK == w_idx) %>%
                         dplyr::pull(prop_harv) %>% max(., 0)})
                 })

        ap_prop_all_1[[d_idx]][[W*(s_idx-1)+w_idx]] <-
          lapply(ap_prop[[d_idx]][[w_idx]][[s_idx]],
                 function(apa) {
                   apa %>%
                     dplyr::mutate(value = value* {harvest %>%
                         dplyr::group_by(DISTRICT) %>%
                         dplyr::mutate(prop_harv = prop.table(HARVEST)) %>%
                         tidyr::replace_na(list(prop_harv = 0)) %>%
                         dplyr::filter(DISTRICT == d_idx,
                                       SUBDISTRICT == s_idx,
                                       STAT_WEEK == w_idx) %>%
                         dplyr::pull(prop_harv) %>% max(., 0)})
                 })

      } # s_idx

      for (chain in seq(nchains)) {
        ap_prop2c[[d_idx]][[w_idx]][[chain]] <-
          lapply(ap_prop2b[[d_idx]][[w_idx]],
                 function(ap) {
                   ap[[chain]]
                 }) %>% dplyr::bind_rows()
      } # chain

    } # w_idx

    for (chain in seq(nchains)) {
      ap_prop_all_2[[d_idx]][[chain]] <-
        lapply(ap_prop_all_1[[d_idx]],
               function(apa) {
                 apa[[chain]]
               }) %>% dplyr::bind_rows()
    } # chain

    rm(idx)

  } # d_idx

  popout_a <- list(ap_prop2c = ap_prop2c,
                   ap_prop_all_2 = ap_prop_all_2)

  return(popout_a)

}

#' Format subdistrict step a
#'
#' @param outraw A subset of MAGMA output.
#' @param dat_in A subset of MAGMA data.
#' @param nreps The same as *nreps* in MAGMA model run.
#' @param nburn The same as *nburn* in MAGMA model run.
#' @param thin The same as *thin* in MAGMA model run.
#' @param nchains The same as *nchains* in MAGMA model run.
#' @param keep_burn The same as *keep_burn* in MAGMA model run.
#' @param C Number of age classes.
#' @param D Number of districts.
#' @param S Number of subdistricts.
#' @param W Number of weeks.
#'
#' @importFrom magrittr %>%
#'
#' @noRd
#'
format_subdistrict_pop_a <- function(outraw, dat_in, nreps, nburn, thin, nchains, keep_burn, C, D, S, W) {

  ### data input ### ----
  metadat <- dat_in$metadat # age and iden info
  harvest <- dat_in$harvest %>% # for calculating p0
    dplyr::mutate(n = apply(.[, 2:4], 1, function(aa) {
      dplyr::filter(metadat,
                    district == aa[1],
                    subdist == aa[2],
                    week == aa[3]) %>% nrow()
    }) ) %>%
    dplyr::mutate(HARVEST = HARVEST * (n != 0)) %>%
    dplyr::select(-n)

  groups <- dat_in$groups # vector id for reporting groups (aka groupvec)
  group_names <- dat_in$group_names # reporting groups

  if (isFALSE(keep_burn)|keep_burn == "false") keep_burn = FALSE

  ### prepare output ### ----
  # organization of out_list:
  # [[chain]][[yr]][[dist]][[sub]][[week]][age, pop, itr]

  ap_prop <- ap_combo_subdis <- list()
  ap_prop_all_1 <- ap_prop_all_2 <- list()

  for (d_idx in 1:D) {
    ap_prop[[d_idx]] <- list()
    # ap_combo_subdis[[d_idx]] <- list()
    ap_prop_all_1[[d_idx]] <- ap_prop_all_2[[d_idx]] <- list()

    idx <- 0
    for (s_idx in 1:S[d_idx]) {
      # ap_combo_subdis[[d_idx]][[s_idx]] <- list()
      ap_prop[[d_idx]][[s_idx]] <- list()
      ap_prop_all_1[[d_idx]][[s_idx]] <- ap_prop_all_2[[d_idx]][[s_idx]] <- list()

      for (w_idx in 1:W) {

        # go in each chain and grab age-pop output according to specified sampling period
        # new! ap_prop: [[d]][[s]][[w]][[chain]][age*itr, c(pop, agevec)]
        # pivot_longer %>% [pop*age*itr, c(itr, age, grp, value)]

        # ap_combo_dis or ap_combo_subdis
        # combine age classes and pop groups
        # weighted by harvest
        # this is done for each week in each district
        # new! [[chain]][[d]][[w]][group*itr, c(itr, grpvec, ages)]

        ap_prop[[d_idx]][[s_idx]][[w_idx]] <-
          apply(outraw, MARGIN = 6,
                function(ol) ol[ , , w_idx, s_idx, d_idx] %>%
                  as.data.frame() %>%
                  tidyr::pivot_longer(cols = 1:(ncol(.)-2)) %>%
                  dplyr::mutate(
                    grpvec = rep(group_names[groups[, d_idx], d_idx],
                                 C* (nreps- nburn* isFALSE(keep_burn))/ thin)
                  ) %>%
                  dplyr::select(-name) %>%
                  dplyr::rename(itr = 1, agevec = 2))

        ap_prop_all_1[[d_idx]][[s_idx]][[w_idx]] <-
          lapply(ap_prop[[d_idx]][[s_idx]][[w_idx]],
                 function(apa) {
                   apa %>%
                     dplyr::mutate(value = value* {harvest %>%
                         dplyr::group_by(DISTRICT, SUBDISTRICT) %>%
                         dplyr::mutate(prop_harv = prop.table(HARVEST)) %>%
                         tidyr::replace_na(list(prop_harv = 0)) %>%
                         dplyr::filter(DISTRICT == d_idx,
                                       SUBDISTRICT == s_idx,
                                       STAT_WEEK == w_idx) %>%
                         dplyr::pull(prop_harv) %>% max(., 0)})
                 })

      } # w_idx

      for (chain in seq(nchains)) {
        ap_prop_all_2[[d_idx]][[s_idx]][[chain]] <-
          lapply(ap_prop_all_1[[d_idx]][[s_idx]],
                 function(apa) {
                   apa[[chain]]
                 }) %>% dplyr::bind_rows()
      } # chain

    } # s_idx

    rm(idx)

  } # d_idx

  popout_a <- list(ap_prop = ap_prop,
                   ap_prop_all_2 = ap_prop_all_2)

  return(popout_a)

}

# end pop ----


# age ----

# Summarize age part of MAGMA output (part a)

#' Format district age class part a
#'
#' @param outraw A subset of MAGMA output.
#' @param dat_in A subset of MAGMA data.
#' @param nreps The same as *nreps* in MAGMA model run.
#' @param nburn The same as *nburn* in MAGMA model run.
#' @param thin The same as *thin* in MAGMA model run.
#' @param nchains The same as *nchains* in MAGMA model run.
#' @param keep_burn The same as *keep_burn* in MAGMA model run.
#' @param C Number of age classes.
#' @param D Number of districts.
#' @param S Number of subdistricts.
#' @param W Number of weeks.
#'
#' @importFrom magrittr %>%
#'
#' @noRd
#'
format_district_age_a <- function(outraw, dat_in, nreps, nburn, thin, nchains, keep_burn, C, D, S, W) {

  ### data input ### ----
  metadat <- dat_in$metadat # age and iden info
  harvest <- dat_in$harvest %>% # for calculating p0
    dplyr::mutate(n = apply(.[, 2:4], 1, function(aa) {
      dplyr::filter(metadat,
                    district==aa[1],
                    subdist==aa[2],
                    week==aa[3]) %>% nrow
    }) ) %>%
    dplyr::mutate(HARVEST = HARVEST * (n != 0)) %>%
    dplyr::select(-n)

  groups <- dat_in$groups # vector id for reporting groups (aka groupvec)
  group_names <- dat_in$group_names # reporting groups
  age_class <- dat_in$age_class # vector id for age classes
  age_classes <- dat_in$age_classes # reporting age classes
  age_classes <- age_classes[sort(unique(age_class))]

  wildpops <- dat_in$wildpops
  K <- length(wildpops)

  if (is.null(dat_in$hatcheries)) {
    hatcheries <- NULL
    H <- 0
  } else {
    hatcheries <- dat_in$hatcheries
    H <- length(hatcheries)
  }

  if (isFALSE(keep_burn)|keep_burn == "false") keep_burn = FALSE

  ### prepare output ### ----
  # organization of out_list:
  # [[chain]][[yr]][[dist]][[sub]][[week]][age, pop, itr]

  ap_prop <- ap_combo_dis <- list()

  for (d_idx in 1:D) {
    ap_prop[[d_idx]] <- list()
    ap_combo_dis[[d_idx]] <- list()

    idx <- 0
    for (w_idx in 1:W) {
      ap_prop[[d_idx]][[w_idx]] <- list()

      for (s_idx in 1:S[d_idx]) {

        # go in each chain and grab age-pop output according to specified sampling period
        # new! ap_prop: [[d]][[s]][[w]][[chain]][age*itr, c(pop, agevec)]
        # pivot_longer %>% [pop*age*itr, c(itr, age, grp, value)]

        # ap_combo_dis or ap_combo_subdis
        # combine age classes and pop groups
        # weighted by harvest
        # this is done for each week in each district
        # new! [[chain]][[d]][[w]][group*itr, c(itr, grpvec, ages)]

        idx <- idx + 1
        ap_prop[[d_idx]][[w_idx]][[s_idx]] <-
          apply(outraw, MARGIN = 6,
                function(ol) ol[ , , w_idx, s_idx, d_idx] %>%
                  as.data.frame() %>%
                  tidyr::pivot_longer(cols = 1:(ncol(.)-2)) %>%
                  dplyr::mutate(
                    grpvec = rep(group_names[groups[, d_idx], d_idx],
                                 C* (nreps- nburn* isFALSE(keep_burn))/ thin)
                  ) %>%
                  dplyr::select(-name) %>%
                  dplyr::rename(itr = 1, agevec = 2))

        ap_combo_dis[[d_idx]][[idx]] <-
          lapply(ap_prop[[d_idx]][[w_idx]][[s_idx]], function(apfile) {
            apfile %>%
              dplyr::group_by(itr, agevec, grpvec) %>%
              dplyr::summarise(ppi = sum(value), .groups = "drop") %>%
              dplyr::mutate(ppi = ppi* {harvest %>%
                  dplyr::group_by(DISTRICT) %>%
                  # proportional within a district
                  dplyr::mutate(prop_harv = prop.table(HARVEST)) %>%
                  dplyr::filter(DISTRICT == d_idx,
                                SUBDISTRICT == s_idx,
                                STAT_WEEK == w_idx) %>%
                  dplyr::pull(prop_harv) %>% max(., 0)}) %>%
              tidyr::pivot_wider(names_from = agevec, values_from = ppi)
          })

      } # s_idx
    } # w_idx

    rm(idx)

  } # d_idx

  return(ap_combo_dis)

}

#' Format subdistrict age classes part a
#'
#' @param outraw MAGMA output
#' @param dat_in MAGMA data
#' @param nreps same as MAGMA model run
#' @param nburn same as MAGMA model run
#' @param thin same as MAGMA model run
#' @param nchains same as MAGMA model run
#' @param keep_burn same as MAGMA model run
#' @param C Age classes
#' @param D Districts
#' @param S Subdistricts
#' @param W Weeks
#'
#' @importFrom magrittr %>%
#'
#' @noRd
#'
format_subdistrict_age_a <- function(outraw, dat_in, nreps, nburn, thin, nchains, keep_burn, C, D, S, W) {

  ### data input ### ----
  metadat <- dat_in$metadat # age and iden info
  harvest <- dat_in$harvest %>% # for calculating p0
    dplyr::mutate(n = apply(.[, 2:4], 1, function(aa) {
      dplyr::filter(metadat,
                    district==aa[1],
                    subdist==aa[2],
                    week==aa[3]) %>% nrow
    }) ) %>%
    dplyr::mutate(HARVEST = HARVEST * (n != 0)) %>%
    dplyr::select(-n)

  groups <- dat_in$groups # vector id for reporting groups (aka groupvec)
  group_names <- dat_in$group_names # reporting groups
  age_class <- dat_in$age_class # vector id for age classes
  age_classes <- dat_in$age_classes # reporting age classes
  age_classes <- age_classes[sort(unique(age_class))]

  wildpops <- dat_in$wildpops
  K <- length(wildpops)

  if (is.null(dat_in$hatcheries)) {
    hatcheries <- NULL
    H <- 0
  } else {
    hatcheries <- dat_in$hatcheries
    H <- length(hatcheries)
  }

  if (isFALSE(keep_burn)|keep_burn == "false") keep_burn = FALSE

  ### prepare output ### ----
  # organization of out_list:
  # [[chain]][[yr]][[dist]][[sub]][[week]][age, pop, itr]

  ap_prop <- ap_combo_subdis <- list()

  for (d_idx in 1:D) {
    ap_prop[[d_idx]] <- list()
    ap_combo_subdis[[d_idx]] <- list()

    idx <- 0
    for (s_idx in 1:S[d_idx]) {
      ap_prop[[d_idx]][[s_idx]] <- list()
      ap_combo_subdis[[d_idx]][[s_idx]] <- list()

      for (w_idx in 1:W) {

        # go in each chain and grab age-pop output according to specified sampling period
        # new! ap_prop: [[yr]][[d]][[s]][[w]][[chain]][age*itr, c(pop, agevec)]
        # pivot_longer %>% [pop*age*itr, c(itr, age, grp, value)]

        # ap_combo_dis or ap_combo_subdis
        # combine age classes and pop groups
        # weighted by harvest
        # this is done for each week in each district
        # new! [[chain]][[yr]][[d]][[w]][group*itr, c(itr, grpvec, ages)]

        ap_prop[[d_idx]][[s_idx]][[w_idx]] <-
          apply(outraw, MARGIN = 6,
                function(ol) ol[ , , w_idx, s_idx, d_idx] %>%
                  as.data.frame() %>%
                  tidyr::pivot_longer(cols = 1:(ncol(.)-2)) %>%
                  dplyr::mutate(
                    grpvec = rep(group_names[groups[, d_idx], d_idx],
                                 C* (nreps- nburn* isFALSE(keep_burn))/ thin)
                  ) %>%
                  dplyr::select(-name) %>%
                  dplyr::rename(itr = 1, agevec = 2))

        ap_combo_subdis[[d_idx]][[s_idx]][[w_idx]] <-
          lapply(ap_prop[[d_idx]][[s_idx]][[w_idx]],
                 function(apfile) {

                   apfile %>%
                     dplyr::group_by(itr, agevec, grpvec) %>%
                     dplyr::summarise(ppi = sum(value), .groups = "drop") %>%
                     dplyr::mutate(ppi = ppi* {harvest %>%
                         dplyr::group_by(DISTRICT, SUBDISTRICT) %>%
                         # proportional within a subdistrict
                         dplyr::mutate(prop_harv = prop.table(HARVEST)) %>%
                         dplyr::filter(DISTRICT == d_idx,
                                       SUBDISTRICT == s_idx,
                                       STAT_WEEK == w_idx) %>%
                         dplyr::pull(prop_harv) %>% max(., 0)}) %>%
                     tidyr::pivot_wider(names_from = agevec, values_from = ppi)

                 })

      } # w_idx

    } # s_idx

    rm(idx)

  } # d_idx

  return(ap_combo_subdis)

}


# end age ----


















