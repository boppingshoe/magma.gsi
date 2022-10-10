

## file contains functions for formatting tbr magma output

## summarize pop and age at the same time ----

#' Format district
#'
#' @param outraw
#' @param dat_in
#' @param nreps
#' @param nburn
#' @param thin
#' @param nchains
#' @param keep_burn
#' @param C
#' @param D
#' @param S
#' @param W
#'
#' @export
#' @importFrom magrittr %>%
#'
#' @noRd
#'
format_district <- function(outraw, dat_in, nreps, nburn, thin, nchains, keep_burn, C, D, S, W) {

  ### data input ### ----
  metadat <- dat_in$metadat # age and iden info
  harvest <- dat_in$harvest %>% # for calculating p0
    dplyr::mutate(n = apply(.[, 2:4], 1, function(aa) {
      dplyr::filter(metadat,
                    # year == aa[1],
                    district == aa[1],
                    subdist == aa[2],
                    week == aa[3]) %>% nrow()
    })) %>%
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

  p_zero <- array(
    apply(
      table(metadat[, c("district", "subdist", "week", "iden")], exclude = NULL),
      seq(3),
      function(alf) any(alf!= 0) # prop = 0 if no sample
    ), c(D, max(S), W, K+ H))

  if (isFALSE(keep_burn)|keep_burn == "false") keep_burn = FALSE

  ### prepare output ### ----
  # organization of out_list:
  # [[chain]][[yr]][[dist]][[sub]][[week]][age, pop, itr]

  ap_prop <- ap_prop2 <- ap_prop2b <- ap_prop2c <- ap_combo_dis <- list()
  ap_prop_all_1 <- ap_prop_all_2 <- list()

  for (d_idx in 1:D) {
    ap_prop[[d_idx]] <- list()
    ap_prop2b[[d_idx]] <- list()
    ap_prop2c[[d_idx]] <- list()
    ap_combo_dis[[d_idx]] <- list()
    ap_prop_all_1[[d_idx]] <- ap_prop_all_2[[d_idx]] <- list()

    idx <- 0
    for (w_idx in 1:W) {
      ap_prop[[d_idx]][[w_idx]] <- list()
      ap_prop2b[[d_idx]][[w_idx]] <- list()
      ap_prop2c[[d_idx]][[w_idx]] <- list()

      for (s_idx in 1:S[d_idx]) {

        # go in each chain and grab age-pop output according to specified sampling period
        # pivot_longer %>% [pop*age*itr, c(itr, age, grp, value)]

        # ap_combo_dis or ap_combo_subdis
        # combine age classes and pop groups
        # weighted by harvest
        # this is done for each week in each district

        idx <- idx + 1
        ap_prop[[d_idx]][[w_idx]][[s_idx]] <-
          apply(outraw, MARGIN = 6,
                function(ol) ol[ , , w_idx, s_idx, d_idx] %>%
                  tibble::as_tibble() %>%
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

  # exclude burn-in while calculate r hat
  keep_list <- ((nburn*keep_burn + 1):(nreps - nburn*isFALSE(keep_burn)))[!((nburn*keep_burn + 1):(nreps - nburn*isFALSE(keep_burn))) %% thin] / thin

  # place holders/empty objects for age/pop summaries
  p_combo <- mc_pop <- summ_pop <- list()
  p_combo_all <- mc_pop_all <- summ_pop_all <- list()
  p_idx <- 0
  ap_prop_grp <- mc_age <- summ_age <- list()
  a_idx <- 0

  for (d_idx in 1:D) {
    ## pops all dist ##
    # sum over ages and combine pop groups

    p_combo_all[[d_idx]] <-
      lapply(ap_prop_all_2[[d_idx]],
             function(apchain) {
               apchain %>%
                 dplyr::select(-agevec) %>%
                 dplyr::group_by(itr, grpvec) %>%
                 dplyr::summarise(p = sum(value), .groups = "drop") %>%
                 tidyr::pivot_wider(names_from = grpvec, values_from = p) %>%
                 dplyr::select(-itr)
             })

    mc_pop_all[[d_idx]] <- coda::as.mcmc.list(
      lapply(p_combo_all[[d_idx]],
             function(rlist) coda::mcmc(rlist[keep_list,])) )

    summ_pop_all[[d_idx]] <-
      lapply(p_combo_all[[d_idx]], function(rlist) rlist[keep_list,]) %>%
      dplyr::bind_rows %>%
      tidyr::pivot_longer(cols = 1:ncol(.)) %>%
      dplyr::group_by(name) %>%
      dplyr::summarise(
        mean = mean(value),
        median = median(value),
        sd = sd(value),
        ci.05 = quantile(value, 0.05),
        ci.95 = quantile(value, 0.95),
        p0 = mean(value < (0.5/ max(1, ( any(p_zero[1, d_idx, , , ])* harvest %>% dplyr::group_by(DISTRICT) %>% dplyr::summarise(sum_harv = sum(HARVEST), .groups = "drop") %>% dplyr::filter(DISTRICT == d_idx) %>% dplyr::pull(sum_harv) )))),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        GR = {if (nchains > 1) {
          coda::gelman.diag(mc_pop_all[[d_idx]],
                            transform = TRUE,
                            autoburnin = FALSE,
                            multivariate = FALSE)$psrf[,"Point est."] %>%
            .[order(group_names[seq(ncol(p_combo_all[[d_idx]][[1]])), d_idx])]
        } else {NA}},
        n_eff = coda::effectiveSize(mc_pop_all[[d_idx]]) %>%
          .[order(group_names[seq(ncol(p_combo_all[[d_idx]][[1]])), d_idx])]
      ) %>%
      dplyr::rename(group = name)

    for (w_idx in 1:W) {
      ## pop ##
      # sum over ages and combine pop groups
      p_idx <- p_idx + 1 # index individual sampling period

      p_combo[[p_idx]] <-
        lapply(ap_prop2c[[d_idx]][[w_idx]],
               function(apfile) {
                 apfile %>%
                   dplyr::select(-agevec) %>%
                   dplyr::group_by(itr, grpvec) %>%
                   dplyr::summarise(p = sum(value), .groups = "drop") %>%
                   tidyr::pivot_wider(names_from = grpvec, values_from = p) %>%
                   dplyr::select(-itr)
               })

      mc_pop[[p_idx]] <- coda::as.mcmc.list(
        lapply(p_combo[[p_idx]],
               function(rlist) coda::mcmc(rlist[keep_list,])) )

      summ_pop[[p_idx]] <-
        lapply(p_combo[[p_idx]], function(rlist) rlist[keep_list,]) %>%
        dplyr::bind_rows %>%
        tidyr::pivot_longer(cols = 1:ncol(.)) %>%
        dplyr::group_by(name) %>%
        dplyr::summarise(
          mean = mean(value),
          median = median(value),
          sd = sd(value),
          ci.05 = quantile(value, 0.05),
          ci.95 = quantile(value, 0.95),
          p0 = mean(value < (0.5/ max(1, ( any(p_zero[1, d_idx, , w_idx, ])* harvest %>% dplyr::group_by(DISTRICT, STAT_WEEK) %>% dplyr::summarise(sum_harv = sum(HARVEST), .groups = "drop") %>% dplyr::filter(DISTRICT == d_idx, STAT_WEEK== w_idx) %>% dplyr::pull(sum_harv) )))),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          GR = {if (nchains > 1) {
            coda::gelman.diag(mc_pop[[p_idx]],
                              transform = TRUE,
                              autoburnin = FALSE,
                              multivariate = FALSE)$psrf[,"Point est."] %>%
              .[order(group_names[seq(ncol(p_combo[[p_idx]][[1]])), d_idx])]
          } else {NA}},
          n_eff = coda::effectiveSize(mc_pop[[p_idx]]) %>%
            .[order(group_names[seq(ncol(p_combo[[p_idx]][[1]])), d_idx])]
        ) %>%
        dplyr::rename(group = name)
    } # W

    ## age (loop through each district) ##
    ap_combo_all = list()
    for(chain in seq(nchains)) {
      ap_combo_all[[chain]] <-
        lapply(ap_combo_dis[[d_idx]],
               function(apc) {
                 apc[[chain]]
               }) %>%
        dplyr::bind_rows() %>%
        dplyr::group_by(itr, grpvec) %>%
        dplyr::summarise(across(.fns = sum), .groups = "drop") %>%
        tidyr::pivot_longer(-c(itr, grpvec)) %>%
        dplyr::group_by(itr, grpvec) %>%
        dplyr::mutate(value = prop.table(value)) %>%
        dplyr::ungroup() %>%
        tidyr::pivot_wider() %>%
        dplyr::select(-itr) %>%
        stats::setNames(., c("grpvec", age_classes))
    } # combine all individual districts

    for (grp in group_names[, d_idx] %>% stats::na.omit) {
      a_idx <- a_idx + 1

      ap_prop_grp[[a_idx]] <-
        lapply(ap_combo_all, function(aoc) {
          ar_temp <- aoc %>%
            dplyr::filter(grpvec == grp) %>%
            stats::setNames(., c("grpvec", age_classes))
          return(ar_temp)# %>% select(-grpvec))
        }) # separate rep groups into its own output

      mc_age[[a_idx]] <-
        coda::as.mcmc.list(
          lapply(ap_prop_grp[[a_idx]], function(rlist) {
            coda::mcmc(rlist[keep_list,] %>% dplyr::select(-grpvec))
          }) )

      harv_dis <- harvest %>%
        dplyr::filter(DISTRICT == d_idx) %>%
        dplyr::pull(HARVEST) %>%
        mean()

      summ_age[[a_idx]] <-
        lapply(ap_prop_grp[[a_idx]], function(rlist) rlist[keep_list,]) %>%
        dplyr::bind_rows %>%
        dplyr::mutate(stock_prop = rowSums(.[,-1])) %>%
        tidyr::pivot_longer(-c(stock_prop, grpvec)) %>%
        dplyr::group_by(name, grpvec) %>%
        dplyr::summarise(
          mean = mean(value),
          median = median(value),
          sd = sd(value),
          ci.05 = quantile(value, 0.05),
          ci.95 = quantile(value, 0.95),
          p0 = mean( value < (0.5/ max(1, harv_dis* stock_prop)) ),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          GR = {if (nchains > 1) {
            coda::gelman.diag(mc_age[[a_idx]],
                              transform = TRUE,
                              autoburnin = FALSE,
                              multivariate = FALSE)$psrf[,"Point est."] %>%
              .[order(age_classes)]
          } else {NA}},
          n_eff = coda::effectiveSize(mc_age[[a_idx]]) %>%
            .[order(age_classes)]
        ) %>%
        dplyr::rename(age = name,
                      group = grpvec) %>%
        dplyr::relocate(group, .before = age)
    } # grp

  } # D

  out <- list()

  out$age_prop <-
    lapply(ap_prop_grp, function(rot) {
      dplyr::bind_rows(rot) %>%
        dplyr::mutate(
          chain = rep(1:nchains,
                      each = (nreps - nburn*isFALSE(keep_burn)) / thin),
          itr = rep(1:((nreps - nburn*isFALSE(keep_burn)) / thin),
                    times = nchains)
        )
    }) # add id for chain and iteration

  out$age_summ <- summ_age

  out$pop_prop <-
    lapply(p_combo, function(rot) {
      dplyr::bind_rows(rot) %>%
        dplyr::mutate(
          chain = rep(1:nchains,
                      each = (nreps - nburn*isFALSE(keep_burn)) / thin),
          itr = rep(1:((nreps - nburn*isFALSE(keep_burn)) / thin),
                    times = nchains)
        )
    }) # add id for chain and iteration

  out$pop_summ <- summ_pop

  out$pop_prop_all <-
    lapply(p_combo_all, function(rot) {
      dplyr::bind_rows(rot) %>%
        dplyr::mutate(
          chain = rep(1:nchains,
                      each = (nreps - nburn*isFALSE(keep_burn)) / thin),
          itr = rep(1:((nreps - nburn*isFALSE(keep_burn)) / thin),
                    times = nchains)
        )
    }) # add id for chain and iteration

  out$pop_summ_all <- summ_pop_all

  return(out)

}


#' Format subdistrict
#'
#' @param outraw
#' @param dat_in
#' @param nreps
#' @param nburn
#' @param thin
#' @param nchains
#' @param keep_burn
#' @param C
#' @param D
#' @param S
#' @param W
#'
#' @export
#' @importFrom magrittr %>%
#'
#' @noRd
#'
format_subdistrict <- function(outraw, dat_in, nreps, nburn, thin, nchains, keep_burn, C, D, S, W) {

  ### data input ### ----
  metadat <- dat_in$metadat # age and iden info
  harvest <- dat_in$harvest %>% # for calculating p0
    dplyr::mutate(n = apply(.[, 2:4], 1, function(aa) {
      dplyr::filter(metadat,
                    # year==aa[1],
                    district==aa[1],
                    subdist==aa[2],
                    week==aa[3]) %>% nrow
    }) ) %>%
    mutate(HARVEST = HARVEST * (n != 0)) %>%
    select(-n)

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

  p_zero <- array(
    apply(
      table(metadat[, c("district", "subdist", "week", "iden")], exclude = NULL),
      seq(3),
      function(alf) any(alf!= 0) # prop = 0 if no sample
    ), c(D, max(S), W, K+ H))

  if (isFALSE(keep_burn)|keep_burn == "false") keep_burn = FALSE

  ### prepare output ### ----
  # organization of out_list:
  # [[chain]][[yr]][[dist]][[sub]][[week]][age, pop, itr]

  ap_prop <- ap_prop2 <- ap_prop2b <- ap_prop2c <- ap_combo_subdis <- list()
  ap_prop_all_1 <- ap_prop_all_2 <- list()

  for (d_idx in 1:D) {
    ap_prop[[d_idx]] <- list()
    ap_combo_subdis[[d_idx]] <- list()
    ap_prop_all_1[[d_idx]] <- ap_prop_all_2[[d_idx]] <- list()

    idx <- 0
    for (s_idx in 1:S[d_idx]) {
      ap_combo_subdis[[d_idx]][[s_idx]] <- list()
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
                  tibble::as_tibble() %>%
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

  # exclude burn-in while calculate r hat
  keep_list <- ((nburn*keep_burn + 1):(nreps - nburn*isFALSE(keep_burn)))[!((nburn*keep_burn + 1):(nreps - nburn*isFALSE(keep_burn))) %% thin] / thin

  # place holders/empty objects for age/pop summaries
  p_combo <- mc_pop <- summ_pop <- list()
  p_combo_all <- mc_pop_all <- summ_pop_all <- list()
  p_idx <- 0
  ap_prop_grp <- mc_age <- summ_age <- list()
  a_idx <- 0

  for (d_idx in 1:D) {
    for (s_idx in 1:S[d_idx]) {

      ## pops all dist ##
      # sum over ages and combine pop groups

      p_combo_all[[max(S)*(d_idx-1)+s_idx]] <-
        lapply(ap_prop_all_2[[d_idx]][[s_idx]],
               function(apchain) {
                 apchain %>%
                   dplyr::select(-agevec) %>%
                   dplyr::group_by(itr, grpvec) %>%
                   dplyr::summarise(p = sum(value), .groups = "drop") %>%
                   tidyr::pivot_wider(names_from = grpvec, values_from = p) %>%
                   dplyr::select(-itr)
               })

      mc_pop_all[[max(S)*(d_idx-1)+s_idx]] <- coda::as.mcmc.list(
        lapply(p_combo_all[[max(S)*(d_idx-1)+s_idx]],
               function(rlist) coda::mcmc(rlist[keep_list,])) )

      summ_pop_all[[max(S)*(d_idx-1)+s_idx]] <-
        lapply(p_combo_all[[max(S)*(d_idx-1)+s_idx]], function(rlist) rlist[keep_list,]) %>%
        dplyr::bind_rows %>%
        tidyr::pivot_longer(cols = 1:ncol(.)) %>%
        dplyr::group_by(name) %>%
        dplyr::summarise(
          mean = mean(value),
          median = median(value),
          sd = sd(value),
          ci.05 = quantile(value, 0.05),
          ci.95 = quantile(value, 0.95),
          p0 = mean(value < (0.5/ max(1, ( any(p_zero[1, d_idx, , , ])* harvest %>% dplyr::group_by(DISTRICT, SUBDISTRICT) %>% dplyr::summarise(sum_harv = sum(HARVEST), .groups = "drop") %>% dplyr::filter(DISTRICT == d_idx, SUBDISTRICT == s_idx) %>% dplyr::pull(sum_harv) )))),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          GR = {if (nchains > 1) {
            coda::gelman.diag(mc_pop_all[[max(S)*(d_idx-1)+s_idx]],
                              transform = TRUE,
                              autoburnin = FALSE,
                              multivariate = FALSE)$psrf[,"Point est."] %>%
              .[order(group_names[seq(ncol(p_combo_all[[max(S)*(d_idx-1)+s_idx]][[1]])), d_idx])]
          } else {NA}},
          n_eff = coda::effectiveSize(mc_pop_all[[max(S)*(d_idx-1)+s_idx]]) %>%
            .[order(group_names[seq(ncol(p_combo_all[[max(S)*(d_idx-1)+s_idx]][[1]])), d_idx])]
        ) %>%
        dplyr::rename(group = name)

      for (w_idx in 1:W) {

        ## pop ##
        # sum over ages and combine pop groups
        p_idx <- p_idx + 1 # index individual sampling period

        p_combo[[p_idx]] <-
          lapply(ap_prop[[d_idx]][[s_idx]][[w_idx]],
                 function(apfile) {
                   apfile %>%
                     dplyr::select(-agevec) %>%
                     dplyr::group_by(itr, grpvec) %>%
                     dplyr::summarise(p = sum(value), .groups = "drop") %>%
                     tidyr::pivot_wider(names_from = grpvec, values_from = p) %>%
                     dplyr::select(-itr) #%>%
                   # stats::setNames(group_names[seq(ncol(.)), d_idx])
                 })

        mc_pop[[p_idx]] <- coda::as.mcmc.list(
          lapply(p_combo[[p_idx]],
                 function(rlist) coda::mcmc(rlist[keep_list,])) )

        summ_pop[[p_idx]] <-
          lapply(p_combo[[p_idx]], function(rlist) rlist[keep_list,]) %>%
          dplyr::bind_rows %>%
          tidyr::pivot_longer(cols = 1:ncol(.)) %>%
          dplyr::group_by(name) %>%
          dplyr::summarise(
            mean = mean(value),
            median = median(value),
            sd = sd(value),
            ci.05 = quantile(value, 0.05),
            ci.95 = quantile(value, 0.95),
            p0 = mean(value < (0.5/ max(1, (any(p_zero[1, d_idx, s_idx, w_idx, ])* dplyr::filter(harvest, DISTRICT== d_idx, SUBDISTRICT== s_idx, STAT_WEEK== w_idx)$HARVEST)))),
            .groups = "drop"
          ) %>%
          mutate(
            GR = {if (nchains > 1) {
              coda::gelman.diag(mc_pop[[p_idx]],
                                transform = TRUE,
                                autoburnin = FALSE,
                                multivariate = FALSE)$psrf[,"Point est."] %>%
                .[order(group_names[seq(ncol(p_combo[[p_idx]][[1]])), d_idx])]
            } else {NA}},
            n_eff = coda::effectiveSize(mc_pop[[p_idx]]) %>%
              .[order(group_names[seq(ncol(p_combo[[p_idx]][[1]])), d_idx])]
          ) %>%
          dplyr::rename(group = name)

      } # W

      ## age (loop through each subdistrict) ##
      ap_combo_all = list()
      for(chain in seq(nchains)) {
        ap_combo_all[[chain]] <-
          lapply(ap_combo_subdis[[d_idx]][[s_idx]],
                 function(apc) {
                   apc[[chain]]
                 }) %>%
          dplyr::bind_rows() %>%
          dplyr::group_by(itr, grpvec) %>%
          dplyr::summarise(across(.fns = sum), .groups = "drop") %>%
          tidyr::pivot_longer(-c(itr, grpvec)) %>%
          dplyr::group_by(itr, grpvec) %>%
          dplyr::mutate(value = prop.table(value)) %>%
          dplyr::ungroup() %>%
          tidyr::pivot_wider() %>%
          dplyr::select(-itr) %>%
          stats::setNames(., c("grpvec", age_classes))
      } # combine all individual districts

      for (grp in group_names[, d_idx] %>% na.omit) {
        a_idx <- a_idx + 1

        ap_prop_grp[[a_idx]] <-
          lapply(ap_combo_all, function(aoc) {
            ar_temp <- aoc %>%
              dplyr::filter(grpvec == grp) %>%
              stats::setNames(., c("grpvec", age_classes))
            return(ar_temp)# %>% select(-grpvec))
          }) # separate rep groups into its own output

        mc_age[[a_idx]] <-
          coda::as.mcmc.list(
            lapply(ap_prop_grp[[a_idx]], function(rlist) {
              coda::mcmc(rlist[keep_list,] %>% dplyr::select(-grpvec))
            }) )

        harv_subdis <- harvest %>%
          dplyr::filter(DISTRICT == d_idx, SUBDISTRICT == s_idx) %>%
          dplyr::pull(HARVEST) %>%
          mean()

        summ_age[[a_idx]] <-
          lapply(ap_prop_grp[[a_idx]], function(rlist) rlist[keep_list,]) %>%
          dplyr::bind_rows %>%
          dplyr::mutate(stock_prop = rowSums(.[,-1])) %>%
          tidyr::pivot_longer(-c(stock_prop, grpvec)) %>%
          dplyr::group_by(name, grpvec) %>%
          dplyr::summarise(
            mean = mean(value),
            median = median(value),
            sd = sd(value),
            ci.05 = quantile(value, 0.05),
            ci.95 = quantile(value, 0.95),
            p0 = mean( value < (0.5/ max(1, harv_subdis* stock_prop)) ),
            .groups = "drop"
          ) %>%
          dplyr::mutate(
            GR = {if (nchains > 1) {
              coda::gelman.diag(mc_age[[a_idx]],
                                transform = TRUE,
                                autoburnin = FALSE,
                                multivariate = FALSE)$psrf[,"Point est."] %>%
                .[order(age_classes)]
            } else {NA}},
            n_eff = coda::effectiveSize(mc_age[[a_idx]]) %>%
              .[order(age_classes)]
          ) %>%
          dplyr::rename(age = name,
                 group = grpvec) %>%
          dplyr::relocate(group, .before = age)
      } # grp

    } # S
  } # D

  out <- list()

  out$age_prop <-
    lapply(ap_prop_grp, function(rot) {
      dplyr::bind_rows(rot) %>%
        dplyr::mutate(
          chain = rep(1:nchains,
                      each = (nreps - nburn*isFALSE(keep_burn)) / thin),
          itr = rep(1:((nreps - nburn*isFALSE(keep_burn)) / thin),
                    times = nchains)
        )
    }) # add id for chain and iteration

  out$age_summ <- summ_age

  out$pop_prop <-
    lapply(p_combo, function(rot) {
      dplyr::bind_rows(rot) %>%
        dplyr::mutate(
          chain = rep(1:nchains,
                      each = (nreps - nburn*isFALSE(keep_burn)) / thin),
          itr = rep(1:((nreps - nburn*isFALSE(keep_burn)) / thin),
                    times = nchains)
        )
    }) # add id for chain and iteration

  out$pop_summ <- summ_pop

  out$pop_prop_all <-
    lapply(p_combo_all, function(rot) {
      dplyr::bind_rows(rot) %>%
        dplyr::mutate(
          chain = rep(1:nchains,
                      each = (nreps - nburn*isFALSE(keep_burn)) / thin),
          itr = rep(1:((nreps - nburn*isFALSE(keep_burn)) / thin),
                    times = nchains)
        )
    }) # add id for chain and iteration

  out$pop_summ_all <- summ_pop_all

  return(out)

}


#' Summarize MAGMA output pop and age at the same time
#'
#' @param outraw MAGMA output
#' @param dat_in MAGMA input data
#' @param nreps The same *nreps* as in MAGMA model run
#' @param nburn The same *nburn* as in MAGMA model run
#' @param thin The same *thin* as in MAGMA model run
#' @param nchains The same *nchains* as in MAGMA model run
#' @param keep_burn The same *keep_burn* as in MAGMA model run
#' @param summ_level Summarize at district or subdistrict level
#'
#' @return Summary tables for reporting groups and age classes
#'
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' # format data
#' wd <- "D:/bobby_adfg/backup_013122/projects/magma/test_TBR" # path to data folder
#' magma_data <- magmatize_data(wd = wd, save_data = FALSE)
#'
#' # model run
#' magma_out <- magmatize_mdl(magma_data, nreps = 50, nburn = 25, thin = 1, nchains = 3)
#'
#' # summary
#' magma_summ <- magmatize_summ(magma_out, nreps = 50, nburn = 25, thin = 1, nchains = 3, summ_level = "district")
#'
magmatize_summ <- function(outraw, dat_in, nreps, nburn, thin, nchains, keep_burn = FALSE, summ_level) {

  ### info needed ### ----
  C <- dat_in$C # number of age classes
  D <- dplyr::n_distinct(dat_in$metadat$district) # number of districts
  S <- dat_in$metadat %>%
    dplyr::group_by(district) %>%
    dplyr::summarise(S = n_distinct(subdist), .groups = "drop") %>%
    dplyr::pull(S) # number of subdistricts
  W <- max(dplyr::n_distinct(dat_in$metadat$week), length(dat_in$stat_weeks)) # number of stats weeks

  ng <- apply(dat_in$groups, 2, max) # number of groups

  if (!is.null(dat_in$districts)) dist_names <- dat_in$districts
  if (!is.null(dat_in$subdistricts)) subdist_names <- dat_in$subdistricts
  if (!is.null(dat_in$stat_weeks)) week_names <- dat_in$stat_weeks

  ### prepare output ### ----

  message("Preparing output (patience grasshopper...)")
  prep_time <- Sys.time()

  if (summ_level == "district") {
    out <- format_district(outraw, dat_in, nreps, nburn, thin, nchains, keep_burn, C, D, S, W)
  } else if (summ_level == "subdistrict") {
    out <- format_subdistrict(outraw, dat_in, nreps, nburn, thin, nchains, keep_burn, C, D, S, W)
  }

  # id for age outputs
  if (summ_level == "district") {
    aout_names <- rep(NA, sum(ng))
    i <- 0
    for (d_i in 1:D) {
      for (g_i in 1:ng[d_i]) {
        i <- i + 1
        aout_names[i] <-
          if (exists("dist_names")) {
            paste(paste0("(", i, ")"),
                  paste0("D", dist_names[d_i]),
                  paste0(dat_in$group_names[g_i, d_i]),
                  sep = "_")
          } else {
            paste(paste0("(", i, ")"),
                  paste0("D", d_i),
                  paste0(dat_in$group_names[g_i, d_i]),
                  sep = "_")
          }
      } # g_i
    } # d_i
  } else if (summ_level == "subdistrict") {
    aout_names <- rep(NA, sum(ng* S))
    i <- 0
    for (d_i in 1:D) {
      for(s_i in 1:S[d_i]) {
        for (g_i in 1:ng[d_i]) {
          i <- i + 1
          aout_names[i] <-
            if (all(exists("dist_names"), exists("subdist_names"))) {
              paste(paste0("(", i, ")"),
                    paste0("D", dist_names[d_i]),
                    paste0("SubD", subdist_names[[d_i]][s_i]),
                    paste0(dat_in$group_names[g_i, d_i]),
                    sep = "_")
            } else {
              paste(paste0("(", i, ")"),
                    paste0("D", d_i),
                    paste0("SubD", s_i),
                    paste0(dat_in$group_names[g_i, d_i]),
                    sep = "_")
            }
        } # g_i
      } # s_i
    } # d_i
  } # end if

  for(li in 1:2) {
    names(out[[li]]) <- aout_names
  } # id for pop output

  if (summ_level == "district") {

    pout_names <- rep(NA, D*W)
    j <- 1
    for (d_i in 1:D) {
      for (w_i in 1:W) {
        pout_names[j] <-
          if (all(exists("dist_names"), exists("week_names"))) {
            paste(paste0("(", j, ")"),
                  paste0("D", dist_names[d_i]),
                  paste0("Wk", week_names[w_i]),
                  sep = "_")
          } else {
            paste(paste0("(", j, ")"),
                  paste0("D", d_i),
                  paste0("Wk", w_i),
                  sep = "_")
          }
        j <- j + 1
      } # w_i
    } # d_i

  } else if (summ_level == "subdistrict") {

    pout_names <- rep(NA, sum(S*W))
    j <- 1
    for (d_i in 1:D) {
      for (s_i in 1:S[d_i]) {
        for (w_i in 1:W) {
          pout_names[j] <-
            if (all(exists("dist_names"), exists("subdist_names"), exists("week_names"))) {
              paste(paste0("(", j, ")"),
                    paste0("D", dist_names[d_i]),
                    paste0("Sub", subdist_names[[d_i]][s_i]),
                    paste0("Wk", week_names[w_i]),
                    sep = "_")
            } else {
              paste(paste0("(", j, ")"),
                    paste0("D", d_i),
                    paste0("Sub", s_i),
                    paste0("Wk", w_i),
                    sep = "_")
            }
          j <- j + 1
        } # w_i
      } # s_i
    } # d_i

  } # end if dist/subdist

  for(li in 3:4) {
    names(out[[li]]) <- pout_names
  } # id for pop output

  if (summ_level == "district") {
    for(li in 5:6) {
      names(out[[li]]) <- dist_names
    }
  } else if (summ_level == "subdistrict") {
    for(li in 5:6) {
      names(out[[li]]) <- subdist_names %>% unlist() %>% names()
    }
  }

  print(Sys.time() - prep_time)
  message(Sys.time())

  return(out)

}


## pop and age separately ----

# district

format_district_pop <- function(outraw, dat_in, nreps, nburn, thin, nchains, keep_burn, C, D, S, W) {

  ### data input ### ----
  metadat <- dat_in$metadat # age and iden info
  harvest <- dat_in$harvest %>% # for calculating p0
    mutate(n = apply(.[,1:4], 1, function(aa) {
      filter(metadat,
             year==aa[1],
             district==aa[2],
             subdist==aa[3],
             week==aa[4]) %>% nrow
    }) ) %>%
    mutate(HARVEST = HARVEST * (n != 0)) %>%
    select(-n)

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

  p_zero <- array(
    apply(
      table(metadat[, c("year","district","subdist","week","iden")], exclude = NULL),
      seq(4),
      function(alf) any(alf!= 0) # prop = 0 if no sample
    ), c(1, D, max(S), W, K+ H))

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
        # new! ap_prop: [[yr]][[d]][[s]][[w]][[chain]][age*itr, c(pop, agevec)]
        # pivot_longer %>% [pop*age*itr, c(itr, age, grp, value)]

        # ap_combo_dis or ap_combo_subdis
        # combine age classes and pop groups
        # weighted by harvest
        # this is done for each week in each district
        # new! [[chain]][[yr]][[d]][[w]][group*itr, c(itr, grpvec, ages)]

        idx <- idx + 1
        ap_prop[[d_idx]][[w_idx]][[s_idx]] <-
          apply(outraw, MARGIN = 7,
                function(ol) ol[ , , w_idx, s_idx, d_idx, 1] %>%
                  as_tibble() %>%
                  pivot_longer(cols = 1:(ncol(.)-2)) %>%
                  mutate(
                    grpvec = rep(group_names[groups[, d_idx], d_idx],
                                 C* (nreps- nburn* isFALSE(keep_burn))/ thin)
                  ) %>%
                  select(-name) %>%
                  rename(itr = 1, agevec = 2))

        ap_prop2b[[d_idx]][[w_idx]][[s_idx]] <-
          lapply(ap_prop[[d_idx]][[w_idx]][[s_idx]],
                 function(ap) {
                   ap %>%
                     mutate(value = value* {harvest %>%
                         group_by(DISTRICT, STAT_WEEK) %>%
                         mutate(prop_harv = prop.table(HARVEST)) %>%
                         replace_na(list(prop_harv = 0)) %>%
                         filter(DISTRICT == d_idx,
                                SUBDISTRICT == s_idx,
                                STAT_WEEK == w_idx) %>%
                         pull(prop_harv) %>% max(., 0)})
                 })

        ap_prop_all_1[[d_idx]][[W*(s_idx-1)+w_idx]] <-
          lapply(ap_prop[[d_idx]][[w_idx]][[s_idx]],
                 function(apa) {
                   apa %>%
                     mutate(value = value* {harvest %>%
                         group_by(DISTRICT) %>%
                         mutate(prop_harv = prop.table(HARVEST)) %>%
                         replace_na(list(prop_harv = 0)) %>%
                         filter(DISTRICT == d_idx,
                                SUBDISTRICT == s_idx,
                                STAT_WEEK == w_idx) %>%
                         pull(prop_harv) %>% max(., 0)})
                 })

      } # s_idx

      for (chain in seq(nchains)) {
        ap_prop2c[[d_idx]][[w_idx]][[chain]] <-
          lapply(ap_prop2b[[d_idx]][[w_idx]],
                 function(ap) {
                   ap[[chain]]
                 }) %>% bind_rows()
      } # chain

    } # w_idx

    for (chain in seq(nchains)) {
      ap_prop_all_2[[d_idx]][[chain]] <-
        lapply(ap_prop_all_1[[d_idx]],
               function(apa) {
                 apa[[chain]]
               }) %>% bind_rows()
    } # chain

    rm(idx)

  } # d_idx

  # exclude burn-in while calculate r hat
  keep_list <- ((nburn*keep_burn + 1):(nreps - nburn*isFALSE(keep_burn)))[!((nburn*keep_burn + 1):(nreps - nburn*isFALSE(keep_burn))) %% thin] / thin

  # place holders/empty objects for age/pop summaries
  p_combo <- mc_pop <- summ_pop <- list()
  p_combo_all <- mc_pop_all <- summ_pop_all <- list()
  p_idx <- 0

  for (d_idx in 1:D) {
    ## pops all dist ##
    # sum over ages and combine pop groups

    p_combo_all[[d_idx]] <-
      lapply(ap_prop_all_2[[d_idx]],
             function(apchain) {
               apchain %>%
                 select(-agevec) %>%
                 group_by(itr, grpvec) %>%
                 summarise(p = sum(value), .groups = "drop") %>%
                 pivot_wider(names_from = grpvec, values_from = p) %>%
                 select(-itr)
             })

    mc_pop_all[[d_idx]] <- as.mcmc.list(
      lapply(p_combo_all[[d_idx]],
             function(rlist) mcmc(rlist[keep_list,])) )

    summ_pop_all[[d_idx]] <-
      lapply(p_combo_all[[d_idx]], function(rlist) rlist[keep_list,]) %>%
      bind_rows %>%
      pivot_longer(cols = 1:ncol(.)) %>%
      group_by(name) %>%
      summarise(
        mean = mean(value),
        median = median(value),
        sd = sd(value),
        ci.05 = quantile(value, 0.05),
        ci.95 = quantile(value, 0.95),
        p0 = mean(value < (0.5/ max(1, ( any(p_zero[1, d_idx, , , ])* harvest %>% group_by(DISTRICT) %>% summarise(sum_harv = sum(HARVEST), .groups = "drop") %>% filter(DISTRICT == d_idx) %>% pull(sum_harv) )))),
        .groups = "drop"
      ) %>%
      mutate(
        GR = {if (nchains > 1) {
          coda::gelman.diag(mc_pop_all[[d_idx]],
                            transform = TRUE,
                            autoburnin = FALSE,
                            multivariate = FALSE)$psrf[,"Point est."] %>%
            .[order(group_names[seq(ncol(p_combo_all[[d_idx]][[1]])), d_idx])]
        } else {NA}},
        n_eff = coda::effectiveSize(mc_pop_all[[d_idx]]) %>%
          .[order(group_names[seq(ncol(p_combo_all[[d_idx]][[1]])), d_idx])]
      ) %>%
      rename(group = name)

    for (w_idx in 1:W) {
      ## pop individual weeks ##
      # sum over ages and combine pop groups
      p_idx <- p_idx + 1 # index individual sampling period

      p_combo[[p_idx]] <-
        lapply(ap_prop2c[[d_idx]][[w_idx]],
               function(apfile) {
                 apfile %>%
                   select(-agevec) %>%
                   group_by(itr, grpvec) %>%
                   summarise(p = sum(value), .groups = "drop") %>%
                   pivot_wider(names_from = grpvec, values_from = p) %>%
                   select(-itr)
               })

      mc_pop[[p_idx]] <- as.mcmc.list(
        lapply(p_combo[[p_idx]],
               function(rlist) mcmc(rlist[keep_list,])) )

      summ_pop[[p_idx]] <-
        lapply(p_combo[[p_idx]], function(rlist) rlist[keep_list,]) %>%
        bind_rows %>%
        pivot_longer(cols = 1:ncol(.)) %>%
        group_by(name) %>%
        summarise(
          mean = mean(value),
          median = median(value),
          sd = sd(value),
          ci.05 = quantile(value, 0.05),
          ci.95 = quantile(value, 0.95),
          p0 = mean(value < (0.5/ max(1, ( any(p_zero[1, d_idx, , w_idx, ])* harvest %>% group_by(DISTRICT, STAT_WEEK) %>% summarise(sum_harv = sum(HARVEST), .groups = "drop") %>% filter(DISTRICT == d_idx, STAT_WEEK== w_idx) %>% pull(sum_harv) )))),
          .groups = "drop"
        ) %>%
        mutate(
          GR = {if (nchains > 1) {
            coda::gelman.diag(mc_pop[[p_idx]],
                              transform = TRUE,
                              autoburnin = FALSE,
                              multivariate = FALSE)$psrf[,"Point est."] %>%
              .[order(group_names[seq(ncol(p_combo[[p_idx]][[1]])), d_idx])]
          } else {NA}},
          n_eff = coda::effectiveSize(mc_pop[[p_idx]]) %>%
            .[order(group_names[seq(ncol(p_combo[[p_idx]][[1]])), d_idx])]
        ) %>%
        rename(group = name)
    } # W

  } # D

  out <- list()

  out$pop_prop <-
    lapply(p_combo, function(rot) {
      bind_rows(rot) %>%
        mutate(
          chain = rep(1:nchains,
                      each = (nreps - nburn*isFALSE(keep_burn)) / thin),
          itr = rep(1:((nreps - nburn*isFALSE(keep_burn)) / thin),
                    times = nchains)
        )
    }) # add id for chain and iteration

  out$pop_summ <- summ_pop

  out$pop_prop_all <-
    lapply(p_combo_all, function(rot) {
      bind_rows(rot) %>%
        mutate(
          chain = rep(1:nchains,
                      each = (nreps - nburn*isFALSE(keep_burn)) / thin),
          itr = rep(1:((nreps - nburn*isFALSE(keep_burn)) / thin),
                    times = nchains)
        )
    }) # add id for chain and iteration

  out$pop_summ_all <- summ_pop_all

  return(out)

}


format_district_age <- function(outraw, dat_in, nreps, nburn, thin, nchains, keep_burn, C, D, S, W) {

  ### data input ### ----
  metadat <- dat_in$metadat # age and iden info
  harvest <- dat_in$harvest %>% # for calculating p0
    mutate(n = apply(.[,1:4], 1, function(aa) {
      filter(metadat,
             year==aa[1],
             district==aa[2],
             subdist==aa[3],
             week==aa[4]) %>% nrow
    }) ) %>%
    mutate(HARVEST = HARVEST * (n != 0)) %>%
    select(-n)

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

  p_zero <- array(
    apply(
      table(metadat[, c("year","district","subdist","week","iden")], exclude = NULL),
      seq(4),
      function(alf) any(alf!= 0) # prop = 0 if no sample
    ), c(1, D, max(S), W, K+ H))

  if (isFALSE(keep_burn)|keep_burn == "false") keep_burn = FALSE

  ### prepare output ### ----
  # organization of out_list:
  # [[chain]][[yr]][[dist]][[sub]][[week]][age, pop, itr]

  ap_prop <- ap_prop2 <- ap_prop2b <- ap_prop2c <- ap_combo_subdis <- ap_combo_dis <- list()
  ap_prop_all_1 <- ap_prop_all_2 <- list()

  for (d_idx in 1:D) {
    ap_prop[[d_idx]] <- list()
    # ap_prop2b[[d_idx]] <- list()
    # ap_prop2c[[d_idx]] <- list()
    # ap_prop_all_1[[d_idx]] <- ap_prop_all_2[[d_idx]] <- list()
    ap_combo_dis[[d_idx]] <- list()

    idx <- 0
    for (w_idx in 1:W) {
      ap_prop[[d_idx]][[w_idx]] <- list()
      # ap_prop2b[[d_idx]][[w_idx]] <- list()
      # ap_prop2c[[d_idx]][[w_idx]] <- list()

      for (s_idx in 1:S[d_idx]) {

        # go in each chain and grab age-pop output according to specified sampling period
        # new! ap_prop: [[yr]][[d]][[s]][[w]][[chain]][age*itr, c(pop, agevec)]
        # pivot_longer %>% [pop*age*itr, c(itr, age, grp, value)]

        # ap_combo_dis or ap_combo_subdis
        # combine age classes and pop groups
        # weighted by harvest
        # this is done for each week in each district
        # new! [[chain]][[yr]][[d]][[w]][group*itr, c(itr, grpvec, ages)]

        idx <- idx + 1
        ap_prop[[d_idx]][[w_idx]][[s_idx]] <-
          apply(outraw, MARGIN = 7,
                function(ol) ol[ , , w_idx, s_idx, d_idx, 1] %>%
                  as_tibble() %>%
                  pivot_longer(cols = 1:(ncol(.)-2)) %>%
                  mutate(
                    grpvec = rep(group_names[groups[, d_idx], d_idx],
                                 C* (nreps- nburn* isFALSE(keep_burn))/ thin)
                  ) %>%
                  select(-name) %>%
                  rename(itr = 1, agevec = 2))

        ap_combo_dis[[d_idx]][[idx]] <-
          lapply(ap_prop[[d_idx]][[w_idx]][[s_idx]], function(apfile) {
            apfile %>%
              group_by(itr, agevec, grpvec) %>%
              summarise(ppi = sum(value), .groups = "drop") %>%
              mutate(ppi = ppi* {harvest %>%
                  group_by(DISTRICT) %>%
                  # proportional within a district
                  mutate(prop_harv = prop.table(HARVEST)) %>%
                  filter(DISTRICT == d_idx,
                         SUBDISTRICT == s_idx,
                         STAT_WEEK == w_idx) %>%
                  pull(prop_harv) %>% max(., 0)}) %>%
              pivot_wider(names_from = agevec, values_from = ppi)
          })

      } # s_idx
    } # w_idx

    rm(idx)

  } # d_idx

  # exclude burn-in while calculate r hat
  keep_list <- ((nburn*keep_burn + 1):(nreps - nburn*isFALSE(keep_burn)))[!((nburn*keep_burn + 1):(nreps - nburn*isFALSE(keep_burn))) %% thin] / thin

  # place holders/empty objects for age/pop summaries
  # p_combo <- mc_pop <- summ_pop <- list()
  # p_combo_all <- mc_pop_all <- summ_pop_all <- list()
  # p_idx <- 0
  ap_prop_grp <- mc_age <- summ_age <- list()
  a_idx <- 0

  for (d_idx in 1:D) {
    ## age (loop through each district) ##
    ap_combo_all = list()
    for(chain in seq(nchains)) {
      ap_combo_all[[chain]] <-
        lapply(ap_combo_dis[[d_idx]],
               function(apc) {
                 apc[[chain]]
               }) %>%
        bind_rows() %>%
        group_by(itr, grpvec) %>%
        summarise(across(.fns = sum), .groups = "drop") %>%
        pivot_longer(-c(itr, grpvec)) %>%
        group_by(itr, grpvec) %>%
        mutate(value = prop.table(value)) %>%
        ungroup() %>%
        pivot_wider() %>%
        select(-itr) %>%
        setNames(., c("grpvec", age_classes))
    } # combine all individual districts

    for (grp in group_names[, d_idx] %>% na.omit) {
      a_idx <- a_idx + 1

      ap_prop_grp[[a_idx]] <-
        lapply(ap_combo_all, function(aoc) {
          ar_temp <- aoc %>%
            filter(grpvec == grp) %>%
            setNames(., c("grpvec", age_classes))
          return(ar_temp)# %>% select(-grpvec))
        }) # separate rep groups into its own output

      mc_age[[a_idx]] <-
        as.mcmc.list(
          lapply(ap_prop_grp[[a_idx]], function(rlist) {
            mcmc(rlist[keep_list,] %>% select(-grpvec))
          }) )

      harv_dis <- harvest %>%
        filter(DISTRICT == d_idx) %>%
        pull(HARVEST) %>%
        mean

      summ_age[[a_idx]] <-
        lapply(ap_prop_grp[[a_idx]], function(rlist) rlist[keep_list,]) %>%
        bind_rows %>%
        mutate(stock_prop = rowSums(.[,-1])) %>%
        pivot_longer(-c(stock_prop, grpvec)) %>%
        group_by(name, grpvec) %>%
        summarise(
          mean = mean(value),
          median = median(value),
          sd = sd(value),
          ci.05 = quantile(value, 0.05),
          ci.95 = quantile(value, 0.95),
          p0 = mean( value < (0.5/ max(1, harv_dis* stock_prop)) ),
          .groups = "drop"
        ) %>%
        mutate(
          GR = {if (nchains > 1) {
            coda::gelman.diag(mc_age[[a_idx]],
                              transform = TRUE,
                              autoburnin = FALSE,
                              multivariate = FALSE)$psrf[,"Point est."] %>%
              .[order(age_classes)]
          } else {NA}},
          n_eff = coda::effectiveSize(mc_age[[a_idx]]) %>%
            .[order(age_classes)]
        ) %>%
        rename(age = name,
               group = grpvec) %>%
        relocate(group, .before = age)
    } # grp

  } # D

  out <- list()

  out$age_prop <-
    lapply(ap_prop_grp, function(rot) {
      bind_rows(rot) %>%
        mutate(
          chain = rep(1:nchains,
                      each = (nreps - nburn*isFALSE(keep_burn)) / thin),
          itr = rep(1:((nreps - nburn*isFALSE(keep_burn)) / thin),
                    times = nchains)
        )
    }) # add id for chain and iteration

  out$age_summ <- summ_age

  return(out)

}


# subdistrict

format_subdistrict_pop <- function(outraw, dat_in, nreps, nburn, thin, nchains, keep_burn, C, D, S, W) {

  ### data input ### ----
  metadat <- dat_in$metadat # age and iden info
  harvest <- dat_in$harvest %>% # for calculating p0
    mutate(n = apply(.[,1:4], 1, function(aa) {
      filter(metadat,
             year==aa[1],
             district==aa[2],
             subdist==aa[3],
             week==aa[4]) %>% nrow
    }) ) %>%
    mutate(HARVEST = HARVEST * (n != 0)) %>%
    select(-n)

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

  p_zero <- array(
    apply(
      table(metadat[, c("year","district","subdist","week","iden")], exclude = NULL),
      seq(4),
      function(alf) any(alf!= 0) # prop = 0 if no sample
    ), c(1, D, max(S), W, K+ H))

  if (isFALSE(keep_burn)|keep_burn == "false") keep_burn = FALSE

  ### prepare output ### ----
  # organization of out_list:
  # [[chain]][[yr]][[dist]][[sub]][[week]][age, pop, itr]

  ap_prop <- ap_prop2 <- ap_prop2b <- ap_prop2c <- ap_combo_subdis <- list()
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
        # new! ap_prop: [[yr]][[d]][[s]][[w]][[chain]][age*itr, c(pop, agevec)]
        # pivot_longer %>% [pop*age*itr, c(itr, age, grp, value)]

        # ap_combo_dis or ap_combo_subdis
        # combine age classes and pop groups
        # weighted by harvest
        # this is done for each week in each district
        # new! [[chain]][[yr]][[d]][[w]][group*itr, c(itr, grpvec, ages)]

        ap_prop[[d_idx]][[s_idx]][[w_idx]] <-
          apply(outraw, MARGIN = 7,
                function(ol) ol[ , , w_idx, s_idx, d_idx, 1] %>%
                  as_tibble() %>%
                  pivot_longer(cols = 1:(ncol(.)-2)) %>%
                  mutate(
                    grpvec = rep(group_names[groups[, d_idx], d_idx],
                                 C* (nreps- nburn* isFALSE(keep_burn))/ thin)
                  ) %>%
                  select(-name) %>%
                  rename(itr = 1, agevec = 2))

        ap_prop_all_1[[d_idx]][[s_idx]][[w_idx]] <-
          lapply(ap_prop[[d_idx]][[s_idx]][[w_idx]],
                 function(apa) {
                   apa %>%
                     mutate(value = value* {harvest %>%
                         group_by(DISTRICT, SUBDISTRICT) %>%
                         mutate(prop_harv = prop.table(HARVEST)) %>%
                         replace_na(list(prop_harv = 0)) %>%
                         filter(DISTRICT == d_idx,
                                SUBDISTRICT == s_idx,
                                STAT_WEEK == w_idx) %>%
                         pull(prop_harv) %>% max(., 0)})
                 })

      } # w_idx

      for (chain in seq(nchains)) {
        ap_prop_all_2[[d_idx]][[s_idx]][[chain]] <-
          lapply(ap_prop_all_1[[d_idx]][[s_idx]],
                 function(apa) {
                   apa[[chain]]
                 }) %>% bind_rows()
      } # chain

    } # s_idx

    rm(idx)

  } # d_idx

  # exclude burn-in while calculate r hat
  keep_list <- ((nburn*keep_burn + 1):(nreps - nburn*isFALSE(keep_burn)))[!((nburn*keep_burn + 1):(nreps - nburn*isFALSE(keep_burn))) %% thin] / thin

  # place holders/empty objects for age/pop summaries
  p_combo <- mc_pop <- summ_pop <- list()
  p_combo_all <- mc_pop_all <- summ_pop_all <- list()
  p_idx <- 0

  for (d_idx in 1:D) {
    for (s_idx in 1:S[d_idx]) {

      ## pops all dist ##
      # sum over ages and combine pop groups

      p_combo_all[[max(S)*(d_idx-1)+s_idx]] <-
        lapply(ap_prop_all_2[[d_idx]][[s_idx]],
               function(apchain) {
                 apchain %>%
                   select(-agevec) %>%
                   group_by(itr, grpvec) %>%
                   summarise(p = sum(value), .groups = "drop") %>%
                   pivot_wider(names_from = grpvec, values_from = p) %>%
                   select(-itr)
               })

      mc_pop_all[[max(S)*(d_idx-1)+s_idx]] <- as.mcmc.list(
        lapply(p_combo_all[[max(S)*(d_idx-1)+s_idx]],
               function(rlist) mcmc(rlist[keep_list,])) )

      summ_pop_all[[max(S)*(d_idx-1)+s_idx]] <-
        lapply(p_combo_all[[max(S)*(d_idx-1)+s_idx]], function(rlist) rlist[keep_list,]) %>%
        bind_rows %>%
        pivot_longer(cols = 1:ncol(.)) %>%
        group_by(name) %>%
        summarise(
          mean = mean(value),
          median = median(value),
          sd = sd(value),
          ci.05 = quantile(value, 0.05),
          ci.95 = quantile(value, 0.95),
          p0 = mean(value < (0.5/ max(1, ( any(p_zero[1, d_idx, , , ])* harvest %>% group_by(DISTRICT, SUBDISTRICT) %>% summarise(sum_harv = sum(HARVEST), .groups = "drop") %>% filter(DISTRICT == d_idx, SUBDISTRICT == s_idx) %>% pull(sum_harv) )))),
          .groups = "drop"
        ) %>%
        mutate(
          GR = {if (nchains > 1) {
            coda::gelman.diag(mc_pop_all[[max(S)*(d_idx-1)+s_idx]],
                              transform = TRUE,
                              autoburnin = FALSE,
                              multivariate = FALSE)$psrf[,"Point est."] %>%
              .[order(group_names[seq(ncol(p_combo_all[[max(S)*(d_idx-1)+s_idx]][[1]])), d_idx])]
          } else {NA}},
          n_eff = coda::effectiveSize(mc_pop_all[[max(S)*(d_idx-1)+s_idx]]) %>%
            .[order(group_names[seq(ncol(p_combo_all[[max(S)*(d_idx-1)+s_idx]][[1]])), d_idx])]
        ) %>%
        rename(group = name)

      for (w_idx in 1:W) {

        ## indivi pop ##
        # sum over ages and combine pop groups
        p_idx <- p_idx + 1 # index individual sampling period

        p_combo[[p_idx]] <-
          lapply(ap_prop[[d_idx]][[s_idx]][[w_idx]],
                 function(apfile) {
                   apfile %>%
                     select(-agevec) %>%
                     group_by(itr, grpvec) %>%
                     summarise(p = sum(value), .groups = "drop") %>%
                     pivot_wider(names_from = grpvec, values_from = p) %>%
                     select(-itr) #%>%
                   # setNames(group_names[seq(ncol(.)), d_idx])
                 })

        mc_pop[[p_idx]] <- as.mcmc.list(
          lapply(p_combo[[p_idx]],
                 function(rlist) mcmc(rlist[keep_list,])) )

        summ_pop[[p_idx]] <-
          lapply(p_combo[[p_idx]], function(rlist) rlist[keep_list,]) %>%
          bind_rows %>%
          pivot_longer(cols = 1:ncol(.)) %>%
          group_by(name) %>%
          summarise(
            mean = mean(value),
            median = median(value),
            sd = sd(value),
            ci.05 = quantile(value, 0.05),
            ci.95 = quantile(value, 0.95),
            p0 = mean(value < (0.5/ max(1, (any(p_zero[1, d_idx, s_idx, w_idx, ])* filter(harvest, DISTRICT== d_idx, SUBDISTRICT== s_idx, STAT_WEEK== w_idx)$HARVEST)))),
            .groups = "drop"
          ) %>%
          mutate(
            GR = {if (nchains > 1) {
              coda::gelman.diag(mc_pop[[p_idx]],
                                transform = TRUE,
                                autoburnin = FALSE,
                                multivariate = FALSE)$psrf[,"Point est."] %>%
                .[order(group_names[seq(ncol(p_combo[[p_idx]][[1]])), d_idx])]
            } else {NA}},
            n_eff = coda::effectiveSize(mc_pop[[p_idx]]) %>%
              .[order(group_names[seq(ncol(p_combo[[p_idx]][[1]])), d_idx])]
          ) %>%
          rename(group = name)

      } # W

    } # S
  } # D

  out <- list()

  out$pop_prop <-
    lapply(p_combo, function(rot) {
      bind_rows(rot) %>%
        mutate(
          chain = rep(1:nchains,
                      each = (nreps - nburn*isFALSE(keep_burn)) / thin),
          itr = rep(1:((nreps - nburn*isFALSE(keep_burn)) / thin),
                    times = nchains)
        )
    }) # add id for chain and iteration

  out$pop_summ <- summ_pop

  out$pop_prop_all <-
    lapply(p_combo_all, function(rot) {
      bind_rows(rot) %>%
        mutate(
          chain = rep(1:nchains,
                      each = (nreps - nburn*isFALSE(keep_burn)) / thin),
          itr = rep(1:((nreps - nburn*isFALSE(keep_burn)) / thin),
                    times = nchains)
        )
    }) # add id for chain and iteration

  out$pop_summ_all <- summ_pop_all

  return(out)

}


format_subdistrict_age <- function(outraw, dat_in, nreps, nburn, thin, nchains, keep_burn, C, D, S, W) {

  ### data input ### ----
  metadat <- dat_in$metadat # age and iden info
  harvest <- dat_in$harvest %>% # for calculating p0
    mutate(n = apply(.[,1:4], 1, function(aa) {
      filter(metadat,
             year==aa[1],
             district==aa[2],
             subdist==aa[3],
             week==aa[4]) %>% nrow
    }) ) %>%
    mutate(HARVEST = HARVEST * (n != 0)) %>%
    select(-n)

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

  p_zero <- array(
    apply(
      table(metadat[, c("year","district","subdist","week","iden")], exclude = NULL),
      seq(4),
      function(alf) any(alf!= 0) # prop = 0 if no sample
    ), c(1, D, max(S), W, K+ H))

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
          apply(outraw, MARGIN = 7,
                function(ol) ol[ , , w_idx, s_idx, d_idx, 1] %>%
                  as_tibble() %>%
                  pivot_longer(cols = 1:(ncol(.)-2)) %>%
                  mutate(
                    grpvec = rep(group_names[groups[, d_idx], d_idx],
                                 C* (nreps- nburn* isFALSE(keep_burn))/ thin)
                  ) %>%
                  select(-name) %>%
                  rename(itr = 1, agevec = 2))

        ap_combo_subdis[[d_idx]][[s_idx]][[w_idx]] <-
          lapply(ap_prop[[d_idx]][[s_idx]][[w_idx]],
                 function(apfile) {

                   apfile %>%
                     group_by(itr, agevec, grpvec) %>%
                     summarise(ppi = sum(value), .groups = "drop") %>%
                     mutate(ppi = ppi* {harvest %>%
                         group_by(DISTRICT, SUBDISTRICT) %>%
                         # proportional within a subdistrict
                         mutate(prop_harv = prop.table(HARVEST)) %>%
                         filter(DISTRICT == d_idx,
                                SUBDISTRICT == s_idx,
                                STAT_WEEK == w_idx) %>%
                         pull(prop_harv) %>% max(., 0)}) %>%
                     pivot_wider(names_from = agevec, values_from = ppi)

                 })

      } # w_idx

    } # s_idx

    rm(idx)

  } # d_idx

  # exclude burn-in while calculate r hat
  keep_list <- ((nburn*keep_burn + 1):(nreps - nburn*isFALSE(keep_burn)))[!((nburn*keep_burn + 1):(nreps - nburn*isFALSE(keep_burn))) %% thin] / thin

  # place holders/empty objects for age/pop summaries
  ap_prop_grp <- mc_age <- summ_age <- list()
  a_idx <- 0

  for (d_idx in 1:D) {
    for (s_idx in 1:S[d_idx]) {

      ## age (loop through each subdistrict) ##

      ap_combo_all = list()
      for(chain in seq(nchains)) {
        ap_combo_all[[chain]] <-
          lapply(ap_combo_subdis[[d_idx]][[s_idx]],
                 function(apc) {
                   apc[[chain]]
                 }) %>%
          bind_rows() %>%
          group_by(itr, grpvec) %>%
          summarise(across(.fns = sum), .groups = "drop") %>%
          pivot_longer(-c(itr, grpvec)) %>%
          group_by(itr, grpvec) %>%
          mutate(value = prop.table(value)) %>%
          ungroup() %>%
          pivot_wider() %>%
          select(-itr) %>%
          setNames(., c("grpvec", age_classes))
      } # combine all individual districts

      for (grp in group_names[, d_idx] %>% na.omit) {
        a_idx <- a_idx + 1

        ap_prop_grp[[a_idx]] <-
          lapply(ap_combo_all, function(aoc) {
            ar_temp <- aoc %>%
              filter(grpvec == grp) %>%
              setNames(., c("grpvec", age_classes))
            return(ar_temp)# %>% select(-grpvec))
          }) # separate rep groups into its own output

        mc_age[[a_idx]] <-
          as.mcmc.list(
            lapply(ap_prop_grp[[a_idx]], function(rlist) {
              mcmc(rlist[keep_list,] %>% select(-grpvec))
            }) )

        harv_subdis <- harvest %>%
          filter(DISTRICT == d_idx, SUBDISTRICT == s_idx) %>%
          pull(HARVEST) %>%
          mean

        summ_age[[a_idx]] <-
          lapply(ap_prop_grp[[a_idx]], function(rlist) rlist[keep_list,]) %>%
          bind_rows %>%
          mutate(stock_prop = rowSums(.[,-1])) %>%
          pivot_longer(-c(stock_prop, grpvec)) %>%
          group_by(name, grpvec) %>%
          summarise(
            mean = mean(value),
            median = median(value),
            sd = sd(value),
            ci.05 = quantile(value, 0.05),
            ci.95 = quantile(value, 0.95),
            p0 = mean( value < (0.5/ max(1, harv_subdis* stock_prop)) ),
            .groups = "drop"
          ) %>%
          mutate(
            GR = {if (nchains > 1) {
              coda::gelman.diag(mc_age[[a_idx]],
                                transform = TRUE,
                                autoburnin = FALSE,
                                multivariate = FALSE)$psrf[,"Point est."] %>%
                .[order(age_classes)]
            } else {NA}},
            n_eff = coda::effectiveSize(mc_age[[a_idx]]) %>%
              .[order(age_classes)]
          ) %>%
          rename(age = name,
                 group = grpvec) %>%
          relocate(group, .before = age)
      } # grp

    } # S
  } # D

  out <- list()

  out$age_prop <-
    lapply(ap_prop_grp, function(rot) {
      bind_rows(rot) %>%
        mutate(
          chain = rep(1:nchains,
                      each = (nreps - nburn*isFALSE(keep_burn)) / thin),
          itr = rep(1:((nreps - nburn*isFALSE(keep_burn)) / thin),
                    times = nchains)
        )
    }) # add id for chain and iteration

  out$age_summ <- summ_age

  return(out)

}


# runners ----
magma_format_pop <- function(outraw, dat_in, nreps, nburn, thin, nchains, keep_burn, summ_level) {

  ### load required packages ### ----
  if(!require("pacman")) install.packages("pacman")
  library(pacman)
  pacman::p_load(coda, tidyverse)

  ### info needed ### ----
  C <- dat_in$C # number of age classes
  D <- dplyr::n_distinct(dat_in$metadat$district) # number of districts
  S <- dat_in$metadat %>%
    group_by(district) %>%
    summarise(S = n_distinct(subdist),
              .groups = "drop") %>%
    pull(S) # number of subdistricts
  W <- max(dplyr::n_distinct(dat_in$metadat$week), length(dat_in$stat_weeks)) # number of stats weeks

  ng <- apply(dat_in$groups, 2, max) # number of groups

  if (!is.null(dat_in$districts)) dist_names <- dat_in$districts
  if (!is.null(dat_in$subdistricts)) subdist_names <- dat_in$subdistricts
  if (!is.null(dat_in$stat_weeks)) week_names <- dat_in$stat_weeks

  ### prepare output ### ----

  message("Preparing output (patience grasshopper...)")
  prep_time <- Sys.time()

  if (summ_level == "district") {
    out <- format_district_pop(outraw, dat_in, nreps, nburn, thin, nchains, keep_burn, C, D, S, W)
  } else if (summ_level == "subdistrict") {
    out <- format_subdistrict_pop(outraw, dat_in, nreps, nburn, thin, nchains, keep_burn, C, D, S, W)
  }

  if (summ_level == "district") {

    pout_names <- rep(NA, D*W)
    j <- 1
    for (d_i in 1:D) {
      for (w_i in 1:W) {
        pout_names[j] <-
          if (all(exists("dist_names"), exists("week_names"))) {
            paste(paste0("(", j, ")"),
                  paste0("D", dist_names[d_i]),
                  paste0("Wk", week_names[w_i]),
                  sep = "_")
          } else {
            paste(paste0("(", j, ")"),
                  paste0("D", d_i),
                  paste0("Wk", w_i),
                  sep = "_")
          }
        j <- j + 1
      } # w_i
    } # d_i

  } else if (summ_level == "subdistrict") {

    pout_names <- rep(NA, sum(S*W))
    j <- 1
    for (d_i in 1:D) {
      for (s_i in 1:S[d_i]) {
        for (w_i in 1:W) {
          pout_names[j] <-
            if (all(exists("dist_names"), exists("subdist_names"), exists("week_names"))) {
              paste(paste0("(", j, ")"),
                    paste0("D", dist_names[d_i]),
                    paste0("Sub", subdist_names[[d_i]][s_i]),
                    paste0("Wk", week_names[w_i]),
                    sep = "_")
            } else {
              paste(paste0("(", j, ")"),
                    paste0("D", d_i),
                    paste0("Sub", s_i),
                    paste0("Wk", w_i),
                    sep = "_")
            }
          j <- j + 1
        } # w_i
      } # s_i
    } # d_i

  } # end if dist/subdist

  for(li in 1:2) {
    names(out[[li]]) <- pout_names
  } # id for pop output

  if (summ_level == "district") {
    for(li in 3:4) {
      names(out[[li]]) <- dist_names
    }
  } else if (summ_level == "subdistrict") {
    for(li in 3:4) {
      names(out[[li]]) <- subdist_names %>% unlist() %>% names
    }
  }

  print(Sys.time() - prep_time)
  message(Sys.time())

  return(out)

}


magma_format_age <- function(outraw, dat_in, nreps, nburn, thin, nchains, keep_burn, summ_level) {

  ### load required packages ### ----
  if(!require("pacman")) install.packages("pacman")
  library(pacman)
  pacman::p_load(coda, tidyverse)

  ### info needed ### ----
  C <- dat_in$C # number of age classes
  D <- dplyr::n_distinct(dat_in$metadat$district) # number of districts
  S <- dat_in$metadat %>%
    group_by(district) %>%
    summarise(S = n_distinct(subdist),
              .groups = "drop") %>%
    pull(S) # number of subdistricts
  W <- max(dplyr::n_distinct(dat_in$metadat$week), length(dat_in$stat_weeks)) # number of stats weeks

  ng <- apply(dat_in$groups, 2, max) # number of groups

  if (!is.null(dat_in$districts)) dist_names <- dat_in$districts
  if (!is.null(dat_in$subdistricts)) subdist_names <- dat_in$subdistricts
  if (!is.null(dat_in$stat_weeks)) week_names <- dat_in$stat_weeks

  ### prepare output ### ----

  message("Preparing output (patience grasshopper...)")
  prep_time <- Sys.time()

  if (summ_level == "district") {
    out <- format_district_age(outraw, dat_in, nreps, nburn, thin, nchains, keep_burn, C, D, S, W)
  } else if (summ_level == "subdistrict") {
    out <- format_subdistrict_age(outraw, dat_in, nreps, nburn, thin, nchains, keep_burn, C, D, S, W)
  }

  # id for age outputs
  if (summ_level == "district") {

    aout_names <- rep(NA, sum(ng))
    i <- 0
    for (d_i in 1:D) {
      for (g_i in 1:ng[d_i]) {
        i <- i + 1
        aout_names[i] <-
          if (exists("dist_names")) {
            paste(paste0("(", i, ")"),
                  paste0("D", dist_names[d_i]),
                  paste0(dat_in$group_names[g_i, d_i]),
                  sep = "_")
          } else {
            paste(paste0("(", i, ")"),
                  paste0("D", d_i),
                  paste0(dat_in$group_names[g_i, d_i]),
                  sep = "_")
          }
      } # g_i
    } # d_i

  } else if (summ_level == "subdistrict") {

    aout_names <- rep(NA, sum(ng* S))
    i <- 0
    for (d_i in 1:D) {
      for(s_i in 1:S[d_i]) {
        for (g_i in 1:ng[d_i]) {
          i <- i + 1
          aout_names[i] <-
            if (all(exists("dist_names"), exists("subdist_names"))) {
              paste(paste0("(", i, ")"),
                    paste0("D", dist_names[d_i]),
                    paste0("SubD", subdist_names[[d_i]][s_i]),
                    paste0(dat_in$group_names[g_i, d_i]),
                    sep = "_")
            } else {
              paste(paste0("(", i, ")"),
                    paste0("D", d_i),
                    paste0("SubD", s_i),
                    paste0(dat_in$group_names[g_i, d_i]),
                    sep = "_")
            }
        } # g_i
      } # s_i
    } # d_i

  } # end if

  for(li in 1:2) {
    names(out[[li]]) <- aout_names
  } # id for pop output

  print(Sys.time() - prep_time)
  message(Sys.time())

  return(out)

}


## format wrapper for TBR ----
tbr_format_wraper <- function(which_dist, outraw, ma_dat, n_itr, n_burn, n_thin, n_chain, keepburn, summlevel, type = NULL) {

  sub_dat <- list(
    x = ma_dat$x[which(ma_dat$metadat$district == which_dist), ],
    y = ma_dat$y,

    metadat = ma_dat$metadat %>%
      filter(district == which_dist) %>%
      mutate(district = 1), # chase
    harvest = ma_dat$harvest %>%
      filter(DISTRICT == which_dist) %>%
      mutate(DISTRICT = 1), # chase

    nstates = ma_dat$nstates,
    nalleles = ma_dat$nalleles,
    C = ma_dat$C,
    groups = cbind(ma_dat$groups[, which_dist, drop = FALSE]),
    group_names = cbind(ma_dat$group_names[, which_dist, drop = FALSE]),

    age_class = ma_dat$age_class,
    age_classes = ma_dat$age_classes,
    wildpops = ma_dat$wildpops,
    hatcheries = ma_dat$hatcheries,

    districts = ma_dat$districts[which_dist, drop = FALSE],
    subdistricts = ma_dat$subdistricts[which_dist, drop = FALSE],
    stat_weeks = ma_dat$stat_weeks
  )

  S <- ma_dat$metadat %>%
    group_by(district) %>%
    summarise(S = n_distinct(subdist),
              .groups = "drop") %>%
    pull(S) # number of subdistricts
  W <- dplyr::n_distinct(ma_dat$metadat$week) # number of weeks
  KH <- nrow(ma_dat$y)

  holder <-
    lapply(outraw, function(dis) {
      lapply(dis[[which_dist]], function(subdis) {
        lapply(subdis, function(wk) {
          wk %>% bind_rows() %>% unlist()
        }) %>% unlist()
      }) %>% unlist()
    })  %>% unlist %>%
    array(., dim = c((ma_dat$C)* (n_itr- n_burn)/ n_thin, KH+2, W, S[which_dist], 1, 1, n_chain))

  if (is.null(type)) {
    sub_out <- magma_format_out(holder, sub_dat, n_itr, n_burn, n_thin, n_chain, keepburn, summlevel)
  } else if (type == "pop") {
    sub_out <- magma_format_pop(holder, sub_dat, n_itr, n_burn, n_thin, n_chain, keepburn, summlevel)
  } else if (type == "age") {
    sub_out <- magma_format_age(holder, sub_dat, n_itr, n_burn, n_thin, n_chain, keepburn, summlevel)
  }

  return(sub_out)

  # wraper to format raw magma output one district at a time

  # Example:
  #
  # nreps = 100
  # nburn = 50
  # thin = 5
  # nchains = 5
  # keep_burn = TRUE
  #
  # out2 <- format_wraper(which_dist = 2,
  #                       outraw = out_raw,
  #                       ma_dat = magmadat,
  #                       n_itr = nreps,
  #                       n_burn = nburn,
  #                       n_thin = thin,
  #                       n_chain = nchains,
  #                       keepburn = keep_burn,
  #                       summlevel = "district",
  #                       type = "pop")

  ## if you don't specify "type", it will summarize both pop and age at the same time

}

# format_wraper <- function(which_dist, outraw, ma_dat, n_itr, n_burn, n_thin, n_chain, keepburn, summlevel, type = NULL) {
#
#   sub_dat <- list(
#     x = ma_dat$x[which(ma_dat$metadat$district == which_dist), ],
#     y = ma_dat$y,
#
#     metadat = ma_dat$metadat %>%
#       filter(district == which_dist) %>%
#       mutate(district = 1), # chase
#     harvest = ma_dat$harvest %>%
#       filter(DISTRICT == which_dist) %>%
#       mutate(DISTRICT = 1), # chase
#
#     nstates = ma_dat$nstates,
#     nalleles = ma_dat$nalleles,
#     C = ma_dat$C,
#     groups = cbind(ma_dat$groups[, which_dist, drop = FALSE]),
#     group_names = cbind(ma_dat$group_names[, which_dist, drop = FALSE]),
#
#     age_class = ma_dat$age_class,
#     age_classes = ma_dat$age_classes,
#     wildpops = ma_dat$wildpops,
#     hatcheries = ma_dat$hatcheries,
#
#     districts = ma_dat$districts[which_dist, drop = FALSE],
#     subdistricts = ma_dat$subdistricts[which_dist, drop = FALSE],
#     stat_weeks = ma_dat$stat_weeks
#   )
#
#   holder <- array(NA, c(dim(outraw)[1:4], 1, 1, n_chain))
#   holder[seq(dim(outraw)[1]), seq(dim(outraw)[2]), seq(dim(outraw)[3]), seq(dim(outraw)[4]), 1, 1, seq(n_chain)] <- outraw[seq(dim(outraw)[1]), seq(dim(outraw)[2]), seq(dim(outraw)[3]), seq(dim(outraw)[4]), which_dist, 1, seq(n_chain)]
#
#   if (is.null(type)) {
#     sub_out <- magma_format_out(holder, sub_dat, n_itr, n_burn, n_thin, n_chain, keepburn, summlevel)
#   } else if (type == "pop") {
#     sub_out <- magma_format_pop(holder, sub_dat, n_itr, n_burn, n_thin, n_chain, keepburn, summlevel)
#   } else if (type == "age") {
#     sub_out <- magma_format_age(holder, sub_dat, n_itr, n_burn, n_thin, n_chain, keepburn, summlevel)
#   }
#
#   return(sub_out)
#
#   # wraper to format raw magma output one district at a time
#
#   # Example:
#   #
#   # nreps = 100
#   # nburn = 50
#   # thin = 5
#   # nchains = 5
#   # keep_burn = TRUE
#   #
#   # out2 <- format_wraper(which_dist = 2,
#   #                       outraw = out_raw,
#   #                       ma_dat = magmadat,
#   #                       n_itr = nreps,
#   #                       n_burn = nburn,
#   #                       n_thin = thin,
#   #                       n_chain = nchains,
#   #                       keepburn = keep_burn,
#   #                       summlevel = "district",
#   #                       type = "pop")
#
#   ## if you don't specify "type", it will summarize both pop and age at the same time
#
# }


## traceplot ----
ggtr_plt <- function (obj, nburn = 0, thin = 1) {

  if ("grpvec" %in% names(obj)) obj <- select(obj, -grpvec)

  pivot_longer({{obj}}, cols = -c(chain, itr)) %>%
    ggplot() +
    geom_line(aes(x = itr, y = value, color = factor(chain))) +
    {if (nburn > 0) {
      annotate("rect", fill = "red", alpha = 0.15,
               xmin = 0, xmax = nburn/thin, ymin = -Inf, ymax = Inf)
    }} +
    facet_grid(name ~ ., scales = "free") +
    labs(color = "MC chain")

} # nburn = 0 if keep_burn = FALSE













