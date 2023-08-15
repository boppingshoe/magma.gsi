#' Summarize big daddy model output ste 3
#'
#' @param out2 Output from step 2
#'
#' @return Model output in process and metadata as a list object.
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' # format data
#' wd <- getwd() # path to data folder
#' magma_data <- magmatize_data(wd = wd, save_data = FALSE)
#'
#' # model run
#' magma_out <- magmatize_mdl(magma_data,
#'   nreps = 50, nburn = 25, thin = 1, nchains = 2)
#'
#' # summary steps 1 to 3
#' tbr1 <- magmatize_summ_bd1(which_dist = 1,
#'   ma_out = magma_out,
#'   ma_dat = magma_data)
#'
#' tbr2 <- magmatize_summ_bd2(tbr1, summ_level = "district", type = "pop")
#'
#' tbr3 <- magmatize_summ_bd3(tbr2)
#' }
#'
#' @export
magmatize_summ_bd3 <- function(out2) {

  nreps <- out2$sub_dat$nreps
  nburn <- out2$sub_dat$nburn
  thin <- out2$sub_dat$thin
  nchains <- out2$sub_dat$nchains
  keep_burn <- out2$sub_dat$keep_burn
  summ_level <- out2$sub_dat$summ_level
  type <- out2$sub_dat$type

  if (type == "pop") {
    sub_out_b <- magmatize_pop_b(out2$holder, out2$sub_dat, nreps, nburn, thin, nchains, keep_burn, summ_level)
  } else if (type == "age") {
    sub_out_b <- magmatize_age_b(out2$holder, out2$sub_dat, nreps, nburn, thin, nchains, keep_burn, summ_level)
  }

  return(sub_out_b)

}

# pop ----

#' Summarize population part of MAGMA output (part b)
#'
#' @param outraw MAGMA output
#' @param dat_in MAGMA data
#' @param nreps same as MAGMA model run
#' @param nburn same as MAGMA model run
#' @param thin same as MAGMA model run
#' @param nchains same as MAGMA model run
#' @param keep_burn same as MAGMA model run
#' @param summ_level Summarize to district or subdistrict
#'
#' @importFrom magrittr %>%
#'
#' @noRd
#'
magmatize_pop_b <- function(out2, dat_in, nreps, nburn, thin, nchains, keep_burn, summ_level) {

  ### info needed ### ----
  C <- dat_in$C # number of age classes
  D <- dplyr::n_distinct(dat_in$metadat$district) # number of districts
  S <- dat_in$metadat %>%
    dplyr::group_by(district) %>%
    dplyr::summarise(S = dplyr::n_distinct(subdist), .groups = "drop") %>%
    dplyr::pull(S) # number of subdistricts
  W <- max(dplyr::n_distinct(dat_in$metadat$week), length(dat_in$stat_weeks)) # number of stats weeks

  if (!is.null(dat_in$districts)) dist_names <- dat_in$districts
  if (!is.null(dat_in$subdistricts)) subdist_names <- dat_in$subdistricts
  if (!is.null(dat_in$stat_weeks)) week_names <- dat_in$stat_weeks

  ### prepare output ### ----

  if (summ_level == "district") {
    out <- format_district_pop_b(out2, dat_in, nreps, nburn, thin, nchains, keep_burn, C, D, S, W)
  } else if (summ_level == "subdistrict") {
    out <- format_subdistrict_pop_b(out2, dat_in, nreps, nburn, thin, nchains, keep_burn, C, D, S, W)
  } else stop("Invalid summ_level.")

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
                    paste0("SubD", subdist_names[[d_i]][s_i]),
                    paste0("Wk", week_names[w_i]),
                    sep = "_")
            } else {
              paste(paste0("(", j, ")"),
                    paste0("D", d_i),
                    paste0("SubD", s_i),
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
      names(out[[li]]) <- subdist_names %>% unlist() %>% names()
    }
  }

  return(out)

}

#' Format district populations step b
#'
#' @param popout_a Output from staep a
#' @param dat_in Subset of MAGMA data
#' @param nreps Same nreps as teh model run
#' @param nburn Same nburn as teh model run
#' @param thin Same thin as teh model run
#' @param nchains Same nchains as teh model run
#' @param keep_burn Same keep_burn as teh model run
#' @param C Age classes
#' @param D Districts
#' @param S Subdistricts
#' @param W Weeks
#'
#' @importFrom magrittr %>%
#'
#' @noRd
#'
format_district_pop_b <- function(popout_a, dat_in, nreps, nburn, thin, nchains, keep_burn, C, D, S, W) {

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

  # groups <- dat_in$groups # vector id for reporting groups (aka groupvec)
  group_names <- dat_in$group_names # reporting groups

  wildpops <- dat_in$wildpops
  K <- length(wildpops)

  H <- ifelse(is.null(dat_in$hatcheries), 0, length(dat_in$hatcheries))

  p_zero <- array(
    apply(
      table(metadat[, c("district", "subdist", "week", "iden")], exclude = NULL),
      seq(3),
      function(alf) any(alf!= 0) # prop = 0 if no sample
    ), c(D, max(S), W, K+ H))

  if (isFALSE(keep_burn) | keep_burn == "false") keep_burn = FALSE

  # exclude burn-in while calculate r hat
  keep_list <- ((nburn*keep_burn + 1):(nreps - nburn*isFALSE(keep_burn)))[!((nburn*keep_burn + 1):(nreps - nburn*isFALSE(keep_burn))) %% thin] / thin

  # place holders/empty objects for pop summaries
  p_combo <- mc_pop <- summ_pop <- list()
  p_combo_all <- mc_pop_all <- summ_pop_all <- list()
  p_idx <- 0

  for (d_idx in 1:D) {
    ## pops all dist ##
    # sum over ages and combine pop groups

    p_combo_all[[d_idx]] <-
      lapply(popout_a$ap_prop_all_2[[d_idx]],
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
      dplyr::bind_rows() %>%
      tidyr::pivot_longer(cols = dplyr::everything()) %>%
      dplyr::group_by(name) %>%
      dplyr::summarise(
        mean = mean(value),
        median = stats::median(value),
        sd = stats::sd(value),
        ci.05 = stats::quantile(value, 0.05),
        ci.95 = stats::quantile(value, 0.95),
        p0 = mean(value < (0.5/ max(1, ( any(p_zero[d_idx, , , ])* harvest %>% dplyr::group_by(DISTRICT) %>% dplyr::summarise(sum_harv = sum(HARVEST), .groups = "drop") %>% dplyr::filter(DISTRICT == d_idx) %>% dplyr::pull(sum_harv) )))),
        .groups = "drop"
      ) %>%
      dplyr::mutate(name = factor(name, levels = c(group_names[seq(ncol(p_combo_all[[d_idx]][[1]])), d_idx]))) %>%
      dplyr::arrange(name) %>%
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
      ## pop individual weeks ##
      # sum over ages and combine pop groups
      p_idx <- p_idx + 1 # index individual sampling period

      p_combo[[p_idx]] <-
        lapply(popout_a$ap_prop2c[[d_idx]][[w_idx]],
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
        dplyr::bind_rows() %>%
        tidyr::pivot_longer(cols = dplyr::everything()) %>%
        dplyr::group_by(name) %>%
        dplyr::summarise(
          mean = mean(value),
          median = stats::median(value),
          sd = stats::sd(value),
          ci.05 = stats::quantile(value, 0.05),
          ci.95 = stats::quantile(value, 0.95),
          p0 = mean(value < (0.5/ max(1, ( any(p_zero[d_idx, , w_idx, ])* harvest %>% dplyr::group_by(DISTRICT, STAT_WEEK) %>% dplyr::summarise(sum_harv = sum(HARVEST), .groups = "drop") %>% dplyr::filter(DISTRICT == d_idx, STAT_WEEK== w_idx) %>% dplyr::pull(sum_harv) )))),
          .groups = "drop"
        ) %>%
        dplyr::mutate(name = factor(name, levels = c(group_names[seq(ncol(p_combo_all[[d_idx]][[1]])), d_idx]))) %>%
        dplyr::arrange(name) %>%
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

  } # D

  out <- list()

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

#' Format subdistrict populations step b
#'
#' @param popout_a Output from staep a
#' @param dat_in Subset of MAGMA data
#' @param nreps Same nreps as teh model run
#' @param nburn Same nburn as teh model run
#' @param thin Same thin as teh model run
#' @param nchains Same nchains as teh model run
#' @param keep_burn Same keep_burn as teh model run
#' @param C Age classes
#' @param D Districts
#' @param S Subdistricts
#' @param W Weeks
#'
#' @importFrom magrittr %>%
#'
#' @noRd
#'
format_subdistrict_pop_b <- function(popout_a, dat_in, nreps, nburn, thin, nchains, keep_burn, C, D, S, W) {

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

  # groups <- dat_in$groups # vector id for reporting groups (aka groupvec)
  group_names <- dat_in$group_names # reporting groups

  wildpops <- dat_in$wildpops
  K <- length(wildpops)

  H <- ifelse(is.null(dat_in$hatcheries), 0, length(dat_in$hatcheries))

  p_zero <- array(
    apply(
      table(metadat[, c("district", "subdist", "week", "iden")], exclude = NULL),
      seq(3),
      function(alf) any(alf!= 0) # prop = 0 if no sample
    ), c(D, max(S), W, K+ H))

  if (isFALSE(keep_burn) | keep_burn == "false") keep_burn = FALSE

  # exclude burn-in while calculate r hat
  keep_list <- ((nburn*keep_burn + 1):(nreps - nburn*isFALSE(keep_burn)))[!((nburn*keep_burn + 1):(nreps - nburn*isFALSE(keep_burn))) %% thin] / thin

  ### prepare output ### ----
  # organization of out_list:
  # [[chain]][[yr]][[dist]][[sub]][[week]][age, pop, itr]

  # place holders/empty objects for age/pop summaries
  p_combo <- mc_pop <- summ_pop <- list()
  p_combo_all <- mc_pop_all <- summ_pop_all <- list()
  p_idx <- 0

  for (d_idx in 1:D) {
    for (s_idx in 1:S[d_idx]) {

      ## pops all dist ##
      # sum over ages and combine pop groups

      p_combo_all[[max(S)*(d_idx-1)+s_idx]] <-
        lapply(popout_a$ap_prop_all_2[[d_idx]][[s_idx]],
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
        dplyr::bind_rows() %>%
        tidyr::pivot_longer(cols = dplyr::everything()) %>%
        dplyr::group_by(name) %>%
        dplyr::summarise(
          mean = mean(value),
          median = stats::median(value),
          sd = stats::sd(value),
          ci.05 = stats::quantile(value, 0.05),
          ci.95 = stats::quantile(value, 0.95),
          p0 = mean(value < (0.5/ max(1, ( any(p_zero[d_idx, , , ])* harvest %>% dplyr::group_by(DISTRICT, SUBDISTRICT) %>% dplyr::summarise(sum_harv = sum(HARVEST), .groups = "drop") %>% dplyr::filter(DISTRICT == d_idx, SUBDISTRICT == s_idx) %>% dplyr::pull(sum_harv) )))),
          .groups = "drop"
        ) %>%
        dplyr::mutate(name = factor(name, levels = c(group_names[seq(ncol(p_combo_all[[max(S)*(d_idx-1)+s_idx]][[1]])), d_idx]))) %>%
        dplyr::arrange(name) %>%
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

        ## indivi pop ##
        # sum over ages and combine pop groups
        p_idx <- p_idx + 1 # index individual sampling period

        p_combo[[p_idx]] <-
          lapply(popout_a$ap_prop[[d_idx]][[s_idx]][[w_idx]],
                 function(apfile) {
                   apfile %>%
                     dplyr::select(-agevec) %>%
                     dplyr::group_by(itr, grpvec) %>%
                     dplyr::summarise(p = sum(value), .groups = "drop") %>%
                     tidyr::pivot_wider(names_from = grpvec, values_from = p) %>%
                     dplyr::select(-itr) #%>%
                   # setNames(group_names[seq(ncol(.)), d_idx])
                 })

        mc_pop[[p_idx]] <- coda::as.mcmc.list(
          lapply(p_combo[[p_idx]],
                 function(rlist) coda::mcmc(rlist[keep_list,])) )

        summ_pop[[p_idx]] <-
          lapply(p_combo[[p_idx]], function(rlist) rlist[keep_list,]) %>%
          dplyr::bind_rows() %>%
          tidyr::pivot_longer(cols = dplyr::everything()) %>%
          dplyr::group_by(name) %>%
          dplyr::summarise(
            mean = mean(value),
            median = stats::median(value),
            sd = stats::sd(value),
            ci.05 = stats::quantile(value, 0.05),
            ci.95 = stats::quantile(value, 0.95),
            p0 = mean(value < (0.5/ max(1, (any(p_zero[d_idx, s_idx, w_idx, ])* dplyr::filter(harvest, DISTRICT== d_idx, SUBDISTRICT== s_idx, STAT_WEEK== w_idx)$HARVEST)))),
            .groups = "drop"
          ) %>%
          dplyr::mutate(name = factor(name, levels = c(group_names[seq(ncol(p_combo_all[[max(S)*(d_idx-1)+s_idx]][[1]])), d_idx]))) %>%
          dplyr::arrange(name) %>%
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

    } # S
  } # D

  out <- list()

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

# end pop ----


# age ----

#' Title
#'
#' @param outraw MAGMA output
#' @param dat_in MAGMA data
#' @param nreps same as MAGMA model run
#' @param nburn same as MAGMA model run
#' @param thin same as MAGMA model run
#' @param nchains same as MAGMA model run
#' @param keep_burn same as MAGMA model run
#' @param summ_level Summarize to district or subdistrict
#'
#' @importFrom magrittr %>%
#'
#' @noRd
#'
magmatize_age_b <- function(out2, dat_in, nreps, nburn, thin, nchains, keep_burn, summ_level) {

  ### info needed ### ----
  C <- dat_in$C # number of age classes
  D <- dplyr::n_distinct(dat_in$metadat$district) # number of districts
  S <- dat_in$metadat %>%
    dplyr::group_by(district) %>%
    dplyr::summarise(S = dplyr::n_distinct(subdist), .groups = "drop") %>%
    dplyr::pull(S) # number of subdistricts
  W <- max(dplyr::n_distinct(dat_in$metadat$week), length(dat_in$stat_weeks)) # number of stats weeks

  ng <- apply(dat_in$groups, 2, max) # number of groups

  if (!is.null(dat_in$districts)) dist_names <- dat_in$districts
  if (!is.null(dat_in$subdistricts)) subdist_names <- dat_in$subdistricts
  if (!is.null(dat_in$stat_weeks)) week_names <- dat_in$stat_weeks

  ### prepare output ### ----

  if (summ_level == "district") {
    out <- format_district_age_b(out2, dat_in, nreps, nburn, thin, nchains, keep_burn, C, D, S, W)
  } else if (summ_level == "subdistrict") {
    out <- format_subdistrict_age_b(out2, dat_in, nreps, nburn, thin, nchains, keep_burn, C, D, S, W)
  } else stop("Invalid summ_level.")

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

  return(out)

}


#' Format district age classes part b
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
format_district_age_b <- function(ageout_a, dat_in, nreps, nburn, thin, nchains, keep_burn, C, D, S, W) {

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

  # groups <- dat_in$groups # vector id for reporting groups (aka groupvec)
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

  # exclude burn-in while calculate r hat
  keep_list <- ((nburn*keep_burn + 1):(nreps - nburn*isFALSE(keep_burn)))[!((nburn*keep_burn + 1):(nreps - nburn*isFALSE(keep_burn))) %% thin] / thin

  # place holders/empty objects for age/pop summaries
  ap_prop_grp <- mc_age <- summ_age <- list()
  a_idx <- 0

  for (d_idx in 1:D) {
    ## age (loop through each district) ##
    ap_combo_all = list()
    for(chain in seq(nchains)) {
      ap_combo_all[[chain]] <-
        lapply(ageout_a[[d_idx]],
               function(apc) {
                 apc[[chain]]
               }) %>%
        dplyr::bind_rows() %>%
        dplyr::group_by(itr, grpvec) %>%
        dplyr::summarise(dplyr::across(.cols = dplyr::everything(),
                                .fns = sum),
                         .groups = "drop") %>%
        tidyr::pivot_longer(-c(itr, grpvec)) %>%
        dplyr::group_by(itr, grpvec) %>%
        dplyr::mutate(value = prop.table(value)) %>%
        dplyr::ungroup() %>%
        tidyr::pivot_wider() %>%
        dplyr::select(-itr) %>%
        stats::setNames(., c("grpvec", age_classes))
    } # combine all individual districts

    for (grp in stats::na.omit(group_names[, d_idx])) {
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
        dplyr::bind_rows() %>%
        dplyr::mutate(stock_prop = rowSums(.[,-1])) %>%
        tidyr::pivot_longer(-c(stock_prop, grpvec)) %>%
        dplyr::group_by(name, grpvec) %>%
        dplyr::summarise(
          mean = mean(value),
          median = stats::median(value),
          sd = stats::sd(value),
          ci.05 = stats::quantile(value, 0.05),
          ci.95 = stats::quantile(value, 0.95),
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

  return(out)

}

#' Format subdistrict age classes part b
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
format_subdistrict_age_b <- function(ageout_a, dat_in, nreps, nburn, thin, nchains, keep_burn, C, D, S, W) {

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

  # groups <- dat_in$groups # vector id for reporting groups (aka groupvec)
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

  # p_zero <- array(
  #   apply(
  #     table(metadat[, c("district", "subdist", "week", "iden")], exclude = NULL),
  #     seq(3),
  #     function(alf) any(alf!= 0) # prop = 0 if no sample
  #   ), c(D, max(S), W, K+ H))

  if (isFALSE(keep_burn)|keep_burn == "false") keep_burn = FALSE

  ### prepare output ### ----
  # organization of out_list:
  # [[chain]][[yr]][[dist]][[sub]][[week]][age, pop, itr]

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
          lapply(ageout_a[[d_idx]][[s_idx]],
                 function(apc) {
                   apc[[chain]]
                 }) %>%
          dplyr::bind_rows() %>%
          dplyr::group_by(itr, grpvec) %>%
          dplyr::summarise(dplyr::across(.cols = dplyr::everything(),
                                  .fns = sum),
                           .groups = "drop") %>%
          tidyr::pivot_longer(-c(itr, grpvec)) %>%
          dplyr::group_by(itr, grpvec) %>%
          dplyr::mutate(value = prop.table(value)) %>%
          dplyr::ungroup() %>%
          tidyr::pivot_wider() %>%
          dplyr::select(-itr) %>%
          stats::setNames(., c("grpvec", age_classes))
      } # combine all individual districts

      for (grp in stats::na.omit(group_names[, d_idx])) {
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
          dplyr::bind_rows() %>%
          dplyr::mutate(stock_prop = rowSums(.[,-1])) %>%
          tidyr::pivot_longer(-c(stock_prop, grpvec)) %>%
          dplyr::group_by(name, grpvec) %>%
          dplyr::summarise(
            mean = mean(value),
            median = stats::median(value),
            sd = stats::sd(value),
            ci.05 = stats::quantile(value, 0.05),
            ci.95 = stats::quantile(value, 0.95),
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

  return(out)

}

# end age ----


utils::globalVariables(c("district", "subdist"))














