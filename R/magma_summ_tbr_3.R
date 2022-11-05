
#' TBR output format wrapper step 3
#'
#' @param out2 Output from step 2
#'
#' @return Model output in process and metadata as a list object.
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
#' # format data
#' wd <- "D:/bobby_adfg/backup_013122/projects/magma/test_TBR" # path to data folder
#' magma_data <- magmatize_data(wd = wd, save_data = FALSE)
#'
#' # model run
#' magma_out <- magmatize_mdl(magma_data,
#'   nreps = 50, nburn = 25, thin = 1, nchains = 3)
#'
#' # summary steps 1 and 2
#' tbr1 <- magmatize_summ_tbr1(which_dist = 2,
#'   outraw = magma_out,
#'   ma_dat = magma_data,
#'   nreps = 50, nburn = 25, thin = 1, nchains = 3,
#'   summ_level = "district", type = "pop")
#'
#' tbr2 <- magmatize_summ_tbr2(tbr1)
#'
#' tbr3 <- magmatize_summ_tbr3(tbr2)
#'
magmatize_summ_tbr3 <- function(out2) {

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
    sub_out_b <- magmatize_age_b(out2$holder, out1$sub_dat, nreps, nburn, thin, nchains, keep_burn, summ_level)
  }

  return(sub_out_b)

}

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

  message("Preparing output (patience grasshopper...)")
  prep_time <- Sys.time()

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
      names(out[[li]]) <- subdist_names %>% unlist() %>% names()
    }
  }

  print(Sys.time() - prep_time)
  message(Sys.time())

  return(out)

}

# pop ----
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
      tidyr::pivot_longer(cols = everything()) %>%
      dplyr::group_by(name) %>%
      dplyr::summarise(
        mean = mean(value),
        median = median(value),
        sd = sd(value),
        ci.05 = quantile(value, 0.05),
        ci.95 = quantile(value, 0.95),
        p0 = mean(value < (0.5/ max(1, ( any(p_zero[d_idx, , , ])* harvest %>% dplyr::group_by(DISTRICT) %>% dplyr::summarise(sum_harv = sum(HARVEST), .groups = "drop") %>% dplyr::filter(DISTRICT == d_idx) %>% dplyr::pull(sum_harv) )))),
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
        tidyr::pivot_longer(cols = everything()) %>%
        dplyr::group_by(name) %>%
        dplyr::summarise(
          mean = mean(value),
          median = median(value),
          sd = sd(value),
          ci.05 = quantile(value, 0.05),
          ci.95 = quantile(value, 0.95),
          p0 = mean(value < (0.5/ max(1, ( any(p_zero[d_idx, , w_idx, ])* harvest %>% dplyr::group_by(DISTRICT, STAT_WEEK) %>% dplyr::summarise(sum_harv = sum(HARVEST), .groups = "drop") %>% dplyr::filter(DISTRICT == d_idx, STAT_WEEK== w_idx) %>% dplyr::pull(sum_harv) )))),
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
        tidyr::pivot_longer(cols = everything()) %>%
        dplyr::group_by(name) %>%
        dplyr::summarise(
          mean = mean(value),
          median = median(value),
          sd = sd(value),
          ci.05 = quantile(value, 0.05),
          ci.95 = quantile(value, 0.95),
          p0 = mean(value < (0.5/ max(1, ( any(p_zero[d_idx, , , ])* harvest %>% dplyr::group_by(DISTRICT, SUBDISTRICT) %>% dplyr::summarise(sum_harv = sum(HARVEST), .groups = "drop") %>% dplyr::filter(DISTRICT == d_idx, SUBDISTRICT == s_idx) %>% dplyr::pull(sum_harv) )))),
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
          tidyr::pivot_longer(cols = everything()) %>%
          dplyr::group_by(name) %>%
          dplyr::summarise(
            mean = mean(value),
            median = median(value),
            sd = sd(value),
            ci.05 = quantile(value, 0.05),
            ci.95 = quantile(value, 0.95),
            p0 = mean(value < (0.5/ max(1, (any(p_zero[d_idx, s_idx, w_idx, ])* dplyr::filter(harvest, DISTRICT== d_idx, SUBDISTRICT== s_idx, STAT_WEEK== w_idx)$HARVEST)))),
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


















