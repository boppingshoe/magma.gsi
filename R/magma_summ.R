## summarizing district/subdistrict ----

#' Format district
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
#' @param which_dist district to summarize
#' @param fst_files fst files location
#' @param save_trace "in_memory" or file path to save trace history
#'
#' @importFrom magrittr %>%
#'
#' @noRd
#'
format_district <- function(outraw, dat_in, keep_list, nrows_ap_prop, nchains, C, D, S, W, which_dist, fst_files, save_trace) {

  ### data input ### ----
  metadat <- dat_in$metadat # age and iden info
  harvest <- dat_in$harvest %>% # for calculating p0
    dplyr::mutate(n = apply(.[, 2:4], 1, function(aa) {
      dplyr::filter(metadat,
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

  ### prepare output ### ----
  # place holders/empty objects for age/pop summaries
  if (save_trace == "in_memory") {
    out <- list(age_prop = list(), age_summ = list(), pop_prop = list(), pop_summ = list(), pop_prop_all = list(), pop_summ_all = list())
  } else {
    out <- list(age_summ = list(), pop_summ = list(), pop_summ_all = list())
    tr_folder <- paste0(save_trace, "/trace_district")
    if (file.exists(paste0(tr_folder, "/p_d", which_dist[1], "all.fst"))) stop("Old trace files already exist. Please move or delete old files before saving new ones.")
    dir.create(tr_folder)
  }

  p_idx <- a_idx <- 0

  # organization of outraw:
  # [[chain]][age * ((nreps - nburn) / thin) * week * sub * dist, pop + 6]

  for (d_idx in 1:D) {
    if (is.null(outraw)) {
      ap_prop <- # set up for age-pop comp and district pop prop
        lapply(seq.int(nchains), function(ch) {
          fst::read.fst(path = paste0(fst_files, "/magma_raw_ch", ch, ".fst"),
                        from = 1 + nrows_ap_prop * (which_dist[d_idx] - 1) * max(S) * W,
                        to = nrows_ap_prop * which_dist[d_idx] * max(S) * W) %>%
            dplyr::left_join({
              harvest %>%
                dplyr::filter(DISTRICT == which_dist[d_idx]) %>%
                # proportional by district
                dplyr::mutate(prop_harv_d = prop.table(HARVEST)) %>%
                # proportional by week within a district (across subdistricts)
                dplyr::mutate(prop_harv_wk = prop.table(HARVEST), .by = STAT_WEEK) %>%
                dplyr::select(-c(YEAR, HARVEST))
            }, dplyr::join_by(d == DISTRICT, s == SUBDISTRICT, w == STAT_WEEK)) %>%
            tidyr::replace_na(list(prop_harv_d = 0, prop_harv_wk = 0)) %>%
            tidyr::pivot_longer(cols = 1:(ncol(.) - 8),
                                names_to = "collection",
                                values_to = "ppi") %>%
            dplyr::mutate( # faster, but rely on the order of collections == groups
              grpname = rep(group_names[groups[, d_idx], d_idx], nrows_ap_prop * max(S) * W)
            ) %>%
            # dplyr::left_join({ # slower, match names with collection
            #   groups %>%
            #     dplyr::as_tibble(rownames = "collection") %>%
            #     dplyr::mutate(grpname = group_names[groups[, d_idx], d_idx]) %>%
            #     dplyr::select(collection, grpname)
            # }, by = "collection") %>%
            dplyr::mutate(ppi_d = ppi * prop_harv_d,
                          ppi_wk = ppi * prop_harv_wk,
                          .keep = "unused")
        })
    } else {
      ap_prop <-
        lapply(outraw, function(o) {
          o %>%
            dplyr::filter(d == which_dist[d_idx]) %>%
            dplyr::left_join({
              harvest %>%
                dplyr::filter(DISTRICT == which_dist[d_idx]) %>%
                dplyr::mutate(prop_harv_d = prop.table(HARVEST)) %>% # proportional by district
                # proportional by week within a district (across subdistricts)
                dplyr::mutate(prop_harv_wk = prop.table(HARVEST), .by = STAT_WEEK) %>%
                dplyr::select(-c(YEAR, HARVEST))
            }, dplyr::join_by(d == DISTRICT, s == SUBDISTRICT, w == STAT_WEEK)) %>%
            tidyr::replace_na(list(prop_harv_d = 0, prop_harv_wk = 0)) %>%
            tidyr::pivot_longer(cols = 1:(ncol(.) - 8),
                                names_to = "collection",
                                values_to = "ppi") %>%
            dplyr::mutate( # rely on the order of collections == groups
              grpname = rep(group_names[groups[, d_idx], d_idx], nrows_ap_prop * max(S) * W),
              ppi_d = ppi * prop_harv_d,
              ppi_wk = ppi * prop_harv_wk,
              .keep = "unused"
            )
        })
    } # else

    ap_combo_dis <-
      lapply(ap_prop, function(ap) {
        ap %>%
          dplyr::summarise(ppi = sum(ppi_d, na.rm = TRUE),
                           .by = c(itr, agevec, grpname, chain)) %>%
          dplyr::mutate(ppi = prop.table(ppi), .by = c(itr, grpname)) %>%
          tidyr::pivot_wider(names_from = agevec, values_from = ppi)
      })

    p_combo_all <-
      lapply(ap_prop, function(ap) {
        ap %>%
          dplyr::summarise(p = sum(ppi_d, na.rm = TRUE),
                           .by = c(itr, grpname, chain)) %>%
          tidyr::pivot_wider(names_from = grpname, values_from = p)
      })

    if (save_trace == "in_memory") {
      out$pop_prop_all[[d_idx]] <- dplyr::bind_rows(p_combo_all)
    } else {
      tidyfst::export_fst(dplyr::bind_rows(p_combo_all),
                          path = paste0(tr_folder, "/p_d", which_dist[d_idx], "all.fst"))
    }

    out$pop_summ_all[[d_idx]] <-
      lapply(p_combo_all, function(rlist) rlist[keep_list,]) %>%
      dplyr::bind_rows() %>%
      dplyr::select(-c(itr, chain)) %>%
      tidyr::pivot_longer(cols = 1:ncol(.), names_to = "group") %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(
        mean = mean(value),
        median = stats::median(value),
        sd = stats::sd(value),
        ci.05 = stats::quantile(value, 0.05),
        ci.95 = stats::quantile(value, 0.95),
        p0 = mean(value < (0.5/ max(1, ( any(p_zero[d_idx, , , ])* {
          harvest %>%
            dplyr::filter(DISTRICT == which_dist[d_idx]) %>%
            dplyr::summarise(sum_harv = sum(HARVEST)) %>%
            dplyr::pull(sum_harv)
        } )))),
        .groups = "drop"
      ) %>% # group is alphabetically ordered
      dplyr::left_join({
        gr_diag(coda::as.mcmc.list(
          lapply(p_combo_all,
                 function(rlist) coda::mcmc(dplyr::select(rlist, -c(itr, chain))[keep_list,])) ), nchains)
      }, by = "group") %>%
      dplyr::mutate(
        group = factor(group,
                       levels = c(group_names[!is.na(group_names[, d_idx]), d_idx]))
      ) %>%
      dplyr::arrange(group)

    for (w_idx in 1:W) {
      p_idx <- p_idx + 1 # index individual sampling period/week

      p_combo <-
        lapply(ap_prop, function(ap) {
          ap %>%
            dplyr::filter(w == w_idx) %>%
            dplyr::summarise(p = sum(ppi_wk, na.rm = TRUE),
                             .by = c(itr, grpname, chain)) %>%
            tidyr::pivot_wider(names_from = grpname, values_from = p)
        })

      if (save_trace == "in_memory") {
        out$pop_prop[[p_idx]] <- dplyr::bind_rows(p_combo)
      } else {
        tidyfst::export_fst(dplyr::bind_rows(p_combo),
                            path = paste0(tr_folder, "/p_d", which_dist[d_idx], "w", w_idx, ".fst"))
      }

      out$pop_summ[[p_idx]] <-
        lapply(p_combo, function(rlist) rlist[keep_list,]) %>%
        dplyr::bind_rows() %>%
        dplyr::select(-c(itr, chain)) %>%
        tidyr::pivot_longer(cols = 1:ncol(.), names_to = "group") %>%
        dplyr::group_by(group) %>%
        dplyr::summarise(
          mean = mean(value),
          median = stats::median(value),
          sd = stats::sd(value),
          ci.05 = stats::quantile(value, 0.05),
          ci.95 = stats::quantile(value, 0.95),
          p0 = mean(value < (0.5/ max(1, ( any(p_zero[d_idx, , w_idx, ])* {
            harvest %>%
              dplyr::filter(DISTRICT == which_dist[d_idx], STAT_WEEK== w_idx) %>%
              dplyr::summarise(sum_harv = sum(HARVEST), .groups = "drop") %>%
              dplyr::pull(sum_harv)
          } )))),
          .groups = "drop"
        ) %>%
        dplyr::left_join({
          gr_diag(coda::as.mcmc.list(
            lapply(p_combo,
                   function(rlist) coda::mcmc(dplyr::select(rlist, -c(itr, chain))[keep_list,])) ), nchains)
        }, by = "group") %>%
        dplyr::mutate(
          group = factor(group,
                         levels = c(group_names[!is.na(group_names[, d_idx]), d_idx]))
        ) %>%
        dplyr::arrange(group)

    } # w_idx

    for (grp in stats::na.omit(group_names[, d_idx])) {
      a_idx <- a_idx + 1

      ap_prop_grp <-
        lapply(ap_combo_dis, function(ap) {
          dplyr::filter(ap, grpname == grp)
        }) # separate rep groups into its own output

      if (save_trace == "in_memory") {
        out$age_prop[[a_idx]] <-
          dplyr::bind_rows(ap_prop_grp) %>%
          dplyr::select(-grpname) # trace plot function doesn't take group names
      } else {
        tidyfst::export_fst(dplyr::bind_rows(ap_prop_grp) %>% dplyr::select(-grpname),
                            path = paste0(tr_folder, "/ap_d", which_dist[d_idx], grp, ".fst"))
      }

      harv_dis <- harvest %>%
        dplyr::filter(DISTRICT == which_dist[d_idx]) %>%
        dplyr::pull(HARVEST) %>%
        mean()

      out$age_summ[[a_idx]] <-
        lapply(ap_prop_grp, function(rlist) rlist[keep_list,]) %>%
        dplyr::bind_rows() %>%
        dplyr::select(-c(itr, chain)) %>%
        dplyr::mutate(stock_prop = rowSums(.[, -1])) %>%
        tidyr::pivot_longer(-c(stock_prop, grpname), names_to = "age") %>%
        dplyr::group_by(age, grpname) %>%
        dplyr::summarise(
          mean = mean(value),
          median = stats::median(value),
          sd = stats::sd(value),
          ci.05 = stats::quantile(value, 0.05),
          ci.95 = stats::quantile(value, 0.95),
          p0 = mean( value < (0.5/ max(1, harv_dis* stock_prop)) ),
          .groups = "drop"
        ) %>%
        dplyr::left_join({
          gr_diag(
            coda::as.mcmc.list(
              lapply(ap_prop_grp, function(rlist) {
              coda::mcmc(rlist[keep_list,] %>% dplyr::select(-c(itr, grpname, chain)))
              }) )
            , nchains)
        }, dplyr::join_by(age == group)) %>%
        dplyr::mutate(
          age = factor(age, levels = age_classes)
        ) %>%
        dplyr::rename(group = grpname) %>%
        dplyr::relocate(group, .before = age) %>%
        dplyr::arrange(age)

    } # grp
  } # d_idx

  return(out)

}


#' Format subdistrict
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
#' @param which_dist district to summarize
#' @param fst_files fst files location
#' @param save_trace "in_memory" or file path to save trace history
#'
#' @importFrom magrittr %>%
#'
#' @noRd
#'
format_subdistrict <- function(outraw, dat_in, keep_list, nrows_ap_prop, nchains, C, D, S, W, which_dist, fst_files, save_trace) {

  ### data input ### ----
  metadat <- dat_in$metadat # age and iden info
  harvest <- dat_in$harvest %>% # for calculating p0
    dplyr::mutate(n = apply(.[, 2:4], 1, function(aa) {
      dplyr::filter(metadat,
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

  ### prepare output ### ----
  # place holders/empty objects for age/pop summaries
  if (save_trace == "in_memory") {
    out <- list(age_prop = list(), age_summ = list(), pop_prop = list(), pop_summ = list(), pop_prop_all = list(), pop_summ_all = list())
  } else {
    out <- list(age_summ = list(), pop_summ = list(), pop_summ_all = list())
    tr_folder <- paste0(save_trace, "/trace_subdistrict")
    if (file.exists(paste0(tr_folder, "/p_d", which_dist[1], "s1all.fst"))) stop("Old trace files already exist. Please move or delete old files before saving new ones.")
    dir.create(tr_folder)
  }

  p_idx <- a_idx <- 0

  # organization of outraw:
  # [[chain]][age * ((nreps - nburn) / thin) * week * sub * dist, pop + 6]

  for (d_idx in 1:D) {
    for (s_idx in 1:max(S)) {
      if (is.null(outraw)) {
        ap_prop <- # set up for age-pop comp and subdistrict pop props
          lapply(seq.int(nchains), function(ch) {
            fst::read.fst(path = paste0(fst_files, "/magma_raw_ch", ch, ".fst"),
                          from = 1 + nrows_ap_prop * (which_dist[d_idx] - 1) * max(S) * W +
                            nrows_ap_prop * (s_idx - 1) * W,
                          to = nrows_ap_prop * (which_dist[d_idx] - 1) * max(S) * W +
                            nrows_ap_prop * s_idx * W) %>%
              dplyr::left_join({
                harvest %>%
                  dplyr::filter(DISTRICT == which_dist[d_idx],
                                SUBDISTRICT == s_idx) %>%
                  # proportional within a subdistrict
                  dplyr::mutate(prop_harv_s = prop.table(HARVEST)) %>%
                  dplyr::select(-c(YEAR, HARVEST))
              }, dplyr::join_by(d == DISTRICT, s == SUBDISTRICT, w == STAT_WEEK)) %>%
              tidyr::replace_na(list(prop_harv_s = 0)) %>%
              tidyr::pivot_longer(cols = 1:(ncol(.) - 7),
                                  names_to = "collection",
                                  values_to = "ppi") %>%
              dplyr::mutate( # rely on the order of collections == groups
                grpname = rep(group_names[groups[, d_idx], d_idx], nrows_ap_prop * W),
                ppi_s = ppi * prop_harv_s
              )
          })
      } else {
        ap_prop <-
          lapply(outraw, function(o) {
            o %>%
              dplyr::filter(d == which_dist[d_idx],
                            s == s_idx) %>%
              dplyr::left_join({
                harvest %>%
                  dplyr::filter(DISTRICT == which_dist[d_idx],
                                SUBDISTRICT == s_idx) %>%
                  # proportional within a subdistrict
                  dplyr::mutate(prop_harv_s = prop.table(HARVEST)) %>%
                  dplyr::select(-c(YEAR, HARVEST))
              }, dplyr::join_by(d == DISTRICT, s == SUBDISTRICT, w == STAT_WEEK)) %>%
              tidyr::replace_na(list(prop_harv_s = 0)) %>%
              tidyr::pivot_longer(cols = 1:(ncol(.) - 7),
                                  names_to = "collection",
                                  values_to = "ppi") %>%
              dplyr::mutate( # rely on the order of collections == groups
                grpname = rep(group_names[groups[, d_idx], d_idx], nrows_ap_prop * W),
                ppi_s = ppi * prop_harv_s
              )
          })
      } # else

      ap_combo_subdis <-
        lapply(ap_prop, function(ap) {
          ap %>%
            dplyr::summarise(ppi = sum(ppi_s, na.rm = TRUE),
                             .by = c(itr, agevec, grpname, chain)) %>%
            dplyr::mutate(ppi = prop.table(ppi), .by = c(itr, grpname)) %>%
            tidyr::pivot_wider(names_from = agevec, values_from = ppi)
        })

      p_combo_all <-
        lapply(ap_prop, function(ap) {
          ap %>%
            dplyr::summarise(p = sum(ppi_s, na.rm = TRUE),
                             .by = c(itr, grpname, chain)) %>%
            tidyr::pivot_wider(names_from = grpname, values_from = p)
        })

      if (save_trace == "in_memory") {
        out$pop_prop_all[[max(S)*(d_idx-1)+s_idx]] <- dplyr::bind_rows(p_combo_all)
      } else {
        tidyfst::export_fst(dplyr::bind_rows(p_combo_all),
                            path = paste0(tr_folder, "/p_d", which_dist[d_idx], "s", s_idx, "all.fst"))
      }

      out$pop_summ_all[[max(S)*(d_idx-1)+s_idx]] <-
        lapply(p_combo_all, function(rlist) rlist[keep_list,]) %>%
        dplyr::bind_rows() %>%
        dplyr::select(-c(itr, chain)) %>%
        tidyr::pivot_longer(cols = 1:ncol(.), names_to = "group") %>%
        dplyr::group_by(group) %>%
        dplyr::summarise(
          mean = mean(value),
          median = stats::median(value),
          sd = stats::sd(value),
          ci.05 = stats::quantile(value, 0.05),
          ci.95 = stats::quantile(value, 0.95),
          p0 = mean(value < (0.5/ max(1, ( any(p_zero[d_idx, s_idx, , ])* {
            harvest %>%
              dplyr::filter(DISTRICT == which_dist[d_idx],
                            SUBDISTRICT == s_idx) %>%
              dplyr::summarise(sum_harv = sum(HARVEST), .groups = "drop") %>%
              dplyr::pull(sum_harv)
          } )))),
          .groups = "drop"
        ) %>% # group is alphabetically ordered
        dplyr::left_join({
          gr_diag(coda::as.mcmc.list(
            lapply(p_combo_all,
                   function(rlist) coda::mcmc(dplyr::select(rlist, -c(itr, chain))[keep_list,])) ), nchains)
        }, by = "group") %>%
        dplyr::mutate(
          group = factor(group,
                         levels = c(group_names[!is.na(group_names[, d_idx]), d_idx]))
        ) %>%
        dplyr::arrange(group)

      for (w_idx in 1:W) {
        p_idx <- p_idx + 1 # index individual sampling period/week

        p_combo <-
          lapply(ap_prop, function(ap) {
            ap %>%
              dplyr::filter(w == w_idx) %>%
              dplyr::summarise(p = sum(ppi, na.rm = TRUE),
                               .by = c(itr, grpname, chain)) %>%
              tidyr::pivot_wider(names_from = grpname, values_from = p)
          })

        if (save_trace == "in_memory") {
          out$pop_prop[[p_idx]] <- dplyr::bind_rows(p_combo)
        } else {
          tidyfst::export_fst(dplyr::bind_rows(p_combo),
                              path = paste0(tr_folder, "/p_d", which_dist[d_idx], "s", s_idx, "w", w_idx, ".fst"))
        }

        out$pop_summ[[p_idx]] <-
          lapply(p_combo, function(rlist) rlist[keep_list,]) %>%
          dplyr::bind_rows() %>%
          dplyr::select(-c(itr, chain)) %>%
          tidyr::pivot_longer(cols = 1:ncol(.), names_to = "group") %>%
          dplyr::group_by(group) %>%
          dplyr::summarise(
            mean = mean(value),
            median = stats::median(value),
            sd = stats::sd(value),
            ci.05 = stats::quantile(value, 0.05),
            ci.95 = stats::quantile(value, 0.95),
            p0 = mean(value < (0.5/ max(1, ( any(p_zero[d_idx, s_idx, w_idx, ])* {
              harvest %>%
                dplyr::filter(DISTRICT == which_dist[d_idx],
                              SUBDISTRICT == s_idx,
                              STAT_WEEK== w_idx) %>%
                dplyr::summarise(sum_harv = sum(HARVEST), .groups = "drop") %>%
                dplyr::pull(sum_harv)
            } )))),
            .groups = "drop"
          ) %>%
          dplyr::left_join({
            gr_diag(coda::as.mcmc.list(
              lapply(p_combo,
                     function(rlist) coda::mcmc(dplyr::select(rlist, -c(itr, chain))[keep_list,])) ), nchains)
          }, by = "group") %>%
          dplyr::mutate(
            group = factor(group,
                           levels = c(group_names[!is.na(group_names[, d_idx]), d_idx]))
          ) %>%
          dplyr::arrange(group)

      } # w_idx

      # loop through each subdistrict
      for (grp in stats::na.omit(group_names[, d_idx])) {
        a_idx <- a_idx + 1

        ap_prop_grp <-
          lapply(ap_combo_subdis, function(ap) {
              dplyr::filter(ap, grpname == grp)
          }) # separate rep groups into its own output

        if (save_trace == "in_memory") {
          out$age_prop[[a_idx]] <-
            dplyr::bind_rows(ap_prop_grp) %>%
            dplyr::select(-grpname) # trace plot function doesn't take group names
        } else {
          tidyfst::export_fst(dplyr::bind_rows(ap_prop_grp) %>% dplyr::select(-grpname),
                              path = paste0(tr_folder, "/ap_d", which_dist[d_idx], "s", s_idx, grp, ".fst"))
        }

        harv_subdis <- harvest %>%
          dplyr::filter(DISTRICT == which_dist[d_idx], SUBDISTRICT == s_idx) %>%
          dplyr::pull(HARVEST) %>%
          mean()

        out$age_summ[[a_idx]] <-
          lapply(ap_prop_grp, function(rlist) rlist[keep_list,]) %>%
          dplyr::bind_rows() %>%
          dplyr::select(-c(itr, chain)) %>%
          dplyr::mutate(stock_prop = rowSums(.[, -1])) %>%
          tidyr::pivot_longer(-c(stock_prop, grpname), names_to = "age") %>%
          dplyr::group_by(age, grpname) %>%
          dplyr::summarise(
            mean = mean(value),
            median = stats::median(value),
            sd = stats::sd(value),
            ci.05 = stats::quantile(value, 0.05),
            ci.95 = stats::quantile(value, 0.95),
            p0 = mean( value < (0.5/ max(1, harv_subdis* stock_prop)) ),
            .groups = "drop"
          ) %>%
          dplyr::left_join({
            gr_diag(
              coda::as.mcmc.list(
                lapply(ap_prop_grp, function(rlist) {
                  coda::mcmc(rlist[keep_list,] %>% dplyr::select(-c(itr, grpname, chain)))
                }) )
              , nchains)
          }, dplyr::join_by(age == group)) %>%
          dplyr::mutate(
            age = factor(age, levels = age_classes)
          ) %>%
          dplyr::rename(group = grpname) %>%
          dplyr::relocate(group, .before = age) %>%
          dplyr::arrange(age)

      } # grp
    } # s_idx
  } # d_idx

  return(out)

}


## Wrapper for model output summarizers (only one that is exported) ----

#' Summarize model output
#'
#' @param ma_out MAGMA output
#' @param ma_dat MAGMA input data
#' @param summ_level Summarize at district or subdistrict level
#' @param which_dist Function format raw magma output one district at a time.
#'   Identify district as 1, 2, ... Default = NULL will summarize all districts.
#' @param fst_files Fst files location if MAGMA model output was saved.
#' @param save_trace default = "in_memory" to have trace history as a part of summary. Or specify the path of a directory to save trace history as fst files.
#'
#' @return Summary tables for reporting groups and/or age classes.
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' # format data
#' wd <- getwd() # path to data folder
#' magma_data <- magmatize_data(wd = wd, save_data = FALSE)
#'
#' # model run
#' magma_out <- magmatize_mdl(magma_data, nreps = 50, nburn = 25, thin = 1, nchains = 2)
#'
#' # summary using output as an object
#' magma_summ <- magmatize_summ(ma_out = magma_out, ma_dat = magma_data, summ_level = "district", which_dist = 1, save_trace = wd)
#'
#' # summary using output saved as Fst files
#' magma_summ <- magmatize_summ(ma_dat = magma_data, summ_level = "district", which_dist = 1, fst_files = wd)
#' }
#'
#' @export
magmatize_summ <- function(ma_out = NULL, ma_dat, summ_level, which_dist = NULL, fst_files = NULL, save_trace = "in_memory") {

  if (is.null(ma_out) & is.null(fst_files)) stop("There's no MAGMA output file. You need to provide the MAGMA output as an object or provide the file path for the saved Fst files.")
  if (!is.null(ma_out) & !is.null(fst_files)) message("You provided MAGMA output both as an object and Fst file path. Summary is done using the saved Fst files, just so you know.")
  if (save_trace != "in_memory" & !dir.exists(save_trace)) stop("wrong specification for `save_trace`. It has to be `in_memory` or a directory path to a folder to save trace output as Fst files.")

  if (is.null(which_dist)) which_dist <- unique(ma_dat$metadat$district)

  sub_dat <- list(
    x = ma_dat$x[which(ma_dat$metadat$district %in% which_dist), ],
    y = ma_dat$y,

    metadat = ma_dat$metadat %>%
      dplyr::filter(district %in% which_dist),
    harvest = ma_dat$harvest %>%
      dplyr::filter(DISTRICT %in% which_dist),

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

  ### info needed ----
  C <- sub_dat$C # number of age classes
  D <- dplyr::n_distinct(sub_dat$metadat$district) # number of districts
  S <- sub_dat$metadat %>%
    dplyr::group_by(district) %>%
    dplyr::summarise(S = dplyr::n_distinct(subdist), .groups = "drop") %>%
    dplyr::pull(S) # number of subdistricts
  W <- max(dplyr::n_distinct(sub_dat$metadat$week), length(sub_dat$stat_weeks)) # number of stats weeks
  KH <- nrow(sub_dat$y)
  ng <- apply(sub_dat$groups, 2, max) # number of groups

  dist_names <- sub_dat$districts
  subdist_names <- sub_dat$subdistricts
  week_names <- sub_dat$stat_weeks

  age_classes <- sub_dat$age_classes

  # organization of outraw:
  # [[chain]][dist*sub*week*age_class*itr, pop] and ordered by district
  if (!is.null(fst_files)) {
    ma_out <- list()
    ma_out$specs <- readRDS(paste0(fst_files, "/magma_specs.rds"))
  }

  nreps <- ma_out$specs["nreps"]
  nburn <- ma_out$specs["nburn"]
  thin <- ma_out$specs["thin"]
  nchains <- ma_out$specs["nchains"]
  keep_burn <- ma_out$specs["keep_burn"] == 1

  # if (isFALSE(keep_burn)|keep_burn == "false") keep_burn <- FALSE

  # exclude burn-in while calculate r hat in formatting functions
  keep_list <- ((nburn*keep_burn + 1):(nreps - nburn*isFALSE(keep_burn)))[!((nburn*keep_burn + 1):(nreps - nburn*isFALSE(keep_burn))) %% thin] / thin

  nrows_ap_prop <- C * (nreps - nburn*isFALSE(keep_burn)) / thin

  ### detect output type and convert ----
  if (!is.null(fst_files)) { # output is saved as Fst files
    holder <- NULL
  } else if (is.list(ma_out$outraw)) { # output as object
    if (all(lapply(ma_out$outraw, is.data.frame) %>% unlist())) {
      # output in 2024 new format
      holder <- lapply(ma_out$outraw, function(o) dplyr::filter(o, d %in% which_dist))
    } else {
      # output is in old (2023) format
      holder <- old_to_new(ma_out$outraw, nrows_ap_prop, nchains, C, W, S, D, age_classes, which_dist)
    }
  } else if (is.array(ma_out$outraw)) { # malia
    if (length(dim(ma_out$outraw)) == 3) { # 2024
      holder <- lapply(seq.int(nchains), function(ch) {
        ma_out$outraw[ , , ch] %>%
          data.table::as.data.table() %>%
          dplyr::rename_with(~c(row.names(sub_dat$y), "itr", "agevec", "d", "s", "w", "chain")) %>%
          dplyr::mutate(agevec = age_classes[agevec])
      })
    } else if (length(dim(ma_out$outraw == 6))) { # 2023
      holder <-
        lapply(1:nchains, function(ch) {
          lapply(1:D, function(d_i) {
            lapply(1:S[d_i], function(s_i) {
              lapply(1:W, function (w_i) {
                data.table::as.data.table(out_jl[ , , w_i, s_i, d_i, ch]) %>%
                  dplyr::rename_with(~c(row.names(sub_dat$y), "itr", "agevec")) %>%
                  dplyr::mutate(d = d_i, s = s_i, w = w_i, chain = ch,
                                agevec = age_classes[agevec])
              }) %>% dplyr::bind_rows()
            }) %>% dplyr::bind_rows()
          }) %>% dplyr::bind_rows()
        })
    }
  } else stop("Format of MAGMA output is not recognized and cannot be summarized using current version of magmatize_summ().")

  ### prepare output ----
  message("Preparing output (patience grasshopper...)")
  prep_time <- Sys.time()

  if (summ_level == "district") {
    out <- format_district(holder, sub_dat, keep_list, nrows_ap_prop, nchains, C, D, S, W, which_dist, fst_files, save_trace)
  } else if (summ_level == "subdistrict") {
    out <- format_subdistrict(holder, sub_dat, keep_list, nrows_ap_prop, nchains, C, D, S, W, which_dist, fst_files, save_trace)
  } else stop("Invalid summ_level.")

  # id for age-by-group output
  if (summ_level == "district") {
    aout_names <- rep(NA, sum(ng))
    i <- 0
    for (d_i in 1:D) {
      for (g_i in 1:ng[d_i]) {
        i <- i + 1
        aout_names[i] <-
          paste(
            paste0("D", dist_names[d_i]),
            paste0(sub_dat$group_names[g_i, d_i]),
            sep = "_")
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
            paste(
              paste0("D", dist_names[d_i]),
              paste0("S", subdist_names[[d_i]][s_i]),
              paste0(sub_dat$group_names[g_i, d_i]),
              sep = "_")
        } # g_i
      } # s_i
    } # d_i
  } # end if

  if (save_trace == "in_memory") {
    names(out[[1]]) <- names(out[[2]]) <- aout_names
  } else names(out[[1]]) <- aout_names

  # id for pop output (weekly)
  if (summ_level == "district") {
    pout_names <- rep(NA, D*W)
    j <- 1
    for (d_i in 1:D) {
      for (w_i in 1:W) {
        pout_names[j] <-
          paste(
            paste0("D", dist_names[d_i]),
            paste0("W", week_names[w_i]),
            sep = "_")
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
            paste(
              paste0("D", dist_names[d_i]),
              paste0("S", subdist_names[[d_i]][s_i]),
              paste0("W", week_names[w_i]),
              sep = "_")
          j <- j + 1
        } # w_i
      } # s_i
    } # d_i
  } # end if dist/subdist

  if (save_trace == "in_memory") {
    names(out[[3]]) <- names(out[[4]]) <- pout_names
  } else names(out[[2]]) <- pout_names

  # id for pop output (combined)
  if (summ_level == "district") {
    if (save_trace == "in_memory") {
      names(out[[5]]) <- names(out[[6]]) <- paste0("D", dist_names)
    } else names(out[[3]]) <- paste0("D", dist_names)
  } else if (summ_level == "subdistrict") {
    if (save_trace == "in_memory") {
      names(out[[5]]) <- names(out[[6]]) <-
        subdist_names %>%
        dplyr::bind_cols() %>%
        tidyr::pivot_longer(dplyr::everything()) %>%
        dplyr::arrange(name) %>%
        dplyr::mutate(name = paste0("D", name)) %>%
        tidyr::unite("subdist", dplyr::everything(), sep = "_S") %>%
        dplyr::pull(subdist)
    } else {
      names(out[[3]]) <- subdist_names %>%
        dplyr::bind_cols() %>%
        tidyr::pivot_longer(dplyr::everything()) %>%
        dplyr::arrange(name) %>%
        dplyr::mutate(name = paste0("D", name)) %>%
        tidyr::unite("subdist", dplyr::everything(), sep = "_S") %>%
        dplyr::pull(subdist)
    }
  }

  print(Sys.time() - prep_time)
  message(Sys.time())

  return(out)

}


## utility functions ----

#' Summarize Gelman-Rubin diagnostics
#'
#' @param mc_list Coda MCMC object
#' @param nchians Number of MC chains
#'
#' @noRd
#'
gr_diag <- function(mc_list, nchains) {
  GR <- {if (nchains > 1) {
    coda::gelman.diag(mc_list,
                      transform = TRUE,
                      autoburnin = FALSE,
                      multivariate = FALSE)$psrf[,"Point est."]
  } else {NA}}
  n_eff <- coda::effectiveSize(mc_list)

  dplyr::as_tibble(cbind(GR, n_eff), rownames = "group")
}


#' 2022 output format compatible
#'
#' @param oldraw Old format output
#' @param nchains Number of MC chains
#' @param nreps Number of MC reps
#' @param C Number of age classes
#' @param W Number of weeks
#' @param S Number of subdistricts
#' @param D Number of districts
#' @param which_dist Subset by district
#'
#' @noRd
#'
old_to_new <- function(oldraw, nrows_ap_prop, nchains, C, W, S, D, age_classes, which_dist) {
  lapply(1:nchains, function(ch) {
    lapply(oldraw[[ch]], function(dis) {
      lapply(dis, function(subdis) {
        lapply(subdis, function(wk) {
          wk
        }) |> dplyr::bind_rows()
      }) |> dplyr::bind_rows()
    }) |> dplyr::bind_rows() |>
      dplyr::mutate(agevec = age_classes[agevec],
                    w = rep(rep(rep(1:W, each = nrows_ap_prop), S), D),
                    s = rep(rep(1:S, each = W*nrows_ap_prop), D),
                    d = rep(1:D, each = S*W*nrows_ap_prop),
                    chain = ch) |>
      dplyr::filter(d %in% which_dist)
  })
}


utils::globalVariables(c(".", "district", "subdist", "week", "HARVEST", "n", "agevec",
                         "grpvec", "ppi", "YEAR", "DISTRICT", "SUBDISTRICT", "STAT_WEEK",
                         "prop_harv", "p", "sum_harv", "stock_prop", "group", "age",
                         "d", "s", "w", "prop_harv_d", "prop_harv_wk", "ppi_d", "ppi_wk",
                         "grpname", "prop_harv_s", "ppi_s", "ch", "out_jl"))










