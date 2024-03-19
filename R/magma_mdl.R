
#' Run MAGMA model
#'
#' @param dat_in Input data list.
#' @param nreps Amount of simulations to run.
#' @param nburn Number of burn-in simulations to discard.
#' @param thin At what interval to keep the simulations.
#' @param nchains Number of independent MCMC chains run in the simulation.
#' @param nadapt Amount of warm-up/adapt runs before the simulation (only for fully Bayesian mode).
#' @param keep_burn Logical (default = `FALSE`). To keep the burn-ins in the output or not.
#' @param age_prior Option to adjust prior weight on the age class proportions.
#'  * `flat` (default): conventional setup that puts a flat prior on the age proportions.
#'  * `zero_out`: prior weight will concentrate on the major age groups that are observed in metadata (i.e., "zero out" the unobserved age classes).
#'  * `weak_flat`: a less influential flat prior than the conventional one
#' @param cond_gsi Logical (default = `TRUE`). Option to use conditional GSI model. See vignette for details.
#'  * `TRUE`: run MAGMA with a hybrid algorithm of conditional GSI and fully Bayesian.
#'  * `FALSE`: run MAGMA with fully Bayesian algorithm.
#' @param file File path for saving the output.
#'   The default is `NULL` for not saving the output.
#' @param seed Option to initialize a pseudo-random number generator (set random seed)
#'   so the output can be reproduced exactly.
#'   Just pick a seed number and make note of it for future reference.
#'   The default is `NULL`.
#' @param iden_output Option to have trace history for individual group membership assignments included in the final output. Default is FALSE.
#'
#' @return A list object contains:
#'   - the raw output of MAGMA as a list/multi-way array that need to be further summarized using summary functions,
#'   - specifications for the model run (information needed for summary),
#'   - and individual group membership assignment history (optional).
#'
#' @importFrom magrittr %>%
#' @importFrom doRNG %dorng%
#' @importFrom foreach %dopar%
#'
#' @examples
#' \dontrun{
#' # format data
#' wd <- getwd() # path to data folder
#' magma_data <- magmatize_data(wd = wd, save_data = FALSE)
#'
#' # model run
#' magma_out <- magmatize_mdl(magma_data, nreps = 50, nburn = 25, thin = 1, nchains = 2)
#' }
#'
#' @export
magmatize_mdl <- function(dat_in, nreps, nburn, thin, nchains, nadapt = 50, keep_burn = FALSE, age_prior = "flat", cond_gsi = TRUE, file = NULL, seed = NULL, iden_output = FALSE) {

  ### save test file (specs) ----
  specs <- c(nreps, nburn, thin, nchains, keep_burn) %>%
    stats::setNames(c("nreps", "nburn", "thin", "nchains", "keep_burn"))

  if(!is.null(file)) {
    message(paste0("MAGMA specifications saved in ", file, "/magma_specs.rds"))
    saveRDS(specs, file = paste0(file, "/magma_specs.rds"))
  }

  ### data input ----
  x <- as.matrix(dat_in$x) # mixture
  y <- as.matrix(dat_in$y) # base
  metadat <- dat_in$metadat # age and iden info
  nstates <- dat_in$nstates # number of allele types and age classes
  nalleles <- dat_in$nalleles # number of allele types
  C <- dat_in$C # number of age classes
  D <- dplyr::n_distinct(metadat$district) # number of districts
  S <- metadat %>%
    dplyr::group_by(district) %>%
    dplyr::summarise(S = dplyr::n_distinct(subdist), .groups = "drop") %>%
    dplyr::pull(S) # number of subdistricts
  W <- dplyr::n_distinct(metadat$week) # number of stats weeks
  groups <- dat_in$groups # vector id for reporting groups (aka groupvec)

  age_class <- dat_in$age_class # vector id for age classes

  wildpops <- dat_in$wildpops
  K <- length(wildpops)

  if (is.null(dat_in$hatcheries)) {
    if (nrow(groups) > length(wildpops)) {
      stop("Reporting groups had hatcheries but no hatchery was identified in the input data.")
    }
    hatcheries <- NULL
    H <- 0
  } else {
    hatcheries <- dat_in$hatcheries
    H <- length(hatcheries)
  }

  allpops <- c(wildpops, hatcheries)

  if (any(metadat$iden > (K + H) | metadat$iden > nrow(groups), na.rm = TRUE)) stop("Unidentified populations in metadat$iden. Maybe there are hatcheries in the data that are not listed in the reporting groups?")

  na_a <- which(is.na(metadat$age))
  na_i <- which(is.na(metadat$iden))

  metadat$age <- factor(metadat$age, levels = seq(C))
  metadat$iden <- factor(metadat$iden, levels = seq(K + H))

  trait_fac <- factor(rep(names(nstates), nstates), levels = names(nstates))
  states <- paste0(paste0(trait_fac, "_"), unlist(lapply(nstates, seq.int)))

  alleles <- states[seq.int(sum(nalleles))] # allele types
  ages <- states[-seq.int(sum(nalleles))] # age classes

  ### specifications ----
  rdirich <- function(alpha0) {
    if (sum(alpha0) > 0) {
      vec = stats::rgamma(length(alpha0), alpha0, 1)
      vec = vec / sum(vec)
      vec[vec == 0] = .Machine$double.xmin
      vec
    } else{
      rep(0, length(alpha0))
    }
  } # og random dirichlet by jj

  if (cond_gsi) nadapt = 0
  if (keep_burn) nburn = 0

  ### ballroom categories ----
  categories <- c("Live, Werk, Pose", "Bring It Like Royalty", "Face", "Best Mother", "Best Dressed", "High Class In A Fur Coat", "Snow Ball", "Butch Queen Body", "Weather Girl", "Labels", "Mother-Daughter Realness", "Working Girl", "Linen Vs. Silk", "Perfect Tens", "Modele Effet", "Stone Cold Face", "Realness", "Intergalatic Best Dressed", "House Vs. House", "Femme Queen Vogue", "High Fashion In Feathers", "Femme Queen Runway", "Lofting", "Higher Than Heaven", "Once Upon A Time")

  encouragements <- c("I'm proud of ", "keep your ", "be a ", "find strength in ", "worry will never change ", "work for a cause, not for ", "gradtitude turns what we have into ", "good things come to ", "your attitude determines your ", "the only limits in life are ", "find joy in ", "surround yourself with only ", "if oppotunity doesn't knock, build ")

  message(paste0("Running model... and ", sample(encouragements, 1), sample(categories, 1), "!"))

  run_time <- Sys.time()

  ### initial values ----
  # hyper-param for relative freq q (allele) and pi (age class)
  beta <- # actually beta and gamma
    matrix(0,
           nrow = nrow(y),
           ncol = ncol(y),
           dimnames = dimnames(y))

  beta[wildpops, alleles] <-
    matrix(
      rep(1 / nalleles, nalleles),
      nrow = K, # number of wildpops (i.e. collection)
      ncol = sum(nalleles),
      byrow = TRUE,
      dimnames = list(wildpops, alleles)
    ) # genetic part of prior (beta)

  beta[allpops, ages] <-
    matrix(
      if (age_prior == "flat") {
        rep(1/ table(age_class), table(age_class)) / length(unique(age_class))
      } else if (age_prior == "zero_out") {
        rep(1/ table(age_class), table(age_class)) / length(unique(age_class)) * (seq.int(length(age_class)) %in% metadat$age) %>% ifelse(. == 0, 1e-5, .) %>% {. / sum(.)}
      } else if (age_prior == "weak_flat") {
        rep(1/ table(age_class), table(age_class)) / (length(unique(age_class)) * 10)
      } else {
        stop("Wrong specification for argument 'age_prior'. Choices are 'flat', 'zero_out', or 'weak_flat'")
      },
      nrow = K + H,
      ncol = C,
      byrow = TRUE,
      dimnames = list(allpops, ages)
    ) # age part of prior (gamma)

  t_q <- apply(y + beta, 1, function(rw) {
    unlist(tapply(rw, trait_fac, function(betty) {
      if (sum(betty)) {
        betty / sum(betty)
      } else {
        rep(1, length(betty))
      }
    }, simplify = FALSE)[names(nstates)])
  }) # transposed (allele freq)

  t_pi <- t_q[-seq.int(sum(nalleles)), ] # initial pi = gamma

  freq <- matrix(
    0,
    nrow = nrow(x),
    ncol = K + H,
    dimnames = list(row.names(x), allpops)
  )

  if (H > 0 & length(na_i) < nrow(x)) {
    freq[-na_i, hatcheries] <-
      t(sapply(as.integer(metadat$iden[-na_i]) - K, function(m) {
        ans = rep(0L, H)
        ans[m] = 1L
        ans
      }))
  } # only when both reporting groups and sample include hatcheries

  # genotype freq prod f(x_m|q_k)h(a_m|pi_k)
  # rows = indiv, cols = pops
  freq[na_i, wildpops] <- exp(x[na_i,] %*% log(t_q[, wildpops]))

  pPrior <- # alpha, hyper-param for p (pop props)
    array(
      abind::abind(
        lapply(seq(D),
               function(d_idx) {
                 abind::abind(
                   lapply(seq(max(S)),
                          function(s_idx) {
                            matrix(
                              (1/ table(groups[, d_idx])/
                                 max(groups[, d_idx]))[groups[, d_idx]],
                              nrow = W,
                              ncol = K + H,
                              byrow = TRUE
                            )
                          }), along = 0)
               }), along = 0), c(D, max(S), W, K + H))

  metadat$iden[na_i] <- unlist(lapply(na_i, function(m) {
    sample(K, 1, TRUE,
           pPrior[metadat$district[m], metadat$subdist[m], metadat$week[m], seq.int(K)] * freq[m, seq.int(K)])
  }))

  p_zero <- array(
    apply(
      table(metadat[, c("district", "subdist", "week", "iden")]),
      seq(3),
      function(alf) any(alf!= 0) # prop = 0 if no sample
      ), c(D, max(S), W, K+ H))

  p <- aperm(
    apply(
      table(metadat[, c("district" ,"subdist" ,"week" ,"iden")]) + p_zero* pPrior,
      seq(3), rdirich
      ), c(2:4, 1)
    )

  ### parallel chains ----
  chains <- paste0("chain", seq(nchains))
  cl <- parallel::makePSOCKcluster(nchains)
  doParallel::registerDoParallel(cl, cores = nchains)
  if (!is.null(seed)) doRNG::registerDoRNG(seed, once = TRUE)

  outraw <- foreach::foreach(
    ch = chains, .packages = c("magrittr", "tidyr", "dplyr")
    ) %dorng% {

    ppi_out <- vector("list", D)
    for (d_idx in 1:D) {
      ppi_out[[d_idx]] <- vector("list", max(S))
      for (s_idx in 1:S[d_idx]) {
        ppi_out[[d_idx]][[s_idx]] <- vector("list", W)
        for (w_idx in 1:W) {
          ppi_out[[d_idx]][[s_idx]][[w_idx]] <- vector("list", (nreps- nburn) / thin)
        } # w
      } # s
    } # d

    if (iden_output == TRUE) iden_out <- list()
    i <- 1

    ## gibbs loop ##
    for (rep in seq(nreps + nadapt)) {

      metadat$age[na_a] <- unlist(lapply(na_a, function(m) {
        sample(C, 1, TRUE, t_pi[, metadat$iden[m]])
      })) # assign ages

      metadat$iden[na_i] <- unlist(lapply(na_i, function(m) {
        sample(K, 1, TRUE,
              p[metadat$district[m], metadat$subdist[m], metadat$week[m], seq.int(K)] * freq[m, seq.int(K)])
      })) # assign pop iden (only for wild pops)

      p <- aperm(
        apply(
          table(metadat[, c("district", "subdist", "week", "iden")]) + p_zero* pPrior,
          seq(3), rdirich),
        c(2:4, 1))

      if ((cond_gsi & rep %% 10 != 0) | rep <= nadapt) { # cond gsi or adapt stage

        t_pi <- apply(
          beta[, ages] + table(metadat$iden, metadat$age),
          1, rdirich) # pi|a,z,gamma ~ dirich(gamma+ sum(z*a))

      } else {

        x_sum <-
          matrix(
            0L,
            nrow = nrow(y),
            ncol = ncol(y),
            dimnames = dimnames(y)
          )

        x[na_a, seq.int(sum(nalleles) + 1, sum(nstates))] <-
          t(sapply(seq.int(length(na_a)),
                   function(mm) {
                     ag <- metadat$age[na_a][mm]
                     ags <- rep(0L, C)
                     ags[+ag] <- 1L
                     ags
                   })) # updating age info in x

        x_sum[as.integer(sort(unique(metadat$iden))),] <-
          rowsum(x, group = metadat$iden) # colsums for new assignment
        # (this is the original code for tallying. It does not account for new age assignments
        # because x is not updated during sampling. Note: x is updated in Feb 2024)

        # edited version because earlier version did not update x
        # x_sum[as.integer(levels(metadat$iden)), seq.int(sum(nalleles))] <-
        #   apply(x[, seq.int(sum(nalleles))], 2, function(x_t) tapply(x_t, metadat$iden, sum)) %>%
        #   tidyr::replace_na(0) # colsums for new assignment
        #
        # x_sum[, -seq.int(sum(nalleles))] <-
        #   table(metadat$iden, metadat$age) # update age identity

        beta_prm <- y + beta + x_sum # posterior q ~ dirich(b')

        t_q <- apply(beta_prm, 1, function(rw) {
          unlist(tapply(rw, INDEX = trait_fac, FUN = rdirich))
        })

        t_pi <- t_q[-seq.int(sum(nalleles)), ] # pi|a,z,gamma ~ dirich(gamma+ sum(z*a))

        freq[na_i, wildpops] <- exp(x[na_i,] %*% log(t_q[, wildpops]))

      }

      # record output based on keep or not keep burn-ins
      if (rep > nadapt) { # after adaptation stage
        if ((rep-nadapt) > nburn & (rep-nadapt - nburn) %% thin == 0) {

          it <- (rep - nadapt - nburn) / thin
          for (d_idx in 1:D) {
            for (s_idx in 1:S[d_idx]) {
              for (w_idx in 1:W) {
                ppi_out[[i]] <-
                  t_pi %*% diag(p[d_idx, s_idx, w_idx, ]) %>%
                  data.table::as.data.table() %>% # x = ages, y = pops
                  dplyr::mutate(itr = it,
                                agevec = age_class,
                                d = d_idx,
                                s = s_idx,
                                w = w_idx,
                                chain = ch)
                i <- i + 1
              } # w
            } # s
          } # d

          if (iden_output == TRUE) {
            iden_out[[it]] <- metadat$iden
          }

        } # if rep > nburn & (rep-nburn) %% thin == 0
      } # if rep > nadapt

    } # end gibbs loop

    if (iden_output == TRUE) {
      out <- list(data.table::rbindlist(ppi_out) %>% dplyr::arrange(d, s, w),
                  {sapply(iden_out, rbind) %>% t() %>% dplyr::as_tibble()})
      colnames(out[[1]])[seq.int(nrow(groups))] <- row.names(groups)
      return(out)
    } else {
      out <- data.table::rbindlist(ppi_out) %>% dplyr::arrange(d, s, w)
      colnames(out)[seq.int(nrow(groups))] <- row.names(groups)
      return(out)
    }

    } # end parallel chains

  parallel::stopCluster(cl)

  ### prepare output ----

  if (iden_output == TRUE) {
    outraw1 <- lapply(outraw, function(ol) ol[[1]])

    outraw_ids <-
      lapply(outraw, function(ol) ol[[2]]) %>%
      dplyr::bind_rows()

    magma_out <- list(outraw = outraw1, idens = outraw_ids, specs = specs)
  } else {
    magma_out <- list(outraw = outraw, specs = specs)
  }

  if (!is.null(file)) {
    message(paste("Raw MAGMA output saved in", file))
    sapply(seq.int(length(magma_out$outraw)), function(i) tidyfst::export_fst(magma_out$outraw[[i]], path = paste0(file, "/magma_raw_ch", i, ".fst")))
    if (iden_output == TRUE) {
      message(paste("IA posteriors saved in", file, "as magma_idens.fst"))
      tidyfst::export_fst(magma_out$idens, path = paste0(file, "/magmma_iden.fst"))
    }
  }

  print(Sys.time() - run_time)
  message(Sys.time())

  return(magma_out)

}


utils::globalVariables(c("district", "subdist"))














