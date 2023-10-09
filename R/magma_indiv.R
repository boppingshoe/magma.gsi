
#' Individual group membership assignment likelihood
#'
#' @param magma_dat MAGMA data formatted using `magmatize_data()`.
#' @param out_repunit Logical (default = `FALSE`). Option to output the likelihood in populations (`FALSE`) or combined in reporting groups (`TRUE`).
#'  * NOTE: combining populations into reporting group only works for single district or multi-district with the same reporting groups.
#'
#' @return A tibble contains likelihood of reporting group memberships for each individual of MAGMA metadata.
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' wd <- getwd() # path to data folder
#' magma_data <- magmatize_data(wd = wd, save_data = FALSE)
#' scaled_like <- magmatize_indiv(magma_data)
#' }
#'
#' @export
magmatize_indiv <- function(magma_dat, out_repunit = FALSE) {
  message("Calculating likelihood...")

  x <- as.matrix(magma_dat$x) # mixture
  y <- as.matrix(magma_dat$y) # base
  metadat <- magma_dat$metadat # age and iden info
  nstates <- magma_dat$nstates # number of allele types and age classes
  nalleles <- magma_dat$nalleles
  wildpops <- magma_dat$wildpops
  C <- magma_dat$C # number of age classes
  K <- length(wildpops)
  if (is.null(magma_dat$hatcheries)) {
    hatcheries <- NULL
    H <- 0
  } else {
    hatcheries <- magma_dat$hatcheries
    H <- length(hatcheries)
  }

  allpops <- c(wildpops, hatcheries)
  age_class <- magma_dat$age_class # vector id for age classes

  na_i <- which(is.na(metadat$iden))

  trait_fac <- factor(rep(names(nstates), nstates), levels = names(nstates))
  states <- paste0(paste0(trait_fac, "_"), unlist(lapply(nstates, seq.int)))
  alleles <- states[seq.int(sum(nalleles))] # allele types
  ages <- states[-seq.int(sum(nalleles))] # age classes

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
      rep(1/ table(age_class), table(age_class)) / max(age_class),
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

  freq <- matrix(
    0,
    nrow = nrow(x),
    ncol = K + H,
    dimnames = list(rownames(x), allpops)
  )

  if (H > 0 & length(na_i) < nrow(x)) {
    freq[-na_i, hatcheries] <-
      t(sapply(as.integer(metadat$iden[-na_i]) - K, function(m) {
        ans = rep(0L, H)
        ans[m] = 1L
        ans
      }))
  } # only when both reporting groups and sample include hatcheries

  freq[na_i, wildpops] <- exp(x[na_i,] %*% log(t_q[, wildpops]))

  scaled_like <-
    apply(freq, 1, prop.table)

  if (isTRUE(out_repunit)) {
    message("Combining populations using reporting groups of district 1.")
    grp_nms <- magma_dat$group_names[, 1]

    scaled_like <-
      scaled_like %>%
      rowsum(., magma_dat$groups[, 1]) %>%
      t() %>%
      dplyr::as_tibble(., .name_repair = "minimal") %>%
      purrr::set_names(grp_nms) %>%
      dplyr::mutate(indiv = rownames(magma_dat$metadat))
  } else {
    scaled_like <-
      t(scaled_like) %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(indiv = rownames(magma_dat$metadat))
  }

  return(scaled_like %>% dplyr::relocate(indiv))

}


utils::globalVariables(c("indiv"))







