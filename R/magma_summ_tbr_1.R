
#' TBR output format wrapper step 1
#'
#' @param which_dist Function format raw magma output one district at a time.
#'   Identify district as 1, 2, ... Default = NULL will summarize all districts.
#' @param outraw MAGMA output.
#' @param ma_dat The same MAGMA input data for model run.
#' @param nreps The same as *nreps* in MAGMA model run.
#' @param nburn The same as *nburn* in MAGMA model run.
#' @param thin The same as *thin* in MAGMA model run.
#' @param nchains The same as *nchains* in MAGMA model run.
#' @param keep_burn The same as *keep_burn* in MAGMA model run.
#' @param summ_level Summarize at district or subdistrict level.
#' @param type Identify "pop" or "age" to summarize populations or age class.
#' @param malia = TRUE if the output is from program malia.
#'
#' @return Model output in multiway array and subset of metadata as a list object.
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
#' # summary step 1
#' tbr1 <- magmatize_summ_tbr1(which_dist = 3,
#'   outraw = magma_out,
#'   ma_dat = magma_data,
#'   nreps = 50, nburn = 25, thin = 1, nchains = 3)
#'
magmatize_summ_tbr1 <- function(which_dist = NULL, outraw, ma_dat, nreps, nburn, thin, nchains, keep_burn = FALSE, malia = FALSE) {#, summ_level, type = NULL) {

  # if (is.null(which_dist) | is.null(type)) stop("Must declare which_dist and a type.")
  if (is.null(which_dist)) stop("Must declare which_dist.")

  if (is.null(which_dist)) which_dist <- unique(ma_dat$metadat$district)

  sub_dat <- list(
    # x = ma_dat$x[which(ma_dat$metadat$district %in% which_dist), ],
    # y = ma_dat$y,

    metadat = ma_dat$metadat %>%
      dplyr::filter(district %in% which_dist) %>%
      { if (length(which_dist) == 1) dplyr::mutate(., district = 1) else . },
    harvest = ma_dat$harvest %>%
      dplyr::filter(DISTRICT %in% which_dist) %>%
      { if (length(which_dist) == 1) dplyr::mutate(., DISTRICT = 1) else . },

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
    stat_weeks = ma_dat$stat_weeks,

    nreps = nreps,
    nburn = nburn,
    thin = thin,
    nchains = nchains,
    keep_burn = keep_burn#,
    # summ_level = summ_level,
    # type = type
  )

  S <- ma_dat$metadat %>%
    dplyr::group_by(district) %>%
    dplyr::summarise(S = dplyr::n_distinct(subdist), .groups = "drop") %>%
    dplyr::pull(S) # number of subdistricts
  W <- dplyr::n_distinct(ma_dat$metadat$week) # number of weeks
  # KH <- nrow(ma_dat$y)
  KH <- length(c(ma_dat$wildpops, ma_dat$hatcheries))

  # organization of outraw:
  # [[chain]][[dist]][[sub]][[week]][age, pop, itr]

  if (isTRUE(malia)) {
    holder <- array(outraw[ , , , , which_dist, ], dim = c(ma_dat$C* (nreps- nburn*isFALSE(keep_burn))/ thin, KH + 2, W, S[which_dist], 1, nchains))
  } else {
    holder <-
      lapply(outraw, function(ch) {
        lapply(ch[[which_dist]], function(subdis) {
          lapply(subdis, function(wk) {
            wk %>% dplyr::bind_rows() %>% unlist()
          }) %>% unlist()
        }) %>% unlist()
      })  %>% unlist() %>%
      array(., dim = c(ma_dat$C* (nreps- nburn*isFALSE(keep_burn))/ thin, KH + 2, W, S[which_dist], 1, nchains))
  }

  out1 <- list(holder = holder,
               sub_dat = sub_dat)

  return(out1)

}
























