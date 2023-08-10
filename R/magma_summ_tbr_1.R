#' Summarize big daddy model output step 1
#'
#' Use this function for summarize big daddy model output (e.g., TBR)
#'
#' @param which_dist Function format raw magma output one district at a time.
#'   Identify district as 1, 2, ... Default = NULL will summarize all districts.
#' @param ma_out MAGMA output.
#' @param ma_dat The same MAGMA input data for model run.
#'
#' @return Model output in multiway array and subset of metadata as a list object.
#' @importFrom magrittr %>%
#'
#' @examples
#' # format data
#' wd <- getwd() # path to data folder
#' magma_data <- magmatize_data(wd = wd, save_data = FALSE)
#'
#' # model run
#' magma_out <- magmatize_mdl(magma_data,
#'   nreps = 50, nburn = 25, thin = 1, nchains = 2)
#'
#' # summary step 1
#' tbr1 <- magmatize_summ_bd1(which_dist = 1,
#'   ma_out = magma_out,
#'   ma_dat = magma_data)
#'
#' @export
magmatize_summ_bd1 <- function(which_dist = NULL, ma_out, ma_dat) {

  if (is.null(which_dist)) stop("Must declare which_dist.")
  # if (is.null(which_dist)) which_dist <- unique(ma_dat$metadat$district)

  outraw <- ma_out$outraw
  nreps <- ma_out$specs["nreps"]
  nburn <- ma_out$specs["nburn"]
  thin <- ma_out$specs["thin"]
  nchains <- ma_out$specs["nchains"]
  keep_burn <- ma_out$specs["keep_burn"] == 1

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

  if (is.array(outraw)) {
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
























