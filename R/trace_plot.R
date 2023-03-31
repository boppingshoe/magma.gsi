
#' Plot MCMC trace
#'
#' @param obj Trace from the model output.
#' @param nburn Number of burn-in you set up when you ran the model.
#'   Default is 0 if you didn't save the burn-ins (keep_burn = FALSE).
#' @param thin Number of thinning you set up when you ran the model.
#'   Default is 1 (no thinning).
#' @param name_order Arrange the reporting groups as you wish. Leave it empty
#'   if you want to accept the default.
#'
#' @return Trace plot in ggplot

#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' # format data
#' wd <- "D:/bobby_adfg/projects/magma/test_TBR" # path to data folder
#' magma_data <- magmatize_data(wd = wd, save_data = FALSE)
#'
#' # model run
#' magma_out <- magmatize_mdl(magma_data,
#'   nreps = 50, nburn = 25, thin = 1, nchains = 3, tbr = TRUE)
#'
#' # summary
#' magma_summ <- magmatize_summ_tbr(which_dist = 2,
#'   outraw = magma_out,
#'   ma_dat = magma_data,
#'   nreps = 50, nburn = 25, thin = 1, nchains = 3,
#'   summ_level = "district")
#'
#' # trace plot
#' tr_plot(obj = magma_summ$pop_prop[[1]])
#'
tr_plot <- function (obj, nburn = 0, thin = 1, name_order = NULL) {

  if (is.null(name_order)) {
    name_order <- dplyr::select(obj, -c(itr, chain)) %>% colnames()
  }

  if ("grpvec" %in% names(obj)) obj <- dplyr::select(obj, -grpvec)

  tidyr::pivot_longer({{ obj }}, cols = -c(chain, itr)) %>%
    dplyr::mutate(name = factor(name, levels = name_order)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(x = itr, y = value, color = factor(chain))) +
    {if (nburn > 0) {
      ggplot2::annotate("rect", fill = "red", alpha = 0.15,
                        xmin = 0, xmax = nburn/thin, ymin = -Inf, ymax = Inf)
    }} +
    ggplot2::facet_grid(name ~ ., scales = "free") +
    ggplot2::labs(color = "MC chain")

} # nburn = 0 if keep_burn = FALSE

utils::globalVariables(c("chain", "itr", "name", "value"))




