
#' Individual group membership assignment summary
#'
#' @param ma_out MAGMA model output object name.
#' @param ma_dat MAGMA input data object name.
#' @param out_repunit Logical (default = `FALSE`). Option to output the summaries in populations (`FALSE`) or combined in reporting groups (`TRUE`).
#'  * NOTE: combining populations into reporting group recommend only for single district or multi-district with the same reporting groups.
#'
#' @return A tibble contains probabilities for reporting group memberships of each individual of MAGMA metadata.
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#' wd <- getwd() # path to data folder
#' magma_data <- magmatize_data(wd = wd, save_data = FALSE)
#' magma_out <- magmatize_mdl(magma_data, nreps = 50, nburn = 25, thin = 1, nchains = 2)
#' magma_assn <- magmatize_indiv(magma_out, magma_datt, out_repunit = TRUE)
#' }
#'
magmatize_indiv <- function(ma_out, ma_dat, out_repunit = FALSE) {

  nreps <- ma_out$specs["nreps"]
  nburn <- ma_out$specs["nburn"]
  thin <- ma_out$specs["thin"]
  # nchains <- ma_out$specs["nchains"]
  keep_burn <- ma_out$specs["keep_burn"] == 1

  keep_list <- ((nburn*keep_burn + 1):(nreps - nburn*isFALSE(keep_burn)))[!((nburn*keep_burn + 1):(nreps - nburn*isFALSE(keep_burn))) %% thin] / thin

  idens <- ma_out$idens[keep_list, ]

  if (isTRUE(out_repunit)) {
    message("Combining populations using reporting groups of district 1.")

    pi <- apply(ma_out$idens, 2,
                function (idens) {
                  factor(idens, levels = seq(length(ma_dat$groups[, 1]))) %>%
                    table(.) %>%
                    tapply(., ma_dat$groups[, 1], sum) %>%
                    prop.table(.)
                })

    tidyr::tibble(indiv = rownames(ma_dat$x)) %>%
      dplyr::bind_cols({
        t(pi) %>%
          as.data.frame() %>%
          stats::setNames(ma_dat$group_names[, 1])
      })
  } else {
    pi <- apply(ma_out$idens, 2,
                function (idens) {
                  factor(idens, levels = seq(length(ma_dat$groups[, 1]))) %>%
                    table(.) %>%
                    prop.table(.)
                })

    tidyr::tibble(indiv = rownames(ma_dat$x)) %>%
      dplyr::bind_cols({
        t(pi) %>%
          as.data.frame() %>%
          stats::setNames(c(ma_dat$wildpops, ma_dat$hatcheries))
      })
  }

}


utils::globalVariables(c("indiv"))


