## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----setup--------------------------------------------------------------------
library(magma.gsi)
# devtools::load_all()


## ----groupvec, message=FALSE--------------------------------------------------
readr::read_table("data/group_namesFAKE.txt")

readr::read_table("data/groupsFAKE.txt")


## ----harvest, message=FALSE---------------------------------------------------
readr::read_table("data/harvestFAKE.txt")


## ----metadat, message=FALSE---------------------------------------------------
head(readr::read_table("data/metadata.txt"))


## ----run-magma-data-----------------------------------------------------------
yomamafat <-
  magmatize_data(wd = getwd(),
                 age_classes = c(13, 21, 22, 23, 31, 32, "other"),
                 fishery = NULL,
                 loci_names = NULL,
                 save_data = FALSE)


## ----message=FALSE------------------------------------------------------------
str(yomamafat)


## ----run-magama---------------------------------------------------------------
freak_out <-
  magmatize_mdl(dat_in = yomamafat,
                nreps = 50, nburn = 25, thin = 1, nchains = 2, nadapt = 0,
                keep_burn = TRUE, age_prior = "zero_out",
                cond_gsi = TRUE, file = NULL, seed = NULL, iden_output = TRUE)


## ----format-magama------------------------------------------------------------
magma_summ <-
  magmatize_summ(which_dist = 1,
                 ma_out = freak_out,
                 ma_dat = yomamafat,
                 summ_level = "district")


## ----age-summ-----------------------------------------------------------------
magma_summ$age_summ


## ----age-prop-----------------------------------------------------------------
magma_summ$age_prop


## ----trace-plot, fig.height=12, fig.width=8, out.width="100%", fig.cap="Trace plots for age composition."----
magmatize_tr_plot(magma_summ$age_prop$D1_Koyukuk, nburn = 25)


## ----ia-----------------------------------------------------------------------
magma_indiv <- magmatize_indiv(ma_out = freak_out, ma_dat = yomamafat, out_repunit = TRUE)

magma_indiv


