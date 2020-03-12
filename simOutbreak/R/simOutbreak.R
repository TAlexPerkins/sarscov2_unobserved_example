# -------------------------------------------------------------------------------- $
#
#   Branching process simulator
#   From Alex Perkins (@TAlexPerkins) and lab (http://perkinslab.weebly.com) on March 10, 2020
#   Translated by Sean Wu (slwu89@berkeley.edu)
#   March 2020
#
# -------------------------------------------------------------------------------- #

#' Branching Process Outbreak Simulator
#'
#'
#' @useDynLib simOutbreak simOutbreak_C
#' @export
simOutbreak <- function(
  timeImport,
  R = 1.97,
  k = 1e3,
  si_mean = 4.56,
  si_sd = 0.95,
  inc_shape = 1.88,
  inc_scale = 7.97,
  symp_to_death_mean = 22.3,
  symp_to_death_sd = 0.42 * 22.3,
  report_delay_rate = 3,
  stopSimulationDay = 68,
  repSims = 1e2,
  asympProp = 0.35,
  asympRFraction = 0.5,
  lnormFlag = TRUE
){
  stopifnot(is.logical(lnormFlag))
  .Call(
    simOutbreak_C,
    as.numeric(timeImport),
    as.numeric(R),
    as.numeric(k),
    as.numeric(si_mean),
    as.numeric(si_sd),
    as.numeric(inc_shape),
    as.numeric(inc_scale),
    as.numeric(symp_to_death_mean),
    as.numeric(symp_to_death_sd),
    as.numeric(report_delay_rate),
    as.integer(stopSimulationDay),
    as.integer(repSims),
    as.numeric(asympProp),
    as.numeric(asympRFraction),
    as.integer(lnormFlag)
  )
}
