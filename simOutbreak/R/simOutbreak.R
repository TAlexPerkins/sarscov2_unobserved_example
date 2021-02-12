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
#' This function works exactly the same as the R version.
#' If both are loaded in the same R session, use this one via \code{simOutbreak::simOutbreak}
#'
#' @param timeImport timing of each imported infection
#' @param R reproduction number
#' @param k dispersion parameter
#' @param si_mean mean of serial interval distribution
#' @param si_sd standard deviation of serial interval distribution
#' @param inc_shape  hape parameter of incubation period distribution (Weibull)
#' @param inc_scale scale parameter of incubation period distribution (Weibull)
#' @param symp_to_death_mean  mean of time between symptom onset and death
#' @param symp_to_death_sd std. dev. of time between symptom onset and death
#' @param report_delay_rate mean delay between symptoms and reporting
#' @param stopSimulationDay day of year since Jan 1 when simulation stops
#' @param repSims number of replicates to produce for these parameters
#' @param asympProp proportion of infections that are asymptomatic
#' @param asympRFractionrelative infectiousness of asymptomatics
#' @param lnormFlag toggles whether serial interval distribution is lognormal
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
  stopifnot(length(timeImport)>=1)
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
