/* --------------------------------------------------------------------------------
#
#   Branching process simulator
#   From Alex Perkins (@TAlexPerkins) and lab (http://perkinslab.weebly.com) on March 10, 2020
#   Translated by Sean Wu (slwu89@berkeley.edu)
#   March 2020
#
-------------------------------------------------------------------------------- */

#ifndef BPSIM_H
#define BPSIM_H

#include <stdlib.h>
#include <string.h>

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include <R_ext/Utils.h> // for user interrupt checking

// branching process simulator
// all pars are doubles unless otherwise specified
SEXP simOutbreak_C(
  SEXP timeImport, // vector of ints
  SEXP R_r,
  SEXP k_r,
  SEXP si_mean_r,
  SEXP si_sd_r,
  SEXP inc_shape_r,
  SEXP inc_scale_r,
  SEXP symp_to_death_mean_r,
  SEXP symp_to_death_sd_r,
  SEXP report_delay_rate_r,
  SEXP stopSimulationDay_r, // int
  SEXP repSims_r, // int
  SEXP asympProp_r,
  SEXP asympRFraction_r,
  SEXP lnormFlag // int (coerce in R)
);

#endif
