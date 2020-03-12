/* --------------------------------------------------------------------------------
#
#   Branching process simulator
#   From Alex Perkins (@TAlexPerkins) and lab (http://perkinslab.weebly.com) on March 10, 2020
#   Translated by Sean Wu (slwu89@berkeley.edu)
#   March 2020
#
-------------------------------------------------------------------------------- */

#include "simOutbreak.h"

// branching process simulator
SEXP simOutbreak_C(
  SEXP timeImport,
  SEXP R_r,
  SEXP k_r,
  SEXP si_mean_r,
  SEXP si_sd_r,
  SEXP inc_shape_r,
  SEXP inc_scale_r,
  SEXP symp_to_death_mean_r,
  SEXP symp_to_death_sd_r,
  SEXP report_delay_rate_r,
  SEXP stopSimulationDay_r,
  SEXP repSims_r,
  SEXP asympProp_r,
  SEXP asympRFraction_r,
  SEXP lnormFlag_r
){

  GetRNGstate();

  // stuff about timeImport
  int timeImport_len = Rf_length(timeImport);
  double* timeImport_ptr = REAL(timeImport);

  // scalar parameters
  double R = Rf_asReal(R_r);
  double k = Rf_asReal(k_r);
  double si_mean = Rf_asReal(si_mean_r);
  double si_sd = Rf_asReal(si_sd_r);
  double inc_shape = Rf_asReal(inc_shape_r);
  double inc_scale = Rf_asReal(inc_scale_r);
  double symp_to_death_mean = Rf_asReal(symp_to_death_mean_r);
  double symp_to_death_sd = Rf_asReal(symp_to_death_sd_r);
  double report_delay_rate = Rf_asReal(report_delay_rate_r);
  int    stopSimulationDay = Rf_asInteger(stopSimulationDay_r);
  int    repSims = Rf_asInteger(repSims_r);
  double asympProp = Rf_asReal(asympProp_r);
  double asympRFraction = Rf_asReal(asympRFraction_r);

  int    lnormFlag = Rf_asInteger(lnormFlag_r);

  // decide what to keep track of across replicate simulations
  SEXP daily_mat = Rf_allocMatrix(INTSXP, repSims, stopSimulationDay);
  SEXP case_mat = Rf_allocMatrix(INTSXP, repSims, stopSimulationDay);
  SEXP death_mat = Rf_allocMatrix(INTSXP, repSims, stopSimulationDay);

  SEXP cum_vec = Rf_allocVector(INTSXP, repSims);

  // calculate meanlog and sdlogs for lnorm distributions
  double si_sdlog, si_meanlog;
  if(lnormFlag==1){
    si_sdlog = sqrt(log( pow(si_sd/si_mean,2) + 1. ));
    si_meanlog = log(si_mean) - 0.5*pow(si_sdlog,2);
  }

  double symp_to_death_sdlog = sqrt(log( pow(symp_to_death_sd/symp_to_death_mean,2) + 1.));
  double symp_to_death_meanlog = log(symp_to_death_mean) - 0.5*pow(symp_to_death_sdlog,2);

  // initiate vector to store daily incidence of autochthonous infections
  int* dailyIncidence = calloc(stopSimulationDay,sizeof(int));

  // initiate vector to store daily reporting of autochthonous infections
  int* dailyCases = calloc(stopSimulationDay,sizeof(int));

  // initiate vector to store daily mortality of autochthonous infections
  int* dailyMortality = calloc(stopSimulationDay,sizeof(int));

  // replicate repSims number of times
  for(int rr=0; rr<repSims; rr++){

    // initiate vector to store timing of parent infections
    int* time_exp = calloc(timeImport_len,sizeof(int));


  }


  PutRNGstate();

};
