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
  int* dailyIncidence = (int*)malloc(stopSimulationDay*sizeof(int));

  // initiate vector to store daily reporting of autochthonous infections
  int* dailyCases = (int*)malloc(stopSimulationDay*sizeof(int));

  // initiate vector to store daily mortality of autochthonous infections
  int* dailyMortality = (int*)malloc(stopSimulationDay*sizeof(int));

  // replicate repSims number of times
  for(int rr=0; rr<repSims; rr++){

    // zero out arrays for this MC rep
    for(int i=0; i<stopSimulationDay; i++){
      dailyIncidence[i] = 0;
      dailyCases[i] = 0;
      dailyMortality[i] = 0;
    }

    // initiate vector to store timing of parent infections
    double* time_exp = (double*)malloc(timeImport_len*sizeof(double));
    memcpy(time_exp,timeImport_ptr,timeImport_len*sizeof(double));
    int time_exp_len = timeImport_len;

    // use this when we replace the variable of parents with offspring
    double* time_exp_new = (double*)malloc(timeImport_len*sizeof(double));

    // loop through generations of transmission until extinct or time is up
    while(timeImport_len > 0){

      // draw a number of offspring for each parent
      double* R_vec = (double*)malloc(sizeof(double)*time_exp_len);
      int* number_offspring = (int*)malloc(sizeof(int)*time_exp_len);
      int number_offspring_tot = 0; // total num off spring from all parents this gen
      for(int j=0; j<time_exp_len; j++){
        // get parent j's R0
        R_vec[j] = unif_rand();
        if(R_vec[j] < asympProp){
          R_vec[j] = R_vec[j]*asympRFraction;
        }
        // sample number of offspring from j
        number_offspring[j] = (int)rnbinom_mu(k,R_vec[j]);
        number_offspring_tot += number_offspring[j]; // increment tot offspring
      }

      if(number_offspring_tot > 0){

        // remove parents with zero offspring
        // we're not gonna do that in C due to too much memory management, it's slow
        // we'll just check the logic below

        // determine the time of exposure of each offspring
        time_exp_new = (double*)realloc(time_exp_new,number_offspring_tot*sizeof(double));
        int time_exp_new_i = 0;
        if(lnormFlag==1){

          // go over each parent
          for(int j=0; j<time_exp_len; j++){
            // if produced offspring cases, loop over them
            if(number_offspring[j] > 0){
              for(int k=0; k<number_offspring[j]; k++){
                // k'th offspring of parent j has exposure time equal to time_exp[j] plus a random variate
                time_exp_new[time_exp_new_i] = time_exp[j] + rnorm(si_mean,si_sd);
                time_exp_new_i++;
              }
            }
          }

        } else {

          // go over each parent
          for(int j=0; j<time_exp_len; j++){
            // if produced offspring cases, loop over them
            if(number_offspring[j] > 0){
              for(int k=0; k<number_offspring[j]; k++){
                // k'th offspring of parent j has exposure time equal to time_exp[j] plus a random variate
                time_exp_new[time_exp_new_i] = time_exp[j] + rlnorm(si_meanlog,si_sdlog);
                time_exp_new_i++;
              }
            }
          }

        }

        // C only: update time_exp
        time_exp = (double*)realloc(time_exp,number_offspring_tot*sizeof(double));
        memcpy(time_exp,time_exp_new,number_offspring_tot*sizeof(double));
        time_exp_len = number_offspring_tot;

        // retain only those infection that occur before time is up
        int inf_within_bound = 0;
        double* time_exp_within_bound = (double*)malloc(sizeof(double)*time_exp_len);
        for(int j=0; j<time_exp_len; j++){
          if( (int)floor(time_exp[j]) <= stopSimulationDay ){
            time_exp_within_bound[inf_within_bound] = time_exp[j];
            inf_within_bound++;
          }
        }
        time_exp = (double*)realloc(time_exp,sizeof(double)*inf_within_bound);
        memcpy(time_exp,time_exp_within_bound,sizeof(double)*inf_within_bound);
        time_exp_len = inf_within_bound;

        // determine incubation period and detection in a single loop
        double* incubation_periods = (double*)malloc(time_exp_len*sizeof(double));
        double* time_det = (double*)malloc(time_exp_len*sizeof(double));
        for(int j=0; j<time_exp_len; j++){
          incubation_periods[j] = rweibull(inc_shape,inc_scale);
          time_det[j] = time_exp[j] + incubation_periods[j] + rpois(report_delay_rate);
        }

        // retain only those detected infections that occur before time is up
        int det_within_bound_i = 0;
        double* det_within_bound = (double*)malloc(sizeof(double)*time_exp_len);
        for(int j=0; j<time_exp_len; j++){
          if( (int)floor(time_det[j]) <= stopSimulationDay ){
            det_within_bound[det_within_bound_i] = time_det[j];
            det_within_bound_i++;
          }
        }
        time_det = (double*)realloc(time_det,sizeof(double)*det_within_bound_i);
        memcpy(time_det,det_within_bound,sizeof(double)*det_within_bound_i);

        // update vector of daily incidence of detected infections






      }

      // free memory (make more efficient later)
      free(R_vec);
      free(number_offspring);
    }





    free(time_exp);
    free(time_exp_new);
  } // end monte carlo reps


  PutRNGstate();

  free(dailyIncidence);
  free(dailyCases);
  free(dailyMortality);


};
