/* --------------------------------------------------------------------------------
#
#   Branching process simulator
#   From Alex Perkins (@TAlexPerkins) and lab (http://perkinslab.weebly.com) on March 10, 2020
#   Translated by Sean Wu (slwu89@berkeley.edu)
#   March 2020
#
-------------------------------------------------------------------------------- */

#include "simOutbreak.h"


/* --------------------------------------------------------------------------------
# one iter
-------------------------------------------------------------------------------- */

SEXP onesim_C(
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
  SEXP lnormFlag_r // int (coerce in R)
){

  /* --------------------------------------------------------------------------------
  # initialize variables
  -------------------------------------------------------------------------------- */

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

  // output (zero out)
  SEXP daily_mat = PROTECT(Rf_allocVector(INTSXP, stopSimulationDay));
  int* daily_mat_ptr = INTEGER(daily_mat);

  SEXP case_mat = PROTECT(Rf_allocVector(INTSXP, stopSimulationDay));
  int* case_mat_ptr = INTEGER(case_mat);

  SEXP death_mat = PROTECT(Rf_allocVector(INTSXP, stopSimulationDay));
  int* death_mat_ptr = INTEGER(death_mat);

  for(int i=0; i<stopSimulationDay; i++){
    daily_mat_ptr[i] = 0;
    case_mat_ptr[i] = 0;
    death_mat_ptr[i] = 0;
  }

  SEXP cum_case = Rf_ScalarInteger(0);
  int* cum_case_ptr = INTEGER(cum_case);
  cum_case_ptr[0] = 0;

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

  // zero out arrays
  for(int i=0; i<stopSimulationDay; i++){
    dailyIncidence[i] = 0;
    dailyCases[i] = 0;
    dailyMortality[i] = 0;
  }

  // initiate vector to store timing of parent infections (this will get resized a lot)
  double* time_exp = (double*)malloc(timeImport_len*sizeof(double));
  memcpy(time_exp,timeImport_ptr,timeImport_len*sizeof(double));
  int parents = timeImport_len;

  /* --------------------------------------------------------------------------------
  #   pointers to dynamically allocated memory
  -------------------------------------------------------------------------------- */

  int* number_offspring = NULL;

  double* time_exp_next = NULL;
  int*    time_exp_floor = NULL;

  double* incubation_periods = NULL;

  double* time_det = NULL;
  int*    time_det_floor = NULL;

  double* time_death = NULL;
  int* time_death_floor = NULL;

  int* time_exp_uniq = NULL;
  int* time_exp_dup = NULL;

  int* time_det_uniq = NULL;
  int* time_det_dup = NULL;

  int* time_death_uniq = NULL;
  int* time_death_dup = NULL;


  /* --------------------------------------------------------------------------------
  #   simulation loop
  -------------------------------------------------------------------------------- */
  int sim=0;
  while(parents > 0){
    sim++;
    Rprintf(" --- SIMULATING GENERATION %d --- \n",sim);

    R_CheckUserInterrupt();

    // draw number of offspring from each case (parent)
    number_offspring = (int*)realloc(number_offspring, parents * sizeof(int));
    int tot_offspring = 0;
    for(int i=0; i<parents; i++){
      double R0;
      double u = unif_rand();
      if(u < asympProp){
        R0 = R * asympRFraction;
      } else {
        R0 = R;
      }
      number_offspring[i] = (int)rnbinom_mu(k,R0);
      tot_offspring += number_offspring[i];
    }

    if(tot_offspring > 0){


      /* --------------------------------------------------------------------------------
      #   determine the time of exposure of each offspring and only for parents with offspring
      -------------------------------------------------------------------------------- */

      time_exp_next = (double*)realloc(time_exp_next,tot_offspring);

      int k=0;
      for(int i=0; i<parents; i++){
        // only for parents with offspring
        if(number_offspring[i] > 0){
          for(int j=0; j<number_offspring[i]; j++){
            double offspring_time_exp = time_exp[i] + rlnorm(si_meanlog,si_sdlog);
            // retain only those infections that occur before time is up
            if( (int)floor(offspring_time_exp) <= stopSimulationDay ){
              time_exp_next[k] = offspring_time_exp;
              k++;
            }
          }
        }
      }

      // new vector time_exp
      parents = k;
      time_exp = (double*)realloc(time_exp,parents * sizeof(double));
      memcpy(time_exp,time_exp_next,parents * sizeof(double));


      /* --------------------------------------------------------------------------------
      #   determine incubation period for each offspring
      -------------------------------------------------------------------------------- */

      incubation_periods = (double*)realloc(incubation_periods,parents * sizeof(double));
      for(int i=0; i<parents; i++){
        incubation_periods[i] = rweibull(inc_shape,inc_scale);
      }


      /* --------------------------------------------------------------------------------
      #   determine the time of infection detection of each offspring
      -------------------------------------------------------------------------------- */

      time_det = (double*)realloc(time_det,parents * sizeof(double));

      k=0;
      for(int i=0; i<parents; i++){
        double time_det_indiv = time_exp[i] + incubation_periods[i] + (double)rpois(report_delay_rate);
        // retain only those detected infections that occur before time is up
        if( (int)floor(time_det_indiv) <= stopSimulationDay ){
          time_det[k] = time_det_indiv;
          k++;
        }
      }
      // time_det will have nonsense after position [0,...,time_det_k-1]
      int time_det_k = k;


      /* --------------------------------------------------------------------------------
      #   determine the time of death of each offspring
      -------------------------------------------------------------------------------- */

      time_death = (double*)realloc(time_death,parents * sizeof(double));

      k=0;
      for(int i=0; i<parents; i++){
        double time_death_indiv = time_exp[i] + incubation_periods[i] + rlnorm(symp_to_death_meanlog,symp_to_death_sdlog);
        // retain only deaths that occur before time is up
        if( (int)floor(time_death_indiv) <= stopSimulationDay ){
          time_death[k] = time_death_indiv;
          k++;
        }
      }
      // time_death will have nonsense after position [0,...,time_death_k-1]
      int time_death_k = k;


      /* --------------------------------------------------------------------------------
      #   update vector of daily incidence of detected infections
      -------------------------------------------------------------------------------- */

      time_exp_floor = (int*)realloc(time_exp_floor,parents * sizeof(int));
      for(int i=0; i<parents; i++){
        time_exp_floor[i] = (int)floor(time_exp[i]);
      }

      time_exp_uniq = (int*)realloc(time_exp_uniq, parents * sizeof(int));
      time_exp_dup = (int*)realloc(time_exp_dup, parents * sizeof(int));
      int j, n, flag;

      time_exp_uniq[0] = time_exp_floor[0];
      int count_exp = 1;

      for(j=0; j <= parents-1; ++j) {
          time_exp_dup[j] = 1;
          flag = 1;;
          /*the unique array will always have 'count' elements*/
          for(n=0; n < count_exp; ++n) {
              if(time_exp_floor[j] == time_exp_uniq[n]) {
                  if(j != n){
                    time_exp_dup[n] += 1;
                  }
                  flag = 0;
              }
          }
          if(flag == 1) {
              ++count_exp;
              time_exp_uniq[count_exp-1] = time_exp_floor[j];
          }
      }

      for(int j=0; j<count_exp; j++){
        dailyIncidence[time_exp_uniq[j]] += time_exp_dup[j];
      }


      /* --------------------------------------------------------------------------------
      #   update vector of daily case reporting
      -------------------------------------------------------------------------------- */

      time_det_floor = (int*)realloc(time_det_floor,time_det_k * sizeof(int));
      for(int i=0; i<time_det_k; i++){
        time_det_floor[i] = (int)floor(time_det[i]);
      }

      time_det_uniq = (int*)realloc(time_det_uniq, time_det_k * sizeof(int));
      time_det_dup = (int*)realloc(time_det_dup, time_det_k * sizeof(int));

      time_det_uniq[0] = time_det_floor[0];
      int count_det = 1;

      for(j=0; j <= time_det_k-1; ++j) {
          time_det_dup[j] = 1;
          flag = 1;;
          /*the unique array will always have 'count' elements*/
          for(n=0; n < count_det; ++n) {
              if(time_det_floor[j] == time_det_uniq[n]) {
                  if(j != n){
                    time_det_dup[n] += 1;
                  }
                  flag = 0;
              }
          }
          if(flag == 1) {
              ++count_det;
              time_det_uniq[count_det-1] = time_det_floor[j];
          }
      }

      for(int j=0; j<count_det; j++){
        dailyCases[time_det_uniq[j]] += time_det_dup[j];
      }


      /* --------------------------------------------------------------------------------
      #   update vector of daily mortality
      -------------------------------------------------------------------------------- */

      time_death_floor = (int*)realloc(time_death_floor,time_death_k * sizeof(int));
      for(int i=0; i<time_death_k; i++){
        time_death_floor[i] = (int)floor(time_death[i]);
      }

      time_death_uniq = (int*)realloc(time_death_uniq, time_death_k * sizeof(int));
      time_death_dup = (int*)realloc(time_death_dup, time_death_k * sizeof(int));

      time_death_uniq[0] = time_death_floor[0];
      int count_death = 1;

      for(j=0; j <= time_death_k-1; ++j) {
          time_death_dup[j] = 1;
          flag = 1;;
          /*the unique array will always have 'count' elements*/
          for(n=0; n < count_death; ++n) {
              if(time_death_floor[j] == time_death_uniq[n]) {
                  if(j != n){
                    time_death_dup[n] += 1;
                  }
                  flag = 0;
              }
          }
          if(flag == 1) {
              ++count_death;
              time_death_uniq[count_death-1] = time_death_floor[j];
          }
      }

      for(int j=0; j<count_death; j++){
        dailyMortality[time_death_uniq[j]] += time_death_dup[j];
      }

    } // end if

  } // end while


  // keep track of things that need to be kept track of
  memcpy(daily_mat_ptr,dailyIncidence,stopSimulationDay * sizeof(int));
  memcpy(case_mat_ptr,dailyCases,stopSimulationDay * sizeof(int));
  memcpy(death_mat_ptr,dailyMortality,stopSimulationDay * sizeof(int));

  for(int i=0; i<stopSimulationDay; i++){
    cum_case_ptr[0] += dailyIncidence[i];
  }

  // return a list to R
  SEXP out = PROTECT(Rf_allocVector(VECSXP,4));
  SET_VECTOR_ELT(out,0,daily_mat);
  SET_VECTOR_ELT(out,1,death_mat);
  SET_VECTOR_ELT(out,2,case_mat);
  SET_VECTOR_ELT(out,3,cum_case);

  SEXP names = PROTECT(Rf_allocVector(STRSXP,4));
  SET_STRING_ELT(names,0,Rf_mkChar("daily"));
  SET_STRING_ELT(names,1,Rf_mkChar("death"));
  SET_STRING_ELT(names,2,Rf_mkChar("cases"));
  SET_STRING_ELT(names,3,Rf_mkChar("cum"));

  Rf_namesgets(out,names);

  PutRNGstate();

  free(dailyIncidence);
  free(dailyCases);
  free(dailyMortality);
  free(time_exp);

  free(number_offspring);
  free(time_exp_next);
  free(time_exp_floor);
  free(incubation_periods);
  free(time_det);
  free(time_det_floor);
  free(time_death);
  free(time_death_floor);
  free(time_exp_uniq);
  free(time_exp_dup);
  free(time_det_uniq);
  free(time_det_dup);
  free(time_death_uniq);
  free(time_death_dup);

  UNPROTECT(5);
  return out;
};






















/* --------------------------------------------------------------------------------
# the whole sim
-------------------------------------------------------------------------------- */

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

  /* --------------------------------------------------------------------------------
  # initialize variables
  -------------------------------------------------------------------------------- */

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
  SEXP daily_mat = PROTECT(Rf_allocMatrix(INTSXP, repSims, stopSimulationDay));
  int* daily_mat_ptr = INTEGER(daily_mat);

  SEXP case_mat = PROTECT(Rf_allocMatrix(INTSXP, repSims, stopSimulationDay));
  int* case_mat_ptr = INTEGER(case_mat);

  SEXP death_mat = PROTECT(Rf_allocMatrix(INTSXP, repSims, stopSimulationDay));
  int* death_mat_ptr = INTEGER(death_mat);

  SEXP cum_vec = PROTECT(Rf_allocVector(INTSXP, repSims));
  int* cum_vec_ptr = INTEGER(cum_vec);

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


  /* --------------------------------------------------------------------------------
  #   allocate memory for simulation
  -------------------------------------------------------------------------------- */

  // initiate vector to store timing of parent infections
  double* time_exp = (double*)malloc(timeImport_len*sizeof(double));
  memcpy(time_exp,timeImport_ptr,timeImport_len*sizeof(double));

  // // use this when we replace the variable of parents with offspring
  // double* time_exp_new = (double*)malloc(timeImport_len*sizeof(double));
  // int time_exp_len = timeImport_len;
  //
  // other variables
  // double* R_vec                   = (double*)malloc(time_exp_len*sizeof(double));
  // int* number_offspring           = (int*)malloc(time_exp_len*sizeof(int));
  // double* time_exp_within_bound   = (double*)malloc(sizeof(double)*50);
  // double* incubation_periods      = (double*)malloc(time_exp_len*sizeof(double));
  // double* time_det                = (double*)malloc(time_exp_len*sizeof(double));
  // int* det_within_bound           = (int*)malloc(time_exp_len*sizeof(int));
  // int* time_exp_floor             = (int*)malloc(time_exp_len*sizeof(int));
  // int* time_exp_unique            = (int*)malloc(time_exp_len*sizeof(int));
  // int* time_exp_duplicated        = (int*)malloc(time_exp_len*sizeof(int));
  // int* time_det_unique            = (int*)malloc(50*sizeof(int));
  // int* time_det_duplicated        = (int*)malloc(50*sizeof(int));
  // double* time_death              = (double*)malloc(time_exp_len*sizeof(double));
  // int* death_within_bound         = (int*)malloc(time_exp_len*sizeof(int));
  // int* time_death_unique          = (int*)malloc(50*sizeof(int));
  // int* time_death_duplicated      = (int*)malloc(50*sizeof(int));

  double* time_exp_new = NULL;
  int time_exp_len = timeImport_len;

  double* R_vec                   = NULL;
  int* number_offspring           = NULL;
  double* time_exp_within_bound   = NULL;
  double* incubation_periods      = NULL;
  double* time_det                = NULL;
  int* det_within_bound           = NULL;
  int* time_exp_floor             = NULL;
  int* time_exp_unique            = NULL;
  int* time_exp_duplicated        = NULL;
  int* time_det_unique            = NULL;
  int* time_det_duplicated        = NULL;
  double* time_death              = NULL;
  int* death_within_bound         = NULL;
  int* time_death_unique          = NULL;
  int* time_death_duplicated      = NULL;

  printf("finished allocating memory\n");

  /* --------------------------------------------------------------------------------
  #   Monte Carlo iterations
  -------------------------------------------------------------------------------- */

  // replicate repSims number of times
  for(int rr=0; rr<repSims; rr++){

    // zero out arrays for this MC rep
    for(int i=0; i<stopSimulationDay; i++){
      dailyIncidence[i] = 0;
      dailyCases[i] = 0;
      dailyMortality[i] = 0;
    }

    printf("arrays zeroed out\n");

    /* --------------------------------------------------------------------------------
    #   simulate through generations of transmission until extinct or time is up
    -------------------------------------------------------------------------------- */

    while(time_exp_len > 0){

      printf(" --- SIM ONE GENERATION --- \n");

      // draw a number of offspring for each parent
      R_vec = (double*)realloc(R_vec,sizeof(double)*time_exp_len);
      if(!R_vec){
        printf(" --- --- --- realloc for R_vec failed --- --- --- \n");
      }
      number_offspring = (int*)realloc(number_offspring,sizeof(int)*time_exp_len);
      if(!number_offspring){
        printf(" --- --- --- realloc for number_offspring failed --- --- --- \n");
      }
      int number_offspring_tot = 0; // total num off spring from all parents this gen
      for(int j=0; j<time_exp_len; j++){
        // get parent j's R0
        R_vec[j] = unif_rand();
        if(R_vec[j] < asympProp){
          R_vec[j] = R*asympRFraction;
        } else {
          R_vec[j] = R;
        }
        // sample number of offspring from j
        number_offspring[j] = (int)rnbinom_mu(k,R_vec[j]);
        printf("sampling num offspring of j: %d",j);
        printf(" number_offspring[j]: %d\n",number_offspring[j]);
        number_offspring_tot += number_offspring[j]; // increment tot offspring
      }

      printf("number_offspring_tot %d\n",number_offspring_tot);

      if(number_offspring_tot > 0){

        // remove parents with zero offspring
        // we're not gonna do that in C due to too much memory management, it's slow
        // we'll just check the logic below

        // determine the time of exposure of each offspring
        time_exp_new = (double*)realloc(time_exp_new,number_offspring_tot*sizeof(double));
        if(!time_exp_new){
          printf(" --- --- --- realloc for time_exp_new failed --- --- --- \n");
        }
        int time_exp_new_i = 0;
        if(lnormFlag!=1){

          // go over each parent
          for(int j=0; j<time_exp_len; j++){
            // if produced offspring cases, loop over them
            if(number_offspring[j] > 0){
              for(int k=0; k<number_offspring[j]; k++){
                // k'th offspring of parent j has exposure time equal to time_exp[j] plus a random variate
                time_exp_new[time_exp_new_i] = time_exp[j] + rnorm(si_mean,si_sd);
                printf("sampling time of exposure of each offspring time_exp_new[time_exp_new_i]: %f\n",time_exp_new[time_exp_new_i]);
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
                printf("sampling time of exposure of each offspring time_exp_new[time_exp_new_i]: %f\n",time_exp_new[time_exp_new_i]);
                time_exp_new_i++;
              }
            }
          }

        }
        printf("time_exp_new_i (should equal number_offspring_tot): %d\n",time_exp_new_i);

        // retain only those infection that occur before time is up
        int inf_within_bound = 0;
        time_exp_within_bound = (double*)realloc(time_exp_within_bound,number_offspring_tot*sizeof(double));
        if(!time_exp_within_bound){
          printf(" --- --- --- realloc for time_exp_within_bound failed --- --- --- \n");
        }
        for(int j=0; j<number_offspring_tot; j++){
          if( (int)floor(time_exp_new[j]) <= stopSimulationDay ){
            time_exp_within_bound[inf_within_bound] = time_exp_new[j];
            printf("bound time_exp_new, this one passed: %f\n",time_exp_within_bound[inf_within_bound]);
            inf_within_bound++;
          }
        }
        printf("total number of ones that passed: %d\n",inf_within_bound);
        time_exp = (double*)realloc(time_exp,sizeof(double)*inf_within_bound);
        if(!time_exp){
          printf(" --- --- --- realloc for time_exp failed --- --- --- \n");
        }
        memcpy(time_exp,time_exp_within_bound,sizeof(double)*inf_within_bound);
        time_exp_len = inf_within_bound;

        printf("retain only those infection that occur before time is up, time_exp_len: %d\n",time_exp_len);

        // // C only: update time_exp
        // time_exp = (double*)realloc(time_exp,number_offspring_tot*sizeof(double));
        // memcpy(time_exp,time_exp_new,number_offspring_tot*sizeof(double));
        // time_exp_len = number_offspring_tot;
        //
        // // retain only those infection that occur before time is up
        // int inf_within_bound = 0;
        // double* time_exp_within_bound = (double*)malloc(sizeof(double)*time_exp_len);
        // for(int j=0; j<time_exp_len; j++){
        //   if( (int)floor(time_exp[j]) <= stopSimulationDay ){
        //     time_exp_within_bound[inf_within_bound] = time_exp[j];
        //     inf_within_bound++;
        //   }
        // }
        // time_exp = (double*)realloc(time_exp,sizeof(double)*inf_within_bound);
        // memcpy(time_exp,time_exp_within_bound,sizeof(double)*inf_within_bound);
        // time_exp_len = inf_within_bound;

        // determine incubation period and detection in a single loop
        incubation_periods = (double*)realloc(incubation_periods,time_exp_len*sizeof(double));
        if(!incubation_periods){
          printf(" --- --- --- realloc for incubation_periods failed --- --- --- \n");
        }
        time_det = (double*)realloc(time_det,time_exp_len*sizeof(double));
        if(!time_det){
          printf(" --- --- --- realloc for time_det failed --- --- --- \n");
        }
        for(int j=0; j<time_exp_len; j++){
          incubation_periods[j] = (double)rweibull(inc_shape,inc_scale);
          time_det[j] = time_exp[j] + incubation_periods[j] + rpois(report_delay_rate);
          printf("sampling incubation_periods[j]: %f\n",incubation_periods[j]);
          printf("sampling time_det[j]: %f\n",time_det[j]);
        }

        // retain only those detected infections that occur before time is up
        // put into new int array det_within_bound_i; because after this point in the code we only use the floor version
        int det_within_bound_i = 0;
        det_within_bound = (int*)realloc(det_within_bound,time_exp_len*sizeof(int));
        if(!det_within_bound){
          printf(" --- --- --- realloc for det_within_bound failed --- --- --- \n");
        }
        for(int j=0; j<time_exp_len; j++){
          if( (int)floor(time_det[j]) <= stopSimulationDay ){
            det_within_bound[det_within_bound_i] = (int)floor(time_det[j]);
            printf("bound time_det, this one passed det_within_bound[det_within_bound_i]: %d\n",det_within_bound[det_within_bound_i]);
            det_within_bound_i++;
          }
        }

        printf("retain only those detected infections that occur before time is up, det_within_bound_i: %d\n",det_within_bound_i);

        // integer vector of time_exp
        time_exp_floor = (int*)realloc(time_exp_floor,time_exp_len*sizeof(int));
        if(!time_exp_floor){
          printf(" --- --- --- realloc for time_exp_floor failed --- --- --- \n");
        }
        for(int j=0; j<time_exp_len; j++){
          time_exp_floor[j] = (int)floor(time_exp[j]);
          printf("time_exp_floor[j]: %d\n",time_exp_floor[j]);
        }

        printf("--- grab unique vals of time_exp --- \n");
        // first, need to get unique values in time_exp
        time_exp_unique      = (int*)realloc(time_exp_unique,time_exp_len*sizeof(int));
        if(!time_exp_unique){
          printf(" --- --- --- realloc for time_exp_unique failed --- --- --- \n");
        }
        time_exp_duplicated  = (int*)realloc(time_exp_duplicated,time_exp_len*sizeof(int));
        if(!time_exp_duplicated){
          printf(" --- --- --- realloc for time_exp_duplicated failed --- --- --- \n");
        }
        int flag;
        time_exp_unique[0] = time_exp_floor[0];
        int count_time_exp = 1; /*one element is counted*/

        for(int j=0; j <= time_exp_len-1; ++j) {
          time_exp_duplicated[j] = 1;
          flag = 1;
          /*the unique array will always have 'count' elements*/
          for(int n=0; n < count_time_exp; ++n) {
            if(time_exp_floor[j] == time_exp_unique[n]) {
              if(j != n){
                time_exp_duplicated[n] += 1;
              }
              flag = 0;
            }
          }
          if(flag == 1) {
            ++count_time_exp;
            time_exp_unique[count_time_exp-1] = time_exp_floor[j];
          }
        }
        printf("there are this many unique values count_time_exp: %d\n",count_time_exp);
        for(int j=0; j<count_time_exp; j++){
          printf("time_exp_unique[j]: %d\n",time_exp_unique[j]);
          printf("time_exp_duplicated[j]: %d\n",time_exp_duplicated[j]);
        }

        printf("--- grab unique vals of time_det --- \n");
        // also need to get unique values for time_det
        time_det_unique      = (int*)realloc(time_det_unique,det_within_bound_i*sizeof(int));
        if(!time_det_unique){
          printf(" --- --- --- realloc for time_det_unique failed --- --- --- \n");
        }
        time_det_duplicated  = (int*)realloc(time_det_duplicated,det_within_bound_i*sizeof(int));
        if(!time_det_duplicated){
          printf(" --- --- --- realloc for time_det_duplicated failed --- --- --- \n");
        }
        time_det_unique[0] = det_within_bound[0];
        int count_time_det = 1; /*one element is counted*/

        for(int j=0; j <= det_within_bound_i-1; ++j) {
          time_det_duplicated[j] = 1;
          flag = 1;
          /*the unique array will always have 'count' elements*/
          for(int n=0; n < count_time_det; ++n) {
            if(det_within_bound[j] == time_det_unique[n]) {
              if(j != n){
                time_det_duplicated[n] += 1;
              }
              flag = 0;
            }
          }
          if(flag == 1) {
            ++count_time_det;
            time_det_unique[count_time_det-1] = det_within_bound[j];
          }
        }
        printf("there are this many unique values count_time_det: %d\n",count_time_det);
        for(int j=0; j<count_time_det; j++){
          printf("time_det_unique[j]: %d\n",time_det_unique[j]);
          printf("time_det_duplicated[j]: %d\n",time_det_duplicated[j]);
        }

        // update vector of daily incidence of detected infections
        printf(" --- update vector of daily incidence of detected infections --- \n");
        for(int j=0; j<count_time_exp; j++){
          dailyIncidence[time_exp_unique[j]] += time_exp_duplicated[j];
        }

        // update vector of daily case reporting
        printf(" --- update vector of daily case reporting --- \n");
        for(int j=0; j<count_time_det; j++){
          dailyCases[time_det_unique[j]] += time_det_duplicated[j];
        }

        // determine the time of death of each offspring
        printf(" --- determine the time of death of each offspring --- \n");
        time_death = (double*)realloc(time_death,time_exp_len*sizeof(double));
        if(!time_death){
          printf(" --- --- --- realloc for time_death failed --- --- --- \n");
        }
        for(int j=0; j<time_exp_len; j++){
          time_death[j] = time_exp[j] + incubation_periods[j] + rlnorm(symp_to_death_meanlog,symp_to_death_sdlog);
        }

        // retain only deaths that occur before time is up
        printf(" --- retain only deaths that occur before time is up --- \n");
        int death_within_bound_i = 0;
        death_within_bound = (int*)realloc(death_within_bound,time_exp_len*sizeof(int));
        if(!death_within_bound){
          printf(" --- --- --- realloc for death_within_bound failed --- --- --- \n");
        }
        for(int j=0; j<time_exp_len; j++){
          if( (int)floor(time_death[j]) <= stopSimulationDay ){
            death_within_bound[death_within_bound_i] = (int)floor(time_death[j]);
            printf("bound time_death, this one passed death_within_bound[death_within_bound_i]: %d\n",death_within_bound[death_within_bound_i]);
            death_within_bound_i++;
          }
        }
        printf(" total num deaths retained: %d\n",death_within_bound_i);

        printf("--- grab unique vals of death_within_bound --- \n");
        // need to get unique values for death_within_bound
        time_death_unique      = (int*)realloc(time_death_unique,death_within_bound_i*sizeof(int));
        if(!time_death_unique){
          printf(" --- --- --- realloc for time_death_unique failed --- --- --- \n");
        }
        time_death_duplicated  = (int*)realloc(time_death_duplicated,death_within_bound_i*sizeof(int));
        if(!time_death_duplicated){
          printf(" --- --- --- realloc for time_death_duplicated failed --- --- --- \n");
        }
        time_death_unique[0] = death_within_bound[0];
        int count_time_death = 1; /*one element is counted*/

        for(int j=0; j <= death_within_bound_i-1; ++j) {
          time_death_duplicated[j] = 1;
          flag = 1;
          /*the unique array will always have 'count' elements*/
          for(int n=0; n < count_time_death; ++n) {
            if(death_within_bound[j] == time_death_unique[n]) {
              if(j != n){
                time_death_duplicated[n] += 1;
              }
              flag = 0;
            }
          }
          if(flag == 1) {
            ++count_time_death;
            time_death_unique[count_time_death-1] = death_within_bound[j];
          }
        }
        printf("there are this many unique values count_time_death: %d\n",count_time_death);
        for(int j=0; j<count_time_death; j++){
          printf("time_death_unique[j]: %d\n",time_death_unique[j]);
          printf("time_death_duplicated[j]: %d\n",time_death_duplicated[j]);
        }

        printf(" --- update vector of daily mortality --- \n");
        // update vector of daily mortality
        for(int j=0; j<count_time_death; j++){
          dailyMortality[time_death_unique[j]] += time_death_duplicated[j];
        }
      } // if statement for offspring >0
    } // while loop over sim

    // keep track of things that need to be kept track of
    int cum_inc_tot = 0;
    for(int i=0; i<stopSimulationDay; i++){
      // daily_mat
      daily_mat_ptr[rr + repSims*i] = dailyIncidence[i];

      // case_mat
      case_mat_ptr[rr + repSims*i] = dailyCases[i];

      // dailyMortality
      death_mat_ptr[rr + repSims*i] = dailyMortality[i];

      // cumulative incidence
      cum_inc_tot += dailyIncidence[i];
    }
    cum_vec_ptr[rr] = cum_inc_tot;


  } // end monte carlo reps

  printf(" --- SIMULATION DONE RETURNING TO R --- \n");

  // return a list to R
  SEXP out = PROTECT(Rf_allocVector(VECSXP,4));
  SET_VECTOR_ELT(out,0,daily_mat);
  SET_VECTOR_ELT(out,1,death_mat);
  SET_VECTOR_ELT(out,2,case_mat);
  SET_VECTOR_ELT(out,3,cum_vec);

  SEXP names = PROTECT(Rf_allocVector(STRSXP,4));
  SET_STRING_ELT(names,0,Rf_mkChar("daily"));
  SET_STRING_ELT(names,1,Rf_mkChar("death"));
  SET_STRING_ELT(names,2,Rf_mkChar("cases"));
  SET_STRING_ELT(names,3,Rf_mkChar("cum"));

  Rf_namesgets(out,names);

  PutRNGstate();

  free(dailyIncidence);
  free(dailyCases);
  free(dailyMortality);
  free(time_exp);
  free(time_exp_new);
  free(R_vec);
  free(number_offspring);
  free(time_exp_within_bound);
  free(incubation_periods);
  free(time_det);
  free(det_within_bound);
  free(time_exp_floor);
  free(time_exp_unique);
  free(time_exp_duplicated);
  free(time_det_unique);
  free(time_det_duplicated);
  free(time_death);
  free(death_within_bound);
  free(time_death_unique);
  free(time_death_duplicated);

  UNPROTECT(6);
  return out;
};
