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
# branching process simulator
-------------------------------------------------------------------------------- */

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

  // output (zero out)
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
  si_sdlog = sqrt(log( pow(si_sd/si_mean,2.) + 1. ));
  si_meanlog = log(si_mean) - 0.5*pow(si_sdlog,2.);

  double symp_to_death_sdlog = sqrt(log( pow(symp_to_death_sd/symp_to_death_mean,2.) + 1.));
  double symp_to_death_meanlog = log(symp_to_death_mean) - 0.5*pow(symp_to_death_sdlog,2.);

  // initiate vector to store daily incidence of autochthonous infections
  int* dailyIncidence = (int*)R_alloc(stopSimulationDay,sizeof(int));

  // initiate vector to store daily reporting of autochthonous infections
  int* dailyCases = (int*)R_alloc(stopSimulationDay, sizeof(int));

  // initiate vector to store daily mortality of autochthonous infections
  int* dailyMortality = (int*)R_alloc(stopSimulationDay, sizeof(int));


  /* --------------------------------------------------------------------------------
  #   Monte Carlo simulation
  -------------------------------------------------------------------------------- */

  for(int rr=0; rr<repSims; rr++){

    // zero out arrays
    for(int i=0; i<stopSimulationDay; i++){
      dailyIncidence[i] = 0;
      dailyCases[i] = 0;
      dailyMortality[i] = 0;
    }

    // initiate vector to store timing of parent infections (this will get resized a lot)
    double* time_exp = (double*)Calloc(timeImport_len,double);
    memmove(time_exp,timeImport_ptr,timeImport_len*sizeof(double));
    int parents = timeImport_len;

    /* --------------------------------------------------------------------------------
    #   simulation loop
    -------------------------------------------------------------------------------- */

    while(parents > 0){

      R_CheckUserInterrupt();

      // draw number of offspring from each case (parent)
      int* number_offspring = (int*)Calloc(parents,int);
      int tot_offspring = 0;
      for(int i=0; i<parents; i++){
        double R0;
        double u = runif(0.,1.);
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

        double* time_exp_next = (double*)Calloc(tot_offspring,double);

        int k;
        if(lnormFlag==0){

          k=0;
          for(int i=0; i<parents; i++){
            // only for parents with offspring
            if(number_offspring[i] > 0){
              for(int j=0; j<number_offspring[i]; j++){
                double offspring_time_exp = time_exp[i] + (double)rnorm(si_mean,si_sd);
                // retain only those infections that occur before time is up
                if( (int)floor(offspring_time_exp) <= stopSimulationDay ){
                  time_exp_next[k] = offspring_time_exp;
                  k++;
                }
              }
            }
          }

        } else {

          k=0;
          for(int i=0; i<parents; i++){
            // only for parents with offspring
            if(number_offspring[i] > 0){
              for(int j=0; j<number_offspring[i]; j++){
                double offspring_time_exp = time_exp[i] + (double)rlnorm(si_meanlog,si_sdlog);
                // retain only those infections that occur before time is up
                if( (int)floor(offspring_time_exp) <= stopSimulationDay ){
                  time_exp_next[k] = offspring_time_exp;
                  k++;
                }
              }
            }
          }

        }

        // new vector time_exp
        parents = k;
        time_exp = (double*)Realloc(time_exp,parents,double);
        memmove(time_exp,time_exp_next,parents * sizeof(double));


        /* --------------------------------------------------------------------------------
        #   determine incubation period for each offspring
        -------------------------------------------------------------------------------- */

        double* incubation_periods = (double*)Calloc(parents,double);
        for(int i=0; i<parents; i++){
          incubation_periods[i] = (double)rweibull(inc_shape,inc_scale);
        }


        /* --------------------------------------------------------------------------------
        #   determine the time of infection detection of each offspring
        -------------------------------------------------------------------------------- */

        double* time_det = (double*)Calloc(parents,double);

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

        double* time_death = (double*)Calloc(parents,double);

        k=0;
        for(int i=0; i<parents; i++){
          double time_death_indiv = time_exp[i] + incubation_periods[i] + (double)rlnorm(symp_to_death_meanlog,symp_to_death_sdlog);
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

        int* time_exp_floor = (int*)Calloc(parents,int);

        for(int i=0; i<parents; i++){
          time_exp_floor[i] = (int)floor(time_exp[i]);
        }

        int* time_exp_uniq = (int*)Calloc(parents,int);
        int* time_exp_dup = (int*)Calloc(parents,int);

        time_exp_uniq[0] = time_exp_floor[0];
        int count_exp = 1;

        for(int j=0; j <= parents-1; ++j) {
            time_exp_dup[j] = 1;
            int flag = 1;
            /*the unique array will always have 'count' elements*/
            for(int n=0; n < count_exp; ++n) {
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

        int* time_det_floor = (int*)Calloc(time_det_k,int);

        for(int i=0; i<time_det_k; i++){
          time_det_floor[i] = (int)floor(time_det[i]);
        }

        int* time_det_uniq = (int*)Calloc(time_det_k,int);
        int* time_det_dup = (int*)Calloc(time_det_k,int);

        time_det_uniq[0] = (int)floor(time_det[0]);
        int count_det = 1;

        for(int j=0; j <= time_det_k-1; ++j) {
            time_det_dup[j] = 1;
            int flag = 1;
            /*the unique array will always have 'count' elements*/
            for(int n=0; n < count_det; ++n) {
                  if((int)floor(time_det[j]) == time_det_uniq[n]) {
                    if(j != n){
                      time_det_dup[n] += 1;
                    }
                    flag = 0;
                }
            }
            if(flag == 1) {
                ++count_det;
                time_det_uniq[count_det-1] = (int)floor(time_det[j]);
            }
        }

        for(int j=0; j<count_det; j++){
          dailyCases[time_det_uniq[j]] += time_det_dup[j];
        }


        /* --------------------------------------------------------------------------------
        #   update vector of daily mortality
        -------------------------------------------------------------------------------- */

        int* time_death_floor = (int*)Calloc(time_death_k,int);

        for(int i=0; i<time_death_k; i++){
          time_death_floor[i] = (int)floor(time_death[i]);
        }

        int* time_death_uniq = (int*)Calloc(time_death_k,int);
        int* time_death_dup = (int*)Calloc(time_death_k,int);

        time_death_uniq[0] = time_death_floor[0];
        int count_death = 1;

        for(int j=0; j <= time_death_k-1; ++j) {
            time_death_dup[j] = 1;
            int flag = 1;
            /*the unique array will always have 'count' elements*/
            for(int n=0; n < count_death; ++n) {
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

        Free(time_exp_next);
        Free(incubation_periods);
        Free(time_det);
        Free(time_death);

        Free(time_exp_floor);
        Free(time_exp_uniq);
        Free(time_exp_dup);

        Free(time_det_floor);
        Free(time_det_uniq);
        Free(time_det_dup);

        Free(time_death_floor);
        Free(time_death_uniq);
        Free(time_death_dup);

      } // end if

      Free(number_offspring);

    } // end while

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

    Free(time_exp);

  } // end MC

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

  UNPROTECT(6);
  return out;
};





// /* --------------------------------------------------------------------------------
// # one iter
// -------------------------------------------------------------------------------- */
//
// SEXP onesim_C(
//   SEXP timeImport, // vector of ints
//   SEXP R_r,
//   SEXP k_r,
//   SEXP si_mean_r,
//   SEXP si_sd_r,
//   SEXP inc_shape_r,
//   SEXP inc_scale_r,
//   SEXP symp_to_death_mean_r,
//   SEXP symp_to_death_sd_r,
//   SEXP report_delay_rate_r,
//   SEXP stopSimulationDay_r, // int
//   SEXP repSims_r, // int
//   SEXP asympProp_r,
//   SEXP asympRFraction_r,
//   SEXP lnormFlag_r // int (coerce in R)
// ){
//
//   /* --------------------------------------------------------------------------------
//   # initialize variables
//   -------------------------------------------------------------------------------- */
//
//   GetRNGstate();
//
//   // stuff about timeImport
//   int timeImport_len = Rf_length(timeImport);
//   double* timeImport_ptr = REAL(timeImport);
//
//   // scalar parameters
//   double R = Rf_asReal(R_r);
//   double k = Rf_asReal(k_r);
//   double si_mean = Rf_asReal(si_mean_r);
//   double si_sd = Rf_asReal(si_sd_r);
//   double inc_shape = Rf_asReal(inc_shape_r);
//   double inc_scale = Rf_asReal(inc_scale_r);
//   double symp_to_death_mean = Rf_asReal(symp_to_death_mean_r);
//   double symp_to_death_sd = Rf_asReal(symp_to_death_sd_r);
//   double report_delay_rate = Rf_asReal(report_delay_rate_r);
//   int    stopSimulationDay = Rf_asInteger(stopSimulationDay_r);
//   int    repSims = Rf_asInteger(repSims_r);
//   double asympProp = Rf_asReal(asympProp_r);
//   double asympRFraction = Rf_asReal(asympRFraction_r);
//
//   int    lnormFlag = Rf_asInteger(lnormFlag_r);
//
//   // output (zero out)
//   SEXP daily_mat = PROTECT(Rf_allocVector(INTSXP, stopSimulationDay));
//   int* daily_mat_ptr = INTEGER(daily_mat);
//
//   SEXP case_mat = PROTECT(Rf_allocVector(INTSXP, stopSimulationDay));
//   int* case_mat_ptr = INTEGER(case_mat);
//
//   SEXP death_mat = PROTECT(Rf_allocVector(INTSXP, stopSimulationDay));
//   int* death_mat_ptr = INTEGER(death_mat);
//
//   SEXP cum_case = Rf_ScalarInteger(0);
//   int* cum_case_ptr = INTEGER(cum_case);
//   cum_case_ptr[0] = 0;
//
//   // calculate meanlog and sdlogs for lnorm distributions
//   double si_sdlog, si_meanlog;
//   if(lnormFlag==1){
//     si_sdlog = sqrt(log( pow(si_sd/si_mean,2) + 1. ));
//     si_meanlog = log(si_mean) - 0.5*pow(si_sdlog,2);
//   }
//
//   double symp_to_death_sdlog = sqrt(log( pow(symp_to_death_sd/symp_to_death_mean,2) + 1.));
//   double symp_to_death_meanlog = log(symp_to_death_mean) - 0.5*pow(symp_to_death_sdlog,2);
//
//   // initiate vector to store daily incidence of autochthonous infections
//   // int* dailyIncidence = malloc(stopSimulationDay*sizeof(int));
//   int* dailyIncidence = (int*)R_alloc(stopSimulationDay,sizeof(int));
//
//   // initiate vector to store daily reporting of autochthonous infections
//   int* dailyCases = (int*)R_alloc(stopSimulationDay, sizeof(int));
//
//   // initiate vector to store daily mortality of autochthonous infections
//   int* dailyMortality = (int*)R_alloc(stopSimulationDay, sizeof(int));
//
//   // zero out arrays
//   for(int i=0; i<stopSimulationDay; i++){
//     dailyIncidence[i] = 0;
//     dailyCases[i] = 0;
//     dailyMortality[i] = 0;
//   }
//
//   // initiate vector to store timing of parent infections (this will get resized a lot)
//   // double* time_exp = (double*)malloc(timeImport_len*sizeof(double));
//   double* time_exp = (double*)Calloc(timeImport_len,double);
//   Rprintf(" *** memmove(time_exp,timeImport_ptr,timeImport_len*sizeof(double)); *** \n");
//   memmove(time_exp,timeImport_ptr,timeImport_len*sizeof(double));
//   int parents = timeImport_len;
//
//   /* --------------------------------------------------------------------------------
//   #   pointers to dynamically allocated memory
//   -------------------------------------------------------------------------------- */
//
//   /* --------------------------------------------------------------------------------
//   #   simulation loop
//   -------------------------------------------------------------------------------- */
//   int sim=0;
//   while(parents > 0){
//     sim++;
//     Rprintf("\n\n\n --- SIMULATING GENERATION %d --- \n",sim);
//     Rprintf(" --- PARENTS: %d --- \n",parents);
//
//     R_CheckUserInterrupt();
//
//     // draw number of offspring from each case (parent)
//     // Rprintf(" *** calling malloc with 'number_offspring' *** \n");
//     // int* number_offspring = malloc(parents * sizeof(int));
//     Rprintf(" *** calling Calloc with 'number_offspring' *** \n");
//     int* number_offspring = (int*)Calloc(parents,int);
//     int tot_offspring = 0;
//     for(int i=0; i<parents; i++){
//       double R0;
//       double u = unif_rand();
//       if(u < asympProp){
//         R0 = R * asympRFraction;
//       } else {
//         R0 = R;
//       }
//       number_offspring[i] = (int)rnbinom_mu(k,R0);
//       tot_offspring += number_offspring[i];
//     }
//
//     Rprintf(" --- TOTAL OFFSPRING: %d --- \n",tot_offspring);
//     if(tot_offspring > 0){
//
//       /* --------------------------------------------------------------------------------
//       #   determine the time of exposure of each offspring and only for parents with offspring
//       -------------------------------------------------------------------------------- */
//
//       // Rprintf(" *** calling malloc with 'time_exp_next' *** \n");
//       // double* time_exp_next = malloc(tot_offspring * sizeof(double));
//       Rprintf(" *** calling Calloc with 'time_exp_next' *** \n");
//       double* time_exp_next = (double*)Calloc(tot_offspring,double);
//
//       int k=0;
//       for(int i=0; i<parents; i++){
//         // only for parents with offspring
//         if(number_offspring[i] > 0){
//           for(int j=0; j<number_offspring[i]; j++){
//             double offspring_time_exp = time_exp[i] + rlnorm(si_meanlog,si_sdlog);
//             // retain only those infections that occur before time is up
//             if( (int)floor(offspring_time_exp) <= stopSimulationDay ){
//               time_exp_next[k] = offspring_time_exp;
//               k++;
//             }
//           }
//         }
//       }
//
//       // new vector time_exp
//       parents = k;
//       // Rprintf(" *** calling realloc with 'time_exp' *** \n");
//       // time_exp = realloc(time_exp,parents * sizeof(double));
//       Rprintf(" *** calling Realloc with 'time_exp' *** \n");
//       time_exp = (double*)Realloc(time_exp,parents,double);
//       Rprintf(" *** memmove(time_exp,time_exp_next,parents * sizeof(double)); *** \n");
//       memmove(time_exp,time_exp_next,parents * sizeof(double));
//
//
//       /* --------------------------------------------------------------------------------
//       #   determine incubation period for each offspring
//       -------------------------------------------------------------------------------- */
//
//       // Rprintf(" *** calling malloc with 'incubation_periods' *** \n");
//       // double* incubation_periods = malloc(parents * sizeof(double));
//       Rprintf(" *** calling Calloc with 'incubation_periods' *** \n");
//       double* incubation_periods = (double*)Calloc(parents,double);
//       for(int i=0; i<parents; i++){
//         incubation_periods[i] = rweibull(inc_shape,inc_scale);
//       }
//
//
//       /* --------------------------------------------------------------------------------
//       #   determine the time of infection detection of each offspring
//       -------------------------------------------------------------------------------- */
//
//       // Rprintf(" *** calling malloc with 'time_det' *** \n");
//       // double* time_det = malloc(parents * sizeof(double));
//       Rprintf(" *** calling Calloc with 'time_det' *** \n");
//       double* time_det = (double*)Calloc(parents,double);
//
//       k=0;
//       for(int i=0; i<parents; i++){
//         double time_det_indiv = time_exp[i] + incubation_periods[i] + (double)rpois(report_delay_rate);
//         // retain only those detected infections that occur before time is up
//         if( (int)floor(time_det_indiv) <= stopSimulationDay ){
//           time_det[k] = time_det_indiv;
//           k++;
//         }
//       }
//       // time_det will have nonsense after position [0,...,time_det_k-1]
//       int time_det_k = k;
//
//
//       /* --------------------------------------------------------------------------------
//       #   determine the time of death of each offspring
//       -------------------------------------------------------------------------------- */
//
//       // Rprintf(" *** calling malloc with 'time_death' *** \n");
//       // double* time_death = malloc(parents * sizeof(double));
//       Rprintf(" *** calling Calloc with 'time_death' *** \n");
//       double* time_death = (double*)Calloc(parents,double);
//
//       k=0;
//       for(int i=0; i<parents; i++){
//         double time_death_indiv = time_exp[i] + incubation_periods[i] + rlnorm(symp_to_death_meanlog,symp_to_death_sdlog);
//         // retain only deaths that occur before time is up
//         if( (int)floor(time_death_indiv) <= stopSimulationDay ){
//           time_death[k] = time_death_indiv;
//           k++;
//         }
//       }
//       // time_death will have nonsense after position [0,...,time_death_k-1]
//       int time_death_k = k;
//
//
//       /* --------------------------------------------------------------------------------
//       #   update vector of daily incidence of detected infections
//       -------------------------------------------------------------------------------- */
//
//       // Rprintf(" *** calling malloc with 'time_exp_floor' *** \n");
//       // int* time_exp_floor = malloc(parents * sizeof(int));
//       Rprintf(" *** calling Calloc with 'time_exp_floor' *** \n");
//       int* time_exp_floor = (int*)Calloc(parents,int);
//       for(int i=0; i<parents; i++){
//         time_exp_floor[i] = (int)floor(time_exp[i]);
//       }
//
//       // Rprintf(" *** calling malloc with 'time_exp_uniq' *** \n");
//       // int* time_exp_uniq = malloc(parents * sizeof(int));
//       // Rprintf(" *** calling malloc with 'time_exp_dup' *** \n");
//       // int* time_exp_dup = malloc(parents * sizeof(int));
//
//       Rprintf(" *** calling Calloc with 'time_exp_uniq' *** \n");
//       int* time_exp_uniq = (int*)Calloc(parents,int);
//       Rprintf(" *** calling Calloc with 'time_exp_dup' *** \n");
//       int* time_exp_dup = (int*)Calloc(parents,int);
//
//       int j, n, flag;
//
//       time_exp_uniq[0] = time_exp_floor[0];
//       int count_exp = 1;
//
//       for(j=0; j <= parents-1; ++j) {
//           time_exp_dup[j] = 1;
//           flag = 1;;
//           /*the unique array will always have 'count' elements*/
//           for(n=0; n < count_exp; ++n) {
//               if(time_exp_floor[j] == time_exp_uniq[n]) {
//                   if(j != n){
//                     time_exp_dup[n] += 1;
//                   }
//                   flag = 0;
//               }
//           }
//           if(flag == 1) {
//               ++count_exp;
//               time_exp_uniq[count_exp-1] = time_exp_floor[j];
//           }
//       }
//
//       for(int j=0; j<count_exp; j++){
//         dailyIncidence[time_exp_uniq[j]] += time_exp_dup[j];
//       }
//
//
//       /* --------------------------------------------------------------------------------
//       #   update vector of daily case reporting
//       -------------------------------------------------------------------------------- */
//
//       // Rprintf(" *** calling malloc with 'time_det_floor' *** \n");
//       // int* time_det_floor = malloc(time_det_k * sizeof(int));
//       Rprintf(" *** calling Calloc with 'time_det_floor' *** \n");
//       int* time_det_floor = (int*)Calloc(time_det_k,int);
//
//       for(int i=0; i<time_det_k; i++){
//         time_det_floor[i] = (int)floor(time_det[i]);
//       }
//
//       // Rprintf(" *** calling malloc with 'time_det_uniq' *** \n");
//       // int* time_det_uniq = malloc(time_det_k * sizeof(int));
//       // Rprintf(" *** calling malloc with 'time_det_dup' *** \n");
//       // int* time_det_dup = malloc(time_det_k * sizeof(int));
//
//       Rprintf(" *** calling Calloc with 'time_det_uniq' *** \n");
//       int* time_det_uniq = (int*)Calloc(time_det_k,int);
//       Rprintf(" *** calling Calloc with 'time_det_dup' *** \n");
//       int* time_det_dup = (int*)Calloc(time_det_k,int);
//
//       time_det_uniq[0] = time_det_floor[0];
//       int count_det = 1;
//
//       for(j=0; j <= time_det_k-1; ++j) {
//           time_det_dup[j] = 1;
//           flag = 1;;
//           /*the unique array will always have 'count' elements*/
//           for(n=0; n < count_det; ++n) {
//               if(time_det_floor[j] == time_det_uniq[n]) {
//                   if(j != n){
//                     time_det_dup[n] += 1;
//                   }
//                   flag = 0;
//               }
//           }
//           if(flag == 1) {
//               ++count_det;
//               time_det_uniq[count_det-1] = time_det_floor[j];
//           }
//       }
//
//       for(int j=0; j<count_det; j++){
//         dailyCases[time_det_uniq[j]] += time_det_dup[j];
//       }
//
//
//       /* --------------------------------------------------------------------------------
//       #   update vector of daily mortality
//       -------------------------------------------------------------------------------- */
//
//       // Rprintf(" *** calling malloc with 'time_death_floor' *** \n");
//       // int* time_death_floor = malloc(time_death_k * sizeof(int));
//
//       Rprintf(" *** calling Calloc with 'time_death_floor' *** \n");
//       int* time_death_floor = (int*)Calloc(time_death_k,int);
//
//       for(int i=0; i<time_death_k; i++){
//         time_death_floor[i] = (int)floor(time_death[i]);
//       }
//
//       // Rprintf(" *** calling malloc with 'time_death_uniq' *** \n");
//       // int* time_death_uniq = malloc(time_death_k * sizeof(int));
//       // Rprintf(" *** calling malloc with 'time_death_dup' *** \n");
//       // int* time_death_dup = malloc(time_death_k * sizeof(int));
//
//       Rprintf(" *** calling Calloc with 'time_death_uniq' *** \n");
//       int* time_death_uniq = (int*)Calloc(time_death_k,int);
//       Rprintf(" *** calling Calloc with 'time_death_dup' *** \n");
//       int* time_death_dup = (int*)Calloc(time_death_k,int);
//
//       time_death_uniq[0] = time_death_floor[0];
//       int count_death = 1;
//
//       for(j=0; j <= time_death_k-1; ++j) {
//           time_death_dup[j] = 1;
//           flag = 1;;
//           /*the unique array will always have 'count' elements*/
//           for(n=0; n < count_death; ++n) {
//               if(time_death_floor[j] == time_death_uniq[n]) {
//                   if(j != n){
//                     time_death_dup[n] += 1;
//                   }
//                   flag = 0;
//               }
//           }
//           if(flag == 1) {
//               ++count_death;
//               time_death_uniq[count_death-1] = time_death_floor[j];
//           }
//       }
//
//       for(int j=0; j<count_death; j++){
//         dailyMortality[time_death_uniq[j]] += time_death_dup[j];
//       }
//
//       // free(time_exp_next);
//       // free(incubation_periods);
//       // free(time_det);
//       // free(time_death);
//       //
//       // free(time_exp_floor);
//       // free(time_exp_uniq);
//       // free(time_exp_dup);
//       //
//       // free(time_det_floor);
//       // free(time_det_uniq);
//       // free(time_det_dup);
//       //
//       // free(time_death_floor);
//       // free(time_death_uniq);
//       // free(time_death_dup);
//
//       Free(time_exp_next);
//       Free(incubation_periods);
//       Free(time_det);
//       Free(time_death);
//
//       Free(time_exp_floor);
//       Free(time_exp_uniq);
//       Free(time_exp_dup);
//
//       Free(time_det_floor);
//       Free(time_det_uniq);
//       Free(time_det_dup);
//
//       Free(time_death_floor);
//       Free(time_death_uniq);
//       Free(time_death_dup);
//
//     } // end if
//
//     // free(number_offspring);
//
//     Free(number_offspring);
//
//   } // end while
//
//   Rprintf(" --- ENDING SIMULATION --- \n");
//   // keep track of things that need to be kept track of
//
//
//   Rprintf(" *** memmove(daily_mat_ptr,dailyIncidence,stopSimulationDay * sizeof(int)); *** \n");
//   memmove(daily_mat_ptr,dailyIncidence,stopSimulationDay * sizeof(int));
//   Rprintf(" *** memmove(case_mat_ptr,dailyCases,stopSimulationDay * sizeof(int)); *** \n");
//   memmove(case_mat_ptr,dailyCases,stopSimulationDay * sizeof(int));
//   Rprintf(" *** memmove(death_mat_ptr,dailyMortality,stopSimulationDay * sizeof(int)); *** \n");
//   memmove(death_mat_ptr,dailyMortality,stopSimulationDay * sizeof(int));
//
//   for(int i=0; i<stopSimulationDay; i++){
//     cum_case_ptr[0] += dailyIncidence[i];
//   }
//
//   // Rprintf(" *** free(dailyIncidence); *** \n");
//   // free(dailyIncidence);
//   // Rprintf(" *** free(dailyCases); *** \n");
//   // free(dailyCases);
//   // Rprintf(" *** free(dailyMortality); *** \n");
//   // free(dailyMortality);
//   Rprintf(" *** free(time_exp); *** \n");
//   Free(time_exp);
//
//   // return a list to R
//   SEXP out = PROTECT(Rf_allocVector(VECSXP,4));
//   SET_VECTOR_ELT(out,0,daily_mat);
//   SET_VECTOR_ELT(out,1,death_mat);
//   SET_VECTOR_ELT(out,2,case_mat);
//   SET_VECTOR_ELT(out,3,cum_case);
//
//   SEXP names = PROTECT(Rf_allocVector(STRSXP,4));
//   SET_STRING_ELT(names,0,Rf_mkChar("daily"));
//   SET_STRING_ELT(names,1,Rf_mkChar("death"));
//   SET_STRING_ELT(names,2,Rf_mkChar("cases"));
//   SET_STRING_ELT(names,3,Rf_mkChar("cum"));
//
//   Rf_namesgets(out,names);
//
//   PutRNGstate();
//
//   UNPROTECT(5);
//   return out;
// };
