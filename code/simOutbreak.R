# This is a script to generate results posted in a thread on
# Twitter by Alex Perkins (@TAlexPerkins) on March 10, 2020 regarding
# the COVID-19 pandemic and its status up to March 8, 2020 in the
# United States.
#
# https://twitter.com/TAlexPerkins/status/1237513140578865153
#
# Contributors to the development of this code include:
# Alex Perkins, taperkins@nd.edu
# Sean Cavany, scavany@nd.edu
# Sean Moore, smoore15@nd.edu
# Anita Lerch, alerch2@nd.edu
#
# Please visit our lab's website at http://perkinslab.weebly.com
# for more information about our research, including future updates
# on our COVID-19 work.

# This function is called by script.R and simulates local transmission.
simOutbreak = function(
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
  lnormFlag = T)
{
  # decide what to keep track of across replicate simulations
  daily.mat = matrix(0,repSims,stopSimulationDay)

  case.mat = matrix(0,repSims,stopSimulationDay)
  death.mat = matrix(0,repSims,stopSimulationDay)

  cum.vec = rep(0,repSims)


  # calculate meanlog and sdlogs for lnorm distributions
  if(lnormFlag){
    si_sdlog = sqrt(log((si_sd/si_mean)^2 + 1))
    si_meanlog = log(si_mean) - 0.5*si_sdlog^2
  }

  symp_to_death_sdlog = sqrt(log((symp_to_death_sd/symp_to_death_mean)^2 + 1))
  symp_to_death_meanlog = log(symp_to_death_mean) - 0.5*symp_to_death_sdlog^2


  # replicate repSims number of times
  for(rr in 1:repSims){
    # initiate vector to store daily incidence of autochthonous infections
    dailyIncidence = rep(0,stopSimulationDay)

    # initiate vector to store daily reporting of autochthonous infections
    dailyCases = rep(0,stopSimulationDay)

    # initiate vector to store daily mortality of autochthonous infections
    dailyMortality = rep(0,stopSimulationDay)

    # initiate vector to store timing of parent infections
    time.exp = timeImport

    # loop through generations of transmission until extinct or time is up
    while(length(time.exp) > 0){

      # draw a number of offspring for each parent
      R_vec = ifelse(runif(length(time.exp),0,1) < asympProp, R*asympRFraction, R)
      number.offspring = rnbinom(n=length(time.exp), mu=R_vec, size=k)

      if(sum(number.offspring > 0)){

        # remove parents with zero offspring
        time.exp = time.exp[number.offspring>0]
        number.offspring = number.offspring[number.offspring>0]

        # determine the time of exposure of each offspring
        if(!lnormFlag) {
          time.exp =
            rep(time.exp,times=number.offspring) +
            rnorm(sum(number.offspring), mean=si_mean, sd=si_sd)
        } else {
          time.exp =
            rep(time.exp,times=number.offspring) +
            rlnorm(sum(number.offspring), meanlog=si_meanlog, sdlog=si_sdlog)
        }

        # retain only those infections that occur before time is up
        time.exp =
          time.exp[floor(time.exp) <= stopSimulationDay]

        # determine incubation period for each offspring
        incubation_periods = rweibull(length(time.exp), shape=inc_shape, scale=inc_scale)

        # determine the time of infection detection of each offspring
        time.det =
          time.exp +
          incubation_periods +
          rpois(length(time.exp), lambda=report_delay_rate)

        # retain only those detected infections that occur before time is up
        time.det =
          time.det[floor(time.det) <= stopSimulationDay]

        # update vector of daily incidence of detected infections
        dailyIncidence[as.numeric(names(table(floor(time.exp))))] =
          dailyIncidence[as.numeric(names(table(floor(time.exp))))] +
          table(floor(time.exp))

        # update vector of daily case reporting
        dailyCases[as.numeric(names(table(floor(time.det))))] =
          dailyCases[as.numeric(names(table(floor(time.det))))] +
          table(floor(time.det))

        # determine the time of death of each offspring
        time.death =
          time.exp +
          incubation_periods +
          rlnorm(length(time.exp),meanlog=symp_to_death_meanlog,
                 sdlog=symp_to_death_sdlog)

        # retain only deaths that occur before time is up
        time.death =
          time.death[floor(time.death) <= stopSimulationDay]

        # update vector of daily mortality
        dailyMortality[as.numeric(names(table(floor(time.death))))] =
          dailyMortality[as.numeric(names(table(floor(time.death))))] +
            table(floor(time.death))
      }
    }
    # keep track of things that need to be kept track of
    daily.mat[rr,] = dailyIncidence
    case.mat[rr,] = dailyCases
    death.mat[rr,] = dailyMortality
    cum.vec[rr] = sum(dailyIncidence)
  }

  # return things that need to be returned
  return(list(daily=daily.mat,death=death.mat,cases=case.mat,cum=cum.vec))
}
