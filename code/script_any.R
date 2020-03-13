# This is a script to generate results posted in a thread on
# Twitter by Alex Perkins (@TAlexPerkins) on March 10, 2020 regarding
# the COVID-19 pandemic and its statCountry up to March 8, 2020 in the
# United States.
#
# https://twitter.com/TAlexPerkins/statUS/1237513140578865153
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

# load libraries
library(extraDistr)
library(doParallel)
library(mc2d)

# load function to simulate autochthonoCountry transmission
source('simOutbreak.R')

# This code has assumptions that might not be true for every country.
# This is just a quick check
run_sims_DANGER <- function(cntry, cntry2) {
  
  # updated 20200313
  datasrc <- list(
    # https://github.com/midas-network/COVID-19/tree/master/data/cases/global/line_listings_nihfogarty
    # linelist = '../data/2020_03_12_1800EST_linelist_NIHFogarty.csv',
    linelist = '../data/2020_03_07_1800EST_linelist_NIHFogarty.csv',
    
    # https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data/csse_covid_19_time_series
    confirmed = '../data/time_series_19-covid-Confirmed.csv', 

    # https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data/csse_covid_19_time_series
    deaths = '../data/time_series_19-covid-Deaths.csv'
  )

  # set random number seed
  set.seed(1234)
  
  # read in line list data for CNTRY
  # data from https://github.com/midas-network/COVID-19/tree/master/data/cases/global/line_listings_nihfogarty
  linelist = read.csv(datasrc$linelist)
  yesCNTRY = subset(linelist, country==cntry)
  
  # remove Diamond Princess repatriated cases
  yesCNTRY = yesCNTRY[grep("Diamond",yesCNTRY$summary,invert=T),]
  
  # number of travelers that were cases or died
  num.CF = c(
    nrow(subset(yesCNTRY,traveler>0))-sum(subset(yesCNTRY,traveler>0)$death>0),
    sum(subset(yesCNTRY,traveler>0)$death>0))
  
  # read in case data internationally
  # data from https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data/csse_covid_19_time_series
  ts = read.csv(datasrc$confirmed)
  ts.natl = matrix(0,length(unique(ts$Country.Region)),ncol(ts)-4)
  for(ii in 1:ncol(ts.natl)){
    ts.natl[,ii] = aggregate(ts[,4+ii],by=list(ts$Country.Region),FUN=sum)[,2]
  }
  row.names(ts.natl) = aggregate(ts[,4+ii],by=list(ts$Country.Region),FUN=sum)[,1]
  for(ii in 1:nrow(ts.natl)){
    ts.natl[ii,-1] = pmax(0,diff(ts.natl[ii,]))
  }
  
  # read in death data internationally
  # updated 20200307
  # data from https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data/csse_covid_19_time_series
  tsd = read.csv(datasrc$deaths)
  tsd.natl = matrix(0,length(unique(tsd$Country.Region)),ncol(tsd)-4)
  for(ii in 1:ncol(ts.natl)){
    tsd.natl[,ii] = aggregate(tsd[,4+ii],by=list(tsd$Country.Region),FUN=sum)[,2]
  }
  row.names(tsd.natl) = aggregate(tsd[,4+ii],by=list(tsd$Country.Region),FUN=sum)[,1]
  for(ii in 1:nrow(tsd.natl)){
    tsd.natl[ii,-1] = pmax(0,diff(tsd.natl[ii,]))
  }
  colnames(tsd.natl) = 22:(ncol(tsd.natl)+21)
  
  # count up local cases by day in the CNTRY
  cases.CNTRY.total = c(rep(0,21),ts.natl[cntry2,])
  cases.CNTRY.imported =
    table(
      as.Date(as.character(subset(yesCNTRY,traveler>0)$reporting.date)) -
        as.Date('2019-12-31'))
  tmp = rep(0,length(cases.CNTRY.total))
  tmp[as.numeric(names(cases.CNTRY.imported))] = cases.CNTRY.imported
  cases.CNTRY.imported = tmp
  rm(tmp)
  cases.CNTRY.local = pmax(0, cases.CNTRY.total - cases.CNTRY.imported)
  
  # count up local deaths by day in the CNTRY
  deaths.CNTRY.total = c(rep(0,21),tsd.natl[cntry2,])
  deaths.CNTRY.imported =
    table(
      as.Date(as.character(subset(yesCNTRY,traveler>0&death>0)$reporting.date)) -
        as.Date('2019-12-31'))
  tmp = rep(0,length(deaths.CNTRY.total))
  tmp[as.numeric(names(deaths.CNTRY.imported))] = deaths.CNTRY.imported
  deaths.CNTRY.imported = tmp
  rm(tmp)
  deaths.CNTRY.local = pmax(0, deaths.CNTRY.total - deaths.CNTRY.imported)
  
  # sample replicates of how many infections have been imported into the CNTRY
  maxCNTRY = 1e4
  # Alex, is this na.rm fine?
  rangeCNTRY = sum(yesCNTRY$traveler==1, na.rm=TRUE):maxCNTRY
  # estimate for asymptomatic proportion based on
  # https://www.medrxiv.org/content/10.1101/2020.02.20.20025866v2
  PrAsymptomatic = exp(optim(par=c(0,0),fn=function(par){
    sum((
      qbeta(c(0.5,0.025,0.975),exp(par[1]),exp(par[2])) -
        c(0.179,0.155,0.202)) ^ 2)})$par)
  # estimate for proportion of symptomatic infections resulting in death based on
  # http://weekly.chinacdc.cn/en/article/id/e53946e2-c6c4-41e9-9a9b-fea8db1a8f51
  PrDeathSymptom = c(1+1023,1+44672-1023)
  
  # set values of unknown parameters
  # note that these values seem to maximize the probability of the cumulative
  # deaths in the CNTRY as of March 8, 2020 predicted by the model
  replicates = 200
  
  # these parameter values are provisional guesses
  # with these values, the model does a decent job predicting deaths in the CNTRY
  PrCaseSymptom.trav = rep(0.38,replicates)
  asympRFraction = rep(1,replicates)
  
  # sample from uncertainty about proportions of infection outcomes
  propns.ASCF = cbind(
    rbeta(replicates,PrAsymptomatic[1],PrAsymptomatic[2]),
    rbeta(replicates,PrDeathSymptom[1],PrDeathSymptom[2]))
  propns.ASCF = cbind(
    propns.ASCF[,1],
    (1-propns.ASCF[,1]) * (1-PrCaseSymptom.trav) * (1-propns.ASCF[,2]),
    (1-propns.ASCF[,1]) * PrCaseSymptom.trav * (1-propns.ASCF[,2]),
    (1-propns.ASCF[,1]) * propns.ASCF[,2])
  
  # draw samples of the number of imported infections
  imports = numeric(length=replicates)
  for(ii in 1:replicates){
    PrImportedInfections =
      dmultinomial(
        x = cbind(
          0:(maxCNTRY-sum(num.CF)),
          num.CF[1],num.CF[2]),
        prob = c(sum(propns.ASCF[ii,1:2]),propns.ASCF[ii,3:4]))
    imports[ii] =
      sample(
        sum(num.CF):maxCNTRY,
        1,
        prob=PrImportedInfections,
        replace=T)
  }
  
  # draw samples of the day on which imported infections arrived
  case.days=vector()
  for(i in 1:length(cases.CNTRY.imported)){
    if(cases.CNTRY.imported[i]>0){
      if(length(case.days)==0){
        case.days=rep(i,cases.CNTRY.imported[i])
      }else{
        case.days=c(case.days,rep(i,cases.CNTRY.imported[i]))      
      }
    }
  }
  import.case.density = density(
    case.days,
    from = 1,
    to = length(cases.CNTRY.imported),
    n = length(cases.CNTRY.imported))$y
  import.doy = list()
  for(ii in 1:replicates){
    import.doy[[ii]] = sample(
      1:length(cases.CNTRY.imported),
      imports[ii],
      prob=import.case.density,
      replace=T)
  }
  
  # simulate local transmission for each draw of imported infections
  local = foreach(ii = 1:replicates) %do% {
    simOutbreak(
      timeImport = import.doy[[ii]], # timing of each imported infection
      R = 1.97, # reproduction number
      k = 1e3, # dispersion parameter
      si_mean = 4.56, # mean of serial interval distribution
      si_sd = 0.95, # standard deviation of serial interval distribution
      inc_shape = 1.88, # shape parameter of incubation period distribution
      inc_scale = 7.97, # scale parameter of incubation period distribution
      symp_to_death_mean = 22.3, # mean of time between symptom onset and death
      symp_to_death_sd = 0.42 * 22.3, # std. dev. of time between symptom onset and death
      report_delay_rate = 3, # mean delay between symptoms and reporting
      stopSimulationDay = 68, # day of year since Jan 1 when simulation stops
      repSims = 1, # number of replicates to produce for these parameters
      asympProp = propns.ASCF[ii,1], # proportion of infections that are asymptomatic
      asympRFraction = 1, # relative infectiousness of asymptomatics
      lnormFlag = T # toggles whether serial interval distribution is lognormal
    )
  }
  
  # set figure margins
  par(mar=c(4,5,1,1))
  
  # plot distribution of cumulative infections
  pdf(paste0('../figures/cumulative_infections_', cntry, '.pdf'),width=5,height=4)
  hist(
    unlist(lapply(local,function(ll)ll$cum)),
    col='gray',xlab='Cumulative infections',
    ylab='Number of simulations',main='',las=1)
  dev.off()
  
  # plot distribution of deaths expected to eventually result
  # from infections through the end of the simulation
  pdf(paste0('../figures/eventual_deaths_', cntry, '.pdf'), width=5,height=4)
  hist(
    unlist(lapply(local,function(ll)ll$cum)) * propns.ASCF[,4],
    col='gray',xlab='Infections so far that will result in death',
    ylab='Number of simulations',main='',las=1)
  dev.off()
  
  # plot all locally acquired infections over time
  pdf(paste0('../figures/infections_daily_', cntry, '.pdf'),width=5,height=4)
  local.mat = t(matrix(
    unlist(lapply(local, function(x) x$daily)),
    length(local[[1]]$daily),
    replicates))
  plot(
    as.Date('2019-12-31') + 1:ncol(local.mat),
    apply(local.mat,2,function(ii)median(ii,na.rm=T)),
    ylim=c(0,max(local.mat)),col=1,lwd=2,type='l',xaxs='i',yaxs='i',las=1,
    xlim=as.Date('2019-12-31') + c(31,ncol(local.mat)),
    xlab='Date',ylab='Daily infections',main='')
  polygon(
    c(as.Date('2019-12-31') + 1:ncol(local.mat),
      rev(as.Date('2019-12-31') + 1:ncol(local.mat))),
    c(apply(local.mat,2,function(ii)quantile(ii,0.025,na.rm=T)),
      rev(apply(local.mat,2,function(ii)quantile(ii,0.975,na.rm=T)))),
    border=NA,col=rgb(0,0,0,0.25))
  dev.off()
  
  # plot locally acquired symptomatic infections over time
  pdf(paste0('../figures/symptomatic_daily_', cntry, '.pdf'),width=5,height=4)
  cases.mat = t(matrix(
    unlist(lapply(local, function(x) x$cases)),
    length(local[[1]]$cases),
    replicates))
  cases.mat.obs = rbinom(length(cases.mat), as.vector(cases.mat), rowSums(propns.ASCF[,2:3]))
  cases.mat = matrix(cases.mat.obs, replicates, ncol(cases.mat))
  plot(
    as.Date('2019-12-31') + 1:ncol(cases.mat),
    apply(cases.mat,2,function(ii)median(ii,na.rm=T)),
    ylim=c(0,max(cases.mat)),col=1,lwd=2,type='l',xaxs='i',yaxs='i',las=1,
    xlim=as.Date('2019-12-31') + c(31,ncol(cases.mat)),
    xlab='Date',ylab='Daily symptomatic infections',main='')
  polygon(
    c(as.Date('2019-12-31') + 1:ncol(cases.mat),
      rev(as.Date('2019-12-31') + 1:ncol(cases.mat))),
    c(apply(cases.mat,2,function(ii)quantile(ii,0.025,na.rm=T)),
      rev(apply(cases.mat,2,function(ii)quantile(ii,0.975,na.rm=T)))),
    border=NA,col=rgb(0,0,0,0.25))
  dev.off()
  
  # plot proportion of locally acquired symptomatic infections reported over time
  pdf(paste0('../figures/symptomatic_detected_', cntry, '.pdf'),width=5,height=4)
  cases.mat = t(matrix(
    unlist(lapply(local, function(x) x$cases)),
    length(local[[1]]$cases),
    replicates))
  p.mat = matrix(NA,nrow(cases.mat),ncol(cases.mat))
  for(ii in 1:nrow(cases.mat)){
    for(jj in 1:ncol(cases.mat)){
      if(cases.mat[ii,jj]){
        p.mat[ii,jj] =
          rbeta(1,1+cases.CNTRY.local[jj],max(1,1+cases.mat[ii,jj]-cases.CNTRY.local[jj])) /
          sum(propns.ASCF[2:3])
      }
    }
  }
  plot(
    as.Date('2019-12-31') + 1:ncol(p.mat),
    apply(p.mat,2,function(ii)median(ii,na.rm=T)),
    ylim=c(0,1),col=1,lwd=2,type='l',xaxs='i',yaxs='i',las=1,
    xlim=as.Date('2019-12-31') + c(31,ncol(p.mat)),
    xlab='Date',ylab='Symptomatics reporting',
    main='')
  polygon(
    c(as.Date('2019-12-31') + 1:ncol(p.mat),
      rev(as.Date('2019-12-31') + 1:ncol(p.mat))),
    c(apply(p.mat,2,function(ii)quantile(ii,0.025,na.rm=T)),
      rev(apply(p.mat,2,function(ii)quantile(ii,0.975,na.rm=T)))),
    border=NA,col=rgb(0,0,0,0.25))
  dev.off()
  
  # plot deaths resulting from locally acquired infections over time
  # note that this does not show deaths that are assumed to happen in
  # the future as a result of recent infection 
  pdf(paste0('../figures/deaths_daily_', cntry, '.pdf'),width=5,height=4)
  death.mat = t(matrix(
    unlist(lapply(local, function(x) x$death)),
    length(local[[1]]$death),
    replicates))
  death.mat.obs = rbinom(length(death.mat), as.vector(death.mat), propns.ASCF[,4])
  death.mat = matrix(death.mat.obs, replicates, ncol(death.mat))
  plot(
    as.Date('2019-12-31') + 1:ncol(death.mat),
    apply(death.mat,2,function(ii)median(ii,na.rm=T)),
    ylim=c(0,max(c(deaths.CNTRY.local,death.mat))),col=1,lwd=2,type='l',xaxs='i',yaxs='i',las=1,
    xlim=as.Date('2019-12-31') + c(31,ncol(death.mat)),
    xlab='Date',ylab='Daily deaths',main='')
  polygon(
    c(as.Date('2019-12-31') + 1:ncol(death.mat),
      rev(as.Date('2019-12-31') + 1:ncol(death.mat))),
    c(apply(death.mat,2,function(ii)quantile(ii,0.025,na.rm=T)),
      rev(apply(death.mat,2,function(ii)quantile(ii,0.975,na.rm=T)))),
    border=NA,col=rgb(0,0,0,0.25))
  legend("topleft",lty=rep("solid",2),lwd=2,
         legend=c("Data", "Model"),col=c("red","black"),
         bty='n') 
  par(new=T)
  plot(deaths.CNTRY.local[31:length(deaths.CNTRY.local)],
       type='l',lwd=2,col=2,ylim=c(0,6),xaxs='i',yaxs='i',
       xaxt='n',yaxt='n',xlab='',ylab='')
  dev.off()
  
}




# main --------------------------------------------------------------------
run_sims_DANGER('USA', 'US')
# run_sims_DANGER('France', 'France')
# run_sims_DANGER('Italy', 'Italy')
# run_sims_DANGER('India', 'India')

