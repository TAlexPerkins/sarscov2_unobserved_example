# This is a script to generate results posted in a thread on
# Twitter by Alex Perkins (@TAlexPerkins) on March 10, 2020 regarding
# the COVID-19 pandemic and its statIndia up to March 8, 2020 in the
# United States.
#
# https://twitter.com/TAlexPerkins/statIndia/1237513140578865153
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

# load function to simulate autochthonoIndia transmission
source('simOutbreak.R')

# set random number seed
set.seed(1234)

# read in line list data for India
# updated 20200313
# data from https://github.com/midas-network/COVID-19/tree/master/data/cases/global/line_listings_nihfogarty
linelist = read.csv('../data/2020_03_12_1800EST_linelist_NIHFogarty.csv')
yesIndia = subset(linelist, country=='India')

# Only 3 cases! The news says 70 as of 3/13

# remove Diamond Princess repatriated cases
# yesIndia = yesIndia[grep("Diamond",yesIndia$summary,invert=T),]

# number of travelers that were cases or died
num.CF = c(
  nrow(subset(yesIndia,traveler>0))-sum(subset(yesIndia,traveler>0)$death>0),
  sum(subset(yesIndia,traveler>0)$death>0))

# read in case data internationally
# updated 20200313
# data from https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data/csse_covid_19_time_series
ts = read.csv('../data/time_series_19-covid-Confirmed.csv')
ts.natl = matrix(0,length(unique(ts$Country.Region)),ncol(ts)-4)
for(ii in 1:ncol(ts.natl)){
  ts.natl[,ii] = aggregate(ts[,4+ii],by=list(ts$Country.Region),FUN=sum)[,2]
}
row.names(ts.natl) = aggregate(ts[,4+ii],by=list(ts$Country.Region),FUN=sum)[,1]
for(ii in 1:nrow(ts.natl)){
  ts.natl[ii,-1] = pmax(0,diff(ts.natl[ii,]))
}

# read in death data internationally
# updated 20200313
# data from https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data/csse_covid_19_time_series
tsd = read.csv('../data/time_series_19-covid-Deaths.csv')
tsd.natl = matrix(0,length(unique(tsd$Country.Region)),ncol(tsd)-4)
for(ii in 1:ncol(ts.natl)){
  tsd.natl[,ii] = aggregate(tsd[,4+ii],by=list(tsd$Country.Region),FUN=sum)[,2]
}
row.names(tsd.natl) = aggregate(tsd[,4+ii],by=list(tsd$Country.Region),FUN=sum)[,1]
for(ii in 1:nrow(tsd.natl)){
  tsd.natl[ii,-1] = pmax(0,diff(tsd.natl[ii,]))
}
colnames(tsd.natl) = 22:(ncol(tsd.natl)+21)

# count up local cases by day in the India
cases.India.total = c(rep(0,21),ts.natl['India',])
cases.India.imported =
  table(
    as.Date(as.character(subset(yesIndia,traveler>0)$reporting.date)) -
      as.Date('2019-12-31'))
tmp = rep(0,length(cases.India.total))
tmp[as.numeric(names(cases.India.imported))] = cases.India.imported
cases.India.imported = tmp
rm(tmp)
cases.India.local = pmax(0, cases.India.total - cases.India.imported)

# count up local deaths by day in the India
deaths.India.total = c(rep(0,21),tsd.natl['India',])
deaths.India.imported =
  table(
    as.Date(as.character(subset(yesIndia,traveler>0&death>0)$reporting.date)) -
      as.Date('2019-12-31'))
tmp = rep(0,length(deaths.India.total))
tmp[as.numeric(names(deaths.India.imported))] = deaths.India.imported
deaths.India.imported = tmp
rm(tmp)
deaths.India.local = pmax(0, deaths.India.total - deaths.India.imported)

# sample replicates of how many infections have been imported into the India
maxIndia = 1e4
rangeIndia = sum(yesIndia$traveler==1):maxIndia
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
# deaths in the India as of March 8, 2020 predicted by the model
replicates = 200

# these parameter values are provisional guesses
# with these values, the model does a decent job predicting deaths in the India
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
        0:(maxIndia-sum(num.CF)),
        num.CF[1],num.CF[2]),
      prob = c(sum(propns.ASCF[ii,1:2]),propns.ASCF[ii,3:4]))
  imports[ii] =
    sample(
      sum(num.CF):maxIndia,
      1,
      prob=PrImportedInfections,
      replace=T)
}

# draw samples of the day on which imported infections arrived
case.days=vector()
for(i in 1:length(cases.India.imported)){
  if(cases.India.imported[i]>0){
    if(length(case.days)==0){
      case.days=rep(i,cases.India.imported[i])
    }else{
      case.days=c(case.days,rep(i,cases.India.imported[i]))      
    }
  }
}
import.case.density = density(
  case.days,
  from = 1,
  to = length(cases.India.imported),
  n = length(cases.India.imported))$y
import.doy = list()
for(ii in 1:replicates){
  import.doy[[ii]] = sample(
    1:length(cases.India.imported),
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
    asympRFraction = 1, # relative infectioIndianess of asymptomatics
    lnormFlag = T # toggles whether serial interval distribution is lognormal
  )
}

# set figure margins
par(mar=c(4,5,1,1))

# plot distribution of cumulative infections
pdf('../figures/cumulative_infections_India.pdf',width=5,height=4)
hist(
  unlist(lapply(local,function(ll)ll$cum)),
  col='gray',xlab='Cumulative infections',
  ylab='Number of simulations',main='',las=1)
dev.off()

# plot distribution of deaths expected to eventually result
# from infections through the end of the simulation
pdf('../figures/eventual_deaths_India.pdf',width=5,height=4)
hist(
  unlist(lapply(local,function(ll)ll$cum)) * propns.ASCF[,4],
  col='gray',xlab='Infections so far that will result in death',
  ylab='Number of simulations',main='',las=1)
dev.off()

# plot all locally acquired infections over time
pdf('../figures/infections_daily_India.pdf',width=5,height=4)
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
pdf('../figures/symptomatic_daily_India.pdf',width=5,height=4)
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
pdf('../figures/symptomatic_detected_India.pdf',width=5,height=4)
cases.mat = t(matrix(
  unlist(lapply(local, function(x) x$cases)),
  length(local[[1]]$cases),
  replicates))
p.mat = matrix(NA,nrow(cases.mat),ncol(cases.mat))
for(ii in 1:nrow(cases.mat)){
  for(jj in 1:ncol(cases.mat)){
    if(cases.mat[ii,jj]){
      p.mat[ii,jj] =
        rbeta(1,1+cases.India.local[jj],max(1,1+cases.mat[ii,jj]-cases.India.local[jj])) /
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
pdf('../figures/deaths_daily_India.pdf',width=5,height=4)
death.mat = t(matrix(
  unlist(lapply(local, function(x) x$death)),
  length(local[[1]]$death),
  replicates))
death.mat.obs = rbinom(length(death.mat), as.vector(death.mat), propns.ASCF[,4])
death.mat = matrix(death.mat.obs, replicates, ncol(death.mat))
plot(
  as.Date('2019-12-31') + 1:ncol(death.mat),
  apply(death.mat,2,function(ii)median(ii,na.rm=T)),
  ylim=c(0,max(c(deaths.India.local,death.mat))),col=1,lwd=2,type='l',xaxs='i',yaxs='i',las=1,
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
plot(deaths.India.local[31:length(deaths.India.local)],
     type='l',lwd=2,col=2,ylim=c(0,6),xaxs='i',yaxs='i',
     xaxt='n',yaxt='n',xlab='',ylab='')
dev.off()
