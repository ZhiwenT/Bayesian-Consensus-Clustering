rm(list = ls())

require(devtools)
#install_version("RcppArmadillo", version = "0.11.2.3.1", repos = "http://cran.us.r-project.org")
library(Rcpp)
library(RcppArmadillo)
library(inline)

sourceCpp("BCC.cpp")
source("BCC function (GLM)_edited.R")
#install.packages(c("mixAK","MCMCpack","label.switching","nnet","mclust","truncdist","MCMCpack","lme4","Rmpfr","mvtnorm"), 
#    repos='http://cran.us.r-project.org', dependencies=TRUE)

#--------------------------------------------------------------#
# PBC Data - Setting Up
#--------------------------------------------------------------#
# install.packages("mixAK")
library(mixAK)
data(PBCseq);head(PBCseq)
# patients known to be alive  and without liver transplantation at 910 days of follow-up
idx <- unique(PBCseq[PBCseq$alive>910,]$id); length(idx)
dnew910 <- PBCseq[PBCseq$id %in% idx,]; length(unique(dnew910$id))
dnew910$time <-  dnew910$month - mean(dnew910$month,na.rm=TRUE)
dnew910$time2 <-  dnew910$time^2

# use only data before 910 days (2.5 years)
dnew910.before <- dnew910[dnew910$day<=910,]; length(unique(dnew910.before$id))

# create a new ID variable with values from 1 to N;
subj <- unique(dnew910.before$id)
N <- length(subj)
id.new <- NULL
for (i in 1:N) {id.new   <- c(id.new, 
                              rep(i,length(dnew910.before[dnew910.before$id==subj[i],]$id)))}
dnew910.before$id.new <- id.new

head(dnew910.before)
#-------------------------------------------------------------------------------------------------------------------#
set.seed(995)
length(unique(dnew910.before[is.na(dnew910.before$bili)==FALSE,]$id.new))
length(unique(dnew910.before[is.na(dnew910.before$platelet)==FALSE,]$id.new))
length(unique(dnew910.before[is.na(dnew910.before$spiders)==FALSE,]$id.new))

dnew910.before$lbili_scale <- scale(dnew910.before$lbili)
  
  
fit.BCC <- BCC.multi (
  mydat = list(dnew910.before$lbili_scale,dnew910.before$lplatelet,dnew910.before$spiders),
  dist = c("gaussian","gaussian","binomial"),
  id = list(dnew910.before$id.new,dnew910.before$id.new,dnew910.before$id.new),
  time = list(dnew910.before$time,dnew910.before$time,dnew910.before$time),
  formula =list(y ~ time +  (1|id),
                y ~ time +  (1|id),
                y ~ time +  (1|id)),
  num.cluster = 2,
  hyper.par  = list(delta=1,a.star=1,b.star=1,aa0=0.001, bb0=0.001, ww0=0,vv0=9, cc0=0.001, dd0=0.001,rr0=4,RR0=3),
  sigma.sq.e.common = 1, 	
  c.ga.tuning = list(1,1,1),		    # tuning parameter for MH algorithm (fixed effect parameters), each parameter corresponds to an outcome/marker		
  c.theta.tuning = list(1,1,1),		# tuning parameter for MH algorithm (random effect), each parameter corresponds to an outcome/marker
  adaptive.tuning = 0,	      		# adaptive tuning parameters, 1 - yes, 0 - no
  align.clusters=0,			# assign clusters	 
  alpha.common=0,				# 1 - common alpha, 0 - separate alphas for each outcome
  sig.var = 0,				# 1 - unstructure random effect variance, 0 - diagonal random effect variance structure
  initials= NULL,			# initial values for model parameters
  initial.cluster.membership = "mixAK", # "mixAK" or "random"
  print.info="FALSE",
  burn.in = 200, 			# number of samples discarded
  thin = 10, 				# thinning
  per = 1, 				# output information every "per" iteration 
  max.iter = 1200) 			# maximum number of iteration 
fit.BCC$summary.stat

fit.BCC$alpha.adjust


















