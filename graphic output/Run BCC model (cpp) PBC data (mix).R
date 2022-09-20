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
dnew910$time <-  dnew910$month
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


#----------------------------------------------------------------------------------------------#
# Individual Trajectories for each marker
#----------------------------------------------------------------------------------------------#
library(ggplot2)
gp1 <- ggplot(data = dnew910.before, aes(x =month, y =lbili))+
  geom_point(size=2,alpha=0.5) +
  geom_line(aes(x = month, y = lbili,group=id ),size=1.5,alpha=0.5)+ 
  geom_smooth(method = "loess", size = 3,se = FALSE,span=2)+
  theme_bw() + 
  theme(legend.position =c(0.15,0.88),legend.title=element_blank(),
        plot.title = element_text(size = 12, face = "bold"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text=element_text(size=10),
        axis.text.x = element_text(angle = 0 ),
        strip.text.x = element_text(size = 16, angle = 0),
        strip.text.y = element_text(size = 16,face="bold")) +
  xlab("months") + ylab("lbili")
gp2 <- ggplot(data = dnew910.before, aes(x =month, y =lsgot))+
  geom_point(size=2,alpha=0.5) +
  geom_line(aes(x = month, y = lsgot,group=id ),size=1.5,alpha=0.5)+ 
  geom_smooth(method = "loess", size = 3,se = FALSE,span=2)+
  theme_bw() + 
  theme(legend.position =c(0.15,0.88),legend.title=element_blank(),
        plot.title = element_text(size = 12, face = "bold"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text=element_text(size=10),
        axis.text.x = element_text(angle = 0 ),
        strip.text.x = element_text(size = 16, angle = 0),
        strip.text.y = element_text(size = 16,face="bold")) +
  xlab("months") + ylab("lsgot")
gp3 <- ggplot(data = dnew910.before, aes(x =month, y = spiders))+
  geom_point(size=3,alpha=0.5) +
  geom_line(aes(x = month, y = spiders,group=id ),size=1.5,alpha=0.5)+ 
  geom_smooth(method = "loess", size = 3,se = FALSE,span=2)+
  theme_bw() + 
  theme(legend.position =c(0.15,0.88),legend.title=element_blank(),
        plot.title = element_text(size = 12, face = "bold"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text=element_text(size=10),
        axis.text.x = element_text(angle = 0 ),
        strip.text.x = element_text(size = 16, angle = 0),
        strip.text.y = element_text(size = 16,face="bold")) +
  xlab("months") + ylab("spiders")

dev.new(width=180, height=90)
library(cowplot)
plot_grid(gp1,gp2, gp3,labels=c("(A)", "(B)", "(C)"), nrow = 1,   align = "v" )


#-------------------------------------------------------------------------------------------------------------------#
# fit the BCC model with 2 clusters
#-------------------------------------------------------------------------------------------------------------------#
set.seed(2002)
fit.BCC <- BCC.multi (
  mydat = list(dnew910.before$lbili,dnew910.before$lsgot,dnew910.before$spiders),
  dist = c("gaussian","gaussian","binomial"),
  id = list(dnew910.before$id.new,dnew910.before$id.new,dnew910.before$id.new),
  time = list(dnew910.before$time,dnew910.before$time,dnew910.before$time),
  formula =list(y ~ time +  (1|id),
                y ~ time +  (1|id),
                y ~ time +  (1|id)),
  num.cluster = 2,
  hyper.par  = list(delta=1,a.star=1,b.star=1,aa0=0.001, bb0=0.001, ww0=0,vv0=25, cc0=0.001, dd0=0.001,rr0=4,RR0=3),
  sigma.sq.e.common = 1, 	
  c.ga.tuning = list(1,1,1),		    # tuning parameter for MH algorithm (fixed effect parameters), each parameter corresponds to an outcome/marker		
  c.theta.tuning = list(1,1,1),		# tuning parameter for MH algorithm (random effect), each parameter corresponds to an outcome/marker
  adaptive.tuning = 0,	      		# adaptive tuning parameters, 1 - yes, 0 - no
  align.clusters=0,			# assign clusters	 
  alpha.common=0,				# 1 - common alpha, 0 - separate alphas for each outcome
  sig.var = 0,				# 1 - unstructured random effect variance, 0 - diagonal random effect variance structure
  initials= NULL,			# initial values for model parameters
  initial.cluster.membership = "mixAK", # "mixAK" or "random"
  print.info="FALSE",
  burn.in = 2000, 			# number of samples discarded
  thin = 10, 				# thinning
  per = 100, 				# output information every "per" iteration 
  max.iter = 12000) 			# maximum number of iteration 

# print the run time in seconds
as.numeric(fit.BCC$run.time) 

# print the summary statistics for all parameters 
fit.BCC$summary.stat

# print the proportion \pi for each cluster (mean, sd, 2.5% and 97.5% percentile)
# geweke statistics (geweke.stat) between -2 and 2 suggests the parameters converge
fit.BCC$summary.stat$PPI

# making a trace plot for PPI to inspect convergence 
num.sample <- (fit.BCC$max.iter -  fit.BCC$burn.in)/fit.BCC$thin # number of MCMC sample
dev.new(width=180, height=90)
par(mfrow=c(1,2))
plot(1:num.sample,fit.BCC$PPI[,1],type="l",col="blue" ,
     xlab="MCMC iterations",ylab=expression(pi[1]),main= expression(paste("Trace of ",pi[1],sep="")),
     cex.lab=1.6,cex.axis=1.6,cex.main=2)
plot(1:num.sample,fit.BCC$PPI[,2],type="l",col="blue" ,
     xlab="MCMC iterations",ylab=expression(pi[2]),main= expression(paste("Trace of ",pi[2],sep="")),
     cex.lab=1.6,cex.axis=1.6,cex.main=2)

# print the adherence parameters alpha for each marker
# geweke statistics (geweke.stat) between -2 and 2 suggests the parameters converge
fit.BCC$summary.stat$ALPHA

# making a trace plot for alpha to inspect convergence 
dev.new(width=180, height=90)
par(mfrow=c(1,3))
plot(1:num.sample,fit.BCC$ALPHA[,1],type="l",col="blue" ,
     xlab="MCMC iterations",ylab=expression(alpha[1]),main= expression(paste("Trace of ",alpha[1],sep="")),
     cex.lab=1.6,cex.axis=1.6,cex.main=2)
plot(1:num.sample,fit.BCC$ALPHA[,2],type="l",col="blue" ,
     xlab="MCMC iterations",ylab=expression(alpha[2]),main= expression(paste("Trace of ",alpha[2],sep="")),
     cex.lab=1.6,cex.axis=1.6,cex.main=2)
plot(1:num.sample,fit.BCC$ALPHA[,3],type="l",col="blue" ,
     xlab="MCMC iterations",ylab=expression(alpha[2]),main= expression(paste("Trace of ",alpha[3],sep="")),
     cex.lab=1.6,cex.axis=1.6,cex.main=2)

# print the fixed effect coefficients \gamma for the first marker (lbili) (with 2 clusters)
fit.BCC$summary.stat$GA[[1]]

# print the fixed effect coefficients \gamma for the second marker (lsgot) (with 2 clusters)
fit.BCC$summary.stat$GA[[2]]

# print the fixed effect coefficients \gamma for the second marker (spider) (with 2 clusters)
fit.BCC$summary.stat$GA[[3]]


#-------------------------------------------------------------------------------------------------------------------#
# Making plot for each marker by clusters
#-------------------------------------------------------------------------------------------------------------------#
library(ggplot2)
dat <- fit.BCC$dat
dat[[1]]$time.org <- dat[[1]]$time + mean(dnew910$month,na.rm=TRUE) # transform the time back to the orignal scale

# For the first marker by local clustering
dat[[1]]$cluster.local <- factor(dat[[1]]$cluster.local,label=c("Cluster 1","Cluster 2"))
gp1 <- ggplot(data = dat[[1]], aes(x =time.org, y =y, color=cluster.local,linetype=cluster.local,fill=cluster.local))+
  geom_point(size=2,alpha=0.2) +
  geom_line(aes(x = time.org, y = y,group=id,color=cluster.local),size=1.5,alpha=0.2)+ 
  geom_smooth(method = "loess", size = 3,se = FALSE,span=2) +
  ggtitle(expression(paste("lbili (",hat(alpha), " = 0.98 )", sep=""))) + 
  theme_bw() + 
  ylab("lbili") + xlab("months")+
  theme(legend.position =c(0.8,0.1),legend.title=element_blank(),
        plot.title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text=element_text(size=10),
        axis.text.x = element_text(angle = 0 ),
        strip.text.x = element_text(size = 16, angle = 0),
        strip.text.y = element_text(size = 16,face="bold")) +
  guides(color=guide_legend(nrow=2,byrow=TRUE), linetype=guide_legend(nrow=2,byrow=TRUE),fill=guide_legend(nrow=2,byrow=TRUE)) 


# For the second marker by local clustering
dat[[2]]$cluster.local <- factor(dat[[2]]$cluster.local,label=c("Cluster 1","Cluster 2"))
dat[[2]]$time.org <- dat[[2]]$time + mean(dnew910$month,na.rm=TRUE) # transform the time back to the original scale

gp2 <- ggplot(data = dat[[2]], aes(x =time.org, y =y, color=cluster.local,linetype=cluster.local,fill=cluster.local))+
  geom_point(size=2,alpha=0.2) +
  geom_line(aes(x = time.org, y = y,group=id,color=cluster.local),size=1.5,alpha=0.2)+ 
  geom_smooth(method = "loess", size = 3,se = FALSE,span=2) +
  ggtitle(expression(paste("lsgot (",hat(alpha), " = 0.79 )", sep=""))) + 
  theme_bw() + 
  ylab("lsgot") + xlab("months")+
  theme(legend.position =c(0.8,0.1),legend.title=element_blank(),
        plot.title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text=element_text(size=10),
        axis.text.x = element_text(angle = 0 ),
        strip.text.x = element_text(size = 16, angle = 0),
        strip.text.y = element_text(size = 16,face="bold")) +
  guides(color=guide_legend(nrow=2,byrow=TRUE), linetype=guide_legend(nrow=2,byrow=TRUE),fill=guide_legend(nrow=2,byrow=TRUE)) 

# For the third marker by local clustering
dat[[3]]$cluster.local <- factor(dat[[3]]$cluster.local,label=c("Cluster 1","Cluster 2"))
dat[[3]]$time.org <- dat[[3]]$time + mean(dnew910$month,na.rm=TRUE) # transform the time back to the original scale

gp3 <- ggplot(data = dat[[3]], aes(x =time.org, y =y, color=cluster.local,linetype=cluster.local,fill=cluster.local))+
  geom_point(size=2,alpha=0.5) +
  geom_smooth(method = "loess", size = 3,se = FALSE,span=2) +
  geom_line(aes(x = time.org, y = y,group=id,color=cluster.local),size=1.5,alpha=0.2)+ 
  ggtitle(expression(paste("spiders (",hat(alpha), " = 0.75 )", sep=""))) + 
  theme_bw() + 
  ylab("spiders") + xlab("months")+
  theme(legend.position =c(0.9,0.4),legend.title=element_blank(),
        plot.title = element_text(size = 16, face = "bold"),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text=element_text(size=10),
        axis.text.x = element_text(angle = 0 ),
        strip.text.x = element_text(size = 16, angle = 0),
        strip.text.y = element_text(size = 16,face="bold")) +
  guides(color=guide_legend(nrow=2,byrow=TRUE), linetype=guide_legend(nrow=2,byrow=TRUE),fill=guide_legend(nrow=2,byrow=TRUE)) 


# For the first marker by global(consensus) clustering
dat[[1]]$cluster.global <- factor(dat[[1]]$cluster.global)
gp4 <- ggplot(data = dat[[1]], aes(x =time.org, y =y, color=cluster.global,linetype=cluster.global,fill=cluster.global))+
  geom_point(size=2,alpha=0.2) +
  geom_line(aes(x = time.org, y = y,group=id,color=cluster.global),size=1.5,alpha=0.2)+ 
  geom_smooth(method = "loess", size = 3,se = FALSE,span=2) +
  ggtitle("lbili") + 
  theme_bw() + 
  ylab("lbili") + xlab("months")+
  theme(legend.position = "none",legend.title=element_blank(),
        plot.title = element_text(size = 16),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text=element_text(size=10),
        axis.text.x = element_text(angle = 0 ),
        strip.text.x = element_text(size = 16, angle = 0),
        strip.text.y = element_text(size = 16,face="bold")) +
  guides(color=guide_legend(nrow=2,byrow=TRUE), linetype=guide_legend(nrow=2,byrow=TRUE),fill=guide_legend(nrow=2,byrow=TRUE)) 


# For the second marker by global(consensus) clustering
dat[[2]]$cluster.global <- factor(dat[[2]]$cluster.global)
gp5 <- ggplot(data = dat[[2]], aes(x =time.org, y = y, color=cluster.global,linetype=cluster.global,fill=cluster.global))+
  geom_point(size=2,alpha=0.2) +
  geom_line(aes(x = time.org, y = y,group=id,color=cluster.global),size=1.5,alpha=0.2)+ 
  geom_smooth(method = "loess", size = 3,se = FALSE,span=2)+
  ggtitle("lsgot") + 
  theme_bw() + 
  ylab("lsgot") + xlab("months")+
  theme(legend.position ="none",legend.title=element_blank(),
        plot.title = element_text(size = 16 ),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text=element_text(size=10),
        axis.text.x = element_text(angle = 0 ),
        strip.text.x = element_text(size = 16, angle = 0),
        strip.text.y = element_text(size = 16))  

# For the third marker by global(consensus) clustering
dat[[3]]$cluster.global <- factor(dat[[3]]$cluster.global)
gp6 <- ggplot(data = dat[[3]], aes(x =time.org, y = y, color=cluster.global,linetype=cluster.global,fill=cluster.global))+
  geom_point(size=2,alpha=0.2) +
  geom_line(aes(x = time.org, y = y,group=id,color=cluster.global),size=1.5,alpha=0.2)+ 
  geom_smooth(method = "loess", size = 3,se = FALSE,span=2)+
  ggtitle("spiders") + 
  theme_bw() + 
  ylab("spiders") + xlab("months")+
  theme(legend.position ="none",legend.title=element_blank(),
        plot.title = element_text(size = 16),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text=element_text(size=10),
        axis.text.x = element_text(angle = 0 ),
        strip.text.x = element_text(size = 16, angle = 0),
        strip.text.y = element_text(size = 16,face="bold")) 

dev.new(width=180, height=90)
library(cowplot)
plot_grid(gp1,gp2,gp3, gp4,gp5,gp6, labels=c("(A)", "(B)", "(C)", "(D)", "(E)", "(F)" ), nrow = 2,   align = "v" )


# Making Kaplan–Meier (KM) for survival probability after 2.5 years
library(survminer)
library(survival)

# create new ID from 1 to N - for dnew910;
subj <- unique(dnew910$id)
N <- length(subj)
id.new <- NULL
for (i in 1:N) {id.new   <- c(id.new, rep(i,length(dnew910[dnew910$id==subj[i],]$id)))}
dnew910$id.new <- id.new
dnew910_uq <- dnew910[!duplicated(dnew910$id, fromLast=TRUE),] # Keep last observation per ID


myplot <- list()
# KM Curve by local clustering for the first marker
dat[[1]]$id.new <- dat[[1]]$id
dat1 <- data.frame(id.new= dat[[1]]$id,cluster.local1 = factor(dat[[1]]$cluster.local,label=c(1,2)))
dat1_uq <- dat1[!duplicated(dat1$id.new, fromLast=TRUE),]; 
dnew910_uq <- merge(dat1_uq,dnew910_uq,by="id.new")
dnew910_uq$cluster <- dnew910_uq$cluster.local1

fit <- survfit(Surv(month, delta.death) ~  cluster,data = dnew910_uq, start.time=30.08)
myplot[[1]] <- ggsurvplot(fit, data = dnew910_uq, risk.table = TRUE,pval=TRUE,conf.int = TRUE,
                          pval.coord=c(40,0.2), xlab="months",title = expression(paste("lbili (",hat(alpha), " = 0.98 )", sep=""))) + 
  guides(colour = guide_legend(nrow = 2))

# KM Curve by local clustering for the second marker
dat[[2]]$id.new <- dat[[2]]$id
dat2 <- data.frame(id.new= dat[[2]]$id,cluster.local2 = factor(dat[[2]]$cluster.local,label=c(1,2)))
dat2_uq <- dat2[!duplicated(dat2$id.new, fromLast=TRUE),]; 
dnew910_uq <- merge(dat2_uq,dnew910_uq,by="id.new")
dnew910_uq$cluster <- dnew910_uq$cluster.local2

fit <- survfit(Surv(month, delta.death) ~  cluster,data = dnew910_uq,start.time=30.08)
myplot[[2]] <- ggsurvplot(fit, data = dnew910_uq, risk.table = TRUE,pval=TRUE,conf.int = TRUE,
                          pval.coord=c(40,0.2), xlab="months",title=expression(paste("lsgot (",hat(alpha), " = 0.79)", sep=""))) + 
  guides(colour = guide_legend(nrow = 2))


#-------------------------------------------------------------------#
dat[[3]]$id.new <- dat[[3]]$id
dat3 <- data.frame(id.new= dat[[3]]$id,cluster.local3 = factor(dat[[3]]$cluster.local,label=c(1,2)))
dat3_uq <- dat3[!duplicated(dat3$id.new, fromLast=TRUE),]; 
dnew910_uq <- merge(dat3_uq,dnew910_uq,by="id.new")
dnew910_uq$cluster <- dnew910_uq$cluster.local3

fit <- survfit(Surv(month, delta.death) ~  cluster,data = dnew910_uq,start.time=30.08)
myplot[[3]] <- ggsurvplot(fit, data = dnew910_uq, risk.table = TRUE,pval=TRUE,
                          conf.int = TRUE,pval.coord=c(40,0.2), xlab="months",title=expression(paste("spiders (",hat(alpha), " = 0.75)", sep=""))) + 
  guides(colour = guide_legend(nrow = 2))

#-------------------------------------------------------------------#
dat <- data.frame(id.new= dat[[1]]$id,cluster.global = factor(dat[[1]]$cluster.global,label=c(1,2)))
dat_uq <- dat[!duplicated(dat$id.new, fromLast=TRUE),]; 
dnew910_uq <- merge(dat_uq,dnew910_uq,by="id.new")
dnew910_uq$cluster <- dnew910_uq$cluster.global

fit <- survfit(Surv(month, delta.death) ~  cluster,data = dnew910_uq,start.time=30.08)
myplot[[4]] <- ggsurvplot(fit, data = dnew910_uq, risk.table = TRUE,pval=TRUE,conf.int = TRUE,
                          pval.coord=c(40,0.2), xlab="months",title="overall clusters") + 
  guides(colour = guide_legend(nrow = 2))


dev.new(width=180, height=90)
# Arrange multiple ggsurvplots and print the output
arrange_ggsurvplots(myplot, print = TRUE, ncol = 4, nrow = 1 )


#--------------------------------------------------------------------------------#
# Cox model for the association between clusters and time to death
#--------------------------------------------------------------------------------#
# association between the local clustering of the first marker (lbili) and time to death
SurvObj <- with(dnew910_uq[dnew910_uq$alive>910,], Surv(month, delta.death == 1))
fit.coxph <- coxph(SurvObj~cluster.local1,data=dnew910_uq[dnew910_uq$alive>910,])
fit.coxph

# association between the local clustering of the second marker (lsgot) and time to death
SurvObj <- with(dnew910_uq[dnew910_uq$alive>910,], Surv(month, delta.death == 1))
fit.coxph <- coxph(SurvObj~cluster.local2 + age + sex,data=dnew910_uq[dnew910_uq$alive>910,])
fit.coxph

# association between the local clustering of the third marker (lsgot) and time to death
SurvObj <- with(dnew910_uq[dnew910_uq$alive>910,], Surv(month, delta.death == 1))
fit.coxph <- coxph(SurvObj~cluster.local3 + age + sex,data=dnew910_uq[dnew910_uq$alive>910,])
fit.coxph

# association between the global(consensus) clustering and time to death
SurvObj <- with(dnew910_uq[dnew910_uq$alive>910,], Surv(month, delta.death == 1))
fit.coxph <- coxph(SurvObj~cluster.global,data=dnew910_uq[dnew910_uq$alive>910,])
fit.coxph











