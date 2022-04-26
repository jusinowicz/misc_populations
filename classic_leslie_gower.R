###############################################################################
#
# Population dynamics of the SPATIALLY IMPLICIT Leslie-Gower 
# (aka annual plant model) for 2 species, with calculation of 
# invasion growth rates and coexistence mechanisms
#
###############################################################################

#=========================================================================
## Load these libraries
#=========================================================================
library(MASS)
library(nlstools) #For prediction intervals
library(tidyverse)
library(gridExtra)
require(grid)
source("./eq_tests.R")

#=========================================================================
# Population dynamics implemented numerically with reciprocal invasions
#=========================================================================
#########
#Tunable
#########
ngens=16000 #Number of generations 
mFr=matrix(c(5,5))   #Mean reproduction rates
mFr2=matrix(c(10,10))   #Mean reproduction rates
sig_Fr= matrix(c(1, -1, -1, 1),2,2) #Interspecific covariance of reproduction
sig_Fr2= matrix(c(1, 0, 0, 1),2,2) #Interspecific covariance of reproduction

alphas=matrix( c(1,0.8,0.8,1),2,2) #Alpha coefficients, competitive interactions
sr=0.9 #Survival
#Invasion times (divide into quarters is the easiest)
invasions =c(1, floor(ngens/4), floor(ngens/2), floor(ngens*3/4) )

#########
#Internal
#########
Fr=abs(mvrnorm(ngens,mFr,sig_Fr)) #Take the abs value of correlated normal
Fr2=abs(mvrnorm(ngens,mFr2,sig_Fr2)) #Take the abs value of correlated normal

#Standardize values between 0 and 1 based on maximum.
# Fr=Fr/(matrix(apply(Fr,2,max),ngens,2,byrow=T)) 
#Population matrixes
nrns1=matrix(0,ngens,1)
nrns2=matrix(0,ngens,1)
nrns3=matrix(0,ngens,1)
nrns4=matrix(0,ngens,1)
lggr1=matrix(0,ngens,1)
lvgr1=matrix(0,ngens,1)

nrns2[1]=0.1 #Seed species 2 and allow it to establish as resident
nrns4[1]=0.1 #Seed species 2 and allow it to establish as resident


for (n in 1:(ngens-1)) {

	#Invasion: Seed or set to near-zero at the appropriate time steps:
	if (n==invasions[2]) { nrns1[n]=0.0001;nrns3[n]=0.0001 }

	#Switch the roles of resident and invader between species 1 and 2
	if (n==invasions[3]) { nrns1[n]=0.1;nrns3[n]=0.1 
							nrns2[n]=0;nrns4[n]=0}
	#Second invasion
	if (n==invasions[4]) { nrns2[n]=0.0001;nrns4[n]=0.0001 }


	#LG model
	nrns1[n+1] = sr*nrns1[n]+Fr[n,1]*nrns1[n]/(1+alphas[1,1]*nrns1[n]+alphas[1,2]*nrns2[n])
	nrns2[n+1] = sr*nrns2[n]+Fr[n,2]*nrns2[n]/(1+alphas[2,1]*nrns1[n]+alphas[2,2]*nrns2[n])
	lggr1[n+1] = sr+Fr[n,1]/(1+alphas[1,1]*nrns1[n]+alphas[1,2]*nrns2[n])
	#LV model
	# nrns3[n+1] = nrns3[n]+0.5*nrns3[n]*(1-alphas[1,1]/Fr2[n,1]*nrns3[n]-alphas[1,2]/Fr2[n,1]*nrns4[n])
	# nrns4[n+1] = nrns4[n]+0.5*nrns4[n]*(1-alphas[2,1]/Fr2[n,2]*nrns3[n]-alphas[2,2]/Fr2[n,2]*nrns4[n])

	nrns3[n+1] = nrns3[n]+0.5/(Fr2[n,1])*nrns3[n]*(Fr2[n,1]-alphas[1,1]*nrns3[n]-alphas[1,2]*nrns4[n])
	nrns4[n+1] = nrns4[n]+0.5/(Fr2[n,2])*nrns4[n]*(Fr2[n,2]-alphas[2,1]*nrns3[n]-alphas[2,2]*nrns4[n])
	lvgr1[n+1] = 1+0.5/(Fr2[n,1])*(Fr2[n,1]-alphas[1,1]*nrns3[n]-alphas[1,2]*nrns4[n])

	# nrns3[n+1] = nrns3[n]+1.5*nrns3[n]*(1-alphas[1,1]*nrns3[n]-alphas[1,2]*nrns4[n])+nrns3[n]*Fr2[n,1]
	# nrns4[n+1] = nrns4[n]+1.5*nrns4[n]*(1-alphas[2,1]*nrns3[n]-alphas[2,2]*nrns4[n])+nrns4[n]*Fr2[n,2]
	# lvgr1[n+1] = 1+1.5*(1-alphas[1,1]*nrns3[n]-alphas[1,2]*nrns4[n]) + Fr2[n,1]
}

tu1 = 1
tu2 = ngens
# tu1 = 4000
# tu2 = 4800
plot(nrns2[tu1:tu2],t="l",ylim=c(0,60))
lines(nrns1[tu1:tu2],col="red") 
plot(nrns3[tu1:tu2],t="l",ylim=c(0,12))
lines(nrns4[tu1:tu2],col="red") 

tu1 = 4000
tu2 = 4800
d1_tmp = data.frame(lg1 =nrns1[tu1:tu2], lg2 =nrns2[tu1:tu2], time = 1:length(nrns2[tu1:tu2]) )
d1 = gather(d1_tmp, "lg1","lg2", key=species, value=Population  )
d2_tmp = data.frame(lv1 =nrns3[tu1:tu2], lv2 =nrns4[tu1:tu2], time = 1:length(nrns3[tu1:tu2]) )
d2 = gather(d2_tmp, "lv1","lv2", key=species, value=Population  )

c_use = c("#440154FF","#35B779FF" )
p1 = ggplot() + geom_line(data=d1, aes(x=time, y = Population,color=species))+
	scale_color_manual(values=c_use) +
	theme_bw() + theme(
	text = element_text(size=20),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none")
p2 = ggplot() + geom_line(data=d2, aes(x=time, y =Population,color=species))+
	scale_color_manual(values=c_use) +
	theme_bw() + theme(
	text = element_text(size=20),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none")


g=grid.arrange(p2, p1,nrow=1)

ggsave(file="populations.jpg", g)


#Growth rate functions: 
tu1 = 8000
tu2 = 9000
grlg1 = nrns1[2:ngens]/nrns1[1:(ngens-1)] 
grlv1 = nrns3[2:ngens]/nrns3[1:(ngens-1)] 
x1 = 0:200
x1a = sort(nrns1[tu1:tu2])
x1b = sort(nrns3[tu1:tu2])
grm_lg = sr+mean(Fr[,1])/(1+alphas[1,1]*x1a+alphas[1,2]*mean(nrns2[tu1:tu2]) )
grm_lv = 1+mean(0.5/(Fr2[,1]) )*(mean(Fr2[,1])-alphas[1,1]*x1b-alphas[1,2]*mean(nrns4[tu1:tu2]) )

gr1a = data.frame(gr1= grlg1[tu1:tu2], Population = nrns1[(tu1-1):(tu2-1)]  )
#gr1a = data.frame(gr1= grlg1, Population = x2a  )
gr1b = data.frame(gr1m= grm_lg, Population = x1a)
gr2a = data.frame(gr1 = grlv1[tu1:tu2], Population = nrns3[(tu1-1):(tu2-1)] )
#gr2a = data.frame(gr1 = grlv1, Population = x2b )
gr2b =  data.frame(gr1m= grm_lv, Population =x1b )

grm_lgM = sr+mean(Fr[,1])/(1+alphas[1,1]*mean(x1a)+alphas[1,2]*mean(nrns2[tu1:tu2]) )
grm_lvM = 1+mean(0.5/(Fr2[,1]))*(mean(Fr2[,1])-alphas[1,1]*mean(x1b)-alphas[1,2]*mean(nrns4[tu1:tu2]) )
grm_lgMM = mean(grlg1[tu1:tu2])
grm_lvMM = mean(grlv1[tu1:tu2])

p1 = ggplot() + geom_point(data=gr1a, aes(x=Population, y = gr1))+
	geom_line(data=gr1b, aes(x=Population, y = gr1m))+
	scale_color_manual(values=c_use[1])  +
	xlim(c(0,60))+
	theme_bw() + theme(
	text = element_text(size=20),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none")
p2 = ggplot() + geom_point(data=gr2a, aes(x=Population, y = gr1))+
	geom_line(data=gr2b, aes(x=Population, y = gr1m))+
	xlim(c(0,10))+
	ylim(c(0.5,1.5))+
	scale_color_manual(values=c_use[1])  +
	theme_bw() + theme(
	text = element_text(size=20),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none")


g=grid.arrange(p1, p2,nrow=1)

ggsave(file="populations.jpg", g)
#=========================================================================
# Calculate invasion growth rates 
#=========================================================================

#Fit a linear model to get invasion growth rates
#These starting and ending times are just being tuned by eye: 
a1=4000
a2=4100
m1=log(nrns1[a1:a2])
xx= a1:a2
m1.lm=lm(m1~xx)

#or take the log(gr1[inv+1]/gr1[inv]) to get invasion growth rates
ni = a1
nj=a1+1
m1.d = log( (sr*nrns1[nj]+Fr[nj,1]*nrns1[nj]/(1+alphas[1,1]*nrns1[nj]+alphas[1,2]*nrns2[nj]))/
(sr*nrns1[ni]+Fr[ni,1]*nrns1[ni]/(1+alphas[1,1]*nrns1[ni]+alphas[1,2]*nrns2[ni])) )

#=========================================================================
# Try recovering original parameters by fitting with NLS
#=========================================================================
#Try all parameters at once, all data at once: 
dat1 = nrns1[1:(ngens-1)]
dat2 = nrns1[2:ngens] #One time step forward
#Fit the model using the growth rate del1 = (change in time)
delta_dat = data.frame( del1 = dat2/dat1, dat1=dat1, datb =nrns2[1:(ngens-1)]  ) 
delta_dat = delta_dat[is.finite(delta_dat[,1]), ] #Only keep finite rows! 
start1 = c(r=1.1, a1=0.5, a2=0.5) #The starting values of parameters that NLS will fit
lg_fit = nls(del1 ~ 0.9+r/(1+a1*dat1+a2*datb), data = data.frame(delta_dat), 
	start=start1)
summary(lg_fit)

### Get prediction intervals for NLS fit. 
lg_fit_boot = nlsBoot(lg_fit) #Use nlstools to repeatadly sample (bootstrap)
#Use the nlsBoot object to make intervals
lg_pred1 = nlsBootPredict (lg_fit_boot, newdata=delta_dat[,2:3], interval = "prediction")
lg_dat =data.frame(
	t = 1:length (delta_dat$del1), del1= delta_dat$del1, pred = lg_pred1[,1], 
	lower = lg_pred1[,2], upper = lg_pred1[,3]) #Combine the data frame with predictions

#An ugly plot just to look at the results
ggplot( )+geom_line(data = lg_dat, aes(x =t, y =del1,col = "Simulation" )) +
			geom_line(data = lg_dat, aes(x =t, y =pred,col="Prediction" )) +
			geom_line(data = lg_dat, aes(x =t, y =lower,col ="PI" )) +
			geom_line(data = lg_dat, aes(x =t, y =upper,col="PI" )) +
			ylab("Population growth rate")+ xlab("Time") + 
			scale_y_continuous(limits = c(0.5, 1.5)) 

#Try only interspecific competition, in the invasion regions: 
#These starting and ending times are just being tuned by eye: 
a1=4000
a2=4100
dat1 = nrns1[a1:a2]
dat2 = nrns1[(a1+1):(a2+1)] #One time step forward
#Fit the model using the growth rate del1 = (change in time)
delta_dat = data.frame( del1 = dat2/dat1, dat1=dat1, datb =nrns2[a1:a2]  ) 
delta_dat = delta_dat[is.finite(delta_dat[,1]), ] #Only keep finite rows! 
start1 = c(r=1.1, a2=0.5) #The starting values of parameters that NLS will fit
lg_fit2 = nls(del1 ~ 0.9+r/(1+a2*datb), data = delta_dat, 
	start1)
summary(lg_fit2) #This does poorly, perhaps because of too little data? 

#Try only interspecific competition, in the invasion region, and assume R is known: 
dat1 = nrns1[a1:a2]
dat2 = nrns1[(a1+1):(a2+1)] #One time step forward
#Fit the model using the growth rate del1 = (change in time)
delta_dat = data.frame( del1 = dat2/dat1, dat1=dat1, datb =nrns2[a1:a2]  ) 
delta_dat = delta_dat[is.finite(delta_dat[,1]), ] #Only keep finite rows! 
start1 = c(a2=0.5) #The starting values of parameters that NLS will fit
lg_fit3 = nls(del1 ~ 0.9+mean(Fr[,1])/(1+a2*datb), data = delta_dat, 
	start1)
summary(lg_fit3) #This does well 
### Get prediction intervals for NLS fits. 
lg_fit_boot3 = nlsBoot(lg_fit3)
lg_pred3 = nlsBootPredict (lg_fit_boot3, newdata=delta_dat, interval = "prediction")

#=========================================================================
# Experimental section: How well do stats for equilibrium work?  
# Use the time series of population growth, nrns1, to examine the 
# equilibrium properties of species 1. If you wanted to run any of this on
# real data, just replace "nrns1" with your matrix/data.frame column of 
# population time series data. 
#=========================================================================

#This version smooths the data with a moving average and then returns the
#1st derivative at each point in the smoothed curve. Equilibrium should happen 
#when the 1st derivative = 0
#This is sensitive to the window size when the data are noisy. 
#Try different window sizes for the moving average. 
ed1 = eq_diff(nrns1[invasions[2]:(invasions[2]+1000)],window=3)
ed1 = eq_diff(nrns1[invasions[2]:(invasions[2]+1000)],window=100)
plot(ed1[1:500,2])

#This is a standard test for stationarity in a time series. The test is 
#actually for non-stationarity, so low p-values indicate non-stationarity. 
#This test does not actually work very well when trying to identify a region
#of stationarity (i.e. equilibrium) in an otherwise stationary series.  
eLB1=eq_LB(nrns1[invasions[2]:(invasions[2]+1000)], window = 10, max.lag = 10 )
eLB1=eq_LB(nrns1[invasions[2]:(invasions[2]+1000)], window = 100, max.lag = 100 )
plot(eLB1[1:500])

#Fit a logistic growth model to the data instead of testing for any sort 
#of equilibrium conditions. Use NLS to find K 
dat1 = nrns1[invasions[2]:(invasions[2]+1000)]
dat2 = nrns1[(invasions[2]+1):(invasions[2]+1000+1)] #One time step forward
#Fit the model using the growth rate del1 = (change in time)
delta_dat = data.frame( del1 = dat2/dat1, dat1=dat1 ) 
start1 = c(r=1.1, K=1.1) #The starting values of parameters that NLS will fit
logistic_fit = nls(del1 ~ 1+r*(1-dat1/K),data = delta_dat, start1)
summary(logistic_fit)

#If you want to get confidence intervals and plot them as well.
###This plot is actually kind of neat. It plots the predicted
#growth rate (del1) as a function of the density of species 1 (dat1)
#and not as a function of time. By relating the data this way, it shows how 
#most of the spread in the data that drives the SE/CI width happens at the 
#low densities and that there is a lot less SE around the equilibrium density
#itself. This might be a useful visual tool for real data to gauge if the 
#population growth has convgerged on an equilibrium. 

delta_dat$pred = predict(logistic_fit) #Get the fits from the model
se = summary(logistic_fit)$sigma #Sigma
ci = outer(delta_dat$pred, c(outer(se, c(-1,1), '*'))*1.96, '+') #CIs
ii = order(delta_dat$dat1)

# shaded area plot
low = ci[ii,1]; high = ci[ii,2]; base = delta_dat[ii,'dat1']
polygon(c(base,rev(base)), c(low,rev(high)), col='grey')
with(delta_dat[ii,], lines(dat1, pred, col='blue'))
with(delta_dat, points(dat1, del1))