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
source("./eq_tests.R")

#=========================================================================
# Population dynamics implemented numerically with reciprocal invasions
#=========================================================================
#########
#Tunable
#########
ngens=16000 #Number of generations 
mFr=matrix(c(5,5))   #Mean reproduction rates
sig_Fr= matrix(c(1, -1, -1, 1),2,2) #Covariance matrix
alphas=matrix( c(0.5,0.5),2,1) #competitive interactions
sr=0.9 #Survival
#Invasion times (divide into quarters is the easiest)
invasions =c(1, floor(ngens/4), floor(ngens/2), floor(ngens*3/4) )

#########
#Internal
#########
Fr=abs(mvrnorm(ngens,mFr,sig_Fr)) #Take the abs value of correlated normal
#Standardize values between 0 and 1 based on maximum.
Fr=Fr/(matrix(apply(Fr,2,max),ngens,2,byrow=T)) 
#Population matrixes
nrns1=matrix(0,ngens,1)
nrns2=matrix(0,ngens,1)
#Exponential form of the model: 
nrns1E=matrix(0,ngens,1)
nrns2E=matrix(0,ngens,1)

nrns2[1]=0.1 #Seed species 2 and allow it to establish as resident
nrns2E[1]=0.1

for (n in 1:(ngens-1)) {

	#Invasion: Seed or set to near-zero at the appropriate time steps:
	if (n==invasions[2]) { nrns1[n]=0.0001; nrns1E[n]=0.0001  }
	#Switch the roles of resident and invader between species 1 and 2
	if (n==invasions[3]) { nrns1[n]=0.1; nrns1E[n]=0.1
							nrns2[n]=0; nrns2E[n]=0   }
	#Second invasion
	if (n==invasions[4]) { nrns2[n]=0.0001;nrns2E[n]=0.0001 }


	#Spatially implicit annual plant model
	nrns1[n+1] = sr*(1-Fr[n,1])*nrns1[n]+Fr[n,1]*mFr[1]*nrns1[n]/(1+alphas[1]*(nrns1[n]*Fr[n,1]+Fr[n,2]*nrns2[n]))
	nrns2[n+1] = sr*(1-Fr[n,2])*nrns2[n]+Fr[n,2]*mFr[2]*nrns2[n]/(1+alphas[2]*(nrns1[n]*Fr[n,1]+Fr[n,2]*nrns2[n]))

	#Exponential form of the model
	nrns1E[n+1] = sr*(1-exp(log(Fr[n,1])))*nrns1E[n]+nrns1E[n]*mFr[1]*exp( log(Fr[n,1])- log(1+alphas[1]*(nrns1E[n]*Fr[n,1]+Fr[n,2]*nrns2E[n])))
	nrns2E[n+1] = sr*(1-exp(log(Fr[n,2])))*nrns2E[n]+nrns2E[n]*mFr[2]*exp( log(Fr[n,2])-log(1+alphas[2]*(nrns1E[n]*Fr[n,1]+Fr[n,2]*nrns2E[n])))


}

plot(nrns2,t="l")
lines(nrns1,col="red") 

#=========================================================================
# Calculate invasion growth rates 
#=========================================================================

#Fit a linear model to get invasion growth rates
a1=4000
a2=4100
m1=log(nrns1[a1:a2])
xx= a1:a2
m1.lm=lm(m1~xx)

#or take the log(gr1[inv+1]/gr1[inv]) to get invasion growth rates
ni = a1
nj=a1+1
m1.d = log( (sr*(1-Fr[nj,1])*nrns1[nj]+Fr[nj,1]*mFr[1]*nrns1[nj]/(1+alphas[1]*(nrns1[nj]*Fr[nj,1]+Fr[nj,2]*nrns2[nj])))/
(sr*(1-Fr[ni,1])*nrns1[ni]+Fr[ni,1]*mFr[1]*nrns1[ni]/(1+alphas[1]*(nrns1[ni]*Fr[ni,1]+Fr[ni,2]*nrns2[ni]))) )

#=========================================================================
# How well do stats for equilibrium work?  
# Use the time series of population growth, nrns1, to examine the 
# equilibrium properties of species 1. 
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