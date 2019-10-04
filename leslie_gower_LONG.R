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
# Calculate invasion growth rates and coexistence mechanisms
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
# Analytical calculations of coexistence mechanisms
#=========================================================================

#Approximation 1, without the log(E) and log(C) or log gr
r1 = 3000:4000 #Region where resident is equilibrated
Ei_f = (Fr[r1,1])
Ei = mean(Ei_f)
n2_eq = mean (nrns2[r1])
Ci_f = (1+alphas[1]*(Fr[r1,2]*nrns2[r1])) 
Ci = mean(Ci_f)
Yi = mFr[1]

gr1=log(-(Yi/Ci^2)*cov(Ei_f,Ci_f)+(Ei*Yi)/Ci^3 *var(Ci_f) + (Ei*Yi)/Ci+(1 - Ei)*sr)

#Approximation 2, without the log(E) and log(C) but log(gr)
r1 = 3000:4000 #Region where resident is equilibrated
Ei_f = (Fr[r1,1])
Ei = mean(Ei_f)
n2_eq = mean (nrns2[r1])
Ci_f = (1+alphas[1]*(Fr[r1,2]*nrns2[r1])) 
Ci = mean(Ci_f)
Yi = mFr[1]
D1 = Yi/Ci-sr
D2 = (Ei*Yi)/Ci+(1-Ei)*sr
gr1l=-(D1^2)/(2*(D2^2))*var(Ei_f)+((Ei*Yi*D1)/(Ci^2*D2^2)-Yi/(Ci^2*D2))*cov(Ei_f,Ci_f)+
	((2*Ei*Yi)/(Ci^3*D2-(Ei^2*Yi^2)/(Ci^4*D2^2)))/2 *var(Ci_f) + log(D2)

#Approximation 3, using logs
Eil_f = log(Fr[r1,1])
Eil = mean(Eil_f)
Cil_f = log(1+alphas[1]*(Fr[r1,2]*nrns2[r1])) 
#Cil_f = Eil_f - log(sr)
#Cil = Eil - log(sr)
Cil = mean(Cil_f)
D1l = exp(Eil-Cil)*Yi-exp(Eil*sr)
D2l = exp(Eil-Cil)*Yi+(1-exp(Eil))*sr
D3l = exp(Eil-Cil)*Yi

gr2 = (D1l/D2l-D1l^2/D2l^2)*var(Eil_f)/2+(D3l*D1l/D2l^2-D3l/D2l)*cov(Eil_f,Cil_f)+
(D3l/D2l-D3l^2/D2l^2)*var(Cil_f)/2+ log(D2l);

#Approximation 3 Resident, using logs
Yr = mFr[2]
Erl_f = log(Fr[r1,2])
Erl = mean(Erl_f)
Crl_f = log(1+alphas[1]*(Fr[r1,2]*nrns2[r1])) 
#Cil_f = Eil_f - log(sr)
#Cil = Eil - log(sr)
Crl = mean(Crl_f)
D1rl = exp(Erl-Crl)*Yr-exp(Erl*sr)
D2rl = exp(Erl-Crl)*Yr+(1-exp(Erl))*sr
D3rl = exp(Erl-Crl)*Yr

gr2r = (D1rl/D2rl-D1rl^2/D2rl^2)*var(Erl_f)/2+(D3rl*D1rl/D2rl^2-D3rl/D2rl)*cov(Erl_f,Crl_f)+
(D3rl/D2rl-D3rl^2/D2rl^2)*var(Crl_f)/2+ log(D2rl);

#Approximation 4, no logs etc, original parameters 
Do = (1+alphas[1]*(mean(Fr[r1,2])*mean(nrns2[r1]))) #Same as Ci, or exp(Ci) depeding on model

# gr3 = (alphas[1]*mFr[1])/Do^2* ( (mean(Fr[r1,1])*mean(Fr[r1,2])^2*alphas[1]*var(nrns2[r1]))/Do+(mean(Fr[r1,1])*mean(nrns2[r1])^2*alphas[1]*var(Fr[r1,2]))/
# 		Do-(mean(Fr[r1,2])*cov(Fr[r1,1],Fr[r1,2] )) )+
# 		mean(Fr[r1,1])*mFr[1]/Do+(1-mean(Fr[r1,1]))*sr

gr3=((mean(Fr[r1,1])*mean(Fr[r1,2])^2*alphas[1]^2*mFr[1]*var(nrns2[r1]))/(Do)^3+
	(mean(Fr[r1,1])*mean(nrns2[r1])^2*alphas[1]^2*mFr[1]*var(Fr[r1,2]))/Do^3)-
	(mean(nrns2[r1])*alphas[1]*mFr[1]*cov(Fr[r1,1],Fr[r1,2] ) )/Do^2+
	mean(Fr[r1,1])*mFr[1]/Do+(1-mean(Fr[r1,1]))*sr



###############################################################################
#Specific components, i.e. the storage effect
###############################################################################


#Invader
#Fr= Fr%*%diag(1/colMeans(Fr))
Esi_p=colMeans(Fr)
Esi=colMeans(log(Fr))
Cs_p=(Esi_p[1]*mFr[1])/((Esi_p[1]-1)*sr+1)
Cs = Esi[1]-log((exp(Esi[1])*sr)/mFr[1]-sr/mFr[1]+1/mFr[1])

Cis_p=(Fr[,1]*mFr[1])/((Fr[,1]-1)*sr+1)
Cis = log(Fr[,1])-log((exp(log(Fr[,1]))*sr)/mFr[1]-sr/mFr[1]+1/mFr[1])

Cis_p2= (1+alphas[1]*Fr[n,2]*nrns2[a1:a2]) 
Cis2 = log((1+alphas[1]*Fr[n,2]*nrns2[a1:a2]) 

EEI = exp(log(Fr)-Esi[1])
EI = mean(EEI)

EEIp = ((Esi_p[1]-Fr)*(sr-1))/Esi_p[1]
EIp = mean(EEIp)

c1=(exp(Esi[1])*sr-sr+1)
gr1=(1-sr)*EI+((1-sr)*var(log(Fr[,1])))/2-c1*cov(log(Fr))[1,2]+1/2*c1*var(Cis)

c2=((Esi_p[1])*sr+1)
gr2=(c2-sr)*EIp-c2^2/(Esi_p[1])^2*1/mFr[1]*cov((Fr))[1,2]+c2^3/(Esi_p[1])^2*(1/mFr[1]^2)*var(Cis_p)

m1=log(nrns1[4000:4100])
xx= a1:a2
m1.lm=lm(m1~xx)





#The resident quantities
mFr2=c(mFr[2],mFr[1])
Es = c(Esi[2],Esi[1])
#nhat= ((Es*mFr2)/ (1-sr*(1-Es))-1)/(alphas*Es)
nhat=(Es*mFr2+(1-Es)*sr-1)/((Es^2-Es)*alphas*sr+Es*alphas)
############
#
# This is the key to the two different answers I was getting! 
# This term, whether the "yield" or "intrinsic growth" rate is 
# on the bottom. This reflects where the term appears in the model
# itself, i.e. whether it is in the competition term. 
#
#################################
#nhat= ((Es*mFr2)/ (1-sr*(1-Es))-1)/(alphas[2]*Es*mFr2)
Cs = Es*alphas*nhat+1
vars=apply(Fr,2,var)
vars=c(vars[2],vars[1])
covars = c(cov(Fr)[1,2],cov(Fr)[2,1])
vars_X = c(var(nrns2[1000:2000]), var(nrns1[10000:12000]) )

#Effect of resident term
se_res = Esi*alphas^2*nhat^2*mFr/Cs^3 * vars
se_resb = (Esi*Es^2*alphas^2*mFr)/(Cs)^3 *vars_X

#Effect if invader term
se_inv = alphas*nhat*mFr/Cs^2*covars

#delta_I
delta_I = se_res+se_resb-se_inv

##################
#Chesson terms: Note, these are only equivalent under the condition
#that fluctuations in the environmental variable have been standardized
#by the mean. So Es^2(se_res+se_resb)-Es*Esi(se_inv) = delta_IC
##################
Beta_i = (1-sr*(1-Esi))/Cs

#Effect of resident term
se_resC = Beta_i*(Cs-1)* vars
#Effect if invader term
se_invC = Beta_i*(Cs-1)*covars

#delta_I
delta_IC = se_resC-se_invC

