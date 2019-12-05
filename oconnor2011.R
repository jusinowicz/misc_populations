#==============================================================================
# Simulate the models from O'Connor et al. 2011, which are essentially just 
# several variations on the MacArthur consumer-resource model. 
# This uses the library deSolve for the ODEs
# This code contains the following parts so far:
# 1. Basic ODE simulation with constant temp-dependent parameters
# 	 Right now only one of the models is implemented. This is model 4 in the 
#	 table. 
#==============================================================================
#==============================================================================
#Load Libraries
#==============================================================================
library(deSolve)

#==============================================================================
###Time to simulate over: 
tend = 1000 #Length of numerical simulation
times  = seq(from = 0, to = tend, by = 0.01)
tl = length(times)

### Numbers of species 
nPsp = 1
nHsp = 1
nspp = nPsp + nHsp

###Activation energies as per table 1
rEp = 0.32
KEp = -0.32
aEh = 0.65
mEh = 0.65
kb = 8.617e-5# 1.38064852*10^(-23)
K = 273

###Initial values? 
rp0 = 0.62
Kp0 = 5680
ah0 = 2
mh0 = 0.1
b0  = 100
e0 = 0.1

#Range of temperatures to explore: 
###For now, while I'm still troubleshooting, just use 1 temp.

tmin = 20 +K
tmax = tmin+1#30 +K
tby = 1
ntemps = (tmax-tmin)/tby
temps = seq(tmin,tmax,tby)

#Output of each temp
out1 = list(matrix(0,ntemps,1))



for (t in 1:ntemps ) { 
	print(t)
	T = temps[t]

	###Calculate model parameters with activation energies
	rp = rp0#*exp(-rEp/(T*kb))  #These parameter values are crazy small! wtf? 
	Kp = Kp0#*exp(-KEp/(T*kb))  
	ah = ah0#*exp(-aEh/(T*kb))  
	mh = mh0#*exp(-mEh/(T*kb))  
	b0 = b0
	e0 = e0

	#=============================================================================
	# Inner loop. Run the model
	#=============================================================================

	parms = list(nspp=nspp, nPsp = nPsp, nHsp = nHsp,
			rp = rp, Kp =Kp,
			ah = ah, mh = mh, b0 = b0, e0 = e0
		 )

	###Define the model for deSolve
	###This is model 4 in Table 2 of OConnor et al. 2011
	oconnor_model4 = function(times,sp,parms){
			     P = matrix(sp[1:nPsp],nPsp, 1)
			     H = matrix(sp[(nPsp+1):(nPsp+nHsp)],nHsp, 1)

				###Resource dynamics: Logistic growth, reduced by consumption
				dP = P
				for( i in 1:nHsp){
					dP[i] = P[i]*( (rp[i]) * (1 - P[i]/Kp[i]) - (t(ah[i])%*%H))
				}

				###Consumer dynamics: LV consumption
				dH = H 
				for( i in 1:nHsp){
					dH[i] = H[i] * ( ah[i]%*%P - H[i]*mh[i] )
				}

			return(list(c(dP,dH)))
		}

	###Run the ODE and store the output
	winit = c(matrix(3,nspp,1))
	out=NULL
	out_temp = ode(y=winit,times=times,func=oconnor_model4,parms=parms)
	out$out = out_temp
	out$parms=parms
	out1[t]= list(out)


	plot(out$out[,"1"],t="l",col="red",ylim = c(0,max(out$out[ ,2:(nspp+1)],na.rm=T)))
	for( n in 2:(nspp) ) {
		lines(out$out[,paste(n)],t="l",col="blue")
	}

}