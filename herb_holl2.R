#==============================================================================
#A model of herbivory with Holling 2 handling time.
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

#######
#######Scenario 1: Keep sea urchins (Y) fixed and change starting pop of seagrass (P)
#######

#model parameters: one value for each iteraction of the model 
#make sure these are all the same length
nPsp = 1 #Just one plant species
b = c(matrix(1.5,100,1) ) #c( seq(1.5,2.5,0.5) ) #intrinsic growth -- this is varied in paper
N = c(matrix(100, length(b),1) ) #carrying capacity
m = c(matrix(0.01, length(b),1) ) #mortality rate
a = c(matrix(0.01, length(b),1) ) #encounter rate
tau = c(matrix(0.23, length(b),1) ) #handling time -- this is varied in paper
Y = c(matrix(19, length(b),1) ) #c(seq(20,40,1)) #herbivore density -- this is varied in paper

nruns = length(b)

#Output of each parameter iteration 
out1 = vector("list",nruns)

#Could add more loops to go over all of the varied parameters. 
for (n in 1:nruns) {


	#=============================================================================
	# Inner loop. Run the model
	#=============================================================================

	parms = list(nPsp = nPsp, b = b[n], N=N[n], m = m[n], a=a[n], Y=Y[n], tau=tau[n])
		

		###Define the model for deSolve
		plant_mod1 = function(times,y,parms){

			with(as.list(c(y,parms)),{ 
			 	###Seagrass dynamics: Logistic growth, reduced by Holling 2 consumption
				dP = P*b * (1 - P/N) - P*m - a*Y/(1+a*tau*P) *P  
				list(dP)
				})

		}


	###Run the ODE and store the output
	#winit = c(matrix(3,nspp,1))
	winit = c(P=n)

	out=NULL
	out_temp = ode(y=winit,times=times,func=plant_mod1,parms=parms)
	out$out = out_temp
	out$parms=parms
	out1[n]= list(out)

	###Plot population densities. Red is the primary producer
	plot(out$out[,"P"],t="l",col="red", 
		xlab="Time",ylab = "Population")


}

#Plot all of the trajectories
plot(out1[[1]]$out[1:1000,"P"],t="l", ylim = c(0,110),
		xlab="Time",ylab = "Population")
for (n in 2:100) {

	lines(out1[[n]]$out[1:1000,"P"])
}

#For tau = 0.23, they al seem to go to the same final density, no matter what P_0 is. 
#Is it a matter of finding the right tau? 

#######
#######Scenario 2: Keep P_0 and Y fixed, and vary tau to look for an ASS? (lol)

#model parameters: one value for each iteraction of the model 
#make sure these are all the same length
nPsp = 1 #Just one plant species
tau = c(seq(0.05, .43, .01) ) #handling time -- this is varied in paper
Y = c(matrix(20,length(tau),1)) #herbivore density -- this is varied in paper
b = c(matrix(1.5,length(Y),1) ) #c( seq(1.5,2.5,0.5) ) #intrinsic growth -- this is varied in paper
N = c(matrix(100, length(b),1) ) #carrying capacity
m = c(matrix(0.01, length(b),1) ) #mortality rate
a = c(matrix(0.01, length(b),1) ) #encounter rate

nruns = length(b)

#Output of each parameter iteration 
out2 = vector("list",nruns)

#Could add more loops to go over all of the varied parameters. 
for (n in 1:nruns) {


	#=============================================================================
	# Inner loop. Run the model
	#=============================================================================

	parms = list(nPsp = nPsp, b = b[n], N=N[n], m = m[n], a=a[n], Y=Y[n], tau=tau[n])
		

		###Define the model for deSolve
		plant_mod1 = function(times,y,parms){

			with(as.list(c(y,parms)),{ 
			 	###Seagrass dynamics: Logistic growth, reduced by Holling 2 consumption
				dP = P*b * (1 - P/N) - P*m - a*Y/(1+a*tau*P) *P  
				list(dP)
				})

		}


	###Run the ODE and store the output
	#winit = c(matrix(3,nspp,1))
	winit = c(P=10) # Based on figure, at this value if ASS is true then some should go to 0

	out=NULL
	out_temp = ode(y=winit,times=times,func=plant_mod1,parms=parms)
	out$out = out_temp
	out$parms=parms
	out2[n]= list(out)

	###Plot population densities. Red is the primary producer
	plot(out$out[,"P"],t="l",col="red", 
		xlab="Time",ylab = "Population")


}


#######
#######Scenario 2: Keep seagrass fixed (P) and let (Y) change 
#######

#model parameters: one value for each iteraction of the model 
#make sure these are all the same length
nPsp = 1 #Just one plant species
Y = c(seq(20,40,1)) #herbivore density -- this is varied in paper
b = c(matrix(2,length(Y),1) ) #c( seq(1.5,2.5,0.5) ) #intrinsic growth -- this is varied in paper
N = c(matrix(100, length(b),1) ) #carrying capacity
m = c(matrix(0.01, length(b),1) ) #mortality rate
a = c(matrix(0.01, length(b),1) ) #encounter rate
tau = c(matrix(0.37, length(b),1) ) #handling time -- this is varied in paper

nruns = length(b)

#Output of each parameter iteration 
out2 = vector("list",nruns)

#Could add more loops to go over all of the varied parameters. 
for (n in 1:nruns) {


	#=============================================================================
	# Inner loop. Run the model
	#=============================================================================

	parms = list(nPsp = nPsp, b = b[n], N=N[n], m = m[n], a=a[n], Y=Y[n], tau=tau[n])
		

		###Define the model for deSolve
		plant_mod1 = function(times,y,parms){

			with(as.list(c(y,parms)),{ 
			 	###Seagrass dynamics: Logistic growth, reduced by Holling 2 consumption
				dP = P*b * (1 - P/N) - P*m - a*Y/(1+a*tau*P) *P  
				list(dP)
				})

		}


	###Run the ODE and store the output
	#winit = c(matrix(3,nspp,1))
	winit = c(P=10)

	out=NULL
	out_temp = ode(y=winit,times=times,func=plant_mod1,parms=parms)
	out$out = out_temp
	out$parms=parms
	out2[n]= list(out)

	###Plot population densities. Red is the primary producer
	plot(out$out[,"P"],t="l",col="red", 
		xlab="Time",ylab = "Population")


}

#Plot all of these trajectories:
par(mfrow =c (1,2) )

plot(out2[[1]]$out[1:1000,"P"],t="l", ylim = c(0,110),
		xlab="Time",ylab = "Population")
for (n in 2:100) {

	lines(out2[[n]]$out[1:1000,"P"])
}

plot(out2a[[1]]$out[1:1000,"P"],t="l", ylim = c(0,110),
		xlab="Time",ylab = "Population")
for (n in 2:100) {

	lines(out2a[[n]]$out[1:1000,"P"])
}

T = tau*N
1/2 *(1 - (1+m/(b*a*T))/(a*T) ) + 1/2*sqrt(1+( (1-m/(b*a*T))/(a*T) )^2 - 4*Y/(b*T) )

1/2 *(1 - (1+m/(b*a*T))/(a*T) ) - 1/2*sqrt(1+( (1-m/(b*a*T))/(a*T) )^2 - 4*Y/(b*T) )

