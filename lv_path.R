#=============================================================================
#Implementation of a lotka-volterra competition model with a density-dependent 
# pathogen mortality term (from Schnitzer et al. 2011)
#=============================================================================
# load libraries
#=============================================================================
library(tidyverse)
library(deSolve)
library(gridExtra)

#=============================================================================
# Define the population dynamics through the following function
#=============================================================================
lv_path= function(times,sp,parms){
	with( as.list(c(parms, sp )),
		{		
		
			nspp = parms$nspp 
			Ni = matrix(sp[1:nspp], nspp, 1)

			dNi = Ni 
			for( i in 1:nspp){
				#Lotka-volterra consumption with a saturating mortality term 
				dNi[i] = Ni[i] * ( ri[i] *(1 - t(alphas[i,])%*%Ni )  
					- Ni[i] / (Ni[i]+Hp[i] ) )
			}

	  	list( c(dNi) )
		})	

}  

#=============================================================================
# Set values of parameters
#=============================================================================
#How many species to start? 
nspp = 20

spp_prms = NULL
spp_prms$ri = matrix(rpois(nspp,10), nspp, 1) #intrinsic growth
spp_prms$Hp = matrix(rpois(nspp,100), nspp, 1) #intrinsic growth

#Competition coefficients: 
spp_prms$alphas = matrix( runif(nspp^2,  0.1, 0.9),nspp,nspp )
diag(spp_prms$alphas) = matrix(rnorm(nspp, 0.99,0.01),nspp,1) #Set the intraspecific alphas = 1

#Pass all of these parameters as a list
parms = list(
	nspp=nspp, ri = spp_prms$ri, Hp = spp_prms$Hp, alphas = spp_prms$alphas
 )

#=============================================================================
# Run the model with initial conditions 
#=============================================================================
tend = 100
delta1 = 0.1
times  = seq(from = 0, to = tend, by = delta1)
tl = length(times)
minit = c( matrix(0.001,nspp,1) )

lv_out = dede(y=minit, times=times, func=lv_path, parms=parms, atol = 1e-9)
lv_out = as.data.frame(lv_out)

#=============================================================================
# Plot
#=============================================================================
lv_long =lv_out %>% gather( species, N, 2:(nspp+1) )

#Total number of species still around at each time: 
lv_long = lv_long %>% group_by(time) %>% mutate(spp_tot = sum(N>0.001) )

#Their collective biomass: 
lv_long = lv_long %>% group_by(time) %>% mutate(biomass = sum(N[N>0.001] ) )

#Each species' time trajectory
p1=ggplot()+ geom_line( data = lv_long, aes ( x = time, y = N, color = species)  )+ 
ylab("Population")+
theme(axis.text.x=element_blank(), axis.title.x=element_blank(), legend.position = "none") 

#Total number of species
p2=ggplot()+ geom_line( data = lv_long, aes ( x = time, y = spp_tot)  )+
ylab("Number of species")+
theme(axis.text.x=element_blank(), axis.title.x=element_blank(), legend.position = "none") 

#Total biomass of community
p3=ggplot()+ geom_line( data = lv_long, aes ( x = time, y = biomass)  )+
ylab("Biomass")

p4 = grid.arrange(p1,p2,p3, nrow = 3 )

ggsave("./feedbacks1.pdf", p4, width = 8, height = 10)



