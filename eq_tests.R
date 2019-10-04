# Tests for equilibrium in stochastic population models. 
# git rebase origin/master
# git push -u origin master
# git pull git@github.com:jusinowicz/info_theory_eco.git master

#=========================================================================
## Load these libraries
#=========================================================================

#=========================================================================
# Find stationary parts of a time series based on the Ljung-Box test for
# autocorrelation. If the p-value is significant, then the signal is 
# non-stationary.  
# x  		data 
# window	window size of the time series to look at 
# max.lag	largest lag size
#=========================================================================

eq_LB = function (x, window = 10, max.lag = 10 ){

	nt = length(x)
	pvals = matrix(0,nt,1) #Store the p-values
	for (t in 1:(nt-window)) { 
		pvals[t] = unlist((Box.test(x[t:(t+window)],lag=max.lag, type="Ljung-Box"))$p.value)
	}

	return(pvals)
}

#=========================================================================
# Find the equilibrium based on where the derivative of the smoothed 
# function reaches 0. Smoothing is done by implementing a moving average.
# x  		data 
# window	window size of the time series for MA
#=========================================================================
eq_diff = function (x, window = 10){
	nt = length(x)
	dat_temp = matrix(0,nt,2) #Store the smooth and the difference
	for (t in 1:(nt-window)) {
		dat_temp[t,1] = mean(x[t:(t+window)]) 
	}
	dat_temp[,2]=c(diff(dat_temp[,1]),0)
	return(dat_temp)
}
