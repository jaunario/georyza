# Author: Jorrel Khalil S. Aunario, j.aunario@irri.org
#   
#
#

# Based on 

# TODO add function Using SinExpFlattened TAir function
ttime <- function(t, tmin, tmax, t.sunrise=6, t.sunset=18, DL=12, P=1.5, noc.coeff=4){ 
	# Using SinExp TAir Function
	if (t>=t.sunrise & t<=t.sunset){
		Temp <- tmin+((tmax-tmin) * sin(pi*(t-t.sunrise)/(DL+2*P)))
	} else {
		Tset <-ttime(t=t.sunset,tmin=tmin, tmax=tmax, t.sunrise=t.sunrise, t.sunset=t.sunset, DL=DL, P=P, noc.coeff=noc.coeff)
		Temp <- (tmin-Tset*exp(-(24-DL)/noc.coeff)+(Tset-tmin)*exp(-(t+ifelse(t<t.sunrise,24,0)-t.sunset)/noc.coeff))/(1-exp(-(24-DL)/noc.coeff))
	}
	return(Temp)
}

tpan <- function(tdew,tair, no.cooling=FALSE){
	if(no.cooling) tp <- tair else tp <- tair*0.78+0.073*RH(tdew,tair) 
	return(tp)
}

#TODO
tsuns <- function(time=t.sunset()){
	
}

#TODO
tsunr <- function(time=t.sunrise()){
	
}