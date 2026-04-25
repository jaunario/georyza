# Author: Jorrel Khalil Aunario j.aunario@irri.org
# Title
# 
###############################################################################

spikelet.fertility <- function(temp, method="Horie"){
	value <- switch(method, 
				Horie=1/(1+exp(0.853*(temp-36.6))),
				Abeysiriwarden=pmin(1,pmax(0,2.7121-0.0782*temp)), #TODO for further investigation
				Jagadish=exp(14.3-0.408*temp)/(1+exp(14.3-0.408*temp)), #TODO for further investigation
				Weerakoon=pmax(0,1-exp(0.3904*(temp-36.18))), #TODO for further investigation
				Yan=pmin(1, pmax(0,1.89-0.0428*temp)), #TODO for further investigation
				Julia1=pmin(1,pmax(0,3.08-0.0761*temp)), #TODO for further investigation
				Julia2=pmin(1, pmax(0,1-0.2*(temp-29.7))),
				NA) #TODO for further investigation
	if(sum(is.na(value))>0){
		warning("Unrecognized method. Kindly check the spelling of the method you specified.")
	}
	return(value)
}
