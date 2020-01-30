# Authors: Raffaele Argiento, Francesco Nutini, Federico Filipponi, Giacinto Manfron, Mirco Boschetti
# Editor: Jorrel Khalil Aunario
# Consiglio Nazionale delle Ricerche (CNR-IMATI and CNR-IREA) - Italy
# Date :  Dec 2012
# Version 1.1
# Licence GPL v3
# Maintainer: Mirco Boschetti <boschetti.m@irea.cnr.it>, Federico Filipponi <federico.filipponi@gmail.com>
# Description: Compute PhenoRice algorithm for rice detection from MODIS data (product MOD09A1)
# PhenoRice Algorithm Reference: Manfron, G., Crema, A., Boschetti, M., Confalonieri, R., 2012. 
#Testing automatic procedures to map rice area and detect phenological crop information exploiting 
#time series analysis of remote sensed MODIS data 8531, 85311E-85311E-11.

#smoothlib <- paste(ifelse(Sys.info()["machine"] %in% c("x86-64","x86_64"),"x86_64","i386"),"/smoothTS",ifelse(Sys.info()["sysname"]=="Windows",".dll",".so"),sep="")

PhenoRice <- function(evi, ndfi, noise, evismoothopt=FALSE, quiet=TRUE){
	TT <- length(evi)
	
	if (evismoothopt) {
		#dyn.load(system.file(smoothlib, package="PhenoRiceMap"))

		evismoothing <- function(evi,noise) {
			pesi <- vector(length=length(evi))
			pesi[which(noise==3)] <- 1
			pesi[which(noise==2)] <- 0.5
			pesi[which(noise==1)] <- 0.00001
			out <- .C("smoothTS", y=as.double(evi), w=as.double(pesi))
			return(out$y)
		}
		
		evi <- evismoothing(evi=evi,noise=noise)
		
		#dyn.unload(system.file(smoothlib, package="PhenoRiceMap"))
	} else {
		
		# compute EVI smoothing using R functions
		
		##############################################################
		##############################################################
		# Smoothing 
		##############################################################
		
		##############################################################
		# Spike detection
		##############################################################
		#----------------------------------------------------------------
		#Polynom for spike filling:
		
		liscia1 <- function(x,griglia=c(1,2,3,5,6,7)){
			if(length(x)!=7){
				stop("Hai sbagliato in fun liscia1")
			}
			x <- x[-4]
			coeff <- lm(x~griglia+I(griglia^2))$coefficients  #poliniomio per interpolazione spike
			#print(coeff)
			return(coeff[1]+coeff[2]*4+coeff[3]*4^2  )
		}
		#----------------------------------------------------------------
		
		#window for spike detection:
		w <- 2
		#posiz <- vector("numeric")
		
		#search the spikes
		for( i in 4:(TT-4)){                  #starts from 4th because: w=2
			app <- evi[ ((i-w):(i+w))[-(w+1)] ] #take the window +/- 2 without the central position
			med <- mean(app)  #calculate median
			sta <- sd(app)  #calculate st.dev.
			if(abs(evi[i]-med)>2*sta){  #condition to have a spike (if evi != mediana+2*st)
				#abs because the condition is ? +/- (spike pos e neg)
				
				evi[i] <- liscia1(evi[(i-3):(i+3)]) #in case of spike detection, fill with polynom 
			}                                         #defined in a window = 3
			
		}
		
		#evi2 = output from spike detection
		#evi2 <- evi
		
		##############################################################
		# Small drops
		##############################################################
		for(i in 2:(TT-2)){
			if( ((evi[i-1]-evi[i])>0) & ( (evi[i+1]-evi[i])>0) ){
				evi[i] <- (evi[i-1]+evi[i+1])/2
			}
		}
		#evi3 <- evi
		
		#repate the small drops filter
		for(i in 2:(TT-2)){
			if( ((evi[i-1]-evi[i])>0) & ( (evi[i+1]-evi[i])>0) ){
				evi[i] <- (evi[i-1]+evi[i+1])/2
			}
		}
		
		#evi4 output from small drops filter
		#evi4 <- evi
		
		
		##############################################################
		# weighted smoothing
		##############################################################
		#weights for smoothing, from "noise" input
		pesi <- vector(length=TT)
		pesi[which(noise==3)] <- 1
		pesi[which(noise==2)] <- 0.5
		pesi[which(noise==1)] <- 0.00001
		
		#----------------------------------------------------------------
		#Polynom for weighted smoothing:
		liscia2 <- function(y,griglia=c(-3,-2,-1,0,1,2,3),weight){
			if(length(y)!=7){
				stop("Hai sbagliato in fun liscia2")
			}
			coeff <- lm(y~griglia+I(griglia^2),weights=weight)$coefficients 
			return(coeff[1]+coeff[2]*griglia[4]+coeff[3]*griglia[4]^2)
		}
		#----------------------------------------------------------------
		
		evi_new <- vector(length=TT) #build vector with length TT
		for(i in 4:(TT-3)){           #weighted smoothing starts from 4th till fourth from the last
			gri <- (i-3):(i+3)        #window for weighted smoothing
			evi_new[i] <- liscia2(y=evi[gri],weight=pesi[gri])
		}
		
		#evi output from weighted smoothing
		evi <- evi_new
		
	}
    
    
  ##############################################################
  ##############################################################
  # Signal analysis
  ##############################################################
    
  ##############################################################
  # derivative analysis
  ##############################################################
  #cut the un-smoothed date (first three and last three)
    evi5 <- evi
    evi <- evi[4:(TT-3)]
    noise <- noise[4:(TT-3)]
    ndfi <- ndfi[4:(TT-3)]
    TT <- length(evi)
  
  # skip pixel if evi max doesn't reach value of 2000
  ricemap <- as.numeric(max(evi)>2000)
  uscita <- data.frame(ricemap=rep(0,3),season=0,riceseason=0,riceseason1=0,riceseason2=0,riceseason3=0,ricedet=0,min=0,startseason=0,max=0,endseason=0)
  uscita_evi <- data.frame(mini_evi=rep(0,3),massi_evi=0,delta_evi=0)
  if(ricemap==0){
    uscita$ricemap[1] <- ricemap
  } else {
    
    #Derivate: calculate relatives max and min after smoothing process
    deriva <- vector(length=TT)
    
    deriva[1] <- -3/2*evi[1]+2*evi[2]-1/2*evi[3] #derivate of first data
    deriva[TT] <- 1/2*evi[TT-2]-2*evi[TT-1]+3/2*evi[TT] #derivate of last data
    deriva[2:(TT-1)] <- diff(evi,lag=2)/2 #doesn't take first and last data because lag = 2
    
    #if the multiplication is negative there is a min or a max
    #i.e. where "minus times plus = minus"
    #                   shift1      shift2
    critici <- which(deriva[2:TT]*deriva[1:(TT-1)]<0) 
    quanti <- length(critici)  #count if there are min or max 
    
    #create empty vectors for min and max, length TT
    massi <- rep(0,TT) 
    mini <- rep(0,TT)
    
    #if there aren't max or min skip
    if(quanti==0){
      #print("no max or min")
      next
    } else {
      #for every min and max (i.e. vector quanti)
      #fill the vectors for min and max
    
      for(i in 1:quanti){
        if(deriva[critici[i]]>0){massi[critici[i]+1]=1 }
        if(deriva[critici[i]]<0){mini[critici[i]+1]=1 }
      }
    }
    
  ##############################################################
  ##############################################################
  # Rice phenology detection
  ##############################################################
  
  ##############################################################
  #Crop max detection 
  ##############################################################
    loc_massi <- which(massi==1)
    for(i in loc_massi){       
      #we consider the Crop max (heading) if occour in the solar year (jan-dec)
      if(i<19 | i> (TT-6)){ #cannot apply criteria aoutside this interval
        massi[i] <- 0                                                   
      } else{
        app_giu <- sum(deriva[(i-5):i]>0)<3 #1?criteria= 3 positive derivate before max (in window of 5)
        app_su <- sum(deriva[i:(i+5)]<0)<3  #1?criteria= 3 negative derivate after max (in window of 5)
        app_tre <- min(evi[i:(i+6)])>= (2/3*evi[i]) #2?criteria= EVI decrease of 1/3 after max (in window of 6)
        app_quattro <- evi[i]<2000 #3?criteria= Crop max have EVI > 2000
        if( app_su | app_giu  | app_tre | app_quattro ){ #all conditions must be verify
          massi[i] <- 0
        }
      }
    }
    
    #loc_massi have the detected Crop max
    loc_massi <- which(massi==1)
    
    
    ##############################################################
    #Veg min detection  
    ##############################################################
    
	loc_mini <- which(mini==1)
  
    for(i in loc_mini){
      if(i<3 | i> (TT-17)){
        mini[i] <- 0
      } else{
        app_su <- sum(deriva[i:(i+5)]>0)<3 #1?criteria= 3 positive derivate after min (in window of 5)
        terzo_criterio <- sum(massi[(i+6):(i+17)])==0 #3? criteria= after a min there is a max (in window of +6:+17)
        
        if(app_su | terzo_criterio){ #all conditions must be verified
          mini[i]=0
        }                
      }
    }
    
    #loc_mini have the detected Rice min
    loc_mini <- which(mini==1)
    loc_minipos <- which(mini>=0)
    loc_minipos <- loc_minipos*mini
    loc_minipos <- replace(loc_minipos, which(loc_minipos==0), 999)
    
    
   ##############################################################
   #Clean absolute min and max   
   ############################################################## 
    for (i in loc_massi) {
      if(sum(mini[(i-17):(i-6)])>1) {mini[min(loc_minipos[(i-17):(i-6)])]=0} # remove first Rice min if there are two consecutive min before each Crop max
      if(sum(mini[(i-17):(i-6)])<1) {massi[i]=0} #remove Crop max if there is not a Rice min in the window max-17:max-6
    }
    
    # clean Crop  max
    for (i in loc_massi) {
      if(sum(mini[(i-17):(i-6)])<1) {massi[i]=0} #remove Crop max if there is not a Rice min in the window max-17:max-6
    }
    
    #loc_mini e loc_massi have the detected Rice seasons
    loc_mini <- which(mini==1)
    loc_massi <- which(massi==1)
    
    #create vector with EVI values in correspondance of max
    loc_massipos <- which(massi>=0)
    loc_massipos <- loc_massipos*massi
    evi_massi <- evi*massi
    loc_evi_massi <- which(evi_massi>0)
    
    
    #if there are more than one max after each min remove the one with lower evi value
    for (i in loc_mini) {
      if(sum(massi[(i+6):(i+17)])>1) {massi[which(evi_massi==min(evi_massi[loc_massipos[(i+6):(i+17)]]))]=0}
    } 
    
    # update loc_massi
    loc_massi <- which(massi==1)
    
    # repeat cleaning filter for mini, if there are two consecutive loc_mini
    # create mini variables for cleaning
    loc_minipos <- which(mini>=0)
    loc_minipos <- loc_minipos*mini
    #loc_miniposfin <- loc_minipos
    loc_minipos <- replace(loc_minipos, which(loc_minipos==0), 999)
    # cleaning mini
    for (i in loc_massi) { #massi[which(evi_massi==min(evi_massi[loc_massipos[(i+6):(i+17)]]))]=0
      if(sum(mini[(i-17):(i-6)])>1) {mini[min(loc_minipos[(i-17):(i-6)])]=0} # remove first absolute min if there are two consecutive min before each max
    }
    
    # update loc_mini
    loc_mini <- which(mini==1)
    
	if (length(loc_mini)!=length(loc_massi)) {
		#print(paste("pixel ", p, ": have a different number of min and max")) # for debug purpose
		uscita$ricemap[1] <- 0
		mini[1:length(mini)] <- 0
		massi[1:length(massi)] <- 0
	}
	
	#loc_miniveg <- which(mini==1)
	#loc_massiveg <- which(massi==1)
	
	##############################################################
	#Rice min detection  
	##############################################################
	#reject flood condition in case of noise data
	ndfi[which(noise<2)]=0
	
	mini_rice <- mini
	loc_mini_rice <- which(mini==1)
	
	#loc_massi_rice <- which(massi==1)
	
	for(i in loc_mini_rice){
		if(i<3 | i> (TT-17)){
			mini_rice[i] <- 0
		} else {
			seconda_riga <- ndfi[(i-2):(i+3)]>0 #2? criteria= flood condition in a time window of -2 +3
			secondo_criterio <- sum(seconda_riga)==0 
			terzo_criterio <- sum(massi[(i+6):(i+17)])==0 #3? criteria= after a min there is a max (in window of +6:+17)
			
			if(secondo_criterio | terzo_criterio){ #all conditions must be verified
				mini_rice[i]=0
			}                
		}
	}
	
	#loc_mini have the detected Rice min
	loc_mini_rice <- which(mini_rice==1)
	
	# compute if detected seasons are rice seasons (preceded by flood event)
	isrice <- as.numeric((as.numeric(loc_mini %in% loc_mini_rice))==1)
	
	if(sum(isrice)>0){
		# compute which detected seasons are rice seasons (preceded by flood event)
		ricedetseason <- paste(which((as.numeric(loc_mini %in% loc_mini_rice))==1))
		ricedet <- paste(which((as.numeric(loc_mini %in% loc_mini_rice))==1), collapse="")
	} else {
		ricedet <- 0
		ricedetseason <- 0
	}
	
	# compute ricemap for the 3 seasons
	if(sum(as.numeric(ricedetseason))==0){
		riceseason1 <- 0
		riceseason2 <- 0
		riceseason3 <- 0
	} else {
		if(sum(as.numeric(ricedetseason %in% 1))>0){
			riceseason1 <- 1
		} else {
			riceseason1 <- 0
		}
		if(sum(as.numeric(ricedetseason %in% 2))>0){
			riceseason2 <- 1
		} else {
			riceseason2 <- 0
		}
		if(sum(as.numeric(ricedetseason %in% 3))>0){
			riceseason3 <- 1
		} else {
			riceseason3 <- 0
		}
	}
	
    ##############################################################
    #calculate SoS EoS and save outputs
    ##############################################################
    #print results into output dataframe
    ricemap <- as.numeric(sum(massi)>0)

	vegmap <- as.numeric(sum(massi)>0)
	
	if(length(loc_mini)!=length(loc_massi)){
		ricemap <- 0
		vegmap <- 0
	}
	
    uscita$ricemap[1] <- ricemap
    
	if (vegmap==1){
      
      nvegseason <- length(loc_massi)
	  nriceseason <- length(loc_mini_rice)
	  ricedetected <- ricedet
	  
      if(nvegseason>3){
        uscita$ricemap[1] <- ricemap
		#nriceseason <- nriceseason
        nvegseason <- 3
        #if(!quiet) message("Computing phenorice found more than 3 seasons: ", tile, " pixel: ", p)
      }else{
      
      # set ricemap to zero if min number differs from max number
      if (length(loc_mini)!=length(loc_massi)) {
        #print(paste("pixel ", p, ": have a different number of min and max")) # for debug purpose
        uscita$ricemap[1] <- 0
      } else {
        #
        # calculate evi delta between min and max
        #delta_evicheck <- vector('numeric')
        #for (i in 1:nriceseason){
        #  delta_evicheck[i] <- evi[loc_massi[i]]-evi[loc_mini[i]]
        #}
        #
        ## remove min and max if evi delta is lower than specified delta (e.g 1400)
        #
        #delta <- 1400
        #
        #for (i in 1:nriceseason){
        #  if(delta_evicheck[i]<delta) {
        #    loc_mini <- loc_mini[-i]
        #    loc_massi <- loc_massi[-i]
        #  }
        #}
        
        ## update number of rice seasons
        #nriceseason <- length(loc_massi)
        #three is the maximum number of season per year  
        #nriceseason <- ifelse(nriceseason>3,3,nriceseason)
        
        #if (nriceseason==0){
        #  uscita$ricemap[1] <- 0
        #} else {
          # calculate SoS and EoS
          inizio_stag <- vector('numeric')
          fine_stag <- vector('numeric')
          app <- evi[loc_mini]+0.1*evi[loc_mini]
          app <- evi[loc_mini]+0.1*(evi[loc_massi]-evi[loc_mini]) #Sos= 10% of difference between Rice min and Rice max
          app1 <- evi[loc_massi]-1/2*evi[loc_massi] #EoS= 50% of difference between Rice min and Rice max  
          
          
          for(i in 1:nvegseason){
            inizio_stag[i] <- which(evi[loc_mini[i]:TT]>app[i])[1] 
            inizio_stag[i] <- loc_mini[i]+ inizio_stag[i]-1
            
            fine_stag[i] <- which(evi[loc_massi[i]:TT]<app1[i])[1] 
            fine_stag[i] <- loc_massi[i]+ fine_stag[i]-1
          }
          
		  uscita$season[1] <- nvegseason
		  uscita$riceseason[1] <- nriceseason
		  uscita$ricedet[1] <- ricedetected
		  uscita$riceseason1[1] <- riceseason1
		  uscita$riceseason2[1] <- riceseason2
		  uscita$riceseason3[1] <- riceseason3
		  
          for (i in 1:nvegseason){            
                  
            uscita$min[i] <- c(loc_mini[i]+3)
            uscita$startseason[i] <- c(inizio_stag[i]+3)
            uscita$max[i] <- c(loc_massi[i]+3)
            uscita$endseason[i] <- c(fine_stag[i]+3)
          }  
          
          if (evismoothopt) {
            # calculate flood evi and peak evi        
            mini_evi <- vector('numeric')
            massi_evi <- vector('numeric')
            delta_evi <- vector('numeric')
            for (i in 1:nvegseason){
              mini_evi[i] <- evi[loc_mini[i]]
              massi_evi[i] <- evi[loc_massi[i]]
              delta_evi[i] <- massi_evi[i]-mini_evi[i]
            }
            # update uscita_evi with results
            
            for (i in 1:nvegseason){
              uscita_evi$mini_evi[i] <- c(mini_evi[i])
              uscita_evi$massi_evi[i] <- c(massi_evi[i])
              uscita_evi$delta_evi[i] <- c(delta_evi[i])
            }
          }
        #}
      }
    }
  }
 } 
  gc(verbose=FALSE, reset=TRUE)
  if (evismoothopt) {
    return(list(uscita,uscita_evi,evi5))
  } else {
    return(list(uscita))
  }
}
