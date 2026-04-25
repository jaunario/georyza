# Author:(translator) Jorrel Khalil S. Aunario, j.aunario@irri.org
# Date :  25 February 2014
# Version 0.0.1
# Licence GPL v3
###############################################################################


hourly.srad <- function(doy,tod,lat,tavg,rh,srad,I00=1365,beta_dust=0.1,pre=1013){
	
	delta <- -23.4 * cos(360 * (doy + 10) / 365 * pi / 180) * pi / 180              # the angle between the sun's ray and the equatorial plane of the earth, solar declination (radian)
	Hour  <-  15 * (tod - 12) * pi / 180                                            # the hour angle of the sun (radian)
	lat_rad <- lat * pi / 180                                                       # latitude (radian)
	cos_theta <- sin(lat_rad) * sin(delta) + cos(lat_rad) * cos(delta) * cos(Hour)  # cosine of theta, the zenith angle of the sun
	
	# Daily mean of solar radiation at the top of the atmosphere, S0d (W/m2)
	
	eta <- 2 * pi / 365 * doy
	dd2 <- 1.00011 + 0.034221 * cos(eta) + 0.00128 * sin(eta) + 0.000719 * cos(2 * eta) + 0.000077 * sin(2 * eta)
	x <-  -tan(lat_rad) * tan(delta)
	H <-  atan(-x / sqrt(-x * x + 1)) + 2 * atan(1)
	S0d <- I00 / pi * dd2 * (H * sin(lat_rad) * sin(delta) + cos(lat_rad) * cos(delta) * sin(H))
	
	# Daily mean of solar radiation at the ground surface on a clear day, Sdf (W/m2)
	C1 <-  0.21 - 0.2 * beta_dust
	mNoon <- (Pre / 1013) * 1 / cos(lat_rad - delta)
	#k3 <- 1.402 - 0.06 * log(beta_dust + 0.02) / log(10#) - 0.1 * (1 / Cos(lat_rad - delta) - 0.91) ^ 0.5
	xx <- log10(beta_dust + 0.02)
	k3 <- 1.402 - 0.06 * log10(beta_dust + 0.02) - 0.1 * (1 / cos(lat_rad - delta) - 0.91) ^ 0.5
	md <- k3 * mNoon
	F1 <- 0.056 + 0.16 * beta_dust ^ 0.5
	e_hPa_mean <- svp(tavg) * rh / 100 # Daily mean of water vapor pressure (hPa) calculated from tavg and RH_mean
	TDEW <- 237.3 * log10(e_hPa_mean / 6.11) / (7.5 - log10(e_hPa_mean / 6.11))       # Dew point temperature (C)
	log10w <- 0.0312 * TDEW - 0.0963
	i3 <- 0.014 * (md + 7 + 2 * log10w) * log10w
	ref <- 0.2   #Albedo of the ground surface
	j1 <- (0.066 + 0.34 * beta_dust ^ 0.5) * (ref - 0.15)
	Sdf <- S0d * (C1 + 0.7 * 10 ^ (-md * F1)) * (1 - i3) * (1 + j1)
	
	
	#Ratio of observed daily mean solar radiation to Sdf on a clear day, Sd_ratio
	Sd_ratio <- Sd_mean / Sdf
	
	#Instantaneous solar radiation at the top of the atmosphere, S0 (W/m2)
	S0 <- I00 * dd2 * cos_theta
	
	if (cos_theta <= 0) {
		Rsd <- 0 
	} else{
		#Instantaneous solar radiation at the ground surface on a clear day, Rsf (W/m2)
		m <- (Pre / 1013) * 1 / cos_theta 
		i1 <- 0.014 * (m + 7 + 2 * log10w) * log10w
		Rsf <- S0 * (C1 + 0.7 * 10 ^ (-m * F1)) * (1 - i1) * (1 + j1)
		
		#Instantaneous solar radiation at the ground surface, RsdiW/m2)
		Rsd <- Rsf * Sd_ratio        
	}
	return(Rsd)    
}

hourly.temp <- function(tavg, tmin, tmax, tod,m1=2.09,m2=0.2,m3=45-pi/180){
	#Instantaneous air temperature, TaiC)
	return(tavg + (tmax - tmin) / m1 * (-cos(pi / 12 * tod - m3) + m2 * cos(pi / 6 * tod - m3)))
}

hourly.wind <- function(wind,tod,m1=2.77,m2=0.17,m3=49.5*pi/180){
	
	wind_max <- 2.2 * wind #Daily maximum of wind speed (m/s) empirically estimated from U_mean
	
	#Instantaneous wind speed, Uim/s)
	return(wind + (wind_max - 0) / m1 * (-cos(pi / 12 * tod - m3) + m2 * cos(pi / 6 * tod - m3)))
}

