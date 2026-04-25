# TODO: See NOAA Solar Calculations
nooa.params <- function(date, timezone, time=0.1){
	date <- as.numeric(date)-as.numeric(as.Date("1900-1-1"))+2
	julianday <- date+2415018.5+(time-timezone)/24
	juliancentury <- (julianday-2451545)/36525
	
	GMLS <- (280.46646+juliancentury*(36000.76983 + juliancentury*0.0003032)) %% 360
	GMAS <- 357.52911+juliancentury*(35999.05029 - 0.0001537*juliancentury)
	EEO <- 0.016708634-juliancentury*(0.000042037+0.0000001267*juliancentury)
	
	SEC <- sin(radians(GMAS))*(1.914602-juliancentury*(0.004817+0.000014*juliancentury))+sin(radians(2*GMAS))*(0.019993-0.000101*juliancentury)+sin(radians(3*GMAS))*0.000289
	
	STL <- GMLS+SEC
	STA <- GMAS+SEC
	#SRA <- degrees(atan2(cos(radians(SAL)),cos(radians(OC))*sin(radians(SAL))))
	
	MOE <- 23+(26+((21.448-juliancentury*(46.815+juliancentury*(0.00059-juliancentury*0.001813))))/60)/60
	OC <- MOE+0.00256*cos(radians(125.04-1934.136*juliancentury))
	SAL <- STL-0.00569-0.00478*sin(radians(125.04-1934.136*juliancentury))
	SD <- degrees(asin(sin(radians(OC))*sin(radians(SAL))))
	
	var.y <- tan(radians(OC/2))*tan(radians(OC/2))
	eq.of.time <- 4*degrees(var.y*sin(2*radians(GMLS))-2*EEO*sin(radians(GMAS))+4*EEO*var.y*sin(radians(GMAS))*cos(2*radians(GMLS))-0.5*var.y^2*sin(4*radians(GMLS))-1.25*EEO^2*sin(2*radians(GMAS)))
	
	return(list(sun.dec=SD,eq.t=eq.of.time))
}


# TODO: See NOAA Solar Calculations for now is 12
solar.time <- function(lat, lon, ...){
	tz <- lon%/%15
	nooa.vals <- nooa.params(..., timezone=tz)
	ha.sunrise <- degrees(acos(cos(radians(90.833))/(cos(radians(lat))*cos(radians(nooa.vals$sun.dec)))-tan(radians(lat))*tan(radians(nooa.vals$sun.dec))))
	solar.noon <- ((720-4*lon-nooa.vals$eq.t+tz*60)/1440)
	sunrise <- solar.noon-ha.sunrise*4/1440
	sunset <- solar.noon+ha.sunrise*4/1440
	sundur <- ha.sunrise*8
	return(list(noon=solar.noon*24,rise=sunrise*24,set=sunset*24,duration=sundur/60))
}


# TODO Also use NOAA Solar Calculations
t.sunrise <- function(DL, LSH=solar.time()$solar.noon){
	return(LSH-0.5*DL)
}

t.sunset <- function(DL, LSH=solar.time()$solar.noon){
	return(LSH+0.5*DL)
}

t.peaktemp <- function(P, LSH=solar.time()$solar.noon){
	return(LSH+P)
}

p.fltime <- function(t, PRDEL=0.1, mn=16.5, stdev=1.71969045571797){
	return( pnorm(t+0.5*PRDEL, mn, stdev)-pnorm(t-0.5*PRDEL, mn, stdev))
}