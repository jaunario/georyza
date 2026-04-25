
smooth.raster <- function(objraster, maskraster=NULL, rm.zeros=TRUE){
	
	if(!is.null(maskraster)){
	  rs1 <- disaggregate(objraster, fact=res(objraster)[1]/res(maskraster)[1])  
	  rs1 <- mask(rs1, maskraster)
	} else {
	  rs1 <- objraster
	}
	rs1 <- focal(rs1, fun=mean, w=matrix(1,5,5), na.rm=TRUE)
	rs1 <- focal(rs1, fun=mean, w=matrix(1,5,5), na.rm=TRUE)
	if(!is.null(maskraster)) rs1 <- mask(rs1, maskraster)
	if(rm.zeros) rs1[values(rs1)<=0] <- NA 
	rs1 <- disaggregate(rs1, fact=3)
	return(rs1)    
}

stress.map <- function(datamatrix, baseraster, maskraster, writeto=getwd(), verbose=TRUE){
	custom.colors <- colorRampPalette(c('#330066', '#3300AA', '#33EEAA', '#66FF66', '#99FF22', '#FFFF00', '#FF8844', '#CC5500', '#CC1100', '#990000'), interpolate='spline', space='rgb')
	
	if (!force.directories(writeto, recursive=TRUE)) stop("Cannot create output directory.")
	grps <- vector()
	grprange <- vector()
	thisgrp <- ""
	sditems <- grep("sd", colnames(datamatrix))                 
	for (i in 1:ncol(datamatrix)){
		key <- unlist(strsplit(colnames(datamatrix)[i],"_"))
		
		if (verbose) show.message("Writing ", colnames(datamatrix)[i], ".tif to disk.", eol="\r")        
		objraster <- raster(baseraster)
		objraster[as.numeric(rownames(datamatrix))] <- datamatrix[,i]
		dttype <- ifelse(grepl("Cold_",colnames(datamatrix)[i])|grepl("ColdYr_",colnames(datamatrix)[i]),"INT1U","FLT4S")    
		objraster <- writeRaster(objraster, filename=paste(writeto,"/",colnames(datamatrix)[i],".tif",sep=""), format="GTiff", datatype=dttype, overwrite=TRUE)
		
		if (length(grep("diff",colnames(datamatrix)[i]))>0) next
		#if (verbose) show.message("Disaggregating and resampling ", colnames(datamatrix)[i], eol="\r")
		
		# get group also
		#searchstr <- ifelse(key[length(key)]=="sd",paste(key[1],"[[:punct:]][[:graph:]]*sd",sep=""),paste(key[1],"[[:punct:]][[:graph:]]*",sep=""))
		#if (key[length(key)]!="sd") {
		#objraster <- focalDisagg(objraster, maskraster)
		#grp <- setdiff(grep(searchstr, colnames(datamatrix)),sditems)
		#thisgrp <- key[1]         
		#} else {
		#objraster <- focalDisagg(objraster, maskraster, rm.zeros=FALSE)
		#grp <- grep(searchstr, colnames(datamatrix))
		#thisgrp <- paste(key[1], "sd", sep="_")
		#}
		#okgrp <- which(grps==thisgrp)
		#if (length(okgrp)!=1){            
		#	grps <- c(grps,thisgrp)
		
		#	grpmin <- floor(min(datamatrix[,grp],na.rm=TRUE))          
		#		grpmax <- ceiling(max(datamatrix[,grp],na.rm=TRUE))
		#		grprange <- rbind(grprange,c(grpmin,grpmax))
		
		# for legend
		#		rs <- raster(baseraster)
		#		rs[as.numeric(rownames(datamatrix))] <- sort(runif(nrow(datamatrix),min=grpmin, max=grpmax))
		#		png(filename = paste(writeto,paste(thisgrp,"_legend.png",sep=""),sep="/"), width=600, height=600, res=100)
		#		plot(rs, col = custom.colors(255))
		#		text(4,4,"xxxx")
		#		mtext("Legend",4,padj=-8.5,line=1.5,las=1,cex=1.2)
		#		if(key[1] %in% c("tmin", "tavg", "tmax")) mtext("?C",4,padj=-6.5,line=4.1,las=1) else mtext("days",4,padj=-6.5,line=4.1,las=1)
		#		dev.off()
		
		#	} else {
		#		grpmin <- grprange[okgrp,1]
		#		grpmax <- grprange[okgrp,2]
		#	}
		
		#	show.message("Writing ", colnames(datamatrix)[i], ".kml to disk.", eol="\r")
		#	KML(objraster, filename=paste(writeto,"/",colnames(datamatrix)[i],".kml",sep=""), col = custom.colors(255), maxpixels=ncell(objraster), zip=FALSE, zlim=c(grpmin,grpmax))
		show.message(colnames(datamatrix)[i], " done.", eol="\n")
	}    
	return(writeto) 
}
