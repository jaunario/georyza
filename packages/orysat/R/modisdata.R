# Authors: Jorrel Khalil S. Aunario
# Date: 18 May 2011

if (!isGeneric("modis.data")){
	setGeneric("modis.data", function(x,...) standardGeneric("modis.data"))
}

#setMethod('modis.data', signature(x='missing'), 
#	function(nrows=2400, ncols=2400, xmn=-180, xmx=180, ymn=-90, ymx=90, crs) {
#		e <- extent(xmn, xmx, ymn, ymx)
#		if (missing(crs)) {
#			if (e@xmin > -360.1 & e@xmax < 360.1 & e@ymin > -90.1 & e@ymax < 90.1) { 
#				crs ="+proj=longlat +datum=WGS84"
#			} else {
#				crs=NA
#			}
#		}
#		r <- modis.data(e, nrows=nrows, ncols=ncols, crs=crs)
#		return(r)
#	}
#)

setMethod("modis.data", signature(x='missing'), 
		function(x, ...) {
			return(new("modis.data",...))
		}
)

setMethod("modis.data", signature(x="modis.data"),
		function(x){
			m <- new("modis.data", product=x@product, acqdate=x@acqdate, zone=x@zone, version=x@version, 
					proddate=x@proddate, projection=x@projection, extent=x@extent, ncols=x@ncols, nrows=x@nrows, imgvals=data.frame())
			return(m)			
		}

)

setMethod("modis.data", signature(x="data.frame"),
		function(x, cached=FALSE, cache.path="./cache"){
			mstack <- stack(x$filename)
			NAvalue(mstack) <- -26872
			
			m <- new("modis.data", product=x$product[1], acqdate=x$acqdate[1], zone=x$zone[1], version=x$version[1], 
					proddate=x$proddate, projection=projection(mstack), extent=extent(mstack), ncols=ncol(mstack), nrows=nrow(mstack), imgvals=data.frame())
			
			if (cached){
				if(!file.exists(cache.path)){
					dir.create(cache.path, recursive=TRUE)	
				} 
				st <- Sys.time()
				m@imgvals <- data.frame(row=numeric(0),cache.file=character(0))
				for (i in 1:2400){
					message("reading row ", i, "\r", appendLF=FALSE)
					tocache <- as.data.frame(values(mstack,i))
					save(tocache,file=paste(cache.path,"/",deparse(substitute(x)),"_",i,".RData",sep=""))
					m@imgvals[i,] <- NA
					m@imgvals$row[i] <- i
					m@imgvals$cache.file[i] <- paste(cache.path,"/",deparse(substitute(x)),"_",i,".RData",sep="")
				}				
				Sys.time()-st
			} else {
				m@imgvals <- as.data.frame(values(mstack))
			}
			return(m)			
		}

)

setMethod("modis.data", signature(x="RasterStack"),
		function(x){
			m <- new("modis.data", product="", acqdate="", zone="", version="", 
					proddate="", projection=projection(x), extent=extent(x), ncols=ncol(x), nrows=nrow(x), imgvals=data.frame(values(x)))
			return(m)			
		}

)

setMethod("modis.data", signature(x="character"),
		function(x, nodata=FALSE){
			info <- modisFiles(path=dirname(x), pattern=basename(x))
			if (length(x)>1){
				base <- stack(x)
			} else {
				base <- raster(x)
			}
			m <- new("modis.data", product=info$product[1], acqdate=info$acqdate[1], zone=info$zone[1], version=info$version[1], 
					proddate=info$proddate[1], projection=projection(base), extent=extent(base), ncols=ncol(base), nrows=nrow(base),imgvals=data.frame())
			if(!nodata) {
				m@imgvals <- data.frame(values(base))
				colnames(m@imgvals) <- info$band
			}	
			rm(base)
			gc(verbose=FALSE,reset=TRUE)
			return(m)
		}		
)

