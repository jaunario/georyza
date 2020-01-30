# Title: IRRI-GIS Implementation of Xiao (2006) Rice Mapping Algorithm

#     -4 - Forest/Thick Vegetation 
#     -5 - Shrub
#     -6 - persistent water

cell.ref <- data.frame(cell=pixels.cropland, rice=NA)


filelist.index <- dir(INDICES_DIR, pattern = paste0(PRODUCTS, ".", TILE, ".*.tif$"), recursive = TRUE, full.names = TRUE)
inv.idxfiles <-  inventory.modis(filelist.index, modisinfo = c("product", "zone", "acqdate", "band"), file.ext = "tif")
inv.idxfiles <- inv.idxfiles[order(inv.idxfiles$band, inv.idxfiles$acqdate),]

stk.jday <- raster::stack(inv.idxfiles$filename[inv.idxfiles$band=="julianday"])
stk.evi <- raster::stack(inv.idxfiles$filename[inv.idxfiles$band=="EVI"])
stk.ndvi <- raster::stack(inv.idxfiles$filename[inv.idxfiles$band=="NDVI"])
stk.lswi <- raster::stack(inv.idxfiles$filename[inv.idxfiles$band=="LSWI"])

# Extract data on 
mat.lswi <- stk.lswi[cell.ref$cell]/10000
mat.evi <- stk.evi[cell.ref$cell]/10000
mat.ndvi <- t(stk.ndvi[cell.ref$cell]/10000)

mat.lswi <- t(mat.lswi)
mat.evi <- t(mat.evi)
mat.ndvi <- t(mat.ndvi)

# TODO: try Fill NAs using approx funtion

cell.ref$na.count <- apply(mat.evi, 2, function(...){return(sum(!is.na(...)))})
cell.ref$rice[apply(mat.ndvi,2,xiaoflags.forest)] <- -4
cell.ref$rice[apply(mat.lswi,2,xiaoflags.shrub)] <- -5
cell.ref$rice[mapply(xiaoflags.persistentwater, ndvi=as.data.frame(mat.ndvi), lswi=as.data.frame(mat.lswi))] <- -6

# Detect Persistent Water
st <- Sys.time()
fill.evi <- as.data.frame(apply(mat.evi[,cell.ref$na.count>=30], 2, approx, x=1:nrow(mat.evi), xout=1:nrow(mat.evi)))
fill.evi <- fill.evi[,grep("y",colnames(fill.evi))]
rm(mat.evi)


fill.ndvi <- as.data.frame(apply(mat.ndvi[,cell.ref$na.count>=30], 2, approx, x=1:nrow(mat.ndvi), xout=1:nrow(mat.ndvi)))
fill.ndvi <- fill.evi[,grep("y",colnames(fill.ndvi))]
rm(mat.evi)

fill.lswi <- as.data.frame(apply(mat.lswi[,cell.ref$na.count>=30], 2, approx, x=1:nrow(mat.lswi), xout=1:nrow(mat.lswi)))
fill.lswi <- fill.lswi[,grep("y",colnames(fill.lswi))]
rm(mat.lswi)
gc(reset)
en <- Sys.time()
en-st

st <- Sys.time()
for (i in 1:nrow(mat.ndvi)){
  if(!is.na(cell.ref$rice[i])) next
  
  pixel.data <- data.frame(evi=mat.evi[,i], ndvi=mat.ndvi[,i], lswi=mat.lswi[,i])
  if(sum(colSums(!is.na(pixel.data)))<6){
    message("cell-", cell.ref$cell[i], ": Not Enough data.") 
    cell.ref$rice[i] <- -99
    next
  }
  temp.fill <- data.frame(apply(pixel.data, 2, approx, x=1:nrow(pixel.data), xout=1:nrow(pixel.data)))
  temp.fill <- temp.fill[,seq(from=2, by=2, length.out = ncol(pixel.data))]
  colnames(temp.fill) <- colnames(pixel.data)
  pixel.data <- temp.fill
  rm(temp.fill)
  
  cell.ref$rice[i] <- rice.Xiao_v1(pixel.data,evi.ricemax = -1, evi.halfricemax = -1)
  rm(pixel.data)
  message("cell-", cell.ref$cell[i], ": Done.")
}
en <- Sys.time()
en-st


