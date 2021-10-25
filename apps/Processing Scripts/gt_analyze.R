# RUN THIS FOR TESTING ALGORITHMS

shp.points <- shapefile("C:/WORKSPACE/SPIA/Validation/BGD_BrahmaputraRegion_2017_AMAN_RnR_SRSP_VGrid.shp")
shp.sin <- spTransform(shp.points, projection(rst.lc))
cells.gt <- cellFromXY(rst.lc, rgeos::gCentroid(shp.sin, byid = TRUE))
idx.ptp <- which(cells.gt$cells.gt %in% pixels.toprocess)

cells.gt <- data.frame(cells.gt, rgeos::gCentroid(shp.sin, byid = TRUE), shp.sin$RnR)
cells.gt <- cells.gt[idx.ptp,]
smevi.gt <- t(smooth.evi[,match(cells.gt$cells.gt, pixels.toprocess)])



gt.idx.maxevi <- apply(smevi.gt, 1, which.max)

# Remove pixels with (idx.maxevi-13) < 0
psv <- which(gt.idx.maxevi<13)

# Analyze Maturity Signature... Consistent decline
post.maxevi <- mapply('[', as.data.frame(smooth.evi), as.data.frame(sapply(idx.maxevi, seq, length.out=10)))
post.diff <- apply(post.maxevi, 2, diff)
post.sign <- apply(post.diff, 2, sign)

post.consdesc <- which(colSums(post.sign[1:4,]) == -4)
post.wcmin <- apply(post.diff, 2, which.min)



write.csv(cells.gt, paste0("SMEVI_2017_h25v06_", format(min(date.acqdates), "%x"),"_",format(max(date.acqdates), "%x"),".csv"), row.names = FALSE)
