#################################################
ob_path<-"../OUTPUT/20231204/Obs_NCUMG_day1_20231204_Thr40_Object1.nc"
model_path<-"../OUTPUT/20231204/20231204_Thr40_NCUMG_day1RotatedModel_1.nc"
rotated_model<-"../OUTPUT/20231204/20231204_Thr40_NCUMG_day1RotatedModel_1.nc"
plot_file_name<-"CRAplot_20231204_Th40_NCUMG_day1_1.png"
tbl1<-"../OUTPUT/20231204/20231204_Thr40_NCUMG_day1_1_Table1.txt"
tbl2<-"../OUTPUT/20231204/20231204_Thr40_NCUMG_day1_1_Table2.txt"
tbl3<-"../OUTPUT/20231204/20231204_Thr40_NCUMG_day1_1_Table3.txt"
tbl4<-"../OUTPUT/20231204/20231204_Thr40_NCUMG_day1_1_Table4.txt"
#################################################
library(raster)
library("RColorBrewer")
orf <- raster(ob_path)
mrf <- raster(model_path)
orf[orf<40]=NA;
mrf[mrf<40]=NA;
rotated<-raster(rotated_model)
crat=40 ######### Enter Rainfall threshold here (RF in mm/day units)
orft<-orf
mrft<-mrf
png(file=plot_file_name,width = 650, height = 650, units = "px", pointsize = 12, bg = "white",  res = NA)
scatter_data<-"../OUTPUT/20231204/20231204_Thr40_NCUMG_day1_1_ScatterData.txt"
#layout(matrix(c(1, 2, 3, 0,4,0),2,3, nrow = T))
#layout(matrix(c(1, 2, 3, 0, 4, 0), nrow = 3, byrow = TRUE))
#layout.show(n=6)
##################################################################
library(rgdal)
shp <- readOGR(dsn = file.path("../SHP/IND_WHOLE.shp"), stringsAsFactors = F)
mycols=c("darkturquoise","darkseagreen1","darkgreen","yellow","goldenrod1","violetred1","maroon4")
myrange=colorRampPalette(mycols)
mycrange=myrange(500)
mycols = mycrange[c(0:500)]

par(mfrow=c(3,2), family="serif", mai=c(0.4,0.4,0.3,0.3),cex.axis = 1.5,cex.main = 2,font = 2,cex.axis = 1.5)
## Obs plot
plota<-raster::plot(orf,
	xlim = c(70,90),
	ylim = c(5,35),
	zlim=c(40,500),
	col=mycols,
	xlab = "", 
	ylab = "",
	interpolate=TRUE,
	asp=1.00,
	legend=TRUE,
	main="Observation",
        xaxs="r",
        yaxs="r",
        font=2,
        cex.axis = 1.5)
#axis(1, font = 2)
#axis(2, font = 2)
title(ylab="Latitude",mgp=c(2,1,0),cex.lab=1.5,font=2)
title(xlab="Longitude",mgp=c(2,1,0),cex.lab=1.5,font=2)
plot(shp,add=TRUE)

## Scatter Model Vs Obs
par(mai=c(0.5,0.6,0.4,0.4))
pdata=read.delim(scatter_data, header = TRUE, sep = "")
pdata[pdata<40]=NA;
plot(x = pdata[,1],y = pdata[,2],
    xlim = c(40,500),
    ylim = c(40,500),
    ylab="",
    xlab="",
        xaxs="r",
        yaxs="r",
    cex.axis = 1.5,
    font=2,
   main = "Obs vs Model")
abline(coef = c(0,1),lwd = 2)
title(ylab="Model (Shifted+Rotated)",mgp=c(3,1,0),cex.lab=2,font=2)
title(xlab="Observation",mgp=c(2.5,1,0),cex.lab=2,font=2)

##Model 
par(mai=c(0.4,0.4,0.3,0.3))
plotb <-raster::plot(mrf,
	xlim = c(70,90),
	ylim = c(5,35),
	zlim=c(40,500),
	col=mycols,
	xlab = "", 
	ylab = "",
	interpolate=TRUE,
	xaxs="r",
	yaxs="r",
	asp=1.00,
	legend=FALSE,
	main="Model",
        font=2,
        cex.axis = 1.5)
axis(1, font = 2)
axis(2, font = 2)
title(ylab="Latitude",mgp=c(2,1,0),cex.lab=1.5,font=2)
title(xlab="Longitude",mgp=c(2,1,0),cex.lab=1.5,font=2)
plot(shp,add=TRUE)
											
plot.new()
##Shifted + Rotated Model 
par(mai=c(0.4,0.4,0.3,0.3))
plotc <-raster::plot(rotated,
	xlim = c(70,90),
	ylim = c(5,35),
	zlim=c(40,500),
	col=mycols,
	xlab = "", 
	ylab = "",
        xaxs="r",
        yaxs="r",
	interpolate=TRUE,
	asp=1.00,
	legend=FALSE,
	main="Model (Shifted + Rotated)",
        font=2,
        cex.axis = 1.5)
title(ylab="Latitude",mgp=c(2,1,0),cex.lab=1.5,font=2)
title(xlab="Longitude",mgp=c(2,1,0),cex.lab=1.5,font=2)
plot(shp,add=TRUE)

#########  Plotting the Tables ##############

detach("package:raster")
library(gridExtra)
library(grid)
library(formattable)

bg_params=list(fill=c(NA, rep("white",5)), col="white")
tt=ttheme_minimal(core=list(fg_params=list(hjust=0, x=0.10,fontsize = 13)),
                      rowhead=list(fg_params=list(hjust=0, x=0,fontsize = 13)))

rows<-c("CRA Threshold=40mm/day    ")
# Push viewport
pushViewport(viewport(x=0.590,y=0.620,height=0.8,width=1.0))
grid.table(rows,theme=tt)
testdf1<-read.table(tbl1)
# Push viewport
pushViewport(viewport(x=0.610,y=0.40,height=0.8,width=0.5))
grid.table(testdf1,theme=tt)
#
testdf2<-read.table(tbl2)
# Push viewport
pushViewport(viewport(x=0.600,y=0.25,height=1.0,width=1.0))
grid.table(testdf2,theme=tt)
#
rows2<-c(substitute(bold("      Error Decomposition:                       ")))
# Push viewport
pushViewport(viewport(x=0.23,y=0.38,height=0.5,width=1.00))
grid.table(rows2,theme=tt)
#
testdf3<-read.table(tbl3,header = FALSE,sep = ",")
rows3<-c(paste(testdf3[1,1],"%  "))
rows4<-c(paste(testdf3[2,1],"%  "))
rows5<-c(paste(testdf3[3,1],"%  "))
rows6<-c(paste(testdf3[4,1],"%  "))
ab=rbind(rows3,rows4,rows5,rows6)
# Push viewport
pushViewport(viewport(x=0.500,y=0.27,height=0.9,width=1.3))
grid.table(ab,theme=tt,rows=NULL)
#
testdf4<-read.table(tbl4)
rows7<-c(paste(round(testdf4[1,1],digits=1),",",round(testdf4[2,1],digits=1),",",round(testdf4[3,1],digits=1)))
bold_text <- textGrob("Location Error (x, y, Angle)=", gp = gpar(fontface = "bold"))
# Push viewport
pushViewport(viewport(x = 0.51, y = 0.230, height = 0.8, width = 0.10))
grid.draw(bold_text)
table_data <- data.frame(
  Lon = round(testdf4[1, 1], digits = 2),
  Lat = round(testdf4[2, 1], digits = 2),
  Angle = round(testdf4[3, 1], digits = 2)
)
# Push viewport
pushViewport(viewport(x = 3.49, y = 0.520, height = 4.5, width = 4.50))
grid.table(rows7, theme = tt)
##############
