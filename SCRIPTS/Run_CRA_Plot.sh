#!/bin/bash
module load gnu/R/3.4.3
##################################################
export date=$1
export ObjNum=$2
export case=$3
export OUT_PATH=$4
export Threshold=$5
echo "Received arguments: $@"
RUNPATH=`pwd`
OUT_PREFIX=${OUT_PATH}/${date}_Thr${Threshold}_${case}
echo $OUT_PREFIX
shpfile=../SHP/IND_WHOLE.shp
figname=CRAplot_${date}_Th${Threshold}_${case}_${ObjNum}.png
###################################################
INPUT_OBS=${RUNPATH}/Obs_${case}_${date}_Thr${Threshold}_Object${ObjNum}.nc
INPUT_MODEL=${RUNPATH}/Model_${case}_${date}_Thr${Threshold}_Object${ObjNum}.nc
INPUT_RotModel=${OUT_PREFIX}RotatedModel_${ObjNum}.nc

############ From Stat Files #################################
CRA_threshold="${Threshold}mm/day"
################# Obs Vs Actual ############# 
cd $RUNPATH
###
cat > ${date}_PlotCRA_Th${Threshold}_${case}_${ObjNum}.R << EOF
# Read raster data and isolate objects by thresholds
# Obtain the area of objects.
library(raster)
library("RColorBrewer")
#################################################
orf <- raster("${INPUT_OBS}")
mrf <- raster("${INPUT_MODEL}")
orf[orf<${Threshold}]=NA;
mrf[mrf<${Threshold}]=NA;
rotated<-raster("${INPUT_RotModel}")
crat=$Threshold ######### Enter Rainfall threshold here (RF in mm/day units)
orft<-orf
mrft<-mrf
png(file="${figname}",width = 650, height = 650, units = "px", pointsize = 12, bg = "white",  res = NA)

mcoords <- raster::xyFromCell(mrf, which(!is.na(values(mrf))))
mx_min <- min(mcoords[, 1])-5
mx_max <- max(mcoords[, 1])+5
my_min <- min(mcoords[, 2])-5
my_max <- max(mcoords[, 2])+5

ocoords <- raster::xyFromCell(orf, which(!is.na(values(orf))))
ox_min <- min(ocoords[, 1])-5
ox_max <- max(ocoords[, 1])+5
oy_min <- min(ocoords[, 2])-5
oy_max <- max(ocoords[, 2])+5

x_min<-min(mx_min,ox_min)
x_max<-max(mx_max,ox_max)
y_min<-min(my_min,oy_min)
y_max<-max(my_max,oy_max)

mrfactcoord<-extent(mrf)
x1<-mrfactcoord[1]
x2<-mrfactcoord[2]
y1<-mrfactcoord[3]
y2<-mrfactcoord[4]

if (x_min < x1 ) { 
  x_min <- x1
} else if (x_max > x2) { 
  x_max <- x2
} else if (y_min < y1) { 
  y_min <- y1
} else if (y_max > y2) { 
  y_max <- y2
} else {
  x_min <- x_min
  x_max <- x_max
  y_min <- y_min
  y_max <- y_max
}

lon1=x_min
lon2=x_max
lat1=y_min
lat2=y_max
#
orfmax<-max(values(orf),na.rm=TRUE)
mrfmax<-max(values(mrf),na.rm=TRUE)
Rainlimt<-max(orfmax,mrfmax)+10
#
#layout(matrix(c(1, 2, 3, 0,4,0),2,3, nrow = T))
#layout(matrix(c(1, 2, 3, 0, 4, 0), nrow = 3, byrow = TRUE))
#layout.show(n=6)
##################################################################
library(rgdal)
shp <- readOGR(dsn = file.path("$shpfile"), stringsAsFactors = F)
mycols=c("darkturquoise","darkseagreen1","darkgreen","yellow","goldenrod1","violetred1","maroon4")
myrange=colorRampPalette(mycols)
mycrange=myrange(Rainlimt)
mycols = mycrange[c(0:Rainlimt)]

par(mfrow=c(3,2), family="serif", mai=c(0.4,0.4,0.3,0.3),cex.axis = 1.5,cex.main = 2,font = 2,cex.axis = 1.5)
## Obs plot
plota<-raster::plot(orf,
	xlim = c(lon1,lon2),
	ylim = c(lat1,lat2),
	zlim=c(${Threshold},Rainlimt),
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
pdata=read.delim("${OUT_PREFIX}_${ObjNum}_ScatterData.txt", header = TRUE, sep = "")
pdata[pdata<${Threshold}]=NA;
plot(x = pdata[,1],y = pdata[,2],
    xlim = c(${Threshold},Rainlimt),
    ylim = c(${Threshold},Rainlimt),
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
	xlim = c(lon1,lon2),
	ylim = c(lat1,lat2),
	zlim=c(${Threshold},Rainlimt),
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
	xlim = c(lon1,lon2),
	ylim = c(lat1,lat2),
	zlim=c(${Threshold},Rainlimt),
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

rows<-c("CRA Threshold=${CRA_threshold}    ")
# Push viewport
pushViewport(viewport(x=0.590,y=0.620,height=0.8,width=1.0))
grid.table(rows,theme=tt)
testdf1<-read.table("${OUT_PREFIX}_${ObjNum}_Table1.txt")
# Push viewport
pushViewport(viewport(x=0.610,y=0.40,height=0.8,width=0.5))
grid.table(testdf1,theme=tt)
#
testdf2<-read.table("${OUT_PREFIX}_${ObjNum}_Table2.txt")
# Push viewport
pushViewport(viewport(x=0.600,y=0.25,height=1.0,width=1.0))
grid.table(testdf2,theme=tt)
#
rows2<-c(substitute(bold("      Error Decomposition:                       ")))
# Push viewport
pushViewport(viewport(x=0.23,y=0.38,height=0.5,width=1.00))
grid.table(rows2,theme=tt)
#
testdf3<-read.table("${OUT_PREFIX}_${ObjNum}_Table3.txt",header = FALSE,sep = ",")
rows3<-c(paste(testdf3[1,1],"%  "))
rows4<-c(paste(testdf3[2,1],"%  "))
rows5<-c(paste(testdf3[3,1],"%  "))
rows6<-c(paste(testdf3[4,1],"%  "))
ab=rbind(rows3,rows4,rows5,rows6)
# Push viewport
pushViewport(viewport(x=0.500,y=0.27,height=0.9,width=1.3))
grid.table(ab,theme=tt,rows=NULL)
#
testdf4<-read.table("${OUT_PREFIX}_${ObjNum}_Table4.txt")
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
EOF

Rscript ${date}_PlotCRA_Th${Threshold}_${case}_${ObjNum}.R
mv ${date}_PlotCRA_Th${Threshold}_${case}_${ObjNum}.R ${OUT_PATH}
mv ${figname} ${OUT_PATH}

exit
