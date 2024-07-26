#!/bin/bash
module load gnu/R/3.4.3
##################################################
date=$1  ### yyyymmdd
RUNPATH=`pwd`
ActINPUT_MODEL=../INPUT/model_20240527.nc
ActINPUT_OBS=../INPUT/obs_20240527.nc 
OUT_PATH=../OUTPUT/${date}
Threshold=40
mkdir -p ${OUT_PATH}
case=Test
OUT_PREFIX=${OUT_PATH}/${date}_Thr${Threshold}_${case}
Ngrids=6 
minsize=10
sepdist=1
iopt=2  # 2: with rotation error and 1:without rotation error 
##################################################

cd $RUNPATH
cat > Calc_CRA_${date}_${case}_Th${Threshold}.R << EOF
library(raster)
library("SpatialVx")
dev.new()
###
rotate <- function(x) t(apply(x, 2, rev))	
antirotate <- function(x) apply(t(x),2,rev)
maxvalue<-function(rf){
	ptop<-rasterToPoints(rf,na.rm=TRUE)
	orfpts<-ptop[,3]; rfmax<-max(orfpts)
	rfmax<-round(rfmax,digits=2)
	return(rfmax)
	}
################## INPUTS ##############################
ActOrf <- raster("${ActINPUT_OBS}")
ActMrf <- raster("${ActINPUT_MODEL}")
################# Regrid ##########################################
target_resolution <- res(ActOrf)
new_raster <- raster(extent(ActOrf), res = target_resolution)
ActMrf <- resample(ActMrf, new_raster, method = "bilinear")
print(ActOrf)
print(ActMrf)
plot(ActMrf,main="Model Regridded")
Obs_nm='Regrid_Obs_${date}.nc'
Mod_nm='Regrid_Mod_${date}.nc'

Regrid_Obsnc <- writeRaster(ActOrf, filename=Obs_nm, format="CDF", overwrite=TRUE)
Regrid_Modnc <- writeRaster(ActMrf, filename=Mod_nm, format="CDF", overwrite=TRUE)
#######################################################
maxob<-maxvalue(ActOrf);maxmodel<-maxvalue(ActMrf)
if (${Threshold}>maxob & ${Threshold}>maxmodel){
   stop(paste0('Error: The used threshold (',${Threshold},') is more than the maximum rainfall in obs (', maxob, ') or in model (',maxmodel,'). Please change the threshold'))
}
################## Identify the Matched Objects ########
p<-raster::as.matrix(ActOrf)
q<-raster::as.matrix(ActMrf)

p_rot<-rotate(p)
q_rot<-rotate(q)

hold <- make.SpatialVx(p_rot,q_rot)
print(hold)
look1 <- FeatureFinder(hold, thresh = ${Threshold},min.size = c(${minsize},${minsize})) 
plot(look1)
numO<-length(look1\$X.feats)
numM<-length(look1\$Y.feats)
  if (numO<1 | numM<1) {
     stop(paste0('Warning!!! There are no objects identified with this threshold either in model(no. of objects:',numM,') or in observation (no. of objects:',numO,'). Please check the threshold'))
  }

look2 <- minboundmatch( look1, type = "multiple", mindist =${sepdist})
look <- MergeForce( look2 )
#plot(look)

SummFeats=summary(look)
ObsSummary<-SummFeats\$X
ForecastSummary<-SummFeats\$Y
num_featO=dim(ObsSummary)
num_featO<-num_featO[1]

    if (num_featO <1) {
       stop("Error !!! There are no paired objects found with this threshold. Function Isolate_Objects stops here. ")    
    } 

print(paste0('No. of identified features in Obs= ',num_featO[1]))
num_featM=dim(ForecastSummary)
num_featM<-num_featM[1]
print(paste0('No. of identified features in Forecast= ',num_featM[1]))

### get the number of matched features ###
MatchObj<-print(look\$matches)
Num_Match<-MatchObj[,1]
TotalNumFeat=length(Num_Match)
MatchFeatModel<-MatchObj[,1]
MatchFeatObs<-MatchObj[,2]

############### Detaching Packages which stop running ´shift´ function ##########
detach(package:SpatialVx)
detach(package:turboEM)
detach(package:quantreg)
detach(package:SparseM)
detach(package:numDeriv)
detach(package:doParallel)
detach(package:iterators)
detach(package:foreach)
detach(package:smatr)
detach(package:smoothie)
detach(package:fields)
detach(package:maps)
detach(package:spam)
detach(package:dotCall64)
detach(package:spatstat)
detach(package:rpart)
detach(package:nlme)
detach(package:spatstat.data)
library(fields)
##
for (ii in 1:TotalNumFeat){
        nam=paste("Area_ObF",ii,sep="")
        assign(nam,SummFeats\$X[ii,3])
        nam2=paste("Area_ModelF",ii,sep="")
        assign(nam2,SummFeats\$Y[ii,3])
        ind_Obs<-(look\$X.labeled)
        ind_Model<-(look\$Y.labeled)
        ind_Obs[ind_Obs!=MatchFeatObs[ii]]<-NA
        ind_Model[ind_Model!=MatchFeatModel[ii]]<-NA
        #image.plot(ind_Obs,main="Object from Observation")
        image(ind_Obs,col="blue", main = "Object from Observation",axes=FALSE)
        axis(side = 1)  # X-axis
        axis(side = 2)
        axis(side = 3)
        axis(side = 4)
        #image.plot(ind_Model,main="Object from Model")
        image(ind_Model,col="red", main = "Object from Model")

	tempObs=as.matrix(ActOrf)
	tempObs<-rotate(tempObs)
	new_Obs<-(ind_Obs*tempObs)/MatchFeatObs[ii]
#	image.plot(new_Obs)
	newrastO<-ActOrf
	new_Obs<-antirotate(new_Obs)
	values(newrastO)=new_Obs
	plot(newrastO,main="Extracted Observation Data")

	tempModel=as.matrix(ActMrf)
	tempModel<-rotate(tempModel)
	new_Model<-(ind_Model*tempModel)/MatchFeatModel[ii]
	newrastM<-ActMrf
	new_Model<-antirotate(new_Model)
	values(newrastM)=new_Model
	plot(newrastM,main="Extracted Model Data")

orf<-newrastO
mrf<-newrastM

######################
Obj_Obs=paste('Obs_${case}_${date}_Thr${Threshold}','_Object',ii,'.nc',sep="")
Obj_Model=paste('Model_${case}_${date}_Thr${Threshold}','_Object',ii,'.nc',sep="")
ObObjnc <- writeRaster(orf, filename=Obj_Obs, format="CDF", overwrite=TRUE)
MoObjnc <- writeRaster(mrf, filename=Obj_Model, format="CDF", overwrite=TRUE)
print("-----------------------------------------------------------------")
print(paste("Object files used for analysis are",Obj_Obs,"and",Obj_Model))
print("-----------------------------------------------------------------")
##################################################################
source("../SRC/CRA_Err_Decomp_Features.R")
CRA_Err_Decomp(Obs_nm,Mod_nm,Obj_Obs,Obj_Model,thr=${Threshold},objnum=ii,"${OUT_PREFIX}",Ngrids=${Ngrids},iopt=${iopt})
print(".................. Running  for plotting ......................... TESTING...............")
iopt=${iopt}
if (iopt==2) {
args=c("${date}",ii,"${case}","${OUT_PATH}","${Threshold}")
system2("../SCRIPTS/Run_CRA_Plot.sh",args=args)
} else {
args=c("${date}",ii,"${case}","${OUT_PATH}","${Threshold}")
system2("../SCRIPTS/Run_CRA_Plot_NoRot.sh",args=args)
}
}
dev.off()
EOF

Rscript Calc_CRA_${date}_${case}_Th${Threshold}.R
#mv *_Object*.nc ${OUT_PATH}
#mv Calc_CRA_${date}_${case}_Th${Threshold}.R ${OUT_PATH}
#mv Rplots.pdf ${OUT_PATH}/CRA_Features_${date}_${case}.pdf

exit
