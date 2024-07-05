Isolate_Objects<-function(ActOrf,ActMrf,objOBS,objMODEL,thr,minsize,sepdist){

library(raster)
library("SpatialVx")
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
ActOrf <- raster(ActOrf)
ActMrf <- raster(ActMrf)
################# Regrid ##########################################
obs_resolution <- res(ActOrf)
model_resolution <- res(ActMrf)
if ((obs_resolution[1] != model_resolution[1]) | (obs_resolution[2] != model_resolution[2])){
new_raster <- raster(extent(ActOrf), res = obs_resolution)
ActMrf<- resample(ActMrf, new_raster, method = "bilinear")
}
print(ActOrf)
print(ActMrf)
#plot(ActMrf,main="Model Regridded")
Obs_nm='Regrid_Obs.nc'
Mod_nm='Regrid_Model.nc'
#
Regrid_Obsnc <- writeRaster(ActOrf, filename=Obs_nm, format="CDF", overwrite=TRUE)
Regrid_Modnc <- writeRaster(ActMrf, filename=Mod_nm, format="CDF", overwrite=TRUE)
#######################################################
maxob<-maxvalue(ActOrf);maxmodel<-maxvalue(ActMrf)
if (thr>maxob & thr>maxmodel){
   stop(paste0('Error: The used threshold (',thr,') is more than the maximum rainfall in obs (', maxob, ') or in model (',maxmodel,'). Please change the threshold'))
}
else
{
################## Identify the Matched Objects ########
 p<-raster::as.matrix(ActOrf)
 q<-raster::as.matrix(ActMrf)
 
 p_rot<-rotate(p)
 q_rot<-rotate(q)
 
 hold <- make.SpatialVx(p_rot,q_rot)
 #print(hold)
 look1 <- FeatureFinder(hold, thresh = thr,min.size = c(minsize,minsize))  
 #print(look1)
 numO<-length(look1$X.feats)
 numM<-length(look1$Y.feats)
 
  if (numO<1 | numM<1) {
     stop(paste0('Warning!!! There are no objects identified with this threshold either in model(no. of objects:',numM,') or in observation (no. of objects:',numO,'). Please check the threshold'))
  }
  else
  {
   look2 <- minboundmatch( look1, type = "multiple", mindist =sepdist) 
   look <- MergeForce( look2 )
   
   SummFeats=summary(look)
   ObsSummary<-SummFeats$X
   ForecastSummary<-SummFeats$Y
   num_featO=dim(ObsSummary)
   num_featO<-num_featO[1]

    if (num_featO <1) {
       stop("Error !!! There are no paired objects found with this threshold. Function Isolate_Objects stops here. ")    
    } 
    else 
    {
     print(paste0('No. of identified features in Obs= ',num_featO[1]))
     num_featM=dim(ForecastSummary)
     num_featM<-num_featM[1]
     print(paste0('No. of identified features in Forecast= ',num_featM[1]))
     
     ### get the number of matched features ###
     MatchObj<-print(look$matches)
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
              assign(nam,SummFeats$X[ii,3])
              nam2=paste("Area_ModelF",ii,sep="")
              assign(nam2,SummFeats$Y[ii,3])
              ind_Obs<-(look$X.labeled)
              ind_Model<-(look$Y.labeled)
              ind_Obs[ind_Obs!=MatchFeatObs[ii]]<-NA
              ind_Model[ind_Model!=MatchFeatModel[ii]]<-NA
      #        image.plot(ind_Obs)
      #        image.plot(ind_Model)
      
      	tempObs=as.matrix(ActOrf)
      	tempObs<-rotate(tempObs)
      	new_Obs<-(ind_Obs*tempObs)/MatchFeatObs[ii]
     # 	image.plot(new_Obs)
      	newrastO<-ActOrf
      	new_Obs<-antirotate(new_Obs)
      	values(newrastO)=new_Obs
      #	plot(newrastO)
      png(paste0(file='ObsObj',ii,'.png'),bg = "white",  res = NA)
      	tempModel=as.matrix(ActMrf)
      	tempModel<-rotate(tempModel)
      	new_Model<-(ind_Model*tempModel)/MatchFeatModel[ii]
      	newrastM<-ActMrf
      	new_Model<-antirotate(new_Model)
      	values(newrastM)=new_Model
      #	plot(newrastM)
      png(paste0(file='ModelObj',ii,'.png'),bg = "white",  res = NA)
      
      orf<-newrastO
      mrf<-newrastM
      
      ######################
      Obj_Obs=paste(objOBS,'_Object',ii,'.nc',sep="")
      Obj_Model=paste(objMODEL,'_Object',ii,'.nc',sep="")
      ObObjnc <- writeRaster(orf, filename=Obj_Obs, format="CDF", overwrite=TRUE)
      MoObjnc <- writeRaster(mrf, filename=Obj_Model, format="CDF", overwrite=TRUE)
      }         #### for loop for objects ends here
      print(paste("Paired rainfall objects are saved as NetCDF files: ",Obj_Obs,"and",Obj_Model))
     }         #### end of if loop for no paired objects        
   }       #### end of if loop for no objects  
 }     #### end of if loop for threshold > max rain ends
}   
