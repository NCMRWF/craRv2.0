### Read Object files identified by feature finder ###
### and does CRA Analysis                          ###
CRA_Err_Decomp<-function(ActOrf,ActMrf,orf,mrf,thr,objnum,out_prefix,Ngrids,iopt){

library(raster)
############# Functions for calculation of RMSE, COG and Statistics ########
RMSE <- function(x, y) { sqrt(mean((x - y)^2)) }
#######
CalcCor<-function(x,y){
	x1<-rasterToPoints(x)
	y1<-rasterToPoints(y)
	x1[x1<=0]<-NA
	y1[y1<=0]<-NA
	data<-cbind(x1[,3],y1[,3])
	data2<-data[complete.cases(data),]
	coreln<-cor(data2[,1],data2[,2],use='complete.obs')
	return(coreln)  }
#######
cog <- function(r){
	ptop<-rasterToPoints(r,na.rm=TRUE)
	xpt<-ptop[,1] ; ypt<-ptop[,2] ; rf<-ptop[,3]
	xp<-sum(xpt*rf)/sum(rf);yp<-sum(ypt*rf)/sum(rf)
	return(c(xp,yp))
	}
#######
objcentroid <- function(dt){
            non_zero_indices <- which(dt != 0, arr.ind = TRUE)
            return(colMeans(non_zero_indices))
            }
#######
dtrotate <- function(x) t(apply(x, 2, rev))
antirotate <- function(x) apply(t(x),2,rev)
#######
approx_na_2d <- function(x, radius = 1) {
   result <- x
   na_indices <- which(x==0, arr.ind = TRUE)
   for (i in 1:nrow(na_indices)) {
    row_index <- na_indices[i, 1]
    col_index <- na_indices[i, 2]
    row_start <- max(1, row_index - radius)
    row_end <- min(nrow(x)-1, row_index + radius)
    col_start <- max(1, col_index - radius)
    col_end <- min(ncol(x)-1, col_index + radius)
    nearby_values <- x[row_start:row_end, col_start:col_end]
    nearby_values <- nearby_values[!is.na(nearby_values)]
# Calculate the average of nearby non-NA values
    if (length(nearby_values) > 0) {
      result[row_index, col_index] <- mean(nearby_values)
       }
   }
return(result)
}
############# for rotation ###########
obj_rotate<-function(obs_to_rotate,model_to_rotate,rot_ref_rms,rot_ref_corr){
          print("Calculating Error due to Rotation of the Rainfall Object...............")
          #print(paste("reference corr =",rot_ref_corr)) 
          matob=as.matrix(obs_to_rotate)
          matob[is.na(matob)]<-0
          matob<-dtrotate(matob)
#          #image.plot(matob)
          matmodel=as.matrix(model_to_rotate)
          matmodel[is.na(matmodel)]<-0
          matmodel<-dtrotate(matmodel)
#          #image.plot(matmodel)
          num_rows <- nrow(matob)
          num_cols <- ncol(matob)
          rot_centre<-objcentroid(matmodel)
          angle<-0.
          ob_mat=as.vector(matob);
          rot_mod=as.vector(matmodel);
          valid_indices <- !is.na(ob_mat) & !is.na(rot_mod)
          ob_mat_valid <- ob_mat[valid_indices]
          rot_mod_valid <- rot_mod[valid_indices]
          avgob=mean(ob_mat_valid) ; avgrot=mean(rot_mod_valid) 
          sdob=sqrt(mean((ob_mat_valid-avgob)*(ob_mat_valid-avgob)))
          sdmod=sqrt(mean((rot_mod_valid-avgrot)*(rot_mod_valid-avgrot)))
          rot_mat_save<-matmodel
          
          for(angle_degrees in seq(from=0, to=350, by=10)) {
             # Convert angle from degrees to radians
             angle_radians <- angle_degrees * (pi / 180)
             print(paste("Running for angle ",angle_degrees,"........."))
             rotational_matrix <- matrix(c(cos(angle_radians), -sin(angle_radians),
                                           sin(angle_radians), cos(angle_radians)),
                                         nrow = 2)
             
             rotated_data <- matrix(0, nrow = num_rows, ncol = num_cols) ### Added due to out of bounds error 
             non_zero_indices<- which(matmodel != 0 , arr.ind = TRUE)
             non_zero_lines<-dim(non_zero_indices)
             for(nl in 1:non_zero_lines[1]) {
              i<-non_zero_indices[nl,1]
              j<-non_zero_indices[nl,2]
              # Translate the coordinates to the centroid
              translated_i <- i - floor(rot_centre[1]) + 1
              translated_j <- j - floor(rot_centre[2]) + 1
              point <- c(floor(translated_i), floor(translated_j))  # Translated coordinates
              rotated_point <- rotational_matrix %*% point
              rotated_i <- floor(rotated_point[1]) + floor(rot_centre[1]) - 1
              rotated_j <- floor(rotated_point[2]) + floor(rot_centre[2]) - 1
              # Retain only the rotated indices within bounds
              if(rotated_i >= 1 && rotated_i <= num_rows && rotated_j >= 1 && rotated_j <= num_cols ) { 
              rotated_data[rotated_i, rotated_j] <- matmodel[i, j]
                 }
               }
          
            myMatrix <- approx_na_2d(rotated_data)
            matob[matob==0]<-NA; myMatrix[myMatrix==0]<-NA;
            ## plotting Observation base and forecast object overlay 
            image(matob,col=c("white","blue"), main = "Rotated Data over Obs")
            legend("topright", legend = "Obs", fill = c("blue"), bty = "n")
            image(myMatrix, col = c("cyan"),add = TRUE)
            legend("bottomleft", legend="Rotated Forecast Object",fill = c("cyan"), bty = "n")
            # Mark the point of rotation
            points(rot_centre[1],rot_centre[2], col = "black", pch = 16)
           
            ### Calculate correlation and RMS after rotating model
             rot_mod<-as.vector(myMatrix); rot_mod[rot_mod==0]<-NA;   
             ob_mat<-as.vector(matob); ob_mat[ob_mat==0]<-NA;
             valid_indices <- !is.na(ob_mat) & !is.na(rot_mod)
             if (any(valid_indices)) {
               ob_mat_valid <- ob_mat[valid_indices]
               rot_mod_valid <- rot_mod[valid_indices]
               rot_rms <- RMSE(ob_mat_valid, rot_mod_valid)
               rot_corr <- cor(ob_mat_valid, rot_mod_valid, use = 'complete.obs')
             } else {
               ob_mat_valid <- numeric(0)  # or another appropriate default
               rot_mod_valid <- numeric(0) # or another appropriate default
               rot_rms <- NA
               rot_corr <- NA
             }
             rot_corr[is.na(rot_corr)]<-(-99999)
             nangle<-9999
              if((angle_degrees <=350) & (rot_corr > rot_ref_corr)) {
#               print("Got better rotated position.. ")
               print(paste("Got better rotated position for angle ",angle_degrees,"...."))
               ### Statistics ###
               rot_ref_corr<-rot_corr
               rot_ref_rms<-rot_rms
               angle<-angle_degrees
               nangle<-angle_degrees
               avgob=mean(ob_mat_valid) ; avgrot=mean(rot_mod_valid) 
               sdob=sqrt(mean((ob_mat_valid-avgob)*(ob_mat_valid-avgob)))
               sdmod=sqrt(mean((rot_mod_valid-avgrot)*(rot_mod_valid-avgrot)))
               rot_mat_save<-myMatrix
               }
             }   ### rotation loop ends
            return(list(angle=angle,rot_ref_corr=rot_ref_corr,rot_ref_rms=rot_ref_rms,sdob=sdob,sdmod=sdmod,array1=ob_mat_valid,array2=rot_mod_valid,avgob=avgob,avgrot=avgrot,rotated_Object=rot_mat_save))
####### If there is no any angle of rotation could make the model to get the best overlap
            if((angle_degrees >=350) & (nangle > 350)){
               angle<-0
               ob_mat=as.vector(matob);
               rot_mod=as.vector(matmodel);
               valid_indices <- !is.na(ob_mat) & !is.na(rot_mod)
               ob_mat_valid <- ob_mat[valid_indices]
               rot_mod_valid <- rot_mod[valid_indices]
               avgob=mean(ob_mat_valid) ; avgrot=mean(rot_mod_valid) 
               sdob=sqrt(mean((ob_mat_valid-avgob)*(ob_mat_valid-avgob)))
               sdmod=sqrt(mean((rot_mod_valid-avgrot)*(rot_mod_valid-avgrot)))
               rot_mat_save<-matmodel
               return(list(angle=angle,rot_ref_corr=rot_ref_corr,rot_ref_rms=rot_ref_rms,sdob=sdob,sdmod=sdmod,array1=ob_mat_valid,array2=rot_mod_valid,avgob=avgob,avgrot=avgrot,rotated_Object=rot_mat_save))
              }
            #return(list(angle=angle,rot_ref_corr=rot_ref_corr,rot_ref_rms=rot_ref_rms,sdob=sdob,sdmod=sdmod,array1=ob_mat_valid,array2=rot_mod_valid,avgob=avgob,avgrot=avgrot,rotated_Object=rot_mat_save))
            }  ##### obj_rotate function ends
################################################

objstats<-function(rf){
	xgrd=xres(rf)
	ltmin<-rf@extent[3];ltmax<-rf@extent[4]
	ptop<-rasterToPoints(rf,na.rm=TRUE)
	orfpts<-ptop[,3] ;nptso<-length(ptop[,3])
	rfmean<-mean(orfpts) ; rfmax<-max(orfpts) # Getting mean & max from stats
	pi=22/7
	cosfac=cos(0.5*(ltmin+ltmax)*pi/180.)
	factor=xgrd*xgrd*111.*(cosfac*111.)*1.e-6    #to convert to rain volume [km^3]
	rfvol=nptso*rfmean*factor
	rfarea=(nptso*factor)*1.e6
	return(c(nptso, rfmean, rfmax, rfvol, rfarea))
	}
###

################################ Reading input file #########################
ActOrf<-raster(ActOrf)  ## Actual Input Obs
ActMrf<-raster(ActMrf)        ## Actual input model
orf<-raster(orf)        ## Feature from Obs
mrf<-raster(mrf)              ## Feature from Model
#############################################################################
xgrd=xres(orf); ygrd=yres(orf); grd=xgrd
CRAmask<-(sum(orf,mrf,na.rm=TRUE))           # get the overlap of model and observation
CRAmask[is.na(orf) & is.na(mrf)]<-NA
CRAmask_orig<-CRAmask      
plot(CRAmask,col='red',main="Verification Mask without Shift")
######### Calulate all stats over the identified feature  ##########
orfstatsThr<-objstats(orf)
mrfstatsThr<-objstats(mrf)
OrfThrnpts<-orfstatsThr[1]; MrfThrnpts<-mrfstatsThr[1]
OrfThravg<-orfstatsThr[2] ; MrfThravg<-mrfstatsThr[2]
OrfThrmax<-orfstatsThr[3] ; MrfThrmax<-mrfstatsThr[3]
OrfThrvol<-orfstatsThr[4] ; MrfThrvol<-mrfstatsThr[4]
OrfThrarea<-orfstatsThr[5]; MrfThrarea<-mrfstatsThr[5]  
Orfpts<-rasterToPoints(orf,na.rm=TRUE);Orfpts<-Orfpts[,3]
Mrfpts<-rasterToPoints(mrf,na.rm=TRUE);Mrfpts<-Mrfpts[,3]
valid_indices <- !is.na(Orfpts) & !is.na(Mrfpts)
Orfpts_valid <- Orfpts[valid_indices]
Mrfpts_valid <- Mrfpts[valid_indices]
sdF=sqrt(mean((Orfpts_valid-OrfThravg)*(Orfpts_valid-OrfThravg)))
sdO=sqrt(mean((Mrfpts_valid-MrfThravg)*(Mrfpts_valid-MrfThravg)))
Act_orig_corr<-cor(Orfpts_valid,Mrfpts_valid,use='complete.obs')
MSEtotalCorrMeth=((OrfThravg-MrfThravg)^2)+((sdO-Act_orig_corr*sdF)^2)+(1-Act_orig_corr^2)*(sdF)^2
############## Apply CRA area mask on Obs and Model ##################
a<-ActOrf ; b<-ActMrf
a[is.na(CRAmask)]<-NA; b[is.na(CRAmask)]<-NA
####   Calculate COG for Obs and Model over CRA Area ###
COGorf<-cog(a); COGmrf<-cog(b)
## Calulate all stats over the CRA masked area ################
astats<-objstats(a);bstats<-objstats(b)
anpts<-astats[1]; bnpts<-bstats[1]
aavg<-astats[2] ; bavg<-bstats[2]
amax<-astats[3] ; bmax<-bstats[3]
avol<-astats[4] ; bvol<-bstats[4]
aarea<-astats[5]; barea<-bstats[5]  
apts<-rasterToPoints(a,na.rm=TRUE);apts<-apts[,3]
bpts<-rasterToPoints(b,na.rm=TRUE);bpts<-bpts[,3]
numapts<-astats[1]

##### Calculate RMSE & Corr before shifting ######
orig_cor<-CalcCor(a,b); old_cor<-orig_cor    
origrmsdata<-cbind(apts,bpts)
origrmsdata2<-origrmsdata[complete.cases(origrmsdata),]
orig_rms<-RMSE(origrmsdata2[,1],origrmsdata2[,2]); old_rms<-orig_rms
print(paste("RMS and correlation before shift and rotation are: ",orig_rms," and ",orig_cor))
################################################
print("RMS Error Minimization Starts....")
	#iopt=1 # 1. Sq. Err minimization 2. for patcorr maximization
igrd=Ngrids   # MAX NUMBER OF GRIDS TO SCAN IN W TO E & N TO S ; 169 SCANS
idx=0
prev_rms=99999
prev_cor=-99999
no_shift=0
for (c in -igrd:igrd){
 dxgrd=c*grd
  for (r in -igrd:igrd){
   idx=idx+1
   dygrd=r*grd
   newb<-shift(ActMrf,dx=dxgrd, dy=dygrd)
   newmod<-shift(mrf,dx=dxgrd, dy=dygrd)             

#### Updation of verification mask based on shifted Model field #########
   verif_mask<-merge(newmod, orf)
#   verif_mask<-merge(newCRAmask_b, orf, tolerance = 0.5)
### 
   figtitle<-paste("Verification Mask after shifting (",c,",",r,")")
   plot(verif_mask,col='red',main=figtitle,xlim=c((extent(ActMrf)[1])-(igrd+1),(extent(ActMrf)[2]+igrd)),ylim=c((extent(ActMrf)[3])-(igrd+1),(extent(ActMrf)[4])+igrd))
   abline(v = (COGorf[1]+dygrd), col = "green", lty = 2,add=TRUE )  # Vertical grid line
   abline(h = (COGorf[2]+dxgrd), col = "green", lty = 2,add=TRUE ) # Horizontal grid line
##
   av<-ActOrf ; abv<-ActMrf; bv<-newb
   av[is.na(verif_mask)]<-NA; abv[is.na(verif_mask)]<-NA; bv[is.na(verif_mask)]<-NA

   ptav<-rasterToPoints(av, na.rm=TRUE); avpts<-ptav[,3]
   ptabv<-rasterToPoints(abv, na.rm=TRUE); abvpts<-ptabv[,3]
   ptbv<-rasterToPoints(bv, na.rm=TRUE); bvpts<-ptbv[,3]
   lnavpts<-length(avpts)
   chkpts<-(lnavpts/numapts)*100
## Calculate the correlation and RMSE after shifting the model field ##
   shift_cor<-CalcCor(av,bv) ; 
   rmsdata<-cbind(avpts,bvpts)
   rmsdata2<-rmsdata[complete.cases(rmsdata),]
   shift_rms<-RMSE(rmsdata2[,1],rmsdata2[,2])
   shift_rms[is.na(shift_rms)]<-99999
   shift_cor[is.na(shift_cor)]<-(-99999)

##### Calculates COG,Corr and CRA Error Components if shift_rms < old_rms (Opt1) ########
    if((shift_cor > prev_cor)) {
     #print(paste("Got better translation point...","RMSE=",shift_rms,"and","Corr=",shift_cor))
     print(paste("Got better translation point... Translated for (",c,",",r,")"))
     model_to_rotate<-bv
     obs_to_rotate<-av
     xpts<-c*grd
     ypts<-r*grd
#########################################
     prev_rms<-shift_rms
     prev_cor<-shift_cor
     pattcorr=shift_cor ; RMSshift=shift_rms
     avstats<-objstats(av);bvstats<-objstats(abv)
     avnpts<-avstats[1] ; bvnpts<-bvstats[1]
     avavg<-avstats[2] ; bvavg<-bvstats[2]
     avmax<-avstats[3] ; bvmax<-bvstats[3]
     avvol<-avstats[4] ; bvvol<-bvstats[4]
     avarea<-avstats[5] ; bvarea<-bvstats[5]  
     mrfstats_shift<-summary(bvpts) ; mrfmax_shift<-mrfstats_shift[6]   
     mrfmax=bmax 
#     
      if(chkpts < 50.){
       print('MORE THAN 50% OF POINTS SHIFTED OUT') ; outcode1=1
      }else{
       print('MORE THAN 50% OF POINTS USED') ; outcode1=0
       }
      if(mrfmax_shift != mrfmax){
       print('MAX RAINFALL LOST') ; outcode2=1
      }else{
       print('MAX RAINFALL RETAINED') ; outcode2=0
       }
     no_shift=no_shift+1;
    } 
  }   
 } ### for loop for translation ends

if(iopt==2){
    if(no_shift>0){
     Part_rot<-obj_rotate(obs_to_rotate,model_to_rotate,RMSshift,pattcorr)
     print("Rotation ends for this translation point...")
     Err_angle=Part_rot$angle;
     if(Err_angle>180){Err_angle=Err_angle-360}
      print(Err_angle)
      rot_fin_corr=Part_rot$rot_ref_corr;
     # print(rot_fin_corr)
      rot_fin_rms=Part_rot$rot_ref_rms;
     # print(rot_fin_rms)
      sdob_rot=Part_rot$sdob;
      sdmod_rot=Part_rot$sdmod;
      pts_ob=Part_rot$array1;
     # print(length(pts_ob))
      pts_rot_mod=Part_rot$array2;
     # print(length(pts_rot_mod))
      avgob=Part_rot$avgob;
      avgrot=Part_rot$avgrot;
      rotated_model=Part_rot$rotated_Object;
     # print(dim(rotated_model))
## save rotated model data
      rotnewrast <- ActMrf
      rot_model_obj <- antirotate(rotated_model)
      values(rotnewrast) <- rot_model_obj
      rot_rast_nm <- paste(out_prefix,'RotatedModel_',objnum,'.nc', sep="")
      MoObjnc <- writeRaster(rotnewrast, filename = rot_rast_nm, format = "CDF", overwrite = TRUE)
      }
    else {
     print("There is no spatial shift")
     obs_to_rotate=orf;
     model_to_rotate=mrf;
     shift_rms=orig_rms;
     shift_cor=orig_cor;
     Part_rot<-obj_rotate(obs_to_rotate,model_to_rotate,shift_rms,shift_cor)
     print("Rotation ends for this translation point...")
     Err_angle=Part_rot$angle;
     if(Err_angle>180){Err_angle=Err_angle-360}
     rot_fin_corr=Part_rot$rot_ref_corr;
     rot_fin_rms=Part_rot$rot_ref_rms;
     sdob_rot=Part_rot$sdob;
     sdmod_rot=Part_rot$sdmod;
     pts_ob=Part_rot$array1;
     #print(length(pts_ob))
     pts_rot_mod=Part_rot$array2;
     #print(length(pts_rot_mod))
     avgob=Part_rot$avgob;
     avgrot=Part_rot$avgrot;
     pattcorr=orig_cor
     RMSshift=orig_rms
     }
    } ### stats calculation from rotated object ends
##########

 if(iopt==1){
     print('Successfully completed the CRA Analysis using squared error minimization method')
     avstats<-objstats(av);bvstats<-objstats(abv)
     avnpts<-avstats[1] ; bvnpts<-bvstats[1]
     avavg<-avstats[2] ; bvavg<-bvstats[2]
     avmax<-avstats[3] ; bvmax<-bvstats[3]
     avvol<-avstats[4] ; bvvol<-bvstats[4]
     avarea<-avstats[5] ; bvarea<-bvstats[5]  
     mrfstats_shift<-summary(bvpts) ; mrfmax_shift<-mrfstats_shift[6]   
     mrfmax=bmax 
     
     COGmrf_shift<-cog(bv)
     avpts[is.na(bvpts)]<-NA; bvpts[is.na(avpts)]<-NA;
     ll<-avpts ; mm<-abvpts ; kk<-bvpts
     difr=(na.omit(kk)-na.omit(ll)) ; toterrshift=sum(difr*difr)
     avga=mean(na.omit(ll)) ; avgb=mean(na.omit(mm)) ; avgbshift=mean(na.omit(kk))
     sda=sqrt(mean((na.omit(ll)-avga)*(na.omit(ll)-avga)))
     sdb=sqrt(mean((na.omit(mm)-avgb)*(na.omit(mm)-avgb)))  
     sdbshift=sqrt(mean((na.omit(kk)-avgbshift)*(na.omit(kk)-avgbshift)))
     Observed_Rain<-ll
     Forecast_Rain_Shifted<-kk
     rfscat<-cbind(Observed_Rain, Forecast_Rain_Shifted)
     ############## CRA error decomposition & other statistics ##############
     rot_fin_rms=NA
     pattcorr_rot=NA
     Err_angle=NA
     MSErotation=NA
     stats0<-c(OrfThrnpts,MrfThrnpts,OrfThravg,MrfThravg,OrfThrmax,MrfThrmax,OrfThrvol,MrfThrvol,OrfThrarea,MrfThrarea)
     stats0<-t(round(stats0,digits=2))
     stats1<-c(anpts,bnpts,aavg,bavg,amax,bmax,avol,bvol,aarea,barea,orig_rms, orig_cor)
     stats1<-t(round(stats1,digits=2))
     stats2<-c(avnpts,bvnpts,avavg,bvavg,avmax,bvmax,avvol,bvvol,avarea,bvarea,
               orig_rms,RMSshift,orig_cor,pattcorr)
     stats2<-t(round(stats2,digits=2))
     RMSorig=orig_rms
     MSEtotal=RMSorig*RMSorig               #total error in CRA
##
     MSEshift=RMSshift*RMSshift           #error left after two blobs superimposed
     MSEdisplacement=MSEtotal-MSEshift    #difference in RMSE before & after shift
     MSEvolume=(avgbshift-avga)*(avgbshift-avga)        #mean intensity difference
     MSEpattern=MSEshift-MSEvolume        #residual fine scale pattern difference
##
     } #### stats calculation for only translation

#### Error decomposition for Opt2 (Maximizing Pattern Correlation Method)
   if(iopt==2){
    print('Successfully completed the CRA Analysis using Pattern Correlation Method')
    pattcorr_rot=rot_fin_corr
    corrorig=orig_cor
    RMSorig=orig_rms
    COGmrf_shift<-cog(model_to_rotate)
    COGorf<-cog(obs_to_rotate)
    Observed_Rain<-pts_ob
    Forecast_Rain_Shifted<-pts_rot_mod
    rfscat<-cbind(Observed_Rain, Forecast_Rain_Shifted)
    ############## CRA error decomposition & other statistics ##############
    stats0<-c(OrfThrnpts,MrfThrnpts,OrfThravg,MrfThravg,OrfThrmax,MrfThrmax,OrfThrvol,MrfThrvol,OrfThrarea,MrfThrarea)
    stats0<-t(round(stats0,digits=2))
    stats1<-c(anpts,bnpts,aavg,bavg,amax,bmax,avol,bvol,aarea,barea,orig_rms, orig_cor)
    stats1<-t(round(stats1,digits=2))
    stats2<-c(avnpts,bvnpts,avavg,bvavg,avmax,bvmax,avvol,bvvol,avarea,bvarea,
              orig_rms,rot_fin_rms,orig_cor,rot_fin_corr)
    stats2<-t(round(stats2,digits=2))
##########
    MSEdisplacement=2.*sdob_rot*sdmod_rot*(pattcorr-corrorig)        #uses correlation before and after shift
    MSEvolume=(avgob-avgrot)*(avgob-avgrot) # mean intensity difference (note this is not identical to MSEvolume for minERR)
    MSEpattern=2.*sdob_rot*sdmod_rot*(1.-pattcorr_rot) + (sdob_rot-sdmod_rot)^2. # conditional bias + imperfect corr. terms
    MSErotation=2*sdob_rot*sdmod_rot*(pattcorr_rot-pattcorr)
    MSEtotal=MSEdisplacement+MSEvolume+MSEpattern+MSErotation
    }
####### Convert CRA Error components to percent ###########
MSEdisplacement=MSEdisplacement/MSEtotal*100.
print(paste("MSEdisplacement=",MSEdisplacement))
MSEvolume      =MSEvolume      /MSEtotal*100.
print(paste("MSEvolume=",MSEvolume))
MSEpattern     =MSEpattern     /MSEtotal*100.
print(paste("MSEpattern=",MSEpattern))
MSErotation     =MSErotation   /MSEtotal*100.
print(paste("MSErotation=",MSErotation))
MSEcomps<-c(MSEtotal,MSEdisplacement,MSEvolume,MSEpattern,MSErotation)
MSEcomps<-t(round(MSEcomps,digits=2))
##
table1_fig <- matrix(c(OrfThrnpts,MrfThrnpts,
                       OrfThravg,MrfThravg,
                       OrfThrmax,MrfThrmax,
                       OrfThrvol,MrfThrvol), nrow = 4, byrow = TRUE)
table1_fig<-round(table1_fig,digits=1)
colnames(table1_fig) <- c("Obs", "Fcst")
rownames(table1_fig) <- c(paste("Gridpoints>=",thr,"mm/day"), "Average Rain Rate (mm/day)", "Maximum Rain (mm/day)", "Rain Volume (km^3)")
##
table2_fig<- matrix(c(orig_rms,RMSshift,rot_fin_rms,
                       orig_cor,pattcorr,pattcorr_rot),nrow=2,byrow=TRUE)
table2_fig<-round(table2_fig,digits=2)
colnames(table2_fig) <- c("Original", "Shifted","Shift+Rotation")
rownames(table2_fig) <- c("RMSE(mm/day)", "Correlation coefficient")
###
table3_fig<- matrix(c(MSEdisplacement,
                      MSEvolume,
                      MSEpattern,
                      MSErotation),nrow=4,byrow=TRUE)
#table3_fig<-round(table3_fig,digits=1)
numeric_indices <- sapply(table3_fig, is.numeric)
table3_fig[numeric_indices] <- round(as.numeric(table3_fig[numeric_indices]), digits = 2)
rownames(table3_fig) <- c("Displacement Error", "Volume Error","Pattern Error","Rotation Error")

####
cogs<-c(COGorf[1],COGorf[2],COGmrf[1],COGmrf[2],COGmrf_shift[1],COGmrf_shift[2])
cogs<-t(round(cogs,digits=2))
orig_londif=(COGorf[1]-COGmrf[1]);
orig_latdif=(COGorf[2]-COGmrf[2]);
shift_londif=(COGmrf[1]-COGmrf_shift[1]);
shift_latdif=(COGmrf[2]-COGmrf_shift[2]);

table4_fig=matrix(c(xpts,ypts,Err_angle),nrow=3,byrow=TRUE)
rownames(table4_fig) <- c("Lon (Deg.)", "Lat (Deg.)","Angle (Deg.)")
colnames(table4_fig) <- c("Location Error")
###
cogs<-c(COGorf[1],COGorf[2],COGmrf[1],COGmrf[2],COGmrf_shift[1],COGmrf_shift[2])
cogs<-t(round(cogs,digits=2))
xdif_orig=abs(COGorf[1]-COGmrf[1])*111.1 ; ydif_orig=abs(COGorf[2]-COGmrf[2])*111.1
xdif_shift=abs(COGorf[1]-COGmrf_shift[1])*111.1 ; ydif_shift=abs(COGorf[2]-COGmrf_shift[2])*111.1
dist_orig=sqrt((xdif_orig*xdif_orig)+(ydif_orig*ydif_orig))
dist_shift=sqrt((xdif_shift*xdif_shift)+(ydif_shift*ydif_shift))
xydistcog<-c(xdif_orig,ydif_orig,xdif_shift,ydif_shift,dist_orig,dist_shift)
xydistcog<-t(round(xydistcog,digits=2))

######################### Writing Statistics to files ################################################
write.table(rfscat,file=paste(out_prefix,objnum,'ScatterData.txt',sep="_"),sep="   ",col.names=T,row.names=F,append=F)
#write.table(cogs,file=paste(out_prefix,objnum,'cog.txt',sep="_"),sep="   ",col.names=F,row.names=F,append=F)
write.table(table1_fig,file=paste(out_prefix,objnum,'Table1.txt',sep="_"),sep="   ",col.names=T,row.names=T,append=F)
write.table(table2_fig,file=paste(out_prefix,objnum,'Table2.txt',sep="_"),sep="   ",col.names=T,row.names=T,append=F)
write.table(table3_fig,file=paste(out_prefix,objnum,'Table3.txt',sep="_"),sep="   ",col.names=F,row.names=T,append=F)
write.table(table4_fig,file=paste(out_prefix,objnum,'Table4.txt',sep="_"),sep="   ",col.names=T,row.names=T,append=F)
}
