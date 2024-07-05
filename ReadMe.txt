******CRA Analysis Package (V2) ************
This package does the CRA Analysis 
Tested with R version - 3.4.3
Input file format tested: netcdf (Both observation and model forecast need not be having same resolution)

File system of CRA_package_V2.tar.gz:
1)INPUT 
	obs_20240527.nc          : Obs test data
	model_20240527.nc        : Model test data
2)SRC
	Isolate_Objects.R        : Function to regrid obs and model data, identify and pair objects from both model and observation based on rainfall threshold. This will provide each objects in separate NetCDF files. 
	CRA_Err_Decomp_Features.R: Main function for CRA analysis
3)SCRIPTS
	CRA_Calc_Submit.sh : This is a shell script which creates a R script for creating the identified paired object files and does CRA Analysis. It also provides a pdf which shows the paired objects and snapshots of the movement of the model raster over CRA area. It runs Run_CRA_Plot.sh within itself to create plots.
	Run_CRA_Plot.sh : Creates figure including the necessary statistic tables
4)OUTPUT
	All output files of this test data are kept here 
###################
User inputs for CRA_Calc_Submit.sh 

date=????????    #### date in the format of yyyymmdd
RUNPATH=         #### Path for rundirectory
ActINPUT_OBS=    #### path for observation file
ActINPUT_MODEL=  #### path for model file
OUT_PATH=        #### output path 
Threshold=       #### threshold rainfall for CRA Analysis 
case=            #### Experiment name
OUT_PREFIX=      #### output file prefix
Ngrids=6         #### No of grids to be shifted which is required for the RMSE minimization; 
			at present it is kept as 6 which means the model raster will be moving 
			6 grids right, left, down, and up over the CRA area. 


###################
User inputs for Run_CRA_Plot.sh
#### passing from CRA_Calc_Submit.sh and these line have to be change by the user while running separately.
export date=$1
export ObjNum=$2
export case=$3
export OUT_PATH=$4
export Threshold=$5
###

date=????????    #### date in the format of yyyymmdd
ObjNum=          #### serial number of the feature to be analysed 
RUNPATH=         #### Path for rundirectory 
OUT_PATH=        #### output path
Threshold=       #### threshold rainfall for CRA Analysis; should be same as in CRA_Calc_Submit.sh 
case=            #### Experiment name; should be same as in CRA_Calc_Submit.sh
OUT_PREFIX=      #### output file prefix; should be same as in CRA_Calc_Submit.sh
Ngrids=6         #### No of grids to be shifted which is required for the RMSE minimization; 
shpfile=         #### path of shp file 
figname=         #### name of the figure to be saved 
INPUT_OBS=       #### observation file which includes only the rainfall object
INPUT_MODEL=     #### model file  which includes only the rainfall object
INPUT_RotModel=  #### rotated object file 

###################

