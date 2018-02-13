"""
This is the code to further prepare the climate forcing: (They will become climate forcing data as separate files for each year)
- resample climate forcing data based on watershed.nc
- convert units for precipitation and temperature
- compute average wind speed using 10-m above ground zonal wind speed, 10-m above ground meridional wind speed from NLDAS data
- compute vapor pressure using specific humidity and surface pressure from NLDAS data
"""
import os
import netcdfFunctions, climateFunctions_2_ori

#### To-DO:    Elevation adjustment line to be added here    <<<<<<<<<<<<==============

"""*********** Input settings for watershed of interest *****************"""
workingDir = "/Projects/Tian_workspace/22yr_Animas_watershed_ueb/fix_poor_data/final/"
print(os.getcwd())
watershedN = 'Animas'
startYear=1999
endYear=2000
cbrfcPWS = watershedN+'P.nc'
cbrfcTWS = watershedN+'T.nc'
nldas2WS = watershedN+'NLDAS2.nc'
refWSnc = 'watershed.nc'
refWSvar = "watershed"
epsgCbrfc = 4326
epsgNldas = 4326
"""****************intermediate files*******************************************"""
intermP = 'inch3hr_'+cbrfcPWS
intermT = 'f3hr_'+cbrfcTWS
intermQ = watershedN+'Q.nc'
intermU = watershedN+'U10.nc'
intermV = watershedN+'V10.nc'
intermPress = watershedN+'Press.nc'
#"""
os.chdir(workingDir)
for cYear in range(startYear,endYear):
    cbrfcPWS = watershedN+str(cYear+1)+'P.nc'
    cbrfcTWS = watershedN+str(cYear+1)+'T.nc'
    nldas2WS = watershedN+str(cYear+1)+'NLDAS2.nc'
    tuString = 'hours since '+str(cYear)+'-10-01 00:00:00 UTC'
    """**************final outputs****************************************************"""
    precIn = 'prcp'+str(cYear-startYear)+'.nc'
    tempIn = 'temp'+str(cYear-startYear)+'.nc'
    windIn = 'windS10m'+str(cYear-startYear)+'.nc'
    vpIn   = 'VP'+str(cYear-startYear)+'.nc'

#    netcdfFunctions.project_subset_and_resample_netCDF_to_referenceNetCDF2(cbrfcPWS, 'ogrid', refWSnc, refWSvar, intermP, in_epsgCode=epsgCbrfc,
#                                                                               tSampling_interval=1, dT = 3.0, time_unitString=tuString)
    netcdfFunctions.project_subset_and_resample_netCDF_to_referenceNetCDF2(cbrfcTWS, 'otgrid', refWSnc, refWSvar, intermT, in_epsgCode=epsgCbrfc,
                                                                               tSampling_interval=1, dT = 3.0, time_unitString=tuString)
    print("P and Ta resampled")
#    netcdfFunctions.convert_netcdf_Units(intermP, precIn, 'ogrid', new_varUnits='m/hr', multiplier_Factor=0.008467)     #0.00847 => inch/3hr = m/hr
    netcdfFunctions.convert_netcdf_Units(intermT, tempIn, 'otgrid', new_varUnits='oC', multiplier_Factor=0.5556, offset=-17.7778)   #C = 0.5556F + (- 17.7778)
    print("P and Ta unit converted ")

print("done")

