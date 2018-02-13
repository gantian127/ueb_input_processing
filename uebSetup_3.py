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
workingDir = "/Projects/Tian_workspace/Animas_1989_watershed_large"
print(os.getcwd())
watershedN = 'Animas'
startYear=1988
endYear=2010
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

    netcdfFunctions.project_subset_and_resample_netCDF_to_referenceNetCDF2(cbrfcPWS, 'ogrid', refWSnc, refWSvar, intermP, in_epsgCode=epsgCbrfc,
                                                                               tSampling_interval=1, dT = 3.0, time_unitString=tuString)
    netcdfFunctions.project_subset_and_resample_netCDF_to_referenceNetCDF2(cbrfcTWS, 'otgrid', refWSnc, refWSvar, intermT, in_epsgCode=epsgCbrfc,
                                                                               tSampling_interval=1, dT = 3.0, time_unitString=tuString)
    print("P and Ta resampled")
    netcdfFunctions.convert_netcdf_Units(intermP, precIn, 'ogrid', new_varUnits='m/hr', multiplier_Factor=0.008467)     #0.00847 => inch/3hr = m/hr
    netcdfFunctions.convert_netcdf_Units(intermT, tempIn, 'otgrid', new_varUnits='oC', multiplier_Factor=0.5556, offset=-17.7778)   #C = 0.5556F + (- 17.7778)
    print("P and Ta unit converted ")
    #

    netcdfFunctions.project_subset_and_resample_netCDF_to_referenceNetCDF2(nldas2WS, 'SPFH2m_110_HTGL', refWSnc, refWSvar, intermQ, in_epsgCode=epsgNldas,
                                                                               time_unitString=tuString, dT=1.0, tSampling_interval=3, in_Xcoord='lon_110', in_Ycoord='lat_110')
    netcdfFunctions.project_subset_and_resample_netCDF_to_referenceNetCDF2(nldas2WS, 'UGRD10m_110_HTGL', refWSnc, refWSvar, intermU, in_epsgCode=epsgNldas,
                                                                               time_unitString=tuString, dT=1.0, tSampling_interval=3, in_Xcoord='lon_110', in_Ycoord='lat_110')
    netcdfFunctions.project_subset_and_resample_netCDF_to_referenceNetCDF2(nldas2WS, 'VGRD10m_110_HTGL', refWSnc, refWSvar, intermV, in_epsgCode=epsgNldas,
                                                                               time_unitString=tuString, dT=1.0, tSampling_interval=3, in_Xcoord='lon_110', in_Ycoord='lat_110')
    netcdfFunctions.project_subset_and_resample_netCDF_to_referenceNetCDF2(nldas2WS, 'PRESsfc_110_SFC', refWSnc, refWSvar, intermPress, in_epsgCode=epsgNldas,
                                                                               time_unitString=tuString, dT=1.0, tSampling_interval=3, in_Xcoord='lon_110', in_Ycoord='lat_110')
    print("Wind and humidity resampled")
    #
    climateFunctions_2_ori.compute_average_windSpeed(intermU, 'UGRD10m_110_HTGL', intermV, 'VGRD10m_110_HTGL', windIn, 'windS10')
    climateFunctions_2_ori.compute_vaporPressure_from_SpecficHumidity(intermQ, 'SPFH2m_110_HTGL', intermPress, 'PRESsfc_110_SFC', vpIn, 'VP')
    print("Wind resultant computed")

    #elevation adjustment
    """
    watershedFunctions_2.project_and_resample_Raster_to_ReferenceNetcdf(nldasDEM,watershedN,projDEMnlads)
    watershedFunctions_2.project_and_resample_Raster_to_ReferenceNetcdf(nedDEM,watershedN,projDEMned)
    climateFunctions_2.adjust_for_elevation_VaporPressure(vpIn,vpIn,'VP',projDEMnldas,projDEMned,baseDateTime=str(cYear)+'/10/01 0')
    #pressure elevation adj
    #temperature elevation adj
    #wind elevation adjustement--logarithmic profile
    #precip elevatin adj. if required
    #elevation adj spechum
    """
print("done")

