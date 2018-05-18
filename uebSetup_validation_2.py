"""
This is used to prepare the terrain and wind, vp, tmax, tmin, prec, temp hourly data used for running RDHM UEB for validation purpose

This is based on the code of uebinRDHM_InputSetup_2.py.
The difference is that this code won't do units conversion for temp and do units conversion for precp as mm/hr

"""
import os
import netcdfFunctions
import climateFunctions_2
from datetime import datetime
import callSubprocess
import watershedFunctions

## Domain bounding box in geographic coordinates left, top, right, bottom.  Must enclose watershed of interest

# Mcphee
leftX, topY, rightX, bottomY = -108.80, 38.05, -107.66, 37.22  # exact box: -108.601067, 37.857910, -107.863539, 37.428745
watershedN = 'Mcphee_MPHC2'

startYear = 2010  # datetime.strptime(startDateTime,"%Y/%m/%d %H").year
endYear = 2015  # datetime.strptime(endDateTime,"%Y/%m/%d %H").year
startMonthDayHour = "10/01 0"
endMonthDayHour = "10/01 0"

workingDir = "/Projects/Tian_workspace/rdhm_ueb_modeling/McPhee_MPHC2/MPHC2_forcing_validation/"
os.chdir(workingDir)
os.mkdir('Forcing')
##reference xmrg file ---must be located in the workingDir

# proj4_string: see the paper: Reed, S.M., and D.R. Maidment, "Coordinate Transformations for Using NEXRAD Data in GIS-based Hydrologic Modeling," Journal of Hydrologic Engineering, 4, 2, 174-182, April 1999
## http://www.nws.noaa.gov/ohd/hrl/distmodel/hrap.htm#backgroundworkingDir
proj4_string = 'PROJCS["Sphere_ARC_INFO_Stereographic_North_Pole",GEOGCS["GCS_Sphere_ARC_INFO",DATUM["Sphere_ARC_INFO",SPHEROID["Sphere_ARC_INFO",6370997,0]],PRIMEM["Greenwich",0],\
UNIT["degree",0.0174532925199433]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",60.00681388888889],PARAMETER["central_meridian",-105],\
PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]]]'
##proj4_string = '+proj=stere +lat_0=90.0 +lat_ts=60.0 +lon_0=-105.0 +k=1 +x_0=0.0 +y_0=0.0 +a=6371200 +b=6371200 +units=m +no_defs'

inputXmrgRaster = 'we0101198906z.gz'
referenceRasterASCII = watershedN+'_refRaster.asc'
referenceRasterTIF = watershedN+'_refRaster.tif'
referenceRasterNC = watershedN+'_refRaster.nc'
cmdString = 'xmrgtoasc -p ster -i ' + inputXmrgRaster + ' -o ' + referenceRasterASCII
callSubprocess.callSubprocess(cmdString, "xmrg to ascii")
cmdString = "gdal_translate -a_srs \"" + proj4_string + "\" "+referenceRasterASCII+" "+referenceRasterTIF
callSubprocess.callSubprocess(cmdString, "ascii to tif")
watershedFunctions.rasterToNetCDF(referenceRasterTIF, referenceRasterNC)

# cbrfc forcing
#cbrfc forcing is converted to xmrg for the whole upper colorado --similar to what RTI are doing
print 'start tair, prec'
workingDir1 = "/Projects/cbrfcTP/"
os.chdir(workingDir1)
targetDir1 = workingDir +"Forcing/"
for year in range(startYear, endYear + 1):
    cmdString = "for i in mm_" + str(year) + "*.nc.gz; do gunzip -cdfq $i>${i//.gz} && ncea -d x," + str(leftX) + "," + str(rightX) \
                + " -d y,"+str(bottomY) + "," + str(topY) + " -O ${i//.gz} ${i//.gz} && mv -f -t " + targetDir1 + " ${i//.gz}; done"
    callSubprocess.callSubprocess(cmdString, 'subset nc and move')

os.chdir(targetDir1)
cmdString = \
    'for i in mm_*.nc; do ' \
    '   ncks --mk_rec_dmn time   -O $i $i &&  ' \
    '   ncks -d time,0 -v otgrid -O $i Tair_${i//.nc}00.nc; '\
    '   ncks -d time,0 -v otgrid -O $i Tair_${i//.nc}01.nc; '\
    '   ncks -d time,0 -v otgrid -O $i Tair_${i//.nc}02.nc; '\
    '   ncks -d time,1 -v otgrid -O $i Tair_${i//.nc}03.nc; '\
    '   ncks -d time,1 -v otgrid -O $i Tair_${i//.nc}04.nc; '\
    '   ncks -d time,1 -v otgrid -O $i Tair_${i//.nc}05.nc; '\
    '   ncks -d time,2 -v otgrid -O $i Tair_${i//.nc}06.nc; '\
    '   ncks -d time,2 -v otgrid -O $i Tair_${i//.nc}07.nc; '\
    '   ncks -d time,2 -v otgrid -O $i Tair_${i//.nc}08.nc; '\
    '   ncks -d time,3 -v otgrid -O $i Tair_${i//.nc}09.nc; '\
    '   ncks -d time,3 -v otgrid -O $i Tair_${i//.nc}10.nc; '\
    '   ncks -d time,3 -v otgrid -O $i Tair_${i//.nc}11.nc; '\
    '   ncks -d time,4 -v otgrid -O $i Tair_${i//.nc}12.nc; '\
    '   ncks -d time,4 -v otgrid -O $i Tair_${i//.nc}13.nc; '\
    '   ncks -d time,4 -v otgrid -O $i Tair_${i//.nc}14.nc; '\
    '   ncks -d time,5 -v otgrid -O $i Tair_${i//.nc}15.nc; '\
    '   ncks -d time,5 -v otgrid -O $i Tair_${i//.nc}16.nc; '\
    '   ncks -d time,5 -v otgrid -O $i Tair_${i//.nc}17.nc; '\
    '   ncks -d time,6 -v otgrid -O $i Tair_${i//.nc}18.nc; '\
    '   ncks -d time,6 -v otgrid -O $i Tair_${i//.nc}19.nc; '\
    '   ncks -d time,6 -v otgrid -O $i Tair_${i//.nc}20.nc; '\
    '   ncks -d time,7 -v otgrid -O $i Tair_${i//.nc}21.nc; '\
    '   ncks -d time,7 -v otgrid -O $i Tair_${i//.nc}22.nc; '\
    '   ncks -d time,7 -v otgrid -O $i Tair_${i//.nc}23.nc; '\
    '   ncks -d time,0 -v ogrid -O $i Prec_${i//.nc}00.nc; '\
    '   ncks -d time,0 -v ogrid -O $i Prec_${i//.nc}01.nc; '\
    '   ncks -d time,0 -v ogrid -O $i Prec_${i//.nc}02.nc; '\
    '   ncks -d time,1 -v ogrid -O $i Prec_${i//.nc}03.nc; '\
    '   ncks -d time,1 -v ogrid -O $i Prec_${i//.nc}04.nc; '\
    '   ncks -d time,1 -v ogrid -O $i Prec_${i//.nc}05.nc; '\
    '   ncks -d time,2 -v ogrid -O $i Prec_${i//.nc}06.nc; '\
    '   ncks -d time,2 -v ogrid -O $i Prec_${i//.nc}07.nc; '\
    '   ncks -d time,2 -v ogrid -O $i Prec_${i//.nc}08.nc; '\
    '   ncks -d time,3 -v ogrid -O $i Prec_${i//.nc}09.nc; '\
    '   ncks -d time,3 -v ogrid -O $i Prec_${i//.nc}10.nc; '\
    '   ncks -d time,3 -v ogrid -O $i Prec_${i//.nc}11.nc; '\
    '   ncks -d time,4 -v ogrid -O $i Prec_${i//.nc}12.nc; '\
    '   ncks -d time,4 -v ogrid -O $i Prec_${i//.nc}13.nc; '\
    '   ncks -d time,4 -v ogrid -O $i Prec_${i//.nc}14.nc; '\
    '   ncks -d time,5 -v ogrid -O $i Prec_${i//.nc}15.nc; '\
    '   ncks -d time,5 -v ogrid -O $i Prec_${i//.nc}16.nc; '\
    '   ncks -d time,5 -v ogrid -O $i Prec_${i//.nc}17.nc; '\
    '   ncks -d time,6 -v ogrid -O $i Prec_${i//.nc}18.nc; '\
    '   ncks -d time,6 -v ogrid -O $i Prec_${i//.nc}19.nc; '\
    '   ncks -d time,6 -v ogrid -O $i Prec_${i//.nc}20.nc; '\
    '   ncks -d time,7 -v ogrid -O $i Prec_${i//.nc}21.nc; '\
    '   ncks -d time,7 -v ogrid -O $i Prec_${i//.nc}22.nc; '\
    '   ncks -d time,7 -v ogrid -O $i Prec_${i//.nc}23.nc; '\
    '   ncks --mk_rec_dmn time $i R_$i &&  ' \
    '   ncra -v otgrid -y min -d time,,,8,8 R_$i Tamin_${i//.nc}00.nc; ' \
    '   ncra -v otgrid -y max -d time,,,8,8 R_$i Tamax_${i//.nc}00.nc && ' \
    '   rm -f R_$i; ' \
    '   rm -f $i; ' \
    'done'
print(cmdString)
callSubprocess.callSubprocess(cmdString, 'subset nc in time')

cmdString = "for i in Prec*.nc; do " \
                "ncap2 -s'ogrid=8.467f*ogrid' -O $i $i && "\
                "ncatted -a units,ogrid,m,c,'mm/hr' -O $i $i; "\
                "ncatted -a long_name,ogrid,m,c,'Hourly precipitation rate' -O $i $i; " \
                "done"
callSubprocess.callSubprocess(cmdString, 'convert netcdf units')
cmdString = "for i in Tair*.nc; do " \
                "ncatted -a units,otgrid,m,c,'oF' -O $i $i; "\
                "ncatted -a long_name,otgrid,m,c,'Temperature' -O $i $i; " \
                "done"
callSubprocess.callSubprocess(cmdString, 'convert netcdf units')
cmdString = "for i in Tamin*.nc; do " \
                "ncatted -a units,otgrid,m,c,'oF' -O $i $i; "\
                "ncatted -a long_name,otgrid,m,c,'Daily minimum temperature' -O $i $i; " \
                "done"
callSubprocess.callSubprocess(cmdString, 'convert netcdf units')
cmdString = "for i in Tamax*.nc; do " \
                "ncatted -a units,otgrid,m,c,'oF' -O $i $i; "\
                "ncatted -a long_name,otgrid,m,c,'Daily maximum temperature' -O $i $i; " \
                "done"
callSubprocess.callSubprocess(cmdString, 'convert netcdf units')

"""
# this can be used when sampling to watershed domain defined by connectivity file (output from RDHM)
#cbrfc forcing
workingDir1 = "/Projects/cbrfcTP/"
os.chdir(workingDir1)
targetDirP = workingDir+"Prec/"
targetDirT = workingDir+"Tair/"
targetDirTmin = workingDir+"Tmin/"
targetDirTmax = workingDir+"Tmax/"
for cYear in range(startYear,endYear):
    startDateTime=str(cYear)+"/"+startMonthDayHour
    endDateTime=str(cYear+1)+"/"+endMonthDayHour
    cbrfcPWS = watershedN+str(cYear+1)+'P.nc'
    cbrfcTWS = watershedN+str(cYear+1)+'T.nc'
    cbrfcTminWS = watershedN + str(cYear + 1) + 'Tmin.nc'
    cbrfcTmaxWS = watershedN + str(cYear + 1) + 'Tmax.nc'
    input_pref   ix = 'mm'
    #these two functions take long to finish--run only if not have the file already
    climateFunctions_2.create_netCDF_from_multple_nc(input_prefix, watershedN, cbrfcPWS, 'ogrid', leftX, topY, rightX, bottomY, startDateTime, endDateTime, 3, out_timeName = 'time')
    cmdString = 'mv -f -t ' + targetDirP + ' ' + cbrfcPWS
    callSubprocess.callSubprocess(cmdString, "move subset nc")
    climateFunctions_2.create_netCDF_from_multple_nc(input_prefix, watershedN, cbrfcTWS, 'otgrid', leftX, topY, rightX, bottomY, startDateTime, endDateTime, 3, out_timeName = 'time')
    cmdString = 'mv -f -t ' + targetDirT + ' ' + cbrfcTWS
    callSubprocess.callSubprocess(cmdString, "move subset nc")
    # min / max temp
    climateFunctions_2.get_daily_minmax_netCDF_from_multple_nc(input_prefix, watershedN, cbrfcTminWS, cbrfcTmaxWS, 'otgrid', leftX, topY, rightX, bottomY, startDateTime, endDateTime, 3, out_timeName='time')
    cmdString = 'mv -f -t ' + targetDirTmin + ' ' + cbrfcTminWS
    callSubprocess.callSubprocess(cmdString, "move subset nc")
    cmdString = 'mv -f -t ' + targetDirTmax + ' ' + cbrfcTmaxWS
    callSubprocess.callSubprocess(cmdString, "move subset nc")
#endfor
"""

#nldas2
print 'start wind, vp preparation'
workingDir2 = "/Projects/nldasMonthly/"
targetDirVPW = workingDir+"Forcing/"
intermQ = watershedN + 'Q.nc'
intermPress = watershedN + 'Press.nc'
intermU = watershedN + 'U10.nc'
intermV = watershedN + 'V10.nc'
refWSvar = "Band1"
epsgNldas = 4326
for cYear in range(startYear,endYear):
    startDateTime = str(cYear) + "/" + startMonthDayHour
    endDateTime = str(cYear + 1) + "/" + endMonthDayHour
    tuString = 'hours since ' + str(cYear) + '-10-01 00:00:00 UTC'
    nldas2WS = watershedN+str(cYear+1)+'NLDAS2.nc'
    vpIn = watershedN + str(cYear + 1) + 'VP.nc'
    windIn = watershedN+str(cYear+1)+'windS10m.nc'

    os.chdir(workingDir2)
    climateFunctions_2.subset_and_concatenate_netCDF_from_multple_nc_NLDAS('NLDAS_FORA0125_H.A_Monthly', watershedN, nldas2WS, leftX, topY, rightX, bottomY, startDateTime, endDateTime,
                                                                             dT=1, in_Xcoord = 'lon_110', in_Ycoord='lat_110',inout_timeName = 'time')
    cmdString = 'mv -f -t ' + workingDir + ' ' + nldas2WS
    callSubprocess.callSubprocess(cmdString, "move subset nc")

    os.chdir(workingDir)
    netcdfFunctions.project_subset_and_resample_netCDF_to_referenceNetCDF2(nldas2WS,'SPFH2m_110_HTGL',referenceRasterNC,refWSvar,intermQ,in_epsgCode=epsgNldas,
                                                                           time_unitString=tuString,dT=1.0,tSampling_interval=1,in_Xcoord='lon_110',in_Ycoord='lat_110')
    netcdfFunctions.project_subset_and_resample_netCDF_to_referenceNetCDF2(nldas2WS,'UGRD10m_110_HTGL',referenceRasterNC,refWSvar,intermU,in_epsgCode=epsgNldas,
                                                                           time_unitString=tuString,dT=1.0,tSampling_interval=1,in_Xcoord='lon_110',in_Ycoord='lat_110')
    netcdfFunctions.project_subset_and_resample_netCDF_to_referenceNetCDF2(nldas2WS,'VGRD10m_110_HTGL',referenceRasterNC,refWSvar,intermV,in_epsgCode=epsgNldas,
                                                                           time_unitString=tuString,dT=1.0,tSampling_interval=1,in_Xcoord='lon_110',in_Ycoord='lat_110')
    netcdfFunctions.project_subset_and_resample_netCDF_to_referenceNetCDF2(nldas2WS,'PRESsfc_110_SFC',referenceRasterNC,refWSvar,intermPress,in_epsgCode=epsgNldas,
                                                                           time_unitString=tuString,dT=1.0,tSampling_interval=1,in_Xcoord='lon_110',in_Ycoord='lat_110')
    #"""
    climateFunctions_2.compute_vaporPressure_from_SpecficHumidity(intermQ, 'SPFH2m_110_HTGL', intermPress, 'PRESsfc_110_SFC', vpIn, 'VP')
    climateFunctions_2.compute_average_windSpeed(intermU,'UGRD10m_110_HTGL',intermV,'VGRD10m_110_HTGL',windIn,'windS10')

    #### for no-REsampling--takes the NLDAS resolution ~ 13 Km
    """
    cmdString = 'nccopy -V UGRD10m_110_HTGL ' + nldas2WS + ' ' + intermU
    callSubprocess.callSubprocess(cmdString, "intermediate U nc")
    cmdString = 'nccopy -V VGRD10m_110_HTGL ' + nldas2WS + ' ' + intermV
    callSubprocess.callSubprocess(cmdString, "intermediate V nc")
    cmdString = 'nccopy -V SPFH2m_110_HTGL ' + nldas2WS + ' ' + intermQ
    callSubprocess.callSubprocess(cmdString, "intermediate Q nc")
    cmdString = 'nccopy -V PRESsfc_110_SFC ' + nldas2WS + ' ' + intermPr
    callSubprocess.callSubprocess(cmdString, "intermediate Pressure nc")
    windIn = watershedN + str(cYear + 1) + 'windS10m.nc'
    vpIn = watershedN + str(cYear + 1) + 'VP.nc'
    climateFunctions_2.compute_average_windSpeed(intermU, 'UGRD10m_110_HTGL', intermV, 'VGRD10m_110_HTGL', windIn,'windS10')
    climateFunctions_2.compute_vaporPressure_from_SpecficHumidity(intermQ, 'SPFH2m_110_HTGL', intermPr, 'PRESsfc_110_SFC', vpIn, 'VP')
    """
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

    cmdString = 'mv -f -t ' + targetDirVPW + ' ' + vpIn
    callSubprocess.callSubprocess(cmdString, "move subset nc")
    cmdString = 'mv -f -t ' + targetDirVPW + ' ' + windIn
    callSubprocess.callSubprocess(cmdString, "move subset nc")


print("done")


