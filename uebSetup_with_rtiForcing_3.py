import os
import netcdfFunctions, climateFunctions_2
from datetime import datetime
import rdhmFunctions
import callSubprocess
"""*********** Convert UEB NC files to XMRG with HRAP projection *****************"""
workingDir = "/Projects/Tian_workspace/rdhm_ueb_modeling/McPhee_MPHC2/MPHC2_forcing/"
watershedN = 'Mcphee_MPHC2'
startDateTime = "1988/10/01 0"
endDateTime = "2010/10/01"
startYear = 1988  # datetime.strptime(startDateTime,"%Y/%m/%d %H").year
endYear = 2010  # datetime.strptime(endDateTime,"%Y/%m/%d %H").year
time_varName='time'

#proj4_string: see the paper: Reed, S.M., and D.R. Maidment, "Coordinate Transformations for Using NEXRAD Data in GIS-based Hydrologic Modeling," Journal of Hydrologic Engineering, 4, 2, 174-182, April 1999
## http://www.nws.noaa.gov/ohd/hrl/distmodel/hrap.htm#background
proj4_string = 'PROJCS["Sphere_ARC_INFO_Stereographic_North_Pole",GEOGCS["GCS_Sphere_ARC_INFO",DATUM["Sphere_ARC_INFO",SPHEROID["Sphere_ARC_INFO",6370997,0]],PRIMEM["Greenwich",0],\
UNIT["degree",0.0174532925199433]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",60.00681388888889],PARAMETER["central_meridian",-105],\
PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]]]'
##proj4_string = '+proj=stere +lat_0=90.0 +lat_ts=60.0 +lon_0=-105.0 +k=1 +x_0=0.0 +y_0=0.0 +a=6371200 +b=6371200 +units=m +no_defs'

####site variables
targetDirSiteV = workingDir+'SiteV/'
siteVarFiles = ['{}{}'.format(watershedN, var) for var in ['CC30m.nc', 'Hcan30m.nc', 'LAI30m.nc', 'Aspect30m.tif', 'Slope30m.tif']]
ueb_site_vars = ['ueb_cc', 'ueb_hcan', 'ueb_lai', 'ueb_aspect', 'ueb_slope']
#os.chdir(workingDir)
os.chdir(targetDirSiteV)
for indx in range(3):
    # convert nc to tif and project
    cmdString = "gdalwarp -t_srs \"" + proj4_string + "\" NETCDF:\"" + siteVarFiles[indx] + "\":Band1 " + ueb_site_vars[indx]+".tif"
    callSubprocess.callSubprocess(cmdString, 'project to Stereographic')
    cmdString = "gdal_translate -of AAIGrid  " + ueb_site_vars[indx]+".tif " +ueb_site_vars[indx] + ".asc"
    callSubprocess.callSubprocess(cmdString, "nc to ascii")
for indx in range(3, 5):
    # project raster
    cmdString = "gdalwarp -t_srs \"" + proj4_string + "\" " + siteVarFiles[indx] + " " + ueb_site_vars[indx] + ".tif"
    callSubprocess.callSubprocess(cmdString, 'project to Stereographic')
    cmdString = "gdal_translate -of AAIGrid  " + ueb_site_vars[indx] + ".tif " + ueb_site_vars[indx] + ".asc"
    callSubprocess.callSubprocess(cmdString, "tif to ascii")
#to xmrg
#os.chdir(targetDirSiteV)
cmdString = " for i in *.asc; do asctoxmrg -i $i -p ster -f par; done"
callSubprocess.callSubprocess(cmdString, "ascii to xmrg")

### forcing
# #cbrfc forc
# inpVar =  ["Prec", "Tair"]
# forcingTargetDir = workingDir+"Forcing/"
# os.chdir(forcingTargetDir)
# for indx in range(2):
#     cmdString = "for i in "+inpVar[indx]+"*.nc; do gdalwarp -ot Float32 -s_srs EPSG:4326 -t_srs \"" + proj4_string + "\"" + " $i ueb${i:0:4}${i:12:2}${i:14:2}${i:8:4}${i:16:2}z.tif; done"
#     callSubprocess.callSubprocess(cmdString, "project forcing")
#     # ascii
#     cmdString = "for i in *.tif; do gdal_translate -of AAIGrid $i ${i//.tif}.asc; done"
#     callSubprocess.callSubprocess(cmdString, "tif to ascii")
#     #to xmrg
#     cmdString = "for i in *.asc; do asctoxmrg -i $i -f par -p ster; done"
#     callSubprocess.callSubprocess(cmdString, "ascii to xmrg")
#     #del intrm files
#     cmdString = " rm -f *tif"
#     callSubprocess.callSubprocess(cmdString, "delete intermediate files")
#     cmdString = " rm -f *asc"
#     callSubprocess.callSubprocess(cmdString, "delete intermediate files")
#     cmdString = " rm -f *.xml *.prj"
#     callSubprocess.callSubprocess(cmdString, "delete intermediate files")
#     ##delete nc file to save space
#     cmdString = " rm -f *.nc"
#     #callSubprocess.callSinpVar =  ["Prec", "Tair", "Tamin", "Tamax"]

#max/min Ta
print 'start Tmax, Tmin preparation'
forcingTargetDir = workingDir+"Forcing/"
os.chdir(forcingTargetDir)
inpVar =  ["Tamin", "Tamax"]
for indx in range(2):
    cmdString = "for i in "+inpVar[indx]+"*.nc; do gdalwarp -ot Float32 -s_srs EPSG:4326 -t_srs \"" + proj4_string + "\"" + " $i ueb${i:0:5}${i:13:2}${i:15:2}${i:9:4}${i:17:2}z.tif; done"
    callSubprocess.callSubprocess(cmdString, "project forcing")
    # ascii
    cmdString = "for i in *.tif; do gdal_translate -of AAIGrid $i ${i//.tif}.asc; done"
    callSubprocess.callSubprocess(cmdString, "tif to ascii")
    #to xmrg
    cmdString = "for i in *.asc; do asctoxmrg -i $i -f par -p ster; done"
    callSubprocess.callSubprocess(cmdString, "ascii to xmrg")
    #del intrm files
    cmdString = " rm -f *tif"
    callSubprocess.callSubprocess(cmdString, "delete intermediate files")
    cmdString = " rm -f *asc"
    callSubprocess.callSubprocess(cmdString, "delete intermediate files")
    cmdString = " rm -f *.xml *.prj"
    callSubprocess.callSubprocess(cmdString, "delete intermediate files")
    ##delete nc file to save space
    cmdString = " rm -f *.nc"
    #callSubprocess.callSubprocess(cmdString, "delete orig files")

#NLDAS VP and WindS
print 'start wind, vp preparation'
forcingTargetDir = workingDir+"Forcing/"
os.chdir(forcingTargetDir)
startMonthDayHour = "10/01 0"
endMonthDayHour = "10/01 0"
#NLDAS VP and WindS
for cYear in range(startYear,endYear):
    print cYear
    nldasinputForc = [watershedN+str(cYear + 1)+'VP.nc',watershedN+str(cYear + 1)+'windS10m.nc']
    outVars = ["uebVp", "uebWindS"]
    nldasvarForc = ['VP', 'windS10']
    startDateTime = str(cYear) + "/" + startMonthDayHour
    for indx in range(2):
        print nldasinputForc[indx]
        print nldasvarForc[indx]
        print startDateTime
        rdhmFunctions.create_multiple_Ascii_fromNetCDF(nldasinputForc[indx], outVars[indx], nldasvarForc[indx], 1,
                                                       startDateTime=startDateTime, time_varName=time_varName,proj4_string=proj4_string)
        print 'done with create multiple netcdf'
        cmdString = "for i in *.asc; do asctoxmrg -i $i -f par -p ster; done"
        callSubprocess.callSubprocess(cmdString, "ascii to xmrg")
        # delete the ascii files
        # cmdString = " rm -f *.asc *asc.aux.xml *.prj"
        # callSubprocess.callSubprocess(cmdString, "delete ascii")

print("done")


