__author__ = 'tsega'
from hydrods_python_client import HydroDS
import settings
from datetime import datetime, timedelta
"""*********** Input settings for watershed of interest *****************"""
workingDir = "/Projects/rdhmUEB/green2009/"
## Domain bounding box in geographic coordinates left, top, right, bottom.  Must enclose watershed of interest
# use rectangular domain--no WS delineation

# Green River near Daniel at Warren Bridge
leftX, topY, rightX, bottomY =  -110.415, 43.593, -109.492, 42.871
# Animas River WS above Durango
#leftX, topY, rightX, bottomY =  -108.15, 38.06, -107.41, 37.16

# Grid projection
#utmZone = int((180 + 0.5*(xmin+xmax))/6) + 1
epsgCode = 26912    #26912                 #26912 utm 12
dx,dy  = 30, 30  #  Grid cell sizes (m) for reprojection
# Set parameters for watershed delineation
streamThreshold = 1000
watershedName = 'Green'
# Cell spacing for subsampled UEB model (m)
dxRes, dyRes = 800, 800
#### model start and end dates
startDateTime = "2009/10/01 0"
endDateTime = "2010/10/01 0"
"""*************************************************************************"""
HDS = HydroDS(username=settings.USER_NAME, password=settings.PASSWORD)

MyFiles = HDS.list_my_files()
for item in MyFiles:
    print(item)

#### Subset DEM and Delineate Watershed
input_static_DEM  = 'nedWesternUS.tif'
subsetDEM_request = HDS.subset_raster(input_raster=input_static_DEM, left=leftX, top=topY, right=rightX,
                                      bottom=bottomY,output_raster=watershedName + 'DEM84.tif')
#Options for projection with epsg full list at: http://spatialreference.org/ref/epsg/
WatershedDEM = HDS.project_resample_raster(input_raster_url_path=subsetDEM_request['output_raster'],
                                                      cell_size_dx=dx, cell_size_dy=dy, epsg_code=epsgCode,
                                                      output_raster=watershedName+'Proj'+str(dx)+'.tif',resample='bilinear')
####Resample watershed grid to coarser grid (to save time)
domainDEM = watershedName + 'DEM' + str(dxRes) + 'm.tif'
Watershed =  HDS.resample_raster(input_raster_url_path= WatershedDEM['output_raster'],
                cell_size_dx=dxRes, cell_size_dy=dyRes, resample='near', output_raster=domainDEM)
#terrain variables
aspect_hires = HDS.create_raster_aspect(input_raster_url_path=WatershedDEM['output_raster'],
                                output_raster=watershedName + 'Aspect' + str(dx)+ '.tif')
aspectRaster = watershedName + 'Aspect' + str(dxRes) + 'm.tif'
aspect = HDS.resample_raster(input_raster_url_path= aspect_hires['output_raster'], cell_size_dx=dxRes,
                                cell_size_dy=dyRes, resample='near', output_raster=aspectRaster)
slope_hires = HDS.create_raster_slope(input_raster_url_path=WatershedDEM['output_raster'],
                                output_raster=watershedName + 'Slope' + str(dx) + '.tif')
slopeRaster = watershedName + 'Slope' + str(dxRes) + 'm.tif'
slope = HDS.resample_raster(input_raster_url_path= slope_hires['output_raster'], cell_size_dx=dxRes,
                                cell_size_dy=dyRes, resample='near', output_raster=slopeRaster)
#Land cover variables
nlcd_raster_resource = 'nlcd2011CONUS.tif'
subset_NLCD_result = HDS.project_clip_raster(input_raster=nlcd_raster_resource,
                                ref_raster_url_path=Watershed['output_raster'],
                                output_raster=watershedName + 'nlcdProj' + str(dxRes) + '.tif')
ccfracNC = watershedName+'CC'+str(dxRes)+'m.nc'
ccfrac = HDS.get_canopy_variable(input_NLCD_raster_url_path=subset_NLCD_result['output_raster'],
                                variable_name='cc', output_netcdf=ccfracNC)
hcanNC = watershedName+'Hcan'+str(dxRes)+'m.nc'
hcan = HDS.get_canopy_variable(input_NLCD_raster_url_path=subset_NLCD_result['output_raster'],
                                variable_name='hcan', output_netcdf=hcanNC)
laiNC = watershedName+'LAI'+str(dxRes)+'m.nc'
lai = HDS.get_canopy_variable(input_NLCD_raster_url_path=subset_NLCD_result['output_raster'],
                                variable_name='lai', output_netcdf=laiNC)

####climate variables
startYear = datetime.strptime(startDateTime,"%Y/%m/%d %H").year
endYear = datetime.strptime(endDateTime,"%Y/%m/%d %H").year
#### we are using data from Daymet; so data are daily
startDate = datetime.strptime(startDateTime, "%Y/%m/%d %H").date().strftime('%m/%d/%Y')
endDate = datetime.strptime(endDateTime, "%Y/%m/%d %H").date().strftime('%m/%d/%Y')

climate_Vars = ['tmin', 'tmax']         #['vp', 'tmin', 'tmax', 'srad', 'prcp']
####iterate through climate variables
for var in climate_Vars:
    for year in range(startYear, endYear+1):
        climatestaticFile1 = var + "_" + str(year) + ".nc4"
        climateFile1 = watershedName + '_' + var + "_" + str(year) + ".nc"
        Year1sub_request = HDS.subset_netcdf(input_netcdf=climatestaticFile1,ref_raster_url_path=Watershed['output_raster'],
                                             output_netcdf=climateFile1)
        concatFile = "conc_"+climateFile1
        if year == startYear:
            concatFile1_url = Year1sub_request['output_netcdf']
        else:
            concatFile2_url = Year1sub_request['output_netcdf']
            concateNC_request = HDS.concatenate_netcdf(input_netcdf1_url_path=concatFile1_url,
                                                       input_netcdf2_url_path=concatFile2_url,
                                                       output_netcdf=concatFile, inout_timeName = 'time')
            concatFile1_url = concateNC_request['output_netcdf']
    climateWYFile = watershedName +'_' +var +str(startYear) +str(endYear) +".nc"
    subset_NC_by_time_result = HDS.subset_netcdf_by_time(input_netcdf_url_path=concatFile1_url,
                                                         time_dimension_name='time', start_date=startDate,
                                                         end_date=endDate, output_netcdf=climateWYFile)
##End for

# NLDAS wind speed and humidity
NLDASVars = ['UGRD10m_110_HTGL', 'VGRD10m_110_HTGL', 'SPFH2m_110_HTGL', 'PRESsfc_110_SFC']  # wind and humidity---but the other variables are listed below
# NLDASVars = ['APCPsfc_110_SFC_acc1h', 'TMP2m_110_HTGL', 'PRESsfc_110_SFC', 'UGRD10m_110_HTGL', 'VGRD10m_110_HTGL', 'SPFH2m_110_HTGL','DLWRFsfc_110_SFC', 'DSWRFsfc_110_SFC']
inTime = 'time'
dTin = 1             # NLDAS data are hourly
ncProj_resample_file_url =[]
for var in NLDASVars:
    for year in range(startYear, endYear+1):
        input_nc_file = 'NLDAS_FORA0125_H.A_' + var + '_' + str(year) + ".nc"
        ouput_nc_file = watershedName + "_" + 'NLDAS_FORA0125_H_A_' + var + '_' + str(year) + ".nc"
        subsetnc = HDS.subset_netcdf_by_coordinates(input_netcdf= input_nc_file,output_netcdf= ouput_nc_file,
                          leftX=leftX, topY=topY, rightX=rightX, bottomY=bottomY, in_Xcoord = 'lon_110', in_Ycoord='lat_110') #, save_as=workingDir+'loganNLDAS_jan2009.nc')
        if year == startYear:
            concatFile1_url = subsetnc['output_netcdf']
        else:
            concatFile2_url = subsetnc['output_netcdf']
            concateNC_request = HDS.concatenate_netcdf(input_netcdf1_url_path=concatFile1_url, input_netcdf2_url_path=concatFile2_url,
                                               output_netcdf=concatFile, inout_timeName=inTime)
            concatFile1_url = concateNC_request['output_netcdf']
    climateWYFile = watershedName + '_' + var + str(startYear) + str(endYear) + ".nc"
    subset_NC_by_time_result = HDS.subset_netcdf_by_datetime(input_netcdf=concatFile1_url, output_netcdf=climateWYFile,
                                       startDateTime= startDateTime, endDateTime= endDateTime, dT= dTin, inout_timeName= inTime, save_as=None)
##End for
#list files
MyFiles = HDS.list_my_files()
for item in MyFiles:
    print(item)

#put files to zip package --can download and upload to HydroShare if needed
tminNC = watershedName +'_tmin' +str(startYear) +str(endYear) +'.nc'
tmaxNC = watershedName +'_tmax' +str(startYear) +str(endYear) +'.nc'
U10NC = watershedName +'_UGRD10m_110_HTGL' +str(startYear) +str(endYear) +'.nc'
V10NC = watershedName +'_VGRD10m_110_HTGL' +str(startYear) +str(endYear) +'.nc'
qNC = watershedName +'_SPFH2m_110_HTGL' +str(startYear) +str(endYear) +'.nc'
PrNC = watershedName +'_PRESsfc_110_SFC' +str(startYear) +str(endYear) +'.nc'

ueb_inputPackage_dict = [domainDEM, aspectRaster, slopeRaster, ccfracNC, hcanNC, laiNC, tminNC, tmaxNC, U10NC, V10NC, qNC, PrNC]
zipFileName = watershedName+str(dxRes)+'_'+str(startYear)+str(endYear)+'.zip'
zip_files_result = HDS.zip_files(files_to_zip=ueb_inputPackage_dict, zip_file_name=zipFileName)
#download to local disk
HDS.download_file(file_url_path=zip_files_result['zip_file_name'], save_as=workingDir+zipFileName)

"""
#### save UEB input package as HydroShare resource
hs_title = 'UEB input package for the '+watershedName+' watershed'
hs_abstract = hs_title +'. It was created by the CI-WATER HydroDS.' + 'Input variables were re-sampled into '+str(dxRes)+ " m grid cells"
hs_keywords=['HydroShare', 'HydroDS', 'DEM', watershedName+' WS']
HDS.set_hydroshare_account(settings.HS_USER_NAME, settings.HS_PASSWORD)
HDS.create_hydroshare_resource(file_name=zipFileName, resource_type ='GenericResource', title= hs_title,
                               abstract= hs_abstract, keywords= hs_keywords)
"""
print('Done')




