"""
This is used to prepare a rectangle area for terrain and land cover data watershed.nc, slope.nc, aspect.nc, cc.nc, hcan.nc, lai.nc
"""
__author__ = 'tsega'

import os
from hydrogate import HydroDS
"""*********** Input settings for watershed of interest *****************"""
workingDir = "/Projects/Tian_workspace/Animas_one_year_rec/"
## Domain bounding box in geographic coordinates left, top, right, bottom.  Must enclose watershed of interest
# use rectangular domain--no WS delineation

# Animas River WS above Durango
leftX, topY, rightX, bottomY =  -108.15, 38.06, -107.41, 37.16
# Green River near Daniel at Warren Bridge
#leftX, topY, rightX, bottomY =  -110.315, 43.493, -109.592, 42.971
# leftX, topY, rightX, bottomY =  -110.415, 43.593, -109.492, 42.871

# Grid projection
#utmZone = int((180 + 0.5*(xmin+xmax))/6) + 1
epsgCode = 26912    #26912                 #26912 utm 12
dx,dy  = 30, 30  #  Grid cell sizes (m) for reprojection
# Set parameters for watershed delineation
streamThreshold = 1000
watershedName = 'Animas'
# Cell spacing for subsampled UEB model (m)
dxRes, dyRes = 1200, 1200
"""*************************************************************************"""
HDS = HydroDS(username='tianG', password='tianGan_2016')

MyFiles = HDS.list_my_files()
for item in HDS.list_my_files():
    try:
        HDS.delete_my_file(item.split('/')[-1])
    except Exception as e:
        print('failed to delete:')
        print(item)
        continue

#### Subset DEM and Delineate Watershed
input_static_DEM  = 'nedWesternUS.tif'
subsetDEM_request = HDS.subset_raster(input_raster=input_static_DEM, left=leftX, top=topY, right=rightX,
                                      bottom=bottomY,output_raster=watershedName + 'DEM84.tif')

#Options for projection with epsg full list at: http://spatialreference.org/ref/epsg/
myWatershedDEM = watershedName + 'Proj' + str(dx) + '.tif'
# you can get epsg code from this website:  http://spatialreference.org/ref/epsg/
WatershedDEM = HDS.project_resample_raster(input_raster_url_path=subsetDEM_request['output_raster'],
                                                      cell_size_dx=dx, cell_size_dy=dy, epsg_code=epsgCode,
                                                      output_raster=myWatershedDEM,resample='bilinear')

####Resample watershed grid to coarser grid
Watershed =  HDS.resample_raster(input_raster_url_path= WatershedDEM['output_raster'],
                cell_size_dx=dxRes, cell_size_dy=dyRes, resample='near', output_raster=watershedName + str(dxRes) + 'm.tif')
HDS.download_file(file_url_path=Watershed['output_raster'], save_as=workingDir+watershedName+str(dxRes)+'.tif')

##  Convert to netCDF for UEB input
Watershed_temp = HDS.raster_to_netcdf(Watershed['output_raster'], output_netcdf='watershed'+str(dxRes)+'.nc')
##  Convert to netCDF for UEB input
#watershedFunctions.wsNetcdf_from_raster(workingDir+watershedName+str(dxRes)+'.tif', workingDir+watershedName+str(dxRes)+'Int.tif',workingDir+watershedName+str(dxRes)+'Int.nc')
#Watershed_temp = HDS.upload_file(workingDir+watershedName+str(dxRes)+'Int.nc')

# In the netCDF file rename the generic variable "Band1" to "watershed"
Watershed_NC = HDS.netcdf_rename_variable(input_netcdf_url_path=Watershed_temp['output_netcdf'],
                                output_netcdf='watershed.nc', input_variable_name='Band1', output_variable_name='watershed')

#terrain variables
# aspect
aspect_hires = HDS.create_raster_aspect(input_raster_url_path=WatershedDEM['output_raster'],
                                output_raster=watershedName + 'Aspect' + str(dx)+ '.tif')
aspect = HDS.resample_raster(input_raster_url_path= aspect_hires['output_raster'], cell_size_dx=dxRes,
                                cell_size_dy=dyRes, resample='near', output_raster=watershedName + 'Aspect' + str(dxRes) + '.tif')
aspect_temp = HDS.raster_to_netcdf(input_raster_url_path=aspect['output_raster'],output_netcdf='aspect'+str(dxRes)+'.nc')
aspect_nc = HDS.netcdf_rename_variable(input_netcdf_url_path=aspect_temp['output_netcdf'],
                                output_netcdf='aspect.nc', input_variable_name='Band1', output_variable_name='aspect')
HDS.download_file(file_url_path=Watershed['output_raster'], save_as=workingDir+watershedName+str(dxRes)+'.tif')
# slope
slope_hires = HDS.create_raster_slope(input_raster_url_path=WatershedDEM['output_raster'],
                                output_raster=watershedName + 'Slope' + str(dx) + '.tif')
slope = HDS.resample_raster(input_raster_url_path= slope_hires['output_raster'], cell_size_dx=dxRes,
                                cell_size_dy=dyRes, resample='near', output_raster=watershedName + 'Slope' + str(dxRes) + '.tif')
slope_temp = HDS.raster_to_netcdf(input_raster_url_path=slope['output_raster'], output_netcdf='slope'+str(dxRes)+'.nc')
slope_nc = HDS.netcdf_rename_variable(input_netcdf_url_path=slope_temp['output_netcdf'],
                                output_netcdf='slope.nc', input_variable_name='Band1', output_variable_name='slope')

#Land cover variables
nlcd_raster_resource = 'nlcd2011CONUS.tif'
subset_NLCD_result = HDS.project_clip_raster(input_raster=nlcd_raster_resource,
                                ref_raster_url_path=Watershed['output_raster'],
                                output_raster=watershedName + 'nlcdProj' + str(dxRes) + '.tif')
#cc
nlcd_variable_result = HDS.get_canopy_variable(input_NLCD_raster_url_path=subset_NLCD_result['output_raster'],
                                variable_name='cc', output_netcdf=watershedName+str(dxRes)+'cc.nc')
cc_nc = HDS.netcdf_rename_variable(input_netcdf_url_path=nlcd_variable_result['output_netcdf'],
                                output_netcdf='cc.nc', input_variable_name='Band1', output_variable_name='cc')
#hcan
nlcd_variable_result = HDS.get_canopy_variable(input_NLCD_raster_url_path=subset_NLCD_result['output_raster'],
                                variable_name='hcan', output_netcdf=watershedName+str(dxRes)+'hcan.nc')
hcan_nc = HDS.netcdf_rename_variable(input_netcdf_url_path=nlcd_variable_result['output_netcdf'],
                                output_netcdf='hcan.nc', input_variable_name='Band1',output_variable_name='hcan')
#lai
nlcd_variable_result = HDS.get_canopy_variable(input_NLCD_raster_url_path=subset_NLCD_result['output_raster'],
                                variable_name='lai', output_netcdf=watershedName+str(dxRes)+'lai.nc')
lai_nc = HDS.netcdf_rename_variable(input_netcdf_url_path=nlcd_variable_result['output_netcdf'],
                                output_netcdf='lai.nc', input_variable_name='Band1',output_variable_name='lai')


# download result
ueb_inputPackage_dict = ['watershed.nc', 'aspect.nc', 'slope.nc', 'cc.nc', 'hcan.nc', 'lai.nc']
zip_files_result = HDS.zip_files(files_to_zip=ueb_inputPackage_dict, zip_file_name='ueb_input.zip')
HDS.download_file(file_url_path=zip_files_result['zip_file_name'], save_as=os.path.join(workingDir,'ueb_input.zip'))

print('Finished UEB input setup')




