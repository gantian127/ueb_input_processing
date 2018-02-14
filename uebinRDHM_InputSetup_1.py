__author__ = 'tsega'
from hydrods_python_client import HydroDS
from datetime import datetime, timedelta
"""*********** Terrain and Land cover for UEB using HydroDS *****************"""

# # user settings to change
workingDir = "/Projects/Tian_workspace/rdhm_ueb_modeling/animas_2007_rec/"

# Grid projection
#utmZone = int((180 + 0.5*(xmin+xmax))/6) + 1
epsgCode = 26912    #26912                 #26912 utm 12
dx,dy = 30, 30  #  Grid cell sizes (m) for reprojection

# Dolores at Dolores
leftX, topY, rightX, bottomY = -108.525, 37.865, -107.85, 37.415  # exact box: -108.51773, 37.857910, -107.863539, 37.428745
watershedN = 'Mcphee_DOLC2'

# # Cell spacing for subsampled UEB model (m)
# dxRes, dyRes = 1200, 1200
## Domain bounding box in geographic coordinates left, top, right, bottom.  Must enclose watershed of interest
# use rectangular domain--no WS delineation

# Green River near Daniel at Warren Bridge
# leftX, topY, rightX, bottomY =  -110.415, 43.593, -109.492, 42.871
# # Animas River WS above Durango
# leftX, topY, rightX, bottomY =  -108.15, 38.06, -107.41, 37.16
# watershedN = 'Animas'

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
input_static_DEM = 'nedWesternUS.tif'
subsetDEM_request = HDS.subset_raster(input_raster=input_static_DEM, left=leftX, top=topY, right=rightX,
                                      bottom=bottomY,output_raster=watershedN + 'DEM84.tif')
#Options for projection with epsg full list at: http://spatialreference.org/ref/epsg/
Watershed30mDEM = watershedN+str(dx)+'m.tif'
WatershedDEM = HDS.project_resample_raster(input_raster_url_path=subsetDEM_request['output_raster'],
                                                      cell_size_dx=dx, cell_size_dy=dy, epsg_code=epsgCode,
                                                      output_raster=Watershed30mDEM,resample='bilinear')
# #Resample watershed grid to coarser grid (to save time)
#
# WatershedResDEM = watershedN + 'DEM' + str(dxRes) + 'm.tif'
# WatershedRes =  HDS.resample_raster(input_raster_url_path= WatershedDEM['output_raster'],
#                 cell_size_dx=dxRes, cell_size_dy=dyRes, resample='near', output_raster=WatershedResDEM)
# ##  Convert to netCDF
# WatershedResNC = watershedN + 'NC' + str(dxRes) + 'm.nc'
# Watershed_resnc = HDS.raster_to_netcdf(WatershedRes['output_raster'], output_netcdf=WatershedResNC)

## use 'referenceRasterTIF' or 'WatershedResNC' for resampling
""" 
## the following xmrg grid output from RDHM can be used as reference raster to resample input to. 
## not used here...high res. 30 m directly input to RDHM 

#proj4_string: see the paper: Reed, S.M., and D.R. Maidment, "Coordinate Transformations for Using NEXRAD Data in GIS-based Hydrologic Modeling," Journal of Hydrologic Engineering, 4, 2, 174-182, April 1999
## http://www.nws.noaa.gov/ohd/hrl/distmodel/hrap.htm#background
proj4_string = 'PROJCS["Sphere_ARC_INFO_Stereographic_North_Pole",GEOGCS["GCS_Sphere_ARC_INFO",DATUM["Sphere_ARC_INFO",SPHEROID["Sphere_ARC_INFO",6370997,0]],PRIMEM["Greenwich",0],\
UNIT["degree",0.0174532925199433]],PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",60.00681388888889],PARAMETER["central_meridian",-105],\
PARAMETER["scale_factor",1],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]]]'
##proj4_string = '+proj=stere +lat_0=90.0 +lat_ts=60.0 +lon_0=-105.0 +k=1 +x_0=0.0 +y_0=0.0 +a=6371200 +b=6371200 +units=m +no_defs'

##reference xmrg file ---must be located in the workingDir
inputXmrgRaster = 'real_uztwc0830200922z.gz'
referenceRasterASCII = watershedN+'_refRaster.asc'
referenceRasterTIF = watershedN+'_refRaster.tif'
referenceRasterNC = watershedN+'_refRaster.nc'
cmdString = 'xmrgtoasc -p ster -i ' + inputXmrgRaster + ' -o ' + referenceRasterASCII
callSubprocess.callSubprocess(cmdString, "xmrg to ascii")
cmdString = "gdal_translate -a_srs \"" + proj4_string + "\" "+referenceRasterASCII+" "+referenceRasterTIF
callSubprocess.callSubprocess(cmdString, "ascii to tif")
watershedFunctions.rasterToNetCDF(referenceRasterTIF, referenceRasterNC)
"""

#terrain variables
print 'aspect'
aspectRaster = watershedN + 'Aspect' + str(dx) + 'm.tif'
aspect = HDS.create_raster_aspect(input_raster_url_path=WatershedDEM['output_raster'],
                                output_raster=aspectRaster)
slopeRaster = watershedN + 'Slope' + str(dx) + 'm.tif'
slope = HDS.create_raster_slope(input_raster_url_path=WatershedDEM['output_raster'],
                                output_raster=slopeRaster)
#Land cover variables
print 'land cover'
nlcd_raster_resource = 'nlcd2011CONUS.tif'
subset_NLCD_result = HDS.project_clip_raster(input_raster=nlcd_raster_resource,
                                ref_raster_url_path=WatershedDEM['output_raster'],
                                output_raster=watershedN + 'nlcdProj' + str(dx) + '.tif')
ccfracNC = watershedN+'CC'+str(dx)+'m.nc'
ccfrac = HDS.get_canopy_variable(input_NLCD_raster_url_path=subset_NLCD_result['output_raster'],
                                variable_name='cc', output_netcdf=ccfracNC)
hcanNC = watershedN+'Hcan'+str(dx)+'m.nc'
hcan = HDS.get_canopy_variable(input_NLCD_raster_url_path=subset_NLCD_result['output_raster'],
                                variable_name='hcan', output_netcdf=hcanNC)
laiNC = watershedN+'LAI'+str(dx)+'m.nc'
lai = HDS.get_canopy_variable(input_NLCD_raster_url_path=subset_NLCD_result['output_raster'],
                                variable_name='lai', output_netcdf=laiNC)

#download to local disk
print 'download'
ueb_inputPackage_dict = [Watershed30mDEM, aspectRaster, slopeRaster, ccfracNC, hcanNC, laiNC]
zipFileName = watershedN+'SiteVars.zip'
zip_files_result = HDS.zip_files(files_to_zip=ueb_inputPackage_dict, zip_file_name=zipFileName)

HDS.download_file(file_url_path=zip_files_result['zip_file_name'], save_as=workingDir+zipFileName)

import zipfile
import os
site_folder = os.path.join(os.getcwd(),'Site_var')
if not os.path.isdir(site_folder):
    os.mkdir(site_folder)

with zipfile.ZipFile(workingDir+zipFileName, 'r') as zip_file:
    zip_file.extractall(site_folder)


print('Done')

