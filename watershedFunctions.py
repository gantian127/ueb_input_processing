"""
This is used for both latest and older version of modeling data
"""
try:
    from osgeo import gdal
except:
    import gdal
try:
    from osgeo import ogr
except:
    import ogr
try:
    from osgeo import osr
except:
    import osr
from gdalconst import *
import shlex
import subprocess
import os
import sys
import math 
import glob
import string
import shutil


#int raster from (possible float) raster (usually dem)
def wsNetcdf_from_raster(input_raster, output_WS_raster, output_WS_netcdf):

    #watershed grid file
    intermRaster = "wsRastfl.tif"
    cmdString = " python \"C:/Python34/Lib/site-packages/osgeo/gdal_calc.py\" -A "+input_raster+" --outfile="+intermRaster+" --calc=A/A"
    callSubprocess(cmdString, "watershed grid computation")

    cmdString = "gdal_translate"+" -ot Int32 " +intermRaster+" "+output_WS_raster
    callSubprocess(cmdString, "convert raster format")

    rasterToNetCDF(output_WS_raster, output_WS_netcdf)
    #To do: Delete temp raster
    #       Free other objects
    #       Direct TauDEM messages to file



def Subset_Raster(input_raster, output_raster, xmin, ymax, xmax, ymin):
    """ Note: upper left (ul) considered origin, i.e. xmin, ymax
    parameters passed as ulx uly lrx lry (xmin, ymax, xmax, ymin)
    The arguments are in decimal degrees with Datum NAD83
    input_raster is in Geographic CS with datum NAD83
    """
    cmdString = "gdal_translate"+" "+"-projwin"+" "+str(xmin)+" "+str(ymax)+" "\
               +str(xmax)+" "+str(ymin)+" "+input_raster+" "+output_raster
    callSubprocess(cmdString, "get DEM subset")


def delineate_Watershed_atShapeFile(input_DEM_raster, input_Outlet_shpFile, output_WS_raster,
                                       output_Outlet_shpFile, streamThreshold):
    """TauDEM doesn't take compressed file; uncompress file
        ToDO:  Check compression first"""
    temp_raster = 'temp.tif'
    uncompressRaster(input_DEM_raster, temp_raster)
    input_raster = os.path.splitext(input_DEM_raster)[0]      #remove the .tif
    # pit remove
    cmdString = " \"C:\Program Files\TauDEM\TauDEM5Exex64\pitremove\" -z "+temp_raster+" -fel "+input_raster+"fel.tif"
    callSubprocess(cmdString,'pit remove')
    #d8 flow dir
    cmdString = " \"C:\Program Files\TauDEM\TauDEM5Exex64\d8flowdir\" -fel "+input_raster+"fel.tif -sd8 "+input_raster+"sd8.tif -p "\
                +input_raster+"p.tif"
    callSubprocess(cmdString, 'd8 flow direction')
    #d8 contributing area without outlet shape file
    cmdString = " \"C:\Program Files\TauDEM\TauDEM5Exex64\\aread8\" -p "+input_raster+"p.tif -ad8 "+input_raster+"ad81.tif -nc"         #check the effect of -nc
    #-o "\ +input_outletshp
    callSubprocess(cmdString, 'd8 contributing area')
    #Get statistics of ad8 file to determine threshold

    #Stream definition by threshold
    cmdString = " \"C:\Program Files\TauDEM\TauDEM5Exex64\\threshold\" -ssa "+input_raster+"ad81.tif -src "+input_raster+"src.tif -thresh "+str(streamThreshold)
    callSubprocess(cmdString, 'Stream definition by threshold')
    #move outlets to stream
    cmdString = " \"C:\Program Files\TauDEM\TauDEM5Exex64\moveoutletstostreams\" -p "+input_raster+"p.tif -src "+input_raster+"src.tif -o "\
                +input_Outlet_shpFile+ " -om "+output_Outlet_shpFile
    callSubprocess(cmdString, 'move outlet to stream')
    #Add projection to moved outlet ---TauDEM excludes the projection from moved outlet; check
    driverName = "ESRI Shapefile"
    driver = ogr.GetDriverByName(driverName)
    dataset = driver.Open(input_Outlet_shpFile)
    layer = dataset.GetLayer()
    srs = layer.GetSpatialRef()
    baseName = os.path.splitext(output_Outlet_shpFile)[0]
    projFile = baseName+".prj"
    #srsString = "+proj=utm +zone="+str(utmZone)+" +ellps=GRS80 +datum=NAD83 +units=m"
    #srs = osr.SpatialReference()
    #srs.ImportFromEPSG(epsgCode)
    #srs.ImportFromProj4(srsString)
    srs.MorphFromESRI()
    file = open(projFile, "w")
    file.write(srs.ExportToWkt())
    file.close()
    #d8 contributing area with outlet shapefile
    cmdString = " \"C:\Program Files\TauDEM\TauDEM5Exex64\\aread8\" -p "+input_raster+"p.tif -ad8 "+input_raster+"ad8.tif -o "+output_Outlet_shpFile+" -nc"
    callSubprocess(cmdString, 'd8 contributing area with outlet shapefile')
    #watershed grid file
    cmdString = " python \"C:/Python34/Lib/site-packages/osgeo/gdal_calc.py\" -A "+input_raster+"ad8.tif --outfile="+output_WS_raster+" --calc=A/A"
    callSubprocess(cmdString, "watershed grid computation")
    #To do: Delete temp raster
    #       Free other objects
    #       Direct TauDEM messages to file



#to get epsgCode, projection code using the EPSG (get list from http://spatialreference.org/ref/epsg/)
def Project_and_resample_Raster_EPSG(input_raster, output_raster,  dx, dy, epsgCode, resample='near'):
    """
    This function projects and re-grids a raster
    parameters are:
        dx, dy are user selected resolution
        epsgCode: projection code using the EPSG (get list from http://spatialreference.org/ref/epsg/)
        resample is regridding/interpolation method
        For images leave the default nearest neighbor interpolation;
        else pass the method required, e.g 'bilinear'
    """
    #Project utm
    cmdString = "gdalwarp -t_srs EPSG:"+str(epsgCode)+" -tr "\
                +str(dx)+" "+str(dy)+" -r "+resample+" -overwrite "+input_raster+" "+output_raster
    callSubprocess(cmdString, "project and re-grid Raster")
    #Delete temp file
    #os.remove("tempRaster.tif")


def project_shapefile_EPSG(input_shp, output_shp, epsgCode):
    """ This projection tested for when the source shape file is in WGS84 Geographic
    coordinate syste (EPSG:4326), but generally gdal/ogr recognizes the input srs
    """
    cmdString = "ogr2ogr -t_srs EPSG:"+str(epsgCode)+" "+output_shp+" "+input_shp
    callSubprocess(cmdString, "Project Shape File")


def Project_and_Resample_Raster_UTMNAD83(input_raster, output_raster,  dx, dy, utmZone, resample='near'):
    """
    This function projects and re-grids a raster
    parameters are:
        dx, dy are user selected resolution
        resample is regridding/interpolation method
        For images leave the default nearest neighbor interpolation;-
        else pass the method required, e.g 'bilinear'
    """
    #Project utm
    cmdString = "gdalwarp -t_srs '+proj=utm +zone=" +str(utmZone)+ " +datum=NAD83' -tr "\
                +str(dx)+" "+str(dy)+" -r "+resample+" -overwrite "+input_raster+" "+output_raster
    callSubprocess(cmdString, "project and re-grid Raster to UTM")
    #Delete temp file
    #os.remove("tempRaster.tif")


def create_OutletShp(shapefileName, outletPointX, outletPointY):
    #Creates a point shapefile (intended for watershed outlet)
    #using point coordinates in WGS84 decimal degrees
    if os.path.isdir(shapefileName):
        shutil.rmtree(shapefileName)
    srs = osr.SpatialReference()
    #srsString = "+proj=latlong +ellps=GRS80 +datum=NAD83"
    #srs.ImportFromProj4(srsString)
    srs.ImportFromEPSG(4326)
    driverName = "ESRI Shapefile"
    drv = ogr.GetDriverByName(driverName )
    if drv is None:
        print("{} driver not available.\n".format(driverName))
        sys.exit( 1 )

    ds = drv.CreateDataSource(shapefileName)   #, 0, 0, 0, gdal.GDT_Unknown )
    if ds is None:
        print("Creation of shapefile failed.\n")
        sys.exit( 1 )

    lyr = ds.CreateLayer(shapefileName, srs, ogr.wkbPoint ) #Add srs
    if lyr is None:
        print("Layer creation failed.\n")
        sys.exit( 1 )

    field_defn = ogr.FieldDefn( "outletName", ogr.OFTString )
    field_defn.SetWidth( 32 )

    if lyr.CreateField ( field_defn ) != 0:
        print ("Creating Name field failed.\n")
        sys.exit( 1 )

#------- add id field for Taudem Gauge watershed
    field_defn = ogr.FieldDefn( "id", ogr.OFTInteger )
    field_defn.SetWidth(32)
    if lyr.CreateField ( field_defn ) != 0:
        print ("Creating Name field failed.\n")
        sys.exit( 1 )
#-------------------------------------------

    x = float(outletPointX)
    y = float(outletPointY)
    name = shapefileName+"Outlet"

    feat = ogr.Feature( lyr.GetLayerDefn())
    feat.SetField( "outletName", name )

    pt = ogr.Geometry(ogr.wkbPoint)
    pt.SetPoint_2D(0, x, y)

    feat.SetGeometry(pt)

    if lyr.CreateFeature(feat) != 0:
        print ("Failed to create feature in shapefile.\n")
        sys.exit( 1 )

    feat.Destroy()

    ds = None


def project_Shapefile_UTM_NAD83(input_shp, output_shp, utmZone):
    """ This projection tested for when the source shape file is in WGS84 Geographic
    coordinate syste (EPSG:4326), but generally gdal/ogr recognizes the input srs
    """
    cmdString = "ogr2ogr -t_srs '+proj=utm +zone=" +str(utmZone)+ " +ellps=GRS80 +datum=NAD83 units=m' "\
              +output_shp+" "+input_shp
    callSubprocess(cmdString, "Project shapefile")


def rasterToNetCDF(input_raster, output_netcdf):
    cmdString = "gdal_translate -of netCDF "+input_raster+" "+output_netcdf
    callSubprocess(cmdString, 'raster to netcdf')


def subset_raster_to_referenceRaster(input_raster, output_raster, reference_Raster):
    """ This function gets part of the input_raster covering the reference_Raster
    To Do: Boundary check-> check if the bounding box of subset raster is
               within the input_raster's boundary
    Boundary parameters extracted from reference_Raster
    """
    data_set = gdal.Open(reference_Raster, GA_ReadOnly)
    geo_transform = data_set.GetGeoTransform()
    # use ulx uly lrx lry
    xmin = geo_transform[0]
    ymax = geo_transform[3]
    dx = geo_transform[1]
    dy = geo_transform[5]
    xmax = xmin + dx * data_set.RasterXSize
    ymin = ymax + dy* data_set.RasterYSize          # dy is -ve

    cmdString = "gdal_translate"+" "+"-projwin"+" "+str(xmin)+" "+str(ymax)+" "\
               +str(xmax)+" "+str(ymin)+" "+input_raster+" "+output_raster
    callSubprocess(cmdString, "subset raster")

    data_set = None


def subset_project_and_resample_Raster_UTMNAD83(input_raster, output_raster, xmin, ymax, xmax, ymin,  dx, dy, resample='near'):
    """
    This function subsets, projects and re-grids a raster
    parameters are:
        (xmin, ymax, xmax, ymin): boundary extents
        in decimal degrees with Datum NAD83
        input_raster is in Geographic CS with datum NAD83
        dx, dy are user selected resolution
        resample is regridding/interpolation method
        For images leave the default nearest neighbor interpolation;
        else pass the method required, e.g 'bilinear'
    """
    cmdString = "gdal_translate"+" "+"-projwin"+" "+str(xmin)+" "+str(ymax)+" "\
               +str(xmax)+" "+str(ymin)+" "+input_raster+" tempRaster.tif"
    callSubprocess(cmdString, "subset raster")
    utmZone = int((180 + 0.5*(xmin+xmax))/6) + 1
    #Project utm
    cmdString = "gdalwarp -t_srs '+proj=utm +zone=" +str(utmZone)+ " +datum=NAD83' -tr "\
                +str(dx)+" "+str(dy)+" -r "+resample+" -overwrite tempRaster.tif "+output_raster
    callSubprocess(cmdString, "project and re-grid DEM")
    #Delete temp file
    os.remove("tempRaster.tif")

def project_raster_UTM_NAD83(input_raster, output_raster, utmZone):
    """ This projection assumes the source spatial reference is known
        i.e. GDAL can read it and recognize it
    """
    cmdString = "gdalwarp -t_srs \"+proj=utm +zone=" +str(utmZone)+ " +ellps=GRS80 +datum=NAD83 +units=m\" "\
                  +input_raster+" "+output_raster
    callSubprocess(cmdString, "project raster")


def resample_Raster(input_raster, output_raster, dx, dy, resample='near'):
    """
    This function re-grids a raster to dx and dy
    parameters are:
          dx, dy are user selected resolution
          resample is regridding/interpolation method
          For images leave the default nearest neighbor interpolation;
          else pass the method required, e.g 'bilinear'
    """
    cmdString = "gdalwarp -tr "+str(dx)+" "+str(dy)+" -et 0.001 -r "+resample+" -overwrite "+\
                 input_raster +" "+ output_raster
    callSubprocess(cmdString, "re-grid raster")
    #Delete temp file


def combineRasters(input_raster1, input_raster2, output_raster):
    """To  Do: may need to specify output no-data value
    """
    cmdString = "gdal_warp "+input_raster1+" "+input_raster2+" "+output_raster
    callSubprocess(cmdString, 'join (stitch) two raster files')


def uncompressRaster(input_raster, output_raster):
    """TauDEM doesn't take compressed file; uncompress file
        ToDO:  Check compression first"""
    cmdString = "gdal_translate -co COMPRESS=NONE"+" "\
               +input_raster+" "+output_raster
    callSubprocess(cmdString, 'uncompress raster')

def download_USGS_NED_and_Subset(output_raster, xmin, ymax, xmax, ymin):
    """
    This function downloads usgs ned dem for user chosen area using usgs web services,
    then it subsets the DEM projects to NAD83 UTM, and re-samples/re-grids.
    parameters are left, top, right, bottom (lon lat coordinates) in decimal degrees (GCS NAD83)
    """
    serviceString = "wget ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/NED/1/IMG/"
    latTop = int(math.ceil(ymax))
    latBottom = int(math.ceil(ymin))
    lonLeft =  -1*(int(math.floor(xmin)))
    lonRight = -1*(int(math.floor(xmax)))
    #check time taken by the following loop
    for lat in range(latBottom, latTop+1):
        for lon in range(lonRight, lonLeft+1):
            if (lon < 100):
                lonStr = "0"+str(lon)
            else:
                lonStr = str(lon)
            zipFile = "n"+str(lat)+"w"+lonStr+".zip"
            cmdString = serviceString+zipFile
            callSubprocess(cmdString, "download USGS DEM using web services")
            #unzip file
            cmdString = "7z e "+ zipFile+ " *.img -r -y"
            callSubprocess(cmdString, "Extract DEM img from Zip file")
    #add to mosaic nedSub temporary file
    cmdString = "gdalwarp -of GTiff -overwrite *.img nedsubTempMosaic.tif"
    callSubprocess(cmdString, "Mosaic images")
    #subset
    cmdString = "gdal_translate"+" "+"-projwin"+" "+str(xmin)+" "+str(ymax)+" "\
               +str(xmax)+" "+str(ymin)+" nedsubTempMosaic.tif tempRaster.tif"
    callSubprocess(cmdString, "subset USGS DEM")


def download_USGSNED_subset_resample_and_Project(output_raster, xmin, ymax, xmax, ymin, dx, dy, resample='near'):
    """
    This function downloads usgs ned dem for user chosen area using usgs web services,
    then it subsets the DEM projects to NAD83 UTM, and re-samples/re-grids.
    parameters are left, top, right, bottom (lon lat coordinates) in decimal degrees (GCS NAD83)
    dx, dy are user selected resolution
    resample is regridding/interpolation method
    For images leave the default nearest neighbor interpolation; else pass the method required, e.g 'bilinear'
    """
    serviceString = "wget ftp://rockyftp.cr.usgs.gov/vdelivery/Datasets/Staged/NED/1/IMG/"
    latTop = int(math.ceil(ymax))
    latBottom = int(math.ceil(ymin))
    lonLeft =  -1*(int(math.floor(xmin)))
    lonRight = -1*(int(math.floor(xmax)))
    #check time taken by the following loop
    for lat in range(latBottom, latTop+1):
        for lon in range(lonRight, lonLeft+1):
            if (lon < 100):
                lonStr = "0"+str(lon)
            else:
                lonStr = str(lon)
            zipFile = "n"+str(lat)+"w"+lonStr+".zip"
            cmdString = serviceString+zipFile
            callSubprocess(cmdString, "download USGS DEM using web services")
            #unzip file
            cmdString = "7z e "+ zipFile+ " *.img -r -y"
            callSubprocess(cmdString, "Extract DEM img from Zip file")
    #add to mosaic nedSub temporary file
    cmdString = "gdalwarp -of GTiff -overwrite *.img nedsubTempMosaic.tif"
    callSubprocess(cmdString, "Mosaic images")
    #subset
    cmdString = "gdal_translate"+" "+"-projwin"+" "+str(xmin)+" "+str(ymax)+" "\
               +str(xmax)+" "+str(ymin)+" nedsubTempMosaic.tif tempRaster.tif"
    callSubprocess(cmdString, "subset USGS DEM")
    utmZone = int((180 + 0.5*(xmin+xmax))/6) + 1
    #Project utm
    cmdString = "gdalwarp -t_srs '+proj=utm +zone=" +str(utmZone)+ " +datum=NAD83' -tr "\
                +str(dx)+" "+str(dy)+" -r "+resample+" -overwrite tempRaster.tif "+output_raster
    callSubprocess(cmdString, "project and re-grid DEM")
    #Delete temp files
    os.remove("tempRaster.tif")
    os.remove("nedsubTempMosaic.tif")
    #delte the individual .img files
    for fl in glob.glob("*.img"):
        os.remove(fl)
    #   Delete the zip files ? ??OR make record so they are not downloaded next time
    # Comment out the following two lines to keep the downloaded zip files for future use
    for fl in glob.glob("*.zip"):
        os.remove(fl)


def callSubprocess(cmdString, debugString):
    cmdargs = shlex.split(cmdString)
    debFile = open('debugFile.txt', 'w')
    debFile.write('Starting %s \n' % debugString)
    retValue = subprocess.call(cmdargs,stdout=debFile)
    if (retValue==0):
        debFile.write('%s Successful\n' % debugString)
        debFile.close()
    else:
        debFile.write('There was error in %s\n' % debugString)
        debFile.close()