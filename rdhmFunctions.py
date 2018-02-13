"""
This is used for the latest version of model input preparation
"""
__author__ = 'tsega'

try:
    from osgeo import gdal, osr, ogr
except:
    import gdal, osr, ogr
from gdalconst import *
import shlex
import subprocess
import os
import numpy as np
import netCDF4
import math
from scipy import interpolate
from datetime import datetime, timedelta
import sys
import shutil
import os

#for rdhm params use xmrg file from SWITmm
def call_subset_Iteratively(reference_raster, fileSuffix='.asc',subdir='sub'):

    for i in os.listdir(os.getcwd()):
        if i.endswith(fileSuffix):
            output_raster = i.replace(fileSuffix,".tif")
            subset_and_resample_Raster(i,reference_raster,output_raster)
            cmdString = "gdal_translate -of AAIGrid "+output_raster+" "+subdir+"/"+i
            callSubprocess(cmdString, 'save ascii')
            #delte output_raster

def subset_and_resample_Raster(input_raster,  reference_raster, output_raster, driverName='GTiff'):  #, output_array):

    srs_data = gdal.Open(input_raster, GA_ReadOnly)
    srs_proj = srs_data.GetProjection() #osr.SpatialReference(wkt
    srs_geotrs = srs_data.GetGeoTransform()

    ref_data = gdal.Open(reference_raster, GA_ReadOnly)
    ref_proj = ref_data.GetProjection()
    ref_geotrs = ref_data.GetGeoTransform()
    Ncols = ref_data.RasterXSize
    Nrows = ref_data.RasterYSize
    ref_data = None

    out_data = gdal.GetDriverByName(driverName).Create(output_raster, Ncols, Nrows, 1, gdal.GDT_Float32)
    out_data.SetGeoTransform(ref_geotrs)
    out_data.SetProjection(ref_proj)

    gdal.ReprojectImage(srs_data,out_data,ref_proj,ref_proj, gdal.GRA_Bilinear )

    srs_data = None
    out_data = None



# dupblicate asci file
def dublicate_ascFile(input_file, filePrefix, time_length, startDateTime='2008/10/01 6', dThour = 6):
    """
    :param filePrefix:
    :param time_length:
    :param startDateTime:
    :param dThour:
    :return:
    """
    for istep in range (int(time_length)):
        curTime = datetime.strptime(startDateTime,'%Y/%m/%d %H')+timedelta(hours=istep*dThour)
        monthStr = str(curTime.month)
        dayStr = str(curTime.day)
        hourStr = str(curTime.hour)
        if curTime.month < 10: monthStr = '0' + str(curTime.month)
        if curTime.day < 10: dayStr = '0' + str(curTime.day)
        if curTime.hour < 10: hourStr = '0'+str(curTime.hour)
        outfileName = filePrefix+monthStr+dayStr+str(curTime.year)+hourStr+"z.asc"
        shutil.copyfile(input_file,outfileName)


# 9.16.17 for netCDF with ensemble of tseries
# creates multiple ascii grid files from 3D netCDF; projection stere recognized in asctoxmrg to convert the grid files to xmrg
# creates hourly data from (possibly) multi-hour input
def create_multiple_Ascii_fromNetCDF_Ensemble(input_netCDF, output_prefix, input_varName, time_length, ensNumber,
                                              dThour=3,
                                              startDateTime='2008/10/01 0', time_varName='time',
                                              ensemble_varName='ensembleNumber',
                                              proj4_string='+proj=stere +lat_0=90.0 +lat_ts=60.00 +lon_0=-105 +k_0=1 +x_0=401.0 +y_0=1601.0 +ellps=sphere +a=6370997 +b=6370997  +units=m +no_defs'):
    # >> proj4_string this is from the paper by David Maidment etal
    """
    :param input_netCDF:  should be in format (time, y, x) or (time, lat, lon)
    :param output_netCDF:
    :param input_varName:
    :param time_length:
    :param time_varName:
    :return:
    """
    tempNetCDF = "temp_" + input_netCDF
    pdqNetCDF = "pdq_" + input_netCDF
    tempTiff1 = "temp2_" + input_netCDF + ".tif"
    for istep in range(int(time_length)):  # "+ensemble_varName+"
        cmdString = "ncea -O -d " + time_varName + "," + str(istep) + "," + str(istep) + " -d " + ensemble_varName + "," + str(ensNumber) + "," + str(ensNumber) + " " + input_netCDF + " " + tempNetCDF
        callSubprocess(cmdString, 'subset netCDF in time dim and ens number')
        # pdq
        cmdString = "ncpdq -O -a " + time_varName + "," + str(ensemble_varName) + ",y,x " + tempNetCDF + " " + pdqNetCDF
        callSubprocess(cmdString, 'permute dimensions')
        # project to stereographic
        cmdString = "gdalwarp -t_srs \"" + proj4_string + "\" NETCDF:\"" + pdqNetCDF + "\":" + input_varName + " " + tempTiff1
        callSubprocess(cmdString, 'project inermediate file')

        for intermT in range(int(dThour)):  # for time steps > 1 hr uniformly interpolate
            curTime = datetime.strptime(startDateTime, '%Y/%m/%d %H') + timedelta(hours=istep * dThour + intermT)
            monthStr = str(curTime.month)
            dayStr = str(curTime.day)
            hourStr = str(curTime.hour)
            if curTime.month < 10: monthStr = '0' + str(curTime.month)
            if curTime.day < 10: dayStr = '0' + str(curTime.day)
            if curTime.hour < 10: hourStr = '0' + str(curTime.hour)
            cmdString = "gdal_translate -of AAIGrid  " + tempTiff1 + " " \
                        + output_prefix + monthStr + dayStr + str(curTime.year) + hourStr + "z.asc"
            callSubprocess(cmdString, 'convert to asc files')


# creates multiple ascii grid files from 3D netCDF; projection stere recognized in asctoxmrg to convert the grid files to xmrg
#note: it outputs hourly data
def create_multiple_Ascii_fromNetCDF(input_netCDF, output_prefix, input_varName, dThour = 3, startDateTime='2008/10/01 0', time_varName='time', in_epsgCode=None,
        proj4_string = '+proj=stere +lat_0=90.0 +lat_ts=37.60 +lon_0=-105 +k_0=1 +x_0=0 +y_0=0 +ellps=sphere +a=6371200 +b=6371200  +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'):

    """
    :param input_netCDF:  should be in format (time, y, x) or (time, lat, lon)
    :param input_varName:
    :param time_length:
    :param time_varName:
    :return:
    """
    ncIn = netCDF4.Dataset(input_netCDF, "r")  # format='NETCDF4')
    time_length = len(ncIn.dimensions[time_varName])
    ncIn.close()

    tempNetCDF = "temp_"+input_netCDF
    tempTiff1 = "temp2_"+input_netCDF+".tif"
    for istep in range (int(time_length)):
        cmdString = "ncea -O -d "+time_varName+","+str(istep)+","+str(istep)+" "+input_netCDF+" "+tempNetCDF
        callSubprocess(cmdString, 'subset netCDF in time dim')
        # project to stereographic
        if in_epsgCode == None:
            cmdString = "gdalwarp -t_srs \"" + proj4_string + "\" NETCDF:\"" + tempNetCDF + "\":" + input_varName + " " + tempTiff1
        else:
            cmdString = "gdalwarp -s_srs EPSG:"+str(in_epsgCode)+" -t_srs \""+proj4_string+"\" NETCDF:\""+ tempNetCDF+"\":"+input_varName+" "+tempTiff1

        callSubprocess(cmdString, 'project inermediate file')

        for intermT in range(int(dThour)):           # for time steps > 1 hr uniformly interpolate
            curTime = datetime.strptime(startDateTime,'%Y/%m/%d %H')+timedelta(hours=istep*dThour+intermT)
            monthStr = str(curTime.month)
            dayStr = str(curTime.day)
            hourStr = str(curTime.hour)
            if curTime.month < 10: monthStr = '0' + str(curTime.month)
            if curTime.day < 10: dayStr = '0' + str(curTime.day)
            if curTime.hour < 10: hourStr = '0'+str(curTime.hour)
            cmdString = "gdal_translate -of AAIGrid  "+tempTiff1+" "\
                        +output_prefix+monthStr+dayStr+str(curTime.year)+hourStr+"z.asc"
            callSubprocess(cmdString, 'convert to asc files')



#leftX=-112.0, topY=42.3, rightX=-111.0, bottomY=41.6, startYear=2009, endYear=2010
# for cbrfc forcing
"""
def create_multiple_ascii_multple_nc(file_prefix, wsName, output_prefix,  var_name, time_length, dThour = 3, time_varName='time',
                startDateTime='2008/10/01 0', endDateTime='2009/10/01 0',
        proj4_string = '+proj=stere +lat_0=90.0 +lat_ts=37.60 +lon_0=-105 +k_0=1 +x_0=0 +y_0=0 +ellps=sphere +a=6371200 +b=6371200  +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'):
"""
    #Subsets and combines multiple netcdf files
    #written for CBRFC temp/prec data sets
    #for datasets that are commpressed with .gz
"""
    startYear = datetime.strptime(startDateTime,"%Y/%m/%d %H").year
    endYear = datetime.strptime(endDateTime,"%Y/%m/%d %H").year
    startDay =  datetime.strptime(startDateTime,"%Y/%m/%d %H").timetuple().tm_yday        #start date = day of year e.g. for 270 for oct 1 2010
    endDay   =  startDay + (datetime.strptime(endDateTime,"%Y/%m/%d %H") - datetime.strptime(startDateTime,"%Y/%m/%d %H")).days          # end date = number of days b/n end-start + start


    for year in range(startYear, endYear+1): 155G  92% /Projects

        cmdString = 'for i in '+file_prefix+'_'+str(year)+'*.nc.gz; do gunzip -cdfq $i>${i//.gz}; done'     #  && rm -f ${i//.gz}  delete the .nc immediately to save space
        #print(cmdString)
        callSubprocess(cmdString, 'unpack nc files for year'+str(year))
    #cmdString="cp mm_20100731.nc.gz mm_20100731.nc.gz_zcopy"

    #cmdString = "for i in "+file_prefix+"_"+str(year)+"*.nc.gz; do gzip -dfq $i; done"      #+subdir+"\/" //&& rm -f ${i//.gz}   delete the .nc immediately to save space
    #callSubprocess(cmdString, 'subset nc files for year '+str(endYear))

    for istep in range (int(time_length)):
        cmdString = "ncea -O -d "+time_varName+","+str(istep)+","+str(istep)+" "+input_netCDF+" "+tempNetCDF
        callSubprocess(cmdString, 'subset netCDF in time dim')
        # project to stereographic
        cmdString = "gdalwarp -t_srs \""+proj4_string+"\" NETCDF:\""+ tempNetCDF+"\":"+input_varName+" "+tempTiff1
        callSubprocess(cmdString, 'project inermediate file')

        for intermT in range(int(dThour)):           # for time steps > 1 hr uniformly interpolate
            curTime = datetime.strptime(startDateTime,'%Y/%m/%d %H')+timedelta(hours=istep*dThour+intermT)
            monthStr = str(curTime.month)
            dayStr = str(curTime.day)
            hourStr = str(curTime.hour)
            if curTime.month < 10: monthStr = '0' + str(curTime.month)
            if curTime.day < 10: dayStr = '0' + str(curTime.day)
            if curTime.hour < 10: hourStr = '0'+str(curTime.hour)
            cmdString = "gdal_translate -of AAIGrid  "+tempTiff1+" "\
                        +output_prefix+monthStr+dayStr+str(curTime.year)+hourStr+"z.asc"
            callSubprocess(cmdString, 'convert to asc files')


    #delete intermediate files
    cmdString = "rm -f R_*.nc "
    callSubprocess(cmdString, "delete intermediate files")
    cmdString = "rm -f "+wsName+"_"+var_name+"*.nc"
    callSubprocess(cmdString, "delete intermediate files")
    #os.remove("R_*.nc")
"""


def sterographic_toHrapX(in_sterX):
    """
    :param in_sterX: in m
    :return:
    from http://www.nws.noaa.gov/oh/hrl/distmodel/hrap.htm
    xster=hrapx*4762.5 - 401*4762.5
    yster=hrapy*4762.5-1601*4762.5
    """
    return (in_sterX/4762.5) + 401

def sterographic_toHrapY(in_sterY):
    """
    :param in_sterY:
    :return:
    from http://www.nws.noaa.gov/oh/hrl/distmodel/hrap.htm
    xster=hrapx*4762.5 - 401*4762.5
    yster=hrapy*4762.5-1601*4762.5
    """
    return (in_sterY/4762.5) + 1601


def checkGeos():

    point = """{"type":"Point","coordinates":[-90.02001032734194,35.127405562848075]}"""
    geom = ogr.CreateGeometryFromJson(point)
    print(geom.GetX())
    print(geom.GetY())


def proj_Ster(input_DEM_raster, project_raster,
         proj4_string = '+proj=stere +lat_0=90.0 +lat_ts=37.60 +lon_0=-105 +k=1 +x_0=0 +y_0=0 +ellps=sphere +a=6371200 +b=6371200  +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'):
    """
    #60.00681402
    creates RDHM connectivity file from input raster and mask watershed raster file
    projections are done inside the function
    outlet shape file must be in the same coordinate system as d8flow dir raster
    :param in_outletShp:
    :param out_connFile:
    :param in_AreaunitConv: Area is in mi^2 hence need unit coverstion e.g m^2 to sq.miles 0.000000386102
    standard parallel use 60.00681402 if radius is 6371200
    :return:
    """
    """TauDEM doesn't take compressed file; uncompress file
        ToDO:  Check compression first"""
    uncompress_raster = 'uncompress_'+input_DEM_raster
    uncompressRaster(input_DEM_raster, uncompress_raster)
    #project raster
    cmdString = "gdalwarp -t_srs \""+proj4_string+"\" "+uncompress_raster+" "+project_raster
    callSubprocess(cmdString, 'project temp file')
    #project watershed file
    data_set = gdal.Open(project_raster, GA_ReadOnly)
    in_proj = data_set.GetProjection()
    geo_transform = data_set.GetGeoTransform()
    # ulx uly lrx lry
    xmin = geo_transform[0]
    ymax = geo_transform[3]
    dx = geo_transform[1]
    dy = geo_transform[5]
    Ncols = data_set.RasterXSize
    Nrows = data_set.RasterYSize

    print('dx=%d dy=%d'%(dx, dy))
    print('xmin=%f ymax=%f'%(xmin, ymax))
    xmax = xmin + dx * (Ncols-1)         # per convention in the rdhm manual
    ymin = ymax + dy* Nrows              # dy is -ve
    ymaxCon = ymax + dy                  # per convention in the rdhm manual

    srsband = data_set.GetRasterBand(1)
    inarray = srsband.ReadAsArray()
    dType = srsband.DataType
    nodata = srsband.GetNoDataValue()
    print('no data: %d'%nodata)
    arr = np.ma.masked_equal(inarray,nodata)
    print(arr)

    #dataSource = None
    data_set =None


#####ToDO: replace exisintg files
def create_RDHM_connectivity_rawRaster(input_DEM_raster, input_Outlet_shpFile, input_WS_raster, streamThreshold, out_movedOutlet_shpFile, out_connFile, in_hrapDXY = 0.1667, in_AreaunitConv = 0.000000386102,
         proj4_string = '+proj=stere +lat_0=90.0 +lat_ts=37.60 +lon_0=-105 +k=1 +x_0=0 +y_0=0 +ellps=sphere +a=6371200 +b=6371200  +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'):
    """
    #60.00681402
    creates RDHM connectivity file from input raster and mask watershed raster file
    projections are done inside the function
    outlet shape file must be in the same coordinate system as d8flow dir raster
    :param in_outletShp:
    :param out_connFile:
    :param in_AreaunitConv: Area is in mi^2 hence need unit coverstion e.g m^2 to sq.miles 0.000000386102
    standard parallel use 60.00681402 if radius is 6371200
    :return:
    """
    """TauDEM doesn't take compressed file; uncompress file
        ToDO:  Check compression first"""
    uncompress_raster = 'uncompress_'+input_DEM_raster
    uncompressRaster(input_DEM_raster, uncompress_raster)
    #project raster
    project_raster = 'proj_'+input_DEM_raster
    cmdString = "gdalwarp -t_srs \""+proj4_string+"\" "+uncompress_raster+" "+project_raster
    callSubprocess(cmdString, 'project temp file')
    #project watershed file
    proj_WS_raster = "proj_"+input_WS_raster
    cmdString = "gdalwarp -t_srs \""+proj4_string+"\" "+input_WS_raster+" "+proj_WS_raster
    callSubprocess(cmdString, 'project temp file')
    #project shapefile
    prject_shapeFile = 'proj_'+input_Outlet_shpFile
    cmdString = "ogr2ogr -t_srs \""+proj4_string+"\" "+prject_shapeFile+" "+input_Outlet_shpFile
    callSubprocess(cmdString, "Project Shape File")
    #mask raster to watershed (this is required because the cells in connectivity file have to be only 'active' cells)
    mask_raster = 'mask1_'+input_DEM_raster
    subset_and_resample_Raster(project_raster,proj_WS_raster,mask_raster)
    #mask_raster = 'mask_'+input_DEM_raster
    #cmdString = " python \"C:/Python34/Lib/site-packages/osgeo/gdal_calc.py\" -A "+mask_raster1+" -B "+proj_WS_raster+" --outfile="+mask_raster+" --calc=A*B"
    #callSubprocess(cmdString, "watershed grid computation")
    # pit remove
    input_raster = os.path.splitext(input_DEM_raster)[0]      #remove the .tif
    cmdString = "pitremove -z "+mask_raster+" -fel "+input_raster+"fel.tif"
    callSubprocess(cmdString,'pit remove')
    #d8 flow dir
    cmdString = "d8flowdir -fel "+input_raster+"fel.tif -sd8 "+input_raster+"sd8.tif -p "\
                +input_raster+"p.tif"
    callSubprocess(cmdString, 'd8 flow direction')
    #d8 contributing area without outlet shape file
    cmdString = "aread8 -p "+input_raster+"p.tif -ad8 "+input_raster+"ad81.tif -nc"         #check the effect of -nc
    #-o "\ +input_outletshp
    callSubprocess(cmdString, 'd8 contributing area')
    #Get statistics of ad8 file to determine threshold
    #Stream definition by threshold
    cmdString = "threshold -ssa "+input_raster+"ad81.tif -src "+input_raster+"src.tif -thresh "+str(streamThreshold)
    callSubprocess(cmdString, 'Stream definition by threshold')
    #move outlets to stream
    mov_shapeFile = 'mov_'+input_Outlet_shpFile
    cmdString = "moveoutletstostrm -p "+input_raster+"p.tif -src "+input_raster+"src.tif -o "\
                +prject_shapeFile+ " -om "+mov_shapeFile
    callSubprocess(cmdString, 'move outlet to stream')

    """
    on windows may need full path to TauDEM
    # pit remove
    input_raster = os.path.splitext(input_DEM_raster)[0]      #remove the .tif
    cmdString = " \"C:\Program Files\TauDEM\TauDEM5Exex64\pitremove\" -z "+mask_raster+" -fel "+input_raster+"fel.tif"
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
    mov_shapeFile = 'mov_'+input_Outlet_shpFile
    cmdString = " \"C:\Program Files\TauDEM\TauDEM5Exex64\moveoutletstostreams\" -p "+input_raster+"p.tif -src "+input_raster+"src.tif -o "\
                +prject_shapeFile+ " -om "+mov_shapeFile
    callSubprocess(cmdString, 'move outlet to stream')
    """
    #Add projection to moved outlet ---TauDEM excludes the projection from moved outlet; check
    driverName = "ESRI Shapefile"
    driver = ogr.GetDriverByName(driverName )
    if driver is None:
        print("{} driver not available.\n".format(driverName))
        sys.exit( 1 )
    dataset = driver.Open(prject_shapeFile)
    layer = dataset.GetLayer()
    srs = layer.GetSpatialRef()
    baseName = os.path.splitext(mov_shapeFile)[0]
    projFile = baseName+".prj"
    srs.MorphFromESRI()
    file = open(projFile, "w")
    file.write(srs.ExportToWkt())
    file.close()
    dataset = None
    #still the proj file is not being recognized--reprojected here
    cmdString = "ogr2ogr -t_srs \""+proj4_string+"\" "+out_movedOutlet_shpFile+" "+mov_shapeFile
    callSubprocess(cmdString, "Project Shape File")

    # outlets in connectivity
    dataSource = driver.Open(out_movedOutlet_shpFile, 0)
    layer = dataSource.GetLayer()
    fCount  = layer.GetFeatureCount()
    print('feature count: %d \n'%fCount)

    in_d8flowDirfile = input_raster+"p.tif"
    data_set = gdal.Open(in_d8flowDirfile, GA_ReadOnly)
    in_proj = data_set.GetProjection()
    geo_transform = data_set.GetGeoTransform()
    # ulx uly lrx lry
    xmin = geo_transform[0]
    ymax = geo_transform[3]
    dx = geo_transform[1]
    dy = geo_transform[5]
    Ncols = data_set.RasterXSize
    Nrows = data_set.RasterYSize

    print('dx=%d dy=%d'%(dx, dy))
    xmax = xmin + dx * (Ncols-1)         # per convention in the rdhm manual
    ymin = ymax + dy* Nrows              # dy is -ve
    ymaxCon = ymax + dy                  # per convention in the rdhm manual

    srsband = data_set.GetRasterBand(1)
    inarray = srsband.ReadAsArray()
    dType = srsband.DataType
    nodata = srsband.GetNoDataValue()
    print('no data: %f'%nodata)
    arr = np.ma.masked_equal(inarray,nodata)

    maskArray = np.empty((Nrows,Ncols),dtype=int)
    maskArray.fill(-999)
    numId = 0
    for ir in range (Nrows):
        for jc in range (Ncols):
            if inarray[ir,jc] != nodata:             #included only masked area; The DEM should be no-data for those outside of mask (e.g. watershed)
                maskArray[ir,jc] = numId
                numId +=1

    conFile = open(out_connFile, 'w')
    conFile.write('TEXT_SEQ \n')
    conFile.write('NUM_HEADER_REC %d %d \n' %(fCount,numId))
    conFile.write('COL %d \n' %Ncols)
    conFile.write('ROW %d \n' %Nrows)
    conFile.write('LLX %f \n' %sterographic_toHrapX(xmin))
    conFile.write('LLY %f \n' %sterographic_toHrapY(ymin))
    conFile.write('URX %f \n' %sterographic_toHrapX(xmax))
    conFile.write('URY %f \n' %sterographic_toHrapY(ymax))
    conFile.write('DXY %f \n' %in_hrapDXY)
    conFile.write('DATA_HRAP \n')

    for feature in layer:
        geom = feature.GetGeometryRef()
        outName = feature.GetField("outletName")
        #gCentr = geom.Centroid()
        outX = geom.GetX()
        outY = geom.GetY()
        print('outX: %f '%outX)
        print('outY: %f '%outY)
        outIx = int((outX - xmin) / dx)
        outJy = int((ymax - outY)/abs(dy))         # dy is -ve
        """
        print(outJy)
        print(outIx)
        print(outX)
        print(xmin)
        print(ymax)
        print(outY)
        """
        if inarray[outJy,outIx] < 1 or inarray[outJy,outIx] > 8 or inarray[outJy,outIx] ==nodata:
            print("outlet location not in the watershed domain")

        cellId = maskArray[outJy,outIx]   #outJy*Nrows + outIx
        conFile.write('%s    1   %d   %f   %f   %f \n' %(outName,cellId, abs(dx*dy)*in_AreaunitConv, sterographic_toHrapX(outX), sterographic_toHrapY(outY)))  # distances are in units of Hrap
    for ir in range (Nrows):
        for jc in range (Ncols):
            if inarray[ir,jc] == 1:
                dcx = 1
                dcy = 0
            elif inarray[ir,jc] == 2:
                dcx = 1
                dcy = -1
            elif inarray[ir,jc] == 3:
                dcx = 0
                dcy = -1
            elif inarray[ir,jc] == 4:
                dcx = -1
                dcy = -1
            elif inarray[ir,jc] == 5:
                dcx = -1
                dcy = 0
            elif inarray[ir,jc] == 6:
                dcx = -1
                dcy = 1
            elif inarray[ir,jc] == 7:
                dcx = 0
                dcy = 1
            elif inarray[ir,jc] == 8:
                dcx = 1
                dcy = 1
            else:                 # no flow dir? this possibly happens only for the outlet cell
                dcx = 0
                dcy = 0
            if inarray[ir,jc] != nodata:             #included only masked area; The DEM should be no-data for those outside of mask (e.g. watershed)
                cellId = maskArray[ir,jc]          #ir*Nrows + jc
                dcellId = maskArray[ir+dcy, jc+dcx]   #(ir+dcy)*Nrows + jc+dcx
                if inarray[ir+dcy, jc+dcx] == nodata and dcellId >0:         # flowing out of region, make sure it is -ve value
                    dcellId = -1*dcellId
                if(dcx==0 and dcy==0) and dcellId >0:                      #outlet
                    dcellId = -1*dcellId
                xcoord = xmin + dx * jc     # + 0.5*dx            #following the conn. file convention, left border of cell
                ycoord = ymax + dy*(ir+1)   # - 0.5*dy            # bottom border
                conFile.write('     %d  %d      Rv   1      %d      %d      %f      %f      %f\n' %(cellId,dcellId,ir,jc, abs(dx*dy)*in_AreaunitConv, sterographic_toHrapX(xcoord), sterographic_toHrapY(ycoord)) )  #distances are in units of Hrap

    conFile.close()
    dataSource = None
    data_set =None


def create_RDHM_connectivity(in_d8flowDirfile, in_outletShp, in_hrapDXY, in_AreaunitConv, out_connFile):
    """
    creates RDHM connectivity file from d8 flow direction raster
    the d8flowdir raster must be mask-subset to exclude the cells outside of basin of interest
    outlset shape file must be in the same coordinate system as d8flow dir raster
    :param in_d8flowDirfile:
    :param in_outletShp:
    :param out_connFile:
    :param in_AreaunitConv: Area is in mi^2 hence need unit coverstion e.g m^2 to sq.miles 0.000000386102
    :return:
    """
    driverName = "ESRI Shapefile"
    drv = ogr.GetDriverByName(driverName )
    if drv is None:
        print("{} driver not available.\n".format(driverName))
        sys.exit( 1 )

    dataSource = drv.Open(in_outletShp, 0)
    layer = dataSource.GetLayer()
    fCount  = layer.GetFeatureCount()
    print('feature count: %d \n'%fCount)

    data_set = gdal.Open(in_d8flowDirfile, GA_ReadOnly)
    in_proj = data_set.GetProjection()
    geo_transform = data_set.GetGeoTransform()
    # ulx uly lrx lry
    xmin = geo_transform[0]
    ymax = geo_transform[3]
    dx = geo_transform[1]
    dy = geo_transform[5]
    Ncols = data_set.RasterXSize
    Nrows = data_set.RasterYSize

    print('dx=%d dy=%d'%(dx, dy))
    xmax = xmin + dx * (Ncols-1)         # per convention in the rdhm manual
    ymin = ymax + dy* Nrows              # dy is -ve
    ymaxCon = ymax + dy                  # per convention in the rdhm manual

    srsband = data_set.GetRasterBand(1)
    inarray = srsband.ReadAsArray()
    dType = srsband.DataType
    nodata = srsband.GetNoDataValue()
    print('no data: %f'%nodata)
    arr = np.ma.masked_equal(inarray,nodata)

    maskArray = np.empty((Nrows,Ncols),dtype=int)
    maskArray.fill(-999)
    numId = 0
    for ir in range (Nrows):
        for jc in range (Ncols):
            if inarray[ir,jc] != nodata:             #included only masked area; The DEM should be no-data for those outside of mask (e.g. watershed)
                maskArray[ir,jc] = numId
                numId +=1

    conFile = open(out_connFile, 'w')
    conFile.write('TEXT_SEQ \n')
    conFile.write('NUM_HEADER_REC %d %d \n' %(fCount,numId))
    conFile.write('COL %d \n' %Ncols)
    conFile.write('ROW %d \n' %Nrows)
    conFile.write('LLX %f \n' %sterographic_toHrapX(xmin))
    conFile.write('LLY %f \n' %sterographic_toHrapY(ymin))
    conFile.write('URX %f \n' %sterographic_toHrapX(xmax))
    conFile.write('URY %f \n' %sterographic_toHrapY(ymax))
    conFile.write('DXY %f \n' %in_hrapDXY)
    conFile.write('DATA_HRAP \n')

    for feature in layer:
        geom = feature.GetGeometryRef()
        outName = feature.GetField("outletName")
        gCentr = geom.Centroid()
        outX = gCentr.GetX()
        outY = gCentr.GetY()
        print('outX: %f '%outX)
        print('outY: %f '%outY)
        outIx = int((outX - xmin) / dx)
        outJy = int((ymax - outY)/abs(dy))         # dy is -ve
        """
        print(outJy)
        print(outIx)
        print(outX)
        print(xmin)
        print(ymax)
        print(outY)
        """
        cellId = maskArray[outJy,outIx]   #outJy*Nrows + outIx
        conFile.write('%s    1   %d   %f   %f   %f \n' %(outName,cellId, abs(dx*dy)*in_AreaunitConv, sterographic_toHrapX(outX), sterographic_toHrapY(outY)))  # distances are in units of Hrap
    for ir in range (Nrows):
        for jc in range (Ncols):
            if inarray[ir,jc] == 1:
                dcx = 1
                dcy = 0
            elif inarray[ir,jc] == 2:
                dcx = 1
                dcy = -1
            elif inarray[ir,jc] == 3:
                dcx = 0
                dcy = -1
            elif inarray[ir,jc] == 4:
                dcx = -1
                dcy = -1
            elif inarray[ir,jc] == 5:
                dcx = -1
                dcy = 0
            elif inarray[ir,jc] == 6:
                dcx = -1
                dcy = 1
            elif inarray[ir,jc] == 7:
                dcx = 0
                dcy = 1
            elif inarray[ir,jc] == 8:
                dcx = 1
                dcy = 1
            else:                 # no flow dir? this possibly happens only for the outlet cell
                dcx = 0
                dcy = 0
            if inarray[ir,jc] != nodata:             #included only masked area; The DEM should be no-data for those outside of mask (e.g. watershed)
                cellId = maskArray[ir,jc]          #ir*Nrows + jc
                dcellId = maskArray[ir+dcy, jc+dcx]   #(ir+dcy)*Nrows + jc+dcx
                if inarray[ir+dcy, jc+dcx] == nodata and dcellId >0:         # flowing out of region, make sure it is -ve value
                    dcellId = -1*dcellId
                if(dcx==0 and dcy==0) and dcellId >0:                      #outlet
                    dcellId = -1*dcellId
                xcoord = xmin + dx * jc     # + 0.5*dx            #following the conn. file convention, left border of cell
                ycoord = ymax + dy*(ir+1)   # - 0.5*dy            # bottom border
                conFile.write('     %d  %d      Rv   1      %d      %d      %f      %f      %f\n' %(cellId,dcellId,ir,jc, abs(dx*dy)*in_AreaunitConv, sterographic_toHrapX(xcoord), sterographic_toHrapY(ycoord)) )  #distances are in units of Hrap
    conFile.close()


#leftX=-112.0, topY=42.3, rightX=-111.0, bottomY=41.6,
def create_netCDF_from_multple_hdf4(file_prefix, output_netcdf, leftX, topY, rightX, bottomY,out_timeName = 'time'):
    """
    Subsets and combines multiple hdf4 raster files
    file_prefix is something like otgrid for  otgrid_20091001120000.hdf
    Sinusoidal projection:
     SR-ORG:6974
     '+proj=sinu +R=6371007.181 +nadgrids=@null +wktext'
    """
    #1 find the hxxvyy for the watershed in question
    #2 cp *.hxxvyy.* sub    copy only the tiles of interest to different directory(sub)
    #3 gdalinfo to find out the name of the variable of interest & the hdf type (sds, eos etc)
    #4 for i in *.hdf; do gdal_translate -of netCDF -co "FORMAT=NC4" HDF4_EOS:EOS_GRID:"$i":MOD_Grid_MOD15A2H:Lai_500m nc/$i.nc; done
    #5 for i in *.nc; do gdalwarp -s_srs '+proj=sinu +R=6371007.181 +nadgrids=@null +wktext' -t_srs EPSG:4326 -of netCDF -co "FORMAT=NC4" $i proj/$i; done
    #6 for i in *.nc; do gdal_translate -projwin -112.0 42.3 -111.0 41.6 -a_srs EPSG:4326 -of netCDF -co "FORMAT=NC4" $i projsub/$i; done
    #7 ncecat -u "time" *.nc logMODISLAI.nc          #concatenate along time dimension


    cmdString = "for i in *.asc; do  gdal_translate -projwin "+str(leftX)+" "+str(topY)+" "+str(rightX)+" "+str(bottomY)\
                +" -of netCDF -co \"FORMAT=NC4\" $i  $i.nc; done"
    callSubprocess(cmdString, 'subset asc files')

    cmdString = "ncecat -u "+out_timeName+" "+file_prefix+"*.nc "+ output_netcdf
    callSubprocess(cmdString, 'concatenate netCDF files')


#leftX=-112.0, topY=42.3, rightX=-111.0, bottomY=41.6,
def create_netCDF_from_multple_asc(file_prefix, output_netcdf, leftX, topY, rightX, bottomY,out_timeName = 'time'):
    """
    Subsets and combines multiple .asc raster files
    file_prefix is something like otgrid for  otgrid_20091001120000.asc
    """
    cmdString = "for i in *.asc; do  gdal_translate -projwin "+str(leftX)+" "+str(topY)+" "+str(rightX)+" "+str(bottomY)\
                +" -of netCDF -co \"FORMAT=NC4\" $i  $i.nc; done"
    callSubprocess(cmdString, 'subset asc files')

    cmdString = "ncecat -u "+out_timeName+" "+file_prefix+"*.nc "+ output_netcdf
    callSubprocess(cmdString, 'concatenate netCDF files')



def project_and_resample_Array(input_array, srs_geotrs, srs_proj, Nxin, Nyin, reference_netcdf):  #, output_array):

    #srs_data = gdal.Open(input_raster, GA_ReadOnly)
    #srs_proj = srs_data.GetProjection() #osr.SpatialReference(wkt
    srs_data = gdal.GetDriverByName('MEM').Create('', Nxin, Nyin, 1, gdal.GDT_Float32)
    srs_data.SetGeoTransform(srs_geotrs)
    srs_data.SetProjection(srs_proj)
    srsband = srs_data.GetRasterBand(1)
    srsband.WriteArray(input_array)
    srsband.FlushCache()

    ref_data = gdal.Open(reference_netcdf, GA_ReadOnly)
    ref_proj = ref_data.GetProjection()
    ref_geotrs = ref_data.GetGeoTransform()
    Ncols = ref_data.RasterXSize
    Nrows = ref_data.RasterYSize
    ref_data = None

    out_data = gdal.GetDriverByName('MEM').Create('', Ncols, Nrows, 1, gdal.GDT_Float32)
    out_data.SetGeoTransform(ref_geotrs)
    out_data.SetProjection(ref_proj)

    gdal.ReprojectImage(srs_data,out_data,srs_proj,ref_proj, gdal.GRA_Bilinear )
    output_array = out_data.ReadAsArray()

    srs_data = None
    out_data = None
    return output_array



def combineRasters(input_raster1, input_raster2, output_raster):
    """To  Do: may need to specify output no-data value
    """
    cmdString = "gdalwarp "+input_raster1+" "+input_raster2+" "+output_raster
    callSubprocess(cmdString, 'join (stitch) two raster files')


def callSubprocess(cmdString, debugString):
    cmdargs = shlex.split(cmdString)
    debFile = open('debug_file.txt', 'w')
    debFile.write('Starting %s \n' % debugString)
    retValue = subprocess.call(cmdargs,stdout=debFile)
    if (retValue==0):
        debFile.write('%s Successful\n' % debugString)
        debFile.close()
    else:
        debFile.write('There was error in %s\n' % debugString)
        debFile.close()


def rasterToNetCDF(input_raster, output_netcdf):
    cmdString = "gdal_translate -of netCDF "+input_raster+" "+output_netcdf
    callSubprocess(cmdString, 'raster to netcdf')

def uncompressRaster(input_raster, output_raster):
    """TauDEM doesn't take compressed file; uncompress file
        ToDO:  Check compression first"""
    cmdString = "gdal_translate -co COMPRESS=NONE"+" "\
               +input_raster+" "+output_raster
    callSubprocess(cmdString, 'uncompress raster')
