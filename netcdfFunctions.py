"""
This is used for both latest and old version of model processing
"""

try:
    from osgeo import gdal, osr, ogr
except:
    import gdal, osr, ogr
from gdalconst import *
import shlex
import subprocess
import os
import numpy
import netCDF4


def flip_netCDF_yaxis(input_netcdf, output_netcdf, input_varname='Band1'):        #, output_varname='Band1'):
    """
    """
    #temp_netcdf = "temp_"+output_netcdf
    cmdString = "nccopy -4 "+input_netcdf+" "+output_netcdf
    callSubprocess(cmdString, 'copy netcdf with dimensions')

    ncOut = netCDF4.Dataset(output_netcdf,"r+", format='NETCDF4')
    timeLen = len(ncOut.dimensions['time'])
    #yin = ncOut.variables['y'][:]
    #ncOut.variables['y'][:] = yin[::-1]

    for tk in range(timeLen):
        varin = ncOut.variables[input_varname][tk,:,:]
        ncOut.variables[input_varname][tk,:,:] = varin[::-1]

    ncOut.close()


def create_Tiff_from_text(input_textFile, output_raster, Nxin, Nyin, topleftX, topleftY, dx, dy, no_data_value=-9999, epsgCode = 4326, reverse_array = True):
    """
    :param input_textFile:
    :return:
    """
    inp1dArray = numpy.loadtxt(input_textFile)
    input_array = numpy.reshape(inp1dArray,(-1,Nxin))
    #print(input_array)


    srs_data = gdal.GetDriverByName('GTiff').Create(output_raster, Nxin, Nyin, 1, gdal.GDT_Float32)
    srs_data.SetGeoTransform((topleftX,dx,0,topleftY,0,dy))
    srs_proj = osr.SpatialReference()
    srs_proj.ImportFromEPSG(epsgCode)
    srs_data.SetProjection(srs_proj.ExportToWkt())
    srsband = srs_data.GetRasterBand(1)
    if reverse_array:
        output_array = input_array[::-1]    # array from text read starting with bottom-left (0,0)
    else:
        output_array = input_array
    srsband.WriteArray(output_array)
    srsband.SetNoDataValue(no_data_value)
    srsband.FlushCache()

    srs_data = None



def create_3Dnetcdf_from_text(input_textFile, output_netcdf, varName, in_Time = 'time', in_Xcoord = 'x', in_Ycoord='y', vardataType1='float32',fillValue=-9999):
    """

    """
    # copy
    inp2dArray = numpy.loadtxt(input_textFile)
    #input_array = numpy.reshape(inp1dArray,(-1,Nxin))
    print(inp2dArray)

    #re-open file to write re-gridded data
    ncOut = netCDF4.Dataset(output_netcdf,"r+") #, format='NETCDF4')
    timeLen = len(ncOut.dimensions[in_Time])
    xin = ncOut.variables[in_Xcoord][:]
    yin = ncOut.variables[in_Ycoord][:]
    vardataType = numpy.dtype(vardataType1)
    varin = numpy.zeros((timeLen),dtype=vardataType)

    for yi in range(len(yin)):
        for xj in range(len(xin)):
            ncOut.variables[varName][yi,xj,:] = inp2dArray[:,xj]
    ncOut.close()
    #delete temp netcdf file
    #os.remove(temp_netcdf)


def set_2Dnetcdf_values_from_text(input_textFile, output_netcdf, varName, in_Time = 'time', in_Xcoord = 'x', in_Ycoord='y', vardataType1='float32',fillValue=-9999):
    """

    """
    # copy
    in1dArray = numpy.loadtxt(input_textFile)
    #input_array = numpy.reshape(inp1dArray,(-1,Nxin))
    print(in1dArray)

    #re-open file to write re-gridded data
    ncOut = netCDF4.Dataset(output_netcdf,"r+") #, format='NETCDF4')
    #timeLen = len(ncOut.dimensions[in_Time])
    xin = ncOut.variables[in_Xcoord][:]
    yin = ncOut.variables[in_Ycoord][:]
    #vardataType = numpy.dtype(vardataType1)
    #varin = numpy.zeros((timeLen),dtype=vardataType)

    for yi in range(len(yin)):
        for xj in range(len(xin)):
            ncOut.variables[varName][yi,xj] = in1dArray[xj]
    ncOut.close()
    #delete temp netcdf file
    #os.remove(temp_netcdf)


#set specfic cells to a value
def set_specific_cells_Values(output_netcdf, varName, in_Xcoord = 'x', in_Ycoord='y', setValue = 1, dataSize = 3, ycoord=[12,22,27], xcoord=[25,30,17], fillValue=-32768):
    """

    """
    #re-open file to write re-gridded data
    ncOut = netCDF4.Dataset(output_netcdf,"r+") #, format='NETCDF4')
    xin = ncOut.variables[in_Xcoord][:]
    yin = ncOut.variables[in_Ycoord][:]
    dataType = ncOut.variables[varName].datatype               # data type for time variable same as x variable
    varout = numpy.zeros((len(yin),len(xin)),dtype=dataType)
    varout[varout == 0] = fillValue
    for cindx in range(dataSize):
        varout[ycoord[cindx],xcoord[cindx]] = setValue
    ncOut.variables[varName][:,:] = varout[:,:]
    ncOut.close()


# new value = multiplier_factor * old value  +  offset
# C = 0.5556F + (- 17.7778)
# F = 1.8C + 32
def convert_netcdf_Units(input_netcdf, output_netcdf, varName, new_varUnits=" ", multiplier_Factor = 1.0, offset = 0.0):
    """
    does unit conversion for a variable in netcdf file
    :param input_netcdf: input a
    :param output_netcdf: output
    :param varName: name of variable of interest
    :param new_varUnits: name of the new unit after conversion
    :param multiplier_Factor: self explanatory
    :param offset: additive factor
    :return:
    """
    temp_netcdf = "temp_"+output_netcdf
    cmdString = "nccopy -4 "+input_netcdf+" "+temp_netcdf             #output_netcdf
    callSubprocess(cmdString, 'copy netcdf with dimensions')

    ncIn = netCDF4.Dataset(input_netcdf,"r") # format='NETCDF4')
    timeLen = len(ncIn.dimensions['time'])

    ncOut = netCDF4.Dataset(temp_netcdf,"r+", format='NETCDF4')
    for tk in range(timeLen):
        varin = ncIn.variables[varName][tk,:,:]
        varOut = offset + multiplier_Factor * varin
        ncOut.variables[varName][tk,:,:] = varOut[:,:]
    ncIn.close()
    ncOut.close()
    #cmdString = "ncap2 -s \' "+varName+" = "+str(offset)+" + "+str(multiplier_Factor)+"f * "+varName+" \' -O "+input_netcdf+" "+temp_netcdf
    #callSubprocess(cmdString, 'convert netcdf units')
    cmdString = "ncatted -a units,"+varName+",m,c,\'"+new_varUnits+"\' -O "+temp_netcdf+" "+output_netcdf
    callSubprocess(cmdString, 'rename netcdf units')


###TODO: 9.6.17 This is not enough need to ensure watershed.nc no-value cells don't get temp value (curr. getting 0 by interp)
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



def project_subset_and_resample_netCDF_to_referenceNetCDF(input_netcdf, reference_netcdf, varName_ref, varName, output_netcdf):
    """This re-grids a netcdf to target/reference resolution
    Input coordinates are time, y, x
    Warning: Works only if the target boundary is within the input boundary & the coordinates directions are
    the same, i.e. y increasing / decreasing """

    # Read input geo information
    #srs_data = gdal.Open(input_netcdf, GA_ReadOnly)
    srs_data = gdal.Open('NetCDF:"'+input_netcdf+'":'+varName)
    srs_geotrs = srs_data.GetGeoTransform()
    Nxin = srs_data.RasterXSize
    Nyin = srs_data.RasterYSize
    srs_proj = srs_data.GetProjection()
    print(srs_proj)
    #srs_projt = srs_proj.ExportToWkt()
    srs_data = None

    #Add dummy dimensions and variables
    temp_netcdf = "temp"+output_netcdf
    cmdString = "ncrename -O -d x,x_2 -d y,y_2 -v x,x_2 -v y,y_2 -v"+varName+","+varName+"_2 "+\
                 input_netcdf+" "+temp_netcdf
    callSubprocess(cmdString, 'copy netcdf with rename old dimensions')

    ncRef = netCDF4.Dataset(reference_netcdf,"r") # format='NETCDF4')
    xout = ncRef.variables['x'][:]
    yout = ncRef.variables['y'][:]

    ref_grid_mapping = getattr(ncRef.variables[varName_ref],'grid_mapping')
    varAtts_ref = ncRef.variables[ref_grid_mapping].ncattrs()
    attDict_ref = dict.fromkeys(varAtts_ref)
    for attName in varAtts_ref:
        attDict_ref[attName] = getattr(ncRef.variables[ref_grid_mapping],attName)

    ncRef.close()

    ncIn = netCDF4.Dataset(input_netcdf,"r") # format='NETCDF4')
    xin = ncIn.variables['x'][:]
    yin = ncIn.variables['y'][:]

    ncOut = netCDF4.Dataset(temp_netcdf,"r+", format='NETCDF4')
    ncOut.createDimension('y',len(yout))
    ncOut.createDimension('x', len(xout))
    dataType = ncIn.variables['x'].datatype
    vardataType = ncIn.variables[varName].datatype
    ncOut.createVariable('y',dataType,('y',))
    ncOut.createVariable('x',dataType,('x',))
    ncOut.variables['y'][:] = yout[:]
    ncOut.variables['x'][:] = xout[:]
    ncOut.createVariable(varName,vardataType,('time','y','x',))

    #grid mapping
    ncOut.createVariable(ref_grid_mapping,'c',())
    ncOut.variables[ref_grid_mapping].setncatts(attDict_ref)

    #Copy attributes
    varAtts = ncIn.variables[varName].ncattrs()
    attDict = dict.fromkeys(varAtts)
    for attName in varAtts:
        attDict[attName] = getattr(ncIn.variables[varName],attName)
        if attName == 'grid_mapping':
            attDict[attName] = ref_grid_mapping
    ncOut.variables[varName].setncatts(attDict)

    xAtts = ncIn.variables['x'].ncattrs()
    attDict = dict.fromkeys(xAtts)
    for attName in xAtts:
        attDict[attName] = getattr(ncIn.variables['x'],attName)
    ncOut.variables['x'].setncatts(attDict)
    yAtts = ncIn.variables['y'].ncattrs()
    attDict = dict.fromkeys(yAtts)
    for attName in yAtts:
        attDict[attName] = getattr(ncIn.variables['y'],attName)
    ncOut.variables['y'].setncatts(attDict)
    ncOut.close()
    #delete the old variables
    cmdString = "ncks -4 -O -C -x -v x_2,y_2,"+varName+"_2 "+\
                 temp_netcdf+" "+output_netcdf
    callSubprocess(cmdString, 'delete old dimensions')

    #re-open file to write re-gridded data
    ncOut = netCDF4.Dataset(output_netcdf,"r+") #, format='NETCDF4')
    varin = numpy.zeros((len(yin),len(xin)),dtype=vardataType)
    varout = numpy.zeros((len(yout),len(xout)),dtype=vardataType)
    timeLen = len(ncIn.dimensions['time'])
    print(timeLen)
    for tk in range(3): #(timeLen):
        varin[:,:] = ncIn.variables[varName][tk,:,:]
        print(varin)
        #Because gdal tif and Daymet nc y axes directions differ, here array is reversed
        #7.11.15 no need for below
        #varin_rev = varin[::-1]
        varout[:,:] = project_and_resample_Array(varin, srs_geotrs, srs_proj, Nxin, Nyin, reference_netcdf)
        print(varout)
        varout_comp = project_and_resample_Array(varin, srs_geotrs, srs_proj, Nxin, Nyin, reference_netcdf)
        print(varout_comp)
        ncOut.variables[varName][tk,:,:] = varout[:,:]
    ncIn.close()
    ncOut.close()
    #delete temp netcdf file
    #os.remove(temp_netcdf)



def project_subset_and_resample_netCDF_to_referenceNetCDF2(input_netcdf, varName, reference_netcdf, varName_ref, output_netcdf,
            in_epsgCode=None, tSampling_interval=1, start_Time = 0.0,dT = 3.0, time_unitString ='hours since 2008-10-01 00:00:00 UTC',
                                                               in_Time = 'time', in_Xcoord = 'lon', in_Ycoord='lat'):
    """This re-grids a netcdf to target/reference resolution
    Input coordinates are time, y, x
    Warning: Works only if the target boundary is within the input boundary & the coordinates directions are
    the same, i.e. y increasing / decreasing """
    #epsg=4326
    # Read input geo information
    #srs_data = gdal.Open(input_netcdf, GA_ReadOnly)
    srs_data = gdal.Open('NetCDF:"'+input_netcdf+'":'+varName)
    srs_geotrs = srs_data.GetGeoTransform()
    Nxin = srs_data.RasterXSize
    Nyin = srs_data.RasterYSize

    if in_epsgCode == None:
        srs_proj = srs_data.GetProjection()
    else:
        srs_proj = osr.SpatialReference()
        srs_proj.ImportFromEPSG(in_epsgCode)

    srs_projt = srs_proj.ExportToWkt()
    srs_data = None

    #Add dummy dimensions and variables
    temp_netcdf = "temp"+output_netcdf
    cmdString = "nccopy -4 "+reference_netcdf+" "+temp_netcdf             #output_netcdf
    callSubprocess(cmdString, 'copy netcdf with dimensions')

    ncIn = netCDF4.Dataset(input_netcdf,"r") # format='NETCDF4')
    xin = ncIn.variables[in_Xcoord][:]
    yin = ncIn.variables[in_Ycoord][:]
    timeLen = len(ncIn.dimensions[in_Time])
    dataType = ncIn.variables[in_Xcoord].datatype               # data type for time variable same as x variable
    vardataType = ncIn.variables[varName].datatype
    tin = numpy.zeros(int(timeLen/tSampling_interval),dtype=dataType)
    for tk in range(int(timeLen/tSampling_interval)):
        tin[tk] = start_Time + tk*dT*tSampling_interval

    ncOut = netCDF4.Dataset(temp_netcdf,"r+", format='NETCDF4')
    xout = ncOut.variables['x'][:]
    yout = ncOut.variables['y'][:]
    ref_grid_mapping = getattr(ncOut.variables[varName_ref],'grid_mapping')
    ncOut.createDimension(in_Time,int(timeLen/tSampling_interval))
    ncOut.createVariable(in_Time,dataType,(in_Time,))
    ncOut.variables[in_Time][:] = tin[:]
    ncOut.createVariable(varName,vardataType,(in_Time,'y','x',))

    #Copy attributes
    varAtts = ncIn.variables[varName].ncattrs()
    attDict = dict.fromkeys(varAtts)
    grid_map_set = False
    for attName in varAtts:
        attDict[attName] = getattr(ncIn.variables[varName],attName)
        if attName == 'grid_mapping':
            attDict[attName] = ref_grid_mapping
            grid_map_set = True
    if grid_map_set == False:                        #no attribute grid-map in input variable attributes
        attDict['grid_mapping'] = ref_grid_mapping
    ncOut.variables[varName].setncatts(attDict)

    """tAtts = ncIn.variables[in_Time].ncattrs()
    attDict = dict.fromkeys(tAtts)
    for attName in tAtts:
        attDict[attName] = getattr(ncIn.variables[in_Time],attName)
    """
    attDict = {'calendar':'standard', 'long_name':'time'}
    attDict['units'] = time_unitString
    ncOut.variables[in_Time].setncatts(attDict)
    ncOut.close()
    #delete old variables
    cmdString = "ncks -4 -C -h -O -x -v "+varName_ref+" "+temp_netcdf+" "+output_netcdf
    callSubprocess(cmdString, 'delete old/reference variable')

    #re-open file to write re-gridded data
    ncOut = netCDF4.Dataset(output_netcdf,"r+") #, format='NETCDF4')
    varin = numpy.zeros((len(yin),len(xin)),dtype=vardataType)
    varout = numpy.zeros((len(yout),len(xout)),dtype=vardataType)

    for tk in range(int(timeLen/tSampling_interval) ):
        varout[:,:] = 0.0
        for tkin in range(tSampling_interval):
            #print("loop: "+str(tk)+" varin = ")
            varin[:,:] = ncIn.variables[varName][(tk*tSampling_interval+tkin),:,:]
            #print(varin)
            #Because gdal and netCDF4 (and NCO) read the data array in (y-) reverse order, need to adjust orientation
            varin_rev = varin[::-1]
            varoutret = project_and_resample_Array(varin_rev, srs_geotrs, srs_projt, Nxin, Nyin, reference_netcdf)
            #print("out resample = ")
            #print(varoutret)
            varout[:,:] = varout[:,:] + varoutret[:,:]
            #print("var out to write")
        varout[:,:] = varout[:,:] / tSampling_interval
        #print(varout)
        ncOut.variables[varName][tk,:,:] = varout[::-1]         # reverse back the array for netCDF4
    ncIn.close()
    ncOut.close()
    #delete temp netcdf file
    #os.remove(temp_netcdf)



def project_and_resample_Array_Byte(input_array, srs_geotrs, srs_proj, Nxin, Nyin, reference_netcdf):  #, output_array):


    #srs_data = gdal.Open(input_raster, GA_ReadOnly)
    #srs_proj = srs_data.GetProjection() #osr.SpatialReference(wkt
    srs_data = gdal.GetDriverByName('MEM').Create('', Nxin, Nyin, 1, gdal.GDT_Unknown)
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

    out_data = gdal.GetDriverByName('MEM').Create('', Ncols, Nrows, 1, gdal.GDT_Byte)
    out_data.SetGeoTransform(ref_geotrs)
    out_data.SetProjection(ref_proj)

    gdal.ReprojectImage(srs_data,out_data,srs_proj,ref_proj, gdal.GRA_NearestNeighbour )
    output_array = out_data.ReadAsArray()

    srs_data = None
    out_data = None
    return output_array


#replaces no data values with the average of each row
def replace_no_data_with_SpatialAverage(input_netcdf, varName, output_netcdf, in_Time = 'time', in_Xcoord = 'lon_110', in_Ycoord='lat_110', vardataType1='float32',fillValue=-9999):
    """

    """
    # copy
    cmdString = "nccopy -k 4 "+input_netcdf+" "+output_netcdf             #output_netcdf
    callSubprocess(cmdString, 'copy netcdf with dimensions')

    #re-open file to write re-gridded data
    ncOut = netCDF4.Dataset(output_netcdf,"r+") #, format='NETCDF4')
    timeLen = len(ncOut.dimensions[in_Time])
    xin = ncOut.variables[in_Xcoord][:]
    yin = ncOut.variables[in_Ycoord][:]
    vardataType = numpy.dtype(vardataType1)
    varin = numpy.zeros((len(yin),len(xin)),dtype=vardataType)
    #varout = numpy.zeros((len(yout),len(xout)),dtype=vardataType)
    for tk in range(timeLen):
        varin = ncOut.variables[varName][tk,:,:]
        #Because gdal and netCDF4 (and NCO) read the data array in (y-) reverse order, need to adjust orientation
        #varin_rev = (varin[::-1])
        varin[varin == fillValue] = numpy.nan
        rowmean = numpy.nanmean(varin,axis=1)
        nanIndices = numpy.where(numpy.isnan(varin))
        varin[nanIndices] = numpy.take(rowmean,nanIndices[0])
        ncOut.variables[varName][tk,:,:] = varin[:,:]
    ncOut.close()
    #delete temp netcdf file
    #os.remove(temp_netcdf)


#fill in missing values with linearly interploated values in the time dimension
def replace_no_data_along_temporal_Dimension(input_netcdf, varName, output_netcdf, in_Time = 'time', in_Xcoord = 'lon', in_Ycoord='lat', vardataType1='float32',fillValue=-9999):
    """
    fill in missing values with linearly interpolated values in the time dimension
    """
    # copy
    cmdString = "nccopy -k 4 "+input_netcdf+" "+output_netcdf             #output_netcdf
    callSubprocess(cmdString, 'copy netcdf with dimensions')

    #re-open file to write re-gridded data
    ncOut = netCDF4.Dataset(output_netcdf,"r+") #, format='NETCDF4')
    timeLen = len(ncOut.dimensions[in_Time])
    xin = ncOut.variables[in_Xcoord][:]
    yin = ncOut.variables[in_Ycoord][:]
    tin = ncOut.variables[in_Time][:]
    vardataType = numpy.dtype(vardataType1)
    varin = numpy.zeros(timeLen,dtype=vardataType)
    #varout = numpy.zeros((len(yout),len(xout)),dtype=vardataType)
    for yj in range(len(yin)):
        for xi in range(len(xin)):
            varin = ncOut.variables[varName][:,yj,xi]
            #Because gdal and netCDF4 (and NCO) read the data array in (y-) reverse order, need to adjust orientation
            #varin_rev = (varin[::-1])
            varin[varin == fillValue] = numpy.nan
            nans = numpy.isnan(varin)
            varin[nans] = numpy.interp(tin[nans],tin[~nans],varin[~nans])
            ncOut.variables[varName][:,yj, xi] = varin[:]
    ncOut.close()
    #delete temp netcdf file
    #os.remove(temp_netcdf)


#find fSCA, grid cells with >=Wthreshold for
def find_values_atTimestep(input_netcdf, varName, output_txtFile, timeStep,
                                   vardataType1='float32', in_Time = 'time', in_Xcoord = 'x', in_Ycoord='y'):
    """
    find the number of grid cells above threshold
    :param input_netcdf:
    :param varName:
    :param output_txtFile:
    :param valThreshold:
    :return:
    """

    ncIn = netCDF4.Dataset(input_netcdf,"r") # format='NETCDF4')
    xin = ncIn.variables[in_Xcoord][:]
    yin = ncIn.variables[in_Ycoord][:]
    timeLen = len(ncIn.dimensions[in_Time])
    dataType = ncIn.variables[in_Xcoord].datatype               # data type for time variable same as x variable
    vardataType = numpy.dtype(vardataType1)

    varin = numpy.zeros((len(yin),len(xin)),dtype=vardataType)
    varinW = numpy.zeros((timeLen),dtype=vardataType)

    #ncInW = netCDF4.Dataset(Winput_netcdf,"r") # format='NETCDF4')
    #varinW[:] = ncInW.variables[WvarName][0,:]

    afFile = open(output_txtFile, 'w')
    afFile.write('cell#  Value \n')
    #afInt = 0

    varin[:,:] = ncIn.variables[varName][timeStep,:,:]
    for yj in range(len(yin)):
        for xi in range(len(xin)):
            afFile.write('%d  %f\n' %((yj*len(xin) + xi), varin[yj,xi]))
    #
    ncIn.close()
    afFile.close()



#find fSCA, grid cells with >=Wthreshold for
def find_gridcells_above_Threshold(input_netcdf, varName, Winput_netcdf, WvarName, output_txtFile, valThreshold,
                                   vardataType1='float32', in_Time = 'time', in_Xcoord = 'x', in_Ycoord='y'):
    """
    find the number of grid cells above threshold
    :param input_netcdf:
    :param varName:
    :param output_txtFile:
    :param valThreshold:
    :return:
    """

    ncIn = netCDF4.Dataset(input_netcdf,"r") # format='NETCDF4')
    xin = ncIn.variables[in_Xcoord][:]
    yin = ncIn.variables[in_Ycoord][:]
    timeLen = len(ncIn.dimensions[in_Time])
    dataType = ncIn.variables[in_Xcoord].datatype               # data type for time variable same as x variable
    vardataType = numpy.dtype(vardataType1)

    varin = numpy.zeros((len(yin),len(xin)),dtype=vardataType)
    varinW = numpy.zeros((timeLen),dtype=vardataType)

    ncInW = netCDF4.Dataset(Winput_netcdf,"r") # format='NETCDF4')
    varinW[:] = ncInW.variables[WvarName][0,:]

    afFile = open(output_txtFile, 'w')
    afFile.write('Timestep   #cells>Threshold  Af  Wa \n')
    afInt = 0
    totCells = len(yin) * len(xin)
    for tk in range(timeLen):
        afInt = 0
        varin[:,:] = ncIn.variables[varName][tk,:,:]
        for yj in range(len(yin)):
            for xi in range(len(xin)):
                if(varin[yj,xi] > valThreshold):
                    afInt = afInt + 1
        #
        afFile.write('%d   %d   %f  %f\n' %(tk, afInt, afInt/totCells,  varinW[tk]))
    #
    ncIn.close()
    afFile.close()


#find  grid cells with >= threshold for NDSI
def find_gridcells_above_Threshold2(input_netcdf, varName, output_txtFile, valThreshold,
                                   vardataType1='float32', in_Time = 'time', in_Xcoord = 'x', in_Ycoord='y'):
    """
    find the number of grid cells above threshold
    :param input_netcdf:
    :param varName:
    :param output_txtFile:
    :param valThreshold:
    :return:
    """

    ncIn = netCDF4.Dataset(input_netcdf,"r") # format='NETCDF4')
    xin = ncIn.variables[in_Xcoord][:]
    yin = ncIn.variables[in_Ycoord][:]
    timeLen = len(ncIn.dimensions[in_Time])
    #dataType = ncIn.variables[in_Xcoord].datatype               # data type for time variable same as x variable
    vardataType = numpy.dtype(vardataType1)

    varin = numpy.zeros((len(yin),len(xin)),dtype=vardataType)

    afFile = open(output_txtFile, 'w')
    afFile.write('Timestep   #cells>Threshold  fraction_totCells \n')
    afInt = 0
    totCells = len(yin) * len(xin)
    for tk in range(timeLen):
        afInt = 0
        varin[:,:] = ncIn.variables[varName][tk,:,:]
        for yj in range(len(yin)):
            for xi in range(len(xin)):
                if(varin[yj,xi] > valThreshold):
                    afInt = afInt + 1
        #
        afFile.write('%d  %d  %f\n' %(tk, afInt, afInt/totCells))
    #
    ncIn.close()
    afFile.close()



#for getting normalized indices
def get_normalized_Index(input_netcdf1, varName1, input_netcdf2, varName2, output_netcdf, varNameout,
        vardataType1='float32',valid_min=0, valid_max=100, multiplier_Factor = 0.001, offset = 0, vUnits='-', in_Time = 'time', in_Xcoord = 'x', in_Ycoord='y'):
    """

    """
    #epsg=4326
    # Read input geo information
    #srs_data = gdal.Open(input_netcdf, GA_ReadOnly)


    #Add dummy dimensions and variables
    temp_netcdf = "temp"+output_netcdf
    cmdString = "nccopy -k 4 "+input_netcdf1+" "+temp_netcdf             #output_netcdf
    callSubprocess(cmdString, 'copy netcdf with dimensions')

    ncIn1 = netCDF4.Dataset(input_netcdf1,"r") # format='NETCDF4')
    ncIn2 = netCDF4.Dataset(input_netcdf2,"r") # format='NETCDF4')
    xin = ncIn1.variables[in_Xcoord][:]
    yin = ncIn1.variables[in_Ycoord][:]
    timeLen = len(ncIn1.dimensions[in_Time])

    ncOut = netCDF4.Dataset(temp_netcdf,"r+", format='NETCDF4')
    xout = ncOut.variables[in_Xcoord][:]
    yout = ncOut.variables[in_Ycoord][:]
    ref_grid_mapping = getattr(ncOut.variables[varName1],'grid_mapping')
    vardataType = numpy.dtype(vardataType1)
    ncOut.createVariable(varNameout,vardataType,(in_Time,in_Ycoord,in_Xcoord,))


    attDict = {'name':varNameout, 'long_name':varNameout}
    attDict['units'] = vUnits
    attDict['grid_mapping'] = ref_grid_mapping
    #attDict['_FillValue'] = fillValue
    ncOut.variables[varNameout].setncatts(attDict)
    ncOut.close()
    #delete old variables
    cmdString = "ncks -4 -C -h -O -x -v "+varName1+" "+temp_netcdf+" "+output_netcdf
    callSubprocess(cmdString, 'delete old/reference variable')

    #re-open file to write re-gridded data
    ncOut = netCDF4.Dataset(output_netcdf,"r+") #, format='NETCDF4')
    varin1 = numpy.zeros((len(yin),len(xin)),dtype=vardataType)
    varin2 = numpy.zeros((len(yin),len(xin)),dtype=vardataType)
    varout = numpy.zeros((len(yout),len(xout)),dtype=vardataType)
    for tk in range(timeLen):
        varin1[:,:] = ncIn1.variables[varName1][tk,:,:]
        varin2[:,:] = ncIn2.variables[varName2][tk,:,:]
        varout = (varin1 - varin2) / (varin1 + varin2)
        ncOut.variables[varNameout][tk,:,:] = varout[::]         # reverse back the array for netCDF4
    ncIn1.close()
    ncIn2.close()
    ncOut.close()



#for MODIS fSCA
def project_subset_and_resample_netCDF_to_referenceNetCDF4(input_netcdf, varName, reference_netcdf, varName_ref, output_netcdf,
        vardataType1='float32',fillValue=-9999, valid_min=0, valid_max=100, in_epsgCode=4326, tSampling_interval=1, tInterpolation_frequency=8, start_Time = 0.0,dT = 1.0, multiplier_Factor = 0.001, offset = 0,
            Vunits='-', time_unitString ='hours since 2014-01-01 00:00:00 UTC', in_Time = 'time', in_Xcoord = 'lon', in_Ycoord='lat'):
    """This re-grids a netcdf to target/reference resolution
    Input coordinates are time, y, x
    Warning: Works only if the target boundary is within the input boundary & the coordinates directions are
    the same, i.e. y increasing / decreasing """
    #epsg=4326
    # Read input geo information
    #srs_data = gdal.Open(input_netcdf, GA_ReadOnly)
    srs_data = gdal.Open('NetCDF:"'+input_netcdf+'":'+varName)
    srs_geotrs = srs_data.GetGeoTransform()
    Nxin = srs_data.RasterXSize
    Nyin = srs_data.RasterYSize
    if in_epsgCode == None:
        srs_proj = srs_data.GetProjection()
        srs_projt = srs_proj.ExportToWkt()
    else:
        srs_proj = osr.SpatialReference()
        srs_proj.ImportFromEPSG(in_epsgCode)
        srs_projt = srs_proj.ExportToWkt()
    srs_data = None

    #Add dummy dimensions and variables
    temp_netcdf = "temp"+output_netcdf
    cmdString = "nccopy -k 4 "+reference_netcdf+" "+temp_netcdf             #output_netcdf
    callSubprocess(cmdString, 'copy netcdf with dimensions')

    ncIn = netCDF4.Dataset(input_netcdf,"r") # format='NETCDF4')
    xin = ncIn.variables[in_Xcoord][:]
    yin = ncIn.variables[in_Ycoord][:]
    timeLen = len(ncIn.dimensions[in_Time])
    dataType = ncIn.variables[in_Xcoord].datatype               # data type for time variable same as x variable
    #vardataType = ncIn.variables[varName].datatype
    tin = numpy.zeros(((timeLen*tInterpolation_frequency)/tSampling_interval),dtype=dataType)
    for tk in range(int((timeLen*tInterpolation_frequency)/tSampling_interval)):
        tin[tk] = start_Time + tk*dT*tSampling_interval

    ncOut = netCDF4.Dataset(temp_netcdf,"r+", format='NETCDF4')
    xout = ncOut.variables['x'][:]
    yout = ncOut.variables['y'][:]
    ref_grid_mapping = getattr(ncOut.variables[varName_ref],'grid_mapping')
    ncOut.createDimension(in_Time,(timeLen*tInterpolation_frequency)/tSampling_interval)
    ncOut.createVariable(in_Time,dataType,(in_Time,))
    ncOut.variables[in_Time][:] = tin[:]
    vardataType = numpy.dtype(vardataType1)
    ncOut.createVariable(varName,vardataType,(in_Time,'y','x',))

    #Copy attributes
    """
    varAtts = ncIn.variables[varName].ncattrs()
    attDict = dict.fromkeys(varAtts)
    grid_map_set = False
    fill_value_set = False
    for attName in varAtts:
        attDict[attName] = getattr(ncIn.variables[varName],attName)
        if attName == 'grid_mapping':
            attDict[attName] = ref_grid_mapping
            grid_map_set = True
        if attName == '_FillValue':
            attDict[attName] = float(fill_Value)
            fill_value_set = True
    if grid_map_set == False:                        #no attribute grid-map in input variable attributes
        attDict['grid_mapping'] = ref_grid_mapping
    if fill_value_set == False:
        attDict['_FillValue'] = float(fill_Value)
    """
    attDict = {'name':varName, 'long_name':varName}
    attDict['units'] = Vunits
    attDict['grid_mapping'] = ref_grid_mapping
    #attDict['_FillValue'] = fillValue
    ncOut.variables[varName].setncatts(attDict)

    """tAtts = ncIn.variables[in_Time].ncattrs()
    attDict = dict.fromkeys(tAtts)
    for attName in tAtts:
        attDict[attName] = getattr(ncIn.variables[in_Time],attName)
    """
    attDict = {'calendar':'standard', 'long_name':'time'}
    attDict['units'] = time_unitString
    ncOut.variables[in_Time].setncatts(attDict)
    ncOut.close()
    #delete old variables
    cmdString = "ncks -4 -C -h -O -x -v "+varName_ref+" "+temp_netcdf+" "+output_netcdf
    callSubprocess(cmdString, 'delete old/reference variable')

    #re-open file to write re-gridded data
    ncOut = netCDF4.Dataset(output_netcdf,"r+") #, format='NETCDF4')
    varin = numpy.zeros((len(yin),len(xin)),dtype=vardataType)
    varout = numpy.zeros((len(yout),len(xout)),dtype=vardataType)
    for tk in range(timeLen):
        varin[:,:] = ncIn.variables[varName][tk,:,:]
        #Because gdal and netCDF4 (and NCO) read the data array in (y-) reverse order, need to adjust orientation
        varin_rev = (varin[::-1])
        for yj in range(len(yin)):
            for xi in range(len(xin)):
                if((varin_rev[yj,xi] < valid_min) or (varin_rev[yj,xi] > valid_max)):
                    varin_rev[yj,xi] = fillValue
        varout[:,:] = offset + multiplier_Factor*(project_and_resample_Array_Byte(varin_rev, srs_geotrs, srs_projt, Nxin, Nyin, reference_netcdf))
        for tf in range (tInterpolation_frequency):
            ncOut.variables[varName][tk*tInterpolation_frequency+tf,:,:] = varout[::-1]         # reverse back the array for netCDF4
    ncIn.close()
    ncOut.close()
    #delete temp netcdf file
    #os.remove(temp_netcdf)


#this code written for time varying LAI from MODIS
def project_subset_and_resample_netCDF_to_referenceNetCDF3(input_netcdf, varName, reference_netcdf, varName_ref, output_netcdf,
        vardataType1='float32',fillValue=-9999.9, valid_min=0, valid_max=100, in_epsgCode=None, tSampling_interval=4, tInterpolation_frequency=8, start_Time = 0.0,dT = 1.0, multiplier_Factor = 1, offset = 0,
            Vunits='m^2/m^2', time_unitString ='hours since 2010-01-01 00:00:00 UTC', in_Time = 'time', in_Xcoord = 'lon_110', in_Ycoord='lat_110'):
    """This re-grids a netcdf to target/reference resolution
    Input coordinates are time, y, x
    Warning: Works only if the target boundary is within the input boundary & the coordinates directions are
    the same, i.e. y increasing / decreasing """
    #epsg=4326
    # Read input geo information
    #srs_data = gdal.Open(input_netcdf, GA_ReadOnly)
    srs_data = gdal.Open('NetCDF:"'+input_netcdf+'":'+varName)
    srs_geotrs = srs_data.GetGeoTransform()
    Nxin = srs_data.RasterXSize
    Nyin = srs_data.RasterYSize
    if in_epsgCode == None:
        srs_proj = srs_data.GetProjection()
        srs_projt = srs_proj.ExportToWkt()
    else:
        srs_proj = osr.SpatialReference()
        srs_proj.ImportFromEPSG(in_epsgCode)
        srs_projt = srs_proj.ExportToWkt()
    srs_data = None

    #Add dummy dimensions and variables
    temp_netcdf = "temp"+output_netcdf
    cmdString = "nccopy -k 4 "+reference_netcdf+" "+temp_netcdf             #output_netcdf
    callSubprocess(cmdString, 'copy netcdf with dimensions')

    ncIn = netCDF4.Dataset(input_netcdf,"r") # format='NETCDF4')
    xin = ncIn.variables[in_Xcoord][:]
    yin = ncIn.variables[in_Ycoord][:]
    timeLen = len(ncIn.dimensions[in_Time])
    dataType = ncIn.variables[in_Xcoord].datatype               # data type for time variable same as x variable
    #vardataType = ncIn.variables[varName].datatype
    tin = numpy.zeros(((timeLen*tInterpolation_frequency)/tSampling_interval),dtype=dataType)
    for tk in range(int((timeLen*tInterpolation_frequency)/tSampling_interval)):
        tin[tk] = start_Time + tk*dT*tSampling_interval

    ncOut = netCDF4.Dataset(temp_netcdf,"r+", format='NETCDF4')
    xout = ncOut.variables['x'][:]
    yout = ncOut.variables['y'][:]
    ref_grid_mapping = getattr(ncOut.variables[varName_ref],'grid_mapping')
    ncOut.createDimension(in_Time,(timeLen*tInterpolation_frequency)/tSampling_interval)
    ncOut.createVariable(in_Time,dataType,(in_Time,))
    ncOut.variables[in_Time][:] = tin[:]
    vardataType = numpy.dtype(vardataType1)
    ncOut.createVariable(varName,vardataType,(in_Time,'y','x',))

    #Copy attributes
    """
    varAtts = ncIn.variables[varName].ncattrs()
    attDict = dict.fromkeys(varAtts)
    grid_map_set = False
    fill_value_set = False
    for attName in varAtts:
        attDict[attName] = getattr(ncIn.variables[varName],attName)
        if attName == 'grid_mapping':
            attDict[attName] = ref_grid_mapping
            grid_map_set = True
        if attName == '_FillValue':
            attDict[attName] = float(fill_Value)
            fill_value_set = True
    if grid_map_set == False:                        #no attribute grid-map in input variable attributes
        attDict['grid_mapping'] = ref_grid_mapping
    if fill_value_set == False:
        attDict['_FillValue'] = float(fill_Value)
    """
    attDict = {'name':varName, 'long_name':varName}
    attDict['units'] = Vunits
    attDict['grid_mapping'] = ref_grid_mapping
    #attDict['_FillValue'] = fillValue
    ncOut.variables[varName].setncatts(attDict)

    """tAtts = ncIn.variables[in_Time].ncattrs()
    attDict = dict.fromkeys(tAtts)
    for attName in tAtts:
        attDict[attName] = getattr(ncIn.variables[in_Time],attName)
    """
    attDict = {'calendar':'standard', 'long_name':'time'}
    attDict['units'] = time_unitString
    ncOut.variables[in_Time].setncatts(attDict)
    ncOut.close()
    #delete old variables
    cmdString = "ncks -4 -C -h -O -x -v "+varName_ref+" "+temp_netcdf+" "+output_netcdf
    callSubprocess(cmdString, 'delete old/reference variable')

    #re-open file to write re-gridded data
    ncOut = netCDF4.Dataset(output_netcdf,"r+") #, format='NETCDF4')
    varin = numpy.zeros((len(yin),len(xin)),dtype=vardataType)
    varout = numpy.zeros((len(yout),len(xout)),dtype=vardataType)
    for tk in range(timeLen):
        varin[:,:] = ncIn.variables[varName][tk,:,:]
        #Because gdal and netCDF4 (and NCO) read the data array in (y-) reverse order, need to adjust orientation
        varin_rev = (varin[::-1])
        for yj in range(len(yin)):
            for xi in range(len(xin)):
                if((varin_rev[yj,xi] < valid_min) or (varin_rev[yj,xi] > valid_max)):
                    varin_rev[yj,xi] = 0
        varout[:,:] = offset + multiplier_Factor*(project_and_resample_Array(varin_rev, srs_geotrs, srs_projt, Nxin, Nyin, reference_netcdf))
        for tf in range (tInterpolation_frequency):
            ncOut.variables[varName][tk*tInterpolation_frequency+tf,:,:] = varout[::-1]         # reverse back the array for netCDF4
    ncIn.close()
    ncOut.close()
    #delete temp netcdf file
    #os.remove(temp_netcdf)



"""char transverse_mercator long form:
{'grid_mapping_name' : "transverse_mercator" ,'longitude_of_central_meridian' : -111. ,'false_easting' : 500000. , ' false_northing' : 0. , 'latitude_of_projection_origin' : 0.  , 'scale_factor_at_central_meridian' : 0.9996,
'longitude_of_prime_meridian' : 0. , 'semi_major_axis' : 6378137. , 'inverse_flattening' : 298.257222101 ,
'spatial_ref' : "PROJCS[\"NAD83 / UTM zone 12N  \",GEOGCS[\"NAD83\",DATUM[\"North_American_Datum_1983\",SPHEROID[\"GRS 1980\",63  78137,298.2572221010002,AUTHORITY[\"EPSG\",\"7019\"]],AUTHORITY[\"EPSG\",\"6269\  "]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433],AUTHORITY[\"EPSG  \",\"4269\"]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin  \",0],PARAMETER[\"central_meridian\",-111],PARAMETER[\"scale_factor\",0.9996],PA  RAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",0],UNIT[\"metre\"  ,1,AUTHORITY[\"EPSG\",\"9001\"]],AUTHORITY[\"EPSG\",\"26912\"]]" ,
'GeoTransform' : "432404.019091 30 0 4662392.4  4692 0 -30 " }
"""

tmAttributes = {'grid_mapping_name' : 'transverse_mercator' ,'longitude_of_central_meridian' : -111. ,
                            'false_easting' : 500000. , 'false_northing' : 0. , 'latitude_of_projection_origin' : 0.,
                            'scale_factor_at_central_meridian' : 0.9996, 'longitude_of_prime_meridian' : 0. ,
                            'semi_major_axis' : 6378137. , 'inverse_flattening' : 298.257222101 }

def lesser(x,y):
    if (x<y):
        return x
    else:
        return y

def greater(x,y):
    if (x>y):
        return x
    else:
        return y

"""
def project_netCDF_epsg(input_netcdf,  output_netcdf, varName, epsgCode):
"""
"""This projection assumes the source spatial reference is known
        i.e. GDAL can read it and recognize it
        varName: is the variable of interest in the netCDF file for which the projection is made
"""
"""
    tmAttributes['longitude_of_central_meridian'] = float( 6* (utmZone - 1) + 3 - 180 )

    data_set = gdal.Open(input_netcdf, GA_ReadOnly)
    s_srs = data_set.GetProjection()
    data_set = None
    cmdString = "nccopy -k 3 "+input_netcdf+" "+output_netcdf
    callSubprocess(cmdString, 'copy netcdf')

    ncData = netCDF4.Dataset(output_netcdf,"r+", format='NETCDF4')
    xArray = ncData.variables['x'][:]
    yArray = ncData.variables['y'][:]

    outArrayX = numpy.zeros(len(xArray))
    outArrayY = numpy.zeros(len(yArray))

    for i in range(len(xArray)):
        outArrayX[i], dummyY = project_a_point_UTM(xArray[i],yArray[0],s_srs,utmZone)
    for j in range(len(yArray)):
        dummyX, outArrayY[j] = project_a_point_UTM(xArray[0],yArray[j],s_srs,utmZone)

    ncData.variables['x'][:]=outArrayX[:]
    ncData.variables['y'][:]=outArrayY[:]
    ncData.createVariable('transverse_mercator','c')
    ncData.variables['transverse_mercator'].setncatts(tmAttributes)
    ncData.variables[varName].setncattr('grid_mapping', 'transverse_mercator')

    ncData.close()
"""

def project_netCDF_UTM_NAD83(input_netcdf,  output_netcdf, varName, utmZone):
    """ This projection assumes the source spatial reference is known
        i.e. GDAL can read it and recognize it
        varName: is the variable of interest in the netCDF file for which the projection is made
    """
    tmAttributes['longitude_of_central_meridian'] = float( 6* (utmZone - 1) + 3 - 180 )

    data_set = gdal.Open(input_netcdf, GA_ReadOnly)
    s_srs = data_set.GetProjection()
    data_set = None
    cmdString = "nccopy -k 3 "+input_netcdf+" "+output_netcdf
    callSubprocess(cmdString, 'copy netcdf')

    ncData = netCDF4.Dataset(output_netcdf,"r+", format='NETCDF4')
    xArray = ncData.variables['x'][:]
    yArray = ncData.variables['y'][:]

    outArrayX = numpy.zeros(len(xArray))
    outArrayY = numpy.zeros(len(yArray))

    for i in range(len(xArray)):
        outArrayX[i], dummyY = project_a_point_UTM(xArray[i],yArray[0],s_srs,utmZone)
    for j in range(len(yArray)):
        dummyX, outArrayY[j] = project_a_point_UTM(xArray[0],yArray[j],s_srs,utmZone)

    ncData.variables['x'][:]=outArrayX[:]
    ncData.variables['y'][:]=outArrayY[:]
    ncData.createVariable('transverse_mercator','c')
    ncData.variables['transverse_mercator'].setncatts(tmAttributes)
    ncData.variables[varName].setncattr('grid_mapping', 'transverse_mercator')

    ncData.close()

#Logan utm12 nad83 xminymin xmaxymax 432760.510, 4612686.409, 461700.887, 4662453.522
def subset_netCDF_to_referenceRaster(input_netcdf, reference_Raster, output_netcdf):
    """ this gives netcdf subset for reference_raster; to get the exact boundary of the
        reference_raster, the input and reference must have same resolution
        The coordinates of the bounding box are projected to the netcdf projection
    To Do: Boundary check-> check if the bounding box of subset raster is
               within the input_netcdf's boundary
    Boundary parameters extracted from reference_Raster
    """
    data_set = gdal.Open(reference_Raster, GA_ReadOnly)
    s_srs = data_set.GetProjection()
    geo_transform = data_set.GetGeoTransform()
    # use ulx uly lrx lry
    xmin = geo_transform[0]
    ymax = geo_transform[3]
    dx = geo_transform[1]
    dy = geo_transform[5]
    xmax = xmin + dx * data_set.RasterXSize
    ymin = ymax + dy* data_set.RasterYSize          # dy is -ve
    data_set = None
    data_set = gdal.Open(input_netcdf, GA_ReadOnly)
    t_srs = data_set.GetProjection()
    geo_transform = data_set.GetGeoTransform()
    dxT = geo_transform[1]
    dyT = -1*(geo_transform[5])       #dy is -ve
    data_set = None

    nwX, nwY = project_a_point_srs(xmin,ymax,s_srs,t_srs)
    neX, neY = project_a_point_srs(xmax,ymax,s_srs,t_srs)
    seX, seY = project_a_point_srs(xmax,ymin,s_srs,t_srs)
    swX, swY = project_a_point_srs(xmin, ymax,s_srs,t_srs)

    #take the bigger cell size for buffer
    if(dx > dxT):
        dxT = dx
    if(-1*dy > dyT):
        dyT = -1*dy
    #add a buffer around the boundary
    xmin = lesser(nwX,swX) - 2*dxT
    xmax = greater(seX,neX) + 2*dxT
    ymin = lesser(swY,seY) - 2*dyT
    ymax = greater(nwY,neY) + 2*dyT

    cmdString = "ncea -4 -d y,"+str(ymin)+","+str(ymax)+" -d x,"+str(xmin)+","+str(xmax)+" -O "\
                 +input_netcdf+" "+output_netcdf
    callSubprocess(cmdString, 'subset netcdf')
    data_set = None

#this gives netcdf subset along the time dimension
def get_netCDF_subset_TimeDim(input_netcdf, output_netcdf, timedimName, starttimeIndex, endtimeIndex):
    #Note: the time bounds are given as index of the time array
    cmdString = "ncea -4 -d "+timedimName+","+str(starttimeIndex)+","+str(endtimeIndex)+" "\
                 +input_netcdf+" "+output_netcdf
    callSubprocess(cmdString, 'subset netcdf')


def reverse_netCDF_yaxis_and_Inherit_fillValue(input_netcdf, output_netcdf, varName):
    """
    """
    ref_data = gdal.Open(input_netcdf, GA_ReadOnly)
    ref_proj = ref_data.GetProjection()
    ref_geotrs = ref_data.GetGeoTransform()
    Ncols = ref_data.RasterXSize
    Nrows = ref_data.RasterYSize
    inband = ref_data.GetRasterBand(1)
    nodata = inband.GetNoDataValue()
    array = inband.ReadAsArray()
    dType = inband.DataType

    output_netcdf1 = "Temp"+input_netcdf
    out_data = gdal.GetDriverByName('NetCDF').Create(output_netcdf1, Ncols, Nrows, 1, dType,["FORMAT=NC4"])
    out_data.SetGeoTransform(ref_geotrs)
    out_data.SetProjection(ref_proj)
    outband = out_data.GetRasterBand(1)
    #outband.SetNoDataValue(nodata)
    array_rev = array[::-1]
    outband.WriteArray(array_rev)
    outband.FlushCache()
    ref_data = None
    out_data = None

    ncIn = netCDF4.Dataset(output_netcdf1,"r+")
    yin = ncIn.variables['y'][:]
    yin_rev = yin[::-1]
    ncIn.variables['y'][:] = yin_rev[:]
    ncIn.close()

    #6.23.15 setnodatavalue() results in loss of variable-attributes. The following lines are intended to remedy it
    if dType == 1:
        att_type = 'b'
    elif dType == 2:
        att_type = 'u'
    elif dType == 3:
        att_type = 's'
    elif dType == 4:
        att_type = 'ui'
    elif dType == 5:
        att_type = 'i'
    elif dType == 6:
        att_type = 'f'
    elif dType == 7:
        att_type = 'd'
    else:
        print("error unknown data type encountered. Function exits before completion ")
        return
    print(dType)
    val =  nodata/1
    print (val)
    print(nodata)
    cmdString = "ncatted -a _FillValue,"+varName+",m,"+att_type+","+"-32768.0"+" "+ output_netcdf1+" "+output_netcdf
    callSubprocess(cmdString, 'edit attributed')

def reverse_netCDF_yaxis(input_netcdf, output_netcdf):
    """
    """
    ref_data = gdal.Open(input_netcdf, GA_ReadOnly)
    ref_proj = ref_data.GetProjection()
    ref_geotrs = ref_data.GetGeoTransform()
    Ncols = ref_data.RasterXSize
    Nrows = ref_data.RasterYSize
    inband = ref_data.GetRasterBand(1)
    #nodata = inband.GetNoDataValue()
    array = inband.ReadAsArray()
    dType = inband.DataType

    out_data = gdal.GetDriverByName('NetCDF').Create(output_netcdf, Ncols, Nrows, 1, dType,["FORMAT=NC4"])
    out_data.SetGeoTransform(ref_geotrs)
    out_data.SetProjection(ref_proj)
    outband = out_data.GetRasterBand(1)
    #outband.SetNoDataValue(nodata)
    array_rev = array[::-1]
    outband.WriteArray(array_rev)
    outband.FlushCache()
    ref_data = None
    out_data = None

    ncIn = netCDF4.Dataset(output_netcdf,"r+")
    yin = ncIn.variables['y'][:]
    yin_rev = yin[::-1]
    ncIn.variables['y'][:] = yin_rev[:]
    ncIn.close()


def reverse_netCDF_yaxis_and_rename_Variable(input_netcdf, output_netcdf, inut_varname='Band1', output_varname='Band1'):
    """
    """
    ncIn = netCDF4.Dataset(input_netcdf,"r") # format='NETCDF4')
    xin = ncIn.variables['x'][:]
    yin = ncIn.variables['y'][:]

    ncOut = netCDF4.Dataset(output_netcdf,"w", format='NETCDF4')
    ncOut.createDimension('y',len(yin))
    ncOut.createDimension('x', len(xin))
    dataType = ncIn.variables['x'].datatype
    vardataType = ncIn.variables[inut_varname].datatype
    ncOut.createVariable('y',dataType,('y',))
    ncOut.createVariable('x',dataType,('x',))
    ncOut.variables['y'][:] = yin[::-1]
    ncOut.variables['x'][:] = xin[:]
    ncOut.createVariable(output_varname,vardataType,('y','x',))
    #Copy attributes
    varAtts = ncIn.ncattrs()
    attDict = dict.fromkeys(varAtts)
    for attName in varAtts:
        attDict[attName] = getattr(ncIn,attName)
    ncOut.setncatts(attDict)
    #variable
    varAtts = ncIn.variables[inut_varname].ncattrs()
    attDict = dict.fromkeys(varAtts)
    for attName in varAtts:
        attDict[attName] = getattr(ncIn.variables[inut_varname],attName)
    ncOut.variables[output_varname].setncatts(attDict)
    #grid mapping var
    gridMapping = attDict['grid_mapping']
    ncOut.createVariable(gridMapping,'c',())
    varAtts = ncIn.variables[gridMapping].ncattrs()
    attDict = dict.fromkeys(varAtts)
    for attName in varAtts:
        attDict[attName] = getattr(ncIn.variables[gridMapping],attName)
    ncOut.variables[gridMapping].setncatts(attDict)
    ncOut.variables[gridMapping].setncatts(attDict)
    #dim variables
    xAtts = ncIn.variables['x'].ncattrs()
    attDict = dict.fromkeys(xAtts)
    for attName in xAtts:
        attDict[attName] = getattr(ncIn.variables['x'],attName)
    ncOut.variables['x'].setncatts(attDict)
    yAtts = ncIn.variables['y'].ncattrs()
    attDict = dict.fromkeys(yAtts)
    for attName in yAtts:
        attDict[attName] = getattr(ncIn.variables['y'],attName)
    ncOut.variables['y'].setncatts(attDict)
    array = ncIn.variables[inut_varname][:]
    ncOut.variables[output_varname][:] = array[::-1]

    ncIn.close()
    ncOut.close()


"""
resample_netCDF_to_referenceNetCDF('SpawnProj_2010.nc','Spawn17.nc','prcp','Res2.nc')
"""
def resample_netCDF_to_referenceNetCDF(input_netcdf, reference_netcdf, output_netcdf, varName, varUnits=" ", multiplier_Factor = 1, offset = 0):
    """This re-grids a netcdf to target/reference netcdf resolution
        the extent and cell size of the output_netcdf will be that of the reference netcdf
        the input netcdf must have the same projection as that of the reference netcdf
    Note: unlike in all other functions, the reference is acutally netcdf (not raster)
    Input coordinates are time, y, x

    Warning: Works only if the target boundary is within the input boundary & the coordinates directions are
    the same, i.e. y increasing / decreasing
    (GDAL generated netcdf have y inverted)
    ToDO: Check GDAL netcdf generation
    ToDO: Check boundary
    """
    #Add dummy dimensions and variables
    temp_netcdf = "temp"+output_netcdf
    cmdString = "ncrename -d x,x_2 -d y,y_2 -v x,x_2 -v y,y_2 -v"+varName+","+varName+"_2 "+\
                 input_netcdf+" "+temp_netcdf
    callSubprocess(cmdString, 'copy netcdf with rename old dimensions')

    ncRef = netCDF4.Dataset(reference_netcdf,"r") # format='NETCDF4')
    xout = ncRef.variables['x'][:]
    yout = ncRef.variables['y'][:]
    ncRef.close()

    ncIn = netCDF4.Dataset(input_netcdf,"r") # format='NETCDF4')
    xin = ncIn.variables['x'][:]
    yin = ncIn.variables['y'][:]
    dx = abs(xin[1] - xin[0])       # regular grid so the dx are same over the grid
    dy = abs(yin[1] - yin[0])

    ncOut = netCDF4.Dataset(temp_netcdf,"r+", format='NETCDF4')
    ncOut.createDimension('y',len(yout))
    ncOut.createDimension('x', len(xout))
    dataType = ncRef.variables['x'].datatype
    vardataType = ncIn.variables[varName].datatype
    ncOut.createVariable('y',dataType,('y',))
    ncOut.createVariable('x',dataType,('x',))
    ncOut.variables['y'][:] = yout[:]
    ncOut.variables['x'][:] = xout[:]
    ncOut.createVariable(varName,vardataType,('time','y','x',))
    #Copy attributes
    varAtts = ncIn.variables[varName].ncattrs()
    attDict = dict.fromkeys(varAtts)
    for attName in varAtts:
        attDict[attName] = getattr(ncIn.variables[varName],attName)
        if attName == 'units':
            if varUnits != " ":
                attDict[attName] = varUnits
    ncOut.variables[varName].setncatts(attDict)
    xAtts = ncIn.variables['x'].ncattrs()
    attDict = dict.fromkeys(xAtts)
    for attName in xAtts:
        attDict[attName] = getattr(ncIn.variables['x'],attName)
    ncOut.variables['x'].setncatts(attDict)
    yAtts = ncIn.variables['y'].ncattrs()
    attDict = dict.fromkeys(yAtts)
    for attName in yAtts:
        attDict[attName] = getattr(ncIn.variables['y'],attName)
    ncOut.variables['y'].setncatts(attDict)
    ncOut.close()
    #delete the old variables
    cmdString = "ncks -4 -C -x -v x_2,y_2,"+varName+"_2 "+\
                 temp_netcdf+" "+output_netcdf
    callSubprocess(cmdString, 'delete old dimensions')

    #re-open file to write re-gridded data
    ncOut = netCDF4.Dataset(output_netcdf,"r+") #, format='NETCDF4')
    varin = numpy.zeros((len(yin),len(xin)),dtype=vardataType)
    varout = numpy.zeros((len(yout),len(xout)),dtype=vardataType)
    timeLen = len(ncIn.dimensions['time'])
    yLen = len(yout)
    xLen = len(xout)
    for tk in range(timeLen):
        varin[:,:] = ncIn.variables[varName][tk,:,:]
        for yi in range(yLen):
            y1 = int(numpy.floor(abs(yout[yi]-yin[0])/dy))    # abs to make sure
            y2 = y1+1
            for xj in range(xLen):
                x1 = int(numpy.floor(abs(xout[xj]-xin[0])/dx))
                x2 =x1+1
                points = [(yin[y1],xin[x1],varin[y1,x1]),(yin[y1],xin[x2],varin[y1,x2]),(yin[y2],xin[x2],varin[y2,x2]),
                          (yin[y2],xin[x1],varin[y2,x1])]
                varout[yi,xj] = offset + multiplier_Factor * (bilinear_interpolation_with_points_outside_Rectangle(yout[yi],xout[xj],points)) # bilinear_interpolation(yout[yi],xout[xj],points)
        ncOut.variables[varName][tk,:,:] = varout[:,:]
    ncIn.close()
    ncOut.close()
    #delete temp netcdf file
    #os.remove(temp_netcdf)


#This function from: http://stackoverflow.com/questions/8661537/how-to-perform-bilinear-interpolation-in-python
def bilinear_interpolation(x, y, points):
    '''Interpolate (x,y) from values associated with four points.

    The four points are a list of four triplets:  (x, y, value).
    The four points can be in any order.  They should form a rectangle.

        >>> bilinear_interpolation(12, 5.5,[(10, 4, 100),(20, 4, 200),(10, 6, 150),(20, 6, 300)])
        165.0

    '''
    # See formula at:  http://en.wikipedia.org/wiki/Bilinear_interpolation

    points = sorted(points)               # order points by x, then by y
    (x1, y1, q11), (_x1, y2, q12), (x2, _y1, q21), (_x2, _y2, q22) = points

    if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
        raise ValueError('points do not form a rectangle')
    if not x1 <= x <= x2 or not y1 <= y <= y2:
        #raise ValueError
        warnString = 'warning! point ('+repr(x)+', '+repr(y)+') not within the rectangle: '+repr(x1)+' '+repr(x2)+', '+\
            repr(y1)+' '+repr(y2)
        raise ValueError(warnString)

    return (q11 * (x2 - x) * (y2 - y) +
            q21 * (x - x1) * (y2 - y) +
            q12 * (x2 - x) * (y - y1) +
            q22 * (x - x1) * (y - y1)
           ) / ((x2 - x1) * (y2 - y1) + 0.0)



#logan boundaries for test  -111.818, 42.113, -111.457, 41.662
def get_netCDFSubset_Geographic(input_netcdf, output_netcdf, lonname, latname, lonmin, lonmax, latmin, latmax):
    #similar to get DEM subset this function gets subset of netcdf in
    #geographic coordinate system; this enables dealing with projection differences between the source netcdf file
    #and the target watershed raster
    #subsettting before projecting to the target reference system avoids dealing with large file
    #however it works only if the input_netcdf has lat/lon coordinate dimenstions
    """ Note: upper left (ul) considered origin, i.e. xmin, ymax
    parameters passed as ulx uly lrx lry (xmin, ymax, xmax, ymin)
    The arguments are in decimal degrees  is in Geographic CS
    latname and lonname: name used to refer the geographic coordinates
    """
    cmdString = "ncea -4 -d "+latname+","+str(latmin)+","+str(latmax)+" -d "+lonname+","+str(lonmin)+","+str(lonmax)+" "\
                 +input_netcdf+" "+output_netcdf
    callSubprocess(cmdString, 'subset netcdf')



#gets the subset of netcdf withing the reference raster
#values are resampled to the resolution of the reference raster
def project_and_subset_netCDF2D(input_netcdf, reference_raster, output_netcdf, resample='bilinear'):
    """
    :param input_raster:
    :param reference_raster:
    :param output_raster:
    :return:
    For images use nearest neighbor interpolation; else pass the method required

    The target extent parameters -te xmin ymin xmax ymax  may be needed to provide a
    region where the destination projection is valid
    projecting a large dataset (e.g. CONUS) into local projections (e.g. NAD 12N) may fail
    because gdal checks the validity of the projection method for the entire region
    """
    data_set = gdal.Open(input_netcdf)
    target_srs = data_set.GetProjection() #osr.SpatialReference(wkt
    #target_srs.ImportFromWkt(data_set.GetPrjectionRef())
    data_set = None
    srsprjFile = 'srsprj.prf'
    prjFilep = open(srsprjFile,'w')
    prjFilep.write(target_srs)
    prjFilep.close()

    data_set = gdal.Open(reference_raster)
    target_srs = data_set.GetProjection() #osr.SpatialReference(wkt
    geo_transform = data_set.GetGeoTransform()
    dx = geo_transform[1]
    dy = geo_transform[5]
    xmin = geo_transform[0]
    ymax = geo_transform[3]
    xmax = xmin + dx * data_set.RasterXSize
    ymin = ymax + dy* data_set.RasterYSize
    data_set = None
    tprjFile = 'destprj.prf'
    prjFilep = open(tprjFile,'w')
    prjFilep.write(target_srs)
    prjFilep.close()

    cmdString = "gdalwarp -s_srs "+srsprjFile+" -t_srs "+tprjFile+" -te "\
                +str(xmin)+" "+str(ymin)+" "+str(xmax)+" "+str(ymax)+" -tr "\
                +str(dx)+" "+str(-1*dy)+" -r "+resample+" -overwrite "+input_netcdf+" tempraster.tif"        #+output_netcdf
    callSubprocess(cmdString, "create intermediate tiff file ")

    cmdString = "gdal_translate -of NetCDF tempraster.tif "+output_netcdf
    callSubprocess(cmdString, "project and clip NetCDF")
     #delete intermediate file
    os.remove('tempRaster.tif')


#This one works only for 2D; with 3D netcdf inputs it images of 2D, band1, band2,...
def project_and_subset_netCDF2D_Image(input_netcdf, reference_raster, output_netcdf):    # resample='near'):
    """
    :param input_raster:
    :param reference_raster:
    :param output_raster:
    :return:
    For images leave the default nearest neighbor interpolation; else pass the method required
    """
    srs_data = gdal.Open(input_netcdf, GA_ReadOnly)
    srs_proj = srs_data.GetProjection() #osr.SpatialReference(wkt
    srs_geotrans = srs_data.GetGeoTransform()

    ref_data = gdal.Open(reference_raster, GA_ReadOnly)
    ref_proj = ref_data.GetProjection()
    ref_geotrs = ref_data.GetGeoTransform()
    Ncols = ref_data.RasterXSize
    Nrows = ref_data.RasterYSize

    out_data = gdal.GetDriverByName('NetCDF').Create(output_netcdf, Ncols, Nrows, 1, GDT_Byte)
    out_data.SetGeoTransform(ref_geotrs)
    out_data.SetProjection(ref_proj)

    gdal.ReprojectImage(srs_data,out_data,srs_proj,ref_proj, GRA_NearestNeighbour)
    out_data = None

def project_a_point_UTM(xcoord, ycoord, s_srs, utmZone):
    s_srsT = osr.SpatialReference()
    s_srsT.ImportFromWkt(s_srs)
    t_srsT = osr.SpatialReference()
    srsString = "+proj=utm +zone=" +str(utmZone)+ " +ellps=GRS80 +datum=NAD83 +units=m "
    t_srsT.ImportFromProj4(srsString)
    #t_srsT.ImportFromWkt(t_srs)
    transform = osr.CoordinateTransformation(s_srsT, t_srsT)
    pointC = ogr.Geometry(ogr.wkbPoint)
    pointC.SetPoint_2D(0,float(xcoord), float(ycoord))
    pointC.Transform(transform)
    xproj = pointC.GetX()
    yproj = pointC.GetY()
    return xproj, yproj

def project_a_point_srs(xcoord, ycoord, s_srs, t_srs):
    s_srsT = osr.SpatialReference()
    s_srsT.ImportFromWkt(s_srs)
    t_srsT = osr.SpatialReference()
    t_srsT.ImportFromWkt(t_srs)
    transform = osr.CoordinateTransformation(s_srsT, t_srsT)
    pointC = ogr.Geometry(ogr.wkbPoint)
    pointC.SetPoint_2D(0,float(xcoord), float(ycoord))
    pointC.Transform(transform)
    xproj = pointC.GetX()
    yproj = pointC.GetY()
    return xproj, yproj

def raster_to_netCDF(input_raster, output_netcdf):
    cmdString = "gdal_translate -of netCDF -co \"FORMAT=NC4\" "+input_raster+" "+output_netcdf
    callSubprocess(cmdString, 'raster to netcdf')

#This concatenates netcdf files along the time dimension
def concatenate_netCDF(input_netcdf1, input_netcdf2, output_netcdf):
    """To  Do: may need to specify output no-data value
    """
    cmdString = "ncks --mk_rec_dmn time "+input_netcdf1+" tempNetCDF1.nc"
    callSubprocess(cmdString, "intermediate netcdf with record dimension")
    cmdString = "ncks --mk_rec_dmn time "+input_netcdf2+" tempNetCDF2.nc"
    callSubprocess(cmdString, "intermediate netcdf with record dimension")
    #
    cmdString = "ncrcat -4 tempNetCDF1.nc tempNetCDF2.nc "+output_netcdf
    callSubprocess(cmdString, "concatenate netcdf files")
    #delete intermediate files
    os.remove('tempNetCDF1.nc')
    os.remove('tempNetCDF2.nc')

#This combines (stitches) (spatially adjacent) netcdf files accross the spatial/horizontal dimensions
def combineNetCDFs(input_netcdf1, input_netcdf2, output_netcdf):
    """To  Do: may need to specify output no-data value
    """
    cmdString = "gdalwarp -of GTiff -overwrite "+input_netcdf1+" "+input_netcdf2+" tempRaster.tif"  #+output_raster
    callSubprocess(cmdString, "create intermediate raster file")
    #print 'done concatenating netcdfs'
    cmdString = "gdal_translate -of NetCDF tempRaster.tif "+output_netcdf
    callSubprocess(cmdString, "combine two netcdf files")
    #print 'done function'
    #delete intermediate file
    os.remove('tempRaster.tif')

def callSubprocess(cmdString, debugString):
    cmdargs = shlex.split(cmdString)
    debFile = open('debug_file.txt', 'w')
    debFile.write('Starting %s \n' % debugString)
    retValue = subprocess.call(cmdString,stdout=debFile,shell=True)               # use shell=True with a single string of commands; shell=False with list of strings with first element the executable
    if (retValue==0):
        debFile.write('%s Successful\n' % debugString)
        debFile.close()
    else:
        debFile.write('There was error in %s\n' % debugString)
        debFile.close()


#This function from: http://stackoverflow.com/questions/8661537/how-to-perform-bilinear-interpolation-in-python
#_added a line for points outside rectangle to take boundary values
def bilinear_interpolation_with_points_outside_Rectangle(x, y, points):
    '''Interpolate (x,y) from values associated with four points.
    The four points are a list of four triplets:  (x, y, value).
    The four points can be in any order.  They should form a rectangle.
        >>> bilinear_interpolation(12, 5.5,[(10, 4, 100),(20, 4, 200),(10, 6, 150),(20, 6, 300)])
        165.0
    '''
    # See formula at:  http://en.wikipedia.org/wiki/Bilinear_interpolation
    points = sorted(points)               # order points by x, then by y
    (x1, y1, q11), (_x1, y2, q12), (x2, _y1, q21), (_x2, _y2, q22) = points

    if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
        raise ValueError('points do not form a rectangle')
    if not x1 <= x <= x2 or not y1 <= y <= y2:
        #raise ValueError
        warnString = 'warning! point ('+repr(x)+', '+repr(y)+') not within the rectangle: '+repr(x1)+' '+repr(x2)+', '+\
            repr(y1)+' '+repr(y2)
        #raise ValueError(warnString)
        """TZG added this 12.5.14 for inspection """
        print (warnString)
        #use boundary values
        if (x < x1):
            x = x1
        if(x > x2):
            x = x2
        if(y < y1 ):
            y = y1
        if(y > y2):
            y = y2

    return (q11 * (x2 - x) * (y2 - y) +
            q21 * (x - x1) * (y2 - y) +
            q12 * (x2 - x) * (y - y1) +
            q22 * (x - x1) * (y - y1)
           ) / ((x2 - x1) * (y2 - y1) + 0.0)


