"""
This is used for the older version of processing the ueb model input
"""
try:
    from osgeo import gdal, osr, ogr
except:
    import gdal, osr, ogr
from gdalconst import *
import shlex
import subprocess
import numpy
import netCDF4
import math
from datetime import datetime, timedelta


def subset_and_concatenate_netCDF_from_multple_nc_NLDAS(file_prefix, wsName, output_netcdf, leftX, topY, rightX, bottomY,
                      startDateTime, endDateTime, dT=1, in_Xcoord = 'lon_110', in_Ycoord='lat_110',inout_timeName = 'time'):

    """
    Subsets and combines multiple netcdf files
    for nldas forcing, with multiple time steps (e.g., organized in monthly files)
    should already have time dim. for ncrcat, made record dim by ncks
    e.g.:
    Logan leftX=-112.0, topY=42.3, rightX=-111.0, bottomY=41.6, startYear=2009, endYear=2010
    for nldas data with time dim (e.g., previously concatenated in time dim)
    """
    startYear = datetime.strptime(startDateTime,"%Y/%m/%d %H").year
    endYear = datetime.strptime(endDateTime,"%Y/%m/%d %H").year
    startMonth = datetime.strptime(startDateTime,"%Y/%m/%d %H").month
    endMonth = datetime.strptime(endDateTime,"%Y/%m/%d %H").month
    startDay =  datetime.strptime(startDateTime,"%Y/%m/%d %H").timetuple().tm_yday        #start date = day of year for 2010
    endDay   =  startDay + (datetime.strptime(endDateTime,"%Y/%m/%d %H") - datetime.strptime(startDateTime,"%Y/%m/%d %H")).days          # end date = day of year for 2011 + 365

    #print(startYear)
    #print(endYear)

    for year in range(startYear, endYear+1):
        for month in range(1, 13):
            if month < 10:
                monthS = '0'+str(month)
            else:
                monthS = str(month)
            cmdString = "for i in "+file_prefix+"*"+str(year)+monthS+"*.nc; do ncea -d "+in_Xcoord+","+str(leftX)+","+str(rightX)\
                    +" -d "+in_Ycoord+","+str(bottomY)+","+str(topY)+" -O $i "+wsName+"_WS_$i; done"      #+subdir+"\/"
            callSubprocess(cmdString, 'subset nc files for year '+str(year))

    cmdString = "for i in "+wsName+"_WS_"+"*.nc; do ncks --mk_rec_dmn "+inout_timeName+" -O $i R_$i; done"
    callSubprocess(cmdString, "intermediate netcdf with record dimension")


    cmdString = "ncrcat -4 -H -h -O  R_"+wsName+"*.nc -o concat_"+output_netcdf                     #-H don't append input file list -h don't append history
    callSubprocess(cmdString, "concatenate netcdf files")

    hD = int(24/dT)
    starttimeIndex = (startDay-1) * hD                #4.24.17 correct index for e.g., oct 1 0 hrs     #starttimeIndex = (startDay-1) * hD + 1
    endtimeIndex = (endDay-1) * hD - 1                   # corr. indx for exclusive of last step e.g., oct 1 0 hrs

    #starttimeIndex = startDay * hD
    #endtimeIndex = endDay * hD
    #starttimeIndex = (startDay-1) * hD + 1               #4.24.17 correct index for e.g., oct 1 0 hrs
    #endtimeIndex = (endDay-1) * hD                   # corr. indx for exclusive of last step e.g., oct 1 0 hrs

    print(starttimeIndex)
    print(endtimeIndex)
    cmdString = "ncea -4 -H -O -d "+inout_timeName+","+str(starttimeIndex)+","+str(endtimeIndex)+" concat_"\
                 +output_netcdf+" "+output_netcdf
    callSubprocess(cmdString, 'subset netcdf in time')


    #delete intermediate files
    cmdString = "rm -f R_*.nc "
    callSubprocess(cmdString, "delete intermediate files")
    cmdString = "rm -f "+wsName+"_WS_"+"*.nc "
    callSubprocess(cmdString, "delete intermediate files")

    #cmdString = "DEL "+"R_"+wsName+"*.nc"
    #callSubprocess(cmdString, "delete intermediate files")
    #cmdString = "DEL "+wsName+"*.nc"
    #callSubprocess(cmdString, "delete intermediate files")
    #os.remove("R_*.nc")


#Logan leftX=-112.0, topY=42.3, rightX=-111.0, bottomY=41.6, startYear=2009, endYear=2010
# RBC leftX, topY, rightX, bottomY =  -111.819, 40.829, -111.736, 40.767
#for nldas data
def create_netCDF_from_multple_nc2(file_prefix, wsName, output_netcdf, leftX, topY, rightX, bottomY,
                      startDateTime, endDateTime, dT=1, in_Xcoord = 'lon_110', in_Ycoord='lat_110',out_timeName = 'time'):

    """
    Subsets and combines multiple netcdf files
    written for CBRFC temp/prec data sets that are commpressed with .gz
    before calling this funcion do gunzip *.nc.gz
    """  
    startYear = datetime.strptime(startDateTime,"%Y/%m/%d %H").year
    endYear = datetime.strptime(endDateTime,"%Y/%m/%d %H").year
    startMonth = datetime.strptime(startDateTime,"%Y/%m/%d %H").month
    endMonth = datetime.strptime(endDateTime,"%Y/%m/%d %H").month
    startDay =  datetime.strptime(startDateTime,"%Y/%m/%d %H").timetuple().tm_yday        #start date = day of year for 2010
    endDay   =  startDay + (datetime.strptime(endDateTime,"%Y/%m/%d %H") - datetime.strptime(startDateTime,"%Y/%m/%d %H")).days          # end date = day of year for 2011 + 365

    #print(startYear)  print(endYear)
    for year in range(startYear, endYear+1):
        for month in range(1, 13):
            if month < 10:
                monthS = '0'+str(month)
            else:
                monthS = str(month)
            cmdString = "for i in "+file_prefix+"*"+str(year)+monthS+"*.nc; do ncea -d "+in_Xcoord+","+str(leftX)+","+str(rightX)\
                    +" -d "+in_Ycoord+","+str(bottomY)+","+str(topY)+" $i -O "+"R_"+wsName+"_$i ;done"      #+subdir+"\/"
            callSubprocess(cmdString, 'subset nc files for year '+str(year))

    #cmdString = "for %i in ("+wsName+"*.nc) do  (ncks --mk_rec_dmn "+out_timeName+" %i R_%i)"
    #callSubprocess(cmdString, "intermediate netcdf with record dimension")

    cmdString = "ncecat  -u \""+out_timeName+"\" -H -h -O R_"+wsName+"*.nc -o concat_"+output_netcdf                #-H don't append input file list -h don't append history
    callSubprocess(cmdString, "concatenate netcdf files")

    hD = int(24/dT)
    starttimeIndex = (startDay-1) * hD                #4.24.17 correct index for e.g., oct 1 0 hrs     #starttimeIndex = (startDay-1) * hD + 1
    endtimeIndex = (endDay-1) * hD - 1                   # corr. indx for exclusive of last step e.g., oct 1 0 hrs
    #print(starttimeIndex)     print(endtimeIndex)
    cmdString = "ncea -4 -d "+out_timeName+","+str(starttimeIndex)+","+str(endtimeIndex)+" -H -h -O concat_"\
                 +output_netcdf+" -o "+output_netcdf
    callSubprocess(cmdString, 'subset netcdf in time')
    #delete intermediate files
    cmdString = "rm -f R_"+wsName+"*.nc"
    callSubprocess(cmdString, "delete intermediate files")
    #os.remove("R_*.nc")


#leftX=-112.0, topY=42.3, rightX=-111.0, bottomY=41.6, startYear=2009, endYear=2010
# for cbrfc forcing
def create_netCDF_from_multple_nc(file_prefix, wsName, output_netcdf, var_name, leftX, topY, rightX, bottomY, startDateTime, endDateTime, dT=3, out_timeName = 'time'):
    """
    Subsets and combines multiple netcdf files
    written for CBRFC temp/prec data sets
    for datasets that are commpressed with .gz
    """
    startYear = datetime.strptime(startDateTime,"%Y/%m/%d %H").year
    endYear = datetime.strptime(endDateTime,"%Y/%m/%d %H").year
    startDay =  datetime.strptime(startDateTime,"%Y/%m/%d %H").timetuple().tm_yday        #start date = day of year e.g. for 274 for oct 1 2010
    endDay   =  startDay + (datetime.strptime(endDateTime,"%Y/%m/%d %H") - datetime.strptime(startDateTime,"%Y/%m/%d %H")).days          # end date = number of days b/n end-start + start


    for year in range(startYear, endYear+1):
        cmdString = 'for i in '+file_prefix+'_'+str(year)+'*.nc.gz; do gunzip -cdfq $i>${i//.gz} && gdal_translate -projwin '+str(leftX)+' '+str(topY)+' '+str(rightX)+' '+str(bottomY)\
                   +' -of netCDF -co "FORMAT=NC4" NETCDF:\"${i//.gz}\":'+var_name+' '+wsName+'_'+var_name+'_${i//.gz} && rm -f ${i//.gz}; done'     #delete the .nc immediately to save space
        #print(cmdString)
        callSubprocess(cmdString, 'subset nc files for year'+str(year))

    cmdString = 'for i in '+wsName+'_'+var_name+'*.nc; do  ncks --mk_rec_dmn '+out_timeName+' $i R_$i; done'
    callSubprocess(cmdString, "intermediate netcdf with record dimension")

    cmdString = "ncrcat -H -h -O R_*.nc -o concat_"+output_netcdf                     #-H don't append input file list -h don't append history
    callSubprocess(cmdString, "concatenate netcdf files")

    hD = int(24/dT)                  #number of time-steps in a day
    starttimeIndex = (startDay-1) * hD                #4.24.17 correct index for e.g., oct 1 0 hrs     #starttimeIndex = (startDay-1) * hD + 1
    endtimeIndex = (endDay-1) * hD - 1                   # corr. indx for exclusive of last step e.g., oct 1 0 hrs
    cmdString = "ncea -4 -d "+out_timeName+","+str(starttimeIndex)+","+str(endtimeIndex)+" -H -h -O concat_"\
                 +output_netcdf+" -o "+output_netcdf
    callSubprocess(cmdString, 'subset netcdf in time')


    #delete intermediate files
    cmdString = "rm -f R_*.nc "
    callSubprocess(cmdString, "delete intermediate files")
    cmdString = "rm -f "+wsName+"_"+var_name+"*.nc"
    callSubprocess(cmdString, "delete intermediate files")
    #os.remove("R_*.nc")


#leftX=-112.0, topY=42.3, rightX=-111.0, bottomY=41.6,
#for LAI
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
#for data from cbrfc
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


def compute_average_windSpeed(input_netcdfU, varNameU, input_netcdfV, varNameV, output_netcdfW, varNameW ):            #output_netcdfVP, varNameVP,
    """This re-grids a netcdf to target/reference resolution
    Input coordinates are time, y, x
    Warning: Works only if the target boundary is within the input boundary & the coordinates directions are
    the same, i.e. y increasing / decreasing """

    #Copy dimensions and variables
    """temp_netcdf = 'temp'+output_netcdfRH
    cmdString = "nccopy -k 4 "+input_netcdfT+" "+temp_netcdf             #output_netcdf
    callSubprocess(cmdString, 'copy netcdf with dimensions')"""

    #delete old variable
    cmdString = "ncks -4 -C -h -O -x -v "+varNameU+" "+input_netcdfU+" "+output_netcdfW
    callSubprocess(cmdString, 'delete old/reference variable')

    ncInU = netCDF4.Dataset(input_netcdfU,"r") # format='NETCDF4')
    vardataType = ncInU.variables[varNameU].datatype
    ref_grid_mapping = 'grid mmapping'    #getattr(ncInU.variables[varNameU],'grid_mapping')
    timeLen = len(ncInU.dimensions['time'])

    ncInV = netCDF4.Dataset(input_netcdfV,"r") # format='NETCDF4')

    ncOut = netCDF4.Dataset(output_netcdfW,"r+", format='NETCDF4')
    ncOut.createVariable(varNameW,vardataType,('time','y','x',))
    attDict = {'name':varNameW, 'long_name':'Average wind speed at height 10 m'}
    attDict['units'] = 'm/s'
    attDict['grid_mapping'] = ref_grid_mapping
    ncOut.variables[varNameW].setncatts(attDict)

    #varin = numpy.zeros((len(yin),len(xin)),dtype=vardataType)
    #varout = numpy.zeros((len(yout),len(xout)),dtype=vardataType)
    for tk in range(timeLen):
        varinU = ncInU.variables[varNameU][tk,:,:]
        varinV = ncInV.variables[varNameV][tk,:,:]
        Wave1 = varinU*varinU
        Wave2 = varinV*varinV
        Wave3 = Wave1+Wave2
        Wave = numpy.sqrt(Wave3)
        ncOut.variables[varNameW][tk,:,:] = Wave[:,:]

    ncInU.close()
    ncInV.close()
    ncOut.close()
    #delete temp netcdf file
    #os.remove(temp_netcdf)


def compute_relative_humidity_from_SpecficHumidity(input_netcdfT, varNameT, input_netcdfP, varNameP, input_netcdf_q, varNameq, output_netcdfRH, varNameRH ):            #output_netcdfVP, varNameVP,
    """This re-grids a netcdf to target/reference resolution
    Input coordinates are time, y, x
    Warning: Works only if the target boundary is within the input boundary & the coordinates directions are
    the same, i.e. y increasing / decreasing """

    #Copy dimensions and variables
    """temp_netcdf = 'temp'+output_netcdfRH
    cmdString = "nccopy -k 4 "+input_netcdfT+" "+temp_netcdf             #output_netcdf
    callSubprocess(cmdString, 'copy netcdf with dimensions')"""

    #delete old variable
    cmdString = "ncks -4 -C -h -O -x -v "+varNameT+" "+input_netcdfT+" "+output_netcdfRH
    callSubprocess(cmdString, 'delete old/reference variable')

    ncInT = netCDF4.Dataset(input_netcdfT,"r") # format='NETCDF4')
    vardataType = ncInT.variables[varNameT].datatype
    ref_grid_mapping = getattr(ncInT.variables[varNameT],'grid_mapping')
    timeLen = len(ncInT.dimensions['time'])

    ncInP = netCDF4.Dataset(input_netcdfP,"r") # format='NETCDF4')
    ncInq = netCDF4.Dataset(input_netcdf_q,"r") # format='NETCDF4')

    ncOut = netCDF4.Dataset(output_netcdfRH,"r+", format='NETCDF4')
    ncOut.createVariable(varNameRH,vardataType,('time','y','x',))
    attDict = {'name':varNameRH, 'long_name':'Relative humidity'}
    attDict['units'] = '-'
    attDict['grid_mapping'] = ref_grid_mapping
    ncOut.variables[varNameRH].setncatts(attDict)

    #varin = numpy.zeros((len(yin),len(xin)),dtype=vardataType)
    #varout = numpy.zeros((len(yout),len(xout)),dtype=vardataType)

    for tk in range(int(timeLen-1)):
        varinT = ncInT.variables[varNameT][tk,:,:]
        varinP = ncInP.variables[varNameP][tk,:,:]
        varinq = ncInq.variables[varNameq][tk,:,:]
        eSat = 611*(17.3*varinT/(varinT+237.3))
        print(eSat)
        RH = varinP*varinq /(eSat*(varinq+0.622))
        ncOut.variables[varNameRH][tk,:,:] = RH[:,:]

    ncInT.close()
    ncInP.close()
    ncInq.close()
    ncOut.close()
    #delete temp netcdf file
    #os.remove(temp_netcdf)


def compute_relative_humidity_from_VaporPressure(input_netcdfT, varNameT, input_netcdfVP, varNameVP, output_netcdfRH, varNameRH ):            #output_netcdfVP, varNameVP,
    """This re-grids a netcdf to target/reference resolution
    Input coordinates are time, y, x
    Warning: Works only if the target boundary is within the input boundary & the coordinates directions are
    the same, i.e. y increasing / decreasing """

    #Copy dimensions and variables
    """temp_netcdf = 'temp'+output_netcdfRH
    cmdString = "nccopy -k 4 "+input_netcdfT+" "+temp_netcdf             #output_netcdf
    callSubprocess(cmdString, 'copy netcdf with dimensions')"""

    #delete old variable
    cmdString = "ncks -4 -C -h -O -x -v "+varNameT+" "+input_netcdfT+" "+output_netcdfRH
    callSubprocess(cmdString, 'delete old/reference variable')

    ncInT = netCDF4.Dataset(input_netcdfT,"r") # format='NETCDF4')
    vardataType = ncInT.variables[varNameT].datatype
    ref_grid_mapping = getattr(ncInT.variables[varNameT],'grid_mapping')
    timeLen = len(ncInT.dimensions['time'])

    ncInq = netCDF4.Dataset(input_netcdfVP,"r") # format='NETCDF4')

    ncOut = netCDF4.Dataset(output_netcdfRH,"r+", format='NETCDF4')
    ncOut.createVariable(varNameRH,vardataType,('time','y','x',))
    attDict = {'name':varNameRH, 'long_name':'Relative humidity'}
    attDict['units'] = '-'
    attDict['grid_mapping'] = ref_grid_mapping
    ncOut.variables[varNameRH].setncatts(attDict)

    #varin = numpy.zeros((len(yin),len(xin)),dtype=vardataType)
    #varout = numpy.zeros((len(yout),len(xout)),dtype=vardataType)

    for tk in range(int(timeLen-1)):
        varinT = ncInT.variables[varNameT][tk,:,:]
        varinVP = ncInq.variables[varNameVP][tk,:,:]
        eSat = 611*(17.3*varinT/(varinT+237.3))
        print(eSat)
        RH = varinVP /eSat
        ncOut.variables[varNameRH][tk,:,:] = RH[:,:]

    ncInT.close()
    ncInq.close()
    ncOut.close()
    #delete temp netcdf file
    #os.remove(temp_netcdf)


def compute_vaporPressure_from_SpecficHumidity(input_netcdf_q, varNameq, input_netcdfP, varNameP,  output_netcdfVP, varNameVP ):            #output_netcdfVP, varNameVP,
    """This re-grids a netcdf to target/reference resolution
    Input coordinates are time, y, x
    Warning: Works only if the target boundary is within the input boundary & the coordinates directions are
    the same, i.e. y increasing / decreasing """

    #Copy dimensions and variables
    """temp_netcdf = 'temp'+output_netcdfRH
    cmdString = "nccopy -k 4 "+input_netcdfT+" "+temp_netcdf             #output_netcdf
    callSubprocess(cmdString, 'copy netcdf with dimensions')"""

    #delete old variable
    cmdString = "ncks -4 -C -h -O -x -v "+varNameq+" "+input_netcdf_q+" "+output_netcdfVP
    callSubprocess(cmdString, 'delete old/reference variable')

    ncInq = netCDF4.Dataset(input_netcdf_q,"r") # format='NETCDF4')
    vardataType = ncInq.variables[varNameq].datatype
    ref_grid_mapping = getattr(ncInq.variables[varNameq],'grid_mapping')
    timeLen = len(ncInq.dimensions['time'])

    ncInP = netCDF4.Dataset(input_netcdfP,"r") # format='NETCDF4')

    ncOut = netCDF4.Dataset(output_netcdfVP,"r+", format='NETCDF4')
    ncOut.createVariable(varNameVP,vardataType,('time','y','x',))
    attDict = {'name':varNameVP, 'long_name':'Vapor pressure'}
    attDict['units'] = '-'
    attDict['grid_mapping'] = ref_grid_mapping
    ncOut.variables[varNameVP].setncatts(attDict)
    #varin = numpy.zeros((len(yin),len(xin)),dtype=vardataType)
    #varout = numpy.zeros((len(yout),len(xout)),dtype=vardataType)

    for tk in range(timeLen):
        varinP = ncInP.variables[varNameP][tk,:,:]
        varinq = ncInq.variables[varNameq][tk,:,:]

        VP = varinP*varinq /(varinq+0.622)
        ncOut.variables[varNameVP][tk,:,:] = VP[:,:]

    ncInP.close()
    ncInq.close()
    ncOut.close()
    #delete temp netcdf file
    #os.remove(temp_netcdf)


#This uses P = P0*((T0 + lambda*dz)/T0)^(-g/lambda*R)
def adjust_for_elevation_SurfaceAirPressure(input_netcdfP,  varNameP, input_netcdfT,  varNameT, output_netcdf, input_raster, target_raster, time_var_name='time', baseDateTime='1980/01/01 0', timeUnits='days'):
    """
    ajusts temerature values with lapse rate based on elevaton
    :param input_netcdf: input
    :param output_netcdf: output
    :param varName: name of variable of interest
    :param input_raster: raster containing the DEM values used when generating the original climate data (in input_netcdf)
    :param target_raster: temperature values adjusted to this DEM
    :return:
    """
    # lapse rates, monthly   ToDo: check values
    lapserateT = [0.0044,0.0059,0.0071,0.0078,0.0081,0.0082,0.0081,0.0081,0.0077,0.0068,0.0055,0.0047]
    grav = 9.81        # m/s^2
    Rgas = 287.04      # J/kg-k

    #copy variables and attriutes to output netcdf
    cmdString = "nccopy -k 3 "+input_netcdfP+" "+output_netcdf
    callSubprocess(cmdString, 'copy netcdf')

    #create elevation difference dem raster
    dzDEM = 'dzDEM_temp.tif'
    cmdString =  " python \"C:/Python34/Lib/site-packages/osgeo/scripts/gdal_calc.py\" -A "+input_raster+" -B "+target_raster +" --outfile="+dzDEM+" --calc=B-A"
    callSubprocess(cmdString, "compute elevation difference")
    #get value of dzDEM array
    ref_data = gdal.Open(dzDEM, GA_ReadOnly)
    inband = ref_data.GetRasterBand(1)
    #nodata = inband.GetNoDataValue()
    dzDEMArray = inband.ReadAsArray()
    dzDEMArrayout = dzDEMArray[::-1]         # 8.12.15 gdal and netcdf4 have arrays in opposite directions
    #dType = inband.DataType
    ref_data = None

    #Open file to write re-gridded data
    ncOut = netCDF4.Dataset(output_netcdf,"r+") #, format='NETCDF4')
    xout = ncOut.variables['x'][:]
    yout = ncOut.variables['y'][:]
    vardataType = ncOut.variables[varNameP].datatype
    varin = numpy.zeros((len(yout),len(xout)),dtype=vardataType)
    varout = numpy.zeros((len(yout),len(xout)),dtype=vardataType)
    timeLen = len(ncOut.dimensions['time'])
    tArray = ncOut.variables[time_var_name][:]
    yLen = len(yout)
    xLen = len(xout)

    ncInT = netCDF4.Dataset(input_netcdfT,"r") # format='NETCDF4')
    for tk in range(timeLen):
        varin[:,:] = ncOut.variables[varNameP][tk,:,:]
        varinT = ncInT.variables[varNameT][tk,:,:]
        if timeUnits == 'days':
            dateTimetk = datetime.strptime(baseDateTime,"%Y/%m/%d %H") + timedelta(days=int(tArray[tk]))
        elif timeUnits == 'hours':
            dateTimetk = datetime.strptime(baseDateTime,"%Y/%m/%d %H") + timedelta(hours=int(tArray[tk]))
        else:
            print("Erroneous time unit")
            return
        monthtk = int(dateTimetk.month - 1)
        #for yi in range(yLen):
        #    for xj in range(xLen):
        #       varout[yi,xj] = varin[yi,xj] - (dzDEMArrayout[yi,xj] * lapserateT[monthtk])  # bilinear_interpolation(yout[yi],xout[xj],points)
        varout = varin* math.pow((((varinT+273.16) + dzDEMArrayout*-1*lapserateT[monthtk])/(varinT+273.16)),(-1*grav/(Rgas*-1*lapserateT[monthtk])))     #T in oC and lapserrate +ve
        ncOut.variables[varNameP][tk,:,:] = varout[:,:]

    ncOut.close()


def adjust_for_elevation_VaporPressure(input_netcdf, output_netcdf, varName, input_raster, target_raster, time_var_name='time', baseDateTime='2008/10/01 0', timeUnits='hours'):
    """
    adjusts values with lapse rate based on elevaton
    :param input_netcdf: input
    :param output_netcdf: output
    :param varName: name of variable of interest
    :param input_raster: raster containing the DEM values used when generating the original climate data (in input_netcdf)
    :param target_raster: temperature values adjusted to this DEM
    :return:
    """
    # lapse rates, monthly ToDo: check values
    lapserateTd = [0.00041,0.00042,0.00040,0.00039,0.00038,0.00036,0.00033,0.00033,0.00036,0.00037,0.00040,0.00040]
    a=611.21
    b=22.452
    c=240.97

    #copy variables and attriutes to output netcdf
    cmdString = "ncea -4 -h -O "+input_netcdf+" "+output_netcdf
    callSubprocess(cmdString, 'copy netcdf')

    #create elevation difference dem raster
    dzDEM = 'dzDEM_temp.tif'
    cmdString =  " python \"C:/Python34/Lib/site-packages/osgeo/scripts/gdal_calc.py\" -A "+input_raster+" -B "+target_raster +" --outfile="+dzDEM+" --calc=B-A"
    callSubprocess(cmdString, "comute elevation difference")
    #get value of dzDEM array
    ref_data = gdal.Open(dzDEM, GA_ReadOnly)
    inband = ref_data.GetRasterBand(1)
    #nodata = inband.GetNoDataValue()
    dzDEMArray = inband.ReadAsArray()
    dzDEMArrayout = dzDEMArray[::-1]         # 8.12.15 gdal and netcdf4 have arrays in opposite directions
    #dType = inband.DataType
    ref_data = None

    #Open file to write re-gridded data
    ncOut = netCDF4.Dataset(output_netcdf,"r+") #, format='NETCDF4')
    xout = ncOut.variables['x'][:]
    yout = ncOut.variables['y'][:]
    vardataType = ncOut.variables[varName].datatype
    varin = numpy.zeros((len(yout),len(xout)),dtype=vardataType)
    varout = numpy.zeros((len(yout),len(xout)),dtype=vardataType)
    timeLen = len(ncOut.dimensions['time'])
    tArray = ncOut.variables[time_var_name][:]
    yLen = len(yout)
    xLen = len(xout)
    counterErr = 0
    for tk in range(timeLen):
        varin[:,:] = ncOut.variables[varName][tk,:,:]
        if timeUnits == 'days':
            dateTimetk = datetime.strptime(baseDateTime,"%Y/%m/%d %H") + timedelta(days=int(tArray[tk]))
        elif timeUnits == 'hours':
            dateTimetk = datetime.strptime(baseDateTime,"%Y/%m/%d %H") + timedelta(hours=int(tArray[tk]))
        else:
            print("Erroneous time unit")
            return
        monthtk = int(dateTimetk.month - 1)
        for yi in range(yLen):
            for xj in range(xLen):
                if varin[yi,xj] <= 0:
                    Tdold = -40                       # assume the lowest possible temp -40 # ToDO: may be it is better to use missing-value?
                    counterErr =+1
                    print(" 0 or negative VP at time step %d ,cell %d %d  ,value = %f"%(tk, yi, xj, varin[yi,xj]))
                else:
                    Tdold = (c * math.log(varin[yi,xj]/a)) / ( b - math.log(varin[yi,xj]/a) )
                Tdnew  = Tdold - ((c * dzDEMArrayout[yi,xj] * lapserateTd[monthtk]) / b )
                varout[yi,xj] = a * math.exp((b * Tdnew)/(c + Tdnew))
        ncOut.variables[varName][tk,:,:] = varout[:,:]

    ncOut.close()
    print("0 or negative VP occured %d times" %counterErr)


def adjust_for_elevation_Precipitation(input_netcdf, output_netcdf, varName, input_raster, target_raster, time_var_name='time', baseDateTime='2008/10/01 0', timeUnits='hours'):
    """
    ajusts temerature values with lapse rate based on elevaton
    :param input_netcdf: input
    :param output_netcdf: output
    :param varName: name of variable of interest
    :param input_raster: raster containing the DEM values used when generating the original climate data (in input_netcdf)
    :param target_raster: temperature values adjusted to this DEM
    :return:
    """
    # precipitaiton adjustment factor, monthly, from Listen and Eleder 2006   ToDo: check values
    Xp = [0.00035,0.00035,0.00035,0.0003,0.00025,0.0002,0.0002,0.0002,0.0002,0.00025,0.0003,0.00035]

    #copy variables and attriutes to output netcdf
    cmdString = "nccopy -k 3 "+input_netcdf+" "+output_netcdf
    callSubprocess(cmdString, 'copy netcdf')

    #create elevation difference dem raster
    dzDEM = 'dzDEM_temp.tif'
    cmdString =  " python \"C:/Python34/Lib/site-packages/osgeo/scripts/gdal_calc.py\" -A "+input_raster+" -B "+target_raster +" --outfile="+dzDEM+" --calc=B-A"
    callSubprocess(cmdString, "comute elevation difference")
    #get value of dzDEM array
    ref_data = gdal.Open(dzDEM, GA_ReadOnly)
    inband = ref_data.GetRasterBand(1)
    #nodata = inband.GetNoDataValue()
    dzDEMArray = inband.ReadAsArray()
    dzDEMArrayout = dzDEMArray[::-1]         # 8.12.15 gdal and netcdf4 have arrays in opposite directions
    #dType = inband.DataType
    ref_data = None

    #Open file to write re-gridded data
    ncOut = netCDF4.Dataset(output_netcdf,"r+") #, format='NETCDF4')
    xout = ncOut.variables['x'][:]
    yout = ncOut.variables['y'][:]
    vardataType = ncOut.variables[varName].datatype
    varin = numpy.zeros((len(yout),len(xout)),dtype=vardataType)
    varout = numpy.zeros((len(yout),len(xout)),dtype=vardataType)
    timeLen = len(ncOut.dimensions['time'])
    tArray = ncOut.variables[time_var_name][:]
    yLen = len(yout)
    xLen = len(xout)
    for tk in range(timeLen):
        varin[:,:] = ncOut.variables[varName][tk,:,:]
        if timeUnits == 'days':
            dateTimetk = datetime.strptime(baseDateTime,"%Y/%m/%d %H") + timedelta(days=int(tArray[tk]))
        elif timeUnits == 'hours':
            dateTimetk = datetime.strptime(baseDateTime,"%Y/%m/%d %H") + timedelta(hours=int(tArray[tk]))
        else:
            print("Erroneous time unit")
            return
        monthtk = int(dateTimetk.month - 1)

        varout = varin*(1 + Xp[monthtk]*dzDEMArrayout)/(1 - Xp[monthtk]*dzDEMArray)  # bilinear_interpolation(yout[yi],xout[xj],points)
        ncOut.variables[varName][tk,:,:] = varout[:,:]

    ncOut.close()


def adjust_for_elevation_Temperature(input_netcdf, output_netcdf, varName, input_raster, target_raster, time_var_name='time', baseDateTime='1980/01/01 0', timeUnits='days'):
    """
    ajusts temerature values with lapse rate based on elevaton
    :param input_netcdf: input
    :param output_netcdf: output
    :param varName: name of variable of interest
    :param input_raster: raster containing the DEM values used when generating the original climate data (in input_netcdf)
    :param target_raster: temperature values adjusted to this DEM
    :return:
    """
    # lapse rates, monthly   ToDo: check values
    lapserateT = [0.0044,0.0059,0.0071,0.0078,0.0081,0.0082,0.0081,0.0081,0.0077,0.0068,0.0055,0.0047]

    #copy variables and attriutes to output netcdf
    cmdString = "nccopy -k 3 "+input_netcdf+" "+output_netcdf
    callSubprocess(cmdString, 'copy netcdf')

    #create elevation difference dem raster
    dzDEM = 'dzDEM_temp.tif'
    cmdString =  " python \"C:/Python34/Lib/site-packages/osgeo/scripts/gdal_calc.py\" -A "+input_raster+" -B "+target_raster +" --outfile="+dzDEM+" --calc=B-A"
    callSubprocess(cmdString, "comute elevation difference")
    #get value of dzDEM array
    ref_data = gdal.Open(dzDEM, GA_ReadOnly)
    inband = ref_data.GetRasterBand(1)
    #nodata = inband.GetNoDataValue()
    dzDEMArray = inband.ReadAsArray()
    dzDEMArrayout = dzDEMArray[::-1]         # 8.12.15 gdal and netcdf4 have arrays in opposite directions
    #dType = inband.DataType
    ref_data = None

    #Open file to write re-gridded data
    ncOut = netCDF4.Dataset(output_netcdf,"r+") #, format='NETCDF4')
    xout = ncOut.variables['x'][:]
    yout = ncOut.variables['y'][:]
    vardataType = ncOut.variables[varName].datatype
    varin = numpy.zeros((len(yout),len(xout)),dtype=vardataType)
    varout = numpy.zeros((len(yout),len(xout)),dtype=vardataType)
    timeLen = len(ncOut.dimensions['time'])
    tArray = ncOut.variables[time_var_name][:]
    yLen = len(yout)
    xLen = len(xout)
    for tk in range(timeLen):
        varin[:,:] = ncOut.variables[varName][tk,:,:]
        if timeUnits == 'days':
            dateTimetk = datetime.strptime(baseDateTime,"%Y/%m/%d %H") + timedelta(days=int(tArray[tk]))
        elif timeUnits == 'hours':
            dateTimetk = datetime.strptime(baseDateTime,"%Y/%m/%d %H") + timedelta(hours=int(tArray[tk]))
        else:
            print("Erroneous time unit")
            return
        monthtk = int(dateTimetk.month - 1)
        #for yi in range(yLen):
        #    for xj in range(xLen):
        #       varout[yi,xj] = varin[yi,xj] - (dzDEMArrayout[yi,xj] * lapserateT[monthtk])  # bilinear_interpolation(yout[yi],xout[xj],points)
        varout = varin - (dzDEMArrayout* lapserateT[monthtk])  # bilinear_interpolation(yout[yi],xout[xj],points)
        ncOut.variables[varName][tk,:,:] = varout[:,:]

    ncOut.close()



def callSubprocess(cmdString, debugString):
    cmdargs = shlex.split(cmdString)
    #print(cmdString)
    #print(cmdargs)
    debFile = open('debug_file.txt', 'w')
    debFile.write('Starting %s \n' % debugString)
    retValue = subprocess.call(cmdString,stdout=debFile,shell=True)               # use shell=True with a single string of commands; shell=False with list of strings with first element the executable
    if (retValue==0):
        debFile.write('%s Successful\n' % debugString)
        debFile.close()
    else:
        debFile.write('There was error in %s\n' % debugString)
        debFile.close()


#create_netCDF_from_multple_nc('mm','sub','logan','loganTP.nc',-112.0,42.3,-111.0,41.6,'2008/10/01 0','2009/10/01 0',3,'time')