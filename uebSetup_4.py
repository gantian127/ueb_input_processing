import os
import netcdfFunctions, rdhmFunctions_ori
"""*********** Input settings for watershed of interest *****************"""
workingDir = "/Projects/Tian_workspace/SAC_Animas_1989_watershed/SAC_input"
mSWIT = 'SWIT.nc'
mmSWIT = 'SWITmm.nc'
startDateTime='1989/10/01 0'
time_varName='time'
output_prefix = 'xmrg'
input_varName = 'SWIT'
new_varUnit = 'mm'
mult_Factor = 1000
offSet = 0.0
time_length = 2904  #11672    #17512
dThour = 3  #this is determined by the swit.nc file time interval
proj4_string = '+proj=stere +lat_0=90.0 +lat_ts=60.0 +lon_0=-105.0 +k=1 +x_0=0.0 +y_0=0.0 +a=6371200 +b=6371200 +units=m +no_defs'  #proj4_string: see the paper by David Maidment etal
#
#proj4_string = '+proj=stere +lat_0=90.0 +lat_ts=37.60 +lon_0=-105 +k_0=1 +x_0=0 +y_0=0 +ellps=sphere +a=6371200 +b=6371200  +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'  #for Animas to reduce area distortion
#standard parallel use 60.00681402 if radius is not 6371200
#proj4_string = '+proj=stere +lat_0=90.0 +lat_ts=60.00681402 +lon_0=-105 +k_0=1 +x_0=401.0 +y_0=1601.0 +ellps=sphere +a=6371200 +b=6371200  +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'
#"""
print(os.getcwd())
os.chdir(workingDir)
print('start units conversion')
netcdfFunctions.convert_netcdf_Units(mSWIT, mmSWIT, input_varName, new_varUnits=new_varUnit, multiplier_Factor=mult_Factor, offset=offSet)
print('start asci file conversion')
rdhmFunctions_ori.create_multiple_Ascii_fromNetCDF(mmSWIT, output_prefix, input_varName, time_length, dThour, startDateTime=startDateTime, time_varName=time_varName, proj4_string=proj4_string)

print("netcdf to ascii file done")

