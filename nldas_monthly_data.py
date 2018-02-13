"""
This is to use the hourly NLDAS data to be combined as monthly data for input preparation

"""


try:
    from osgeo import gdal, osr, ogr
except:
    import gdal, osr, ogr
import shlex
import subprocess
from datetime import datetime


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


def group_monthly_netCDF_from_multple_nc_NLDAS(file_prefix, startDateTime, endDateTime,
                                            out_timeName = 'time'):
    """
    Subsets and combines multiple netcdf files
    for nldas forcing
    """
    startYear = datetime.strptime(startDateTime,"%Y/%m/%d %H").year
    endYear = datetime.strptime(endDateTime,"%Y/%m/%d %H").year
    startMonth = datetime.strptime(startDateTime,"%Y/%m/%d %H").month
    endMonth = datetime.strptime(endDateTime,"%Y/%m/%d %H").month
    startDay =  datetime.strptime(startDateTime,"%Y/%m/%d %H").timetuple().tm_yday        #start date = day of year for 2010
    endDay   =  startDay + (datetime.strptime(endDateTime,"%Y/%m/%d %H") - datetime.strptime(startDateTime,"%Y/%m/%d %H")).days          # end date = day of year for 2011 + 365

    print(startYear)
    print(endYear)
    for year in range(startYear, endYear+1):
        for month in range(1, 13):
            if month < 10:
                monthS = '0'+str(month)
            else:
                monthS = str(month)
            cmdString = "ncecat -4 -H -h -O -u "+out_timeName+"  "+file_prefix+"*"+str(year)+monthS+"*.nc  -o "+file_prefix+"_Monthly_"+str(year)+monthS+".nc"                      #-H don't append input file list -h don't append history
            #cmdString = "dir /b "+file_prefix+"*"+str(year)+monthS+"*.nc | ncecat -u \""+out_timeName+"\" -H -O -o "+file_prefix+"_Monthly_"+str(year)+monthS+".nc"
            try:
                callSubprocess(cmdString, "concatenate netcdf files into monthly bundle")
                'Done for {}/{}'.format(year, monthS)
            except Exception as e:
                print 'Failed for {}/{}'.format(year, monthS)



group_monthly_netCDF_from_multple_nc_NLDAS(file_prefix='NLDAS_FORA0125_H.A',startDateTime='2016/01/01 0', endDateTime='2017/12/31 0', out_timeName = 'time')

print 'Done'