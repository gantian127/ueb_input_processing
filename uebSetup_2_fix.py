import os
import climateFunctions_2_ori
import netcdfFunctions_ori
"""*********** Input settings for watershed of interest *****************"""
watershedN = 'Green'
startYear=1988
endYear=2009
## Domain bounding box in geographic coordinates left, top, right, bottom.  Must enclose watershed of interest
# Upper Colorado River WS above Lake Powel
#leftX, topY, rightX, bottomY =  -113.860, 43.493, -104.420, 35.950
# Animas River WS above Durango
#leftX, topY, rightX, bottomY =  -108.0505, 37.9695, -107.5150, 37.2626

#4.20 for 22 years Animas simulation
#leftX, topY, rightX, bottomY =  -108.15, 38.06, -107.41, 37.16

# rect Domain: Green River near Daniel at Warren Bridge
leftX, topY, rightX, bottomY =  -110.415, 43.593, -109.492, 42.871

# Little Bear River WS upstream of Cutler reservoir, near Logan UT
# leftX, topY, rightX, bottomY =  -112.05, 41.74, -111.49, 41.36
#Logan River above First Dam
#leftX, topY, rightX, bottomY = -111.8037, 42.1248, -111.4255, 41.6946
# Little Bear River Watershed upstream of Hyrum Reservoir
#leftX, topY, rightX, bottomY = -111.97, 41.629, -111.48, 41.36
# Weber River Watershed upstream of GSL
#leftX, topY, rightX, bottomY = -112.2378, 41.4803, -110.8147, 40.5478
# San Juan River Above Farmington, NW
# leftX, topY, rightX, bottomY =
# Dolores-Dolores, upstream of Mcphee reservoir
#leftX, topY, rightX, bottomY =  -108.71, 37.88, -107.81, 37.41
# Blue river below Green Mountain Reservoir
#leftX, topY, rightX, bottomY =  -106.419, 39.931, -105.772, 39.349
# Blue river below Dillon
#leftX, topY, rightX, bottomY =  -106.282, 39.685, -105.772, 39.349
# Williams Fork below Williams Fork Reservoir
#leftX, topY, rightX, bottomY =  -106.301, 40.045, -105.870, 39.677
# Frazer river at Granby
#leftX, topY, rightX, bottomY =  -106.07, 40.137, -105.663, 39.769
# Frazer river at Winter Park
#leftX, topY, rightX, bottomY =  -105.814, 39.909, -105.680, 39.79
# Green River near Daniel at Warren Bridge
#leftX, topY, rightX, bottomY =  -110.315, 43.493, -109.592, 42.971

# Uncompahgre river at Ridgeway Reservoir
#leftX, topY, rightX, bottomY =  -108.002, 38.276, -107.507, 37.867
# Red Butte Creek at Red Butte Reservoir
#leftX, topY, rightX, bottomY =  -111.819, 40.829, -111.736, 40.767


#cbrfc forcing
"""
workingDir = "/Projects/cbrfcTP/"
os.chdir(workingDir)
for cYear in range(startYear,endYear):
    startDateTime=str(cYear)+'/10/01 0'
    endDateTime=str(cYear+1)+'/07/30 0'
    cbrfcPWS = watershedN+str(cYear+1)+'P.nc'
    cbrfcTWS = watershedN+str(cYear+1)+'T.nc'
    climateFunctions_2.create_netCDF_from_multple_nc('mm', watershedN, cbrfcPWS, 'ogrid', leftX, topY, rightX, bottomY, startDateTime, endDateTime, 3, out_timeName = 'time')
    climateFunctions_2.create_netCDF_from_multple_nc('mm', watershedN, cbrfcTWS, 'otgrid', leftX, topY, rightX, bottomY, startDateTime, endDateTime, 3, out_timeName = 'time')
"""

#fill -9999 cbrfc temp
###TODO: 9.6.17 This is not enough need to ensure watershed.nc no-value cells don't get temp value (curr. getting 0 by interp)
workingDir = "/Projects/Tian_workspace"
os.chdir(workingDir)
for cYear in range(startYear,endYear):
    cbrfcTWS = watershedN+str(cYear+1)+'T.nc'
    cbrfcTWS_Fill = watershedN + str(cYear + 1) + 'Tfill.nc'
    cbrfcTWS_FillT = watershedN + str(cYear + 1) + 'TfillT.nc'
    ##spatial fill
    netcdfFunctions_ori.replace_no_data_with_SpatialAverage(cbrfcTWS, 'otgrid', cbrfcTWS_Fill, in_Xcoord='lon', in_Ycoord='lat')
    ##temporal fill
    netcdfFunctions_ori.replace_no_data_along_temporal_Dimension(cbrfcTWS_Fill, 'otgrid', cbrfcTWS_FillT, in_Xcoord='lon', in_Ycoord='lat')

#nldas2
"""
workingDir = "/Projects/nldasMonthly/"
os.chdir(workingDir)
for cYear in range(startYear,endYear):
    startDateTime=str(cYear)+'/10/01 0'
    endDateTime=str(cYear+1)+'/10/01 0'
    nldas2WS = watershedN+str(cYear+1)+'NLDAS2.nc'
    climateFunctions_2.subset_and_concatenate_netCDF_from_multple_nc_NLDAS('NLDAS_FORA0125_H.A_Monthly', watershedN, nldas2WS, leftX, topY, rightX, bottomY, startDateTime, endDateTime,
                                                                               dT=1, in_Xcoord = 'lon_110', in_Ycoord='lat_110',inout_timeName = 'time')
    #create_netCDF_from_multple_nc2('NLDAS', watershedN, nldas2WS, leftX, topY, rightX, bottomY, startDateTime, endDateTime, dT=1, in_Xcoord = 'lon_110', in_Ycoord='lat_110',out_timeName = 'time')

"""
print("done")

