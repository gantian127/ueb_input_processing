"""
This is to move the original prepared files for 2010,1988, 1989-2009 to a climate forcing folder
"""

import shutil
import os
import re


input_dir_list = [
                    r'D:\3_NASA_Project\Model output and plots\22yr_Animas_watershed\tian_resutls_prcp_temp_wind_vp1988',
                    r'D:\3_NASA_Project\Model output and plots\22yr_Animas_watershed\tian_results_prcp_temp_wind_vp2010',
                    r'D:\3_NASA_Project\Model output and plots\22yr_Animas_watershed\tian_results_prcp_temp_1989_2009',
                    r'D:\3_NASA_Project\Model output and plots\22yr_Animas_watershed\tian_results_vp_wind_1989_2009'
                 ]

add_factor_list = [
                   0,
                   21,
                   1,
                   1,
                ]
output_dir = r'D:\3_NASA_Project\Model output and plots\22yr_Animas_watershed\climate_forcing'


for input_dir, add_factor in zip(input_dir_list, add_factor_list):
    print input_dir
    print os.listdir(input_dir)
    for file_name in os.listdir(input_dir):
        src = os.path.join(input_dir, file_name)
        number = re.search(r'\d+', file_name[-5:]).group()
        new_name = file_name[:-5]+file_name[-5:].replace(number, str(int(number)+add_factor))
        print number, new_name
        des = os.path.join(output_dir, new_name)
        shutil.copy(src, des)

print 'Done!'