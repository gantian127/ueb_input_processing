"""
This is to update the xmrg climate forcing file name to the right UTC time

# 5min/wy for uebVp, uebWindS, uebTamin, uebTamax
# 43min/22yr for uebVp, uebWindS, uebTamin, uebTamax

"""

import os
from datetime import datetime, timedelta
import csv
from shutil import copy2

# step1 user settings ###################################################
ori_forcing_folder = r'/Projects/Tian_workspace/rdhm_ueb_modeling/McPhee_MPHC2/MPHC2_forcing_validation/Forcing'
update_forcing_folder = r'/Projects/Tian_workspace/rdhm_ueb_modeling/McPhee_MPHC2/MPHC2_forcing_validation/Final_forcing'

if not os.path.isdir(update_forcing_folder):
    os.mkdir(update_forcing_folder)

shift_info = {
    'uebTair': [12, -1],  # hour, -1 backward, 1 forward
    'uebPrec': [12, -1],
    'uebTamin': [0, 0],
    'uebTamax': [24, -1],
    'uebVp':  [0, 0],
    'uebWindS': [0, 0],
}

# step2 update file name and copy file to new folder ####################################
for var, shift_setting in shift_info.items():
    print 'start {}'.format(var)
    shift_hr, shift_dire = shift_setting
    ori_file_list = [file_name for file_name in os.listdir(ori_forcing_folder) if var in file_name and 'z.gz' in file_name]

    record_list = [['ori file', 'update file']]
    if ori_file_list:
        for ori_file in ori_file_list:
            try:
                # update file name
                ori_time_str = ori_file.replace(var, '').replace('z.gz', '')
                ori_time = datetime.strptime(ori_time_str, '%m%d%Y%H')
                update_time = ori_time + shift_dire * timedelta(hours=shift_hr)  # remember to specify as hours =
                update_time_str = update_time.strftime('%m%d%Y%H')
                update_file = var+update_time_str+'z.gz'
                record_list.append([ori_file, update_file])

                # copy file to new folder
                ori_path = os.path.join(ori_forcing_folder, ori_file)
                new_path = os.path.join(update_forcing_folder, update_file)
                copy2(ori_path, new_path)
            except Exception as e:
                print 'failed for {}'.format(ori_file)
                continue

        with open('{}.csv'.format(var), 'w') as var_file:
            writer = csv.writer(var_file)
            writer.writerows(record_list)

print 'time shift is done'

# fix missing tair
#  find . -name 'uebTair10122010*'
#  for i in uebTair10122010*; do mv $i ori_$i ; done
#  find . -name 'uebTair10122010*'
#  find . -name 'ori_uebTair10122010*'
#  for i in uebTair10112010*; do cp $i ${i:0:9}12${i:11:10}; done
#  find . -name 'uebTair10122010*'

# fix missing tmax tmin
# find . -name 'uebTamax1012201000z.gz'
# mv uebTamax1012201000z.gz ori_uebTamax1012201000z.gz
# find . -name 'uebTamax1012201000z.gz'
# cp uebTamax1011201000z.gz uebTamax1012201000z.gz
# find . -name 'uebTamax1012201000z.gz'

# find . -name 'uebTamin1012201000z.gz'
# mv uebTamin1012201000z.gz ori_uebTamin1012201000z.gz
# find . -name 'uebTamin1012201000z.gz'
# cp uebTamin1011201000z.gz uebTamin1012201000z.gz
# find . -name 'uebTamin1012201000z.gz'
