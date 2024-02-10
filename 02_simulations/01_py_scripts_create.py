# create python scripts for running opendrift simualations and biology module

# Create a folder for each time period
# Create an Opendrift simulation for each grouping of mpas within each period (9*9)
# Create a biology script for each time period (1*9)

import os
from shutil import copyfile


scripts_dir = 'scripts'
folders = [
    1101, 1105, 1108,
    1401, 1405, 1408,
    1701, 1705, 1708]
folders = [f'sim{str(i)}' for i in folders]
groups = list(range(1,10))
times = [
    '20110101_20110316', '20110501_20110714', '20110801_20111014',
    '20140101_20140316', '20140501_20140714', '20140801_20141014',
    '20170101_20170316', '20170501_20170714', '20170801_20171014',]


# Opendrift scripts
with open('templates/mpa_20210205_cluster.py', 'r') as file:
    od_script = file.read()
for folder, time in zip(folders, times):
    dir = os.path.join(scripts_dir, folder)
    if not os.path.exists(dir):
        os.makedirs(dir)
    for group in groups:
        name = os.path.join(dir, f'od{str(group)}.py')
        with open(name, 'w') as pys:
            pys.write(f'''dates = '{time}'\n''')
            pys.write(f'shp_group = {group}\n')
            pys.write(od_script)


# Biology scripts
with open('templates/biology.py', 'r') as file:
    bi_script = file.read()
for folder in folders:
    dir = os.path.join(scripts_dir, folder)
    if not os.path.exists(dir):
        os.makedirs(dir)
    for group in groups:
        name = os.path.join(dir, f'bi{str(group)}.py')
        with open(name, 'w') as pys:
            pys.write(f'nc_group = {group}\n')
            pys.write(bi_script)    