# join metapopulation persistence tables to patch centroids


import arcpy
import pandas as pd
import numpy as np


# MPA centroids
# Use the ones that are the ends of the connectivity lines
mpa_centroids = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\cluster_results\scripts\COMBINED.gdb\patch_centroids'

# mpa metapopulation persistence table
dir = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\metapop_pers\csv'
timesteps = 1750
metapop = f'metapop_pers_{timesteps}_'

# working gdb
gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\metapop_pers\metapop_pers.gdb'
arcpy.env.workspace = gdb

# mpas to exclude that are deep in inlets or narrow areas
mpas_exclude = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas.gdb\M10_toexcludefromanalysis'



# make list of files that match metapop
files = os.listdir(dir)
fs = []
for file in files:
    if file.startswith(metapop):
        fs.append(file)

# if I want to exclude any, do it manually:
fs = fs[:]

# create df equal to first file in the list
df = pd.read_csv(os.path.join(dir, fs[0]))
df = df.drop(['Unnamed: 0'], axis=1)
# for each file, read in as df
for file in fs[1:]:
    df_metapop = pd.read_csv(os.path.join(dir, file))
    df_metapop = df_metapop.drop(['Unnamed: 0', 'popn'], axis=1)

    # check if there are any duplicates and remove
    col_exist = list(df.columns)
    col_new = list(df_metapop.columns)[1:] # first one will be uID
    matches = list(set(col_exist).intersection(col_new))
    if len(matches) > 0:
        for col in matches:
            df_metapop = df_metapop.drop([col], axis=1)

    # merge
    df = df.merge(df_metapop, on='uID_20201124')



# output to arc
x = np.array(np.rec.fromrecords(df.values))
names = df.dtypes.index.tolist()
x.dtype.names = tuple(names)
out_string = f'metapop_pers_{timesteps}'
arcpy.da.NumPyArrayToTable(x, os.path.join(gdb, out_string))


# join and copy
arcpy.env.qualifiedFieldNames = False
joined_table = arcpy.AddJoin_management(mpa_centroids, 'uID_202011', out_string, 'uID_20201124')

# remove mpas to be excluded
to_include = []
with arcpy.da.SearchCursor(mpas_exclude, ['uID_20201124', 'exclude']) as cursor:
    for row in cursor:
        if not row[1] == 1:
            to_include.append(row[0])
mpas_sel = arcpy.SelectLayerByAttribute_management(joined_table, 'NEW_SELECTION', '"uID_202011" IN {}'.format(str(tuple(to_include))))


arcpy.CopyFeatures_management(mpas_sel, f'metapop_pers_centroids')
