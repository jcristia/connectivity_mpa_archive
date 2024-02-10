# Calculate connectivity metrics for the time-averaged connectivity lines for
# each PLD

import networkx as nx
import arcpy
import numpy as np
import pandas as pd
import os



root = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\cluster_results\scripts'
input_gdb = 'COMBINED.gdb'
output_gdb = r'CONN_METRICS.gdb'
plds = [1,3,7,10,21,30,40,60]
centroids = 'patch_centroids'


if not arcpy.Exists(os.path.join(root, output_gdb)):
    arcpy.CreateFileGDB_management(root, output_gdb)

# copy centroids to new gdb
for pld in plds:
    arcpy.env.workspace = os.path.join(root, output_gdb)
    out_pts = 'conn_metrics_pld{}'.format(str(pld))
    arcpy.CopyFeatures_management(os.path.join(root, input_gdb, centroids), out_pts)

# calc metrics for each PLD
arcpy.env.workspace = os.path.join(root, input_gdb)
fcs = arcpy.ListFeatureClasses(wild_card='conn_avg_pld*')
for fc in fcs:
    arcpy.env.workspace = os.path.join(root, input_gdb)
    pld = fc[12:]
    if int(pld) not in plds:
        continue

    arr = arcpy.da.FeatureClassToNumPyArray(fc, ('from_id', 'to_id', 'prob_avg'))
    df = pd.DataFrame(arr)

    # remove self connections
    df = df[df.from_id != df.to_id]
    # This is also where I could remove specific MPA connections if necessary.

    G = nx.from_pandas_edgelist(df, source='from_id', target='to_id', edge_attr='prob_avg', create_using=nx.DiGraph)
    
    bt = nx.betweenness_centrality(G, k=None, normalized=True, weight='prob_avg', endpoints=False, seed=None)
    dca = nx.degree_centrality(G)
    dci = nx.in_degree_centrality(G)
    dco = nx.out_degree_centrality(G)

    # add metrics as node attributes
    nx.set_node_attributes(G, bt, 'bt')
    nx.set_node_attributes(G, dca, "dca")
    nx.set_node_attributes(G, dci, "dci")
    nx.set_node_attributes(G, dco, "dco")
    

    df_att = pd.DataFrame({
        'node':list(G.nodes), 
        'bt':bt.values(),
        'dca':dca.values(),
        'dci':dci.values(),
        'dco':dco.values()
            })


    arcpy.env.workspace = os.path.join(root, output_gdb)
    pts = 'conn_metrics_pld{}'.format(str(pld))
    metrics = ['bt', 'dca', 'dci', 'dco']
    for field in metrics:
        arcpy.AddField_management(pts, field, 'DOUBLE')

    # arcpy CalculateField is giving me a weird error saying it can't find Python
    # for the sake of time, I will just do this manually with an UpdateCursor
    fields = metrics[:]  # need to be careful how I copy lists
    fields.append('uID_202011')        
    with arcpy.da.UpdateCursor(pts, fields) as cursor:
        for row in cursor:
            if row[-1] not in df_att.node.values:
                continue # I should only need this for testing
            dftemp = df_att[df_att.node==row[-1]]
            for i,m in enumerate(metrics):
                row[i] = dftemp[m].values[0]
            cursor.updateRow(row)



# stop there until I know what I want from these