
# if running on cedar, then py_scripts_create.py adds a variable above this line
# for which group nc file to use


# Script to add bio part of a biophysical model
# This script takes the output netcdf file from an Opendrift simulation and 
# modifies particle end trajectories by consdering precompetency period, 
# settlement, and mortality.

# John Cristiani
# University of British Columbia
# origin: 2019-04-19
# chapter 2 updates: 2021-01-15

# env: biology


import netCDF4 as nc
import numpy as np
from shapely.geometry import shape, Point, LineString, Polygon
import pandas as pd
import geopandas
import logging
import math
import sys
import os
logging.basicConfig(level=logging.INFO)



###################
# paths and variables
###################

# unix version of paths:
cluster_run = True
root = os.getcwd()
base = os.path.dirname(__file__)
####
input_nc_dir = os.path.join(root, base, 'outputs/nc') # where nc files were output
input_npy_dir = os.path.join(root, base, 'outputs/npy') # where npy files were output
####
shp_og = os.path.join(root, 'inputs/mpas/mpa_.shp') # full release_polys dataset
shp_og_buff = os.path.join(root, 'inputs/mpas/mpa_buff.shp')  # buffered for checking settlement
####
crs_input_shp = {'init' :'epsg:3005'} # BC Albers
output_shp_dir = os.path.join(root, base, 'outputs/shp') # connection lines

# # Local version of paths:
# cluster_run = False
# ####
# input_nc_dir = r'outputs/nc' # where nc files were output
# input_npy_dir = r'outputs/npy' # where npy files were output
# ####
# shp_og = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas_shp_release\mpa_.shp' # full release_polys dataset
# shp_og_buff = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas_shp_release\mpa_buff.shp'  # buffered for checking settlement
# ####
# crs_input_shp = {'init' :'epsg:3005'} # BC Albers
# output_shp_dir = r'outputs/shp' # connection lines

precomp = 0
# if I am using 'stranding' in opendrift, then I likely need at least a small 
# precompetency period. Otherwise everything just ends up settling at home. 
# Be careful setting this. It depends on the time_step and time_step_output 
# you used in the run. It is in units of the timestep output. If time step 
# output is 30 minutes, then a precomp of 2 is 1 hour.

# get these values from the simulation script
time_step_output = 0.5 # in hours. It will be in seconds in the opendrift script
interval_of_release = 4 # in hours, interval can't be less than time step output
# if no delayed release then just put same value as time_step_output
num_of_releases = 84 # if no delayed release then just put 1

# allow particles to settle on a poly?
settlement_apply = True

# mortality
mortality_rate = 0.15 # instantaneous daily rate
mort_period = 8 # after how many time_step_outputs to apply mortality rate 
# MAKE THIS A FACTOR OF 24
# The mortality rate will be scaled appropriately so that it still matches the 
# daily rate. This option is given because it seems uncessary to apply it at 
# every time step, but for some species with short PLDs, it will make sense to 
# apply it more often than once per day. If mortality rate is 0, the also set 
# this to 0.

# if it is an Opendrift backwards run
# !!! I noticed some inconsistent behavior with the shape of resulting arrays
# for each particle in the time dimension. If they are not all the same length
# this would create an issue for how I analyze backwards runs. I will need to
# look into this more if I need to do backwards runs.
# So, for now this feature is considered in development. 
backwards_run = False

# PLD limit
# I run opendrift simulations for what I expect the max PLD to be that I am 
# considering. Then in this script I can set a smaller PLD and see which 
# connections are made if I had only run the simulation up to a certain 
# timestep. I need to do PLDs all at once on each run of this script because 
# mortality is random and I want all PLDs done on one random selection of 
# particles instead of on different selections.
# Provide PLDs in a list in units of days
plds = [1, 3, 7, 10, 21, 30, 40, 60]
pld_index = None
if cluster_run:
    pld_index = int(sys.argv[1])
    plds = plds[pld_index:pld_index+1]

# !!! Reduce the range of nc files that get processed?
# see for loop at bottom

###################
# BIOLOGY FUNCTIONS
###################

###################
# assign a polygon unique ID to each particle
###################

def get_particle_originPoly(
    origin_marker, traj, crs_input_shp, lat_np, lon_np, backwards_run):

    logging.info("Getting origin coordinates of each particle")

    # get starting lat/lon for each particle
    lons = np.load(lon_np)
    lats = np.load(lat_np)
    # reverse order for backwards runs (opendrift records them in reverse in 
    # the netcdf file)
    if backwards_run:
        lons = lons[::-1]
        lats = lats[::-1]
    
    # get origin_marker (uID)
    om = []
    for i in origin_marker:
        o = i[0]
        if np.ma.is_masked(o):
            o = i[i.mask == False][0]
        om.append(o)

    # add origin coords , trajid and uID to df
    df = pd.DataFrame()
    df['o_coords'] = list(zip(lons, lats))
    df['o_coords'] = df['o_coords'].apply(Point)
    df['traj_id'] = list(traj)
    df['uID_part'] = list(om)
    points = geopandas.GeoDataFrame(df, geometry='o_coords')
    points.crs = {'init' :'epsg:4326'}
    points = points.to_crs(crs_input_shp)
    origin = pd.DataFrame(data=points)
    
    return origin


###################
# determine precompetency and release intervals
###################

def calc_precomp(
    precomp, time_step_output, particles_per_release, interval_of_release, 
    num_of_releases, traj):

    # the timesteps when we release particles
    timesteps_with_release = []
    for release in range(num_of_releases):
        ts = (float(interval_of_release) / float(time_step_output)) * release
        timesteps_with_release.append(int(ts))

    # when the precompetency period ends for each group of particle releases
    precomp_end_timestep = []
    for release in timesteps_with_release:
        ts_e = release + precomp
        precomp_end_timestep.append(ts_e)

    # the range of time periods of the precomp period for each group of particles
    precomp_range = []
    for p in precomp_end_timestep:
        precomp_range.append([p-precomp, p])

    # the corresponding particle IDs for each release
    particle_range = []
    if num_of_releases == 1:
        particle_range = [[1, len(traj) + 1]]
    else:
        for release in range(1,num_of_releases+1):
            # Opendrift keeps particles in order that they are released. Hopefully this never changes.
            p_range = [1 + ((release-1) * particles_per_release),(release * particles_per_release) +1]
            particle_range.append(p_range)

    return timesteps_with_release, precomp_end_timestep, precomp_range, particle_range


###################
# settlement
# associate a particle with a poly if it strands on the coast on it
# or if it is over it at the end of a PLD period
###################

def settlement(
    settlement_apply, origin, shp_og_buff, timestep, status, lon, lat, traj, 
    crs_input_shp, precomp, precomp_range, particle_range, pld_int,
    timesteps_with_release):

    poly  = geopandas.GeoDataFrame.from_file(shp_og_buff)
    poly.crs = crs_input_shp
    dest_df = pd.DataFrame(columns=['d_coords','traj_id','dest_id','time_int'])

    if settlement_apply: # if this is false, then it will just join the blank 
        # dest_df to origin, and the get_destination_coords function will fill 
        # in the rest. This would be if you just want the points shapefile.

        for ts_i in range(0,len(timesteps_with_release)):
            # I am applying active settlement at the end of each pld period,
            # However, since some particles are on a delay release, I need to account
            # for the pld plus the additional timesteps that they were delayed

            ts = timesteps_with_release[ts_i]
            particles = list(range(particle_range[ts_i][0],particle_range[ts_i][1]))

            logging.info(
                "Processing settlement for particle_range {} out of {}".format(particle_range[ts_i], particle_range[-1][-1]))

            # I'm now only going to check the timestep at the end of the PLD
            # When I am doing it for every PLD, it was too slow.
            # If a future project needs active settlement at each time step then
            # I will just write a separate biology script for it. It is too
            # difficult to maintain this script for multiple uses.
            i = ts+pld_int
            if i >= len(status[0]): # setting this as >= instead of > because I'm using i as an index, so it has to be 1 less than the length
                # check if the pld is greater than the steps written in the nc
                # file. If every particle settles before a certain PLD, then
                # it will be shorter.
                # In that case, just use the last timestep available.
                i = len(status[0])-1
            
            # remove particles that are not part of the current particle range
            mask = np.isin(traj, particles)
            t = traj[mask]       

            # since some particles may be masked at this timestep, find last
            # timestep where a particle is not masked
            lons_dest = []
            lats_dest = []
            time_steps = []
            for particle in t:
                index = np.where(traj[:] ==  particle)[0][0]
                if np.ma.is_masked(lon[index][i]):
                    # if the last value is masked then it must have stranded and we can 
                    # search by where it is 1.
                    # If a particle goes outside of the grid it gets coded as '2 - 
                    # missing data'. Anything above 0 is considered deactivated.
                    j = np.where(status[index] > 0)[0][0]
                    lo = lon[index][j]
                    lons_dest.append(lo)
                    la = lat[index][j]
                    lats_dest.append(la)
                    time_steps.append(j)
                else: # otherwise just use i
                    lons_dest.append(lon[index][i])
                    lats_dest.append(lat[index][i])
                    time_steps.append(i)

            df = pd.DataFrame()
            df['d_coords'] = list(zip(lons_dest, lats_dest))
            df['d_coords'] = df['d_coords'].apply(Point)
            df['traj_id'] = list(t)
            df['time_int'] = time_steps
            points = geopandas.GeoDataFrame(df, geometry='d_coords')
            points.crs = {'init' :'epsg:4326'}
            points = points.to_crs(crs_input_shp)
            # This should be 1:M join. A point can overlap many MPAs:
            pointInPolys = geopandas.tools.sjoin(points, poly, how='inner')
            pointInPolys = pointInPolys.rename(columns={'uID_202011':'dest_id'})
            dest_df = dest_df.append(pointInPolys[['d_coords','traj_id','dest_id','time_int']], ignore_index=True)
    
    # join the two tables
    logging.info("merging destination and origin dataframes")
    # Need to coerce merge. traj_id must be numeric. The dest_df data types were
    # all "object". This was not a problem on windows, but when running on the 
    # cluster it woud give an error
    dest_df = dest_df.infer_objects()
    origin = origin.infer_objects()
    dest_df.traj_id = dest_df.traj_id.astype('float')
    origin.traj_id = origin.traj_id.astype('float')
    # This could be a M:1 relationship if a point overlaps many MPAs.
    # Pandas should duplicate origin
    origin_dest = dest_df.merge(origin, on='traj_id', how='outer')

    return origin_dest

###################
#
#  Add the final destination coordinates to particles that did not settle on a patch
#
###################

def get_destination_coords(
    origin_dest, traj, lon, lat, timestep, crs_input_shp, status, particle_range,
    pld_int, timesteps_with_release):

    dest_df = pd.DataFrame(columns=['Coordinates','traj_id','time_step'])

    for ts_i in range(0,len(timesteps_with_release)):
        # Check the destination of all particles that did not settle at the end
        # of the PLD.
        # However, since some particles are on a delay release, I need to account
        # for the pld plus the additional timesteps that they were delayed

        ts = timesteps_with_release[ts_i]
        particles = list(range(particle_range[ts_i][0],particle_range[ts_i][1]))

        logging.info(
            "Getting destination coordinates for particle_range {} out of {}".format(particle_range[ts_i], particle_range[-1][-1]))

        i = ts+pld_int
        if i >= len(status[0]): # setting this as >= instead of > because I'm using i as an index, so it has to be 1 less than the length
            # check if the pld is greater than the steps written in the nc
            # file. If every particle settles before a certain PLD, then
            # it will be shorter.
            # In that case, just use the last timestep available.
            i = len(status[0])-1
        
        # remove particles that are not part of the current particle range
        mask = np.isin(traj, particles)
        t = traj[mask]

        lons_dest = []
        lats_dest = []
        time_steps = []
        for particle in t:
            index = np.where(traj[:] == particle)[0][0]
            if np.ma.is_masked(lon[index][i]):
                # if the value at that pld is masked then it must have stranded eariler
                # and we can search by where it is 1.
                # If a particle goes outside of the grid it gets coded as '2 - 
                # missing data'. Anything above 0 is considered deactivated.
                j = np.where(status[index] > 0)[0][0]
                lo = lon[index][j]
                lons_dest.append(lo)
                la = lat[index][j]
                lats_dest.append(la)
                time_steps.append(j)
            else: # otherwise just use i
                lons_dest.append(lon[index][i])
                lats_dest.append(lat[index][i])
                time_steps.append(i)

        df = pd.DataFrame()
        df['Coordinates'] = list(zip(lons_dest, lats_dest))
        df['Coordinates'] = df['Coordinates'].apply(Point)
        df['traj_id'] = list(t)
        df['time_step'] = time_steps
        points_dest = geopandas.GeoDataFrame(df, geometry='Coordinates')
        points_dest.crs = {'init' :'epsg:4326'}
        points_dest = points_dest.to_crs(crs_input_shp)
        dest_df = dest_df.append(points_dest[['Coordinates', 'traj_id', 'time_step']], ignore_index=True)

    # join, fill in values where null, remove columns
    logging.info("joining destination coordinates to dataframe")
    dest_df = dest_df.infer_objects()
    dest_df.traj_id = dest_df.traj_id.astype('float')
    origin_dest = origin_dest.merge(dest_df, on='traj_id')
    origin_dest['time_int'].loc[origin_dest['time_int'].isnull()] = origin_dest['time_step']
    origin_dest['d_coords'].loc[origin_dest['d_coords'].isnull()] = origin_dest['Coordinates']
    origin_dest = origin_dest.drop(['time_step', 'Coordinates'], axis=1)

    origin_dest = origin_dest.sort_values(by=['traj_id'])
    origin_dest = origin_dest.reset_index(drop=True)

    return origin_dest


###################
# calculate mortality
###################

def calc_mortality(
    mortality_rate, traj, timestep, origin_dest, time_step_output, mort_period, 
    timesteps_with_release, particle_range):

    logging.info("calculating mortality")

    mortality_p = pd.DataFrame(columns=['traj_id','mortstep'])

    if mortality_rate > 0:
        timestep_days = time_step_output / 24   # proportion of a day for one timestep
        mort_timesteps = np.arange(0, len(timestep)-1, mort_period)   # timesteps to apply mortality
        # instantaneous mortality rate for the mortality application interval
        inst_rate = 1 - (math.exp(math.log(1-mortality_rate) * (timestep_days * mort_period)))
        # so as long as the interval that I calculate mortality stays the same 
        # throughout a simulation, I don't need to worry about a changing rate 
        # or how many new particles enter the system

        # need to not consider particles that are not released yet 
        timesteps_with_release = np.array(timesteps_with_release)
        particle_range = np.array(particle_range)

        mortality_p = pd.DataFrame(columns=['traj_id','mortstep'])

        for i in mort_timesteps[1:]:
            mortality_selection = traj[:]
            # remove ones that have not seeded yet
            if i <= np.max(timesteps_with_release):
                # get indices of all timesteps >= current time step
                periods_exempt = particle_range[np.where(timesteps_with_release[:] >= i)]
                # get particles that should not be considered for mortality
                p_remove = np.arange(np.min(periods_exempt), np.max(periods_exempt))
                # remove from mortality_selection
                mortality_selection = np.setdiff1d(mortality_selection, p_remove)

            # remove ones that have already been added to mortality_p
            mortality_selection = np.setdiff1d(mortality_selection, mortality_p['traj_id'].values)

            # remove ones that have settled before this time step
            # select from origin_dest where time_int is this timestep and where 
            # dest_id is not null (so actually settled somewhere, I can still 
            # killed stranded ones)
            # .values turns the selection into a numpy array
            p_settled = origin_dest['traj_id'].loc[(origin_dest["time_int"] <= i) & (origin_dest['dest_id'].notnull())].values
            mortality_selection = np.setdiff1d(mortality_selection, p_settled)

            # select random particles based on inst_rate
            num_to_kill = int(len(mortality_selection) * inst_rate)
            mortality_selection = np.random.choice(mortality_selection, num_to_kill, replace=False)

            # append this selection to mortality_p with the timestep that they were killed
            df = pd.DataFrame({'traj_id':mortality_selection, 'mortstep':i})
            mortality_p = mortality_p.append(df, ignore_index=True, sort=True)

    mortality_p = mortality_p.infer_objects()
    origin_dest = origin_dest.infer_objects()
    mortality_p.traj_id = mortality_p.traj_id.astype('float')
    origin_dest.traj_id = origin_dest.traj_id.astype('float')
    # join to origin_dest
    # This could be a M:1 join (more than 1 particle with the same traj_id)
    origin_dest_mort = origin_dest.merge(mortality_p, on='traj_id', how='outer')
    # we still want the join to happen even if mortality_p is empty. It will 
    # just make the column NaN which we later turn to -1.

    return origin_dest_mort, mortality_p

###################
# Add in starting time interval so that we know the full time period of 
# particles that may not have been released at the first time step.
# This allows us calculate different connections for different PLDs
###################

def start_time_int(origin_dest_mort, timesteps_with_release, particle_range, traj):
    
    logging.info("adding in particle start time")

    df = pd.DataFrame()
    df['traj_id'] = list(traj)

    # starting time step of each particle
    time_int_start = []
    for t in range(len(timesteps_with_release)):
        for particle in range(particle_range[t][0],particle_range[t][1]):
            time_int_start.append(timesteps_with_release[t])
    
    df['time_int_s'] = time_int_start # name shortened for shapefile

    df = df.infer_objects()
    origin_dest_mort = origin_dest_mort.infer_objects()
    df.traj_id = df.traj_id.astype('float')
    origin_dest_mort.traj_id = origin_dest_mort.traj_id.astype('float')   
    origin_dest_mort = origin_dest_mort.merge(df, on='traj_id')

    return origin_dest_mort

####################
# fill na with -1
# this is to prevent Arc from turning them to 0 on export
####################

def fill_na(origin_dest_mort):
    origin_dest_mort = origin_dest_mort.fillna(
        value={'dcoords': -1, 'dest_id': -1, 'mortstep': -1})
    return origin_dest_mort






###################
# OUTPUTS
###################

#### output destinaiton points to shapefile ####
def out_shp_dest_points(origin_dest_mort, crs_input_shp, shp_out, date_start):
    logging.info("writing points to shapefile")
    # can only have one geometry column
    # remove origin spatial column since for origin I am just concernced about origin poly ID
    od = origin_dest_mort.drop(['o_coords'], axis=1)
    od = geopandas.GeoDataFrame(od, geometry='d_coords')
    od.crs = crs_input_shp
    od['date_start'] = date_start
    od.to_file(filename=shp_out, driver='ESRI Shapefile')

#### create connection lines ####
def connection_lines(origin_dest_mort, shp_og, crs_input_shp, conn_lines_out, date_start, pld_int, pld):

    logging.info("writing connection lines to shapefile")
    od = origin_dest_mort
    sg = geopandas.read_file(shp_og)

    # on od, select particles where time_int_s minus time_int is less than or equal to PLD
    od_pld = od[(od.time_int - od.time_int_s <= pld_int)]

    # get each unique combination of originID and destID and get count of particles that survived
    od_unique = od_pld[(od_pld.dest_id != -1) & (od_pld.mortstep == -1)].groupby(['uID_part','dest_id']).size().reset_index(name='Freq')
    # how to read this:
    # first we select from od the ones that settled and survived
    # then we groupby unique combinations of uID and dest_id
    # then we get the count of those unique combinations
    # this normally makes uID the index and doesn't have a column name for count 
    # (the series we created), so we reset index and give the count a column name

    # df of time interval where first settlement occurred
    df_time_int = od_pld[(od_pld.dest_id != -1) & (od_pld.mortstep == -1)].groupby(['uID_part','dest_id'])['time_int'].min().reset_index(name='time_int')
    
    # set up for creating self connection lines. Size of circle lines based on amount settled and average area of all patches.
    def CircleCoords(xLeft, yCenter, r, n): # credit: MGET. Also, see circle_coords.xlsx for explanation of equation.
        return [(xLeft + r - math.cos(2*math.pi/n*x)*r, math.sin(2*math.pi/n*x)*r + yCenter) for x in range(n+1)]
    # min and max quantities used for normalization
    quantity_min = od_unique[od_unique.dest_id == od_unique.uID_part].Freq.min()
    quantity_max = od_unique[od_unique.dest_id == od_unique.uID_part].Freq.max()
    # get average area
    area_mean = sg.area.mean()
    # get radius of a circle with this area
    radius = math.sqrt(area_mean/math.pi)
    
    # for each unique combinaton create line from centroids
    connection_lines = pd.DataFrame(columns=['from_id','to_id','quantity','totalori','prob','time_int', 'line', 'pld'])
    conn_i = 0
    for row in od_unique.itertuples(index=False):
        # get total amount of particles released from patch
        total = od.uID_part[od.uID_part ==  row[0]].value_counts().values[0]

        # time interval where first settlement occurred
        time_int = df_time_int[(df_time_int.uID_part == row[0]) & (df_time_int.dest_id == row[1])]['time_int'].values[0]
    
        # get centroid of from and to patches
        centroid_origin = sg[sg.uID_202011 == row[0]].centroid
        centroid_dest = sg[sg.uID_202011 == row[1]].centroid

        if row[0] != row[1]:
            geom_line = LineString([centroid_origin.tolist()[0], centroid_dest.tolist()[0]])
        else:
            # normalize the quantites to 0.5 - 1 range (or I can do 0-1 but then the smallest one won't show up)
            #quantity_norm = 0.5 * (row[2] - quantity_min) / float(quantity_max - quantity_min) + 0.5
            quantity_norm = (row[2] - quantity_min) / float(quantity_max - quantity_min)
            radius_adj = radius * quantity_norm
            geom_line = LineString(CircleCoords(centroid_origin.x.tolist()[0], centroid_origin.y.tolist()[0], radius_adj, 90))
    
        connection_lines.loc[conn_i] = [row[0],row[1],float(row[2]),float(total),row[2]/float(total), time_int,geom_line, pld]
        conn_i += 1
    
    connection_lines['date_start'] = date_start   
    connection_lines = geopandas.GeoDataFrame(connection_lines, geometry='line')
    connection_lines.crs = crs_input_shp
    if len(connection_lines) == 0:
        logging.warning('Connection lines df is blank. Cannot write to shapefile')
        return
    connection_lines.to_file(filename=conn_lines_out, driver='ESRI Shapefile')

#### output patch centroids to shapefile (for use in network analysis) ####
def out_shp_patch_centroids(shp_og, patch_centroids_out, crs_input_shp, date_start):
    sg = geopandas.read_file(shp_og)
    # copy poly to new GeoDataFrame
    points = sg.copy()
    # change the geometry
    points.geometry = points['geometry'].centroid
    # same crs
    points.crs = crs_input_shp
    points['date_start'] = date_start
    points.to_file(filename=patch_centroids_out, driver='ESRI Shapefile')





###################
# RUN biology
# run all functions, even if you aren't applying settlement and/or mortality
###################

# commented out for running on Cedar. I just manually created the folders.
# when I run jobs at the same time, sometimes the conflict when checking and
# creating directories
# # create output shp folder
# if not os.path.exists(output_shp_dir):
#     os.makedirs(output_shp_dir)
dest_pts_dir = os.path.join(output_shp_dir, 'dest_pts')
# if not os.path.exists(dest_pts_dir):
#     os.makedirs(dest_pts_dir)
conn_lines_dir = os.path.join(output_shp_dir, 'conn_lines')
# if not os.path.exists(conn_lines_dir):
#     os.makedirs(conn_lines_dir)

# get nc files
nc_files = os.listdir(input_nc_dir)
ncs = []
for file in nc_files:
    ncs.append(os.path.join(input_nc_dir, file))

# reduce that range if if I don't want to run them all
if cluster_run:  # nc_group gets added for cluster runs
    for n in ncs:
        base = os.path.splitext(os.path.basename(n))[0]
        base = base.split('_')[1]
        if int(base) == nc_group:
            nc_index = ncs.index(n)
    ncs = ncs[nc_index:nc_index+1]
else:
    nc_index = None
    ncs = ncs[:]

# run biology functions for each release shapefile
for ncf in ncs:

    # get base name
    base = os.path.splitext(os.path.basename(ncf))[0]
    base = base.split('_')[1]
    logging.info("Processing group {}".format(base))

    # the lat/lon numpy files of starting coords saved from the opendrift run
    lat_np = os.path.join(input_npy_dir, 'lat_' + base + '.npy')
    lon_np = os.path.join(input_npy_dir, 'lon_' + base + '.npy')

    dataset = nc.Dataset(ncf, "r") # no '+'
    lon = dataset.variables["lon"]
    lat = dataset.variables["lat"]
    traj = dataset.variables["trajectory"]
    status = dataset.variables["status"]
    timestep = dataset.variables["time"]
    origin_marker = dataset.variables['origin_marker']
    #age_seconds = dataset.variables['age_seconds']
    # this will give the actual 30 SECOND interval that it was deactivated
    date_start = dataset.time_coverage_start

    particles_per_release = int(len(traj) / num_of_releases)

    origin = get_particle_originPoly(
        origin_marker, traj, crs_input_shp, lat_np, lon_np, backwards_run
        )

    timesteps_with_release, precomp_end_timestep, \
    precomp_range, particle_range = calc_precomp(
        precomp, time_step_output, particles_per_release, 
        interval_of_release, num_of_releases, traj
        )

    # put into for loop here, for each pld
    for pld in plds:

        logging.info("Processing pld {} for nc file {}".format(pld, base))
        # check that pld is not longer than length of timestep
        pld_int = int((pld * 24) / time_step_output)
        if pld_int > len(timestep):
            logging.warning("PLD provided is greater than length of time")

        origin_dest = settlement(
            settlement_apply, origin, shp_og_buff, timestep, status, lon, lat, traj, 
            crs_input_shp, precomp, precomp_range, particle_range, pld_int,
            timesteps_with_release
            )

        origin_dest = get_destination_coords(
            origin_dest, traj, lon, lat, timestep, crs_input_shp, status,
            particle_range, pld_int, timesteps_with_release
            )

        origin_dest_mort, mortality_p = calc_mortality(
            mortality_rate, traj, timestep, origin_dest, time_step_output, 
            mort_period, timesteps_with_release, particle_range
            )

        origin_dest_mort = start_time_int(
            origin_dest_mort, timesteps_with_release, particle_range, traj
            )

        origin_dest_mort = fill_na(origin_dest_mort)

        ### outputs
        shp_out = os.path.join(dest_pts_dir, 'dest_biology_pts_{}_pld{}.shp'.format(base, str(pld)))
        out_shp_dest_points(origin_dest_mort, crs_input_shp, shp_out, date_start)

        if settlement_apply:
            conn_lines_out = os.path.join(conn_lines_dir, 'connectivity_{}_pld{}.shp'.format(base, str(pld)))
            connection_lines(
                origin_dest_mort, shp_og, crs_input_shp, conn_lines_out, date_start, pld_int, pld
                )

# just create this once
if (nc_index == 0 or nc_index == None) and (pld_index == 0 or pld_index == None):
    patch_centroids_out = os.path.join(output_shp_dir, 'patch_centroids.shp')
    out_shp_patch_centroids(
        shp_og, patch_centroids_out, crs_input_shp, date_start
        )