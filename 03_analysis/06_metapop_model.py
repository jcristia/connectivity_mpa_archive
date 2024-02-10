# Metapopulation model with stochasticity

# Based on Patricks metacommunity model:
# https://github.com/plthompson/mcomsimr/blob/master/R/MC_simulate.R

import arcpy
import pandas as pd
import numpy as np
import seaborn as sns
from datetime import datetime



#### Inputs ####

root = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity'
conn = r'cluster_results\scripts\COMBINED.gdb\conn_avg_pld{}' # connectivity probabilities per PLD averaged over time
patches = os.path.join(root, r'spatial\MPA\mpas.gdb\M09_mpa_joined')
out_gdb = os.path.join(root, r'scripts\metapop_pers\metapop_pers.gdb')

timesteps = 1750
#plds = [3,7,10]
plds = [21,30,40]
mort_rates = [0.1] # mortality rate
#prop_reproducing = 0.5 # proportion of the population the reproduces in a time step
larvae_dispersing_per_adult = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35] # doing it this way with 
# fractional amounts allows us to eliminate the prop_reproducing variable. By 
# allowing fractional amounts, it represents that not all of the population 
# reproduces at each time step but still allows us to have numbers greater than 1
remove_retention = False



#### Get starting population based on patch area ####

field_names = [i.name for i in arcpy.ListFields(patches) if i.type != 'OID']
field_names = [i.name for i in arcpy.ListFields(patches) if i.name in ['uID_20201124', 'Shape_Area']]
cursor = arcpy.da.SearchCursor(patches, field_names)
patches = pd.DataFrame(data=[row for row in cursor], columns=field_names)

# remove mpas that are in exclude list
# Connections to/from these mpas are minimal, so there is little impact anyway
mpas_exclude = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas.gdb\M10_toexcludefromanalysis'
to_include = []
with arcpy.da.SearchCursor(mpas_exclude, ['uID_20201124', 'exclude']) as cursor:
    for row in cursor:
        if not row[1] == 1:
            to_include.append(row[0])
patches = patches[patches.uID_20201124.isin(to_include)]

# assign proportionally with largest patch receiving 100 million individuals
# I'll probably want to scale this differently later
max_area = patches.Shape_Area.max()
patches['popn'] = (patches.Shape_Area * 100000000) / max_area
# for any patches with fewer than 1000 individuals, assign them 100
patches.loc[patches.popn < 1000, 'popn'] = 1000
patches = patches.round({'popn':0})

# convert to array
patches = patches.sort_values(by=['uID_20201124'])
patch_pop = np.array(patches.popn)



#### Simulation for each PLD ####

# df to hold output data
out_patches = patches

start_time = datetime.now()
temp_start = start_time

for pld in plds:

    #### read in connectivity data as matrix ####
    print(f'processing pld {pld}')

    # feature class to pandas df
    conns = os.path.join(root, conn.format(pld))
    field_names = [i.name for i in arcpy.ListFields(conns) if i.type != 'OID']
    cursor = arcpy.da.SearchCursor(conns, field_names)
    df = pd.DataFrame(data=[row for row in cursor], columns=field_names)

    # but now remove the mpas to exclude that are deep in inlets or narrow areas
    df = df[(df.from_id.isin(to_include)) & (df.to_id.isin(to_include))]

    # add missing mpas
    # there are a few mpas that have nothing coming into them, which means
    # I don't get a square matrix
    #i = 1
    # while i < 407:
    #     if i not in df.to_id.unique():
    #         df = df.append({'to_id':i, 'from_id':i}, ignore_index=True)
    #     i += 1
    for i in to_include:
        if i not in df.to_id.unique():
            df = df.append({'to_id':i, 'from_id':i}, ignore_index=True)


    # to pandas matrix df
    df_p = df.pivot(index='to_id', columns='from_id', values='prob_avg')

    # pandas df to numpy matrix
    conn_matrix_orig = df_p.to_numpy()
    conn_matrix_orig = np.nan_to_num(conn_matrix_orig)
    if np.shape(conn_matrix_orig)[0] != np.shape(conn_matrix_orig)[1]:
        raise ValueError('Not a square matrix')

    # testing: remove retention:
    if remove_retention:
        np.fill_diagonal(conn_matrix_orig, 0)

    for proportion in larvae_dispersing_per_adult:

        print(f'processing dispersal proportion {proportion}')

        for mort in mort_rates:

            print(f'processing mortality rate {mort}')

            # reset population and conn_matrix to start
            pop_net = patch_pop
            conn_matrix = conn_matrix_orig

            # dataframe for correlations
            df_cc = pd.DataFrame(columns={'timestep', 'cc', 'prop_occ'})

            for i in list(range(timesteps)):

                if i % 50 == 0:
                    print(f'time step {i}')

                ### apply mortality ####
                pop_net = pop_net * (1-mort)
                #### expected population size ###
                # add demographic stochasticity
                # use poisson distribution to ensure I get a whole number of individuals
                pop_adj = np.random.poisson(pop_net)

                ### amount that is reproduced and disperses ###
                # no longer using a binomial distribution on the entire pop
                pop_disp = pop_adj * proportion
                # add stochasticity around the amount that disperses
                # ensure that I get whole number of individuals dispersing
                pop_disp = np.random.poisson(pop_disp)

                ### immigration ###
                # approach if we want to move whole individuals
                if i==0:
                    newrow = 1 - np.sum(conn_matrix, axis=0) # add a new row so columns add to 1
                    conn_matrix = np.vstack([conn_matrix, newrow])
                # assign each particle based on dispersal probabilities
                def randomChoice(col):
                    samp = np.random.choice(a=conn_matrix.shape[0], size=pop_disp[col], p=conn_matrix[:,col])
                    unique, counts = np.unique(samp, return_counts=True)
                    samp = np.zeros(conn_matrix.shape[0])
                    np.put(samp, unique, counts)
                    samp = samp[:-1]
                    return(samp)
                # for each column, figure out where its particles go
                samp_all = map(randomChoice, list(range(conn_matrix.shape[1])))
                samp_all = list(samp_all)
                immigrate = sum(samp_all)

                ### net ###
                net = pop_adj + immigrate
                pop_net = net
                pop_net[pop_net<0] = 0

                # apply carrying capacity, simply reduce it to that amount if over
                # K = patch_pop
                pop_net = np.where(pop_net > patch_pop, patch_pop, pop_net)


                #### track equilibrium ####
                # set any patches with individuals to 1
                # compare matrices between timesteps
                # patches will blink in and out, but eventually the matrix will reach equilibrium
                pers = np.where(pop_net > 0, 1, 0)
                if (i != 0) and (i % 5 == 0): # not the first timestep and only every 5 time steps
                    cc_m = np.corrcoef([pers_prev, pers])
                    cc = cc_m[1,0]
                    if np.isnan(cc):
                        cc = 1 # if there is no variance (e.g. all values are 1), then numpy outputs nan
                    proportion_occupied = np.sum(pers) / len(pers)
                    #print(proportion_occupied)
                    df_cc = df_cc.append({'timestep':i, 'cc': cc, 'prop_occ':proportion_occupied}, ignore_index=True)
                pers_prev = pers


            ## plot correlations ###
            df_cc = df_cc.astype({'timestep':'int64', 'cc':'float32', 'prop_occ':'float32'})
            g = sns.lineplot(data=df_cc, x='timestep', y='cc')
            g.set(ylim=(0.91, 1.01))
            str_prop = str(proportion).replace('.', '')
            str_mort = str(mort).split('.')[1]
            g.figure
            if remove_retention:
                out_string = f'cor_pld{pld}_prop{str_prop}_m{str_mort}_noret.png'
            else:
                out_string = f'cor_pld{pld}_prop{str_prop}_m{str_mort}_ret.png'
            g.figure.savefig(os.path.join(root,'scripts/metapop_pers/figs', out_string))
            g.figure.clf()
            
            ## plot proportion occupied ##
            f = sns.lineplot(data=df_cc, x='timestep', y='prop_occ')
            f.figure
            if remove_retention:
                out_string = f'pro_pld{pld}_prop{str_prop}_m{str_mort}_noret.png'
            else:
                out_string = f'pro_pld{pld}_prop{str_prop}_m{str_mort}_ret.png'
            f.figure.savefig(os.path.join(root,'scripts/metapop_pers/figs', out_string))
            f.figure.clf()

            #### add persistence to out_patches ####
            out_patches[f'pld{pld}_prop{str_prop}_m{str_mort}'] = pers

            # time for 1 run:
            print('Runtime: {}'.format(datetime.now() - temp_start))
            temp_start = datetime.now()



print("Total Runtime: {}".format(datetime.now() - start_time))


### output to csv file ###
# doing this way allows me to do simulations separately and combine them later
out_patches = out_patches.drop(columns=['Shape_Area'])
current_time = datetime.now().strftime("%y%m%d%H%M")
out_string = f'metapop_pers_{timesteps}_{current_time}.csv'
out_patches.to_csv('csv/'+out_string)


### output to feature class table ###
# out_patches = out_patches.drop(columns=['Shape_Area'])
# x = np.array(np.rec.fromrecords(out_patches.values))
# names = out_patches.dtypes.index.tolist()
# x.dtype.names = tuple(names)
# out_string = f'metapop_pers_{timesteps}'
# arcpy.da.NumPyArrayToTable(x, os.path.join(out_gdb, out_string))
