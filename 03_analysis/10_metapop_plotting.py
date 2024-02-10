# plots
# proportion of persistent populations by PLD, dispersal proportion, and r_0


import arcpy
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

gdb = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\metapop_pers\metapop_pers.gdb'
arcpy.env.workspace = gdb
metapop_pers_pts = f'metapop_pers_centroids'
mpas_orig = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\spatial\MPA\mpas.gdb\M09_mpa_joined'

canada_only = False
plds = [1,3,7,10,21,30,40,60]
larvae_dispersing_per_adult = [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4]
mort_rates = [0.05, 0.10, 0.15]


# read in pts as pandas df
field_names = [i.name for i in arcpy.ListFields(metapop_pers_pts)]
cursor = arcpy.da.SearchCursor(metapop_pers_pts, field_names)
persistence_df = pd.DataFrame(data=[row for row in cursor], columns=field_names)

# process only Canadian mpas?
if canada_only:

    # read in mpas_orig dataset as df
    field_names = [i.name for i in arcpy.ListFields(mpas_orig)]
    cursor = arcpy.da.SearchCursor(mpas_orig, field_names)
    mpas_orig = pd.DataFrame(data=[row for row in cursor], columns=field_names)

    # remove US ones from metapop df
    uids_keep = mpas_orig.uID_20201124[mpas_orig.provstate=='British Columbia']
    persistence_df = persistence_df[persistence_df.uID_202011.isin(uids_keep)]

# get total length
mpa_total_count = len(persistence_df)

# create blank df to append to
df_summary = pd.DataFrame(columns={'canada_only', 'pld', 'disp_prop', 'mortality_rate', 'persistent_percent'})

# for each persistence field
# get the pld, prop, r values
# calculate % persistent
# append to df
for pld in plds:
    for prop in larvae_dispersing_per_adult:
        for mort in mort_rates:

            str_prop = str(prop).replace('.', '')
            str_mort = str(mort).split('.')[1]
            field = f'pld{pld}_prop{str_prop}_m{str_mort}'

            # check if field exists. If not, set perc_pers to zero.
            if field in list(persistence_df.columns):
                perc_pers = persistence_df[field].sum() / mpa_total_count * 100
            else:
                perc_pers = pd.to_numeric('')

            df_summary = df_summary.append({
                'canada_only': str(canada_only),
                'pld':pld, 
                'disp_prop':prop, 
                'mortality_rate':mort, 
                'persistent_percent':perc_pers
                }, ignore_index=True)


# plots (multiple plots across mortality rates)
sns.set()
sns.set_style('white')
sns.set_context('paper')
f = sns.relplot(
    data = df_summary,
    x = 'disp_prop',
    y = 'persistent_percent',
    hue = 'pld',
    col='mortality_rate',
    kind='line',
    palette = 'tab10',
    marker='o'
)
#f.set(xscale='log')
f.set(xlabel='Proportion of population dispersing at each timestep', ylabel='% of MPAs persistent')
f.savefig('figs/mpas_persistent.svg')


# Final plot for just 1 mortality rate:
df_final = df_summary[df_summary.mortality_rate == 0.1]
sns.set()
sns.set_style('white')
sns.set_context('paper', font_scale=1.25, rc={"lines.linewidth": 2})
f = sns.relplot(
    data = df_final,
    x = 'disp_prop',
    y = 'persistent_percent',
    hue = 'pld',
    kind='line',
    palette = 'tab10',
    marker='o'
)
f._legend.remove()
plt.legend(title='PLD')
f.set(xlabel='Proportion of population produced and \n dispersing at each timestep ($\it{d}$)', ylabel='% of MPAs persistent')
f.savefig('figs/fig03_mpaspersistent.svg')
