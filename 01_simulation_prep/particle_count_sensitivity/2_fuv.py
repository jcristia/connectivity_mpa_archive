# Access the correlation matrix in the band statistics output file.
# Calculate FUV which is the fraction of unexplained variance.
# Plot FUV vs. Percentage of Points.


mpa_uid = 207

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


bandstats_dir = r'C:\Users\jcristia\Documents\GIS\MSc_Projects\MPA_connectivity\scripts\particle_count\analysis\band_statistics'
bandstats_txt = os.path.join(bandstats_dir, f'bandstats_mpa_{mpa_uid}.txt')


# the way it outputs now, I just need the 3rd to last line
# kinda hokey, but this is probably the easiest way to do it
file = open(bandstats_txt, 'r')
lines = file.readlines()
file.close()
line_list = lines[-2].split()
line_list = line_list[1:]
line_list = [float(i) for i in line_list]

# create dataframe
perc = list(range(5,105, 5))
df = pd.DataFrame({'perc_pts': perc, 'corr_coef': line_list})

# calc FUV
df['fuv'] = 1 - (df.corr_coef**2)

# plot
g = sns.lineplot(data=df, x='perc_pts', y='fuv')
# add line for 0.05
g.axhline(0.05, linewidth=0.5, color='r', linestyle='--')
g.figure.savefig(os.path.join(bandstats_dir, f'fuv_particles_{mpa_uid}.png'))


##### General interpretation

# Essentially, we want to see FUV drop below 0.05. That is the level of adequate
# points for capturing most of the variation in the system.
# However, if we never see this level out, then we likely need to release more
# points. We start by assuming that the total amount we release is more than
# enough, but if we still see large drops at higher percentages then we probably
# need to release more.

# My values may jump around a bit, but the key will be for when it consistently 
# gets below a certain value.
# To get below 0.05, it needs to have a correlation greater than 0.98.
# Why 0.05?
# Simons 2013 just says that it is reminiscint of a 95% confidence interval.
# In a way we are capturing 95% of the potential variance.
# I think a visual assessment is a valuable assessment as well - all we see at 
# some point is single point cells blink in and out.

# Lastly, why did I choose 5km cell size?
# Well if we chose something super small 500m, then things are probably not
# going to line up well even with lots of particles.
# For the purposes of MPA connectivity and successful recruitment, I'm most
# interested in the general area of settlement. Managment is likely done at this
# coarse scale as well.
# Also, the resolution of the hydro models (500m - 2500m) doesn't allow us to
# say much at or below these scales.

# For MPA 207, it looks like it just barely stabilizes before 100%. This is fine.
# I am releasing an adequate number of particles, but I probably wouldn't want
# to go lower.
# THEREFORE, my equation for releasing particles is good.