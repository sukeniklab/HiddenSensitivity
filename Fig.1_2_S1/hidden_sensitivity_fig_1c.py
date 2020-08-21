##########################################################################
# Figure 1C: Bar graph with swarm plot of experimental and simulated chi
#
# Supporting material for "Probing the Hidden Sensitivity of Intrinsically 
# Disordered Proteins to their Chemical Environment"
# https://www.biorxiv.org/content/10.1101/2020.08.17.252478v1
#
# Contact: David Moses dmoses5@ucmerced.edu
##########################################################################

from hidden_sensitivity import get_all_idrs_chi_df

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
import csv

# Parameters for figure
mpl.rc("figure", figsize=(9, 5))
plt.rc('xtick', labelsize=24) 
plt.rc('ytick', labelsize=24) 
plt.rcParams["axes.labelsize"] = 24
plt.rcParams["xtick.major.size"] = 8
plt.rcParams["ytick.major.size"] = 8
plt.rcParams["xtick.minor.size"] = 4
plt.rcParams["ytick.minor.size"] = 4
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['axes.linewidth'] = 2
plt.rcParams['legend.fontsize'] = 22
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
sns.set_style("white")
mypallet = sns.color_palette(['tab:red','tab:purple'])
mypallet2 = sns.color_palette(['tab:pink','fuchsia'])

# Experimental data -- build data frame containing chi values for all IDRs in all solution conditions
all_idrs_chi_df = get_all_idrs_chi_df()
Chis=pd.DataFrame(columns=['idr','exp_sim','chi'])
for index, row in all_idrs_chi_df[((all_idrs_chi_df.sol=='Buffer (standard)') & (all_idrs_chi_df.idr.isin(['PUMA','Ash1','p53','E1A','GS24'])))].iterrows():
    Chis = Chis.append({'idr':row.idr,'exp_sim':'experiment','chi':row.chi}, ignore_index=True)

# Simulation data
with open('sim-idr-chi.csv', newline='') as f:
    reader = csv.reader(f)
    sim_idr_chi = list(reader)
for row in sim_idr_chi:
    Chis = Chis.append({'idr':row[0],'exp_sim':'simulation','chi':float(row[1])}, ignore_index=True)

# Make barplot and swarmplot
fig = plt.figure()
ax = sns.barplot(x="idr", y="chi", hue="exp_sim", data=Chis, alpha=1, errwidth=3, errcolor='k', palette=mypallet)
sns.stripplot(x="idr", y="chi", hue="exp_sim", data=Chis, dodge=True, palette=mypallet2, ax=ax, 
              size=7, alpha=1, jitter=0.3, edgecolor='gray', linewidth=0.5)
# sns.swarmplot(x="idr", y="chi", hue="exp_sim", data=Chis, dodge=True, color='k', ax=ax, size=10, alpha=0.4)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('$\chi$')
ax.set_xlabel('')
ax.set_yticks([-0.4,-0.2,0,0.2,0.4,0.6,0.8])
ax.set_yticklabels(['-0.4','-0.2','0','0.2','0.4','0.6','0.8'])
ax.tick_params(direction='in')
ax.xaxis.tick_bottom()
ax.yaxis.tick_left()

# Dotted line at chi = 0
x_vals = np.array(ax.get_xlim())
sns.lineplot(x_vals, [0,0], ax=ax, color="black")
ax.lines[10].set_linestyle("--")

# Fix the legend
handles, labels = ax.get_legend_handles_labels()
handles = handles[2:4]
labels = labels[2:4]
ax.legend(handles, labels, loc='upper left')

plt.savefig('hidden_sensitivity_fig_1c.png', bbox_inches='tight')
plt.show()
