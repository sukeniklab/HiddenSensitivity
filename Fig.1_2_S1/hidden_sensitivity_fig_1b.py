##############################################################################
# Figure 1B: E_fret in buffer of GS linkers vs. length of amino acid chain (N)
#
# Supporting material for "Probing the Hidden Sensitivity of Intrinsically 
# Disordered Proteins to their Chemical Environment"
# https://www.biorxiv.org/content/10.1101/2020.08.17.252478v1
#
# Contact: David Moses dmoses5@ucmerced.edu
##############################################################################

from hidden_sensitivity import build_base_correction_factor_df
from hidden_sensitivity import get_efret
from hidden_sensitivity import convert_efret_df_to_chi_df
from hidden_sensitivity import preprocess_tq_ng_data

import pandas as pd
import numpy as np
import matplotlib as mpl
from scipy import stats
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
import random
plt.rc('xtick', labelsize=24) 
plt.rc('ytick', labelsize=24) 
plt.rcParams["axes.labelsize"] = 24
mpl.rcParams['axes.linewidth'] = 2
plt.rcParams["xtick.major.size"] = 8
plt.rcParams["ytick.major.size"] = 8
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
sns.set_style("white")

# Read in and preprocess raw data for mTurquoise2 and mNeonGreen base spectra
TQ_data, NG_data = preprocess_tq_ng_data()

# Read in raw data for selected IDRs
GS0_data_1=pd.read_csv("GS0-1-repeat-1.csv",sep=",",skiprows=1).join(pd.read_csv('GS0-2-repeat-1.csv',
                    skiprows=1).iloc[:,3:],rsuffix='2').join(pd.read_csv('GS0-3-repeat-1.csv',
                    skiprows=1).iloc[:,3:],rsuffix='3')
GS0_data_2=pd.read_csv("GS0-1-repeat-2.csv",sep=",",skiprows=1).join(pd.read_csv('GS0-2-repeat-2.csv',
                    skiprows=1).iloc[:,3:],rsuffix='2').join(pd.read_csv('GS0-3-repeat-2.csv',
                    skiprows=1).iloc[:,3:],rsuffix='3')
GS8_data_1=pd.read_csv("GS8-1-repeat-1.csv",sep=",",skiprows=1).join(pd.read_csv('GS8-2-repeat-1.csv',
                    skiprows=1).iloc[:,3:],rsuffix='2').join(pd.read_csv('GS8-3-repeat-1.csv',
                    skiprows=1).iloc[:,3:],rsuffix='3')
GS8_data_2=pd.read_csv("GS8-1-repeat-2.csv",sep=",",skiprows=1).join(pd.read_csv('GS8-2-repeat-2.csv',
                    skiprows=1).iloc[:,3:],rsuffix='2').join(pd.read_csv('GS8-3-repeat-2.csv',
                    skiprows=1).iloc[:,3:],rsuffix='3')
GS16_data_1=pd.read_csv("GS16-1-repeat-1.csv",sep=",",skiprows=1).join(pd.read_csv('GS16-2-repeat-1.csv',
                    skiprows=1).iloc[:,3:],rsuffix='2').join(pd.read_csv('GS16-3-repeat-1.csv',
                    skiprows=1).iloc[:,3:],rsuffix='3')
GS16_data_2=pd.read_csv("GS16-1-repeat-2.csv",sep=",",skiprows=1).join(pd.read_csv('GS16-2-repeat-2.csv',
                    skiprows=1).iloc[:,3:],rsuffix='2').join(pd.read_csv('GS16-3-repeat-2.csv',
                    skiprows=1).iloc[:,3:],rsuffix='3')
GS24_data_1=pd.read_csv("GS24-1-repeat-1.csv",sep=",",skiprows=1).join(pd.read_csv('GS24-2-repeat-1.csv',
                    skiprows=1).iloc[:,3:],rsuffix='2').join(pd.read_csv('GS24-3-repeat-1.csv',
                    skiprows=1).iloc[:,3:],rsuffix='3')
GS24_data_2=pd.read_csv("GS24-1-repeat-2.csv",sep=",",skiprows=1).join(pd.read_csv('GS24-2-repeat-2.csv',
                    skiprows=1).iloc[:,3:],rsuffix='2').join(pd.read_csv('GS24-3-repeat-2.csv',
                    skiprows=1).iloc[:,3:],rsuffix='3')
GS32_data_1=pd.read_csv("GS32-1-repeat-1.csv",sep=",",skiprows=1).join(pd.read_csv('GS32-2-repeat-1.csv',
                    skiprows=1).iloc[:,3:],rsuffix='2').join(pd.read_csv('GS32-3-repeat-1.csv',
                    skiprows=1).iloc[:,3:],rsuffix='3')
GS32_data_2=pd.read_csv("GS32-1-repeat-2.csv",sep=",",skiprows=1).join(pd.read_csv('GS32-2-repeat-2.csv',
                    skiprows=1).iloc[:,3:],rsuffix='2').join(pd.read_csv('GS32-3-repeat-2.csv',
                    skiprows=1).iloc[:,3:],rsuffix='3')
GS48_data_1=pd.read_csv("GS48-1-repeat-1.csv",sep=",",skiprows=1).join(pd.read_csv('GS48-2-repeat-1.csv',
                    skiprows=1).iloc[:,3:],rsuffix='2').join(pd.read_csv('GS48-3-repeat-1.csv',
                    skiprows=1).iloc[:,3:],rsuffix='3')
GS48_data_2=pd.read_csv("GS48-1-repeat-2.csv",sep=",",skiprows=1).join(pd.read_csv('GS48-2-repeat-2.csv',
                    skiprows=1).iloc[:,3:],rsuffix='2').join(pd.read_csv('GS48-3-repeat-2.csv',
                    skiprows=1).iloc[:,3:],rsuffix='3')
TQNG_data_1=pd.read_csv("TQNGfree-1-repeat-1.csv",sep=",",skiprows=1).join(pd.read_csv('TQNGfree-2-repeat-1.csv',
                    skiprows=1).iloc[:,3:],rsuffix='2').join(pd.read_csv('TQNGfree-3-repeat-1.csv',
                    skiprows=1).iloc[:,3:],rsuffix='3')

# Build a data frame with correction factors for the donor and acceptor at each concentration of each solute
base_corr_fact_df=build_base_correction_factor_df(TQ_data, NG_data)

# Build an E_fret data frame for each repeat for each IDR
GS0_efret_1=get_efret('GS0',GS0_data_1,range(3,GS0_data_1.shape[1]),base_corr_fact_df,TQ_data,NG_data,True)
GS0_efret_2=get_efret('GS0',GS0_data_2,range(3,GS0_data_2.shape[1]),base_corr_fact_df,TQ_data,NG_data,True)
GS8_efret_1=get_efret('GS8',GS8_data_1,range(3,GS8_data_1.shape[1]),base_corr_fact_df,TQ_data,NG_data,True)
GS8_efret_2=get_efret('GS8',GS8_data_2,range(3,GS8_data_2.shape[1]),base_corr_fact_df,TQ_data,NG_data,True)
GS16_efret_1=get_efret('GS16',GS16_data_1,range(3,GS16_data_1.shape[1]),base_corr_fact_df,TQ_data,NG_data,True)
GS16_efret_2=get_efret('GS16',GS16_data_2,range(3,GS16_data_2.shape[1]),base_corr_fact_df,TQ_data,NG_data,True)
GS24_efret_1=get_efret('GS24',GS24_data_1,range(3,GS24_data_1.shape[1]),base_corr_fact_df,TQ_data,NG_data,True)
GS24_efret_2=get_efret('GS24',GS24_data_2,range(3,GS24_data_2.shape[1]),base_corr_fact_df,TQ_data,NG_data,True)
GS32_efret_1=get_efret('GS32',GS32_data_1,range(3,GS32_data_1.shape[1]),base_corr_fact_df,TQ_data,NG_data,True)
GS32_efret_2=get_efret('GS32',GS32_data_2,range(3,GS32_data_2.shape[1]),base_corr_fact_df,TQ_data,NG_data,True)
GS48_efret_1=get_efret('GS48',GS48_data_1,range(3,GS48_data_1.shape[1]),base_corr_fact_df,TQ_data,NG_data,True)
GS48_efret_2=get_efret('GS48',GS48_data_2,range(3,GS48_data_2.shape[1]),base_corr_fact_df,TQ_data,NG_data,True)
TQNG_efret_1=get_efret('TQNG',TQNG_data_1,range(3,GS0_data_1.shape[1]),base_corr_fact_df,TQ_data,NG_data,True) 

################
# Make the graph
################

# X values for GS linkers and other IDRs
Ns = np.array([0,16,32,48,64,96,150])
idrs = ['Tethered', 'GS8','GS16','GS24','GS32','GS48','Untethered']

# Set up for broken axis
fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(8,5))
axs[0]=plt.subplot2grid((1, 20), (0, 1), rowspan=1, colspan=15)
axs[1]=fig.add_subplot(1,8,8)

# Hide the spines between ax and ax2
axs[0].spines['right'].set_visible(False)
axs[1].spines['left'].set_visible(False)
axs[0].yaxis.tick_left()
axs[1].set_yticklabels([])

# Divide up the data
axs[0].set_xlim(-5,111)
axs[1].set_xlim(142.5,157.5)

# Labels and ticks
axs[0].set(xlabel='N', ylabel='$E_{f}$')
plt.sca(axs[0])
plt.xticks(Ns[:6], Ns[:6])
plt.sca(axs[1])
plt.xticks(Ns[6:7], '')
axs[0].tick_params(direction='in')
axs[1].tick_params(direction='in')
fig.text(0.856, 0.056, 'UT', ha='center', fontsize=24)

d = .01 # how big to make the diagonal lines in axes coordinates
# Arguments to pass plot, just so we don't keep repeating them
kwargs = dict(transform=axs[0].transAxes, color='k', clip_on=False)
axs[0].plot((1-d,1+d), (-2*d,+1.5*d), **kwargs)
axs[0].plot((1-d,1+d),(1-2*d,1+1.5*d), **kwargs)
kwargs.update(transform=axs[1].transAxes)  # switch to the bottom axes
axs[1].plot((-8*d,+7.5*d), (1-2.1*d,1+1.6*d), **kwargs)
axs[1].plot((-8*d,+7.5*d), (-2.1*d,+1.6*d), **kwargs)

# Get data for swarmplot
data = pd.DataFrame({ "GS0":pd.concat([GS0_efret_1[GS0_efret_1['sol'].isin(['Buffer','Buffer.1','Buffer.2','Buffer.3','Buffer.4','Buffer.5',
                                                             'Buffer2','Buffer.12','Buffer.22','Buffer.32','Buffer.42','Buffer.52'])]['E_f'],
                                       GS0_efret_2[GS0_efret_2['sol'].isin(['Buffer','Buffer.1','Buffer.2','Buffer.3','Buffer.4','Buffer.5',
                                                             'Buffer2','Buffer.12','Buffer.22','Buffer.32','Buffer.42','Buffer.52'])]['E_f']], 
                                      axis=0, ignore_index=True),
                      "GS8":pd.concat([GS8_efret_1[GS8_efret_1['sol'].isin(['Buffer','Buffer.1','Buffer.2','Buffer.3','Buffer.4','Buffer.5',
                                                             'Buffer2','Buffer.12','Buffer.22','Buffer.32','Buffer.42','Buffer.52'])]['E_f'],
                                       GS8_efret_2[GS8_efret_2['sol'].isin(['Buffer','Buffer.1','Buffer.2','Buffer.3','Buffer.4','Buffer.5',
                                                             'Buffer2','Buffer.12','Buffer.22','Buffer.32','Buffer.42','Buffer.52'])]['E_f']], 
                                      axis=0, ignore_index=True),
                      "GS16":pd.concat([GS16_efret_1[GS16_efret_1['sol'].isin(['Buffer','Buffer.1','Buffer.2','Buffer.3','Buffer.4','Buffer.5',
                                                             'Buffer2','Buffer.12','Buffer.22','Buffer.32','Buffer.42','Buffer.52'])]['E_f'],
                                        GS16_efret_2[GS16_efret_2['sol'].isin(['Buffer','Buffer.1','Buffer.2','Buffer.3','Buffer.4','Buffer.5',
                                                             'Buffer2','Buffer.12','Buffer.22','Buffer.32','Buffer.42','Buffer.52'])]['E_f']], 
                                        axis=0, ignore_index=True),
                      "GS24":pd.concat([GS24_efret_1[GS24_efret_1['sol'].isin(['Buffer','Buffer.1','Buffer.2','Buffer.3','Buffer.4','Buffer.5',
                                                             'Buffer2','Buffer.12','Buffer.22','Buffer.32','Buffer.42','Buffer.52'])]['E_f'],
                                        GS24_efret_2[GS24_efret_2['sol'].isin(['Buffer','Buffer.1','Buffer.2','Buffer.3','Buffer.4','Buffer.5',
                                                             'Buffer2','Buffer.12','Buffer.22','Buffer.32','Buffer.42','Buffer.52'])]['E_f']], 
                                        axis=0, ignore_index=True), 
                      "GS32":pd.concat([GS32_efret_1[GS32_efret_1['sol'].isin(['Buffer','Buffer.1','Buffer.2','Buffer.3','Buffer.4','Buffer.5',
                                                             'Buffer2','Buffer.12','Buffer.22','Buffer.32','Buffer.42','Buffer.52'])]['E_f'],
                                        GS32_efret_2[GS32_efret_2['sol'].isin(['Buffer','Buffer.1','Buffer.2','Buffer.3','Buffer.4','Buffer.5',
                                                             'Buffer2','Buffer.12','Buffer.22','Buffer.32','Buffer.42','Buffer.52'])]['E_f']], 
                                        axis=0, ignore_index=True),
                      "GS48":pd.concat([GS48_efret_1[GS48_efret_1['sol'].isin(['Buffer','Buffer.1','Buffer.2','Buffer.3','Buffer.4','Buffer.5',
                                                             'Buffer2','Buffer.12','Buffer.22','Buffer.32','Buffer.42','Buffer.52'])]['E_f'],
                                        GS48_efret_2[GS48_efret_2['sol'].isin(['Buffer','Buffer.1','Buffer.2','Buffer.3','Buffer.4','Buffer.5',
                                                             'Buffer2','Buffer.12','Buffer.22','Buffer.32','Buffer.42','Buffer.52'])]['E_f']], 
                                        axis=0, ignore_index=True),
                      "GSinf":TQNG_efret_1[TQNG_efret_1['sol'].isin(['Buffer','Buffer.1','Buffer.2','Buffer.3','Buffer.4','Buffer.5',
                                                             'Buffer2','Buffer.12','Buffer.22','Buffer.32','Buffer.42','Buffer.52'])]['E_f']}) 

# Average Efret in buffer for GS linkers
efrets_means = []
efrets_means.append(data.GS0.mean())
efrets_means.append(data.GS8.mean())
efrets_means.append(data.GS16.mean())
efrets_means.append(data.GS24.mean())
efrets_means.append(data.GS32.mean())
efrets_means.append(data.GS48.mean())
efrets_means.append(data.GSinf.mean())
efrets_means = np.array(efrets_means)

# Standard deviation of Efret in buffer for GS linkers
efrets_stdevs = []
efrets_stdevs.append(data.GS0.std())
efrets_stdevs.append(data.GS8.std())
efrets_stdevs.append(data.GS16.std())
efrets_stdevs.append(data.GS24.std())
efrets_stdevs.append(data.GS32.std())
efrets_stdevs.append(data.GS48.std())
efrets_stdevs.append(data.GSinf.std())
efrets_stdevs = np.array(efrets_stdevs)

# After linear regression of GS Efret values, Efret(GS) = slope * N + intercept 
# Can plot line if needed
slope, intercept, r_value, p_value, std_err = stats.linregress(Ns[0:5], efrets_means[0:5])
x_vals = np.array(axs[0].get_xlim())
y_vals = intercept + (slope * x_vals)
axs[0].plot(x_vals, y_vals, c='slategray')
# sns.lineplot(x_vals, y_vals, ax=axs[0], c='slategray')
axs[0].lines[2].set_linestyle("--")

# Plot the averages (not normalized) on both axes
group = np.array(idrs)
for g in np.unique(group):
    ix = np.where(group == g)
    axs[0].errorbar(Ns[ix], efrets_means[ix], efrets_stdevs[ix], marker='o', c='black', barsabove=True, ecolor='silver', label=g, ms=10, alpha=1, zorder=100)
    axs[1].errorbar(Ns[ix], efrets_means[ix], efrets_stdevs[ix], marker='o', c='black', barsabove=True, ecolor='silver', label=g, ms=10, alpha=1, zorder=100)

# Arrange data for swarmplot
swarm_df = pd.DataFrame(columns=['idr','N','E_f'])
for x in data.GS0:
    swarm_df = swarm_df.append({'idr':idrs[0], 'N':Ns[0] + random.uniform(-5, 5), 'E_f':x}, ignore_index=True)
for x in data.GS8:
    swarm_df = swarm_df.append({'idr':idrs[1], 'N':Ns[1] + random.uniform(-5, 5), 'E_f':x}, ignore_index=True)
for x in data.GS16:
    swarm_df = swarm_df.append({'idr':idrs[2], 'N':Ns[2] + random.uniform(-5, 5), 'E_f':x}, ignore_index=True)
for x in data.GS24:
    swarm_df = swarm_df.append({'idr':idrs[3], 'N':Ns[3] + random.uniform(-5, 5), 'E_f':x}, ignore_index=True)
for x in data.GS32:
    swarm_df = swarm_df.append({'idr':idrs[4], 'N':Ns[4] + random.uniform(-5, 5), 'E_f':x}, ignore_index=True)
for x in data.GS48:
    swarm_df = swarm_df.append({'idr':idrs[5], 'N':Ns[5] + random.uniform(-5, 5), 'E_f':x}, ignore_index=True)
for x in data.GSinf:
    swarm_df = swarm_df.append({'idr':idrs[6], 'N':Ns[6] + random.uniform(-5, 5), 'E_f':x}, ignore_index=True)
    
# Plot the swarms on both axes
swarm_E_f = swarm_df.E_f.to_numpy()
swarm_Ns = swarm_df.N.to_numpy()
swarm_idrs = swarm_df.idr.to_numpy()
swarm_group = np.array(swarm_idrs)

for g in np.unique(swarm_group):
    ix = np.where(swarm_group == g)
    axs[0].scatter([swarm_Ns[x] for x in ix], [swarm_E_f[x] for x in ix], c = 'darkorchid', s = 40, alpha=0.4, zorder=1)
    axs[1].scatter([swarm_Ns[x] for x in ix], [swarm_E_f[x] for x in ix], c = 'darkorchid', s = 40, alpha=0.4, zorder=1)
    
for ax in axs:
    ax.xaxis.tick_bottom()
axs[0].yaxis.tick_left()

axs[0].text(0.1,0.12,r'$E_f$ = {:.5f} $\times$ N + {:.3f}'.format(slope,intercept),horizontalalignment='left',fontsize=24)
axs[0].text(0.1,0.04,'$R^2$ = {:.3f}'.format(r_value**2),horizontalalignment='left',fontsize=24)

plt.savefig('hidden_sensitivity_fig_1b.png', bbox_inches='tight')
plt.show()
