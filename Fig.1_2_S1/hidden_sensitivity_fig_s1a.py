##################################################################################
# Figure S1A: Heatmap of chi values for GS linkers in various solution conditions
#
# Supporting material for "Probing the Hidden Sensitivity of Intrinsically 
# Disordered Proteins to their Chemical Environment"
# https://www.biorxiv.org/content/10.1101/2020.08.17.252478v1
#
# Contact: David Moses dmoses5@ucmerced.edu
##################################################################################

from hidden_sensitivity import get_all_idrs_chi_df

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.cm as cm
from scipy import stats
from matplotlib.ticker import FormatStrFormatter
plt.rc('xtick', labelsize=50) 
plt.rc('ytick', labelsize=50) 
plt.rcParams["xtick.major.size"] = 12
plt.rcParams["ytick.major.size"] = 12
mpl.rcParams['xtick.major.width'] = 3
mpl.rcParams['ytick.major.width'] = 3
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['axes.linewidth'] = 2
mpl.rc('axes',edgecolor='white')
sns.set_style("white")

###########################
# Make colormap for heatmap
###########################
minl=-1.3
maxl=1.3
cspace=np.linspace(minl,maxl,30)
location = np.argmin(np.abs(cspace))
# Blue-white-red
R=np.hstack([np.ones(location), np.linspace(1,0,30-location)])
G=np.hstack([np.linspace(0,1,location),np.linspace(1,0,30-location)])
B=np.hstack([np.linspace(0,1,location),np.ones(30-location)])
#green-white-magenta
# R=np.hstack([np.linspace(0,1,location),np.ones(30-location)])
# G=np.hstack([np.ones(location),np.linspace(1,0,30-location)])
# B=np.hstack([np.linspace(0,1,location),np.ones(30-location)])
cmap=np.vstack([R,G,B]).T
 
#####################################################################################
# plot_chi_heatmap: plots chi vs. concentration for all specified solutes for one IDR
#
# Arguments:
# - idr: the IDR
# - row: the row of the big heatmap where we plot this IDR
# - df: dataframe containing chi values to plot
# - repeats: number of experimental repeats done for this IDR
# - solutes: list of solutes -- each solute is one column of the heatmap
# - titles: column titles representing the solutes
#####################################################################################
def plot_chi_heatmap(idr, row, df, repeats, solutes, titles):
    # Choose a chi value on which to center the y axis
    center_chi = df[(df.idr==idr) & (df.sol=='Urea') & (df.conc==0)].chi.mean()
    for i in range(len(solutes)):
        # Calculate mean and standard deviation
        chi_series_list = []
        for rpt in range(1,repeats+1): 
            start_idx = df[(df.idr==idr) & (df.sol==solutes[i]) & (df.conc==0) & (df.repeat==rpt)].index.tolist()[0]
            Xs = df['conc'][start_idx:start_idx + 7].to_numpy()
            chi_series = df['chi'][start_idx:start_idx + 7].to_numpy()
            chi_series_list.append(chi_series)
        chis_mean = np.mean(chi_series_list, axis=0)
        chis_stdev = np.std(chi_series_list, axis=0)
        
        # Plot the cell for this IDR and solute
        col = solutes.index(solutes[i])
        axs[row][col].scatter(Xs,chis_mean,c='black',s=300,zorder=50)
        axs[row][col].errorbar(Xs,chis_mean,chis_stdev,c='gainsboro',zorder=100,fmt='none')
        if solutes[i] not in ['NaCl','KCl']:
            polyfit_x = Xs.astype(float)
            a=np.polyfit(polyfit_x,chis_mean,1)
            lloc = np.argmin(np.abs(cspace-a[0]*2*np.max(polyfit_x)));
            axs[row][col].set_facecolor(cmap[lloc,:])
        else:
            axs[row][col].set_facecolor('plum')

        # Set ticks
        top_conc = Xs[-1]
        margin = top_conc * 0.125
        left_lim = -margin
        right_lim = top_conc + margin
        axs[row][col].set_xlim([left_lim, right_lim])
        if top_conc==24:
            axs[row][col].set_xticks([0, 20])
            axs[row][col].set_xticklabels(['0', '20'])
        elif top_conc==12:
            axs[row][col].set_xticks([0, 10])
            axs[row][col].set_xticklabels(['0', '10'])
        elif top_conc==6:
            axs[row][col].set_xticks([0, 5])
            axs[row][col].set_xticklabels(['0', '5'])
        elif top_conc < 2:
            axs[row][col].set_xticks([0, 1])
            axs[row][col].set_xticklabels(['0', '1'])

        # Set y limits
        axs[row][col].set_ylim([center_chi-0.3, center_chi+0.3])

        if row == 0:
            axs[row][col].text(0,1.1,titles[i],horizontalalignment='left',
                               fontsize=60,rotation=30,transform=axs[row][col].transAxes)
        if col != 0:
            axs[row][col].yaxis.set_ticks_position('none')
        if col == 0:
            axs[row][col].set_ylabel('$\chi$', fontsize=60, rotation=0)
            axs[row][col].yaxis.set_label_coords(-0.8,0.38)
            axs[row][col].text(-1.1,.5,idr,verticalalignment='center',
                               fontsize=60,ha='right',transform=axs[row][col].transAxes)
    return

###################################################################################
# Solute and title lists (title list allows shorter names or nicknames for solutes)
###################################################################################
solutes = ['D-Sorbitol', 'Glycerol', 'Xylitol', 'Meso-Erythritol', 'D-Mannitol', 'Myo-Inositol',
           'EG', 'PEG200', 'PEG400', 'PEG1500', 'PEG2000', 'PEG4000', 'Ficoll', 
           'Glycine', 'L-Proline', 'Tricine', 'Sarcosine', 'BetaineH2O', 'L-Tryptophan',
           'Sucrose', 'L(+)-Arabinose', 'D-Galactose', 'Trehalose2H2O',
           'Urea', 'GuHCl', 'NaCl', 'KCl']
titles = ['Sorbitol', 'Glycerol', 'Xylitol', 'Erythritol', 'Mannitol', 'Inositol',
           'EG', 'PEG200', 'PEG400', 'PEG1500', 'PEG2000', 'PEG4000', 'Ficoll', 
           'Glycine', 'Proline', 'Tricine', 'Sarcosine', 'Betaine', 'Tryptophan',
           'Sucrose', 'Arabinose', 'Galactose', 'Trehalose',
           'Urea', 'GuHCl', 'NaCl', 'KCl']

########################################################################################################################
# Choose whether to make figure for naturally occurring IDRs (2A) or GS linkers (S1A) -- keep one, comment out the other
########################################################################################################################
figname = 'linkers' # Fig. S1A
# figname = 'idrs' # Fig. 2A

if figname == 'linkers':
    idrs = ['GS8','GS16','GS24','GS32','GS48']
elif figname == 'idrs':
    idrs = ['PUMA','Ash1','E1A','p53']

# Create and format the figure
ncols = len(solutes)
nrows = len(idrs)
wspace=0.4
hspace=0.4
fig, axs = plt.subplots(nrows=nrows, ncols=ncols, sharex=False, sharey=False, 
                        figsize=(ncols*2.5 + (wspace*(ncols-1)), nrows*2.5 + (hspace*(nrows-1))))
plt.subplots_adjust(wspace=0.04, hspace=0.04)  
fig.text(0.453, 0, 'Solute concentration (conc)', fontsize=60, ha='center')
fig.text(0.843, 0, 'Solute concentration (M)', fontsize=60, ha='center')

for row in range(nrows):
    axs[row][0].set_yticks([-0.2, 0, 0.2])
    axs[row][0].set_yticklabels(['-0.2', '0', '0.2'])
    axs[row][0].yaxis.tick_left()

for row in range(ncols):
    axs[nrows-1][row].xaxis.tick_bottom()
    
# Build data frame containing chi values for all IDRs in all solution conditions
all_idrs_chi_df = get_all_idrs_chi_df()

# Populate the cells in the figure
if figname == 'linkers':
    plot_chi_heatmap(idrs[0], 0, all_idrs_chi_df, 2, solutes, titles)
    plot_chi_heatmap(idrs[1], 1, all_idrs_chi_df, 2, solutes, titles)
    plot_chi_heatmap(idrs[2], 2, all_idrs_chi_df, 2, solutes, titles)
    plot_chi_heatmap(idrs[3], 3, all_idrs_chi_df, 2, solutes, titles)
    plot_chi_heatmap(idrs[4], 4, all_idrs_chi_df, 2, solutes, titles)
elif figname == 'idrs':
    plot_chi_heatmap(idrs[0], 0, all_idrs_chi_df, 2, solutes, titles)
    plot_chi_heatmap(idrs[1], 1, all_idrs_chi_df, 2, solutes, titles)
    plot_chi_heatmap(idrs[2], 2, all_idrs_chi_df, 2, solutes, titles)
    plot_chi_heatmap(idrs[3], 3, all_idrs_chi_df, 2, solutes, titles)

# More formatting
for row in range(nrows):
    for col in range (ncols):
        axs[row][col].label_outer()
        axs[row][col].tick_params(direction='in', length=12, width=5)
        axs[row][col].axhline(0, color='black', lw=4, alpha=1, linestyle='dashed')

# Save the figure
if figname == 'linkers':
    fig.savefig('hidden_sensitivity_fig_s1a.png', bbox_inches='tight')
elif figname == 'idrs':
    fig.savefig('hidden_sensitivity_fig_2a.png', bbox_inches='tight')

plt.show()
