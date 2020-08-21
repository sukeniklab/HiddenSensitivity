##############################################################################
# Figure S1B: Chi vs. concentration for several IDRs focusing on a few solutes
#
# Supporting material for "Probing the Hidden Sensitivity of Intrinsically 
# Disordered Proteins to their Chemical Environment"
# https://www.biorxiv.org/content/10.1101/2020.08.17.252478v1
#
# Contact: David Moses dmoses5@ucmerced.edu
##############################################################################

from hidden_sensitivity import get_all_idrs_chi_df

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import matplotlib.cm as cm
from scipy import stats
from matplotlib.ticker import FormatStrFormatter
plt.rc('xtick', labelsize=28) 
plt.rc('ytick', labelsize=28) 
plt.rcParams["xtick.major.size"] = 8
plt.rcParams["ytick.major.size"] = 8
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
sns.set_style("white")

################################################################################################################
# plot_delta_chi_series_n_solutes_one_idr: plots chi vs. concentration for one IDR for specified solutes
#
# Arguments:
# - idr: the IDR
# - df: dataframe containing chi values to plot
# - repeats: number of experimental repeats done for this IDR
# - solutes: list of solutes, each of which is plotted on a different axis according to its position in the list
# - color: the color of the plotted points for the IDR
################################################################################################################
def plot_delta_chi_series_n_solutes_one_idr(idr, df, repeats, solutes, color):
    for i in range(len(solutes)):
        # Calculate mean and standard deviation
        delta_chi_series_list = []
        for rpt in range(1,repeats+1): 
            start_idx = df[(df.idr==idr) & (df.sol==solutes[i]) & (df.conc==0) & (df.repeat==rpt)].index.tolist()[0]
            Xs = df['conc'][start_idx:start_idx + 7].to_numpy()
            chi_series = df['chi'][start_idx:start_idx + 7].to_numpy()
            delta_chi_series = chi_series - chi_series[0]
            delta_chi_series_list.append(delta_chi_series)
        delta_chis_mean = np.mean(delta_chi_series_list, axis=0)
        delta_chis_stdev = np.std(delta_chi_series_list, axis=0)
        axs[i].scatter(Xs,delta_chis_mean,color=color,s=300,label=idr,edgecolors='black')
        axs[i].errorbar(Xs,delta_chis_mean,delta_chis_stdev,color='black',fmt='none')
        axs[i].text(0.5,1.05,solutes[i],horizontalalignment='center',fontsize=30,transform=axs[i].transAxes)
        axs[i].legend()
        if i != 0:
            axs[i].yaxis.set_ticks_position('none')
        if i == 0:
            axs[i].set_ylabel('$\Delta$$\chi$', fontsize=32)
        # Set ticks
        top_conc = Xs[-1]
        margin = top_conc * 0.125
        left_lim = -margin
        right_lim = top_conc + margin
        axs[i].set_xticks([0, top_conc/3, 2*top_conc/3, top_conc])
        if top_conc == 24:
            axs[i].set_xticklabels(['0', '8', '16', '24'])
        else:
            axs[i].set_xticklabels(['0', str(top_conc/3), str(2*top_conc/3), str(top_conc)])
        axs[i].set_xlim([left_lim, right_lim])

    return

################
# Choose solutes
################
solutes = ['PEG200', 'PEG2000', 'Ficoll', 'Sarcosine', 'Tricine', 'Urea', 'NaCl']

########################################################################################################################
# Choose whether to make figure for naturally occurring IDRs (2B) or GS linkers (S1B) -- keep one, comment out the other
########################################################################################################################
figname = 'linkers' # Fig. S1B
# figname = 'idrs' # Fig. 2B

# Create and format the figure
ncols = len(solutes)
nrows = 1 
fig, axs = plt.subplots(nrows=nrows, ncols=ncols, sharex=False, sharey=True, figsize=(5*ncols + 1.25, 5))
plt.subplots_adjust(wspace=0.05, hspace=0.05)  
fig.text(0.381, -0.03, 'Solute concentration (conc)', fontsize=30, ha='center')
fig.text(0.773, -0.03, 'Solute concentration (M)', fontsize=30, ha='center')
for i in range(ncols):
    axs[i].xaxis.tick_bottom()
axs[0].yaxis.tick_left()

# Build data frame containing chi values for all IDRs in all solution conditions
all_idrs_chi_df = get_all_idrs_chi_df()

# Populate the cells
if figname == 'linkers':
    plot_delta_chi_series_n_solutes_one_idr('GS8', all_idrs_chi_df, 2, solutes, 'black')
    plot_delta_chi_series_n_solutes_one_idr('GS16', all_idrs_chi_df, 2, solutes, 'dimgray')
    plot_delta_chi_series_n_solutes_one_idr('GS24', all_idrs_chi_df, 2, solutes, 'darkgray')
    plot_delta_chi_series_n_solutes_one_idr('GS32', all_idrs_chi_df, 2, solutes, 'silver')
    plot_delta_chi_series_n_solutes_one_idr('GS48', all_idrs_chi_df, 2, solutes, 'lightgray')
elif figname == 'idrs':
    plot_delta_chi_series_n_solutes_one_idr('GS24', all_idrs_chi_df, 2, solutes, 'darkgray')
    plot_delta_chi_series_n_solutes_one_idr('PUMA', all_idrs_chi_df, 2, solutes, 'tab:cyan')
    plot_delta_chi_series_n_solutes_one_idr('Ash1', all_idrs_chi_df, 2, solutes, 'tab:pink')
    plot_delta_chi_series_n_solutes_one_idr('E1A', all_idrs_chi_df, 2, solutes, 'tab:orange')
    plot_delta_chi_series_n_solutes_one_idr('p53', all_idrs_chi_df, 2, solutes, 'limegreen')

for col in range (ncols):
    axs[col].grid(False)
    axs[col].label_outer()
    axs[col].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    axs[col].axhline(0, color='black', lw=3, alpha=1, linestyle='dashed')
    axs[col].tick_params(direction='in', length=8, width=3)

for i in range(ncols):
    axs[i].get_legend().remove()
axs[4].legend(fontsize=14, loc='lower right')

# Save the figure
if figname == 'linkers':
    fig.savefig('hidden_sensitivity_fig_s1b.png', bbox_inches='tight')
elif figname == 'idrs':
    fig.savefig('hidden_sensitivity_fig_2b.png', bbox_inches='tight')

plt.show()
