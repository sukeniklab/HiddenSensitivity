##########################################################################
# Figure 1A: FRET experiments
#
# Supporting material for "Probing the Hidden Sensitivity of Intrinsically 
# Disordered Proteins to their Chemical Environment"
# https://www.biorxiv.org/content/10.1101/2020.08.17.252478v1
#
# Contact: David Moses dmoses5@ucmerced.edu
##########################################################################

from hidden_sensitivity import preprocess_tq_ng_data
from hidden_sensitivity import build_base_correction_factor_df
from hidden_sensitivity import get_efret

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.offsetbox import (TextArea, DrawingArea, OffsetImage, AnnotationBbox)
from matplotlib.cbook import get_sample_data
import matplotlib.image as mpimg
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import math
import seaborn as sns
sns.set_style("white")
plt.rc('xtick', labelsize=24) 
plt.rc('xtick.major', size=8, pad=7)
plt.rc('ytick', labelsize=24) 
plt.rc('ytick.major', size=8, pad=7)
plt.rcParams["axes.labelsize"] = 24
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'

# Read in and preprocess raw data for mTurquoise2 and mNeonGreen base spectra
TQ_data, NG_data = preprocess_tq_ng_data()

# Read in data for one IDR
GS16_data_2=pd.read_csv("GS16-1-repeat-2.csv",sep=",",skiprows=1).join(pd.read_csv('GS16-2-repeat-2.csv',
                    skiprows=1).iloc[:,3:],rsuffix='2').join(pd.read_csv('GS16-3-repeat-2.csv',
                    skiprows=1).iloc[:,3:],rsuffix='3')
                    
# Build a data frame with correction factors for the donor and acceptor at each concentration of each solute
base_corr_fact_df=build_base_correction_factor_df(TQ_data, NG_data)
                    
# Build an E_fret data frame for that IDR
GS16_efret_2=get_efret('GS16',GS16_data_2,range(3,GS16_data_2.shape[1]),base_corr_fact_df,TQ_data,NG_data,True)

# Get example spectra in buffer, in an expanding solution, and in a compacting solution
df = pd.DataFrame({ "wavelengths":GS16_efret_2[GS16_efret_2.sol=='Buffer']['x_range'].item(),
                    "f_d":GS16_efret_2[GS16_efret_2.sol=='Buffer']['f_d'].item(),
                    "f_a":GS16_efret_2[GS16_efret_2.sol=='Buffer']['f_a'].item(),
                    "fit_buffer":GS16_efret_2[GS16_efret_2.sol=='Buffer']['fit'].item(),
                    "fit_expanding":GS16_efret_2[GS16_efret_2.sol=='Tricine.4']['fit'].item(),
                    "fit_compacting":GS16_efret_2[GS16_efret_2.sol=='PEG2000.4']['fit'].item() })

fig, ax = plt.subplots(1, figsize=(8,8))

plt.plot(df.wavelengths, 100 * df.fit_buffer / df.fit_buffer[26], c='k', label='buffer')
plt.plot(df.wavelengths, 100 * df.fit_expanding /  df.fit_expanding[26], c='b', label='expanding solution')
plt.plot(df.wavelengths, 100 * df.fit_compacting /  df.fit_compacting[26], c='r', label='compacting solution')
plt.plot(df.wavelengths, 5000 * df.f_d / df.f_d.sum(), c='cyan', label = 'FRET donor')
plt.plot(df.wavelengths, 3000 * df.f_a / df.f_a.sum(), c='lime', label = 'FRET acceptor')
plt.fill_between(df.wavelengths, 5000 * df.f_d / df.f_d.sum(), color='cyan', alpha=0.5)
plt.fill_between(df.wavelengths, 3000 * df.f_a / df.f_a.sum(), color='lime', alpha=0.3)

img_expand = mpimg.imread("1a1.png")
img_compact = mpimg.imread("1a2.png")
img_buffer = mpimg.imread("1a3.png")
imagebox_expand = OffsetImage(img_expand, zoom=0.059)
imagebox_compact = OffsetImage(img_compact, zoom=0.055)
imagebox_buffer = OffsetImage(img_buffer, zoom=0.053)
ab_expand = AnnotationBbox(imagebox_expand, xy=(570, 112), frameon = False)
ab_buffer = AnnotationBbox(imagebox_buffer, (568, 152), frameon = False)
ab_compact = AnnotationBbox(imagebox_compact, (570, 195), frameon = False)
ax.add_artist(ab_expand)
ax.add_artist(ab_compact)
ax.add_artist(ab_buffer)

plt.xlabel('Wavelength (nm)', fontsize=24)
plt.ylabel('Fluorescence (AU)', fontsize=24)
ax.legend(fontsize=20, loc='center left', bbox_to_anchor=(1, 0.8))
ax.yaxis.tick_left()
ax.xaxis.tick_bottom()
ax.set_xlim([450, 600])

plt.savefig('hidden_sensitivity_fig_1a.png', bbox_inches='tight')
plt.show()
