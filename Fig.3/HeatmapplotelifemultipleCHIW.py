# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 00:49:24 2020
重要 DP02138 OMEGA 值改爲0！！！！！！！！！！
@author: ShaharGroup-fyu
"""
from colormap import Colormap
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['svg.fonttype'] = 'none'
def plot_examples(colormaps,a,b):
    """
    Helper function to plot data with associated colormap.
    """
    np.random.seed(19680801)
    data = np.random.randn(30, 30)
    n = len(colormaps)
    fig, axs = plt.subplots(1, n, figsize=(n * 2 + 2, 3),
                            constrained_layout=True, squeeze=False)
    for [ax, cmap] in zip(axs.flat, colormaps):
        psm = ax.pcolormesh(data, cmap=cmap, rasterized=True, vmin=a, vmax=b)
        fig.colorbar(psm, ax=ax,orientation=0)
    plt.savefig(str(b)+'.svg')
    plt.show()
def convert(columnname,rawvalue):
    maxvalue=data.max()[columnname]
    minvalue=data.min()[columnname]
    output=(rawvalue-minvalue)/(maxvalue-minvalue)
    return(output)
def convertoverwrite(columnname,rawvalue,Overwrite):
    maxvalue=Overwrite
    minvalue=data.min()[columnname]
    output=(rawvalue-minvalue)/(maxvalue-minvalue)
    return(output)
def convertminus(columnname,rawvalue):
    maxvalue=data.max()[columnname]
    minvalue=data.min()[columnname]
    output=(maxvalue-rawvalue)/(maxvalue-minvalue)
    return(output)
def pandatolist(data,MTFE,term):
    newname=data[data['MTFE']==MTFE][term]
    newname=newname.tolist()
    return(newname[0])  
minll=-1.35
maxll=1.35
minl=-1
maxl=1
cspacenew2=np.linspace(minl,maxl,256)
cspacenew3=np.linspace(minl,maxl,256)
cspacenew=np.linspace(minll,maxll,256)
cspace=np.linspace(minl,maxl,30)
location = np.argmin(np.abs(cspace))
#blue-white-red
R=np.hstack([np.ones(location), np.linspace(1,0,30-location)])
G=np.hstack([np.linspace(0,1,location),np.linspace(1,0,30-location)])
B=np.hstack([np.linspace(0,1,location),np.ones(30-location)])
RR=np.hstack([np.ones(location), np.linspace(1,0,30-location)])
GG=np.hstack([np.linspace(0,1,location),np.linspace(0,1,30-location)])
BB=np.hstack([np.linspace(0,1,location),np.ones(30-location)])
#green-white-magenta
# R=np.hstack([np.linspace(0,1,location),np.ones(30-location)])
# G=np.hstack([np.ones(location),np.linspace(1,0,30-location)])
# B=np.hstack([np.linspace(0,1,location),np.ones(30-location)])
cmap=np.vstack([R,G,B]).T
cmaptemp=np.vstack([GG,B,R]).T
c = Colormap()
mycmap = c.cmap_linear('purple(w3c)', 'white', 'green')
mycmap2= c.cmap_linear('purple(w3c)', 'white', 'orange')
mycmap3= c.cmap_linear('coral', 'white', 'deeppink')
mycmap4=c.cmap_linear('red','white','blue')
plot_examples([mycmap],-0.56641,1.33256)
plot_examples([mycmap2],20,158)
plot_examples([mycmap3],0.0372871,0.98424)
plot_examples([mycmap4],-1,1)
dataraw = pd.read_csv('ADATA.csv')
data = dataraw.sort_values('ChiW')
plt.tight_layout()
order=1
plt.figure(figsize=(12,5))
plt.subplots_adjust(left=0, bottom=0.2, right=0.98, top=0.8, wspace=0, hspace=0.15)
for j in range(0, len(data)):
    plt.register_cmap(cmap=mycmap)
    plt.subplot(9,95,order)
    plt.ylim(0,10)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=30,rotation=90)
    colorindexraw=data.iloc[j]['ChiW']
    lloc = np.argmin(np.abs(cspacenew-colorindexraw))
    ax = plt.gca()
    if order==1:
        ax.set_ylabel('\u03C7           ',fontsize=25,rotation=0)
        ax.set_yticklabels([])
    else:
        ax.get_yaxis().set_visible(False)
    ax.set_facecolor(mycmap2(lloc))
    ax.get_xaxis().set_visible(False)
    order+=1
order=96
for j in range(0, len(data)):
    plt.subplot(9,95,order)
    plt.ylim(0,10)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30,rotation=90)
    colorindexraw=data.iloc[j]['Chi(3-0)']
    colorindex=convert('Chi(3-0)',colorindexraw)
    lloc = np.argmin(np.abs(cspace-colorindex))
    ax = plt.gca()
    if order==96:
        ax.set_ylabel('\u0394\u03C7$_{exp}$           ',fontsize=25,rotation=0)
        ax.set_yticklabels([])
    else:
        ax.get_yaxis().set_visible(False)
    ax.set_facecolor(cmap[lloc,:])
    ax.get_xaxis().set_visible(False)
    plt.xlabel('\u03c8%',fontsize=30)
    order+=1
order=191
for j in range(0, len(data)):
    plt.subplot(9,95,order)
    plt.ylim(0,10)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=30)
    print(data.iloc[j]['Protein name'])
    colorindexraw=data.iloc[j]['Chi(-3-0)']
    colorindex=-1*convertminus('Chi(-3-0)',colorindexraw)
    lloc = np.argmin(np.abs(cspace-colorindex))
    ax = plt.gca()
    if order==191:
        ax.set_ylabel('\u0394\u03C7$_{com}$           ',fontsize=25,rotation=0)
        ax.set_yticklabels([])
    else:
        ax.get_yaxis().set_visible(False)
    ax.set_facecolor(cmap[lloc,:])
    ax.get_xaxis().set_visible(False)
    order+=1
order=286
for j in range(0, len(data)):
    plt.subplot(9,95,order)
    plt.ylim(0,10)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=30,rotation=90)
    colorindexraw=data.iloc[j]['Length']
    colorindex=convertoverwrite('Length',colorindexraw,100)
    lloc = np.argmin(np.abs(cspacenew2-colorindex))
    ax = plt.gca()
    if order==286:
        ax.set_yticklabels([])
        ax.set_ylabel('N           ',fontsize=25,rotation=0)
    else:
        ax.get_yaxis().set_visible(False)
    ax.set_facecolor(mycmap2(lloc))
    ax.get_xaxis().set_visible(False)
    order+=1
order=381
for j in range(0, len(data)):
    plt.subplot(9,95,order)
    plt.ylim(0,10)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=30,rotation=90)
    colorindexraw=data.iloc[j]['Kappa']
    colorindex=convert('Kappa',colorindexraw)
    lloc = np.argmin(np.abs(cspacenew2-colorindex))
    ax = plt.gca()
    fig = plt.gcf()
    if order==381:
        ax.set_ylabel('\u03BA           ',fontsize=25,rotation=0)
        ax.set_yticklabels([])
    else:
        ax.get_yaxis().set_visible(False)
    ax.set_facecolor(mycmap2(lloc))
    ax.get_xaxis().set_visible(False)
    order+=1
order=476
for j in range(0, len(data)):
    plt.subplot(9,95,order)
    plt.ylim(0,10)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=30,rotation=90)
    colorindexraw=data.iloc[j]['Omega']
    lloc = np.argmin(np.abs(cspacenew2-colorindexraw))
    ax = plt.gca()
    fig = plt.gcf()
    if order==476:
        ax.set_ylabel('\u03A9           ',fontsize=25,rotation=0)
        ax.set_yticklabels([])
    else:
        ax.get_yaxis().set_visible(False)
    ax.set_facecolor(mycmap2(lloc))
    ax.get_xaxis().set_visible(False)
    order+=1
order=571
for j in range(0, len(data)):
    plt.subplot(9,95,order)
    plt.ylim(0,10)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=30,rotation=90)
    colorindexraw=data.iloc[j]['NCPR']
    lloc = np.argmin(np.abs(cspacenew2-colorindexraw))
    ax = plt.gca()
    fig = plt.gcf()
    if order==571:
        ax.set_ylabel('NCPR          ',fontsize=25,rotation=0)
        ax.set_yticklabels([])
    else:
        ax.get_yaxis().set_visible(False)
    ax.set_facecolor(mycmap2(lloc))
    ax.get_xaxis().set_visible(False)
    order+=1
order=666
for j in range(0, len(data)):
    plt.subplot(9,95,order)
    plt.ylim(0,10)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=30,rotation=90)
    colorindexraw=data.iloc[j]['FCR']
    lloc = np.argmin(np.abs(cspacenew2-colorindexraw))
    ax = plt.gca()
    if order==666:
        ax.set_yticklabels([])
        ax.set_ylabel('FCR          ',fontsize=25,rotation=0)
    else:
        ax.get_yaxis().set_visible(False)
    ax.set_facecolor(mycmap2(lloc))
    ax.get_xaxis().set_visible(False)
    order+=1
plt.savefig('abc.svg')
# order=761
# for j in range(0, len(data)):
#     plt.subplot(9,95,order)
#     plt.ylim(0,10)
#     plt.xticks(fontsize=20)
#     plt.yticks(fontsize=30,rotation=90)
#     colorindexraw=data.iloc[j]['Hydro']
#     lloc = np.argmin(np.abs(cspacenew2-colorindexraw))
#     ax = plt.gca()
#     if order==761:
#         ax.set_yticklabels([])
#         ax.set_ylabel('FHR           ',fontsize=25,rotation=0)
#     else:
#         ax.get_yaxis().set_visible(False)
#     ax.set_facecolor(mycmap2(lloc))
#     ax.get_xaxis().set_visible(False)
#     order+=1
# plt.show()
# plt.figure(figsize=(10,1))
# plt.subplots_adjust(left=0, bottom=0.2, right=1, top=0.8, wspace=0, hspace=0)
# order=1
# orderd=[14,17,21,25,29]
# for kk in orderd:
#     plt.subplot(1,5,order)
#     ax = plt.gca()
#     ax.set_facecolor(cmap[kk])
#     ax.get_yaxis().set_visible(False)
#     ax.get_xaxis().set_visible(False)
#     order+=1