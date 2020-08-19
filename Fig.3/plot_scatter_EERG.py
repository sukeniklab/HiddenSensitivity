# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 12:06:18 2020

@author: ssukenik
"""

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['svg.fonttype'] = 'none'
proteinlist=['E1A','Ash1','p53','puma_wild']
data = pd.read_csv('AData.csv')
#sns.set_palette("Set2",100)
#data=dataraw[(dataraw['Raw0/norm']>0.75)&(dataraw['Raw0/norm']<0.85)]
def plotScatter(data,name,metric,proteins='.*',colorBy='idx',legend=0,xlabel='',ylabel='',ylimit=''):
    '''
    input:
        metric: the features you want to plot (RG0, Kappa, etc.)
        proteins: A regexp string to select which proteins to plot. '.*' selects
        everything (default)
        colorBy: A feature that will determine colors (default is by index)
        legend: boolean definind if you're plotting a legend or not
    '''
    fix,ax=plt.subplots(1,len(metric),figsize=[3*len(metric),5])
    plt.tight_layout()
    B=data[data['Protein name'].str.match(proteins)].reset_index()
    B['idx']=B['Protein name']
    print(B)
    plt.figure(figsize=(20,20))
    # pkmn_type_colors = ['#78C850',  # Grass
    #                 '#F08030',  # Fire
    #                 '#6890F0',  # Water
    #                 '#A8B820',  # Bug
    #                 '#A8A878',  # Normal
    #                 '#A040A0',  # Poison
    #                 '#F8D030',  # Electric
    #                 '#E0C068',  # Ground
    #                 '#EE99AC',  # Fairy
    #                 '#C03028',  # Fighting
    #                 '#F85888',  # Psychic
    #                 '#B8A038',  # Rock
    #                 '#705898',  # Ghost
    #                 '#98D8D8',  # Ice
    #                 '#7038F8',  # Dragon
    #                ]
    order=0
    colortitle=['Attractive solution','Repulsive solution']
    colorref=['b','r']
    standard=['p53','E1A','puma_wildfull','Ash1']
    Highlight=data[data['Protein name'].isin(standard)]
    print(Highlight)
    for j in range(len(name)):
        for i in range(len(metric)):
           # sns.scatterplot(x=B[metric[i]], y=B[name[j]], data=B, color=colorref[order])
            plt.scatter(x=B[metric[i]], y=B[name[j]],s=B['Length'],color=colorref[order],alpha=0.25,linewidths=15,label=colortitle[order])
            if order==0:
                aaa='stdChi(3-0)'
            else:
                aaa='stdChi(-3-0)'
            plt.errorbar(x=B[metric[i]], y=B[name[j]],xerr=B['EEstdnorm'],yerr=B[aaa],color=colorref[order],alpha=0.15,fmt='none')
            plt.scatter(x=Highlight[metric[i]], y=Highlight[name[j]],marker='s',linewidths=20)
            ax=plt.gca()
            plt.xticks([-0.5,0,0.5,1,1.5],fontsize=40)
            ax.tick_params(direction='out', length=12, width=2)
            plt.yticks(fontsize=40)
            plt.xlabel(xlabel,fontsize=60)
            plt.ylabel(ylabel,fontsize=60)
            plt.ylim(-1.5,ylimit)
            plt.legend(fontsize=40)
            order+=1
        plt.savefig('abc.svg')
#    plt.vlines(B[metric],B[name[0]],B[name[1]], linestyles='dashdot')
    
# Plot metrics 3-11 only synthetic peptides by name match (4 capital letters)
#plotSwarm(data.columns[2:10],'[A-Z][A-Z][A-Z][A-Z]')
# Plot metrics 3-11 for only disprot proteins by name match (start with "DP")
#plotSwarm(data.columns[2:10],'DP.*')
# PLot only Kappa and RG_0 metrics for all proteins. color by length
plotScatter(data,['Chi(3-0)','Chi(-3-0)'],['ChiW'],'.*','Hydro',0,'$\u03C7$','$\Delta\u03C7$',1.5)
#plotScatter(data,['absMAXdif1'],['Raw0/norm'],'.*','Hydro',0,'$R_g$/$R^\Theta_g$','$\Delta (\Delta R_g$/$\Delta R_g^M)$',0.5)
#plotScatter(['(RawRGW-RawRGm3)/MAXdif'],['Raw0/norm'],'.*','Hydro',0,'$R_g$/$R^\Theta_g$','$\Delta R_g$/$\Delta R_g^M$')

# plotSwarm(['RawRG(0)','NormRG(0)','Helifraction(0)','Heliframe','HB','(NormRG(+1)-NRG(-1))/NRG(0)','(Heli(+1)-Heli(-1))/Heli(0)'],'.*','NormRG(0)')
# plotSwarm(['RawRG(0)','NormRG(0)','Helifraction(0)','Heliframe','HB','(NormRG(+1)-NRG(-1))/NRG(0)','(Heli(+1)-Heli(-1))/Heli(0)'],'.*','Helifraction(0)')
# plotSwarm(['RawRG(0)','NormRG(0)','Helifraction(0)','Heliframe','HB','(NormRG(+1)-NRG(-1))/NRG(0)','(Heli(+1)-Heli(-1))/Heli(0)'],'.*','Heliframe')
# plotSwarm(['RawRG(0)','NormRG(0)','Helifraction(0)','Heliframe','HB','(NormRG(+1)-NRG(-1))/NRG(0)','(Heli(+1)-Heli(-1))/Heli(0)'],'.*','HB')
# plotSwarm(['RawRG(0)','NormRG(0)','Helifraction(0)','Heliframe','HB','(NormRG(+1)-NRG(-1))/NRG(0)','(Heli(+1)-Heli(-1))/Heli(0)'],'.*','Hydrophobic')
# plotSwarm(['RawRG(0)','NormRG(0)','Helifraction(0)','Heliframe','HB','(NormRG(+1)-NRG(-1))/NRG(0)','(Heli(+1)-Heli(-1))/Heli(0)'],'.*','FCR')
# import seaborn as sns
# sns.set(style="whitegrid")
# tips = sns.load_dataset("tips")
# ax = sns.swarmplot(x=tips["total_bill"])
