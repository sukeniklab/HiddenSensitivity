# -*- coding: utf-8 -*-
"""
This script plots the weak interaction figure. 
To plot other interaction conditions, change the Chi(1-0) Chi(-1-0) and stdChi(1-0) stdChi(-1-0) to corresponding data
"""

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
data = pd.read_csv('DeltaChiData.csv')
def plotScatter(data,name,metric,proteins='.*',colorBy='idx',legend=0,xlabel='',ylabel='',ylimit=''):
    fix,ax=plt.subplots(1,len(metric),figsize=[3*len(metric),5])
    plt.figure(figsize=(20,20))
    order=0
    colorref=['black','r']
    plt.figure(figsize=(20,20))
    order=0
    colorref=['b','r']
    for j in range(len(name)):
        for i in range(len(metric)):
            plt.scatter(x=data[metric[i]], y=data[name[j]],color=colorref[order],alpha=0.25,linewidths=15)
            if order==0:
                aaa='stdChi(1-0)'
            else:
                aaa='stdChi(-1-0)'
            plt.errorbar(x=data[metric[i]], y=data[name[j]],xerr=data['stdChi(0)'],yerr=data[aaa],color=colorref[order],alpha=0.5,linewidths=10,fmt='none')
            ax=plt.gca()
            plt.xticks([-0.5,0,0.5,1,1.5],fontsize=40)
            ax.tick_params(direction='out', length=12, width=2)
            plt.xticks(fontsize=40)
            plt.text(0.3,1,'$Weak$ $Interaction$',fontsize=40)
            plt.yticks(fontsize=40)
            plt.xlabel(xlabel,fontsize=60)
            plt.ylabel(ylabel,fontsize=60)
            plt.ylim(-1.5,ylimit)
            order+=1
    plt.savefig('S4.svg')
plotScatter(data,['Chi(1-0)','Chi(-1-0)'],['Chi(0)'],'.*','Hydro',0,'$\u03C7$','$\Delta\u03C7$',1.5)

