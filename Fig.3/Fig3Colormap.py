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
def colormap(orderset,numberofsubplot,column,symbol,colormapnumber,cspacename,specialcmap=0,convertfunction=0):
    order=orderset
    for j in range(0, len(data)):
        plt.subplot(9,70,order)
        plt.ylim(0,10)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=30,rotation=90)
        ax = plt.gca()
        colorindexraw=data.iloc[j][column]
        if convertfunction==0:
            lloc = np.argmin(np.abs(cspacename-colorindexraw))
            ax.set_facecolor(colormapnumber(lloc))
        elif specialcmap==1:
            colorindex=convertfunction(column,colorindexraw)
            lloc = np.argmin(np.abs(cspacename-colorindex))
            ax.set_facecolor(cmap[lloc,:])
        elif convertfunction==convertoverwrite:
            colorindex=convertfunction(column,colorindexraw,100)
            lloc = np.argmin(np.abs(cspacename-colorindex))
            ax.set_facecolor(colormapnumber(lloc))
        else:
            colorindex=convertfunction(column,colorindexraw)
            lloc = np.argmin(np.abs(cspacename-colorindex))
            ax.set_facecolor(colormapnumber(lloc))
        if order==orderset:
            ax.set_ylabel(symbol,fontsize=25,rotation=0)
            ax.set_yticklabels([])
        else:
            ax.get_yaxis().set_visible(False)
        ax.get_xaxis().set_visible(False)
        order+=1
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
    return(-1*output)
def pandatolist(data,MTFE,term):
    newname=data[data['MTFE']==MTFE][term]
    newname=newname.tolist()
    return(newname[0])  
#Colormap setup
minll=-1.35
maxll=1.35
minl=-1
maxl=1
cspacenew2=np.linspace(minl,maxl,256)
cspacenew3=np.linspace(minl,maxl,256)
cspacenew=np.linspace(minll,maxll,256)
cspace=np.linspace(minl,maxl,30)
location = np.argmin(np.abs(cspace))
R=np.hstack([np.ones(location), np.linspace(1,0,30-location)])
G=np.hstack([np.linspace(0,1,location),np.linspace(1,0,30-location)])
B=np.hstack([np.linspace(0,1,location),np.ones(30-location)])
cmap=np.vstack([R,G,B]).T
c = Colormap()
mycmap2= c.cmap_linear('purple(w3c)', 'white', 'orange')
mycmap4=c.cmap_linear('red','white','blue')
dataraw = pd.read_csv('ADATA_short.csv')
data = dataraw.sort_values('ChiW')
plt.tight_layout()
plt.figure(figsize=(9,5))
plt.subplots_adjust(left=0, bottom=0.2, right=0.98, top=0.8, wspace=0, hspace=0.15)
nofsubplot=70
count=1
colormap(count,nofsubplot,'ChiW','\u03C7           ',mycmap2,cspacenew)
count+=nofsubplot
colormap(count,nofsubplot,'Chi(3-0)','\u0394\u03C7$_{exp}$           ',cmap,cspace,1,convert)
count+=nofsubplot
colormap(count,nofsubplot,'Chi(-3-0)','\u0394\u03C7$_{com}$           ',cmap,cspace,1,convertminus)
count+=nofsubplot
colormap(count,nofsubplot,'Length','N           ',mycmap2,cspacenew2,0,convertoverwrite)
count+=nofsubplot
colormap(count,nofsubplot,'Kappa','\u03BA           ',mycmap2,cspacenew2,0,convert)
count+=nofsubplot
colormap(count,nofsubplot,'Omega','\u03A9           ',mycmap2,cspacenew2,0)
count+=nofsubplot
colormap(count,nofsubplot,'NCPR','NCPR           ',mycmap2,cspacenew2,0)
count+=nofsubplot
colormap(count,nofsubplot,'FCR','FCR           ',mycmap2,cspacenew2,0)
plt.savefig('S3Colormap.svg')
