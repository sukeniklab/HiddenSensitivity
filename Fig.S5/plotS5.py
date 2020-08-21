# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 13:59:08 2020

@author: ssukenik
"""

import pandas as pd
import matplotlib.pyplot as plt
import math

allN=[20, 30, 40, 60, 80, 100] #all lengths
AllData=pd.DataFrame()
selPsi=150
 #The psi value to plot
R='RE' #Can be RE for end-to-end or RG for gyration
#The for loop loops over all lengths
fig,ax=plt.subplots(2,3,figsize=[12,12])
i=0
for N in allN:
    
    Rdata = pd.read_csv(R+'_'+str(N)+'.csv')
    Rdata.rename(columns={'Unnamed: 0':'psi'},inplace=True)
    Rdata.iloc[:,0]=[-300, -200, -150, -100, -50, 0, 50, 100, 150, 200, 300]
    Rdata=Rdata.T
    Rdata.reset_index(inplace=True)
    Rdata.columns=Rdata.loc[0]
    Rdata.drop(0,inplace=True)
    Rdata.rename(columns={'psi':'seq'}, inplace=True)
    Rdata.columns.name=''
    
    print(str(math.floor(i/3))+","+str(i%3))
    ax[math.floor(i/3),i%3].plot(Rdata.T.iloc[2:-1,:],c='black',linewidth=0.1)
    ax[math.floor(i/3),i%3].text(0,2,"N="+str(N))
    i+=1
    plt.savefig('Fig.S5.svg')