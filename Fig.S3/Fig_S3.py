# -*- coding: utf-8 -*-
"""
@author: Feng Yu
Here we load the data from csv and groupby protein name.
We sort the values based on MTFE (we called psi in the paper)
And then make plotting
"""
import pandas as pd
import matplotlib.pyplot as plt
def pandatolist(data,MTFE,term):
    newname=data[data['MTFE']==MTFE][term]
    newname=newname.tolist()
    return(newname[0])  
data = pd.read_csv('HEATDATA_2.csv') # HEATDATA_1 is for S3A and HEATDATA_2 is for S3B
dataname = pd.read_csv('SIData.csv') 
C=data.groupby('protein')
dash_x=[-0.4,3.4]
dash_y=[0,0]
plt.tight_layout()
order=1
plt.figure(figsize=(40,60))
plt.subplots_adjust(left=0, bottom=0.2, right=0.98, top=0.8, wspace=0, hspace=0.15)
for index,j in data.groupby('protein'):
    plt.subplot(7,5,order)
    print(j.size) # We have two different data size. One kind of data have 7 data points. Another one have 9.
    if j.size==45 :
        jt=j.sort_values(["MTFE"],ascending=True).tail(5)
        jh=j.sort_values(["MTFE"],ascending=True).head(5)
    else:
        jt=j.sort_values(["MTFE"],ascending=True).tail(4)
        jh=j.sort_values(["MTFE"],ascending=True).head(4)
    plt.scatter(x=jt['MTFE'], y=jt['Chi'],color='blue',label=index,linewidths=20)
    plt.scatter(x=-jh['MTFE'], y=jh['Chi'],color='red',label=index,linewidths=20)
    plt.plot(dash_x,dash_y, linestyle='dashed',color='black',dashes=(5, 20))
    plt.errorbar(x=jt['MTFE'], y=jt['Chi'],yerr=jt['stdChi'],fmt='none')
    plt.errorbar(x=-jh['MTFE'], y=jh['Chi'],yerr=jh['stdChi'],fmt='none')
    plt.ylim(-1,2.5)
    plt.xlim(-0.4,3.4)
    uniprot=dataname[dataname['simulation name']==index]['Uniprot names']
    uniprottext=uniprot.tolist()[0]
    plt.text(0,2,uniprottext,fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks([-1,0,1,2],fontsize=30)
    if order>=29:                       #Adjust x axis
        plt.xlabel('$\u03C8$',fontsize=40)
    order+=1
    ax = plt.gca()
    orderlist=[2,7,12,17,22,27,32] # Adjust y axis
    if order in orderlist :
        plt.ylabel('$\u03C7$',fontsize=40)
    else:
        ax.get_yaxis().set_visible(False)
plt.savefig('S3B.svg')
