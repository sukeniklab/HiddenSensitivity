
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['svg.fonttype'] = 'none'

data = pd.read_csv('AData_short.csv')
def plotScatter(data,name,metric,xlabel='',ylabel='',ylimit=''):
    plt.figure(figsize=(20,20))
    order=0
    colortitle=['Attractive solution','Repulsive solution']
    colorref=['b','r']
    standard=['p53','E1A','puma_wildfull','Ash1']
    Highlight=data[data['Protein name'].isin(standard)]
    print(Highlight)
    for j in range(len(name)):
        for i in range(len(metric)):
            plt.scatter(x=data[metric[i]], y=data[name[j]],s=data['Length'],color=colorref[order],alpha=0.25,linewidths=15,label=colortitle[order])
            if order==0:
                err='stdChi(3-0)'
            else:
                err='stdChi(-3-0)'
            plt.errorbar(x=data[metric[i]], y=data[name[j]],xerr=data['stdChi(0)'],yerr=data[err],color=colorref[order],alpha=0.15,fmt='none')
            plt.scatter(x=Highlight[metric[i]], y=Highlight[name[j]],marker='s',linewidths=20) #highlight the protein in the experiment
            ax=plt.gca()
            plt.xticks([-0.5,0,0.5,1,1.5],fontsize=40)
            ax.tick_params(direction='out', length=12, width=2)
            plt.yticks(fontsize=40)
            plt.xlabel(xlabel,fontsize=60)
            plt.ylabel(ylabel,fontsize=60)
            plt.ylim(-1.5,ylimit)
            plt.legend(fontsize=40)
            order+=1
        plt.savefig('S3scatter.svg')
plotScatter(data,['Chi(3-0)','Chi(-3-0)'],['ChiW'],'$\u03C7$','$\Delta\u03C7$',1.5)
