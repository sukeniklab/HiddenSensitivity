# -*- coding: utf-8 -*-
"""
@author: Feng Yu
Here we get the data of GSlinker from SIdata.csv and us scipy to fit the curve.
"""

from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
x=[16,32,48,64,80,96,128]
x=np.array(x)
y=[1.985,2.915,3.597,4.115,4.651,4.683,5.700]
yerr=[0.006,0.026,0.043,0.163,0.194,0.456,0.315]#These data got from GS linker's simulation data.It is included in the Table_3.
squared_numbers = [number ** 6 for number in y]
def guess (x,k1,k2):
    return ((k2*x ** k1))
popt, pcov = curve_fit(guess, x, y)
print(popt)
print(popt[0])
perr = np.sqrt(np.diag(pcov))
print(perr)
plt.figure(figsize=(10,8))
plt.xticks(fontsize=30)
ax = plt.gca()
ax.tick_params(width=5,length=10)
plt.yticks(fontsize=30)
plt.ylabel('$R_e$ (nm)',fontsize=30)
plt.xlabel('$N$',fontsize=30)
xplot=np.arange(0,130,1)
plt.plot(xplot,0.55*(xplot)**0.48,linewidth=5)
plt.text(50,2,r'$R_e =0.55 \times N^{0.48}$',fontsize=30)
plt.scatter(x, y,linewidths=20,color='y',marker='.')
plt.errorbar(x, y, yerr=yerr,fmt='none',linewidth=3,barsabove=True,color='black')
plt.savefig('S2.svg')
