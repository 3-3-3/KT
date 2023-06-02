#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 15:52:14 2022

@author: ppurdue
"""

import numpy as np
import matplotlib.pylab as plt        # need for labels and legend
import pandas as pd                   # need for scatter log-log plot
import matplotlib.lines as mlines     # need for creating legend
import matplotlib.patches as mpatches  # need for fill patches in legend

# ---------------  Reading in data & making scatter plot ---------------------

# Read the data from a file:
df1 = pd.read_csv('Taumaline/full_DWD_1.txt', sep=' ', skiprows=11)

# Assign colors based on type of data point and SNR:
btp = df1.loc[:,'WDbintype']
nummax = len(btp)

clr=[]
for jj in range(0,nummax):
    if btp[jj] == 1:
        clr.append('red')
#        mkr.append('x')
    elif btp[jj] == 2:
        clr.append('orange')
#        mkr.append('*')
    elif btp[jj] == 3:
        clr.append('yellow')
#        mkr.append('o')
    elif btp[jj] == 4:
        clr.append('green')
#        mkr.append('1')
    elif btp[jj] == 5:
        clr.append('blue')
#        mkr.append('3')
    else:
        clr.append('magenta')
#        mkr.append('s')

# Create a scatter plot of the data.
fig = df1.plot.scatter(x='Freq', y='dhnorm',loglog=True,marker='8',s=.001,c=clr, zorder=2)

# --------------- Setting the x-values for line sections ------------------

xmin = 5e-4
xmax = 1

jcn1 = .0015
jcn2 = .027
jcn3 = .04


# --------------- Plotting the lines on the graph -------------------------


# 100 linearly spaced numbers (full range)
xfull = np.geomspace(xmin,xmax,100)

# needed subsets
x12 = np.geomspace(jcn1,jcn2,50)
x23 = np.geomspace(jcn2,jcn3,50)
xminto3 = np.geomspace(xmin,jcn3,50)
x3tomax = np.geomspace(jcn3,xmax,50)

# boundary A
bA = pow(151*pow(xfull,2),0.33333)
fig.loglog(xfull,bA,color='black',linewidth=.5)

# boundary B
bBfull = 10**(0.703+0.637*np.log10(xfull)-0.017*pow(np.log10(xfull),2)
          +0.298*pow(np.log10(xfull),3)+0.061*pow(np.log10(xfull),4))
fig.loglog(xfull,bBfull, 'black',linewidth=.5,linestyle='--')

bB3tomax = 10**(0.703+0.637*np.log10(x3tomax)-0.017*pow(np.log10(x3tomax),2)
          +0.298*pow(np.log10(x3tomax),3)+0.061*pow(np.log10(x3tomax),4))

bBminto3 = 10**(0.703+0.637*np.log10(xminto3)-0.017*pow(np.log10(xminto3),2)
          +0.298*pow(np.log10(xminto3),3)+0.061*pow(np.log10(xminto3),4))

# boundary C
bCfull = 10**(0.761+1.005*np.log10(xfull)+0.7*pow(np.log10(xfull),2)
          +0.7*pow(np.log10(xfull),3)+0.214*pow(np.log10(xfull),4)
         +0.023*pow(np.log10(xfull),5))
fig.loglog(xfull,bCfull, 'black',linewidth=.5,linestyle=':')

bC3tomax = 10**(0.761+1.005*np.log10(x3tomax)+0.7*pow(np.log10(x3tomax),2)
          +0.7*pow(np.log10(x3tomax),3)+0.214*pow(np.log10(x3tomax),4)
         +0.023*pow(np.log10(x3tomax),5))

bC12 = 10**(0.761+1.005*np.log10(x12)+0.7*pow(np.log10(x12),2)
          +0.7*pow(np.log10(x12),3)+0.214*pow(np.log10(x12),4)
         +0.023*pow(np.log10(x12),5))

bC23 = 10**(0.761+1.005*np.log10(x23)+0.7*pow(np.log10(x23),2)
          +0.7*pow(np.log10(x23),3)+0.214*pow(np.log10(x23),4)
         +0.023*pow(np.log10(x23),5))

# boundary D
bD = 10**(2.141+1.686*np.log10(xminto3)-0.141*pow(np.log10(xminto3),2)
          +0.007*pow(np.log10(xminto3),3))
fig.loglog(xminto3,bD, 'gray',linewidth=.5,dashes=[4, 2, 1, 2, 1, 2])

bD23 = 10**(2.141+1.686*np.log10(x23)-0.141*pow(np.log10(x23),2)
          +0.007*pow(np.log10(x23),3))

# boundary E
bE = 10**(-1.381-2.108*np.log10(x12)-1.394*pow(np.log10(x12),2)
         -0.167*pow(np.log10(x12),3))
fig.loglog(x12,bE, 'black',linewidth=.5,linestyle='-.')



# --------------- Shading regions on the graph -------------------------

# inspiral vs mass-transfer regions (I & II)
fig.fill_between(xfull,bA,bCfull,color='lightcyan')
fig.fill_between(xfull,bBfull,bCfull,color='lightyellow')

# region III
plt.fill_between(x3tomax,bC3tomax,bB3tomax,color='bisque')
plt.fill_between(xminto3,bBminto3,bD,color='bisque')

#region IV
plt.fill_between(x12,bC12,bE, color='lavender')
plt.fill_between(x23,bC23,bD23,color='lavender')

# no DWD regions
plt.fill_between(xfull,5,bA,color='lightgray')
plt.fill_between(xfull,bBfull,1e-4,color='lightgray')




# -------------------- Labelling the graph ------------------------------

# Label the x and y axes and title the graph.

plt.xlabel(r'$f_{GW}\ \mathrm{(Hz)}$')
plt.ylabel(r'$dh_\mathrm{norm}$')
plt.title('DWD on KT plot (SNR cutoff = 1)')

# Create a custom legend, defining each marker-type explicitly.

bndA = mlines.Line2D([], [], color='black',
                           label=r'A ($M_a=M_d=M_{Ch}$)')
bndB = mlines.Line2D([], [], color='black',linestyle='--',
                           label=r'B ($q=1$)')
bndC = mlines.Line2D([], [], color='black',linestyle=':',
                           label=r'C ($M_a=M_{Ch}$)')
bndD = mlines.Line2D([], [], color='gray',dashes=[4, 2, 1, 2, 1, 2],
                          label=r'D ($q>q_{crit}$)')
bndE = mlines.Line2D([], [], color='black',linestyle='-.',
                          label=r'E ($M_{tot}=M_{Ch}$)')

reddot = mlines.Line2D([], [], color='red', marker='o',linestyle='None',
                          markersize=2, label='He-He')
orangedot = mlines.Line2D([], [], color='orange', marker='o',linestyle='None',
                          markersize=2, label='He-C/O')
yellowdot = mlines.Line2D([], [], color='yellow', marker='o',linestyle='None',
                          markersize=2, label='He-O/Ne')
greendot = mlines.Line2D([], [], color='green', marker='o',linestyle='None',
                          markersize=2, label='C/O-C/O')
bluedot = mlines.Line2D([], [], color='blue', marker='o',linestyle='None',
                          markersize=2, label='C/O-O/Ne')
purpledot = mlines.Line2D([], [], color='magenta', marker='o',linestyle='None',
                          markersize=2, label='O/Ne-O/Ne')

#noDWD_patch = mpatches.Patch(color='lightgray', label='no DWD')
inspiral_patch = mpatches.Patch(color='lightcyan', label='I: inspiral')
masstransfer_patch= mpatches.Patch(color='lightyellow', label='II: mass transfer')
region3_patch = mpatches.Patch(color='bisque', label='III: unstable MT')
region4_patch = mpatches.Patch(color='lavender',label='IV: SN Ia progenitors')

#plt.legend(fontsize=8,handlelength = 3,
#           handles=[bndA,bndB,bndC,bndD,bndE,inspiral_patch,masstransfer_patch,
 #                   region3_patch,region4_patch], loc=4)

#plt.legend(fontsize=8,handlelength = 3,
#          handles=[bndA,bndB,bndC,inspiral_patch,masstransfer_patch], loc=4)

#plt.legend(fontsize=8,handlelength = 3,
 #          handles=[reddot,orangedot,yellowdot,greendot,bluedot,purpledot,
  #                  bndA,bndB,bndC,bndD,bndE], loc=4)

plt.legend(fontsize=8,handlelength = 3,
           handles=[reddot,orangedot,yellowdot,greendot,bluedot,purpledot,
                    inspiral_patch,masstransfer_patch,region3_patch,region4_patch], loc=4)

# Additional labels
fig.text(.003,1,'no DWDs')
fig.text(.01,.001,'no DWDs')

fig.text(.004,.04,'I')
fig.text(.005,.006,'II')
fig.text(.04,.21,'III')
fig.text(.025,.15,'IV')

fig.text(.0011,.08,'A')
fig.text(.004,.0005,'B')
fig.text(.0011,.0015,'C')
fig.text(.002,.0005,'D')

fig.plot([.017,.03],[.08,.02],'black',linewidth=.3)
fig.text(.032,.015,'E')

# --------------------- Adjusting the plot's range ---------------------
#plt.xlim([5e-4,1e-1])
#plt.ylim([3e-4,1])

plt.xlim([.001,.5])
plt.ylim([3e-4,3])

# -------------- Exporting the plot as high-quality pdf  ---------------
plt.savefig("graph.pdf")
