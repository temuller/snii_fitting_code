#!/usr/bin/env python
'''
SNIIP Fitting Code Output Summary

This script plots the output of Ondrej's fitting code, modified by Tomás Müller
'''

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors

sn = input('Enter a SN named(e.g. SN2013am): ')
phot_file = f'{sn}/{sn}_phot.dat'
model_phot_file = f'{sn}/{sn}_model_phot.dat'
velo_file = f'{sn}/{sn}_velo.dat'
model_velo_file = f'{sn}/{sn}_model_velo.dat'

filter_labels = {'U':0,'B':1,'V':2,'R':3,'I':4,'J':5,'H':6,'K':7,'g':9,'r':10,'i':11,
                    'z':12,'uvw2':13,'uvm2':14,'uvw1':15,'u':16,'b':17,'v':18,'Z':19,'Y':20}

filter_labels_inv = {value:key for key, value in filter_labels.items()}  # for plotting purposes

##############################
# photometry
columns = ['MJD', 'mag', 'err', 'filters', 'star_no', 'dy', 'Pi', 'tau', '1-w']
phot_df = pd.read_csv(phot_file, delim_whitespace=True, names=columns)


columns = ['MJD']+list(filter_labels.keys())
model_phot_df = pd.read_csv(model_phot_file, delim_whitespace=True, names=columns, skiprows=1)

# velocity
columns = ['MJD', 'vel', 'err', 'filters', 'star_no', 'dy', 'Pi', '1-w']
velo_df = pd.read_csv(velo_file, delim_whitespace=True, names=columns)

model_velo_df = pd.read_csv(model_velo_file, delim_whitespace=True, names=['MJD', 'vel'])
##############################

colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
color = list(colors.values())[10::5]
uv = ['U', 'uvw2','uvm2','uvw1','u','b','v']
opt = ['B','V','R','I','g','r','i','z']
ir = ['J','H','K','Z','Y']

fig, ax = plt.subplots(2, 2, sharex=True, figsize=(12,12))

# UV
for i, label in enumerate(phot_df.filters.unique()):
    band = filter_labels_inv[label]
    if band in uv:
        filter_df = phot_df[phot_df.filters==label]
        filter_df = phot_df[phot_df.filters==label]
        ax[0,0].errorbar (filter_df.MJD, filter_df.mag, filter_df.err,
                              fmt='.', label=band, color=color[i], mec='k', ms=10)
        ax[0,0].plot(model_phot_df.MJD, model_phot_df[band], color=color[i])
ax[0,0].set_ylabel('Apparent Magnitude', fontsize=16)
ax[0,0].legend()

# Opt
for i, label in enumerate(phot_df.filters.unique()):
    band = filter_labels_inv[label]
    if band in opt:
        filter_df = phot_df[phot_df.filters==label]
        filter_df = phot_df[phot_df.filters==label]
        ax[0,1].errorbar (filter_df.MJD, filter_df.mag, filter_df.err,
                              fmt='.', label=band, color=color[i], mec='k', ms=10)
        ax[0,1].plot(model_phot_df.MJD, model_phot_df[band], color=color[i])
ax[0,1].set_ylabel('Apparent Magnitude', fontsize=16)
ax[0,1].legend()

# IR
for i, label in enumerate(phot_df.filters.unique()):
    band = filter_labels_inv[label]
    if band in ir:
        filter_df = phot_df[phot_df.filters==label]
        filter_df = phot_df[phot_df.filters==label]
        ax[1,0].errorbar (filter_df.MJD, filter_df.mag, filter_df.err,
                              fmt='.', label=band, color=color[i], mec='k', ms=10)
        ax[1,0].plot(model_phot_df.MJD, model_phot_df[band], color=color[i])
ax[1,0].set_xlabel('MJD', fontsize=16)
ax[1,0].set_ylabel('Apparent Magnitude', fontsize=16)
ax[1,0].legend()

# Vel
ax[1,1].errorbar(velo_df.MJD, velo_df.vel, velo_df.err,
                      fmt='.', color='k', ms=10)
ax[1,1].plot(model_velo_df.MJD, model_velo_df.vel, color='k')
ax[1,1].set_xlabel('MJD', fontsize=16)
ax[1,1].set_ylabel('Velocity [km/s]', fontsize=16)

fig.suptitle(sn, fontsize=20)
y_min, y_max = phot_df.mag.min()-0.5, phot_df.mag.max()+0.5
ax[0,0].set_ylim(y_min, y_max)
ax[0,1].set_ylim(y_min, y_max)
ax[1,0].set_ylim(y_min, y_max)
ax[1,1].set_ylim(velo_df.vel.min()-500, velo_df.vel.max()+500)
x_min, x_max = phot_df.MJD.min()-15, phot_df.MJD.max()+20
ax[0,0].set_xlim(x_min, x_max)
ax[0,1].set_xlim(x_min, x_max)
ax[1,0].set_xlim(x_min, x_max)
ax[1,1].set_xlim(x_min, x_max)
ax[0,0].set_facecolor("lightgrey")
ax[0,1].set_facecolor("lightgrey")
ax[1,0].set_facecolor("lightgrey")
ax[1,1].set_facecolor("lightgrey")
ax[0,0].invert_yaxis()
ax[0,1].invert_yaxis()
ax[1,0].invert_yaxis()
plt.show()
