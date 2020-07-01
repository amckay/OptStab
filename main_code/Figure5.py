#!/usr/bin/env python
# coding: utf-8



# Compute the cyclicality of benefit payments shown in Figure 5


import numpy as np
import pandas as pd

from fredapi import Fred
fred = Fred(api_key='55d72a00536ea0269b79e23f26ba6ab0')

from matplotlib import pyplot as plt
from statsmodels.tsa.filters.hp_filter import hpfilter
import statsmodels.formula.api as smf

import os.path



# download data from FRED
CC = fred.get_series('CCSA')['1967-01-01':'2017-07-01'].asfreq(freq='QS',method= 'ffill')['1968-01-01':'2017-07-01']  #continued claims
DPI = fred.get_series('DPI')['1947-01-01':'2017-07-01']
u = fred.get_series('UNRATE')['1968-01-01':'2017-07-01'].asfreq(freq='QS') / 100
bPayments = fred.get_series('W825RC1Q027SBEA')['1968-01-01':'2017-07-01']
Pop = fred.get_series('CNP16OV').asfreq(freq='QS')['1968-01-01':'2017-07-01']*1000


#smooth DPI
DPI_cyc, DPI_trend = hpfilter(np.log(DPI),lamb=1600)
DPI_trend = np.exp(DPI_trend['1968-01-01':'2017-07-01'])


# Calculation
b = (0.5 + 0.5*bPayments/DPI_trend / (CC/Pop))**(1-0.151)


def AddRecBars(ax):
    rec = fred.get_series('USREC')['1967-10-01':'2017-07-01'].asfreq(freq='QS')
    peaks = rec.index[1:][np.where(rec.diff()[1:] == 1.0)]
    troughs = rec.index[1:][np.where(rec.diff()[1:] == -1.0)]

    y1, y2 = ax.get_ylim()
    for pt in zip(peaks,troughs):
        x = pt
        ax.fill_between(x, y1=y1, y2=y2, alpha=.25, color='k')

    ax.set_ylim(y1,y2)
    return ax


# make figure
f = plt.figure(figsize=(12,8))
plt.plot(b, color='k')
ax = f.get_axes()[0]
ax.set_ylim(0.65,0.9)
AddRecBars(ax)

FIGURE_DIR = os.path.join(os.path.pardir,'data','results')
f.savefig(os.path.join(FIGURE_DIR,'cyclical_b.png'),bbox_inches = 'tight')
plt.show()
