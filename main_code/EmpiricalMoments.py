# -*- coding: utf-8 -*-
"""
Compute empirical moments for calibration

Created on Thu Apr 16 13:40:28 2015

@author: iragm01
"""

import pandas as pd
import statsmodels.api as sm
import numpy as np
from scipy.optimize import brentq
import json
import os

START = '1960'
END = '2014'
WRITE_MOMENTS = False

print 'Moments using data from {s} through {e}'.format(s = START, e = str(int(END)-1))

Moments = {}

Y = pd.Series.from_csv('/Users/iragm01/projects/OptStab/data/FRED/GDP.csv')
G = pd.Series.from_csv('/Users/iragm01/projects/OptStab/data/FRED/GCE.csv')
u = pd.Series.from_csv('/Users/iragm01/projects/OptStab/data/FRED/unrateq.csv')
P = pd.Series.from_csv('/Users/iragm01/projects/OptStab/data/FRED/CorePCE.csv')

inflation = np.log(P).diff()[START:END]

GY = G[START:END]/Y[START:END]

gdp_cycle, gdp_trend = sm.tsa.filters.hpfilter(np.log(Y[START:END]),lamb = 1600)

print 'StDev log Y  = {0}'.format(gdp_cycle.std())
Moments['StDevY'] = gdp_cycle.std()

print 'StDev inflation  = {0}'.format(inflation.std())
Moments['StDevinfl'] = inflation.std()

GToPotentialY = G[START:END]/np.exp(gdp_trend)


print 'Mean of G/Y  = {0}'.format(GY.mean())
Moments['meanGY'] = GY.mean()
print 'StDev of G/Y  = {0}'.format(GY.std())
Moments['StDevGY'] = GY.std()
print 'Mean of G/Y*  = {0}'.format(GToPotentialY.mean())
Moments['meanGToPotentialY'] = GToPotentialY.mean()
print 'StDev of G/Y*  = {0}'.format(GToPotentialY.std())
Moments['StDevGtoPotentialY'] = GToPotentialY.std()
print 'Covariance of G/Y* with output gap = {0}'.format(GToPotentialY.cov(gdp_cycle))
Moments['CovGtoPotentialYwithOutputGap'] = GToPotentialY.cov(gdp_cycle)


G_cycle, G_trend = sm.tsa.filters.hpfilter(np.log(G[START:END]),lamb = 1600)

print 'Autocorrelation of detrended log G  = {0}'.format(G_cycle.autocorr())
Moments['AutocorrDetrendedLogG'] = G_cycle.autocorr()

u = u/100
print 'Mean of u = {0}'.format(u[START:END].mean())
#u = sm.tsa.filters.hpfilter(u[START:END],lamb = 100000)[0]
Moments['meanu'] = u[START:END].mean()
print 'StDev of u = {0}'.format(u[START:END].std())
Moments['StDevu'] = u[START:END].std()
print 'Covariance of u and detrended log Y = {0}'.format(u.cov(gdp_cycle))
Moments['CovuOutputGap'] = u.cov(gdp_cycle)



##### Shimer calculations for the job-separation rate
urate = pd.Series.from_csv('/Users/iragm01/projects/OptStab/data/FRED/unrate.csv')
U = pd.Series.from_csv('/Users/iragm01/projects/OptStab/data/FRED/UNEMPLOY.csv')
US = pd.Series.from_csv('/Users/iragm01/projects/OptStab/data/FRED/UEMPLT5.csv')

f = -np.log((U.shift(-1)-US.shift(-1))/U)
f.index = pd.date_range(U.index[0], periods=len(f),freq='M')

def geturatestar(s, f):
    return s/(s+f)

def checkEq4(s_t,f_t,urate_t,urate_tplus1):
    # Returns the residual in equation 4
    uratestar_t = geturatestar(s_t,f_t)
    lamb = 1 - np.exp(-s_t - f_t)
    return -urate_tplus1 + lamb * uratestar_t + (1-lamb)*urate_t


urate = urate.values/100
T = len(urate)
svalues = np.zeros(T)
svalues[-1] = np.NaN
for t in range(0,T-1):
    svalues[t] = brentq(checkEq4,0.001,0.3,args = (f[t],urate[t],urate[t+1]))

index = pd.date_range(U.index[0], periods=T,freq='M')
MyD = pd.DataFrame({'f':f,'s':svalues},index = index)

#MyD has monthly Poisson outflow and inflow rates.  Now we convert to the quarterly rates
MyD = 1.0 - np.exp(-12./4 * MyD)
MyD = MyD.resample('QS')

# print 'Average unemployment inflow rate = {0}'.format(MyD['s'][START:END].mean())
# Moments['meanuInflowRate'] = MyD['s'][START:END].mean()


# hours adjustment
EMRatio = pd.Series.from_csv('/Users/iragm01/projects/OptStab/data/FRED/EMRATIOq.csv')
HoursPerWorkerIndex = pd.Series.from_csv('/Users/iragm01/projects/OptStab/data/FRED/HoursPerWorkerIndex.csv')

TotalHoursPerCapita = EMRatio * HoursPerWorkerIndex / 1000

EMRatio = sm.tsa.filters.hpfilter(np.log(EMRatio[START:END]),lamb = 100000)[0]
HoursPerWorkerIndex = sm.tsa.filters.hpfilter(np.log(HoursPerWorkerIndex[START:END]),lamb = 100000)[0]
TotalHoursPerCapita = sm.tsa.filters.hpfilter(np.log(TotalHoursPerCapita[START:END]),lamb = 100000)[0]

Moments['StDevEMratio'] = EMRatio.std()
print 'St Dev Emp-Pop Ratio = {0}'.format(Moments['StDevEMratio'])

Moments['StDevHoursPerWorkerIndex'] = HoursPerWorkerIndex.std()
print 'St Dev Hours Per Worker Ratio = {0}'.format(Moments['StDevHoursPerWorkerIndex'])

Moments['StDevTotalHours'] = TotalHoursPerCapita.std()
print 'St Dev TotalHours = {0}'.format(Moments['StDevTotalHours'])

print 'Correlation of u and hours per worker = {0}'.format(u.corr(HoursPerWorkerIndex))
Moments['Corr_U_HoursPerWorkerIndex'] = u.corr(HoursPerWorkerIndex)

# import matplotlib.pyplot as plt
# plt.subplot(2,2,1)
# plt.plot(EMRatio+HoursPerWorkerIndex)
# plt.title('Total Hours')
# plt.subplot(2,2,2)
# plt.plot(HoursPerWorkerIndex)
# plt.title('Per Worker')
# plt.subplot(2,2,3)
# plt.plot(EMRatio)
# plt.title('E-P Ratio')
# plt.show()

if WRITE_MOMENTS:
    with open(os.path.join(os.pardir,'data','calibration','Moments.txt'), 'w') as outfile:
        json.dump(Moments, outfile)
else:
    print("Moments not written to file")


# inflow rate = (1-u)*upsilon
# unemployment rate  = upsilon*(1-qM)

# q = ??
# u = 0.0610422322775
# upsilon = 0.0919001812373 / (1-u)
# M = (1 - u / upsilon) / q
