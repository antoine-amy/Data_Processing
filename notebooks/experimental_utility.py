import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import quad
import scipy.special as sc
from scipy.integrate import odeint
import scipy.integrate as integ
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import time
import json
import os

import math
from IPython.display import display, clear_output
import scipy.optimize as spo
import iminuit
import uproot


def read_file(filename, samples=1024):
    data=np.fromfile(filename, dtype=np.int16)
    n=len(data)/samples
    data=np.array(np.array_split(np.array(data),n))
    return data

def area(dataLED):
    n=len(dataLED)
    print("Number of events: ", n)

    nn=len(dataLED[0]) #nn=1024 is the number of samples, each taken in an interval of 10ns
    print('Number of samples (10ns intervals): ', nn)

    inf=int(nn/2)-512  #lower limit of the range in which we search for the maximum
    sup=int(nn/2)+512  #Upper limit of the range in which we search for the maximum
    print("Maximum research interval: [", inf, ",", sup,"]")
    BASE, MAX_POS , MAX, INF, SUP, A, B, AREA, TAU, ENT = [], [], [], [], [], [], [], [], [], []

    for i in range(n):
        wf=(dataLED[i])*(-1)

        max_pos=inf+np.where(wf[inf:sup]==np.max(wf[inf:sup]))[0][0]  #Search of the maximum in the defined interval [inf, sup]

        p=max_pos-20
        bl=np.mean(dataLED[i][p-40:p]) #definition of the baseline as an average over 40 samples just before the selected peak
        std=np.std(dataLED[i][p-40:p]) #Definition of the standard deviation on 40 samples just before the selected peak

        wf=(dataLED[i]-bl)*(-1)

        maxx=wf[max_pos]   #value of the waveform maximum

        #definition of integration limits
        try:
            if wf[max_pos]>3*std :
                idx1= np.where(wf[0:max_pos]<3*std)[0][-1]

            else:
                idx1= np.where(wf[0:max_pos]<0)[0][-1]

        except:
            print("Unable to determine the limits of an integral for waveform: ", i)

        if wf[max_pos]>3*std :
            try:
                idx2= max_pos+np.where(wf[max_pos:nn]<3*std)[0][0]

            except:
                print(f"Unable to determine the integral limit b for waveform: ",i)
        else:
            try:
                idx2= max_pos+np.where(wf[max_pos:nn]<0)[0][0]

            except:
                print(f"Unable to determine the integral limit b for waveform: ",i)

        a=idx1-2
        b=idx2+2

        area=np.sum(wf[a:b]) #integral calculation

        #entropy definition
        if np.sum(wf) > 1:
                norm = np.abs(wf[wf!=0]/np.sum(wf))
                entropy = -np.sum(norm*np.log10(norm))
        else: entropy = 0

        #tau definition
        try:
                tt10 = np.where(wf[max_pos:]<maxx*0.1)[0][0] + max_pos
                tt90 = np.where(wf[max_pos:]<maxx*0.9)[0][0] + max_pos
                tau = tt10 - tt90
        except:
                tau = 0

        BASE.append(bl)
        MAX_POS.append(max_pos)    
        MAX.append(maxx)
        INF.append(inf)
        SUP.append(sup)
        A.append(a)
        B.append(b)
        AREA.append(area)
        TAU.append(tau)
        ENT.append(entropy)

    #dataframe creation    
    data = pd.DataFrame(columns=['pos_max','wf_max', 'inf', 'sup', 'a','b','area', 'tau', 'entropy'])
    data['inf']=INF; data['sup']=SUP
    data['wf_max']=MAX; data['pos_max']=MAX_POS
    data['a']=A; data['b']=B
    data['area']=AREA; data['tau']=TAU; data['entropy']=ENT

    data.to_hdf(f'data_ABALONE_.h5', key='df', mode='w')
    return data