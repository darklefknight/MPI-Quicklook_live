#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Personal Tooblox from Tobias Machnitzki
(tobias.machnitzki@mpimet.mpg.de)
"""

import numpy as np




#%%

def Dewpoint(T,RH):
    '''
    Calculates the Dewpointtemperature from a given Temperature T and Relative Humidity RH
    #z_e = water vapor pressure
    #z_es = saturation water vapor pressure
    #T = temperature
    #tempd = dewpoint temperature
    #RH = relative Humidity '''

    z_es = 610.78*np.exp(17.08085*(T)/(235+(T)))
    z_e = RH*z_es/100
    tempd = (235*np.log(z_e/610.5))/(17.1-np.log(z_e/610.78))
    return tempd


#%%
def nanargmin(a):
    '''
    Calculates the Position of the minimum of a n-dimensional Array with nan-values.
    This function was created with code-segments from Marcus Klingebiel.
    '''
    idx = np.nanargmin(a, axis=None)
    multi_idx = np.unravel_index(idx, a.shape)
    if np.isnan(a[multi_idx]):
        nan_count = np.sum(np.isnan(a))
        # In numpy < 1.8 use idx = np.argsort(a, axis=None)[-nan_count-1]
        idx = np.argpartition(a, -nan_count-1, axis=None)[-nan_count-1]
        multi_idx = np.unravel_index(idx, a.shape)
    return multi_idx


#%%
def nanargmax(a):
    '''
    Calculates the Position of the maximum of a n-dimensional Array with nan-values.
    This function was created with code-segments from Marcus Klingebiel.
    '''
    idx = np.nanargmax(a, axis=None)
    multi_idx = np.unravel_index(idx, a.shape)
    if np.isnan(a[multi_idx]):
        nan_count = np.sum(np.isnan(a))
        # In numpy < 1.8 use idx = np.argsort(a, axis=None)[-nan_count-1]
        idx = np.argpartition(a, -nan_count-1, axis=None)[-nan_count-1]
        multi_idx = np.unravel_index(idx, a.shape)
    return multi_idx
