#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 14:24:01 2023

@author: zhangyu
"""

import os
import pandas as pd
import matplotlib.pyplot as plt

def read_entries_from_csv(file_path):
    try:
        df = pd.read_csv(file_path, comment='#')
        return df['entries']
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None

def plot_entries(folder_path):
    plt.figure(figsize=(12, 6))

    for file_name in os.listdir(folder_path):
        if file_name.endswith('.csv'):
            file_path = os.path.join(folder_path, file_name)
            entries = read_entries_from_csv(file_path)
            if entries is not None:
                plt.plot(entries, label=file_name)

    plt.xlabel(r'$\gamma$ energy[keV]')
    plt.ylabel(r'counts[$s^{-1}$]')
    plt.legend()
    plt.show()

# set file path
folder_path = '/home/zhangyu/Program_test/geant4test/HPGe-calibration/FoilCopper/15cm/cu-ion/' 
plot_entries(folder_path)

########################################

import os
import pandas as pd
import numpy as np
import re  # regulation library

import scipy.optimize as optimize

energy = np.array([121.8, 244.7, 344.3, 444, 661.7, 778.9, 964.1, 1112.1, 1173.2, 1332.5, 1408])*1e-3
eff_exp_10cm = np.array([0.0106475078435655,0.00727785333723451,0.00559899669575573,0.00444161842294093,\
                         0.00327399422074109,0.00278206217497877,0.00232325211815306,0.00204550010717191,\
                             0.00203125818304615,0.00182904446914907,0.00175869453728399])
eff_exp_15cm = np.array([5.35e-3,3.75e-3,2.88e-3,2.27e-3,1.68e-3,1.44e-3,1.23e-3,1.1e-3,1.06e-3,9.75e-4,9.22e-4])
eff_exp_20cm = np.array([3.16984e-3,2.235e-3,1.772e-3,1.438e-3,1.0264e-3,8.732e-4,7.5583e-4,6.681e-4,6.587e-4,5.9451e-4,5.5588e-4])
energy_plot = np.linspace(0.1218, 1.408)


def read_max_entry_from_csv(file_path):
    try:
        df = pd.read_csv(file_path, comment='#')
        return df['entries'].max()
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None

def extract_number_from_filename(filename):

    match = re.search(r'\d+', filename)
    return int(match.group()) if match else None

def find_max_entries(folder_path,angle):
    data = []

    for file_name in os.listdir(folder_path):
        if file_name.endswith('.csv'):
            number = extract_number_from_filename(file_name)
            if number is not None:
                file_path = os.path.join(folder_path, file_name)
                max_entry = read_max_entry_from_csv(file_path)/(1e6/(0.5*(1-np.cos(angle/180*np.pi))))
                if max_entry is not None:
                    data.append([number, max_entry])

    return np.array(data)


folder_path = '/home/zhangyu/Program_test/geant4test/HPGe-calibration/FoilCopper/15cm'
max_entries_array = find_max_entries(folder_path,14)

sorted_array = max_entries_array[max_entries_array[:, 0].argsort()]


print(max_entries_array)


def poly_func(x, *coefs):
    return sum([c * x**i for i, c in enumerate(coefs)])

#########################
compareData = eff_exp_15cm

degree = 6

coefs, _ = optimize.curve_fit(poly_func, energy, compareData, p0=np.ones(degree + 1))
##########################
degree = 6

coefs2, _ = optimize.curve_fit(poly_func, energy, sorted_array[:,1], p0=np.ones(degree + 1))
########################
difference = 1-sorted_array[:,1]/compareData


fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), gridspec_kw={'height_ratios': [5, 1]})
ax1.scatter(energy,compareData*100)
ax1.scatter(energy,sorted_array[:,1]*100)
ax1.plot(energy_plot,poly_func(energy_plot, *coefs)*100)
ax1.plot(energy_plot,poly_func(energy_plot, *coefs2)*100)
ax1.set_xlabel(r'$\gamma$ Energy[MeV]')
ax1.set_ylabel('detection efficiency[%]')
ax1.legend(['experiment data','simulation data','experiment fitting curve','simulation fitting curve'])
for i in range(len(energy)):
    ax1.text(energy[i], compareData[i]*100+0.01, f"{compareData[i]*100:.4f}%", ha='center', va='bottom',color='blue')
for i in range(len(energy)):
    ax1.text(energy[i], sorted_array[i,1]*100-0.025, f"{sorted_array[i,1]*100:.4f}%", ha='center', va='bottom',color='orange')

ax2.plot(energy, difference*100, label='y1 - y2', color='green')
ax2.set_ylabel('Difference')
for i in range(len(energy)):
    ax2.text(energy[i], difference[i]*100, f"{difference[i]*100:.2f}%", ha='center', va='bottom',color='green')
