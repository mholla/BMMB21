#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This code reproduces the analysis of gyral-sulcal thickness ratio in 
human brains from Wang et al. 2021 (DOI: 10.1007/s10237-020-01400-w)

Input files (a text file listing subject names, and the subject 
folder) should be placed into this directory.  The ratio array is 
written to a .asc file in the user-defined folder for each subject.

"""

def t_ratio():

    import nibabel as nib
    import numpy as np
    import os
    
    
    input_txt = 'Yale_TD.txt'
    path = os.path.join(os.getcwd(),input_txt)
    with open(path) as f: subjects = f.read().splitlines()
    
    for subject in subjects:

        print(subject)
        
        subjects_name = 'Yale_TD'
        subjects_dir = os.path.join(os.getcwd(),subjects_name)
        hemi = 'lh'
        surf = 'pial'
        curv = 'H'
        
        c = nib.freesurfer.io.read_morph_data(os.path.join(subjects_dir , subject, 'surf', '{h}.{s}.{c}.crv'.format(h=hemi, s=surf, c=curv)))
        t = nib.freesurfer.io.read_morph_data(os.path.join(subjects_dir , subject, 'surf', '{h}.thickness'.format(h=hemi)))
        
        gyral_sum = np.zeros(len(c))
        sulcal_sum = np.zeros(len(c))
                    
        if curv == 'H':
            
            c = -c
            for i in range(len(c)):
                if c[i] > 0.0 and c[i] <= 0.5:
                    gyral_sum[i] = c[i]
                
            for i in range(len(c)):
                if c[i] < 0.0 and c[i] >= -0.5:
                    sulcal_sum[i] = c[i]
            
            t_gyral = np.zeros(len(t))
            t_sulcal = np.zeros(len(t))
            
            for i in range(len(t)):
                if gyral_sum[i] != 0:
                    t_gyral[i] = t[i]
                    
            for i in range(len(t)):
                if sulcal_sum[i] != 0:
                    t_sulcal[i] = t[i]
        else:
            exit
        
        n_g = np.count_nonzero(t_gyral)
        n_s = np.count_nonzero(t_sulcal)
        
        coverage = np.array([0.01/100, 0.02/100, 0.04/100, 0.06/100, 0.08/100, 0.1/100, 0.2/100, 0.4/100, 0.6/100, 0.8/100, 1/100, 2/100, 4/100,
                             6/100, 8/100, 10/100, 20/100, 30/100, 40/100, 50/100, 60/100, 70/100, 80/100, 90/100, 1])
    
        N_g = n_g * coverage
        N_s = n_s * coverage

        # Cortical thickness ratio calculation for gyral and sulcal regions depending on sampling ratio (coverage)
        ratio = np.zeros(len(N_g))
        gyral_sorted = np.sort(gyral_sum)
        gyral_sorted_rev = gyral_sorted[::-1]
        sulcal_sorted = np.sort(sulcal_sum)
        
        for j in range(len(N_g)):
                
            t_g = np.zeros(len(t))
            t_s = np.zeros(len(t))
            
            gyral_reduced = gyral_sorted_rev[:int(N_g[j])]
            sulcal_reduced = sulcal_sorted[:int(N_s[j])]
            
            for k in range(len(gyral_sum)):
                if gyral_sum[k] >= min(gyral_reduced):
                    t_g[k] = t_gyral[k]
            
            for m in range(len(sulcal_sum)):
                if sulcal_sum[m] <= max(sulcal_reduced):
                    t_s[m] = t_sulcal[m]
            
            t_g = t_g[t_g != 0]
            t_s = t_s[t_s != 0]
                   
            ratio[j] = np.mean(t_g)/np.mean(t_s)   

        array_name = os.path.join(os.getcwd(), 'thickness_ratios_{c}_{s}.asc'.format(c=curv, s=subject))
        np.savetxt(array_name, ratio, fmt='%10.6f', delimiter=' '' ')
        
        return 
    
t_ratio()