#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This code reproduces the analysis of mean gyral and sulcal 
thicknesses in human brains from Wang et al. 2021 
(DOI: 10.1007/s10237-020-01400-w)

Input files (a text file listing subject names, and the subject 
folder) should be placed into this directory.  The ratio array is 
written to a .asc file in the user-defined folder for each subject.

"""

def t_mean():

    import nibabel as nib
    import numpy as np
    import os
    
    input_txt = 'Yale_TD.txt'
    path = os.path.join(os.getcwd(),input_txt)
    with open(path) as f: subjects = f.read().splitlines()
    
    j = 0
    
    mean_t = np.zeros((len(subjects)+1, 2))
    
    for subject in subjects:
        
        print(subject)
        
        subjects_name = 'Yale_TD'
        subjects_dir = os.path.join(os.getcwd(),subjects_name)
        hemi = 'lh'
        surf = 'pial'
        curv = 'H'
        
        c = nib.freesurfer.io.read_morph_data(os.path.join(subjects_dir, subject, 'surf', '{h}.{s}.{c}.crv'.format(h=hemi, s=surf, c=curv)))
        t = nib.freesurfer.io.read_morph_data(os.path.join(subjects_dir, subject, 'surf', '{h}.thickness'.format(h=hemi)))
        
        gyral_sum = np.zeros(len(c))
        sulcal_sum = np.zeros(len(c))
        
        c = -c
        for i in range(len(c)):
            if c[i] > 0.0 and c[i] <= max(c):
                gyral_sum[i] = c[i]
            
        for i in range(len(c)):
            if c[i] < 0.0 and c[i] >= min(c):
                sulcal_sum[i] = c[i]
        
        t_gyral = np.zeros(len(t))
        t_sulcal = np.zeros(len(t))
        
        for i in range(len(t)):
            if gyral_sum[i] != 0:
                t_gyral[i] = t[i]
                
        for i in range(len(t)):
            if sulcal_sum[i] != 0:
                t_sulcal[i] = t[i]
                
        t_gyral = t_gyral[t_gyral != 0]
        t_sulcal = t_sulcal[t_sulcal != 0]
        
        t_gyral_mean = np.mean(t_gyral)
        mean_t[j,0] = t_gyral_mean
        
        t_sulcal_mean = np.mean(t_sulcal)
        mean_t[j,1] = t_sulcal_mean
        
        j = j + 1
        
    mean_t[len(subjects), 0] = np.mean(mean_t[0:len(subjects)-1,0])
    mean_t[len(subjects), 1] = np.mean(mean_t[0:len(subjects)-1,1])
    
    array_name = os.path.join(os.getcwd(), 'mean_thickness.asc')
    np.savetxt(array_name, mean_t, fmt='%10.6f', delimiter=' '' ')
    
    return 

t_mean()
