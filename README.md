# BMMB21

This code reproduces the results of Wang et al. 2021 (https://doi.org/10.1007/s10237-020-01400-w)


## in silico computational simulations

Input files and Abaqus Explicit user subroutines are provided for the reproduction of all the simulations included in the paper. 


## in vivo human brain data

The provided codes calculate the mean thickness in the gyral and sulcal regions of the subject brains, and the gyral-sulcal thickness ratio at different levels of coverage.  They read in the mean curvature (?h.pial.H.crv) and cortical thickness (?h.thickness) output files, created by the automatic neuroimaging pipeline Freesurfer (--recon-all), for a given set of subjects.  In Wang et al. 2021, we used the control subjects of the ABIDE-YALE dataset; the input files Yale_TD.txt and Yale_TD are included here.


## figure plotting

Python scripts are provided for the plotting of Figs. 2 and 9 in the paper.  Fig. 9 uses data compiled from the results of the human_data scripts.