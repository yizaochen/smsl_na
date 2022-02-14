"""
Reference: ./deprecated/1_make_noh_avg_fit_dcd.ipynb
"""
from trajectory.crd_dcd import FoldersBuilder


rootfolder = '/home/ytcdata/traj/'
suffix = 'before_clean'
host_lst = ['pnas_dna']

for host in host_lst:
    builder = FoldersBuilder(rootfolder, suffix, host)
    builder.initialize_folders()
    builder.copy_input_gro_xtc() # This should be replaced when the locations of input gro and xtc change