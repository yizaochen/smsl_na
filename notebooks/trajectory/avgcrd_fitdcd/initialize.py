import imp
from trajectory.crd_dcd import FoldersBuilder

rootfolder = '/home/ytcdata/traj/'
suffix = 'before_clean'
host_lst = ['pnas_dna']

for host in host_lst:
    builder = FoldersBuilder(rootfolder, suffix, host)
    builder.initialize_folders()
    builder.copy_input_gro_xtc()