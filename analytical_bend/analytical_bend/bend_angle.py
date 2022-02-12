from os import path
import numpy as np
import MDAnalysis
from miscell.file_util import check_dir_exist_and_make
from miscell.hd5_util import save_array_to_hd5, read_array_from_hd5
from miscell.na_bp import d_n_bp

class FoldersBuilder:
    def __init__(self, rootfolder, host):
        self.rootfolder = rootfolder
        self.host = host
        self.host_folder = path.join(rootfolder, host)

        self.discre_h5 = path.join(self.host_folder, 'discretized.hdf5')
        self.theta_h5 = path.join(self.host_folder, 'theta.hdf5')

    def initialize_folders(self):
        for folder in [self.rootfolder, self.host_folder]:
            check_dir_exist_and_make(folder)

class Discretizer(FoldersBuilder):
    key = 'discretized'

    def __init__(self, rootfolder, host):
        super().__init__(rootfolder, host)
        self.n_bp = d_n_bp[host]
        self.discre_array = None
    
    def make_discre_h5(self):
        crd, dcd = self.get_crd_dcd()
        mda_u = MDAnalysis.Universe(crd, dcd)
        n_frame = len(mda_u.trajectory)
        discre_array = np.zeros((n_frame, self.n_bp, 3))
        basepairs = self.get_basepairs()
        for frame_id in range(n_frame):
            mda_u.trajectory[frame_id]
            discre_array[frame_id] = self.get_midpoint_array(basepairs)
        self.discre_array = discre_array
        save_array_to_hd5(self.discre_h5, self.key, self.discre_array)

    def read_discre_h5(self):
        self.discre_array = read_array_from_hd5(self.discre_h5, self.key)

    def get_basepairs(self, u, short_segid=False):
        basepairs = {}
        bp_id = 1
        if short_segid:
            segid_i = 'STR1'
            segid_j = 'STR2'
        else:
            segid_i = 'STRAND1'
            segid_j = 'STRAND2'

        for resid_i in range(1, self.n_bp + 1):
            # resid_i: guide strand resid  # resid_j: target strand resid
            resid_j = miscell.get_antistrand_resid(resid_i, n_bp=self.n_bp)
            seq_i = self.sequence['guide']
            seq_j = self.sequence['target']
            resname_i = seq_i[resid_i - 1]
            resname_j = seq_j[resid_j - 1]
            temp = (u.select_atoms('segid {0} and resid {1} and name {2}'.format(segid_i, resid_i, c8_c6[resname_i])),
                    u.select_atoms('segid {0} and resid {1} and name {2}'.format(segid_j, resid_j, c8_c6[resname_j])))
            basepairs[bp_id] = temp
            bp_id += 1
        return basepairs

    def get_midpoint_array(self, basepairs):
        midpoint_array = np.zeros((self.n_bp, 3))
        for bp_id in range(self.n_bp):
            atom1 = basepairs[bp_id+1][0]
            atom2 = basepairs[bp_id+1][1]
            midpoint = (atom1.positions[0] + atom2.positions[0]) / 2
            midpoint_array[bp_id] = midpoint
        return midpoint_array

    def get_crd_dcd(self):
        return 'avg-crd', 'fit-dcd'