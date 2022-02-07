from os import path
from itertools import combinations
import MDAnalysis as mda
import numpy as np
import pandas as pd
import h5py
from curves.curves_main_util import PreliminaryAgent
from miscell.file_util import check_dir_exist_and_make, check_file_exist, copy_verbose
from miscell.na_bp import d_n_bp
from pdb_util.atom import Atom
from pdb_util.pdb import PDBWriter

class FoldersBuilder(PreliminaryAgent):
    def __init__(self, rootfolder, host):
        super().__init__(rootfolder, host)
        self.haxis_discretize_folder = path.join(self.host_folder, 'haxis_discretize')
        self.discre_pdb_dcd_folder = path.join(self.host_folder, 'discre_pdb_dcd')
        self.bend_shape_folder = path.join(self.host_folder, 'bend_shape')

        self.dcre_pdb = path.join(self.discre_pdb_dcd_folder, 'discretized.pdb')
        self.dcre_dcd = path.join(self.discre_pdb_dcd_folder, 'discretized.dcd')
        self.make_haxis_tcl = '/home/yizaochen/codes/smsl_na/curves/tcl_scripts/make_haxis.tcl'

        self.lnorm_theta_h5 = path.join(self.bend_shape_folder, 'l_norm_theta.hdf5')

    def initialize_folders(self):
        for folder in [self.haxis_discretize_folder, self.discre_pdb_dcd_folder, self.bend_shape_folder]:
            check_dir_exist_and_make(folder)

class Discretizer(FoldersBuilder):

    def discretize_haxis_to_pdb(self, start_frame, stop_frame):
        for frame_id in range(start_frame, stop_frame):
            nodes_atg = self.extract_nodes_from_smooth_haxis(frame_id)
            self.write_discretized_pdb(frame_id, nodes_atg)
            self.print_progress(frame_id, 500)

    def extract_nodes_from_smooth_haxis(self, frame_id):
        input_pdb = path.join(self.haxis_smooth_folder, f'{frame_id}.pdb')
        return NodesExtract(input_pdb, self.n_bp).get_atomgroups()

    def write_discretized_pdb(self, frame_id, atom_groups):
        output_pdb = path.join(self.haxis_discretize_folder, f'{frame_id}.pdb')
        PDBWriter(output_pdb, atom_groups).write_pdb()

    def print_progress(self, frame_id, interval):
        if frame_id % interval == 0:
                print(f'Discretize smooth helical axis for {self.host}, Frame-ID: {frame_id}')

    
class NodesExtract:
    def __init__(self, pdb, n_bp):
        self.pdb = pdb
        self.n_bp = n_bp
        self.u = mda.Universe(pdb)

    def get_atomgroups(self):        
        atomgroups = list()
        for resid in range(1, self.n_bp):
            atomgroups.append(self.get_atom(resid, False))
        atomgroups.append(self.get_atom(resid, True))
        return atomgroups

    def get_atom(self, resid, tail):
        if tail:
            atom_data = ['ATOM', resid+1, 'S', 'HAI', resid+1]
        else:
            atom_data = ['ATOM', resid, 'S', 'HAI', resid]
        xyz_list = self.get_xyz_list(resid, tail)
        atom_data += xyz_list
        atom_data += [1.00, 1.00]
        return Atom(atom_data, False)

    def get_xyz_list(self, resid, tail):
        atg_sele = self.u.select_atoms(f'resid {resid}')
        if not tail:
            return list(atg_sele.positions[0])
        else:
            return list(atg_sele.positions[-1])

class Compressor(FoldersBuilder):
    
    def make_ref_pdb(self):
        if not check_file_exist(self.dcre_pdb):
            first_pdb = path.join(self.haxis_discretize_folder, '0.pdb')
            copy_verbose(first_pdb, self.dcre_pdb)

    def compress_allpdbs_to_dcd(self, start_frame, stop_frame):
        print('----Enter the following into terminal----')
        print(f'cd {self.haxis_discretize_folder}')
        print('vmd')
        print('----Enter the following into VMD tkconsole----')
        print(f'source {self.make_haxis_tcl}')
        print(f'read_all_pdb_files {start_frame} {stop_frame-1}')
        print(f'animate write dcd {self.dcre_dcd} beg {start_frame} end {stop_frame-1} waitfor all')
        print('----Finish----')

    def check_result_fit_original_traj(self):
        print('----Enter the following into terminal----')
        print(f'vmd -pdb {self.dcre_pdb} {self.dcre_dcd}')
        print('----Enter the following into VMD tkconsole----')
        print(f'mol new {self.input_pdb}')
        print(f'mol addfile {self.input_xtc} type xtc first 0 last -1 step 1 waitfor all 1' )
        print('----Finish----')


class LNormTheta(FoldersBuilder):
    columns = ['Frame_ID', 'i', 'j', '|l_i|', '|l_j|', 'theta']

    def __init__(self, rootfolder, host):
        super().__init__(rootfolder, host)
        self.u = mda.Universe(self.dcre_pdb, self.dcre_dcd)
        self.n_bead = d_n_bp[host]
        self.n_bead_minus_1 = self.n_bead - 1
        self.d_result = None

    def make_l_norm_theta_hd5(self):
        self.set_d_result_by_iterate_traj()
        self.save_d_result_to_hd5()

    def read_l_modulus_theta_hd5(self):
        h5_handle = h5py.File(self.lnorm_theta_h5, 'r')
        self.d_result = {key: h5_handle.get('key') for key in self.columns}
        h5_handle.close()

    def get_ensemble_average_l(self):
        l_result = list()
        for frame_id in range(len(self.u.trajectory)):
            self.u.trajectory[frame_id]
            vectors = [self.u.atoms.positions[i + 1] - self.u.atoms.positions[i] for i in range(self.n_bead_minus_1)]
            l_result += [np.linalg.norm(vector) for vector in vectors]
        l_result = np.array(l_result)
        return l_result.mean()

    def set_d_result_by_iterate_traj(self):
        d_result = {key: list() for key in self.columns}
        pair_list = list(combinations(range(self.n_bead_minus_1), 2))
        for ts in self.u.trajectory:
            vectors = [self.u.atoms.positions[i + 1] - self.u.atoms.positions[i] for i in range(self.n_bead_minus_1)]
            for i, j in pair_list:
                li_norm, lj_norm, theta = self.get_modulus_angle_between_two_vectors(vectors[i], vectors[j])
                d_result['Frame_ID'].append(ts.frame)
                d_result['i'].append(i)
                d_result['j'].append(j)
                d_result['|l_i|'].append(li_norm)
                d_result['|l_j|'].append(lj_norm)
                d_result['theta'].append(theta)
        self.d_result = d_result

    def save_d_result_to_hd5(self):
        h5_handle = h5py.File(self.lnorm_theta_h5, 'w') 
        for key in self.columns:
            h5_handle.create_dataset(key, data=np.array(self.d_result[key]))
        h5_handle.close()
        print(f'make {self.lnorm_theta_h5}')

    def get_modulus_angle_between_two_vectors(v1, v2):
        v1_modulus = np.linalg.norm(v1)
        v2_modulus = np.linalg.norm(v2)
        v1_u = v1 / v1_modulus
        v2_u = v2 / v2_modulus
        angle = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
        return v1_modulus, v2_modulus, angle