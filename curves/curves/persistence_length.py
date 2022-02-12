from os import path
from itertools import combinations
import MDAnalysis as mda
import numpy as np
import pandas as pd
from miscell.hd5_util import save_d_result_to_hd5, read_d_result_from_hd5
from miscell.pd_util import write_csv
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
        self.fourier_amp_h5 = path.join(self.bend_shape_folder, 'fourier_amplitude.hdf5')

        self.lp_big_wn_csv = path.join(self.bend_shape_folder, 'persistence_length_big_wn.csv')
        self.lp_small_wn_csv = path.join(self.bend_shape_folder, 'persistence_length_small_wn.csv')

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
        save_d_result_to_hd5(self.lnorm_theta_h5, self.columns, self.d_result)

    def read_l_modulus_theta_hd5(self):
        self.d_result = read_d_result_from_hd5(self.lnorm_theta_h5, self.columns)

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
                li_norm, lj_norm, theta = LNormTheta.get_modulus_angle_between_two_vectors(vectors[i], vectors[j])
                d_result['Frame_ID'].append(ts.frame)
                d_result['i'].append(i)
                d_result['j'].append(j)
                d_result['|l_i|'].append(li_norm)
                d_result['|l_j|'].append(lj_norm)
                d_result['theta'].append(theta)
        self.d_result = d_result

    @staticmethod
    def get_modulus_angle_between_two_vectors(v1, v2):
        v1_modulus = np.linalg.norm(v1)
        v2_modulus = np.linalg.norm(v2)
        v1_u = v1 / v1_modulus
        v2_u = v2 / v2_modulus
        angle = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
        return v1_modulus, v2_modulus, angle

class FourierShape(FoldersBuilder):
    n_begin = 0 # the first Fourier-mode 
    n_end = 9   # the last Fourier-mode

    def __init__(self, rootfolder, host):
        super().__init__(rootfolder, host)
        self.mode_lst_int = list(range(self.n_begin, self.n_end+1))
        self.mode_lst_str = [f'{mode_id}' for mode_id in self.mode_lst_int]
        self.bp_id_first, self.bp_id_last = self.get_bp_id_first_last()

        self.df_lnorm_theta = None
        self.d_result = None

    def make_fourier_amp_h5(self):
        last_frame = self.get_last_frame()
        d_result = {mode_id: list() for mode_id in self.mode_lst_str}
        self.set_df_lnorm_theta()
        for frame_id in range(1, last_frame):
            df_filter = self.get_filter_df(frame_id)
            for mode_id in self.mode_lst_int:
                d_result[f'{mode_id}'].append(self.get_an(mode_id, df_filter))
        self.d_result = d_result
        save_d_result_to_hd5(self.fourier_amp_h5, self.mode_lst_str, self.d_result)

    def read_fourier_amp_h5(self):
        self.d_result = read_d_result_from_hd5(self.fourier_amp_h5, self.mode_lst_str)

    def get_last_frame(self):
        u = mda.Universe(self.dcre_pdb, self.dcre_dcd)
        return len(u.trajectory)

    def set_df_lnorm_theta(self):
        if self.df_lnorm_theta is None:
            d_lnorm_theta = read_d_result_from_hd5(self.lnorm_theta_h5, LNormTheta.columns)
            self.df_lnorm_theta = pd.DataFrame(d_lnorm_theta)    

    def get_filter_df(self, frame_id):
        mask = (self.df_lnorm_theta['i'] == self.bp_id_first)
        df0 = self.df_lnorm_theta[mask]
        mask = (df0['Frame_ID'] == frame_id)
        df1 = df0[mask]
        mask = (df1['j'].between(self.bp_id_first+1, self.bp_id_last))
        df2 = df1[mask]
        return df2

    def get_an(self, n, df):
        length_angstrom = df['|l_j|'].sum() # L, unit: angstrom
        length_nm = df['|l_j|'].sum() / 10 # L, unit: nm
        scale_factor = np.sqrt(2/length_nm)
        delta_s_list = df['|l_j|'].tolist()
        theta_list = df['theta'].tolist()
        s_mid_list = FourierShape.get_s_mid_list(df)
        summation = 0
        for delta_s, theta, s_mid in zip(delta_s_list, theta_list, s_mid_list):
            in_cos_term = n * np.pi / length_angstrom
            cos_term = np.cos(in_cos_term * s_mid)
            summation += delta_s * 0.1 * theta * cos_term
        return scale_factor * summation

    def get_bp_id_first_last(self):
        start_end = {'atat_21mer': (3, 17), 'g_tract_21mer': (3, 17), 
                     'a_tract_21mer': (3, 17), 'gcgc_21mer': (3, 17),
                     'ctct_21mer': (3, 17), 'tgtg_21mer': (3, 17),
                     'pnas_dna': (3, 12), 'pnas_rna': (3, 12)}
        return start_end[self.host][0], start_end[self.host][1]
   
    @staticmethod
    def get_s_mid_list(df):
        s_mid_list = np.zeros(df.shape[0])
        delta_s_list = df['|l_j|'].tolist()
        s_total = 0
        for i, delta_s in enumerate(delta_s_list):
            s_mid = s_total + delta_s/2
            s_total += delta_s
            s_mid_list[i] = s_mid
        return s_mid_list


class BigWindowPersistenceLength(FourierShape):

    def __init__(self, rootfolder, host):
        super().__init__(rootfolder, host)
        self.appr_length = 0.34 * (self.bp_id_last - self.bp_id_first) # Unit: nm
        self.columns = [f'n={mode_id}' for mode_id in self.mode_lst_str]
        self.df_lp = None

    def make_lp_big_wn_csv(self):
        df_an = self.get_df_an()
        d_result = self.get_lp_d_result(df_an)
        write_csv(self.lp_big_wn_csv, d_result, self.columns)

    def read_lp_big_wn_csv(self):
        self.df_lp = pd.read_csv(self.lp_big_wn_csv)
        print(f'Read {self.lp_big_wn_csv} into df_lp')

    def get_lp_d_result(self, df_an):
        d_result = {f'{key}': list() for key in self.columns}
        for mode_id in self.mode_lst_int:
            var_an = df_an[f'{mode_id}'].var()
            lp = BigWindowPersistenceLength.get_lp(self.appr_length, mode_id, var_an)
            d_result[f'n={mode_id}'].append(lp)
        return d_result

    def get_df_an(self):
        d_fourier_amp = read_d_result_from_hd5(self.fourier_amp_h5, self.mode_lst_str)
        return pd.DataFrame(d_fourier_amp)

    @staticmethod
    def get_lp(appr_length, mode_id, var_an):
        return np.square(appr_length) / (np.square(mode_id) * np.square(np.pi) * var_an)

class SmallWindowPersistenceLength(BigWindowPersistenceLength):
    def __init__(self, rootfolder, host, n_frames_per_window):
        super().__init__(rootfolder, host)
        self.n_frames_per_window = n_frames_per_window
        self.split_frame_list = None
        self.n_window = None
        self.columns = ['wn-id'] + [f'n={mode_id}' for mode_id in self.mode_lst_str]
        self.df_an_big_wn = self.get_df_an()

    def print_avg_std(self, mode_id):
        if self.df_lp is None:
            self.read_lp_small_wn_csv()
        avg = self.df_lp[f'n={mode_id}'].mean()
        std = self.df_lp[f'n={mode_id}'].std()
        print(f'Mode-ID: {mode_id}')
        print(f'Lp ={avg:.1f}±{std:.1f}')

    def make_lp_small_wn_csv(self):
        if self.split_frame_list is None:
            self.set_split_frame_list()
        d_result = {key: list() for key in self.columns}
        d_result['wn-id'] = list(range(self.n_window))
        for window_id in range(self.n_window):
            df_an_small_wn = self.get_df_an_small_wn(window_id)
            d_result_wn = self.get_lp_d_result(df_an_small_wn)
            for mode_id in self.mode_lst_int:
                d_result[f'n={mode_id}'].append(d_result_wn[f'n={mode_id}'][0])
        write_csv(self.lp_small_wn_csv, d_result, self.columns)

    def read_lp_small_wn_csv(self):
        self.df_lp = pd.read_csv(self.lp_small_wn_csv)
        print(f'Read {self.lp_small_wn_csv} into df_lp')

    def get_df_an_small_wn(self, window_id):
        frame_id_start, frame_id_end = self.split_frame_list[window_id]
        return self.df_an_big_wn.iloc[frame_id_start:frame_id_end+1]

    def set_split_frame_list(self):
        n_total_frames = self.get_last_frame()
        middle_interval = self.n_frames_per_window / 2
        split_frame_list = list()

        frame_id_1 = 0
        frame_id_2 = frame_id_1 + self.n_frames_per_window - 1

        execute_loop = True
        while execute_loop:
            if frame_id_2 > (n_total_frames + 1):
                execute_loop = False
                break
            split_frame_list.append((int(frame_id_1), int(frame_id_2)))
            frame_id_1 += middle_interval
            frame_id_2 = frame_id_1 + self.n_frames_per_window - 1
        self.split_frame_list = split_frame_list
        self.n_window = len(split_frame_list)
