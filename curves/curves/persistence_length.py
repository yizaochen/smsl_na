from os import path
from itertools import combinations
import MDAnalysis as mda
import numpy as np
import pandas as pd
from miscell.hd5_util import save_d_result_to_hd5, read_d_result_from_hd5
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
        self.mode_lst = list(range(self.n_begin, self.n_end+1))
        self.d_result = None
        self.bp_id_first, self.bp_id_last = self.get_bp_id_first_last()


    """
    def __init__(self, workfolder, host, bp_id_first=None, bp_id_last=None):
        self.host = host
        self.rootfolder = path.join(workfolder, host)
        self.df_folder = path.join(self.rootfolder, 'l_theta')
        self.an_folder = path.join(self.rootfolder, 'an_folder')
        self.n_bead = self.d_n_bp[host]
        self.bp_id_first = self.__set_bp_id_first(bp_id_first)
        self.bp_id_last = self.__set_bp_id_last(bp_id_last)
        self.df_name = path.join(self.df_folder, f'l_modulus_theta_{self.n_bead}_beads.csv')
        self.df = None

        self.__check_and_make_folders()
    """

    def make_fourier_amp_h5(self, last_frame):
        last_frame = self.get_last_frame()
        d_result = {mode_id: list() for mode_id in self.mode_lst}

        for frame_id in range(1, last_frame+1):
            df_filter = self.get_filter_df(frame_id)
            for mode_id in self.mode_lst:
                d_result[mode_id].append(self.get_an(mode_id, df_filter))
        self.d_result = d_result

        save_d_result_to_hd5(self.fourier_amp_h5, self.mode_lst, self.d_result)

    def read_df_an(self, n_begin, n_end):
        f_in = path.join(self.an_folder, f'an_{n_begin}_{n_end}_bpfirst{self.bp_id_first}_bplast{self.bp_id_last}.csv')
        df_an = pd.read_csv(f_in)
        return df_an

    def get_last_frame(self):
        u = mda.Universe(self.dcre_pdb, self.dcre_dcd)
        return len(u.trajectory)

    def get_filter_df(self, frame_id):
        mask = (self.df['i'] == self.bp_id_first)
        df0 = self.df[mask]
        mask = (df0['Frame_ID'] == frame_id)
        df1 = df0[mask]
        mask = (df1['j'].between(self.bp_id_first+1, self.bp_id_last))
        df2 = df1[mask]
        return df2

    def get_bp_id_first_last(self):
        start_end = {'atat_21mer': (3, 17), 'g_tract_21mer': (3, 17), 
                     'a_tract_21mer': (3, 17), 'gcgc_21mer': (3, 17),
                     'ctct_21mer': (3, 17), 'tgtg_21mer': (3, 17),
                     'pnas_dna': (3, 12), 'pnas_dna': (3, 12)}
        return start_end[self.host][0], start_end[self.host][1]
        

    """
    def get_appr_L(self):
        # Unit: nm
        return 0.34 * (self.bp_id_last - self.bp_id_first)
        
    def read_l_modulus_theta(self):
        self.df =  pd.read_csv(self.df_name)

    def get_L_of_frameid(self, frame_id):
        df_filter = self.get_filter_df(frame_id)
        L = self.__get_L(df_filter) # unit: angstrom
        return L

    def get_smid_and_interpolation_theta(self, frame_id):
        df_filter = self.get_filter_df(frame_id)
        theta_list = [0]
        theta_list += self.__get_theta_list(df_filter)
        n_theta = len(theta_list)
        s_mid_list = self.__get_s_mid_list(df_filter)
        interpolation_list = list()
        for i in range(n_theta-1):
            theta_inter = (theta_list[i] + theta_list[i+1]) / 2
            interpolation_list.append(theta_inter)
        return s_mid_list, interpolation_list


    def get_an(self, n, df_filter):
        L = self.__get_L(df_filter) # unit: angstrom
        L_nm = L / 10 # unit: nm
        delta_s_list = self.__get_delta_s_list(df_filter)
        theta_list = self.__get_theta_list(df_filter)
        s_mid_list = self.__get_s_mid_list(df_filter)
        scale_factor = np.sqrt(2/L_nm)
        summation = 0
        for delta_s, theta, s_mid in zip(delta_s_list, theta_list, s_mid_list):
            in_cos_term = n * np.pi / L
            cos_term = np.cos(in_cos_term * s_mid)
            summation += delta_s * 0.1 * theta * cos_term
        return scale_factor * summation

    def get_an_simplified(self, L, scale_factor, n, df_filter):
        delta_s_list = self.__get_delta_s_list(df_filter)
        theta_list = self.__get_theta_list(df_filter)
        s_mid_list = self.__get_s_mid_list(df_filter)
        summation = 0
        for delta_s, theta, s_mid in zip(delta_s_list, theta_list, s_mid_list):
            in_cos_term = n * np.pi / L
            cos_term = np.cos(in_cos_term * s_mid)
            summation += delta_s * 0.1 * theta * cos_term
        return scale_factor * summation

    def get_mode_shape_list(self, n, df_filter):
        L = self.__get_L(df_filter) # unit: angstrom
        L_nm = L / 10 # unit: nm
        s_list = np.array(self.__get_s_list(df_filter))
        scale_factor = np.sqrt(2/L_nm)
        in_cos_term = (n * np.pi * s_list) / L
        an = self.get_an(n, df_filter)
        cos_list = an * scale_factor * np.cos(in_cos_term)
        return s_list, cos_list

    def get_cos_list(self, s_list, L, scale_factor, n, df_filter):
        in_cos_term = (n * np.pi * s_list) / L
        an = self.get_an_simplified(L, scale_factor, n, df_filter)
        cos_list = an * scale_factor * np.cos(in_cos_term)
        return cos_list

    def get_cos_list_an(self, s_list, L, scale_factor, n, df_filter):
        in_cos_term = (n * np.pi * s_list) / L
        an = self.get_an_simplified(L, scale_factor, n, df_filter)
        cos_list = an * scale_factor * np.cos(in_cos_term)
        return cos_list, an

    def get_slist_thetalist(self, frame_id):
        df_filter = self.get_filter_df(frame_id)
        s_list = self.__get_s_list(df_filter)
        theta_list = self.__get_theta_list(df_filter)
        s_list = [0] + list(s_list)
        theta_list = [0] + theta_list
        return s_list, theta_list

    def get_approximate_theta(self, frame_id, n_begin, n_end):
        df_filter = self.get_filter_df(frame_id)
        L = self.__get_L(df_filter) # unit: angstrom
        L_nm = L / 10 # unit: nm
        scale_factor = np.sqrt(2/L_nm)
        s_list = self.__get_s_list(df_filter)
        appr_theta_list = np.zeros(len(s_list))
        for n in range(n_begin, n_end+1):
            cos_list = self.get_cos_list(s_list, L, scale_factor, n, df_filter)
            if n == 0:
                appr_theta_list += cos_list / 2
            else:
                appr_theta_list += cos_list
        s_list = [0] + list(s_list)
        appr_theta_list = [0] + list(appr_theta_list)
        return s_list, appr_theta_list

    def get_approximate_theta_singlemode(self, frame_id, n_select):
        df_filter = self.get_filter_df(frame_id)
        L = self.__get_L(df_filter) # unit: angstrom
        L_nm = L / 10 # unit: nm
        scale_factor = np.sqrt(2/L_nm)
        s_list = self.__get_s_list(df_filter)
        appr_theta_list = np.zeros(len(s_list))
        for n in [0, n_select]:
            cos_list, an = self.get_cos_list_an(s_list, L, scale_factor, n, df_filter)
            if n == 0:
                appr_theta_list += cos_list / 2
            else:
                appr_theta_list += cos_list
        s_list = [0] + list(s_list)
        appr_theta_list = [0] + list(appr_theta_list)
        return s_list, appr_theta_list, an

    

    

    def __get_L(self, df):
        return df['|l_j|'].sum()
        
    def __get_theta_list(self, df):
        return df['theta'].tolist()
        
    def __get_delta_s_list(self, df):
        return df['|l_j|'].tolist()
        
    def __get_s_mid_list(self, df):
        s_mid_list = np.zeros(df.shape[0])
        delta_s_list = df['|l_j|'].tolist()
        s_total = 0
        for i, delta_s in enumerate(delta_s_list):
            s_mid = s_total + delta_s/2
            s_total += delta_s
            s_mid_list[i] = s_mid
        return s_mid_list
        
    def __get_s_list(self, df):
        s_list = np.zeros(df.shape[0])
        delta_s_list = df['|l_j|'].tolist()
        s_total = 0
        for i, delta_s in enumerate(delta_s_list):
            s_total += delta_s
            s_list[i] = s_total
        return s_list

    def __check_and_make_folders(self):
        for folder in [self.an_folder]:
            check_dir_exist_and_make(folder)
    """
    