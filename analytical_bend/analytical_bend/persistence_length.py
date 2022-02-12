from os import path
from itertools import combinations
import miscell
import xyz_util
import re
import numpy as np
import pandas as pd
from MDAnalysis import Universe
import crd_util

pycharm_root = '/Users/yizao/PycharmProjects/connect_macro_micro'
ideal_root = '/Users/yizao/IdeaProjects/python/connect_macro_micro'
consist_folder = '/Users/yizao/PycharmProjects/consistency_quasi_fm'
seqeunces = {'pnas': {'arna+arna': {'guide': 'CAAUGGAGUA',
                                    'target': 'UACUCCAUUG'},
                      'bdna+bdna': {'guide': 'CAATGGAGTA',
                                    'target': 'TACTCCATTG'}},
             'pnas_amber': {'arna+arna': {'guide': 'CAAUGGAGUA',
                                          'target': 'UACUCCAUUG'},
                            'bdna+bdna': {'guide': 'CAATGGAGTA',
                                          'target': 'TACTCCATTG'}},
             'pnas_clean': {'arna+arna': {'guide': 'CAAUGGAGUA',
                                          'target': 'UACUCCAUUG'},
                            'bdna+bdna': {'guide': 'CAATGGAGTA',
                                          'target': 'TACTCCATTG'}},
             'pnas_amber_clean': {'arna+arna': {'guide': 'CAAUGGAGUA',
                                                'target': 'UACUCCAUUG'},
                                  'bdna+bdna': {'guide': 'CAATGGAGTA',
                                                'target': 'TACTCCATTG'}},
             'pnas_amber_16mer': {'arna+arna': {'guide': 'GCGCAAUGGAGUACGC', 'target': 'GCGUACUCCAUUGCGC'},
                                   'bdna+bdna': {'guide': 'GCGCAATGGAGTACGC', 'target': 'GCGTACTCCATTGCGC'}}
             }
c8_c6 = {'A': 'C8', 'T': 'C6', 'C': 'C6', 'G': 'C8', 'U': 'C6'}
x_axis = np.array([1, 0, 0])
y_axis = np.array([0, 1, 0])
z_axis = np.array([0, 0, 1])
d_l_avg = {'base': {'pnas_amber': {'arna+arna': 3.71, 'bdna+bdna': 3.48},
                    'pnas': {'arna+arna': 3.88, 'bdna+bdna': 3.52},
                    'pnas_amber_clean': {'arna+arna': 3.71, 'bdna+bdna': 3.48},
                    'pnas_amber_16mer': {'arna+arna': 3.71, 'bdna+bdna': 3.48},
                    'pnas_clean': {'arna+arna': 3.86, 'bdna+bdna': 3.52}},
           'pp': {'pnas_amber': {'arna+arna': 3.44, 'bdna+bdna': 3.45},
                  'pnas': {'arna+arna': 3.64, 'bdna+bdna': 3.52}}
           }
d_segid_abbr = {'STRAND1': 'STR1', 'STRAND2': 'STR2'}
d_mode = {'arna+arna': 1272, 'bdna+bdna': 1230}
d_natoms = {'arna+arna': 424, 'bdna+bdna': 410}


class CustomError(Exception):
    pass


class PersistenLengthAgent:
    def __init__(self, host, type_na, n_bp=10, n_atoms=None):
        self.host = host
        self.type_na = type_na
        if n_atoms is None:
            self.n_atoms = d_natoms[type_na]
        else:
            self.n_atoms = n_atoms
        if re.match('pnas_amber_clean', self.host):
            seq_host = 'pnas_amber_clean'
            self.sequence = seqeunces[seq_host][type_na]
        elif re.match('pnas_amber_16mer', self.host):
            seq_host = 'pnas_amber_16mer'
            self.sequence = seqeunces[seq_host][type_na]
        else:
            self.sequence = seqeunces[self.host][type_na]
        self.n_bp = n_bp  # Number of Base-Pair
        self.pyc_na_folder = path.join(pycharm_root, host, type_na)
        self.ideal_input_folder = path.join(ideal_root, host, type_na, 'input', 'heavyatoms')
        self.avg_crd = path.join(self.ideal_input_folder, '{0}.nohydrogen.avg.crd'.format(type_na))
        self.fit_dcd = path.join(self.ideal_input_folder, '{0}.nohydrogen.fitavg.dcd'.format(type_na))
        # Contour Length and end-to-end length folder defined by phosphate
        self.lc_lee_folder = path.join(self.pyc_na_folder, 'lc_lee_pp')
        self.midpoints_folder = path.join(self.pyc_na_folder, 'midpoints')
        self.ee_align_z_folder = path.join(self.pyc_na_folder, 'ee_align_z')
        self.xyz_traj_folder = path.join(self.pyc_na_folder, 'xyz_traj')
        self.l_theta_phi_folder = path.join(self.pyc_na_folder, 'l_theta_phi')
        self.result_folder = path.join(self.pyc_na_folder, 'pl_result')
        self.cosdelta_folder = path.join(self.pyc_na_folder, 'cos_delta')
        self.smallr_folder = path.join(self.pyc_na_folder, 'small_r')

        self.make_folders()

    def make_folders(self):
        folders = [self.lc_lee_folder, self.midpoints_folder, self.ee_align_z_folder, self.xyz_traj_folder,
                   self.l_theta_phi_folder, self.result_folder, self.cosdelta_folder, self.smallr_folder]
        for folder in folders:
            miscell.check_dir_exist_and_make(folder)

    def get_lc_lee_of_md(self, rtype='base', n_bead=10):
        u = Universe(self.avg_crd, self.fit_dcd)
        if rtype == 'pp':
            basepairs = self.get_basepairs_pp(u)
            f_out = path.join(self.lc_lee_folder, 'md_pp_nbead_{0}.csv'.format(n_bead))
        else:
            basepairs = self.get_basepairs(u)
            f_out = path.join(self.lc_lee_folder, 'md_nbead_{0}.csv'.format(n_bead))
        d = {'Frame': list(), 'Lc': list(), 'Lee': list()}
        for ts in u.trajectory[:]:
            midpoints = self.get_midpoints(basepairs, n_bead=n_bead)
            lc = self.get_distance_sum(midpoints, n_bead=n_bead)
            lee = self.get_end_to_end_phosphate(midpoints, n_bead=n_bead)
            d['Frame'].append(ts.frame + 1)
            d['Lc'].append(lc)
            d['Lee'].append(lee)
        df = pd.DataFrame(d)
        df.to_csv(f_out, index=False)
        return df

    def read_lc_lee_of_md(self, rtype='base', n_bead=10):
        if rtype == 'pp':
            f_in = path.join(self.lc_lee_folder, 'md_pp_nbead_{0}.csv'.format(n_bead))
        else:
            f_in = path.join(self.lc_lee_folder, 'md_nbead_{0}.csv'.format(n_bead))
        df = pd.read_csv(f_in)
        return df

    def get_lc_lee_of_avgcrd(self, rtype='base', n_bead=10):
        """

        :param rtype: pp, base
        :param n_bead: 3,4,5,6,7,8,9,10
        :return:
        """
        u = Universe(self.avg_crd, self.avg_crd)
        if rtype == 'pp':
            basepairs = self.get_basepairs_pp(u)
            f_out = path.join(self.lc_lee_folder, 'avg_pp_nbead_{0}.csv'.format(n_bead))
        else:
            basepairs = self.get_basepairs(u)
            f_out = path.join(self.lc_lee_folder, 'avg_nbead_{0}.csv'.format(n_bead))
        d = {'Frame': list(), 'Lc': list(), 'Lee': list()}
        midpoints = self.get_midpoints(basepairs, n_bead=n_bead)
        lc = self.get_distance_sum(midpoints, n_bead=n_bead)
        lee = self.get_end_to_end_phosphate(midpoints, n_bead=n_bead)
        d['Frame'].append('Avg')
        d['Lc'].append(lc)
        d['Lee'].append(lee)
        df = pd.DataFrame(d)
        df.to_csv(f_out, index=False)
        return df

    def read_lc_lee_of_avgcrd(self, rtype='base', n_bead=10):
        if rtype == 'pp':
            f_in = path.join(self.lc_lee_folder, 'avg_pp_nbead_{0}.csv'.format(n_bead))
        else:
            f_in = path.join(self.lc_lee_folder, 'avg_nbead_{0}.csv'.format(n_bead))
        df = pd.read_csv(f_in)
        return df

    def get_mean_theta_from_avgcrd(self, i=0, j=8, n_bead=10):
        u = Universe(self.avg_crd, self.avg_crd)
        basepairs = self.get_basepairs(u)
        midpoints = self.get_midpoints(basepairs, n_bead=n_bead)
        h_ids = get_h_id_by_i_j(i, j)
        l1 = midpoints[h_ids[1]] - midpoints[h_ids[0]]
        l2 = midpoints[h_ids[3]] - midpoints[h_ids[2]]
        l1_norm = np.linalg.norm(l1)
        l2_norm = np.linalg.norm(l2)
        dot_product = np.dot(l1, l2)
        mean_theta = np.arccos(dot_product / (l1_norm * l2_norm))
        return mean_theta

    def make_lee_lc_table_by_bpids(self, bpid1=1, bpid2=2):
        u = Universe(self.avg_crd, self.fit_dcd)
        basepairs = self.get_basepairs(u)
        d = {'Frame': list(), 'Lc': list(), 'Lee': list()}
        for ts in u.trajectory[:]:
            midpoints = self.get_midpoints(basepairs)
            lc = get_distance_sum_by_bpids(midpoints, bpid1, bpid2)
            lee = get_lee_by_bpids(midpoints, bpid1, bpid2)
            d['Frame'].append(ts.frame + 1)
            d['Lc'].append(lc)
            d['Lee'].append(lee)
        df = pd.DataFrame(d)
        f_out = path.join(self.lc_lee_folder, 'lc_lee_resid{0}_resid{1}.csv'.format(bpid1, bpid2))
        df.to_csv(f_out, index=False)
        return df

    def make_persistence_length_table(self, first_bpid=2, last_bpid=10):
        bpid1 = 1
        i = 1
        d_result = {'i': list(), '<R^2>': list(), '<L>': list(), '<L>^2': list(), 'Dr(L)': list()}
        for bpid2 in range(first_bpid, last_bpid):
            f_in = path.join(self.lc_lee_folder, 'lc_lee_resid{0}_resid{1}.csv'.format(bpid1, bpid2))
            df = pd.read_csv(f_in)
            lees = df['Lee'].tolist()
            lee_square_mean = np.square(lees)
            lee_square_mean = lee_square_mean.mean()
            d_result['<R^2>'].append(lee_square_mean)
            lc_mean = df['Lc'].mean()
            d_result['<L>'].append(lc_mean)
            lc_mean_square = np.square(lc_mean)
            d_result['<L>^2'].append(lc_mean_square)
            dr_l = 3 * (1 - (lee_square_mean/lc_mean_square))
            d_result['Dr(L)'].append(dr_l)
            d_result['i'].append(i)
            i += 1
        df = pd.DataFrame(d_result)
        f_out = path.join(self.lc_lee_folder, 'preliminary_table_pl.csv')
        df.to_csv(f_out, index=False)
        return df

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

    def get_basepairs_pp(self, u):
        basepairs = {}
        bp_id = 1
        for resid_i in range(1, self.n_bp + 1):
            # resid_i: guide strand resid  # resid_j: target strand resid
            resid_j = miscell.get_antistrand_resid(resid_i)
            temp = (u.select_atoms('segid STRAND1 and resid {0} and name P'.format(resid_i)),
                    u.select_atoms('segid STRAND2 and resid {0} and name P'.format(resid_j)))
            basepairs[bp_id] = temp
            bp_id += 1
        return basepairs

    def get_midpoints(self, basepairs, n_bead=10):
        if n_bead == 10:
            n_bp = self.n_bp
        else:
            n_bp = n_bead
        midpoints = dict()
        for bp_id in range(1, n_bp + 1):
            atom1 = basepairs[bp_id][0]
            atom2 = basepairs[bp_id][1]
            midpoint = (atom1.positions[0] + atom2.positions[0]) / 2
            midpoints[bp_id] = midpoint
        return midpoints

    def get_distance_sum(self, midpoints, n_bead=10):
        if n_bead == 10:
            n_bp = self.n_bp
        else:
            n_bp = n_bead
        h_array = [np.linalg.norm(midpoints[bp_id] - midpoints[bp_id + 1])
                   for bp_id in range(1, n_bp)]
        return sum(h_array)

    def get_end_to_end_phosphate(self, points, n_bead=10):
        if n_bead == 10:
            n_bp = self.n_bp
        else:
            n_bp = n_bead
        return np.linalg.norm(points[1] - points[n_bp])

    def make_avgcrd_midpoints_npy(self):
        u = Universe(self.avg_crd, self.avg_crd)
        basepairs = self.get_basepairs(u)
        f_out = path.join(self.midpoints_folder, 'midpoints_avgcrd.npy')
        midpoints = self.get_midpoints(basepairs)
        temp_list = [midpoints[bp_id] for bp_id in range(1, self.n_bp+1)]
        result_list = np.array(temp_list)
        np.save(f_out, result_list)
        return result_list

    def read_avgcrd_midpoints_npy(self):
        f_in = path.join(self.midpoints_folder, 'midpoints_avgcrd.npy')
        midpoints = np.load(f_in)
        return midpoints

    def make_allmidpoints_into_npy(self, rtype='base'):
        u = Universe(self.avg_crd, self.fit_dcd)
        if rtype == 'pp':
            basepairs = self.get_basepairs_pp(u)
            f_out = path.join(self.midpoints_folder, 'midpoints_allframes_pp.npy')
        else:
            basepairs = self.get_basepairs(u)
            f_out = path.join(self.midpoints_folder, 'midpoints_allframes.npy')
        result_list = list()
        for ts in u.trajectory[:]:
            midpoints = self.get_midpoints(basepairs)
            temp_list = [midpoints[bp_id] for bp_id in range(1, self.n_bp+1)]
            result_list.append(temp_list)
        result_list = np.array(result_list)
        np.save(f_out, result_list)
        return result_list

    def make_npy_align_ee_to_zaxis(self, n_bead=3, rtype='base'):
        if rtype == 'pp':
            f_in = path.join(self.midpoints_folder, 'midpoints_allframes_pp.npy')
            f_out = path.join(self.ee_align_z_folder, 'align_points_{0}_beads_pp.npy'.format(n_bead))
        else:
            f_in = path.join(self.midpoints_folder, 'midpoints_allframes.npy')
            f_out = path.join(self.ee_align_z_folder, 'align_points_{0}_beads.npy'.format(n_bead))
        midpoints_allframe = np.load(f_in)
        result_list = list()
        for frame in midpoints_allframe:
            temp = list()
            for bead_idx in range(n_bead):
                temp.append(frame[bead_idx])
            temp = np.array(temp)
            new_points = align_ee_to_zaxis(temp)
            result_list.append(new_points)
        result_list = np.array(result_list)
        np.save(f_out, result_list)
        return result_list

    def read_npy_align_ee_to_zaxis(self, n_bead=3, rtype='base'):
        if rtype == 'pp':
            f_in = path.join(self.ee_align_z_folder, 'align_points_{0}_beads_pp.npy'.format(n_bead))
        else:
            f_in = path.join(self.ee_align_z_folder, 'align_points_{0}_beads.npy'.format(n_bead))
        result_list = np.load(f_in)
        return result_list

    def align_md_allthreepoints_npy(self, v1_id='01', v2_id='89'):
        result_list = list()
        u = Universe(self.avg_crd, self.avg_crd)
        basepairs = self.get_basepairs(u)
        midpoints = self.get_midpoints(basepairs)
        temp_list = [midpoints[bp_id] for bp_id in range(1, self.n_bp + 1)]
        threepoints = translate_two_heads_to_origin(temp_list, v1_id=v1_id, v2_id=v2_id)
        threepoints_xy = align_cross_of_v1v2_to_zaxis(threepoints)
        threepoints_final = align_v1_to_xaxis(threepoints_xy)
        result_list.append(threepoints_final)

        u = Universe(self.avg_crd, self.fit_dcd)
        basepairs = self.get_basepairs(u)
        f_out = path.join(self.midpoints_folder, 'threepoints_l{0}_l{1}_allframes.npy'.format(v1_id, v2_id))
        for ts in u.trajectory[:]:
            midpoints = self.get_midpoints(basepairs)
            temp_list = [midpoints[bp_id] for bp_id in range(1, self.n_bp + 1)]
            threepoints = translate_two_heads_to_origin(temp_list, v1_id=v1_id, v2_id=v2_id)
            threepoints_xy = align_cross_of_v1v2_to_zaxis(threepoints)
            threepoints_final = align_v1_to_xaxis(threepoints_xy)
            result_list.append(threepoints_final)
        result_list = np.array(result_list)
        np.save(f_out, result_list)
        return result_list

    def read_md_allthreepoints_npy(self, v1_id='01', v2_id='89'):
        f_in = path.join(self.midpoints_folder, 'threepoints_l{0}_l{1}_allframes.npy'.format(v1_id, v2_id))
        result_list = np.load(f_in)
        return result_list

    def get_align_threepoints_arbitrary_crd(self, f_crd, v1_id='01', v2_id='89'):
        u = Universe(f_crd, f_crd)
        basepairs = self.get_basepairs(u, short_segid=True)
        midpoints = self.get_midpoints(basepairs)
        temp_list = [midpoints[bp_id] for bp_id in range(1, self.n_bp + 1)]
        threepoints = translate_two_heads_to_origin(temp_list, v1_id=v1_id, v2_id=v2_id)
        threepoints_xy = align_cross_of_v1v2_to_zaxis(threepoints)
        threepoints_final = align_v1_to_xaxis(threepoints_xy)
        return threepoints_final

    def write_traj_xyz(self, n_bead, points, rtype='base'):
        if rtype == 'pp':
            f_out = path.join(self.xyz_traj_folder, '{0}_beads_pp.xyz'.format(n_bead))
        else:
            f_out = path.join(self.xyz_traj_folder, '{0}_beads.xyz'.format(n_bead))
        xyz_util.write_trajectory(points, f_out)
        print('Write coordinates into {0}'.format(f_out))

    def get_l_modulus_theta(self, n_bead=3, rtype='base'):
        points_allframes = self.read_npy_align_ee_to_zaxis(n_bead, rtype=rtype)
        d_result = {'Frame_ID': list(), 'i': list(), 'j': list(), '|l_i|': list(), '|l_j|': list(), 'theta': list()}
        pair_list = get_pair_list(n_bead)
        frame_id = 1
        for points in points_allframes:
            vectors = [points[i + 1] - points[i] for i in range(len(points) - 1)]
            for i, j in pair_list:
                vi_modulus, vj_modulus, theta = miscell.get_modulus_angle_between_two_vectors(vectors[i], vectors[j])
                d_result['Frame_ID'].append(frame_id)
                d_result['i'].append(i)
                d_result['j'].append(j)
                d_result['|l_i|'].append(vi_modulus)
                d_result['|l_j|'].append(vj_modulus)
                d_result['theta'].append(theta)
            frame_id += 1
        df = pd.DataFrame(d_result)
        df = df[['Frame_ID', 'i', 'j', '|l_i|', '|l_j|', 'theta']]
        if rtype == 'pp':
            f_out = path.join(self.l_theta_phi_folder, 'modulus_theta_{0}_beads_pp.csv'.format(n_bead))
        else:
            f_out = path.join(self.l_theta_phi_folder, 'modulus_theta_{0}_beads.csv'.format(n_bead))
        df.to_csv(f_out, index=False)
        return df

    def read_l_modulus_theta(self, n_bead=3, rtype='base'):
        if rtype == 'pp':
            f_in = path.join(self.l_theta_phi_folder, 'modulus_theta_{0}_beads_pp.csv'.format(n_bead))
        else:
            f_in = path.join(self.l_theta_phi_folder, 'modulus_theta_{0}_beads.csv'.format(n_bead))
        df = pd.read_csv(f_in)
        return df

    def get_statistical_parameters_from_df(self, i, j, n_bead=3):
        f_in = path.join(self.l_theta_phi_folder, 'modulus_theta_{0}_beads.csv'.format(n_bead))
        df = pd.read_csv(f_in)
        d_result = {'i': list(), 'j': list(), '<l_i^2>': list(), '<l_j^2>': list(),
                    'cov(|l_i||l_j|,cos(theta+delta))': list(), '<|l_i||l_j|>': list(), 'theta_mean': list(),
                    '<sin(delta)>': list()}
        d_result['i'].append(i)
        d_result['j'].append(j)
        d_result['<l_i^2>'].append(np.square(df['|l_i|']).mean())
        d_result['<l_j^2>'].append(np.square(df['|l_j|']).mean())
        lilj = np.multiply(df['|l_i|'], df['|l_j|'])
        cos_thetaplusdelta = np.cos(df['theta'])
        temp = np.array([lilj, cos_thetaplusdelta])
        cov_mat = np.cov(temp)
        d_result['cov(|l_i||l_j|,cos(theta+delta))'].append(cov_mat[0, 1])
        d_result['<|l_i||l_j|>'].append(lilj.mean())
        theta_mean = df['theta'].mean()
        d_result['theta_mean'].append(theta_mean)
        delta = df['theta'] - theta_mean
        d_result['<sin(delta)>'].append(np.sin(delta).mean())
        df = pd.DataFrame(d_result)
        df = df[['i', 'j', '<l_i^2>', '<l_j^2>', 'cov(|l_i||l_j|,cos(theta+delta))', '<|l_i||l_j|>', 'theta_mean',
                 '<sin(delta)>']]
        f_out = path.join(self.l_theta_phi_folder, 'sta_parameters_{0}_beads.csv'.format(n_bead))
        df.to_csv(f_out, index=False)
        return df

    def get_mean_ee_square(self, n_bead=3):
        """
        Get <R^2>
        :param n_bead:
        :return:
        """
        points_allframes = self.read_npy_align_ee_to_zaxis(n_bead)
        r2_container = [np.linalg.norm(points[0]-points[-1]) for points in points_allframes]
        r2_container = np.square(r2_container)
        return r2_container, r2_container.mean()

    def get_l_theta_phi_allframes(self, n_bead=3):
        points_allframes = self.read_npy_align_ee_to_zaxis(n_bead)
        d_result = {'Frame_ID': list(), 'vector_ID': list(), 'l': list(), 'theta': list(), 'phi': list(),
                    't_x': list(), 't_y': list(), 't_z': list()}
        frame_id = 1
        for points in points_allframes:
            vectors = [points[i + 1] - points[i] for i in range(len(points) - 1)]
            v_id = 1
            for vector in vectors:
                d_result['Frame_ID'].append(frame_id)
                d_result['vector_ID'].append(v_id)
                length = np.linalg.norm(vector)
                d_result['l'].append(length)
                t_x, t_y, t_z, theta, phi = get_tangent_theta_phi(vector, length)
                d_result['theta'].append(theta)
                d_result['phi'].append(phi)
                d_result['t_x'].append(t_x)
                d_result['t_y'].append(t_y)
                d_result['t_z'].append(t_z)
                v_id += 1
            frame_id += 1
        df = pd.DataFrame(d_result)
        df = df[['Frame_ID', 'vector_ID', 'l', 'theta', 'phi', 't_x', 't_y', 't_z']]
        f_out = path.join(self.l_theta_phi_folder, 'l_phi_theta_tangent_{0}_beads.csv'.format(n_bead))
        df.to_csv(f_out, index=False)
        return df

    def make_ee2_vector(self, n_bead=10, rtype='base'):
        points_allframes = self.read_npy_align_ee_to_zaxis(n_bead=n_bead, rtype=rtype)
        d_result = dict()
        for i in range(n_bead):
            d_result[i] = [np.linalg.norm(points[i]-points[0]) for points in points_allframes]
        l_result = [np.square(d_result[i]).mean() for i in range(n_bead)]
        l_result = np.array(l_result)
        if rtype == 'pp':
            f_out = path.join(self.l_theta_phi_folder, 'mean_of_R_square_pp.npy')
        else:
            f_out = path.join(self.l_theta_phi_folder, 'mean_of_R_square.npy')
        np.save(f_out, l_result)
        return l_result

    def get_ee2(self, n_bead=3, rtype='base'):
        if rtype == 'pp':
            f_in = path.join(self.l_theta_phi_folder, 'mean_of_R_square_pp.npy')
        else:
            f_in = path.join(self.l_theta_phi_folder, 'mean_of_R_square.npy')
        result_list = np.load(f_in)
        return result_list[n_bead-1]

    def make_ee2_matrix(self, n_bead=10, rtype='base'):
        """
        Make <R^2>_ij
        :param n_bead:
        :param rtype:
        :return:
        """
        points_allframes = self.read_npy_align_ee_to_zaxis(n_bead=n_bead, rtype=rtype)
        ee2_mat = np.zeros((n_bead, n_bead))
        pair_list = get_pair_list(n_bead+1)
        for i, j in pair_list:
            temp = [np.linalg.norm(points[i]-points[j]) for points in points_allframes]
            ee2_mat[i, j] = np.square(temp).mean()
            ee2_mat[j, i] = ee2_mat[i, j]
        if rtype == 'pp':
            f_out = path.join(self.l_theta_phi_folder, 'mean_of_R_square_mat_pp.npy')
        else:
            f_out = path.join(self.l_theta_phi_folder, 'mean_of_R_square_mat.npy')
        np.save(f_out, ee2_mat)
        return ee2_mat

    def get_ee2_matrix(self, rtype='base'):
        if rtype == 'pp':
            f_in = path.join(self.l_theta_phi_folder, 'mean_of_R_square_mat_pp.npy')
        else:
            f_in = path.join(self.l_theta_phi_folder, 'mean_of_R_square_mat.npy')
        result_list = np.load(f_in)
        return result_list

    def get_ee2_from_ee2mat(self, i, j, rtype='base'):
        ee2_mat = self.get_ee2_matrix(rtype=rtype)
        return ee2_mat[i, j]

    def make_mean_l2_vector(self, n_bead=10, rtype='base'):
        points_allframes = self.read_npy_align_ee_to_zaxis(n_bead=n_bead, rtype=rtype)
        d_result = dict()
        n_vector = n_bead - 1
        for i in range(n_vector):
            d_result[i] = [np.linalg.norm(points[i]-points[i+1]) for points in points_allframes]
        l_result = [np.square(d_result[i]).mean() for i in range(n_vector)]
        l_result = np.array(l_result)
        if rtype == 'pp':
            f_out = path.join(self.l_theta_phi_folder, 'mean_l2_square_pp.npy')
        else:
            f_out = path.join(self.l_theta_phi_folder, 'mean_l2_square.npy')
        np.save(f_out, l_result)
        return l_result

    def get_mean_l2(self, l_id, rtype='base'):
        if rtype == 'pp':
            f_in = path.join(self.l_theta_phi_folder, 'mean_l2_square_pp.npy')
        else:
            f_in = path.join(self.l_theta_phi_folder, 'mean_l2_square.npy')
        result_list = np.load(f_in)
        return result_list[l_id]

    def make_cov_term_matrix(self, n_bead=10, rtype='base'):
        if rtype == 'pp':
            f_in = path.join(self.l_theta_phi_folder, 'modulus_theta_{0}_beads.csv'.format(n_bead))
            f_out = path.join(self.l_theta_phi_folder, 'cov_term_mat_pp.npy')
        else:
            f_in = path.join(self.l_theta_phi_folder, 'modulus_theta_{0}_beads.csv'.format(n_bead))
            f_out = path.join(self.l_theta_phi_folder, 'cov_term_mat.npy')
        df = pd.read_csv(f_in)
        n_vector = n_bead - 1
        cov_result = np.zeros((n_vector, n_vector))
        pair_list = get_pair_list(n_bead)
        for i, j in pair_list:
            mask = ((df['i'] == i) & (df['j'] == j))
            df_temp = df[mask]
            lilj = np.multiply(df_temp['|l_i|'], df_temp['|l_j|'])
            cos_thetaplusdelta = np.cos(df_temp['theta'])
            temp = np.array([lilj, cos_thetaplusdelta])
            cov_mat = np.cov(temp)
            cov_result[i, j] = cov_mat[0, 1]
            cov_result[j, i] = cov_result[i, j]
        np.save(f_out, cov_result)
        return cov_result

    def read_cov_term_matrix(self, rtype='base'):
        if rtype == 'pp':
            f_in = path.join(self.l_theta_phi_folder, 'cov_term_mat_pp.npy')
        else:
            f_in = path.join(self.l_theta_phi_folder, 'cov_term_mat.npy')
        cov_term_mat = np.load(f_in)
        return cov_term_mat

    def make_aij_matrix(self, n_bead=10, rtype='base'):
        if rtype == 'pp':
            f_in = path.join(self.l_theta_phi_folder, 'modulus_theta_{0}_beads_pp.csv'.format(n_bead))
            f_out = path.join(self.l_theta_phi_folder, 'a_mat_pp.npy')
        else:
            f_in = path.join(self.l_theta_phi_folder, 'modulus_theta_{0}_beads.csv'.format(n_bead))
            f_out = path.join(self.l_theta_phi_folder, 'a_mat.npy')
        df = pd.read_csv(f_in)
        n_vector = n_bead - 1
        a_mat = np.zeros((n_vector, n_vector))
        pair_list = get_pair_list(n_bead)
        for i, j in pair_list:
            mask = ((df['i'] == i) & (df['j'] == j))
            df_temp = df[mask]
            lilj = np.multiply(df_temp['|l_i|'], df_temp['|l_j|'])
            a_mat[i, j] = lilj.mean()
            a_mat[j, i] = a_mat[i, j]
        np.save(f_out, a_mat)
        return a_mat

    def read_aij_matrix(self, rtype='base'):
        if rtype == 'pp':
            f_in = path.join(self.l_theta_phi_folder, 'a_mat_pp.npy')
        else:
            f_in = path.join(self.l_theta_phi_folder, 'a_mat.npy')
        result_mat = np.load(f_in)
        return result_mat

    def make_st_sd_ct_matrix(self, n_bead=10, rtype='base'):
        if rtype == 'pp':
            f_in = path.join(self.l_theta_phi_folder, 'modulus_theta_{0}_beads_pp.csv'.format(n_bead))
            st_out = path.join(self.l_theta_phi_folder, 'st_mat_pp.npy')
            sd_out = path.join(self.l_theta_phi_folder, 'sd_mat_pp.npy')
            ct_out = path.join(self.l_theta_phi_folder, 'ct_mat_pp.npy')
        else:
            f_in = path.join(self.l_theta_phi_folder, 'modulus_theta_{0}_beads.csv'.format(n_bead))
            st_out = path.join(self.l_theta_phi_folder, 'st_mat.npy')
            sd_out = path.join(self.l_theta_phi_folder, 'sd_mat.npy')
            ct_out = path.join(self.l_theta_phi_folder, 'ct_mat.npy')
        df = pd.read_csv(f_in)
        n_vector = n_bead - 1
        st = np.zeros((n_vector, n_vector))
        sd = np.zeros((n_vector, n_vector))
        ct = np.zeros((n_vector, n_vector))
        pair_list = get_pair_list(n_bead)
        for i, j in pair_list:
            mask = ((df['i'] == i) & (df['j'] == j))
            df_temp = df[mask]
            theta_mean = df_temp['theta'].mean()
            delta = df_temp['theta'] - theta_mean
            st[i, j] = np.sin(theta_mean)
            st[j, i] = st[i, j]
            sd[i, j] = np.sin(delta).mean()
            sd[j, i] = sd[i, j]
            ct[i, j] = np.cos(theta_mean)
            ct[j, i] = ct[i, j]
        np.save(st_out, st)
        np.save(sd_out, sd)
        np.save(ct_out, ct)
        return st, sd, ct

    def read_st_sd_ct_matrix(self, rtype='base'):
        if rtype == 'pp':
            st_in = path.join(self.l_theta_phi_folder, 'st_mat_pp.npy')
            sd_in = path.join(self.l_theta_phi_folder, 'sd_mat_pp.npy')
            ct_in = path.join(self.l_theta_phi_folder, 'ct_mat_pp.npy')
        else:
            st_in = path.join(self.l_theta_phi_folder, 'st_mat.npy')
            sd_in = path.join(self.l_theta_phi_folder, 'sd_mat.npy')
            ct_in = path.join(self.l_theta_phi_folder, 'ct_mat.npy')
        st = np.load(st_in)
        sd = np.load(sd_in)
        ct = np.load(ct_in)
        return st, sd, ct

    def get_coefficient_and_constant(self, n_bead):
        # Get the first term, <|l_i|^2>
        first_term = 0
        for i in range(n_bead-1):
            first_term += self.get_mean_l2(i)

        # Get the second term
        pair_list = get_pair_list(n_bead)
        cov_mat = self.read_cov_term_matrix()
        a_mat = self.read_aij_matrix()
        st, sd, ct = self.read_st_sd_ct_matrix()
        second_term = 0
        for i, j in pair_list:
            second_term += 2 * cov_mat[i, j]
            second_term -= (2 * a_mat[i, j] * st[i, j] * sd[i, j])

        # Get the constant
        ee_2 = self.get_ee2(n_bead)
        c = first_term + second_term - ee_2

        # Get the coefficients
        coefficients = dict()
        for i in range(n_bead-2):
            coefficients[i+1] = 0
        for i, j in pair_list:
            key = j - i
            temp = 2 * a_mat[i, j] * ct[i, j]
            coefficients[key] += temp
        return coefficients, c

    def get_ensemble_average_l(self, n_bead=10, rtype='base'):
        points_allframes = self.read_npy_align_ee_to_zaxis(n_bead=n_bead, rtype=rtype)
        l_result = list()
        n_vector = n_bead - 1
        for i in range(n_vector):
            for points in points_allframes:
                l_result.append(np.linalg.norm(points[i]-points[i+1]))
        l_result = np.array(l_result)
        return l_result.mean()

    def get_coefficient_and_constant_by_ij(self, i, j, rtype='base'):
        """
        Calculate coefficient from bead_i and bead_j
        :param i: id of bead_i
        :param j: id of bead_j
        :param rtype: pp or base
        :return:
        """
        pair_list, vector_list = get_pair_list_by_beadid(i, j)

        # Get the first term, <|l_i|^2>
        first_term = 0
        for v_id in vector_list:
            first_term += self.get_mean_l2(v_id, rtype=rtype)

        # Get the second term
        cov_mat = self.read_cov_term_matrix(rtype=rtype)
        a_mat = self.read_aij_matrix(rtype=rtype)
        st, sd, ct = self.read_st_sd_ct_matrix(rtype=rtype)
        second_term = 0
        for v_i, v_j in pair_list:
            second_term += 2 * cov_mat[v_i, v_j]
            second_term -= (2 * a_mat[v_i, v_j] * st[v_i, v_j] * sd[v_i, v_j])

        # Get the constant
        ee_2 = self.get_ee2_from_ee2mat(i, j, rtype=rtype)
        c = first_term + second_term - ee_2

        # Get the coefficients
        coefficients = dict()
        n_bead = j - i + 1
        for i in range(n_bead-2):
            coefficients[i+1] = 0
        for v_i, v_j in pair_list:
            key = v_j - v_i
            temp = 2 * a_mat[v_i, v_j] * ct[v_i, v_j]
            coefficients[key] += temp
        return coefficients, c

    def write_pl_mean_err(self, data, n_bead_min, n_bead_max, rtype='base'):
        x_list = range(n_bead_min, n_bead_max + 1)
        d_result = {'n_bead': list(), 'mean': list(), 'std': list()}
        for i in x_list:
            d_result['n_bead'].append(i)
            temp = np.array(data[i])
            d_result['mean'].append(temp.mean())
            if len(temp) >= 3:
                d_result['std'].append(temp.std())
            else:
                d_result['std'].append(np.nan)
        df = pd.DataFrame(d_result)
        df = df[['n_bead', 'mean', 'std']]
        if rtype == 'pp':
            f_out = path.join(self.result_folder, 'pl_result_pp.csv')
        else:
            f_out = path.join(self.result_folder, 'pl_result.csv')
        df.to_csv(f_out, index=False)
        return df

    def get_pl_mean_err(self, rtype='base'):
        if rtype == 'pp':
            f_in = path.join(self.result_folder, 'pl_result_pp.csv')
        else:
            f_in = path.join(self.result_folder, 'pl_result.csv')
        df = pd.read_csv(f_in)
        x_list = df['n_bead'].tolist()
        mean_list = df['mean'].tolist()
        err_list = df['std'].tolist()
        return x_list, mean_list, err_list

    def get_mean_theta_cosdelta_array_by_ij(self, i=0, j=8, n_bead=10):
        df = self.read_l_modulus_theta(n_bead=n_bead)
        mask = (df['i'] == i) & (df['j'] == j)
        df_temp = df[mask]
        #  mean_theta = df_temp['theta'].mean()
        mean_theta = self.get_mean_theta_from_avgcrd(i, j, n_bead)
        delta_theta_array = df_temp['theta'].values - mean_theta
        cos_delta_theta = np.cos(delta_theta_array)
        cos_theta_array = np.cos(df_temp['theta'].values)
        return mean_theta, cos_delta_theta, cos_theta_array, df_temp['theta'].values

    def write_statistics_for_diff_ij(self, n_bead=10):
        df_data = self.read_l_modulus_theta(n_bead=n_bead)
        pair_list = get_pair_list(n_bead=n_bead)
        d_result = {'i': list(), 'j': list(), 'pl': list()}
        d_result_1 = dict()
        for i, j in pair_list:
            mean_theta = self.get_mean_theta_from_avgcrd(i, j, n_bead)
            pl, j_minus_i = get_persistence_length(i, j, df_data, mean_theta, self.host, self.type_na)
            d_result['i'].append(i)
            d_result['j'].append(j)
            d_result['pl'].append(pl)
            if j_minus_i in d_result_1:
                d_result_1[j_minus_i].append(pl)
            else:
                d_result_1[j_minus_i] = [pl]
        f_out_df = path.join(self.result_folder, 'pl_statistic.csv')
        df = pd.DataFrame(d_result)
        df = df[['i', 'j', 'pl']]
        df.to_csv(f_out_df, index=False)
        f_out_dict = path.join(self.result_folder, 'pl_statistic_dict.pkl')
        miscell.save_dict_to_pkl(d_result_1, f_out_dict)
        return df, d_result_1

    def read_statistics_for_diff_ij(self):
        f_in_df = path.join(self.result_folder, 'pl_statistic.csv')
        f_in_dict = path.join(self.result_folder, 'pl_statistic_dict.pkl')
        d_result = miscell.load_pkl_from_dict(f_in_dict)
        df = pd.read_csv(f_in_df)
        return df, d_result


class ThetaVarianceAgent(PersistenLengthAgent):
    def __init__(self, host, type_na, n_bp=10):
        super(ThetaVarianceAgent, self).__init__(host, type_na, n_bp)
        self.derivative_folder = path.join(self.pyc_na_folder, 'derivative_to_c')
        self.nm_single_folder = path.join(self.pyc_na_folder, 'nm_single_mode')
        self.nm_single_fm_folder = path.join(consist_folder, host, type_na, 'nm_single_mode', 'fm')
        self.matrices_folder = path.join(self.pyc_na_folder, 'matrices')
        self.coef_md_folder = path.join(self.pyc_na_folder, 'coef_all_frames')
        self.md_proj_mode_c = path.join(self.coef_md_folder, 'coefficient.csv')

        self.avg_structure = self.get_avg_structure()
        self.eigenvectors = dict()
        self.eigenvectors_fm = dict()
        self.derivative_array = None
        self.mean_theta = None

    def get_avg_structure(self):
        d = miscell.read_structure(self.avg_crd)
        return d

    def get_atomids(self, i=0, j=8):
        u = Universe(self.avg_crd, self.avg_crd)
        basepairs = self.get_basepairs(u)
        h_ids = get_h_id_by_i_j(i, j)
        atomids = list()
        for h_id in h_ids:
            atomids.append(basepairs[h_id][0].indices[0] + 1)
            atomids.append(basepairs[h_id][1].indices[0] + 1)
        return atomids

    def get_md_proj_c_df(self):
        df = pd.read_csv(self.md_proj_mode_c)
        return df

    def read_eigenvector(self, modeid, basis='QHA'):
        if basis == 'FM':
            f_in = path.join(self.nm_single_fm_folder, 'mode.{0}.txt'.format(modeid))
        else:
            f_in = path.join(self.nm_single_folder, 'mode.{0}.txt'.format(modeid))
        data = np.genfromtxt(f_in)
        temp = list()
        for arr in data:
            for element in arr:
                temp.append(element)
        if basis == 'FM':
            self.eigenvectors_fm[modeid] = np.array(temp)
        else:
            self.eigenvectors[modeid] = np.array(temp)

    def get_eigenvector(self, modeid, basis='QHA'):
        try:
            if basis == 'FM':
                return self.eigenvectors_fm[modeid]
            else:
                return self.eigenvectors[modeid]
        except KeyError:
            self.read_eigenvector(modeid, basis=basis)
            if basis == 'FM':
                return self.eigenvectors_fm[modeid]
            else:
                return self.eigenvectors[modeid]

    def get_eigenvector_fm(self, modeid):
        try:
            return self.eigenvectors[modeid]
        except KeyError:
            self.read_eigenvector(modeid)
            return self.eigenvectors[modeid]

    def get_eigenvector_vector(self, mode, atomids, basis='QHA'):
        """

        :param mode:
        :param atomids:
        :param basis: QHA or FM
        :return:
        """
        eigenvector = self.get_eigenvector(mode, basis=basis)
        eigenvector = miscell.get_d_by_vector(eigenvector)
        vector = list()
        for atomid in atomids:
            for label in ['x', 'y', 'z']:
                vector.append(eigenvector[atomid][label])
        return vector

    def get_derivative_matrix(self, i_l=0, j_l=8):
        mat = list()
        atomids = self.get_atomids(i_l, j_l)
        positions = get_positions(atomids, self.avg_structure)
        midpoints, f, h, h_square, g, g_square = get_f_h_g_and_square(positions)

        #  mean_theta = self.get_mean_theta_from_avgcrd(i=i_l, j=j_l)
        #  cos_mean_theta = np.cos(mean_theta)
        cos_theta = f / (g * h)
        cos_theta_square = np.square(cos_theta)
        big_factor = -1 / np.sqrt(1-cos_theta_square)
        for i, atomid in enumerate(atomids):
            if i < 4:
                a = 1 / h
                b = g
                x_factor, y_factor, z_factor = get_xyz_factor(midpoints[3], midpoints[2])
                x_factor_1, y_factor_1, z_factor_1 = get_xyz_factor(midpoints[1], midpoints[0])
                if i < 2:
                    f_sign = -1/2
                    d_sign = -1
                else:
                    f_sign = 1/2
                    d_sign = 1
            else:
                a = 1 / g
                b = h
                x_factor, y_factor, z_factor = get_xyz_factor(midpoints[1], midpoints[0])
                x_factor_1, y_factor_1, z_factor_1 = get_xyz_factor(midpoints[3], midpoints[2])
                if i < 6:
                    f_sign = -1/2
                    d_sign = -1
                else:
                    f_sign = 1/2
                    d_sign = 1
            b_square = np.square(b)
            for j in range(3):
                if j == 0:
                    c = f_sign * x_factor
                    d = (d_sign * x_factor_1) / (2 * b)
                elif j == 1:
                    c = f_sign * y_factor
                    d = (d_sign * y_factor_1) / (2 * b)
                else:
                    c = f_sign * z_factor
                    d = (d_sign * z_factor_1) / (2 * b)
                answer = a * ((c * b) - (f * d)) / b_square
                mat.append(answer)
        mat = np.array(mat)
        mat = mat * big_factor
        #  mat = mat / cos_mean_theta
        return mat

    def get_derivative_matrix_benddistance(self, i_l=0, j_l=8):
        mat = list()
        atomids = self.get_atomids(i_l, j_l)
        positions = get_positions(atomids, self.avg_structure)
        big_x, big_y, big_z, r_i_lprime_norm = get_bigxyz_norm_benddistance(positions)
        for i, atomid in enumerate(atomids):
            if i in [2, 3, 6, 7]:
                factor = 1/2
            else:
                factor = -1/2
            for j in range(3):
                if j == 0:
                    a = big_x
                elif j == 1:
                    a = big_y
                else:
                    a = big_z
                answer = (a / r_i_lprime_norm) * factor
                mat.append(answer)
        mat = np.array(mat)
        return mat

    def get_derivative_matrix_costheta_square(self, i_l=0, j_l=8):
        mat = list()
        atomids = self.get_atomids(i_l, j_l)
        positions = get_positions(atomids, self.avg_structure)
        midpoints, f, h, h_square, g, g_square = get_f_h_g_and_square(positions)
        cos_theta = f / (g * h)
        big_factor = 2 * cos_theta
        for i, atomid in enumerate(atomids):
            if i < 4:
                a = 1 / h
                b = g
                x_factor, y_factor, z_factor = get_xyz_factor(midpoints[3], midpoints[2])
                x_factor_1, y_factor_1, z_factor_1 = get_xyz_factor(midpoints[1], midpoints[0])
                if i < 2:
                    f_sign = -1/2
                    d_sign = -1
                else:
                    f_sign = 1/2
                    d_sign = 1
            else:
                a = 1 / g
                b = h
                x_factor, y_factor, z_factor = get_xyz_factor(midpoints[1], midpoints[0])
                x_factor_1, y_factor_1, z_factor_1 = get_xyz_factor(midpoints[3], midpoints[2])
                if i < 6:
                    f_sign = -1/2
                    d_sign = -1
                else:
                    f_sign = 1/2
                    d_sign = 1
            b_square = np.square(b)
            for j in range(3):
                if j == 0:
                    c = f_sign * x_factor
                    d = (d_sign * x_factor_1) / (2 * b)
                elif j == 1:
                    c = f_sign * y_factor
                    d = (d_sign * y_factor_1) / (2 * b)
                else:
                    c = f_sign * z_factor
                    d = (d_sign * z_factor_1) / (2 * b)
                answer = a * ((c * b) - (f * d)) / b_square
                mat.append(answer)
        mat = np.array(mat)
        mat = big_factor * mat
        return mat

    def make_derivative_df(self, i=0, j=8, mode_range=(1, 1272), basis='QHA'):
        """

        :param i: l_i
        :param j: l_j
        :param mode_range:
        :param basis: QHA or FM
        :return:
        """
        atomids = self.get_atomids(i=i, j=j)
        d_mat = self.get_derivative_matrix(i_l=i, j_l=j)
        d_mat = d_mat.flatten()
        d_result = {'mode': list(), 'derivative': list()}
        for mode in range(mode_range[0], mode_range[1]+1):
            d_result['mode'].append(mode)
            eigvects = self.get_eigenvector_vector(mode, atomids, basis=basis)
            d_result['derivative'].append(np.dot(d_mat, eigvects))
        if basis == 'FM':
            f_out = path.join(self.derivative_folder,
                              'cosdelta.fm.mode{0}.to.mode{1}.i{2}j{3}.csv'.format(mode_range[0], mode_range[1], i, j))
        else:
            f_out = path.join(self.derivative_folder,
                              'cosdelta.mode{0}.to.mode{1}.i{2}j{3}.csv'.format(mode_range[0], mode_range[1], i, j))
        df = pd.DataFrame(d_result)
        df.to_csv(f_out, index=False)
        return df

    def read_derivative_df(self, i=0, j=8, mode_range=(1, 1272), basis='QHA'):
        if basis == 'FM':
            f_in = path.join(self.derivative_folder,
                             'cosdelta.fm.mode{0}.to.mode{1}.i{2}j{3}.csv'.format(mode_range[0], mode_range[1], i, j))
        else:
            f_in = path.join(self.derivative_folder,
                             'cosdelta.mode{0}.to.mode{1}.i{2}j{3}.csv'.format(mode_range[0], mode_range[1], i, j))
        df = pd.read_csv(f_in)
        return df

    def make_derivative_df_theta_ij_from_lstsq(self, beta_vec, v1_id='01', v2_id='89'):
        """

        :param beta_vec: dtheta/dx
        :param v1_id: l_j
        :param v2_id:
        :return:
        """
        i = int(v1_id[0])
        j = int(v2_id[0])
        atomids = self.get_atomids(i=i, j=j)
        d_mat = beta_vec
        d_result = {'mode': list(), 'derivative': list()}
        for mode in range(1, d_mode[self.type_na]+1):
            d_result['mode'].append(mode)
            eigvects = self.get_eigenvector_vector(mode, atomids)
            d_result['derivative'].append(np.dot(d_mat, eigvects))
        first = 'theta.v1{0}.v2{1}.lstsq.deriv.df.csv'.format(v1_id, v2_id)
        f_out = path.join(self.derivative_folder, first)
        df = pd.DataFrame(d_result)
        df.to_csv(f_out, index=False)
        return df

    def read_derivative_df_theta_ij_from_lstsq(self, v1_id='01', v2_id='89'):
        first = 'theta.v1{0}.v2{1}.lstsq.deriv.df.csv'.format(v1_id, v2_id)
        f_in = path.join(self.derivative_folder, first)
        df = pd.read_csv(f_in)
        return df

    def make_derivative_df_benddistance(self, i=0, j=8, mode_range=(1, 1272), basis='QHA'):
        """

        :param i: l_i
        :param j: l_j
        :param mode_range:
        :param basis: QHA or FM
        :return:
        """
        atomids = self.get_atomids(i=i, j=j)
        d_mat = self.get_derivative_matrix_benddistance(i_l=i, j_l=j)
        d_mat = d_mat.flatten()
        d_result = {'mode': list(), 'derivative': list()}
        for mode in range(mode_range[0], mode_range[1]+1):
            d_result['mode'].append(mode)
            eigvects = self.get_eigenvector_vector(mode, atomids, basis=basis)
            d_result['derivative'].append(np.dot(d_mat, eigvects))
        if basis == 'FM':
            f_out = path.join(self.derivative_folder,
                              'benddistance.fm.mode{0}.to.mode{1}.i{2}j{3}.csv'.format(mode_range[0],
                                                                                       mode_range[1], i, j))
        else:
            f_out = path.join(self.derivative_folder,
                              'benddistance.mode{0}.to.mode{1}.i{2}j{3}.csv'.format(mode_range[0], mode_range[1], i, j))
        df = pd.DataFrame(d_result)
        df.to_csv(f_out, index=False)
        return df

    def read_derivative_df_benddistance(self, i=0, j=8, mode_range=(1, 1272), basis='QHA'):
        if basis == 'FM':
            f_in = path.join(self.derivative_folder,
                             'benddistance.fm.mode{0}.to.mode{1}.i{2}j{3}.csv'.format(mode_range[0],
                                                                                      mode_range[1], i, j))
        else:
            f_in = path.join(self.derivative_folder,
                             'benddistance.mode{0}.to.mode{1}.i{2}j{3}.csv'.format(mode_range[0], mode_range[1], i, j))
        df = pd.read_csv(f_in)
        return df

    def make_derivative_df_costheta_square(self, i=0, j=8, mode_range=(1, 1272), basis='QHA'):
        """

        :param i: l_i
        :param j: l_j
        :param mode_range:
        :param basis: QHA or FM
        :return:
        """
        atomids = self.get_atomids(i=i, j=j)
        d_mat = self.get_derivative_matrix_costheta_square(i_l=i, j_l=j)
        d_mat = d_mat.flatten()
        d_result = {'mode': list(), 'derivative': list()}
        for mode in range(mode_range[0], mode_range[1]+1):
            d_result['mode'].append(mode)
            eigvects = self.get_eigenvector_vector(mode, atomids, basis=basis)
            d_result['derivative'].append(np.dot(d_mat, eigvects))
        if basis == 'FM':
            f_out = path.join(self.derivative_folder,
                              'costheta_square.fm.mode{0}.to.mode{1}.i{2}j{3}.csv'.format(mode_range[0],
                                                                                          mode_range[1], i, j))
        else:
            f_out = path.join(self.derivative_folder,
                              'costheta_square.mode{0}.to.mode{1}.i{2}j{3}.csv'.format(mode_range[0],
                                                                                       mode_range[1], i, j))
        df = pd.DataFrame(d_result)
        df.to_csv(f_out, index=False)
        return df

    def read_derivative_df_costheta_square(self, i=0, j=8, mode_range=(1, 1272), basis='QHA'):
        if basis == 'FM':
            f_in = path.join(self.derivative_folder,
                             'costheta_square.fm.mode{0}.to.mode{1}.i{2}j{3}.csv'.format(mode_range[0],
                                                                                         mode_range[1], i, j))
        else:
            f_in = path.join(self.derivative_folder,
                             'costheta_square.mode{0}.to.mode{1}.i{2}j{3}.csv'.format(mode_range[0], mode_range[1],
                                                                                      i, j))
        df = pd.read_csv(f_in)
        return df

    def make_dl_over_dc_2_df(self, i_l=0, j_l=8, mode_range=(1, 1272), basis='QHA'):
        df = self.read_derivative_df(i=i_l, j=j_l, mode_range=mode_range, basis=basis)
        vector = df['derivative'].tolist()
        mat = list()
        for i in range(mode_range[0], mode_range[1]+1):
            temp = list()
            for j in range(mode_range[0], mode_range[1]+1):
                value = vector[i-1] * vector[j-1]
                temp.append(value)
            mat.append(temp)
        mat = np.array(mat)
        if basis == 'FM':
            f_out = path.join(self.matrices_folder,
                              'dcosdelta_over_dc_2.fm.mode{0}.to.mode{1}.i{2}j{3}.npy'.format(mode_range[0],
                                                                                              mode_range[1], i_l, j_l))
        else:
            f_out = path.join(self.matrices_folder,
                              'dcosdelta_over_dc_2.mode{0}.to.mode{1}.i{2}j{3}.npy'.format(mode_range[0],
                                                                                           mode_range[1], i_l, j_l))
        np.save(f_out, mat)
        return mat

    def read_dl_over_dc_2_df(self, i=0, j=8, mode_range=(1, 1272), basis='QHA'):
        if basis == 'FM':
            f_in = path.join(self.matrices_folder,
                             'dcosdelta_over_dc_2.fm.mode{0}.to.mode{1}.i{2}j{3}.npy'.format(mode_range[0],
                                                                                             mode_range[1], i, j))
        else:
            f_in = path.join(self.matrices_folder,
                             'dcosdelta_over_dc_2.mode{0}.to.mode{1}.i{2}j{3}.npy'.format(mode_range[0],
                                                                                          mode_range[1], i, j))
        mat = np.load(f_in)
        return mat

    def make_dtheta_over_dc_2_df_lstsq(self, v1_id='01', v2_id='89'):
        df = self.read_derivative_df_theta_ij_from_lstsq(v1_id=v1_id, v2_id=v2_id)
        vector = df['derivative'].tolist()
        mat = list()
        for i in range(1, d_mode[self.type_na]+1):
            temp = list()
            for j in range(1, d_mode[self.type_na]+1):
                value = vector[i-1] * vector[j-1]
                temp.append(value)
            mat.append(temp)
        mat = np.array(mat)
        first = 'dtheta_dc2_v1{0}v2{1}.npy'.format(v1_id, v2_id)
        f_out = path.join(self.matrices_folder, first)
        np.save(f_out, mat)
        return mat

    def read_dtheta_over_dc_2_df_lstsq(self, v1_id='01', v2_id='89'):
        first = 'dtheta_dc2_v1{0}v2{1}.npy'.format(v1_id, v2_id)
        f_in = path.join(self.matrices_folder, first)
        mat = np.load(f_in)
        return mat

    def make_dtheta_over_dc_2_df_lstsq_version2(self, beta_vec, v1_id='01', v2_id='89'):
        vector = beta_vec
        mat = list()
        for i in range(1, d_mode[self.type_na]+1):
            temp = list()
            for j in range(1, d_mode[self.type_na]+1):
                value = vector[i-1] * vector[j-1]
                temp.append(value)
            mat.append(temp)
        mat = np.array(mat)
        first = 'dtheta_dc2_v1{0}v2{1}.version2.npy'.format(v1_id, v2_id)
        f_out = path.join(self.matrices_folder, first)
        np.save(f_out, mat)
        return mat

    def read_dtheta_over_dc_2_df_lstsq_version2(self, v1_id='01', v2_id='89'):
        first = 'dtheta_dc2_v1{0}v2{1}.version2.npy'.format(v1_id, v2_id)
        f_in = path.join(self.matrices_folder, first)
        mat = np.load(f_in)
        return mat

    def make_dl_over_dc_2_df_pseudo(self, v1_id='01', v2_id='89', mode_range=(1, 1272)):
        df = self.read_df_theta_numeri_c_eigenmode(v1_id=v1_id, v2_id=v2_id)
        vector = df['derivative'].tolist()
        mat = list()
        for i in range(mode_range[0], mode_range[1]+1):
            temp = list()
            for j in range(mode_range[0], mode_range[1]+1):
                value = vector[i-1] * vector[j-1]
                temp.append(value)
            mat.append(temp)
        mat = np.array(mat)
        f_out = path.join(self.matrices_folder,
                          'dtheta_over_dc_2.pseudo.mode{0}.to.mode{1}.npy'.format(mode_range[0], mode_range[1]))
        np.save(f_out, mat)
        return mat

    def read_dl_over_dc_2_df_pseudo(self, mode_range=(1, 1272)):
        f_in = path.join(self.matrices_folder,
                         'dtheta_over_dc_2.pseudo.mode{0}.to.mode{1}.npy'.format(mode_range[0], mode_range[1]))
        mat = np.load(f_in)
        return mat

    def make_dl_over_dc_2_df_benddistance(self, i_l=0, j_l=8, mode_range=(1, 1272), basis='QHA'):
        df = self.read_derivative_df_benddistance(i=i_l, j=j_l, mode_range=mode_range, basis=basis)
        vector = df['derivative'].tolist()
        mat = list()
        for i in range(mode_range[0], mode_range[1]+1):
            temp = list()
            for j in range(mode_range[0], mode_range[1]+1):
                value = vector[i-1] * vector[j-1]
                temp.append(value)
            mat.append(temp)
        mat = np.array(mat)
        if basis == 'FM':
            f_out = path.join(self.matrices_folder,
                              'dbenddist_over_dc_2.fm.mode{0}.to.mode{1}.i{2}j{3}.npy'.format(mode_range[0],
                                                                                              mode_range[1], i_l, j_l))
        else:
            f_out = path.join(self.matrices_folder,
                              'dbenddist_over_dc_2.mode{0}.to.mode{1}.i{2}j{3}.npy'.format(mode_range[0],
                                                                                           mode_range[1], i_l, j_l))
        np.save(f_out, mat)
        return mat

    def read_dl_over_dc_2_df_benddistance(self, i=0, j=8, mode_range=(1, 1272), basis='QHA'):
        if basis == 'FM':
            f_in = path.join(self.matrices_folder,
                             'dbenddist_over_dc_2.fm.mode{0}.to.mode{1}.i{2}j{3}.npy'.format(mode_range[0],
                                                                                             mode_range[1], i, j))
        else:
            f_in = path.join(self.matrices_folder,
                             'dbenddist_over_dc_2.mode{0}.to.mode{1}.i{2}j{3}.npy'.format(mode_range[0],
                                                                                          mode_range[1], i, j))
        mat = np.load(f_in)
        return mat

    def make_dl_over_dc_2_df_costheta_square(self, i_l=0, j_l=8, mode_range=(1, 1272), basis='QHA'):
        df = self.read_derivative_df_costheta_square(i=i_l, j=j_l, mode_range=mode_range, basis=basis)
        vector = df['derivative'].tolist()
        mat = list()
        for i in range(mode_range[0], mode_range[1]+1):
            temp = list()
            for j in range(mode_range[0], mode_range[1]+1):
                value = vector[i-1] * vector[j-1]
                temp.append(value)
            mat.append(temp)
        mat = np.array(mat)
        if basis == 'FM':
            f_out = path.join(self.matrices_folder,
                              'dcostheta_square_over_dc_2.fm.mode{0}.to.mode{1}.i{2}j{3}.npy'.format(mode_range[0],
                                                                                                     mode_range[1],
                                                                                                     i_l, j_l))
        else:
            f_out = path.join(self.matrices_folder,
                              'dcostheta_square_over_dc_2.mode{0}.to.mode{1}.i{2}j{3}.npy'.format(mode_range[0],
                                                                                                  mode_range[1],
                                                                                                  i_l, j_l))
        np.save(f_out, mat)
        return mat

    def read_dl_over_dc_2_df_costheta_square(self, i=0, j=8, mode_range=(1, 1272), basis='QHA'):
        if basis == 'FM':
            f_in = path.join(self.matrices_folder,
                             'dcostheta_square_over_dc_2.fm.mode{0}.to.mode{1}.i{2}j{3}.npy'.format(mode_range[0],
                                                                                                    mode_range[1],
                                                                                                    i, j))
        else:
            f_in = path.join(self.matrices_folder,
                             'dcostheta_square_over_dc_2.mode{0}.to.mode{1}.i{2}j{3}.npy'.format(mode_range[0],
                                                                                                 mode_range[1], i, j))
        mat = np.load(f_in)
        return mat
    
    def make_covar_mat_of_c(self, mode_range=(1, 1272), basis='QHA'):
        if basis == 'FM':
            f_in = path.join(self.coef_md_folder, 'coefficient.fm.csv')
        else:
            f_in = path.join(self.coef_md_folder, 'coefficient.csv')
        df = pd.read_csv(f_in)
        mat = list()
        for i in range(mode_range[0], mode_range[1]+1):
            temp = list()
            for j in range(mode_range[0], mode_range[1]+1):
                value = (df[str(i)] * df[str(j)]).mean()
                temp.append(value)
            mat.append(temp)
        mat = np.array(mat)
        if basis == 'FM':
            f_out = path.join(self.matrices_folder,
                              'covar_mat_of_c.fm.mode{0}.to.mode{1}.npy'.format(mode_range[0], mode_range[1]))
        else:
            f_out = path.join(self.matrices_folder,
                              'covar_mat_of_c.mode{0}.to.mode{1}.npy'.format(mode_range[0], mode_range[1]))
        np.save(f_out, mat)
        return mat    

    def read_covar_mat_of_c(self, mode_range=(1, 1272), basis='QHA'):
        if basis == 'FM':
            f_in = path.join(self.matrices_folder,
                             'covar_mat_of_c.fm.mode{0}.to.mode{1}.npy'.format(mode_range[0], mode_range[1]))
        else:
            f_in = path.join(self.matrices_folder,
                             'covar_mat_of_c.mode{0}.to.mode{1}.npy'.format(mode_range[0], mode_range[1]))
        mat = np.load(f_in)
        return mat

    def write_theta_normalized_variance_matrix(self, normal_var_mat, i=0, j=8, basis='QHA'):
        if basis == 'FM':
            f_out = path.join(self.matrices_folder,
                              'normalized_variance_matrix_theta.fm.i{0}j{1}.npy'.format(i, j))
        else:
            f_out = path.join(self.matrices_folder,
                              'normalized_variance_matrix_theta.i{0}j{1}.npy'.format(i, j))
        np.save(f_out, normal_var_mat)

    def read_theta_normalized_variance_matrix(self, i=0, j=8, basis='QHA'):
        if basis == 'FM':
            f_in = path.join(self.matrices_folder,
                             'normalized_variance_matrix_theta.fm.i{0}j{1}.npy'.format(i, j))
        else:
            f_in = path.join(self.matrices_folder,
                             'normalized_variance_matrix_theta.i{0}j{1}.npy'.format(i, j))
        mat = np.load(f_in)
        return mat

    def write_benddist_normalized_variance_matrix(self, normal_var_mat, i=0, j=8, basis='QHA'):
        if basis == 'FM':
            f_out = path.join(self.matrices_folder,
                              'normalized_variance_matrix_benddist.fm.i{0}j{1}.npy'.format(i, j))
        else:
            f_out = path.join(self.matrices_folder,
                              'normalized_variance_matrix_benddist.i{0}j{1}.npy'.format(i, j))
        np.save(f_out, normal_var_mat)

    def read_benddist_normalized_variance_matrix(self, i=0, j=8, basis='QHA'):
        if basis == 'FM':
            f_in = path.join(self.matrices_folder,
                             'normalized_variance_matrix_benddist.fm.i{0}j{1}.npy'.format(i, j))
        else:
            f_in = path.join(self.matrices_folder,
                             'normalized_variance_matrix_benddist.i{0}j{1}.npy'.format(i, j))
        mat = np.load(f_in)
        return mat

    def write_costheta_square_normalized_variance_matrix(self, normal_var_mat, i=0, j=8, basis='QHA'):
        if basis == 'FM':
            f_out = path.join(self.matrices_folder,
                              'normalized_variance_matrix_costheta_square.fm.i{0}j{1}.npy'.format(i, j))
        else:
            f_out = path.join(self.matrices_folder,
                              'normalized_variance_matrix_costheta_square.i{0}j{1}.npy'.format(i, j))
        np.save(f_out, normal_var_mat)

    def read_costheta_square_normalized_variance_matrix(self, i=0, j=8, basis='QHA'):
        if basis == 'FM':
            f_in = path.join(self.matrices_folder,
                             'normalized_variance_matrix_costheta_square.fm.i{0}j{1}.npy'.format(i, j))
        else:
            f_in = path.join(self.matrices_folder,
                             'normalized_variance_matrix_costheta_square.i{0}j{1}.npy'.format(i, j))
        mat = np.load(f_in)
        return mat

    def read_coefficient_df(self, basis='QHA'):
        if basis == 'FM':
            f_in = path.join(self.coef_md_folder, 'coefficient.fm.csv')
        else:
            f_in = path.join(self.coef_md_folder, 'coefficient.csv')
        df = pd.read_csv(f_in)
        return df

    def get_benddistance_from_avgcrd(self, i=0, j=8, n_bead=10):
        u = Universe(self.avg_crd, self.avg_crd)
        basepairs = self.get_basepairs(u)
        midpoints = self.get_midpoints(basepairs, n_bead=n_bead)
        h_ids = get_h_id_by_i_j(i, j)
        r_i_lprime = midpoints[h_ids[3]] - midpoints[h_ids[2]] + midpoints[h_ids[1]] - midpoints[h_ids[0]]
        return np.linalg.norm(r_i_lprime)

    def get_normalized_benddistance_from_avgcrd(self, i=0, j=8, n_bead=10):
        u = Universe(self.avg_crd, self.avg_crd)
        basepairs = self.get_basepairs(u)
        midpoints = self.get_midpoints(basepairs, n_bead=n_bead)
        h_ids = get_h_id_by_i_j(i, j)
        l1 = midpoints[h_ids[1]] - midpoints[h_ids[0]]
        l1_norm = np.linalg.norm(l1)
        l2 = midpoints[h_ids[3]] - midpoints[h_ids[2]]
        l2_norm = np.linalg.norm(l2)
        r_i_lprime = l1 / l1_norm + l2 / l2_norm
        return np.linalg.norm(r_i_lprime)

    def get_costheta_square_from_avgcrd(self, i=0, j=8, n_bead=10):
        u = Universe(self.avg_crd, self.avg_crd)
        basepairs = self.get_basepairs(u)
        midpoints = self.get_midpoints(basepairs, n_bead=n_bead)
        h_ids = get_h_id_by_i_j(i, j)
        l1 = midpoints[h_ids[1]] - midpoints[h_ids[0]]
        l1_norm = np.linalg.norm(l1)
        l2 = midpoints[h_ids[3]] - midpoints[h_ids[2]]
        l2_norm = np.linalg.norm(l2)
        return np.square(np.dot(l1, l2) / (l1_norm * l2_norm))
        #  return np.dot(l1, l2) / (l1_norm * l2_norm)

    def get_cosdeltatheta_from_avgcrd(self, i=0, j=8, n_bead=10):
        u = Universe(self.avg_crd, self.avg_crd)
        basepairs = self.get_basepairs(u)
        midpoints = self.get_midpoints(basepairs, n_bead=n_bead)
        h_ids = get_h_id_by_i_j(i, j)
        l1 = midpoints[h_ids[1]] - midpoints[h_ids[0]]
        l1_norm = np.linalg.norm(l1)
        l2 = midpoints[h_ids[3]] - midpoints[h_ids[2]]
        l2_norm = np.linalg.norm(l2)
        mean_theta = self.get_mean_theta_from_avgcrd(i, j, n_bead)
        temp = np.square(np.dot(l1, l2) / (l1_norm * l2_norm))
        return temp / np.cos(mean_theta)
        #  return np.dot(l1, l2) / (l1_norm * l2_norm)

    def get_benddistance_from_allatommd(self, i=0, j=8, n_bead=10):
        u = Universe(self.avg_crd, self.fit_dcd)
        basepairs = self.get_basepairs(u)
        h_ids = get_h_id_by_i_j(i, j)
        d = {'Frame': list(), 'BendDistance': list()}
        for ts in u.trajectory:
            midpoints = self.get_midpoints(basepairs, n_bead=n_bead)
            r_i_lprime = midpoints[h_ids[3]] - midpoints[h_ids[2]] + midpoints[h_ids[1]] - midpoints[h_ids[0]]
            d['Frame'].append(ts.frame + 1)
            d['BendDistance'].append(np.linalg.norm(r_i_lprime))
        df = pd.DataFrame(d)
        f_out = path.join(self.cosdelta_folder, 'md.benddistance.i{0}j{1}.csv'.format(i, j))
        df.to_csv(f_out, index=False)
        return df

    def get_normalized_benddistance_from_allatommd(self, i=0, j=8, n_bead=10):
        u = Universe(self.avg_crd, self.fit_dcd)
        basepairs = self.get_basepairs(u)
        h_ids = get_h_id_by_i_j(i, j)
        d = {'Frame': list(), 'BendDistance': list()}
        for ts in u.trajectory:
            midpoints = self.get_midpoints(basepairs, n_bead=n_bead)
            l1 = midpoints[h_ids[1]] - midpoints[h_ids[0]]
            l1_norm = np.linalg.norm(l1)
            l2 = midpoints[h_ids[3]] - midpoints[h_ids[2]]
            l2_norm = np.linalg.norm(l2)
            r_i_lprime = l1 / l1_norm + l2 / l2_norm
            d['Frame'].append(ts.frame + 1)
            d['BendDistance'].append(np.linalg.norm(r_i_lprime))
        df = pd.DataFrame(d)
        f_out = path.join(self.cosdelta_folder, 'md.normalized.benddistance.i{0}j{1}.csv'.format(i, j))
        df.to_csv(f_out, index=False)
        return df

    def get_costheta_square_from_allatommd(self, i=0, j=8, n_bead=10):
        u = Universe(self.avg_crd, self.fit_dcd)
        basepairs = self.get_basepairs(u)
        h_ids = get_h_id_by_i_j(i, j)
        d = {'Frame': list(), 'BendDistance': list()}
        for ts in u.trajectory:
            midpoints = self.get_midpoints(basepairs, n_bead=n_bead)
            l1 = midpoints[h_ids[1]] - midpoints[h_ids[0]]
            l1_norm = np.linalg.norm(l1)
            l2 = midpoints[h_ids[3]] - midpoints[h_ids[2]]
            l2_norm = np.linalg.norm(l2)
            d['Frame'].append(ts.frame + 1)
            d['BendDistance'].append(np.square(np.dot(l1, l2) / (l1_norm * l2_norm)))
            #  d['BendDistance'].append(np.dot(l1, l2) / (l1_norm * l2_norm))
        df = pd.DataFrame(d)
        f_out = path.join(self.cosdelta_folder, 'md.costheta_square.i{0}j{1}.csv'.format(i, j))
        df.to_csv(f_out, index=False)
        return df

    def read_benddistance_from_allatommd(self, i=0, j=8):
        f_in = path.join(self.cosdelta_folder, 'md.benddistance.i{0}j{1}.csv'.format(i, j))
        df = pd.read_csv(f_in)
        return df

    def read_normalized_benddistance_from_allatommd(self, i=0, j=8):
        f_in = path.join(self.cosdelta_folder, 'md.normalized.benddistance.i{0}j{1}.csv'.format(i, j))
        df = pd.read_csv(f_in)
        return df

    def read_costheta_square_from_allatommd(self, i=0, j=8):
        f_in = path.join(self.cosdelta_folder, 'md.costheta_square.i{0}j{1}.csv'.format(i, j))
        df = pd.read_csv(f_in)
        return df

    def read_x_matrix(self):
        f_in = path.join(self.matrices_folder, 'x_matrix_lstsq.npy')
        mat = np.load(f_in)
        return mat

    def write_theta_ij_beta_vec(self, beta_vec, i=0, j=8):
        f_out = path.join(self.matrices_folder, 'betavec_theta_i{0}j{1}.npy'.format(i, j))
        np.save(f_out, beta_vec)
        print('Write beta vector into {0}'.format(f_out))

    def read_theta_ij_beta_vec(self, i=0, j=8):
        f_in = path.join(self.matrices_folder, 'betavec_theta_i{0}j{1}.npy'.format(i, j))
        mat = np.load(f_in)
        return mat

    def get_theta_numerical_eigenmode(self, modeid, c_agent, avg_coords_flat, v1_id='01', v2_id='89', alpha=None):
        tilde_r_array = c_agent.get_tilde_r(modeid, alpha=alpha)
        f_out = c_agent.make_tilde_r_newcrd_change_segid(tilde_r_array, modeid, d_segid_abbr)
        three_points = self.get_align_threepoints_arbitrary_crd(f_out, v1_id=v1_id, v2_id=v2_id)
        theta = get_theta_btw_v1v2(three_points)
        tilde_r_array_flat = tilde_r_array.flatten()
        diff = tilde_r_array_flat - avg_coords_flat
        eigvector = c_agent.get_eigenvector_qha(modeid)
        c = np.dot(diff, eigvector)
        return theta, c

    def get_theta_analytical_eigenmode(self, modeid, c, v1_id='01', v2_id='89'):
        i = int(v1_id[0])
        j = int(v2_id[0])
        if self.derivative_array is None:
            df_deriv = self.read_derivative_df(i=i, j=j, mode_range=(1, d_mode[self.type_na]), basis='QHA')
            self.derivative_array = df_deriv['derivative'].values
        deriv_of_mode = self.derivative_array[modeid-1]
        if self.mean_theta is None:
            self.mean_theta = self.get_mean_theta_from_avgcrd(i, j)
        theta = self.mean_theta + c * deriv_of_mode
        return theta

    def make_df_theta_analy_numeri_eigenmode(self, first_mode, last_mode, c_agent, avg_coords_flat,
                                             v1_id='01', v2_id='89', alpha=None):
        d_result = {'Mode-ID': list(), 'theta_analytical': list(), 'theta_numerical': list()}
        for modeid in range(first_mode, last_mode+1):
            theta_numerical, c = self.get_theta_numerical_eigenmode(modeid, c_agent, avg_coords_flat,
                                                                    v1_id=v1_id, v2_id=v2_id, alpha=alpha)
            theta_analytical = self.get_theta_analytical_eigenmode(modeid, c, v1_id=v1_id, v2_id=v2_id)
            d_result['Mode-ID'].append(modeid)
            d_result['theta_analytical'].append(theta_analytical)
            d_result['theta_numerical'].append(theta_numerical)
        df = pd.DataFrame(d_result)
        df = df[['Mode-ID', 'theta_analytical', 'theta_numerical']]
        if alpha is None:
            f_out = path.join(self.cosdelta_folder, 'theta_numer_analy_mode{0}_mode{1}.csv'.format(first_mode,
                                                                                                   last_mode))
        else:
            first = 'theta_numer_analy_mode{0}_mode{1}_'.format(first_mode, last_mode)
            second = 'alpha{0:.2f}.csv'.format(alpha)
            third = first + second
            f_out = path.join(self.cosdelta_folder, third)
        df.to_csv(f_out, index=False)
        return df

    def read_df_theta_analy_numeri_eigenmode(self, first_mode, last_mode, alpha=None):
        if alpha is None:
            f_in = path.join(self.cosdelta_folder, 'theta_numer_analy_mode{0}_mode{1}.csv'.format(first_mode,
                                                                                                  last_mode))
        else:
            first = 'theta_numer_analy_mode{0}_mode{1}_'.format(first_mode, last_mode)
            second = 'alpha{0:.2f}.csv'.format(alpha)
            third = first + second
            f_in = path.join(self.cosdelta_folder, third)
        df = pd.read_csv(f_in)
        return df

    def make_df_theta_numeri_c_eigenmode(self, c_agent, avg_coords_flat, v1_id='01', v2_id='89'):
        i = int(v1_id[0])
        j = int(v2_id[0])
        d_result = {'Mode-ID': list(), 'theta_numerical': list(), 'c': list(), 'derivative': list()}
        mean_theta = self.get_mean_theta_from_avgcrd(i=i, j=j)
        for modeid in range(1, d_mode[self.type_na]+1):
            theta_numerical, c = self.get_theta_numerical_eigenmode(modeid, c_agent, avg_coords_flat,
                                                                    v1_id=v1_id, v2_id=v2_id)
            derivative = (theta_numerical - mean_theta) / c
            d_result['Mode-ID'].append(modeid)
            d_result['theta_numerical'].append(theta_numerical)
            d_result['c'].append(c)
            d_result['derivative'].append(derivative)
        df = pd.DataFrame(d_result)
        df = df[['Mode-ID', 'theta_numerical', 'c', 'derivative']]
        f_out = path.join(self.cosdelta_folder, 'theta_i{0}j{1}_pseudo_derivative.csv'.format(i, j))
        df.to_csv(f_out, index=False)
        return df

    def read_df_theta_numeri_c_eigenmode(self, v1_id='01', v2_id='89'):
        i = int(v1_id[0])
        j = int(v2_id[0])
        f_in = path.join(self.cosdelta_folder, 'theta_i{0}j{1}_pseudo_derivative.csv'.format(i, j))
        df = pd.read_csv(f_in)
        return df

    def get_r_small_perturbation(self, posi_id, delta_x, xyz, v1_id='01', v2_id='89'):
        """

        :param posi_id: 0-23
        :param delta_x:
        :param v1_id:
        :param v2_id:
        :param xyz: 'x', 'y', 'z'
        :return:
        """
        crd = crd_util.CRDAgent(self.avg_crd)
        tilde_r = np.array(crd.get_xyz_batch())

        i = int(v1_id[0])
        j = int(v2_id[0])
        atomids = self.get_atomids(i=i, j=j)
        atomid = atomids[posi_id]

        displacement = np.zeros(tilde_r.flatten().shape)
        displacement = miscell.get_d_by_vector(displacement)
        displacement[atomid][xyz] = delta_x
        displacement = miscell.get_vector_by_d(displacement)
        displacement = displacement.reshape((self.n_atoms, 3))
        tilde_r = tilde_r + displacement
        return tilde_r

    def make_small_r_newcrd(self, coords, posi_id):
        old_crd = crd_util.CRDAgent(self.avg_crd)
        old_crd.set_xyz_batch(coords)
        old_crd.change_segid(d_segid_abbr)
        f_out = path.join(self.smallr_folder, 'small_r.mode{0}.crd'.format(posi_id))
        old_crd.write_coor(f_out)
        return f_out


def get_distance_sum_by_bpids(midpoints, bpid1=1, bpid2=2):
    h_array = [np.linalg.norm(midpoints[bp_id] - midpoints[bp_id + 1])
               for bp_id in range(bpid1, bpid2)]
    return sum(h_array)


def get_lee_by_bpids(points, bpid1, bpid2):
    return np.linalg.norm(points[bpid1] - points[bpid2])


def align_ee_to_zaxis(points):
    ee_vector = points[-1] - points[0]
    rot_mat = miscell.rotation_mat_between_2vect(ee_vector, z_axis)
    points_after_transform = np.dot(rot_mat, points.T)
    trans_vector = -points_after_transform[:, 0]
    for col_idx in range(points_after_transform.shape[1]):
        points_after_transform[:, col_idx] = points_after_transform[:, col_idx] + trans_vector
    return points_after_transform.T


def align_cross_of_v1v2_to_zaxis(threepoints):
    v1 = threepoints[1] - threepoints[0]
    v2 = threepoints[2] - threepoints[0]
    cross_vec = np.cross(v1, v2)
    rot_mat = miscell.rotation_mat_between_2vect(cross_vec, z_axis)
    points_after_transform = np.dot(rot_mat, threepoints.T)
    return points_after_transform.T


def align_v1_to_xaxis(threepoints):
    v1 = threepoints[1] - threepoints[0]
    rot_mat = miscell.rotation_mat_between_2vect(v1, x_axis)
    points_after_transform = np.dot(rot_mat, threepoints.T)
    return points_after_transform.T


def get_theta_btw_v1v2(threepoints):
    v1 = threepoints[1] - threepoints[0]
    v2 = threepoints[2] - threepoints[0]
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)
    return np.arccos(np.dot(v1, v2))


def get_tangent_theta_phi(v, l):
    v_unit = v / l  # Make a unit vector

    # get v_unit projected onto yz, v_unit_yz
    v_unit_yz = v_unit - (np.dot(v_unit, x_axis) * x_axis)
    print('v_yz: {0}'.format(v_unit_yz))
    # theta = miscell.angle_between(z_axis, v_unit_yz)
    theta = np.arccos(np.dot(z_axis, v_unit_yz))

    # get v_unit projected onto xy, v_unit_xy
    v_unit_xy = v_unit - (np.dot(v_unit, z_axis) * z_axis)
    print('v_xy: {0}'.format(v_unit_xy))
    phi = miscell.angle_between(y_axis, v_unit_xy)
    # phi = np.arccos(np.dot(y_axis, v_unit_xy))

    return v_unit[0], v_unit[1], v_unit[2], theta, phi


def get_pair_list(n_bead):
    return list(combinations(range(n_bead - 1), 2))


def get_pair_list_by_beadid(i, j):
    """
    Basically, id start from 0
    For example: 10 basepairs
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9
    :param i: id of bead_i
    :param j: id of bead_j
    :return:
    """
    n_vector = j - i
    if n_vector < 2:
        raise CustomError("Number of vector must greater or equal than 2")
    temp_list = range(i, j)
    return list(combinations(temp_list, 2)), temp_list


def get_ij_pairs(n_bead):
    max_id = 9
    result_list = list()
    i = 0
    j = i + n_bead - 1
    while j <= max_id:
        result_list.append((i,j))
        i += 1
        j += 1
    return result_list


def polynomial(alpha, coefficients, c):
    result = c
    for i in range(len(coefficients)):
        key = i + 1
        result += coefficients[key] * (alpha ** key)
    return result


def get_lp(alpha, l_avg):
    lp = -l_avg / (2 * np.log(alpha))
    lp = lp / 10
    return lp


def get_h_id_by_i_j(i, j):
    return [i+1, i+2, j+1, j+2]


def get_positions(atomids, struct):
    positions = list()
    for atomid in atomids:
        positions.append(miscell.get_position_array(struct, atomid))
    return np.array(positions)


def get_midpoint_by_posi_posj(posi, posj):
    return (posi + posj) / 2


def get_f_h_g_and_square(positions):
    midpoints = list()
    for midpoint_id in range(4):
        i = midpoint_id * 2
        j = i+1
        midpoints.append(get_midpoint_by_posi_posj(positions[i], positions[j]))
    l1 = midpoints[1] - midpoints[0]
    l2 = midpoints[3] - midpoints[2]
    f = np.dot(l1, l2)
    g = np.linalg.norm(l1)
    h = np.linalg.norm(l2)
    h_square = np.square(h)
    g_square = np.square(g)
    return midpoints, f, h, h_square, g, g_square


def get_xyz_factor(midpoint1, midpoint2):
    diff = midpoint1 - midpoint2
    return diff[0], diff[1], diff[2]


def get_bigxyz_norm_benddistance(positions):
    midpoints = list()
    for midpoint_id in range(4):
        i = midpoint_id * 2
        j = i+1
        midpoints.append(get_midpoint_by_posi_posj(positions[i], positions[j]))
    r_i_lprime = midpoints[3] - midpoints[2] + midpoints[1] - midpoints[0]
    r_i_lprime_norm = np.linalg.norm(r_i_lprime)
    big_x = r_i_lprime[0]
    big_y = r_i_lprime[1]
    big_z = r_i_lprime[2]
    return big_x, big_y, big_z, r_i_lprime_norm


def get_persistence_length(i, j, df_data, mean_theta, host, type_na):
    mask = (df_data['i'] == i) & (df_data['j'] == j)
    df_temp = df_data[mask]
    #  mean_theta = df_temp['theta'].mean()  # Debug
    delta_theta_array = df_temp['theta'].values - mean_theta
    cos_delta_theta = np.cos(delta_theta_array)
    mean_cos_delta_theta = cos_delta_theta.mean()
    j_minus_i = j - i
    l_ensemble_average = d_l_avg['base'][host][type_na]
    log_mean_cos_delta_theta = np.log(mean_cos_delta_theta)
    lp = (-j_minus_i * l_ensemble_average) / (2 * log_mean_cos_delta_theta)
    return lp, j_minus_i


def get_position_ids_from_vect_id(v_id):
    return int(v_id[0]), int(v_id[1])


def translate_two_heads_to_origin(midpoints, v1_id='01', v2_id='89'):
    r1_id, r2_id = get_position_ids_from_vect_id(v1_id)
    r3_id, r4_id = get_position_ids_from_vect_id(v2_id)
    v1 = midpoints[r2_id] - midpoints[r1_id]
    v2 = midpoints[r4_id] - midpoints[r3_id]
    return np.array([[0, 0, 0], v1, v2])

