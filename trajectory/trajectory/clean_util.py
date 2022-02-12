from os import path
from shutil import copyfile
import re
import numpy as np
from MDAnalysis import Universe, Writer
import pandas as pd
import miscell

root_folder = '/Users/yizao/PycharmProjects/connect_macro_micro'
idea_cmm_folder = '/Users/yizao/IdeaProjects/python/connect_macro_micro'
sequences = {'pnas': {'arna+arna': {'guide': 'CAAUGGAGUA',
                                    'target': 'UACUCCAUUG'},
                      'bdna+bdna': {'guide': 'CAATGGAGTA',
                                    'target': 'TACTCCATTG'}},
             'pnas_amber': {'arna+arna': {'guide': 'CAAUGGAGUA',
                                          'target': 'UACUCCAUUG'},
                            'bdna+bdna': {'guide': 'CAATGGAGTA',
                                          'target': 'TACTCCATTG'}}
             }


class CleanTrajAgent:
    def __init__(self, host, type_na):
        self.host = host
        self.type_na = type_na
        self.host_folder = path.join(root_folder, host)
        self.na_folder = path.join(self.host_folder, type_na)
        self.idea_na_folder = path.join(idea_cmm_folder, host, type_na)
        self.clean_traj_folder = path.join(self.na_folder, 'clean_traj')
        self.hb_dat = path.join(self.idea_na_folder, 'x3dna', '0', 'h-bond_g.dat')
        self.hb_csv = path.join(self.clean_traj_folder, 'h-bond_g.csv')
        self.guide_sequence = self.get_guide_sequence()
        self.refhblist = self.get_hb_list()

        self.inp_gro = path.join(self.idea_na_folder, 'input', 'allatoms', '{0}.npt4.all.gro'.format(type_na))
        self.inp_xtc = path.join(self.idea_na_folder, 'input', 'allatoms', '{0}.all.xtc'.format(type_na))
        self.out_gro = path.join(self.clean_traj_folder, '{0}.npt4.all.gro'.format(type_na))
        self.out_xtc = path.join(self.clean_traj_folder, '{0}.all.xtc'.format(type_na))

        self.make_folders()

        if not path.exists(self.hb_csv):
            self.hb_df = make_df_by_hbond(self.hb_dat, self.guide_sequence)
            self.hb_df.to_csv(self.hb_csv, index=False)
        else:
            self.hb_df = pd.read_csv(self.hb_csv)

        self.clean_data_ratio = None

    def get_guide_sequence(self):
        seq = sequences[self.host][self.type_na]['guide']
        return seq

    def get_hb_list(self):
        d = {'A': '2', 'T': '2', 'U': '2', 'C': '3', 'G': '3'}
        return ''.join([d[b] for b in self.guide_sequence])

    def make_folders(self):
        folders = [self.clean_traj_folder]
        for folder in folders:
            miscell.check_dir_exist_and_make(folder)

    def fix_hb_df(self):
        df = self.hb_df.drop(self.hb_df.index[[1]])
        df.to_csv(self.hb_csv, index=False)
        return df

    def get_clean_df_strict_criteria(self):
        mask = (self.hb_df['hb_num'].astype(str) == self.refhblist)
        df = self.hb_df[mask]
        df = df.reset_index()
        self.clean_data_ratio = df.shape[0] * 100. / self.hb_df.shape[0]
        print('After cleaning, {0} % data left.'.format(round(self.clean_data_ratio,2)))
        return df

    def get_clean_df_weak_hb(self):
        self.hb_df['criteria'] = self.hb_df.apply(weak_hb, ref=self.refhblist, axis=1)
        mask = (self.hb_df['criteria'] == 1)
        df = self.hb_df[mask]
        df = df.reset_index()
        self.clean_data_ratio = df.shape[0] * 100. / self.hb_df.shape[0]
        print('After cleaning, {0} % data left.'.format(round(self.clean_data_ratio, 2)))
        return df

    def write_clean_traj(self, sele_indexes):
        t = np.array(sele_indexes)
        t = t.astype('int')
        t = t - 1
        t = list(t)
        u = Universe(self.inp_gro, self.inp_xtc)
        with Writer(self.out_xtc, u.trajectory.n_atoms) as w:
            for ts in u.trajectory:
                if ts.frame in t:
                    w.write(ts)
        print('Write {0} ...'.format(self.out_xtc))
        copyfile(self.inp_gro, self.out_gro)
        print('cp {0} {1}'.format(self.inp_gro, self.out_gro))
        print('vmd -gro {0} {1}'.format(self.out_gro, self.out_xtc))


class SplitTrajAgent:
    def __init__(self, host, type_na):
        self.host = host
        self.type_na = type_na
        self.host_folder = path.join(root_folder, host)
        self.na_folder = path.join(self.host_folder, type_na)
        self.idea_na_folder = path.join(idea_cmm_folder, host, type_na)
        self.split_traj_folder = path.join(self.na_folder, 'split_traj')

        self.inp_gro = path.join(self.idea_na_folder, 'input', 'allatoms', '{0}.npt4.all.gro'.format(type_na))
        self.inp_xtc = path.join(self.idea_na_folder, 'input', 'allatoms', '{0}.all.xtc'.format(type_na))
        #   self.out_gro = path.join(self.split_traj_folder, '{0}.npt4.all.gro'.format(type_na))
        #   self.out_xtc = path.join(self.split_traj_folder, '{0}.all.xtc'.format(type_na))

        self.make_folders()

    def make_folders(self):
        folders = [self.split_traj_folder]
        for folder in folders:
            miscell.check_dir_exist_and_make(folder)

    def get_split_indexes(self, n_split=5):
        u = Universe(self.inp_gro, self.inp_xtc)
        n_frames = len(u.trajectory)
        d_index = dict()
        for i in range(n_split):
            d_index[i] = range(i, n_frames, n_split)
        return d_index

    def write_split_traj(self, sele_indexes, split_id):
        out_gro = path.join(self.split_traj_folder, '{0}.npt4.all.split{1}.gro'.format(self.type_na, split_id))
        out_xtc = path.join(self.split_traj_folder, '{0}.all.split{1}.xtc'.format(self.type_na, split_id))
        t = np.array(sele_indexes)
        t = t.astype('int')
        t = t - 1
        t = list(t)
        u = Universe(self.inp_gro, self.inp_xtc)
        with Writer(out_xtc, u.trajectory.n_atoms) as w:
            for ts in u.trajectory:
                if ts.frame in t:
                    w.write(ts)
        print('Write {0} ...'.format(out_xtc))
        copyfile(self.inp_gro, out_gro)
        print('cp {0} {1} ... done'.format(self.inp_gro, out_gro))
        print('vmd -gro {0} {1}'.format(out_gro, out_xtc))


def make_df_by_hbond(filename, seq):
    d = {'time': [], 'hb_num': []}
    pattern = '# Time ='
    f = open(filename, 'r')
    lines = f.readlines()
    for i, line in enumerate(lines):
        if re.match(pattern, line):
            t = re.sub('# Time = ', '', line)
            d['time'].append(float(t))
            hb_list = list()
            for j in range(len(seq)):
                temp = lines[i+j+2].rstrip()
                temp = temp[0]
                hb_list.append(temp)
            d['hb_num'].append(''.join(hb_list))
    f.close()
    df = pd.DataFrame(d)
    return df


def weak_hb(row, ref):
    x = row['hb_num']
    if not isinstance(x, str):
        x = str(x)
    x = x[:-2]
    for num, num_ref in zip(x, ref):
        if int(num) < (int(num_ref)-1):
            return 0
    return 1

