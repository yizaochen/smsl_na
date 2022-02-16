from os import path, system
from miscell.file_util import check_dir_exist_and_make, copy_verbose
from miscell.na_seq import d_sequences
from pdb_util.atom import Atom


class AvgcrdFitdcd:
    def __init__(self, rootfolder, suffix, host):
        self.rootfolder = rootfolder # '/home/ytcdata/traj/'
        self.suffix = suffix # 'test', 'before_clean', 'after_clean'
        self.host = host
        self.type_na = self.get_type_na()

        self.host_folder = path.join(self.rootfolder, self.host)
        self.suffix_folder = path.join(self.host_folder, self.suffix)

        self.aa_folder = path.join(self.suffix_folder, 'allatoms')
        self.heavy_folder = path.join(self.suffix_folder, 'heavyatoms')
        self.inp_folder = path.join(self.suffix_folder, 'charmm_inp')
        self.dat_folder = path.join(self.suffix_folder, 'charmm_dat')
        self.mkcrd_folder = path.join(self.suffix_folder, 'make_crd')

        self.raw_xtc = path.join(self.aa_folder, 'raw.xtc')
        self.npt4_gro = path.join(self.aa_folder, 'npt4.gro')
        self.npt4_pdb = path.join(self.aa_folder, 'npt4.pdb')

        self.strand1_pdb = path.join(self.mkcrd_folder, 'strand1.pdb')
        self.strand2_pdb = path.join(self.mkcrd_folder, 'strand2.pdb')

    def main(self, exe_build_folder, exe_vmd_check_input_gro_xtc, exe_setup_pdb_dcd, exe_pdb2crd, exe_remove_hydrogen, exe_make_avg_crd, exe_make_fit_dcd, exe_vmd_check_output_crd_dcd):
        self.build_folder(exe_build_folder)
        self.vmd_check_input_gro_xtc(exe_vmd_check_input_gro_xtc)
        self.setup_pdb_dcd(exe_setup_pdb_dcd) # split two strands
        self.pdb2crd(exe_pdb2crd) # Make CRD (combine two strands) 
        self.remove_hydrogen(exe_remove_hydrogen) # Make CRD and DCD without hydrogen atoms
        self.make_avg_crd(exe_make_avg_crd) # Make Average CRD 
        self.make_fit_dcd(exe_make_fit_dcd) # fitting no-H dcd to average crd
        self.vmd_check_output_crd_dcd(exe_vmd_check_output_crd_dcd)

    def build_folder(self, exec):
        if exec:
            agent = FoldersBuilder(self.rootfolder, self.suffix, self.host)
            agent.initialize_folders()
            agent.copy_input_gro_xtc() # This should be replaced when the locations of input gro and xtc change

    def vmd_check_input_gro_xtc(self, exec):
        if exec:
            print(f'vmd -gro {self.npt4_gro} {self.raw_xtc}')

    def setup_pdb_dcd(self, exec):
        if exec:
            agent = Preliminary(self.rootfolder, self.suffix, self.host)
            agent.gro2pdb()
            agent.split_pdb_to_2_strands()

    def pdb2crd(self, exec):
        if exec:
            pass

    def remove_hydrogen(self, exec):
        if exec:
            pass

    def make_avg_crd(self, exec):
        if exec:
            pass

    def make_fit_dcd(self, exec):
        if exec:
            pass

    def vmd_check_output_crd_dcd(self, exec):
        if exec:
            pass

    def get_type_na(self):
        if self.host in ['pnas_rna']:
            return 'arna+arna'
        else:
            return 'bdna+bdna'
    
class FoldersBuilder(AvgcrdFitdcd):

    def initialize_folders(self):
        for folder in [self.rootfolder, self.host_folder, self.suffix_folder, self.aa_folder, self.heavy_folder, self.inp_folder, self.dat_folder, self.mkcrd_folder]:
            check_dir_exist_and_make(folder)

    def copy_input_gro_xtc(self):
        old_rootfolder = '/home/yizaochen/codes/dna_rna/all_systems'
        old_gro = path.join(old_rootfolder, self.host, self.type_na, 'input', f'{self.type_na}.npt4.all.gro')
        old_xtc = path.join(old_rootfolder, self.host, self.type_na, 'input', f'{self.type_na}.all.xtc')
        copy_verbose(old_gro, self.npt4_gro)
        copy_verbose(old_xtc, self.raw_xtc)

class Preliminary(AvgcrdFitdcd):
    gmx = '/usr/bin/gmx'

    def gro2pdb(self):
        system(f'{self.gmx} editconf -f {self.npt4_gro} -o {self.npt4_pdb}')

    def split_pdb_to_2_strands(self):
        n_atom_lst = StrandNatom(self.host, self.type_na).get_n_atom_lst()
        pdb_lst = [self.strand1_pdb, self.strand2_pdb]
        with open(self.npt4_pdb, 'r') as f:
            lines = f.readlines()
        line_id = 4
        for n_atom, pdb in zip(n_atom_lst, pdb_lst):
            FormatProcessor(pdb, lines[line_id:line_id+n_atom]).gromacs_format_to_charmm_format()


class FormatProcessor:
    # gromacs-format pdb to charmm-format pdb
    # based on https://github.com/yizaochen/fluctmatch/blob/master/shell_scripts/pdb_gro2charmm.sh

    def __init__(self, lines, pdb_out):
        self.lines = lines
        self.pdb_out = pdb_out

    def gromacs_format_to_charmm_format(self):
        atomgroup = self.get_atomgroup(self.lines)
        for atom in atomgroup:
            self.reset_resname(atom)
            self.reset_name(atom)
        self.write_pdb(self.pdb_out, atomgroup, verbose=True)

    @staticmethod
    def reset_resname(atom):
        if atom.resname in ['DA5P', 'RA5P', 'DA5', 'RA5', 'DA3', 'RA3', 'DA', 'RA']:
            atom.set_resname('ADE')
        elif atom.resname in ['DT5P', 'RT5P', 'DT5', 'RT5', 'DT3', 'RT3', 'DT', 'RT']:
            atom.set_resname('THY')
        elif atom.resname in ['DC5P', 'RC5P', 'DC5', 'RC5', 'DC3', 'RC3', 'DC', 'RC']:
            atom.set_resname('CYT')
        elif atom.resname in ['DG5P', 'RG5P', 'DG5', 'RG5', 'DG3', 'RG3', 'DG', 'RG']:
            atom.set_resname('GUA')
        elif atom.resname in ['DU5P', 'RU5P', 'DU5', 'RU5', 'DU3', 'RU3', 'DU', 'RU']:
            atom.set_resname('URA')

    @staticmethod
    def reset_name(atom):
        if atom.name in ["1H5'", "2H5'", "1H2'", "2H2'", "2HO'", "1H6", "2H6", "1H4", "2H4"]:
            d_atomname_map = {"1H5'": "H5'", "2H5'": "H5''", "1H2'": "H2'",  "2H2'": "H2''", 
                              "2HO'": "H2''", "1H6": "H61", "2H6": "H62",  "1H4": "H41", "2H4": "H42"}
            atom.set_name(d_atomname_map[atom.name])

    @staticmethod
    def get_atomgroup(lines):
        atomgroup = list()
        for line in lines:
            atomgroup.append(Atom(line, segid_exist=False))
        return atomgroup
    
    @staticmethod
    def write_pdb(pdb, atomgroup, verbose=False):
        with open(pdb, 'w') as f:
            f.write('REMARK\n')
            for atom in atomgroup:
                strout = atom.get_format_str_pdb()
                f.write(f'{strout}\n')
            f.write('END')
        if verbose:
            print(f'Write PDB: {pdb}')
        

class StrandNatom:

    def __init__(self, host, type_na):
        self.host = host
        self.type_na = type_na
        self.n_atom_lst = None

    def get_n_atom_lst(self):
        if self.n_atom_lst is None:
            self.set_n_atom_lst()
        return self.n_atom_lst

    def set_n_atom_lst(self):
        d_5_prime, d_central, d_3_prime = self.get_three_dict()
        n_atom_lst = list()
        for strand_name in ['guide', 'target']:
            seq = d_sequences[self.host][strand_name]
            n_atom = 0
            n_atom += d_5_prime[seq[0]]
            for nt in seq[1:-1]:
                n_atom += d_central[nt]
            n_atom += d_3_prime[-1]
            n_atom_lst.append(n_atom)
        self.n_atom_lst = n_atom_lst

    def get_three_dict(self):
        if self.type_na == 'bdna+bdna':
            d_5_prime = {'A': 30, 'T': 30, 'C': 28, 'G': 31}
            d_central = {'A': 32, 'T': 32, 'C': 30, 'G': 33}
            d_3_prime = {'A': 33, 'T': 33, 'C': 31, 'G': 34}
        else:
            d_5_prime = {'A': 31, 'U': 28, 'C': 29, 'G': 32}
            d_central = {'A': 33, 'U': 30, 'C': 31, 'G': 34}
            d_3_prime = {'A': 34, 'U': 31, 'C': 32, 'G': 35}
        return d_5_prime, d_central, d_3_prime