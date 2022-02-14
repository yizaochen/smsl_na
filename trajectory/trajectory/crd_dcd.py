from os import path
from miscell.file_util import check_dir_exist_and_make, check_file_exist, copy_verbose

class AvgcrdFitdcd:
    def __init__(self, rootfolder, suffix, host):
        self.rootfolder = rootfolder # '/home/ytcdata/traj/'
        self.suffix = suffix # 'test', 'before_clean', 'after_clean'
        self.host = host

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
        self.setup_pdb_dcd(exe_setup_pdb_dcd)
        self.pdb2crd(exe_pdb2crd) # Make CRD (split two strands, then combine) 
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
            pass

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

    
class FoldersBuilder(AvgcrdFitdcd):

    def initialize_folders(self):
        for folder in [self.rootfolder, self.host_folder, self.suffix_folder, self.aa_folder, self.heavy_folder, self.inp_folder, self.dat_folder, self.mkcrd_folder]:
            check_dir_exist_and_make(folder)

    def copy_input_gro_xtc(self):
        old_rootfolder = '/home/yizaochen/codes/dna_rna/all_systems'
        if self.host in ['pnas_rna']:
            type_na = 'arna+arna'
        else:
            type_na = 'bdna+bdna'
        old_gro = path.join(old_rootfolder, self.host, type_na, 'input', f'{type_na}.npt4.all.gro')
        old_xtc = path.join(old_rootfolder, self.host, type_na, 'input', f'{type_na}.all.xtc')
        copy_verbose(old_gro, self.npt4_gro)
        copy_verbose(old_xtc, self.raw_xtc)

    