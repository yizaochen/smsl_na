from os import path
from miscell.file_util import check_dir_exist_and_make, check_file_exist, copy_verbose

class FoldersBuilder:

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

    def vmd_check_gro_xtc(self):
        cmd = f'vmd -gro {self.npt4_gro} {self.raw_xtc}'
        print(cmd)