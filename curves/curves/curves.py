from os import path
from miscell.file_util import check_dir_exist_and_make, check_file_exist, copy_verbose


class PreliminaryAgent:
    def __init__(self, rootfolder, host):
        self.rootfolder = rootfolder
        self.host = host

        self.host_folder = path.join(rootfolder, host)
        self.pdb_xtc_folder = path.join(self.host_folder, 'pdb_xtc')
        self.single_pdb_folder = path.join(self.host_folder, 'single_pdb')
        self.lis_folder = path.join(self.host_folder, 'lis')
        self.haxis_smooth_folder = path.join(self.host_folder, 'haxis_smooth')
        self.workspace_folder = path.join(self.host_folder, 'workspace')

        self.input_pdb = path.join(self.pdb_xtc_folder, f'input.pdb')
        self.input_xtc = path.join(self.pdb_xtc_folder, f'input.xtc')
        self.input_pdb_exist = None
        self.input_xtc_exist = None

    def initialize_folders(self):
        for folder in [self.host_folder, self.pdb_xtc_folder, self.single_pdb_folder, self.lis_folder, 
                       self.haxis_smooth_folder, self.workspace_folder]:
            check_dir_exist_and_make(folder)

    def check_input_pdb_xtc(self):
        self.input_pdb_exist = check_file_exist(self.input_pdb)
        self.input_xtc_exist = check_file_exist(self.input_xtc)

    def copy_input_pdb_xtc(self, old_pdb, old_xtc):
        copy_verbose(old_pdb, self.input_pdb)
        copy_verbose(old_xtc, self.input_xtc)
        self.check_input_pdb_xtc()

    def vmd_check_pdb_xtc(self):
        cmd = f'vmd -pdb {self.input_pdb} {self.input_xtc}'
        print(cmd)

        