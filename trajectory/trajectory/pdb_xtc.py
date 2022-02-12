from os import path
from miscell.file_util import check_dir_exist_and_make, check_file_exist, copy_verbose

class BeforeClean:
    rootfolder = '/home/ytcdata/traj/'
    suffix = 'before_clean'

    def __init__(self, host, initialize=False):
        self.host = host
        self.host_folder = path.join(self.rootfolder, self.host)
        self.pdb_xtc_folder = path.join(self.host_folder, self.suffix)

        self.input_pdb = path.join(self.pdb_xtc_folder, 'input.pdb')
        self.input_xtc = path.join(self.pdb_xtc_folder, 'input.xtc')

        if initialize:
            self.initialize_folders()

    def initialize_folders(self):
        for folder in [self.rootfolder, self.host_folder, self.pdb_xtc_folder]:
            check_dir_exist_and_make(folder)

    def copy_input_pdb_xtc(self, old_pdb, old_xtc):
        copy_verbose(old_pdb, self.input_pdb)
        copy_verbose(old_xtc, self.input_xtc)
        self.check_input_pdb_xtc()

    def vmd_check_pdb_xtc(self):
        cmd = f'vmd -pdb {self.input_pdb} {self.input_xtc}'
        print(cmd)

    def check_input_pdb_xtc(self):
        self.input_pdb_exist = check_file_exist(self.input_pdb)
        self.input_xtc_exist = check_file_exist(self.input_xtc)

class AfterClean(BeforeClean):
    suffix = 'after_clean'

class Test(BeforeClean):
    suffix = 'test'