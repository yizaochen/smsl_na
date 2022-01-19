from os import path
from miscell.file_util import check_dir_exist_and_make, check_file_exist, copy_verbose
from miscell.na_bp import d_n_bp
from pdb_util.pdb import PDBReader, PDBWriter


class PreliminaryAgent:
    def __init__(self, rootfolder, host):
        self.rootfolder = rootfolder
        self.host = host
        self.n_bp = d_n_bp[host]

        self.host_folder = path.join(rootfolder, host)
        self.pdb_xtc_folder = path.join(self.host_folder, 'pdb_xtc')
        self.single_pdb_folder = path.join(self.host_folder, 'single_pdb')
        self.lis_folder = path.join(self.host_folder, 'lis')
        self.haxis_smooth_folder = path.join(self.host_folder, 'haxis_smooth')
        self.workspace_folder = path.join(self.host_folder, 'workspace')

        self.input_pdb_backup = path.join(self.pdb_xtc_folder, f'input.backup.pdb')
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

    def vim_check_pdb(self):
        resid_i = self.n_bp + 1
        resid_j = self.n_bp * 2
        print(f'Please check the resid of second strand is from {resid_i} to {resid_j}')
        print('If not, please add A and B at the end of the lines for the two strands by vim command')
        print(":{line_begin},{line_end}s/$/A/g")
        print(":{line_begin},{line_end}s/$/B/g")
        print('Remember to trim the file becuase of PDBReader skip_header=1, skip_footer=1')
        print(f'vim {self.input_pdb}')

    def change_input_pdb_resid(self, execute=False):
        if execute:
            copy_verbose(self.input_pdb, self.input_pdb_backup)
            reader = PDBReader(self.input_pdb, segid_exist=True)
            atgs = reader.get_atomgroup()

            for atom in atgs:
                if atom.segid == 'B':
                    atom.resid += self.n_bp
                    
            writer = PDBWriter(self.input_pdb, atgs)
            writer.write_pdb()

        