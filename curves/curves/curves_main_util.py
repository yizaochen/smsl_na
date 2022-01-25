from os import path, remove
import subprocess
from glob import glob
from shutil import move
import MDAnalysis as mda
from miscell.file_util import check_dir_exist_and_make, check_file_exist, copy_verbose
from miscell.na_bp import d_n_bp, d_type_na
from pdb_util.pdb import PDBReader, PDBWriter


class PreliminaryAgent:
    d_new_resname = {'RA5': 'A', 'RA3': 'A', 'RA': 'A',
                     'RU5': 'U', 'RU3': 'U', 'RU': 'U',
                     'RC5': 'C', 'RC3': 'C', 'RC': 'C',
                     'RG5': 'G', 'RG3': 'G', 'RG': 'G'}

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

    def check_rna_with_modify_resname(self):
        if self.is_rna():
            copy_verbose(self.input_pdb, self.input_pdb_backup)
            reader = PDBReader(self.input_pdb, segid_exist=False)
            atgs = reader.get_atomgroup()
            for atom in atgs:
                self.modify_rna_resname(atom)
            writer = PDBWriter(self.input_pdb, atgs)
            writer.write_pdb()

    def is_rna(self):
        type_na = d_type_na[self.host]
        return type_na == 'arna+arna'

    def modify_rna_resname(self, atom):
        new_resname = self.d_new_resname[atom.resname]
        atom.set_resname(new_resname)

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

class ExtractPDBAgent(PreliminaryAgent):
    def __init__(self, rootfolder, host):
        super().__init__(rootfolder, host)
        self.u = mda.Universe(self.input_pdb, self.input_xtc)

    def get_n_frames(self):
        return len(self.u.trajectory)

    def print_n_frames(self):
        n_frame = self.get_n_frames()
        print(f'n_frame: {n_frame}')
        
    def extract_pdb_from_xtc(self, start_frame, stop_frame):
        for ts in self.u.trajectory[start_frame:stop_frame]:
            pdb_out = path.join(self.single_pdb_folder, f'{ts.frame}.pdb')
            self.process_single_frame(pdb_out)
            if ts.frame % 500 == 0:
                print(f'Extract PDB for {self.host}, Frame-ID: {ts.frame}')

    def process_single_frame(self, pdb_out):
        with mda.Writer(pdb_out, bonds=None, n_atoms=self.u.atoms.n_atoms) as PDBOUT:
            PDBOUT.write(self.u.atoms)

class ExecCurvesAgent(ExtractPDBAgent):
    lis_name = 'r+bdna'

    def execute_curve_plus(self, start_frame, stop_frame):
        for frame_id in range(start_frame, stop_frame):
            self.clean_files()
            f_single_pdb = path.join(self.single_pdb_folder, f'{frame_id}.pdb')  
            # Start to execute curves
            cmd = self.get_exectue_curve_plus_cmd(f_single_pdb)
            errlog = open(path.join(self.workspace_folder, 'err.log'), 'w')
            outlog = open(path.join(self.workspace_folder, 'out.log'), 'w')
            subprocess.run(cmd, shell=True, stdout=outlog, stderr=errlog,check=False)
            errlog.close()
            outlog.close()
            # Store .lis and _X.pdb files
            workspace_lis = path.join(self.workspace_folder, f'{self.lis_name}.lis')
            workspace_pdb = path.join(self.workspace_folder, f'{self.lis_name}_X.pdb')
            f_lis = path.join(self.lis_folder, f'{frame_id}.lis')
            f_x_pdb = path.join(self.haxis_smooth_folder, f'{frame_id}.pdb')
            move(workspace_lis, f_lis)
            move(workspace_pdb, f_x_pdb)
            if frame_id % 500 == 0:
                print(f'Curves+ for {self.host}, Frame-ID: {frame_id}')

    def clean_files(self):
        pathname = path.join(self.workspace_folder, f'{self.lis_name}*')
        filelist = glob(pathname)
        if len(filelist) != 0:
            for fname in filelist:
                remove(fname)

    def get_exectue_curve_plus_cmd(self, f_single_pdb):
        curve = '/home/yizaochen/opt/curves+/Cur+'
        inp_end_txt = self.get_inp_end(f_single_pdb)
        n1, n2, n3, n4 = self.get_four_numbers()
        cmd1 = f'{curve} <<!\n'
        cmd2 = f'  &inp {inp_end_txt}&end\n'
        cmd3 = '2 1 -1 0 0\n'
        cmd4 = f'{n1}:{n2}\n'
        cmd5 = f'{n3}:{n4}\n'
        cmd6 = '!'
        return cmd1 + cmd2 + cmd3 + cmd4 + cmd5 + cmd6

    def get_inp_end(self, f_single_pdb):
        curve_folder = '/home/yizaochen/opt/curves+/standard'
        lis = path.join(self.workspace_folder, self.lis_name)
        return f'file={f_single_pdb},\n  lis={lis},\n  lib={curve_folder},\n naxlim=3'
    
    def get_four_numbers(self):
        return 1, self.n_bp, 2*self.n_bp, self.n_bp+1

        