from os import path
import MDAnalysis as mda
from curves.curves_main_util import PreliminaryAgent
from miscell.file_util import check_dir_exist_and_make
from pdb_util.atom import Atom
from pdb_util.pdb import PDBWriter
class FoldersBuilder(PreliminaryAgent):
    def __init__(self, rootfolder, host):
        super().__init__(rootfolder, host)
        self.haxis_discretize_folder = path.join(self.host_folder, 'haxis_discretize')
        self.discre_pdb_dcd_folder = path.join(self.host_folder, 'discre_pdb_dcd')
        self.bend_shape_folder = path.join(self.host_folder, 'bend_shape')

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

class Compressor(Discretizer):
    pass

    
