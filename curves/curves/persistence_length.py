from os import path
from curves.curves_main_util import PreliminaryAgent
from miscell.file_util import check_dir_exist_and_make

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

    def discretize_haxis_to_pdb(self):
        pass


class Compressor(Discretizer):
    pass

    
