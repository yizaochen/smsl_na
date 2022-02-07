from curves.curves_main_util import ExtractPDBAgent
from curves.persistence_length import FoldersBuilder, Discretizer

rootfolder = '/home/ytcdata/curves_data'
host_lst = ['pnas_rna', 
            'a_tract_21mer', 'g_tract_21mer', 
            'atat_21mer', 'gcgc_21mer', 
            'ctct_21mer', 'tgtg_21mer']

# Notebook Version: discretize_compress.ipynb
for host in host_lst:
    builder = FoldersBuilder(rootfolder, host)
    builder.initialize_folders()

    extract_agent = ExtractPDBAgent(rootfolder, host)
    start_frame = 0 # start from 0
    stop_frame = extract_agent.get_n_frames()

    discretizer = Discretizer(rootfolder, host)
    discretizer.discretize_haxis_to_pdb(start_frame, stop_frame)