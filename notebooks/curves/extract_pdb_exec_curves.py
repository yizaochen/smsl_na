from curves.curves import ExtractPDBAgent, ExecCurvesAgent

rootfolder = '/home/ytcdata/curves_data'
host_lst = ['pnas_dna', 'pnas_rna', 
            'a_tract_21mer', 'g_tract_21mer', 
            'atat_21mer', 'gcgc_21mer', 
            'ctct_21mer', 'tgtg_21mer']

# Notebook Version: extract_pdb.ipynb
for host in host_lst:
    e_agent = ExtractPDBAgent(rootfolder, host)
    start_frame = 0 # start from 0
    end_frame = e_agent.get_n_frames()
    e_agent.extract_pdb_from_xtc(start_frame, end_frame)

# Notebook Version: exec_curves.ipynb
for host in host_lst:
    exe_agent = ExecCurvesAgent(rootfolder, host)
    start_frame = 0 # start from 0
    end_frame = e_agent.get_n_frames()
    exe_agent.execute_curve_plus(start_frame, end_frame)