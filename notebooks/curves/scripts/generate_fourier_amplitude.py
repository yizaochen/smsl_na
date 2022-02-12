from curves.persistence_length import LNormTheta, FourierShape
rootfolder = '/home/ytcdata/curves_data'
host_lst = ['a_tract_21mer', 'g_tract_21mer', 
            'atat_21mer', 'gcgc_21mer', 
            'ctct_21mer', 'tgtg_21mer']

# Notebook Version: discretize_compress.ipynb
for host in host_lst:
    l_agent = LNormTheta(rootfolder, host)
    l_agent.make_l_norm_theta_hd5()

    shape = FourierShape(rootfolder, host)
    shape.make_fourier_amp_h5()
