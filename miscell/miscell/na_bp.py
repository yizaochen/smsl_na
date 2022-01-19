host_lst = ['pnas_dna', 'pnas_rna', 
            'a_tract_21mer', 'g_tract_21mer', 
            'atat_21mer', 'gcgc_21mer', 
            'ctct_21mer', 'tgtg_21mer']

d_type_na = {'pnas_dna': 'bdna+bdna', 'pnas_rna': 'arna+arna', 
             'a_tract_21mer': 'bdna+bdna', 'g_tract_21mer': 'bdna+bdna', 
             'atat_21mer': 'bdna+bdna', 'gcgc_21mer': 'bdna+bdna', 
             'ctct_21mer': 'bdna+bdna', 'tgtg_21mer': 'bdna+bdna'}

d_n_bp =  {'pnas_dna': 16, 'pnas_rna': 16, 
           'a_tract_21mer': 21, 'g_tract_21mer': 21, 
           'atat_21mer': 21, 'gcgc_21mer': 21, 
           'ctct_21mer': 21, 'tgtg_21mer': 21}

d_time_interval_lst =  {'pnas_dna': ['0_1us'], 'pnas_rna': ['0_1us'], 
           'a_tract_21mer': ['0_1us', '1_2us', '2_3us', '3_4us', '4_5us'], 
           'g_tract_21mer': ['0_1us', '1_2us', '2_3us', '3_4us', '4_5us'], 
           'atat_21mer': ['0_1us', '1_2us', '2_3us', '3_4us', '4_5us'], 
           'gcgc_21mer': ['0_1us', '1_2us', '2_3us', '3_4us', '4_5us'],
           'ctct_21mer': ['0_1us', '1_2us', '2_3us', '3_4us', '4_5us'],
           'tgtg_21mer': ['0_1us', '1_2us', '2_3us', '3_4us', '4_5us']}

d_host_abbr = {'pnas_dna': 'B-DD', 'pnas_rna': 'A-RR', 
               'a_tract_21mer': 'polyA', 'g_tract_21mer': 'polyG', 
               'atat_21mer': 'TpA', 'gcgc_21mer': 'CpG', 
               'ctct_21mer': 'TpC', 'tgtg_21mer': 'GpT'}

d_color = {'pnas_dna': 'blue', 'pnas_rna': 'red', 
           'a_tract_21mer': "#5c8ecb", 'g_tract_21mer': "#ea6251", 
           'atat_21mer': "#8cf8d5", 'gcgc_21mer': "#ef85ed", 
           'ctct_21mer': "#5aed49", 'tgtg_21mer': "#f09f47"}
