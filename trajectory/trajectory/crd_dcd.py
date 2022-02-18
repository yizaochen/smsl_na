from os import path, system
from miscell.file_util import check_dir_exist_and_make, copy_verbose
from miscell.na_seq import d_sequences
from miscell.na_bp import d_n_bp
from pdb_util.atom import Atom


class AvgcrdFitdcd:
    def __init__(self, rootfolder, suffix, host):
        self.rootfolder = rootfolder # '/home/ytcdata/traj/'
        self.suffix = suffix # 'test', 'before_clean', 'after_clean'
        self.host = host
        self.n_bp = d_n_bp[self.host]
        self.type_na = self.get_type_na()
        self.n_frame = self.get_n_frame()

        self.host_folder = path.join(self.rootfolder, self.host)
        self.suffix_folder = path.join(self.host_folder, self.suffix)

        self.aa_folder = path.join(self.suffix_folder, 'allatoms')
        self.heavy_folder = path.join(self.suffix_folder, 'heavyatoms')
        self.inp_folder = path.join(self.suffix_folder, 'charmm_inp')
        self.dat_folder = path.join(self.suffix_folder, 'charmm_dat')
        self.mkcrd_folder = path.join(self.suffix_folder, 'make_crd')

        self.raw_xtc = path.join(self.aa_folder, 'raw.xtc')
        self.raw_dcd = path.join(self.aa_folder, 'raw.dcd')
        self.npt4_gro = path.join(self.aa_folder, 'npt4.gro')
        self.npt4_pdb = path.join(self.aa_folder, 'npt4.pdb')
        self.xtc2dcd_tcl = path.join(self.aa_folder, 'xtc2dcd.tcl')

        self.strand1_pdb = path.join(self.mkcrd_folder, 'strand1.pdb')
        self.strand2_pdb = path.join(self.mkcrd_folder, 'strand2.pdb')

    def main(self, exe_build_folder, exe_vmd_check_input_gro_xtc, exe_setup_pdb_dcd, exe_reset_resid, exe_pdb2crd, exe_remove_hydrogen, exe_make_avg_crd, exe_make_fit_dcd, exe_vmd_check_output_crd_dcd):
        self.build_folder(exe_build_folder)
        self.vmd_check_input_gro_xtc(exe_vmd_check_input_gro_xtc)
        self.setup_pdb_dcd(exe_setup_pdb_dcd) # split two strands
        self.reset_resid(exe_reset_resid)
        self.pdb2crd(exe_pdb2crd) # Make CRD (combine two strands) 
        self.remove_hydrogen(exe_remove_hydrogen) # Make CRD and DCD without hydrogen atoms
        self.make_avg_crd(exe_make_avg_crd) # Make Average CRD 
        self.make_fit_dcd(exe_make_fit_dcd) # fitting no-H dcd to average crd
        self.vmd_check_output_crd_dcd(exe_vmd_check_output_crd_dcd)

    def build_folder(self, exec):
        if exec:
            agent = FoldersBuilder(self.rootfolder, self.suffix, self.host)
            agent.initialize_folders()
            agent.copy_input_gro_xtc() # This should be replaced when the locations of input gro and xtc change

    def vmd_check_input_gro_xtc(self, exec):
        if exec:
            print(f'vmd -gro {self.npt4_gro} {self.raw_xtc}')

    def setup_pdb_dcd(self, exec):
        if exec:
            agent = Preliminary(self.rootfolder, self.suffix, self.host)
            agent.gro2pdb()
            agent.split_pdb_to_2_strands()
            agent.xtc2dcd()

    def reset_resid(self, exec):
        if exec:
            agent = CRDMaker(self.rootfolder, self.suffix, self.host)
            agent.reset_na2_pdb_resid()

    def pdb2crd(self, exec):
        if exec:
            pass

    def remove_hydrogen(self, exec):
        if exec:
            pass

    def make_avg_crd(self, exec):
        if exec:
            pass

    def make_fit_dcd(self, exec):
        if exec:
            pass

    def vmd_check_output_crd_dcd(self, exec):
        if exec:
            pass

    def get_type_na(self):
        if self.host in ['pnas_rna']:
            return 'arna+arna'
        else:
            return 'bdna+bdna'

    def get_n_frame(self):
        # determined by suffix, host
        if self.suffix == 'test':
            return 10
        elif self.suffix == 'before_clean':
            d_n_frame =  {'pnas_dna': 10000, 'pnas_rna': 10000,  'a_tract_21mer': 50000, 'g_tract_21mer': 50000, 
            'atat_21mer': 50000, 'gcgc_21mer': 50000, 'ctct_21mer': 50000, 'tgtg_21mer': 50000}
            return d_n_frame[self.host]
        else: # self.suffix == 'after_clean'
            d_n_frame =  {'pnas_dna': 10000, 'pnas_rna': 10000,  'a_tract_21mer': 50000, 'g_tract_21mer': 50000, 
            'atat_21mer': 50000, 'gcgc_21mer': 50000, 'ctct_21mer': 50000, 'tgtg_21mer': 50000}
            return d_n_frame[self.host]
    
class FoldersBuilder(AvgcrdFitdcd):

    def initialize_folders(self):
        for folder in [self.rootfolder, self.host_folder, self.suffix_folder, self.aa_folder, self.heavy_folder, self.inp_folder, self.dat_folder, self.mkcrd_folder]:
            check_dir_exist_and_make(folder)

    def copy_input_gro_xtc(self):
        old_rootfolder = '/home/yizaochen/codes/dna_rna/all_systems'
        old_gro = path.join(old_rootfolder, self.host, self.type_na, 'input', f'{self.type_na}.npt4.all.gro')
        old_xtc = path.join(old_rootfolder, self.host, self.type_na, 'input', f'{self.type_na}.all.xtc')
        copy_verbose(old_gro, self.npt4_gro)
        copy_verbose(old_xtc, self.raw_xtc)

class Preliminary(AvgcrdFitdcd):
    gmx = '/usr/bin/gmx'
    vmd = '/home/yizaochen/opt/vmd/bin/vmd'

    def gro2pdb(self):
        system(f'{self.gmx} editconf -f {self.npt4_gro} -o {self.npt4_pdb}')

    def split_pdb_to_2_strands(self):
        n_atom_lst = StrandNatom(self.host, self.type_na).get_n_atom_lst()
        pdb_lst = [self.strand1_pdb, self.strand2_pdb]
        with open(self.npt4_pdb, 'r') as f:
            lines = f.readlines()
        line_id = 4
        for n_atom, pdb in zip(n_atom_lst, pdb_lst):
            FormatProcessor(pdb, lines[line_id:line_id+n_atom]).gromacs_format_to_charmm_format()

    def xtc2dcd(self):
        # generate temp script
        with open(self.xtc2dcd_tcl, 'w') as tcl:
            tcl.write(f'mol new {self.npt4_pdb} type pdb\n')
            tcl.write(f'mol addfile {self.raw_xtc} type xtc first 0 last -1 step 1 waitfor all 1')
            tcl.write(f'animate write dcd {self.raw_dcd} beg 1 end {self.n_frame+1} waitfor all\n')
            tcl.write(f'exit\n')
        system(f'{self.vmd} -dispdev text -e {self.xtc2dcd_tcl}')

class FormatProcessor:
    # gromacs-format pdb to charmm-format pdb
    # based on https://github.com/yizaochen/fluctmatch/blob/master/shell_scripts/pdb_gro2charmm.sh

    def __init__(self, lines, pdb_out):
        self.lines = lines
        self.pdb_out = pdb_out

    def gromacs_format_to_charmm_format(self):
        atomgroup = self.get_atomgroup(self.lines)
        for atom in atomgroup:
            self.reset_resname(atom)
            self.reset_name(atom)
        self.write_pdb(self.pdb_out, atomgroup, verbose=True)

    @staticmethod
    def reset_resname(atom):
        if atom.resname in ['DA5P', 'RA5P', 'DA5', 'RA5', 'DA3', 'RA3', 'DA', 'RA']:
            atom.set_resname('ADE')
        elif atom.resname in ['DT5P', 'RT5P', 'DT5', 'RT5', 'DT3', 'RT3', 'DT', 'RT']:
            atom.set_resname('THY')
        elif atom.resname in ['DC5P', 'RC5P', 'DC5', 'RC5', 'DC3', 'RC3', 'DC', 'RC']:
            atom.set_resname('CYT')
        elif atom.resname in ['DG5P', 'RG5P', 'DG5', 'RG5', 'DG3', 'RG3', 'DG', 'RG']:
            atom.set_resname('GUA')
        elif atom.resname in ['DU5P', 'RU5P', 'DU5', 'RU5', 'DU3', 'RU3', 'DU', 'RU']:
            atom.set_resname('URA')

    @staticmethod
    def reset_name(atom):
        if atom.name in ["1H5'", "2H5'", "1H2'", "2H2'", "2HO'", "1H6", "2H6", "1H4", "2H4"]:
            d_atomname_map = {"1H5'": "H5'", "2H5'": "H5''", "1H2'": "H2'",  "2H2'": "H2''", 
                              "2HO'": "H2''", "1H6": "H61", "2H6": "H62",  "1H4": "H41", "2H4": "H42"}
            atom.set_name(d_atomname_map[atom.name])

    @staticmethod
    def get_atomgroup(lines):
        atomgroup = list()
        for line in lines:
            atomgroup.append(Atom(line, segid_exist=False))
        return atomgroup
    
    @staticmethod
    def write_pdb(pdb, atomgroup, verbose=False):
        with open(pdb, 'w') as f:
            f.write('REMARK\n')
            for atom in atomgroup:
                strout = atom.get_format_str_pdb()
                f.write(f'{strout}\n')
            f.write('END')
        if verbose:
            print(f'Write PDB: {pdb}')
        

class StrandNatom:

    def __init__(self, host, type_na):
        self.host = host
        self.type_na = type_na
        self.n_atom_lst = None

    def get_n_atom_lst(self):
        if self.n_atom_lst is None:
            self.set_n_atom_lst()
        return self.n_atom_lst

    def set_n_atom_lst(self):
        d_5_prime, d_central, d_3_prime = self.get_three_dict()
        n_atom_lst = list()
        for strand_name in ['guide', 'target']:
            seq = d_sequences[self.host][strand_name]
            n_atom = 0
            n_atom += d_5_prime[seq[0]]
            for nt in seq[1:-1]:
                n_atom += d_central[nt]
            n_atom += d_3_prime[-1]
            n_atom_lst.append(n_atom)
        self.n_atom_lst = n_atom_lst

    def get_three_dict(self):
        if self.type_na == 'bdna+bdna':
            d_5_prime = {'A': 30, 'T': 30, 'C': 28, 'G': 31}
            d_central = {'A': 32, 'T': 32, 'C': 30, 'G': 33}
            d_3_prime = {'A': 33, 'T': 33, 'C': 31, 'G': 34}
        else:
            d_5_prime = {'A': 31, 'U': 28, 'C': 29, 'G': 32}
            d_central = {'A': 33, 'U': 30, 'C': 31, 'G': 34}
            d_3_prime = {'A': 34, 'U': 31, 'C': 32, 'G': 35}
        return d_5_prime, d_central, d_3_prime

class CRDMaker(AvgcrdFitdcd):
    def make_input(self, amber=False, firstter=None, lastter=None): 
        if self.type_na == 'arna+arna':
            na = 'arna'
            supplement1 = None
            supplement2 = None 
        else:
            na = 'bdna'
            supplement1 = get_patch(self.seq1, 1)
            supplement2 = get_patch(self.seq2, 2)

        crd1 = path.join(self.mkcrd_folder, '{0}1.crd'.format(na))
        inp1 = Script(path.join(self.mkcrd_folder, '{0}1.inp'.format(na)))
        inp1.write_bomlev()
        inp1.initialize_rtf_prm(amber=amber)
        inp1.write_seq(self.seq1, firstter=firstter, lastter=lastter, segid='strand1')
        if supplement1 is not None:
            inp1.write_supplement(supplement1)
        inp1.gen_angle_dihedral()
        inp1.read_pdb(path.join(self.mkcrd_folder, '{0}1.1.pdb'.format(na)))
        inp1.write_crd(crd1)
        inp1.end()

        crd2 = path.join(self.mkcrd_folder, '{0}2.crd'.format(na))
        inp2 = Script(path.join(self.mkcrd_folder, '{0}2.inp'.format(na)))
        inp2.write_bomlev()
        inp2.initialize_rtf_prm(amber=amber)
        inp2.write_seq(self.seq2, firstter=firstter, lastter=lastter, segid='strand2')
        if supplement2 is not None:
            inp2.write_supplement(supplement2)
        inp2.gen_angle_dihedral()
        inp2.read_pdb(path.join(self.mkcrd_folder, '{0}2.1.pdb'.format(na)))
        inp2.write_crd(crd2)
        inp2.end()

        inp3 = Script(path.join(self.mkcrd_folder, '{0}.inp'.format(self.type_na)))
        inp3.write_bomlev()
        inp3.initialize_rtf_prm(amber=amber)
        inp3.write_seq(self.seq1, firstter=firstter, lastter=lastter, segid='strand1')
        if supplement1 is not None:
            inp3.write_supplement(supplement1)
        inp3.write_seq(self.seq2, firstter=firstter, lastter=lastter, segid='strand2')
        if supplement2 is not None:
            inp3.write_supplement(supplement2)
        inp3.gen_angle_dihedral()
        inp3.read_crd(crd1, selection='segid strand1', ignore=True)
        inp3.read_crd(crd2, selection='segid strand2', ignore=True)
        inp3.write_crd(path.join(self.mkcrd_folder, '{0}.crd'.format(self.type_na)))
        inp3.end()

    def make_crd(self):
        if self.type_na == 'arna+arna':
            na = 'arna'
        else:
            na = 'bdna'

        inp1 = path.join(self.mkcrd_folder, '{0}1.inp'.format(na))
        inp1_dat = path.join(self.mkcrd_folder, '{0}1.dat'.format(na))
        exec_charmm(inp1, inp1_dat)

        inp2 = path.join(self.mkcrd_folder, '{0}2.inp'.format(na))
        inp2_dat = path.join(self.mkcrd_folder, '{0}2.dat'.format(na))
        exec_charmm(inp2, inp2_dat)

        inp3 = path.join(self.mkcrd_folder, '{0}.inp'.format(self.type_na))
        inp3_dat = path.join(self.mkcrd_folder, '{0}.dat'.format(self.type_na))
        exec_charmm(inp3, inp3_dat)

    def reset_na2_pdb_resid(self, offset):
        offset = -self.n_bp
        pdb_name = path.join(self.mkcrd_folder, 'bdna2.1.pdb')
        f_backup = path.join(self.mkcrd_folder, 'bdna2.1.backup.pdb')
        copyfile(pdb_name, f_backup)
        print(f'{pdb_name} {f_backup}')
        reader = PDB.PDBReader(pdb_name, skip_header=2, skip_footer=1)
        for atom in reader.atomgroup:
            resid = atom.resid
            atom.set_resid(resid + offset)
        writer = PDB.PDBWriter(pdb_name, reader.atomgroup)
        writer.write_pdb()
        print(f'Reset {pdb_name} resid by offset {offset}!')
        print(f'Check by...\nvim {pdb_name}')