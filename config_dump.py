# Config class for use when reading lammps dump files. Basic atom properties 
# are stored in list of Atom objects. Statistical analysis is performed on 
# whole configuration

# UPDATE 6/01/2014
# Now includes a data object for storing data sets and performing averages,
# variance etc. 

from numpy import matrix
from atom_dump import Atom, Molecule
from vector_operations import *
import sys


class Data_set:
# Class to be used to store data sets (profiles, etc.) to be averaged later
    def __init__(self):
        self.x_index = []
        self.y_data_sets = []
        self.y_data_ave = []
        self.y_data_stderr = []
        self.header = ""
        
    def add_data_set(self, head, x_data, y_data):

        if not self.x_index:
            self.x_index.extend(x_data)  
        if not self.header:
            self.header = head  

        self.y_data_sets.append(y_data)

    def average_data(self):
       
        N_sets = len(self.y_data_sets)
        N_points = len(self.x_index)
        y_sum = N_points*[0.0]
        y_sq_sum = N_points*[0.0]
        self.y_data_ave = N_points*[0.0]
        self.y_data_stderr = N_points*[0.0]

        print N_sets
        print N_points


        for i in range(0, len(self.x_index)):
            for j in range(0, len(self.y_data_sets)):
                y_sum[i] += self.y_data_sets[j][i]
                y_sq_sum[i] += self.y_data_sets[j][i]*self.y_data_sets[j][i]  

            self.y_data_ave[i] = y_sum[i]/N_sets
            diff_sq = (y_sq_sum[i]/N_sets - 
                       self.y_data_ave[i]*self.y_data_ave[i]) 
             
            print "{} {}".format(i, diff_sq)
                  
            self.y_data_stderr[i] = math.sqrt((y_sq_sum[i]/N_sets - 
                               self.y_data_ave[i]*self.y_data_ave[i])/N_sets) 

    def data_to_file(self, fname):
        f_out = open(fname, "w") 
        f_out.write(self.header+"\n")        

        for i in range(0, len(self.x_index)):
            f_out.write("{0:<10.6}    {1:<6.4}    {2:<6.4}\n".format(
                  self.x_index[i], self.y_data_ave[i], self.y_data_stderr[i]))
                     

class VelocityLogMol:
    """
    Class used to calculate velocity autocorellation function.
    contains mol_id and list of velocities with length equal to the
    binning window. Initialised as None
    """

    def __init__(self, MOL_ID, N_BIN):

        self.mol_id = MOL_ID
        self.v_list = [None]*N_BIN 	

    def add_new_velocity(self, v):
        found_flag = False
        i = 0
        while not found_flag and i < len(self.v_list):
            if self.v_list[i] == None:
                self.v_list[i] = v
                found_flag = True
            i += 1 
        return found_flag 

class  Hbond:
    """
    Class used to calculate "autocorellation function" for hydrogen bonds,
    based on the same criteria as for the hbond profile
    """ 
    def __init__(self, MOLI_ID, MOLJ_ID):
        self.moli_id = MOLI_ID
        self.molj_id = MOLJ_ID
        self.del_t = 0

    def increment_t(self):
        self.del_t += 1

# Dictionary of types for CG sites (v4)
class Config:

    def __init__(self):
        self.type_index = {'MPD_B':1, 'MPD_N':2, 'MPD_SA':3, 'MPD_SB':4,
                           'TMC_B':5, 'TMC_C':6, 'TMC_SA':7, 'TMC_SB':8,
                           'WALL':9}

        self.inv_type_index = {j:i for i, j in self.type_index.items()}

        self.mass_index = {1:26.0, 2:3.0, 3:6.0, 4:6.0,
                           5:26.0, 6:12.0, 7:8.0, 8:8.0, 9:1.0}

        pass        

    def new_type_list(self, type_index_new, mass_index_new):
        self.type_index = type_index_new
        self.inv_type_index = {j:i for i, j in self.type_index.items()}
        self.mass_index = mass_index_new

    def open_file(self, fname):
        self.f_dump = open(fname)
        if not self.f_dump:
            print "Unable to open file \"{}\"".format(fname)
        
    def read_config(self):
        # Read timestep
        EOF_flag = False  

        line = self.f_dump.readline()

        if not line: EOF_flag = True
        
        if not EOF_flag:
            self.atomlist = []
            if line == "ITEM: TIMESTEP\n":
                self.timestep = int(self.f_dump.readline())
                print "Timestep: {}".format(self.timestep)
                 
            # Read number of atoms    
            line = self.f_dump.readline()
            if line == "ITEM: NUMBER OF ATOMS\n":
                self.natoms = int(self.f_dump.readline())
                print "Number of atoms: {}".format(self.natoms)
            # Read box extents
            line = self.f_dump.readline()
            if line == "ITEM: BOX BOUNDS pp pp pp\n":
                split = self.f_dump.readline().split() 
                self.xlo = float(split[0])
                self.xhi = float(split[1])
                split = self.f_dump.readline().split()
                self.ylo = float(split[0])
                self.yhi = float(split[1])
                split = self.f_dump.readline().split()
                self.zlo = float(split[0])
                self.zhi = float(split[1])
            print "Cell extents: {} {}".format(self.xlo, self.xhi)
            print "              {} {}".format(self.ylo, self.yhi)
            print "              {} {}".format(self.zlo, self.zhi)
            # Read atom attribute list
            # NOTE: only "default" attributes are currently included
            # Assume same order for attributes - change later
            line = self.f_dump.readline()
            split = line.split()
            if split[1] == "ATOMS":
                # create attribute list
                split.pop(0) # remove "ITEM:"
                split.pop(0) # remove "ATOMS"
                attributes = []
                for i in split:
                    attributes.append(i)
 
                new_atomlist = []
                ## CHANGE 
                for i in range(0, self.natoms):
                    split = self.f_dump.readline().split()
                    new_atom = Atom(split[0], split[1], 
                                self.inv_type_index[int(split[2])], split[3], 
                                split[4], split[5], split[6], split[7],
                                split[8], split[9])
                    new_atomlist.append(new_atom)
            
                self.atomlist = new_atomlist 

        return EOF_flag 

    def print_xyz_from_list(self, fname, mol_id_list):
        """
        Prints CG atoms from list of mol_ids
        """
        xyz_out = open(fname+'.xyz', 'a')
        xyz_out.write(str(len(self.atomlist))+'\n') 
        xyz_out.write('Timestep'+str(self.timestep)+'\n')
        print mol_id_list
        for m in self.mollist:
            for mol_id in mol_id_list:
                if m.mol_id == mol_id:
                    for atom in m.atomlist:
                        if atom.spec == 'MPD_B': string = 'Zn'
                        if atom.spec == 'TMC_B': string = 'Ag'
                        if atom.spec == 'MPD_N': string = 'Zn'
                        if atom.spec == 'TMC_C': string = 'Ag'
                        if atom.spec == 'MPD_SA': string = 'H'
                        if atom.spec == 'MPD_SB': string = 'H'
                        if atom.spec == 'TMC_SA': string = 'H'
                        if atom.spec == 'TMC_SB': string = 'H'
                        if atom.spec == 'WALL': string = 'W'
                        xyz_out.write("{}  {}  {}  {}\n".format(string, 
                          atom.rx, atom.ry, atom.rz))
  

    def print_xyz(self,fname):
        xyz_out = open(fname+'.xyz', 'a')
        xyz_out.write(str(len(self.atomlist))+'\n') 
        xyz_out.write('Timestep '+str(self.timestep)+'\n')
        N = len(self.atomlist)
        atomlist = self.atomlist
        i = 1
        atomlist = sorted(atomlist, key=lambda atom: atom.atom_id)
        while i <= N:
            found_flag = False
            j = 0
            while not found_flag:
                if int(atomlist[j].atom_id) == i:
                    found_flag = True   
                    i += 1 
                    string = atomlist[j].spec
                    if atomlist[j].spec == 'MPD_B': string = 'C'
                    if atomlist[j].spec == 'TMC_B': string = 'O'
                    if atomlist[j].spec == 'MPD_N': string = 'N'
                    if atomlist[j].spec == 'TMC_C': string = 'C'
                    if atomlist[j].spec == 'MPD_SA': string = 'H'
                    if atomlist[j].spec == 'MPD_SB': string = 'H'
                    if atomlist[j].spec == 'TMC_SA': string = 'H'
                    if atomlist[j].spec == 'TMC_SB': string = 'H'
                    if atomlist[j].spec == 'WALL': string = 'W'
                    xyz_out.write("{}  {}  {}  {}\n".format(string, 
                       atomlist[j].rx, atomlist[j].ry, atomlist[j].rz))
                    del atomlist[j]   
                j += 1 
        print 'done' 

    def print_xyz_shift(self,fname, l_off, axis):

        

        xyz_out = open(fname+'.xyz', 'a')
        xyz_out.write(str(len(self.atomlist))+'\n') 
        xyz_out.write('Timestep '+str(self.timestep)+'\n')
        N = len(self.atomlist)
        atomlist = self.atomlist
        i = 1
        atomlist = sorted(atomlist, key=lambda atom: atom.atom_id)
        while i <= N:
            found_flag = False
            j = 0
            while not found_flag:
                if int(atomlist[j].atom_id) == i:
                    found_flag = True   
                    i += 1 
                    string = atomlist[j].spec
                    if atomlist[j].spec == 'MPD_B': string = 'N'
                    if atomlist[j].spec == 'TMC_B': string = 'O'
                    if atomlist[j].spec == 'MPD_N': string = 'N_M'
                    if atomlist[j].spec == 'TMC_C': string = 'O_T'
                    if atomlist[j].spec == 'MPD_SA': string = 'H'
                    if atomlist[j].spec == 'MPD_SB': string = 'H'
                    if atomlist[j].spec == 'TMC_SA': string = 'H'
                    if atomlist[j].spec == 'TMC_SB': string = 'H'
                    if atomlist[j].spec == 'WALL': string = 'W'
                   
                    x_j = atomlist[j].rx
                    y_j = atomlist[j].ry
                    z_j = atomlist[j].rz
                    
                    if axis=='x':
                        x_j -= l_off
                    if axis=='y':
                        y_j -= l_off
                    if axis=='z':
                        z_j -= l_off

                    xyz_out.write("{}  {}  {}  {}\n".format(string, 
                       x_j, y_j, z_j))
                    del atomlist[j]   
                j += 1 
        print 'done' 

    def print_config(self):
        for a in self.atomlist:
            a.print_atom()         

    def print_new_atomlist(self):
        for a in self.new_atomlist:
            a.print_atom()         
    
    def construct_molecule_list(self):
        self.molecule_list = []
        flag = 0
        i = 0
        j = 0
        for a in self.atom_list:
            i += 1
     #      print a.mol_id 
            flag = 0
            for m in self.molecule_list:
                j += 1
                if a.mol_id == m.mol_id:
                   m.add_atom(a)
                   flag = 1
            if not flag:
                new_molecule = Molecule(a.mol_id)
                new_molecule.add_atom(a)
                if a.spec[:3] == "TMC": new_molecule.set_type("TMC")
                if a.spec[:3] == "MPD": new_molecule.set_type("MPD")
                self.mollist.append(new_molecule)
     #      print "i = {}".format(i) 


        sys.stdout.write("There are {} molecules.\n".format(len(self.mollist)))





#   def create_mol_table(self):
#       """
#       Generates table that associates mol_id values with string label for 
#       each molecule type.
#       """
#       for m in self.mollist:
            

    def generate_atomistic_from_list(self, mol_id_list, outname):
        """
        New (temporary) function to generate atomstic configurations from 
        a subset of molecules represented by a list of mol_ids.
        """

        # Populate new sublists of molecules and atoms          
        self.mollist_tmp = []
        self.atomlist_tmp = []
        for mol in self.mollist:
            for mol_id in mol_id_list:
                if mol.mol_id == mol_id:
                    self.mollist_tmp.append(mol)
                    for atom in mol.atomlist:
                        self.atomlist_tmp.append(atom)
         
                  
        self.new_atomlist = []  
        self.new_mollist = []       
 
        self.remove_pbc_mollist(self.mollist_tmp)

        sys.stdout.write("Generating atomistic configuration.\n")
        sys.stdout.flush()

        self.add_benzene(self.mollist_tmp)
 
        self.add_nonbonding_groups(self.mollist_tmp)

        self.add_bond_sites(self.mollist_tmp, self.atomlist_tmp)
 
        R_COM = self.COM_subset(mol_id_list)

        self.print_new_to_file(outname+'.xyz', R_COM) 

    def generate_atomistic(self):
        """
        This will be main function for the generation of atomistic
        representation of coarse-grained TMC/MPD structures.

        All key (global) variables such as atomic list and bond structures will
        be defined here. All additional functions will be called from here. 
        """
       
        self.new_atomlist = []  
        self.new_mollist = []       
 
        self.remove_pbc_mollist(self.mollist)

        sys.stdout.write("Generating atomistic configuration.\n")
        sys.stdout.flush()

        self.add_benzene(self.mollist)
 
        self.add_nonbonding_groups(self.mollist)

        self.add_bond_sites()

        self.print_new_to_file('atomistic_out.xyz') 
           
    def remove_pbc_mollist(self, mollist):
        """
        Removes periodic boundary condition effects of split molecules
        over box extents
        """
        Lx = self.xhi - self.xlo
        Ly = self.yhi - self.ylo
        Lz = self.zhi - self.zlo

        l_x = matrix([[Lx, 0.0, 0.0]])
        l_y = matrix([[0.0, Ly, 0.0]])
        l_z = matrix([[0.0, 0.0, Lz]])

        for m in mollist:
            i = 0
            N = len(m.atomlist)
            while i < N:
                j = i 
                while j < N:
                    r_i = m.atomlist[i].return_r_matrix()
                    r_j = m.atomlist[j].return_r_matrix()
                    r_ij = r_i - r_j
                 #  print "i:{} j:{}".format(i, j)
                    if r_ij[0,0] > Lx/2:
                        m.atomlist[j].set_r(r_j + l_x)
                    elif r_ij[0,0] < -Lx/2:
                        m.atomlist[j].set_r(r_j - l_x)
                    if r_ij[0,1] > Ly/2:
                        m.atomlist[j].set_r(r_j + l_y)
                    elif r_ij[0,1] < -Ly/2:
                        m.atomlist[j].set_r(r_j - l_y)
                    if r_ij[0,2] > Lz/2:
                        m.atomlist[j].set_r(r_j + l_z)
                    elif r_ij[0,2] < -Lz/2:
                        m.atomlist[j].set_r(r_j - l_z)
                    j += 1
                i += 1    
  
    def return_rotation_basis(self, mol_index):
        """ 
        Returns basis vectors for rotation and translation vector 
        corresponding to orientation of CG benzene triangle 
        """
 
        mol_flag = False
        i = 0 
        CG_benzenes = [] 
        while not mol_flag:
            if self.mollist[i].mol_id == mol_index:
            
                #   target molecule found
                mol_flag = True
                m = self.mollist[i]

                for a in m.atomlist:
                    if (a.spec == 'TMC_B' or a.spec == 'MPD_B'):
                        CG_benzenes.append(a)

                # Calculate separation vectors
                r_i = CG_benzenes[0].return_r_matrix()
                r_j = CG_benzenes[1].return_r_matrix()
                r_k = CG_benzenes[2].return_r_matrix()
             
                r_ij = r_i - r_j
                r_jk = r_j - r_k
                r_ki = r_k - r_i 

                e_prime1 = r_i - r_j
                e_prime2 = r_k - (r_i + r_j)/2  
                e_prime = get_basis(e_prime1, e_prime2)    

                # Check orthogonality:
                e1e2 = float(matrix.dot(e_prime[0], e_prime[1].T))  
                e1e3 = float(matrix.dot(e_prime[0], e_prime[2].T))  
                e2e3 = float(matrix.dot(e_prime[1], e_prime[2].T))  
#               print "e1.e2 = {}".format(e1e2)
#               if abs(e1e2) > 0.1: quit()  
#               print "e1.e3 = {}".format(e1e3)
#               if abs(e1e3) > 0.1: quit()  
#               print "e2.e3 = {}".format(e2e3)
#               if abs(e2e3) > 0.1: quit()  
                
 
                r_trans = (r_i + r_j + r_k)/3  

            i += 1 
        return e_prime, r_trans     
                           


    def add_benzene(self, mollist):
        """
        Adds benzene rings to new_atomlist. Called from "generate_atomistic"
        """

        # Identify TMD/MPD benzene representations in each molecule
        for m in mollist:
            flag = False 
            CG_benzenes = []
            for a in m.atomlist:
                if (a.spec == 'TMC_B' 
                   or a.spec == 'MPD_B'):

                   CG_benzenes.append(a)
                   flag = True             
            if flag == True: 

                r_i = CG_benzenes[0].return_r_matrix()
                r_j = CG_benzenes[1].return_r_matrix()
                r_k = CG_benzenes[2].return_r_matrix()
             
                r_ij = r_i - r_j
                r_jk = r_j - r_k
                r_ki = r_k - r_i 

                e_prime1 = r_i - r_j
                e_prime2 = r_k - (r_i + r_j)/2  
                e_prime = get_basis(e_prime1, e_prime2)    
 
                r_trans = (r_i + r_j + r_k)/3       

                benzene_list = self.benzene_ring()
   
 
                benzene_list = transform_translate(e_prime, 
                                                   r_trans, benzene_list)
               
                #fix mol and atom ids later

                new_molecule = Molecule(m.mol_id)

                for b in benzene_list:
                    new_atom = Atom(0, m.mol_id, 'C_b', b[0,0], b[0,1], b[0,2],
                                0.0, 0.0, 0.0, 12.0)                 
                    new_molecule.add_atom(new_atom)
                    self.new_atomlist.append(new_atom) 

                self.new_mollist.append(new_molecule)  

    def add_nonbonding_groups(self, mollist):
        """
        Adds non-bonding groups to atomistic representation based on 
        CG side group positions.
        """ 
        sys.stdout.write("Adding non-bonding groups... ")
        sys.stdout.flush()

        for m in mollist:
            for new_m in self.new_mollist:
                if new_m.mol_id == m.mol_id: # find matching atomistic molecule
                    for a in m.atomlist:
                        group_points = []
	                if a.spec == 'MPD_N':
                            new_atom = Atom(0, m.mol_id, 'N', 
                                            a.rx, a.ry, a.rz, 
                                            0.0, 0.0, 0.0, 14.0)
                            self.new_atomlist.append(new_atom)
                            new_m.add_atom(new_atom)

	                if a.spec == 'TMC_C':
                            new_atom = Atom(0, m.mol_id, 'C', 
                                            a.rx, a.ry, a.rz, 
                                            0.0, 0.0, 0.0, 12.0)
                            self.new_atomlist.append(new_atom)
                            new_m.add_atom(new_atom)

        for m in self.new_mollist: # Add hydrogen to benzene
            for i in m.atomlist:
                l_CH = 1.09 
                neighbourlist = []
                if i.spec == 'C_b': # NOTE: later just check bonds
                    for j in m.atomlist:
                        if j != i:
                            r_ij = (i.return_r_matrix() 
                                      - j.return_r_matrix())
                            mag_r_ij = float(r_ij*r_ij.T)
                            if mag_r_ij < 1.6*1.6:
                                neighbourlist.append(j)
 
                    if len(neighbourlist) == 2:
                        r1 = (i.return_r_matrix() -
                               neighbourlist[0].return_r_matrix()) 
                        r2 = (i.return_r_matrix() -
                               neighbourlist[1].return_r_matrix()) 
                        r_CH = (l_CH*(r1+r2)/
                                math.sqrt((r1+r2)*(r1+r2).T))
                        r_H = i.return_r_matrix() + r_CH
                        new_atom = Atom(0, m.mol_id, 'H', 
                                        r_H[0,0], r_H[0,1], r_H[0,2], 
                                        0.0, 0.0, 0.0, 1.0)
                        self.new_atomlist.append(new_atom) 
                        new_m.add_atom(new_atom)   
                                  
        print "Done."

    def add_bond_sites(self, mollist, atomlist):
        """
        Loops over molecules in CG representation and adds atoms to bonding 
        sites in atomistic representation.
        """
        sys.stdout.write("Adding bonded sites... ")
        sys.stdout.flush() 

        l_bond = 0.35
        l_ab = 1.50
        l_sbN = 1.0
        l_sbC = l_sbN
        l_sbB = 2.5
        l_NH = 1.01
        l_CO = 1.19
        l_CCl = 1.82 
        l_NB = 1.6
        l_CB = l_NB
        rt3 = math.sqrt(3.0)

        Lx = self.xhi - self.xlo
        Ly = self.yhi - self.ylo
        Lz = self.zhi - self.zlo

        l_x = matrix([[Lx, 0.0, 0.0]])
        l_y = matrix([[0.0, Ly, 0.0]])
        l_z = matrix([[0.0, 0.0, Lz]])

        for mi in mollist:
            for ai in mi.atomlist:
                bond_flag = 0
                for mj in mollist:
                    for aj in mj.atomlist:
                    #   print "species i: {} species j: {}".format(
                     #  ai.spec, aj.spec)
                        if ((ai.spec == "MPD_SA" and aj.spec == "TMC_SA") 
                        or (ai.spec == "TMC_SA" and aj.spec == "MPD_SA")):
                            r_i = ai.return_r_matrix()
                            r_j = aj.return_r_matrix()
                            # PBCs here!
                            r_ij = r_i - r_j
                            if r_ij[0,0] > Lx/2: r_ij -= l_x 
                            if r_ij[0,0] < -Lx/2: r_ij += l_x 
                            if r_ij[0,1] > Ly/2: r_ij -= l_y 
                            if r_ij[0,1] < -Ly/2: r_ij += l_y 
                            if r_ij[0,2] > Lz/2: r_ij -= l_z 
                            if r_ij[0,2] < -Lz/2: r_ij += l_z 
                            r_ij2 = float(r_ij*r_ij.T)
                              
                            if r_ij2 < l_bond*l_bond:
                                if (aj.spec == "MPD_SA"):
                                    mj_index = mj.mol_id            
                         #          print "TMC index m = {}".format(mi.mol_id)
                         #          print "MPD index m = {}".format(mj_index)
                                    aj_index = aj.atom_id
                         #          print "MPD index a = {}".format(aj_index)
                                bond_flag = 1

                if (bond_flag and ai.spec ==  "TMC_SA"):
             #      print "TMC index m = {}".format(mi.mol_id)
                    # Find corresponding MPD site
                    j = 0 
                    while (mj_index != mollist[j].mol_id):
                        j += 1
                    mj = mollist[j]
                    mj.mol_id
             #      print "MPD index m = {}".format(mj.mol_id)
                    j = 0
           #        print self.atomlist[j].atom_id 
                    while (aj_index != atomlist[j].atom_id):
                         
                        j += 1
                    aj = atomlist[j]
             #      print "MPD index a = {}".format(aj.mol_id)    
                    r_j = aj.return_r_matrix()
                    
                    for ak in mi.atomlist: 
                        r_k = ak.return_r_matrix()
                        r_ik = r_i - r_k
                        r_ik2 = float(r_ik*r_ik.T)
                        if (ak.spec == "TMC_C" and r_ik2 < l_sbC*l_sbC):
                            r_C = r_k
                      #     print "in TMC loop 1" 
                            for am in mi.atomlist:
                                r_m = am.return_r_matrix()
                                r_Cm = r_C - r_m
                                r_Cm2 = float(r_Cm*r_Cm.T)
                 #              print "r_Cm = {}, l_CB = {}".format(
                 #                                  math.sqrt(r_Cm2), l_CB)
                                if (am.spec == "TMC_B" and r_Cm2 < l_CB*l_CB):
               #                    print "            2" 
                                    r_B = am.return_r_matrix()
                                    r_CB = r_C - r_B
                                                                      
                    for al in mj.atomlist: 
                        r_l = al.return_r_matrix()
                        r_jl = r_j - r_l
                        r_jl2 = float(r_jl*r_jl.T)
                        if (al.spec == "MPD_N" and r_jl2 < l_sbN*l_sbN):
                            r_N = r_l
               #            print "in MPD loop 1" 
                            for an in mj.atomlist:
                                r_n = an.return_r_matrix()
                                r_Nn = r_N - r_n
                                r_Nn2 = float(r_Nn*r_Nn.T)
                                if (an.spec == "MPD_B" and r_Nn2 < l_NB*l_NB):
                #                   print "            2" 
                                    r_B = an.return_r_matrix()
                                    r_NB = r_N - r_B
                     
                    r_CN = r_C - r_N 
                    # Apply pbc
                    if r_CN[0,0] > Lx/2: r_CN -= l_x 
                    if r_CN[0,0] < -Lx/2: r_CN += l_x 
                    if r_CN[0,1] > Ly/2: r_CN -= l_y 
                    if r_CN[0,1] < -Ly/2: r_CN += l_y 
                    if r_CN[0,2] > Lz/2: r_CN -= l_z 
                    if r_CN[0,2] < -Lz/2: r_CN += l_z 
                    # Set up Oxygen
                    e1 = r_CB
                    e1 /= math.sqrt(float(e1*e1.T))
                    e2 = get_normal(e1, r_CN)
                    e_prime = get_basis(e1, e2)
                    r_trans = r_C
                    O = []
                    O.append(matrix([[l_CO/2, 0.0, -l_CO*rt3/2]]))
 
                    new_O = transform_translate(e_prime, r_trans, O)
         
                    new_atom = Atom(0, mi.mol_id, 'O_peptide', 
                                   new_O[0][0,0], new_O[0][0,1], new_O[0][0,2],
                                   0.0, 0.0, 0.0, 1.0)
               
                    for i in range(0, len(self.new_mollist)):
                          if (mi.mol_id == self.new_mollist[i].mol_id):
                              index = i
          #         print "index = {}".format(index) 
                    self.new_mollist[index].add_atom(new_atom)
                    self.new_atomlist.append(new_atom)
                          
                    # Set up Hydrogen
                    e1 = r_NB
                    e1 /= math.sqrt(float(e1*e1.T))
                    e2 = get_normal(e1, r_CN)
                    e_prime = get_basis(e1, e2)
                    r_trans = r_N
                    H = []
                    H.append(matrix([[l_NH/2, 0.0, l_NH*rt3/2]]))
  
                    # transform vectors
               
                    new_H = transform_translate(e_prime, r_trans, H)

        #           print new_H 
                    new_atom = Atom(0, mj.mol_id, 'H_peptide',
                                   new_H[0][0,0], new_H[0][0,1], new_H[0][0,2],
                                   0.0, 0.0, 0.0, 1.0)

                    for i in range(0, len(self.new_mollist)):
                          if (mj.mol_id == self.new_mollist[i].mol_id):
                              index = i
         #          print "index = {}".format(index) 
                    self.new_mollist[index].add_atom(new_atom)
                    self.new_atomlist.append(new_atom)
                    
                        

                if not bond_flag and ai.spec == "TMC_SA":
                    # Add Chlorine and oxygen to TMC Carbon
                    r_sa = ai.return_r_matrix()
                    for aj in mi.atomlist:
                        if aj.spec == "TMC_SB":
                            r_sb = aj.return_r_matrix()
                            r_ab = r_sa - r_sb
                            r_ab2 = float(r_ab*r_ab.T)
                            if r_ab2 < l_ab*l_ab: # locate adjacent bond site
                                r_C = matrix([[0,0,0]])
                                r_B = matrix([[20,20,20]])
                                for ak in mi.atomlist: 
                                    # search for nitrogen and carbon
                                    r_k = ak.return_r_matrix()
                                    r_sbk = r_sb - r_k
                                    r_sbk2 = float(r_sbk*r_sbk.T)

                                    if (ak.spec == "TMC_C" and 
                                    r_sbk2 < l_sbC*l_sbC):
                                        r_C = r_k

                                    if (ak.spec == "TMC_B" and 
                                    r_sbk2 < l_sbB*l_sbB):
                                        r_B = r_k 

                                e1 = r_C - r_B
                                r_sbC = r_sb - r_C
                                e1 /= math.sqrt(float(e1*e1.T))
                                e2 = get_normal(e1, r_sbC)
                                e_prime = get_basis(e1, e2)
                                r_trans = r_C

                                OCl_set = []
                                OCl_set.append(
                                matrix([[l_CO/2, 0.0, l_CO*rt3/2]]))
                                OCl_set.append(
                                matrix([[l_CCl/2, 0.0, -l_CCl*rt3/2]]))
                     #          H_set.append(
                      #         matrix([[l_NH, 0.0, 0.0]]))
  
                                # transform vectors
               
                                new_OCl_set = transform_translate(
                                            e_prime, r_trans, OCl_set)

                                OxygenChlorine = []
                                
                                new_atom = Atom(0, mi.mol_id, 'O',
                                               new_OCl_set[0][0,0], 
                                               new_OCl_set[0][0,1], 
                                               new_OCl_set[0][0,2],
                                               0.0, 0.0, 0.0, 1.0)
                                OxygenChlorine.append(new_atom)
                                new_atom = Atom(0, mi.mol_id, 'Cl',
                                               new_OCl_set[1][0,0], 
                                               new_OCl_set[1][0,1], 
                                               new_OCl_set[1][0,2],
                                               0.0, 0.0, 0.0, 1.0)
                                OxygenChlorine.append(new_atom)

                                for i in range(0, len(self.new_mollist)):
                                    if (mi.mol_id == 
                                    self.new_mollist[i].mol_id):
                                        index = i 
                                for atom in OxygenChlorine:
                                    self.new_mollist[index].add_atom(atom)
                                    self.new_atomlist.append(atom)
                    
                if not bond_flag and ai.spec == "MPD_SA":
                    # Add hydrogens to MPD Nitrogen
                    r_sa = ai.return_r_matrix()
                    for aj in mi.atomlist:
                        if aj.spec == "MPD_SB":
                            r_sb = aj.return_r_matrix()
                            r_ab = r_sa - r_sb
                            r_ab2 = float(r_ab*r_ab.T)
                            if r_ab2 < l_ab*l_ab: # locate adjacent bond site
                                r_N = matrix([[0,0,0]])
                                r_B = matrix([[20,20,20]])
                                for ak in mi.atomlist: 
                                    # search for nitrogen and carbon
                                    r_k = ak.return_r_matrix()
                                    r_sbk = r_sb - r_k
                                    r_sbk2 = float(r_sbk*r_sbk.T)
                                    if (ak.spec == "MPD_N" and 
                                    r_sbk2 < l_sbN*l_sbN):
                                        r_N = r_k
                                        
                                    if (ak.spec == "MPD_B" and 
                                    r_sbk2 < l_sbB*l_sbB):
                                        r_B = r_k 
 
                                e1 = r_N - r_B
                                r_sbN = r_sb - r_N 
                                e1 /= math.sqrt(float(e1*e1.T))
                                e2 = get_normal(e1, r_sbN)
          #                     print "e2 = {}".format(e2)
                                e_prime = get_basis(e1, e2)
                                r_trans = r_N
                                # set up vectors for Hydrogen
                                H_set = []
                                H_set.append(
                                matrix([[l_NH/2, 0.0, l_NH*rt3/2]]))
                                H_set.append(
                                matrix([[l_NH/2, 0.0, -l_NH*rt3/2]]))
  
                                # transform vectors
               
                                new_H_set = transform_translate(
                                            e_prime, r_trans, H_set)

                                Hydrogens = []

                                for r_H in new_H_set:
                                    new_atom = Atom(0, mi.mol_id, 'H',
                                               r_H[0,0], r_H[0,1], r_H[0,2],
                                               0.0, 0.0, 0.0, 1.0)
                                    Hydrogens.append(new_atom)

                                for i in range(0, len(self.new_mollist)):
                                    if (mi.mol_id == 
                                    self.new_mollist[i].mol_id):
                                        index = i 
                                for H in Hydrogens:
                                    self.new_mollist[index].add_atom(H)
                                    self.new_atomlist.append(H)
                                                  
        print "Done."                        
                           
    def benzene_ring(self):
        """
        Returns lists of points (matrices) for benzene ring
        """
        
        a = 1.41
        
        rt3 = math.sqrt(3.0)

        benzene = []

        benzene.append(matrix([[0.0, a, 0.0]]))  
        benzene.append(matrix([[a*rt3/2, a/2, 0.0]]))  
        benzene.append(matrix([[a*rt3/2, -a/2, 0.0]]))  
        benzene.append(matrix([[0.0, -a, 0.0]]))  
        benzene.append(matrix([[-a*rt3/2, -a/2, 0.0]]))  
        benzene.append(matrix([[-a*rt3/2, a/2, 0.0]]))

        return benzene

    def print_new_to_file(self, fname, r_shift):
        print "Printing atomistic configuration to file \"{}\".".format(fname)  
        f_output = open(fname, 'w')
        n_atom = len(self.new_atomlist)
        f_output.write("{}\n\n".format(n_atom))

        for i in self.new_atomlist:
            f_output.write("{0} {1:6f} {2:6f} {3:6f}\n".format(i.spec, 
                     i.rx-r_shift[0,0], i.ry-r_shift[0,1], i.rz-r_shift[0,2]))

    def print_to_file(self, fname):
        f_output = open(fname, 'w')
        n_atom = len(self.atomlist)
        f_output.write("{}\n\n".format(n_atom))

        for i in self.atomlist:
            name = "NULL"
    #       print "species: {}".format(i.spec) 
            if self.type_index[i.spec] == 1: name = "M" 
            if self.type_index[i.spec] == 2: name = "M" 
            if self.type_index[i.spec] == 3: name = "X" 
            if self.type_index[i.spec] == 4: name = "Y" 
            if self.type_index[i.spec] == 5: name = "T" 
            if self.type_index[i.spec] == 6: name = "T" 
            if self.type_index[i.spec] == 7: name = "X" 
            if self.type_index[i.spec] == 8: name = "Y" 
  #         if int(i.spec) == 2: name = "M" 
#           if int(i.spec) == 3: name = "X" 
 #          if int(i.spec) == 4: name = "Y" 
   #        if int(i.spec) == 5: name = "T" 
    #       if int(i.spec) == 6: name = "T" 
     #      if int(i.spec) == 7: name = "X" 
      #     if int(i.spec) == 8: name = "Y" 
            f_output.write("{0} {1:6f} {2:6f} {3:6f}\n"
                           .format(name, i.rx, i.ry, i.rz))

    def COM_velocity(self, molecule):
        """
        Returns the centre of mass velocity for a given molecule
        """
        M = 0.0
        vm = matrix([[0.0, 0.0, 0.0]])

        for atom in molecule.atomlist:
            M += atom.mass
            vm += atom.return_v_matrix()*atom.mass
        V = vm/M
 
        return V

    def COM_vector(self, molecule):
        """
        Returns the centre of mass vector for a given molecule
        """

        Lx = self.xhi - self.xlo
        Ly = self.yhi - self.ylo
        Lz = self.zhi - self.zlo

        l_x = matrix([[Lx, 0.0, 0.0]])
        l_y = matrix([[0.0, Ly, 0.0]])
        l_z = matrix([[0.0, 0.0, Lz]])


        M = 0.0
        rm = matrix([[0.0, 0.0, 0.0]]) 

 
        for atom in molecule.atomlist:
            M += atom.mass
            rm += atom.return_r_matrix()*atom.mass
        R = rm/M

        if R[0,0] < self.xlo: R += l_x 
        if R[0,0] > self.xhi: R -= l_x 
        if R[0,1] < self.ylo: R += l_y 
        if R[0,1] > self.yhi: R -= l_y 
        if R[0,2] < self.zlo: R += l_z 
        if R[0,2] > self.zhi: R -= l_z 
        
        return R

    def COM_cell(self):
        """
        Returns centre of mass vector of entire system
        """
        
        
        Lx = self.xhi - self.xlo
        Ly = self.yhi - self.ylo
        Lz = self.zhi - self.zlo

        l_x = matrix([[Lx, 0.0, 0.0]])
        l_y = matrix([[0.0, Ly, 0.0]])
        l_z = matrix([[0.0, 0.0, Lz]])


        M = 0.0
        rm = matrix([[0.0, 0.0, 0.0]]) 

        for atom in self.atomlist:
            M += atom.mass
            rm += atom.return_r_matrix()*atom.mass 

        R = rm/M 

        if R[0,0] < self.xlo: R += l_x 
        if R[0,0] > self.xhi: R -= l_x 
        if R[0,1] < self.ylo: R += l_y 
        if R[0,1] > self.yhi: R -= l_y 
        if R[0,2] < self.zlo: R += l_z 
        if R[0,2] > self.zhi: R -= l_z 

        return R

    def COM_subset(self, mol_id_list):
        """
        Returns the centre of mass of a given subset of molecules (supplied 
        mol_id_list)
        """
        
        Lx = self.xhi - self.xlo
        Ly = self.yhi - self.ylo
        Lz = self.zhi - self.zlo

        l_x = matrix([[Lx, 0.0, 0.0]])
        l_y = matrix([[0.0, Ly, 0.0]])
        l_z = matrix([[0.0, 0.0, Lz]])


        M = 0.0
        rm = matrix([[0.0, 0.0, 0.0]]) 
       
        for mol in self.mollist:
            for index in mol_id_list:
               if mol.mol_id == index:
                   for atom in mol.atomlist:
                       M += atom.mass
                       rm += atom.return_r_matrix()*atom.mass

        R = rm/M 

        if R[0,0] < self.xlo: R += l_x 
        if R[0,0] > self.xhi: R -= l_x 
        if R[0,1] < self.ylo: R += l_y 
        if R[0,1] > self.yhi: R -= l_y 
        if R[0,2] < self.zlo: R += l_z 
        if R[0,2] > self.zhi: R -= l_z 

        return R

    ### Methods used in ConstructAtomisticConfigurationDFS class        ###
    def return_all_carbons_TMC(self, TMC_mol):

        C_atoms = []
        for a in TMC_mol.atomlist:
            if a.spec == 'TMC_C':
               C_atoms.append(a)

        return C_atoms    

    def return_all_nitrogens_MPD(self, MPD_mol):

        N_atoms = []
        for a in MPD_mol.atomlist:
            if a.spec == 'MPD_N':
               N_atoms.append(a)

        return N_atoms    

    def return_bonded_carbon_TMC(self, TMC_mol, TMC_bond_atom):
        """
        Finds carbon closest to given bond site in TMC molecule.
        Returns position vector and orientation vector.
        """


        tol = 0.1
        l_Cbond = 0.69
        l_Cbond2 = (l_Cbond+tol)*(l_Cbond+tol)
        l_CB = 1.42
        l_CB2 = (l_CB+tol)*(l_CB+tol)
        r_bond = TMC_bond_atom.return_r_matrix()

        for atomi in TMC_mol.atomlist:
            if atomi.spec == 'TMC_C':
                r_i = atomi.return_r_matrix()
                r_bondi = r_i - r_bond
                if float(r_bondi*r_bondi.T) < l_Cbond2:
                    # found relevant carbon atom
                    # generate orientation vector for carbon
                    r_C = r_i


        # Find closest Benzene proxy
        for atomi in TMC_mol.atomlist:
            if atomi.spec == 'TMC_B':
                r_B_temp = atomi.return_r_matrix()
                r_CB_temp = r_C - r_B_temp 
                r_CB2 = float(r_CB_temp*r_CB_temp.T)
                if r_CB2 < l_CB2:
                    l_CB2 = r_CB2 
                    r_CB = r_CB_temp
                    r_B = r_B_temp

        r_CB /= math.sqrt(float(r_CB*r_CB.T))
        return r_CB, r_C, r_B

    def return_bonded_nitrogen_MPD(self, MPD_mol, MPD_bond_atom):
        """
        Finds nitrogen closest to given bond site in MPD molecule.
        Returns position vector and orientation vector.
        """
       
        tol = 0.1 
        l_Nbond = 0.69
        l_Nbond2 = (l_Nbond+tol)*(l_Nbond+tol)
        l_NB = 1.42
        l_NB2 = (l_NB+tol)*(l_NB+tol)
        r_bond = MPD_bond_atom.return_r_matrix()
        l_NB2 = 100000

         
        # Find nitrogen
        for atomi in MPD_mol.atomlist:
            if atomi.spec == 'MPD_N':
                r_i = atomi.return_r_matrix()
                r_bondi = r_i - r_bond
                if float(r_bondi*r_bondi.T) < l_Nbond2:
                    # found relevent nitrogen atom
                    # generate orientation vector for nitrogen
                    r_N = r_i
                   

        # Find closest Benzene proxy
        for atomi in MPD_mol.atomlist:
            if atomi.spec == 'MPD_B':
                r_B_temp = atomi.return_r_matrix()
                r_NB_temp = r_N - r_B_temp 
                r_NB2 = float(r_NB_temp*r_NB_temp.T)
                if r_NB2 < l_NB2:
                    l_NB2 = r_NB2
                    r_NB = r_NB_temp 
                    r_B = r_B_temp

        r_NB /= math.sqrt(float(r_NB*r_NB.T))

        return r_NB, r_N, r_B
 
    def return_atom(self, atom_id):
        """ 
        Returns atom with given id from self.atomlist
        """
        found_flag = 0
        new_atom = 0        

        i = 0
        while not found_flag:
            if self.atomlist[i].atom_id == atom_id:
                found_flag = 1
                new_atom = self.atomlist[i]  
            i += 1
        return new_atom

    def return_molecule(self, mol_id):
        """
        Returns molecule with given id from self.mollist
        """

        found_flag = 0
        new_molecule = 0        
       
        i = 0
        while not found_flag:
            if self.mollist[i].mol_id == mol_id:
                found_flag = 1
                new_molecule = self.mollist[i]  
            i += 1
        return new_molecule
        

    def radial_density_profile(self, mol_id_list, delr, nbin):
        """
        Generates radial density profile of a subset of molecules
        averaged over spherical shells centred at the centre of mass
        of a given subset of molecules supplied by mol_id_list.
        """
        Lx = self.xhi - self.xlo
        Ly = self.yhi - self.ylo
        Lz = self.zhi - self.zlo
 
        L = matrix([[Lx, Ly, Lz]])

        l_x = matrix([[Lx, 0.0, 0.0]])
        l_y = matrix([[0.0, Ly, 0.0]])
        l_z = matrix([[0.0, 0.0, Lz]])

        r_COM = self.COM_subset(mol_id_list)

        rho_bin = [0.0]*nbin
        r_max = nbin*delr 

        for mol in self.mollist:
            for index in mol_id_list:
                if mol.mol_id == index:
                    # Calcualte radial distance
                    for atom in mol.atomlist:
                        r_i = atom.return_r_matrix()
                        r_iCOM = r_i - r_COM
                        # Apply PBCs
                        if r_iCOM[0,0] > Lx/2: r_iCOM -= l_x 
                        if r_iCOM[0,0] < -Lx/2: r_iCOM += l_x 
                        if r_iCOM[0,1] > Ly/2: r_iCOM -= l_y 
                        if r_iCOM[0,1] < -Ly/2: r_iCOM += l_y 
                        if r_iCOM[0,2] > Lz/2: r_iCOM -= l_z 
                        if r_iCOM[0,2] < -Lz/2: r_iCOM += l_z 
                        
                        r_iCOM2 = float(r_iCOM*r_iCOM.T)
                        r_iCOM_mag = math.sqrt(r_iCOM2)
                        if r_iCOM_mag < r_max: 
                            bin_index = int(nbin*(r_iCOM_mag/r_max)) 
                            rho_bin[bin_index] += atom.mass
        # Normailse bins
        for i in range(0, nbin):
            r = i*delr
            dV = (4*math.pi/3)*((delr*delr*delr) + 3*(r*r*delr + r*delr*delr))
            print "shell = {}, volume = {}".format(i, dV)
            rho_bin[i] /= dV
        
        return rho_bin  

    def columnar_density_profile(self, mol_id_grid, gi, gj, 
                                 axis_i, axis_j, axis_k, nbin, l_k):
        """
        Generates gi*gj density profiles along the k axis. 
        Profiles are calculate with respect to the centre of mass of the 
        molecules contained within each column. The lists are supplied by
        mol_id_grid.
        """
        dq_k =  l_k/nbin

        r_com_grid = ([[ matrix([[0.0, 0.0, 0.0]]) 
                      for i in range(0, gi)] 
                      for j in range(0, gj)]) 
        rho_bin_grid = ([[ [0.0]*nbin
                      for i in range(0, gi)] 
                      for j in range(0, gj)]) 
        MPD_bin_grid = ([[ [0.0]*nbin
                      for i in range(0, gi)] 
                      for j in range(0, gj)]) 
        TMC_bin_grid = ([[ [0.0]*nbin
                      for i in range(0, gi)] 
                      for j in range(0, gj)]) 

        Lx = self.xhi - self.xlo
        Ly = self.yhi - self.ylo
        Lz = self.zhi - self.zlo
 
        L = matrix([[Lx, Ly, Lz]])
        dV = L[0, axis_i]*L[0, axis_j]*l_k/(gi*gj*nbin)   

        for i in range(0, gi):          
            for j in range(0, gj):
                print "i = {}, j = {}".format(i,j)
                r_com_grid[i][j] = self.COM_subset(mol_id_grid[i][j])         
                for mol in self.mollist:
                    for mol_id_ij in mol_id_grid[i][j]:
                        if mol.mol_id == mol_id_ij:
                            k_off = r_com_grid[i][j][0,axis_k] - 0.5*l_k
                            r_k = self.COM_vector(mol)[0,axis_k]
                            bin_index = int(nbin*(r_k - k_off)/l_k)
                            if mol.mol_type == "TMC":
                                TMC_bin_grid[i][j][bin_index] += 1/dV
                            if mol.mol_type == "MPD":
                                MPD_bin_grid[i][j][bin_index] += 1/dV
          
                            # distribute mass appropriately
                            for atom in mol.atomlist:
                                k_off = r_com_grid[i][j][0,axis_k] - 0.5*l_k
                                r_k = atom.return_r_matrix()[0,axis_k]
                                bin_index = int(nbin*(r_k - k_off)/l_k)
                                rho_bin_grid[i][j][bin_index] += atom.mass/dV
 
        return rho_bin_grid, TMC_bin_grid, MPD_bin_grid, r_com_grid                    
                              
    def track_molecules(self, axis, qlo, qhi, init_flag, mol_flag,  
                        mol_flag_init, mol_intimes, exc_list):
        """
        Routine that moniters all molecules over a defined distance in one 
        dimension in the simulation cell, giving the time interval between
        entering and leaving.
        """
        # check axis
        if axis == "x":
            L = self.xhi - self.xlo
            axis = 0
        elif axis == "y":
            L = self.yhi - self.ylo
            axis = 1
        elif axis == "z":
            L = self.zhi - self.zlo
            axis = 2
        else:
            print "An axis must be specified. x, y, or z"
            err_flg = 1

        t = self.timestep
 
        duration_list = []                
 
        if init_flag:
            # initialise molecule flags - if molecule is already in region
            # disregard data until it has left the region
            for m in self.mollist:
                index = m.mol_id - 1
  #             print index
  #             print "number of molecules = {}".format(len(self.mollist)) 
                r =  self.COM_vector(m) 
                if r[0,axis] < qhi and r[0,axis] > qlo:
                    mol_flag_init[index] = True 

        # Main loop          
        for m in self.mollist:
            index = m.mol_id - 1
            r =  self.COM_vector(m) 
            if r[0,axis] < qhi and r[0,axis] > qlo:
                # Molecule is in region
                if not mol_flag_init[index]:
                    if mol_flag[index] == False:
                        # Molecule has entered region - start tracking
                        mol_flag[index] = True
                        mol_intimes[index] = t 

            else:
                if mol_flag_init[index]:
                    # molecule initially in region has left - reset init_flag
                    mol_flag_init[index] = False
                if not mol_flag_init[index]:
                    if mol_flag[index] == True:
                        # Molecule has left region, calculate residence time 
                        mol_flag[index] = False  
                        res_time = t - mol_intimes[index]
                        exc_flag = 0 
                        # Test molecule for excluded atoms
                        for a in m.atomlist:
                            for exc_type in exc_list: 
                                if a.spec == exc_type: exc_flag = 1
                        if not exc_flag:
                            # Record time
                            duration_list.append(res_time)   

        print duration_list 
        return mol_flag, mol_flag_init, mol_intimes, duration_list                 
    def velocity_autocorrelation_region(self, xlo, xhi, ylo, yhi, zlo, zhi,
                 xyz_mask, init_flag, mol_v_log_list_old,
                 C_vv, N_vv, N_bin):
        """
        Routine to calculate the velocity autocorrelation function for all
        molecules in a region. 
        xyz_mask is a matrix identifying the components for which 
        v(0).v(t) will be calculated.
        Cvv and Nvv are arrays of length N_bin 
        """ 

        # Initialise velocity logs for particles in region
        if init_flag:
            mol_v_log_list_new = [] 
            for m in self.mollist: 
                r =  self.COM_vector(m)
                if ((r[0,0] > xlo and r[0,0] < xhi) and
                    (r[0,1] > ylo and r[0,1] < yhi) and
                    (r[0,2] > zlo and r[0,2] < zhi)):
                    
                    v_COM = self.COM_velocity(m)
                    mol_v_log = VelocityLogMol(m.mol_id, N_bin)
                    mol_v_log.add_new_velocity(v_COM)
                    mol_v_log_list_new.append(mol_v_log)  

        # Update lists             
        else:                 
            # cross reference current list - update region
            new_mollist = []
            mol_v_log_list_new = [] 
            for m in self.mollist:
                r = self.COM_vector(m)
                if ((r[0,0] > xlo and r[0,0] < xhi) and
                    (r[0,1] > ylo and r[0,1] < yhi) and
                    (r[0,2] > zlo and r[0,2] < zhi)):
 
                    v_COM = self.COM_velocity(m)
                    mol_id_current = m.mol_id 
                    found_flag = 0
                    i = 0
                    N_list = len(mol_v_log_list_old)
                    while not found_flag and i < N_list:
                        mol_id_list = mol_v_log_list_old[i].mol_id
                        if mol_id_current == mol_id_list:
                            if mol_v_log_list_old[i].add_new_velocity(v_COM):
                                mol_v_log_list_new.append(mol_v_log_list_old[i])
                                found_flag = 1

                        i += 1    
 
                    if not found_flag:
                        # Add new molecule
                        # Also adds new molecule if binning window has been 
                        # exceeded
                        mol_v_log = VelocityLogMol(m.mol_id, N_bin)
                        mol_v_log.add_new_velocity(v_COM)
                        mol_v_log_list_new.append(mol_v_log)
               

        print "Logging v(0).v(t) for {} molecules".format(
                                   len(mol_v_log_list_new)) 
        # Calculate C_vv and update N_vv
        for m in mol_v_log_list_new:
            # find latest value
            none_flag = 0 
            i = 0 
            while not none_flag and i < len(m.v_list):
                print m.v_list[i]
                if m.v_list[i] == None: 
                    none_flag = 1
                else: 
                    v_t =  m.v_list[i]*xyz_mask
                    i += 1 
 
            i -= 1

            j = 0
            while i >= 0:
                v_0 = m.v_list[i]*xyz_mask      
                vv = float(v_0*v_t.T)     
                C_vv[j] += vv
                N_vv[j] += 1.0
                j += 1
                i -= 1

        return mol_v_log_list_new, C_vv, N_vv

    def track_molecules_io(self, axis, qlo, qhi, init_flag, mol_flag,  
                        mol_flag_init, mol_flag_left, mol_intimes, exc_list):
        """
        Routine that moniters all molecules over a defined distance in one 
        dimension in the simulation cell, giving the time interval between
        entering and leaving.

        Introduce new flag for tracking from left hand side
        """
        # check axis
        if axis == "x":
            L = self.xhi - self.xlo
            axis = 0
        elif axis == "y":
            L = self.yhi - self.ylo
            axis = 1
        elif axis == "z":
            L = self.zhi - self.zlo
            axis = 2
        else:
            print "An axis must be specified. x, y, or z"
            err_flg = 1

        t = self.timestep
 
        duration_list = []                
 
        if init_flag:
            # initialise molecule flags - if molecule is already in region
            # disregard data until it has left the region
            for m in self.mollist:
                index = m.mol_id - 1
  #             print index
  #             print "number of molecules = {}".format(len(self.mollist)) 
                r =  self.COM_vector(m) 
                if r[0,axis] < qhi and r[0,axis] > qlo:
                    mol_flag_init[index] = True 

        # Main loop          
        for m in self.mollist:
            index = m.mol_id - 1
            r =  self.COM_vector(m)
            if r[0,axis] < qhi and r[0,axis] > qlo:
                # Molecule is in region
                if not mol_flag_init[index]:
                    if mol_flag[index] == False and mol_flag_left[index]==True:
                        # Molecule has entered region from the left - 
                        # start tracking
                        mol_flag[index] = True
                        mol_intimes[index] = t 

            else:
                if mol_flag_init[index]:
                    # molecule initially in region has left - reset init_flag
                    mol_flag_init[index] = False
                if not mol_flag_init[index]:
                    if mol_flag[index] == True:
                        # Molecule has left region, calculate residence time 
                        mol_flag[index] = False  
                        res_time = t - mol_intimes[index]
                        exc_flag = 0 
                        # Test molecule for excluded atoms
                        for a in m.atomlist:
                            for exc_type in exc_list: 
                                if a.spec == exc_type: exc_flag = 1
                        if not exc_flag and r[0,axis] > qhi:
                            # Record time
                            duration_list.append(res_time)   
            # Update side flags
            if r[0,axis] < qhi:
                mol_flag_left[index] = True
            if r[0,axis] > qhi: 
                mol_flag_left[index] = False

        print duration_list 
        return mol_flag, mol_flag_init, mol_flag_left,mol_intimes, duration_list                 
    def velocity_field_grid(self, gx, gy, gz, xlo, xhi, ylo, yhi, zlo, zhi):
        """
        Generates velocity field over grid (gx, gy, gz)
        """
        print "Generating velocity field data..."
        lx = xhi - xlo
        ly = yhi - ylo
        lz = zhi - zlo

        

        v_grid = ([[[ matrix([[0.0, 0.0, 0.0]]) 
                      for k in range(0, gz)]
                      for j in range(0, gy)] 
                      for i in range(0, gx)])

        N_grid = ([[[ 0 
                      for k in range(0, gz)]
                      for j in range(0, gy)] 
                      for i in range(0, gx)])
        # Loop over atoms
        for a in self.atomlist:

            r = a.return_r_matrix()
            if ((r[0,0] > xlo and r[0,0] < xhi) and
                (r[0,1] > ylo and r[0,1] < yhi) and
                (r[0,2] > zlo and r[0,2] < zhi)):
                
                x_index = int(gx*(r[0,0] - xlo)/lx)
                y_index = int(gy*(r[0,1] - ylo)/ly)
                z_index = int(gz*(r[0,2] - zlo)/lz)
                 
        #       print v_grid      
        #       print r
        #       print x_index               
        #       print y_index               
        #       print z_index               

         #      print len(v_grid)
         #      print len(v_grid[0])
         #      print len(v_grid[0][0])

                v_grid[x_index][y_index][z_index] += a.return_v_matrix()

                N_grid[x_index][y_index][z_index] += 1

        # Average over number of atoms
        for i in range(0, gx):
            for j in range(0, gy):
                for k in range(0, gz):
                    if N_grid[i][j][k] > 0: v_grid[i][j][k] /= N_grid[i][j][k] 
        print "Done."

        return v_grid
                
    
       
    def hbond_lifetime_region(self, xlo, xhi, ylo, yhi, zlo, zhi,
                 old_hbond_list, C_hb, N_bin):
        """
        Generates hydrogen bond lifetime distribution. 
        hb_log_list is a list of hydrogen bond objects each containing a
        pair of mol_ids and a counter. 

        If the bond is no longer present it is not added to updated list.        
        A bond only exists between molecules in the same layer
        """

        Lx = self.xhi - self.xlo
        Ly = self.yhi - self.ylo
        Lz = self.zhi - self.zlo
        l_x = matrix([[Lx, 0.0, 0.0]])
        l_y = matrix([[0.0, Ly, 0.0]])
        l_z = matrix([[0.0, 0.0, Lz]])

        # Distance criteria
        pi = 4.0*math.atan(1.0)
        tol = 0.05

        R_OO_c = 3.4
        R_OO_c2 = R_OO_c*R_OO_c
        R_OH_c = 2.425
        R_OH_c2 = R_OH_c*R_OH_c

        phi_c = 30*pi/180

        ## Create new list of all molecules in layer 
        new_mollist = []
        for m in self.mollist:
            r =  self.COM_vector(m)

            if ((r[0,0] > xlo and r[0,0] < xhi) and
                (r[0,1] > ylo and r[0,1] < yhi) and
                (r[0,2] > zlo and r[0,2] < zhi)):

                new_mollist.append(m)

        new_hbond_list = []  
                 
#        ## Loop over mollist to find bonds
#        for m_i in new_mollist:
#            HB_count = 0
#            HB_mollist = []
#            HB_mollist.append(m_i)  
#            for a_i in m_i.atomlist:
#                if a_i.spec == "O":
#                    r_O_i = a_i.return_r_matrix()
#
#                    for m_j in new_mollist:
#                        cont_flag = 0
#                        r_OH_list = []
#                        r_H_list_j = []
#                        if m_i != m_j:
#                            
#                            for a_j in m_j.atomlist:
#
#                                if a_j.spec == "O":
#                                    r_O_j = a_j.return_r_matrix()
#                                    r_OO = r_O_i - r_O_j
#  
#                                    # Apply PBCs 
#                                    if r_OO[0,0] > Lx/2: r_OO -= l_x
#                                    if r_OO[0,0] < -Lx/2: r_OO += l_x
#                                    if r_OO[0,1] > Ly/2: r_OO -= l_y
#                                    if r_OO[0,1] < -Ly/2: r_OO += l_y
#                                    if r_OO[0,2] > Lz/2: r_OO -= l_z
#                                    if r_OO[0,2] < -Lz/2: r_OO += l_z
#
#                            for a_j in m_j.atomlist:
#
#                                if a_j.spec == "H":
#                                    r_H_j = a_j.return_r_matrix()
#                                    r_OH = r_O_j - r_H_j
#  
#                                    # Apply PBCs 
#                                    if r_OH[0,0] > Lx/2: r_OH -= l_x
#                                    if r_OH[0,0] < -Lx/2: r_OH += l_x
#                                    if r_OH[0,1] > Ly/2: r_OH -= l_y
#                                    if r_OH[0,1] < -Ly/2: r_OH += l_y
#                                    if r_OH[0,2] > Lz/2: r_OH -= l_z
#                                    if r_OH[0,2] < -Lz/2: r_OH += l_z
#
#                                    r_H_list_j.append(r_H_j)
#                                    r_OH_list.append(r_OH)
#
#
#                            # Check O-O distance
#                            r_OO2 = float(r_OO*r_OO.T)
#                            if r_OO2 < R_OO_c2: 
#                                r_OH2_temp = 0.0
#                                
#                                for r_H_j in r_H_list_j:
#
#                                    r_OH_ij = r_O_i - r_H_j
#                                    # Apply PBCs 
#                                    if r_OH_ij[0,0] > Lx/2: r_OH_ij -= l_x
#                                    if r_OH_ij[0,0] < -Lx/2: r_OH_ij += l_x
#                                    if r_OH_ij[0,1] > Ly/2: r_OH_ij -= l_y
#                                    if r_OH_ij[0,1] < -Ly/2: r_OH_ij += l_y
#                                    if r_OH_ij[0,2] > Lz/2: r_OH_ij -= l_z
#                                    if r_OH_ij[0,2] < -Lz/2: r_OH_ij += l_z
#                                    r_OH_ij2 = float(r_OH_ij*r_OH_ij.T)
#
#                                    if r_OH_ij2 < R_OH_c2:
#
#                                        for r_OH in r_OH_list:
#                                            r_OH2 = float(r_OH*r_OH.T)     
#                                            r_OOdotr_OH = (float(
#                                                          r_OO*r_OH.T)/(
#                                                          math.sqrt(r_OO2)
#                                                          *math.sqrt(r_OH2)))
#                               
#                                            phi = math.acos(r_OOdotr_OH)
#                                            print phi*180/pi
#                              #             print r_OH  
#                                            if phi < phi_c:
#                                            # found hbond:
#                                                print m_i.mol_id
#                                                print m_j.mol_id
#                                                print "|r_OO| = {}".format(
#                                                        math.sqrt(r_OO2))   
#                                                print "|r_OH| = {}".format(
#                                                        math.sqrt(r_OH2))   
#                                                new_hbond = Hbond(m_i.mol_id, 
#                                                      m_j.mol_id)  
#                                                new_hbond_list.append(new_hbond)
#                                                HB_count += 1.0
#                                                HB_mollist.append(m_j)
#                                                
#                                                if HB_count > 2.0:
#                                                    HB_fout = open(
#                                                    "HB_out.xyz", "w")
#                                                    HB_fout.write("{}\n"
#                                                    .format(3*len(HB_mollist)))
#                                                    HB_fout.write(
#                                      "Many hydrogen bonds!\n")
#                                                    for m in HB_mollist:
#                                                        for a in m.atomlist:
#                                                            r_a = a.return_r_matrix() 
#                                                            HB_fout.write("{} {} {} {}\n".
#                                                            format(a.spec, r_a[0,0],
#                                                               r_a[0,1], 
#                                                               r_a[0,2])) 
#                                            
                                        
#       for hbond_new in new_hbond_list:
#           print "HB: {} {}".format(hbond_new.moli_id, hbond_new.molj_id) 

        for mol in new_mollist:
            for a in mol.atomlist:
                if a.spec == 'O':
                    # Get bin index from COM
              #     print a.spec
                    r_O_a = a.return_r_matrix()
                    r_com = self.COM_vector(mol)
              #     print "r_O_a = {}, r_com_a = {}".format(r_O_a,
               #    r_com) 
                    
                    # find vectors
                              
                    for molb in new_mollist:
                        for b in molb.atomlist:
                            cont_flag = 0  
                            if str(b.spec) == 'O' and b != a:
                               # Identify potential HB candidate
                                r_O_b = b.return_r_matrix()
                                r_OO = r_O_a - r_O_b
                                if r_OO[0,0] > Lx/2: r_OO -= l_x 
                                if r_OO[0,0] < -Lx/2: r_OO += l_x 
                                if r_OO[0,1] > Ly/2: r_OO -= l_y 
                                if r_OO[0,1] < -Ly/2: r_OO += l_y 
                                if r_OO[0,2] > Lz/2: r_OO -= l_z 
                                if r_OO[0,2] < -Lz/2: r_OO += l_z 
                                 
                                # Check O-O distance
                                r_OO2 = float(r_OO*r_OO.T)
                                if r_OO2 < R_OO_c2: cont_flag = 1

                            if cont_flag: # continue
                                cont_flag = 0
                                r_OH_list_b = []
                                r_H_list_b = []
                                      
                                for j in range(0,len(molb.atomlist)):
                                    # set up r_H and r_OH vectors
                                    if molb.atomlist[j].spec == 'H':
                                        r_H = (molb.atomlist[j].
                                               return_r_matrix())
                                        r_H_list_b.append(r_H)
                                        r_OH = r_H - r_O_b 

                                        if r_OH[0,0] > Lx/2: 
                                            r_OH -= l_x 
                                        if r_OH[0,0] < -Lx/2: 
                                            r_OH += l_x 
                                        if r_OH[0,1] > Ly/2: 
                                            r_OH -= l_y 
                                        if r_OH[0,1] < -Ly/2: 
                                            r_OH += l_y 
                                        if r_OH[0,2] > Lz/2: 
                                            r_OH -= l_z 
                                        if r_OH[0,2] < -Lz/2: 
                                            r_OH += l_z 
                                        r_OH_list_b.append(r_OH)
                                
                                for j in range(0, len(r_H_list_b)):
                                    r_OH = r_O_a - r_H_list_b[j]
                                    # PBCs
                                    if r_OH[0,0] > Lx/2: r_OH -= l_x 
                                    if r_OH[0,0] < -Lx/2: r_OH += l_x 
                                    if r_OH[0,1] > Ly/2: r_OH -= l_y 
                                    if r_OH[0,1] < -Ly/2: r_OH += l_y 
                                    if r_OH[0,2] > Lz/2: r_OH -= l_z 
                                    if r_OH[0,2] < -Lz/2: r_OH += l_z 
                                    r_OH2 = r_OH*r_OH.T
                                    if r_OH2 < R_OH_c2: #proceed
             #                          print "|r_OH| = {}".format(
                 #                             math.sqrt(r_OH2))

                                        for r_OH_b in r_OH_list_b:
                                            r_OH2 = r_OH_b*r_OH_b.T
                                            r_OO2 = r_OO*r_OO.T
                                            r_OOr_OHb = r_OO*r_OH_b.T
                                                  
                                            OOdotOHb = (float(r_OOr_OHb)
                                                /(math.sqrt(float(r_OO2))
                                                *math.sqrt(float(r_OH2))))

                                            phi = math.acos(OOdotOHb)
                 #                              print ("phi = {}"
                 #                                     .format(phi*180/pi))
                                            if phi < phi_c:
                 #                              print "HB found!"
                #                               print mol.mol_id
                #                               print molb.mol_id
                #                               print "|r_OO| = {}".format(
                #                                       math.sqrt(r_OO2))   
                #                               print "|r_OH| = {}".format(
                #                                       math.sqrt(r_OH2))   
                                                new_hbond = Hbond(mol.mol_id, 
                                                      molb.mol_id)  
                                                new_hbond_list.append(new_hbond)
                                              
        for hb in old_hbond_list:
            print "{} {} {}".format(hb.moli_id, hb.molj_id, hb.del_t)

        for hbond_old in old_hbond_list:
            in_list = False
             
            for hbond_new in new_hbond_list:

                if ((hbond_new.moli_id == hbond_old.moli_id and 
                     hbond_new.molj_id == hbond_old.molj_id) or 
                    (hbond_new.moli_id == hbond_old.molj_id and 
                     hbond_new.molj_id == hbond_old.moli_id)):
                    in_list = True 
                    hbond_new.del_t = hbond_old.del_t
                    hbond_new.increment_t()
                    print "MATCH"

            if not in_list:
                print "NO MATCH"
                # Bin times: hydrogen bond no longer present
                if hbond_old.del_t < N_bin:
                    C_hb[hbond_old.del_t] += 1.0
                        
                      
        return C_hb, new_hbond_list                    

    def generate_hbond_profile(self, axis, nbin, data, print_flag):
        """
        Routine to generate hydrogen bonding profile where HB is defined 
        from the criteria given in Guardia et al. Journal of Molecular Liquids,
        117, 2005.
        cutoff radii and angles
        R_OO_c = 3.4A
        R_OH_c = 2.425A
        phi_c = 30 degrees            
        """

        R_OH = 1.0 
        tol = 0.05
        pi = 4.0*math.atan(1.0)
 
        HB_bins = [0.0]*nbin
        bins = [0.0]*nbin 
        count_bins = [0]*nbin
        err_flg = 0
  
        # Distance criteria
        R_OO_c = 3.4
        R_OO_c2 = R_OO_c*R_OO_c
        R_OH_c = 2.425
        R_OH_c2 = R_OH_c*R_OH_c
       
        phi_c = 30*pi/180
        Lx = self.xhi - self.xlo
        Ly = self.yhi - self.ylo
        Lz = self.zhi - self.zlo
            
        # check axis
        if axis == "x":
            L = self.xhi - self.xlo
        elif axis == "y":
            L = self.yhi - self.ylo
        elif axis == "z":
            L = self.zhi - self.zlo
        else:
            print "An axis must be specified. x, y, or z"
            err_flg = 1
        l_x = matrix([[Lx, 0.0, 0.0]])
        l_y = matrix([[0.0, Ly, 0.0]])
        l_z = matrix([[0.0, 0.0, Lz]])

        if not err_flg:
            print ("Generating hydrogen bonding profile for water molecules.")
            vol = ((self.xhi - self.xlo)*
                   (self.yhi - self.ylo)*
                   (self.zhi - self.zlo)/nbin)
            # Loop over molecules 
            count = 0
            mol_count = 0 
            for mol in self.mollist:
                mol_count += 1
                for a in mol.atomlist:
                    if a.spec == 'O':
                        # Get bin index from COM
                  #     print a.spec
                        r_O_a = a.return_r_matrix()
                        r_com = self.COM_vector(mol)
                  #     print "r_O_a = {}, r_com_a = {}".format(r_O_a,
                   #    r_com) 
                        if axis == "x":
                            q = r_com[0,0]
                        elif axis == "y":
                            q = r_com[0,1]
                        elif axis == "z":
                            q = r_com[0,2]
                        if q > L: q -= L  
                        if q < 0: q += L 
                        
                        
                        ibin = int(nbin*(q)/L)
                  #     print "q = {}, ibin = {}".format(q, ibin) 
                        
                        count_bins[ibin] += 1
                        # find vectors
                                  
                        for molb in self.mollist:
                            for b in molb.atomlist:
                                cont_flag = 0  
                                if str(b.spec) == 'O' and b != a:
                                   # Identify potential HB candidate
                                    r_O_b = b.return_r_matrix()
                                    r_OO = r_O_a - r_O_b
                                    if r_OO[0,0] > Lx/2: r_OO -= l_x 
                                    if r_OO[0,0] < -Lx/2: r_OO += l_x 
                                    if r_OO[0,1] > Ly/2: r_OO -= l_y 
                                    if r_OO[0,1] < -Ly/2: r_OO += l_y 
                                    if r_OO[0,2] > Lz/2: r_OO -= l_z 
                                    if r_OO[0,2] < -Lz/2: r_OO += l_z 
                                     
                                    # Check O-O distance
                                    r_OO2 = float(r_OO*r_OO.T)
                                    if r_OO2 < R_OO_c2: cont_flag = 1

                                if cont_flag: # continue
                                    cont_flag = 0
                                    r_OH_list_b = []
                                    r_H_list_b = []
                                          
                                    for j in range(0,len(molb.atomlist)):
                                        # set up r_H and r_OH vectors
                                        if molb.atomlist[j].spec == 'H':
                                            r_H = (molb.atomlist[j].
                                                   return_r_matrix())
                                            r_H_list_b.append(r_H)
                                            r_OH = r_H - r_O_b 

                                            if r_OH[0,0] > Lx/2: 
                                                r_OH -= l_x 
                                            if r_OH[0,0] < -Lx/2: 
                                                r_OH += l_x 
                                            if r_OH[0,1] > Ly/2: 
                                                r_OH -= l_y 
                                            if r_OH[0,1] < -Ly/2: 
                                                r_OH += l_y 
                                            if r_OH[0,2] > Lz/2: 
                                                r_OH -= l_z 
                                            if r_OH[0,2] < -Lz/2: 
                                                r_OH += l_z 
                                            r_OH_list_b.append(r_OH)
                                    
                                    for j in range(0, len(r_H_list_b)):
                                        r_OH = r_O_a - r_H_list_b[j]
                                        # PBCs
                                        if r_OH[0,0] > Lx/2: r_OH -= l_x 
                                        if r_OH[0,0] < -Lx/2: r_OH += l_x 
                                        if r_OH[0,1] > Ly/2: r_OH -= l_y 
                                        if r_OH[0,1] < -Ly/2: r_OH += l_y 
                                        if r_OH[0,2] > Lz/2: r_OH -= l_z 
                                        if r_OH[0,2] < -Lz/2: r_OH += l_z 
                                        r_OH2 = r_OH*r_OH.T
                                        if r_OH2 < R_OH_c2: #proceed
                 #                          print "|r_OH| = {}".format(
                     #                             math.sqrt(r_OH2))

                                            for r_OH_b in r_OH_list_b:
                                                r_OH2 = r_OH_b*r_OH_b.T
                                                r_OO2 = r_OO*r_OO.T
                                                r_OOr_OHb = r_OO*r_OH_b.T
                                                      
                                                OOdotOHb = (float(r_OOr_OHb)
                                                    /(math.sqrt(float(r_OO2))
                                                    *math.sqrt(float(r_OH2))))

                                                phi = math.acos(OOdotOHb)
                     #                              print ("phi = {}"
                     #                                     .format(phi*180/pi))
                                                if phi < phi_c:
                     #                              print "HB found!"
                                                    r_com_b = (
                                                    self.COM_vector(molb))
 
                                                    if axis == "x":
                                                        q_b = r_com_b[0,0]
                                                    elif axis == "y":
                                                        q_b = r_com_b[0,1]
                                                    elif axis == "z":
                                                        q_b = r_com_b[0,2]

                                                    if q_b > L: q_b -= L  
                                                    if q_b < 0: q_b += L 
                                                    ibin_b = int(
                                                            nbin*(q_b)/L)
                                                  
                                                    HB_bins[ibin] += 1.0       
                                                    HB_bins[ibin_b] += 1.0
       
        for i in range(0, len(HB_bins)):
            # HB per molecule
            if count_bins[i] == 0: HB_bins[i] = 0.0 
            else:
                HB_bins[i] /= count_bins[i]  
#               print HB_bins[i]
#               print count_bins[i]
        header = "{} (A)      HB per molecule   Err".format(axis)
        axis_label = nbin*[0.0]
        for i in range(0, nbin):
                
            axis_label[i] = (i+0.5)/nbin*L

        if data: 
            data.add_data_set(header, axis_label, HB_bins)  
                                
    
    

    def generate_density_profile(self, axis, nbin, mol_sig, data, print_flag):
        """ 
        Routine to generate "molecular" number density profile from 
        current configuration.
        axis is specified as string "x", "y", or "z"
        nbin is number of bins for histogram
        mol_sig is the first part of the string that prefixes the atom type
        associated with the atom number. e.g. "MPD" in "MPD_N".
        The density profile will only be created for one particular 
        molecule type. This may be altered later.  
        """
        
        bins = [0.0]*nbin
        err_flg = 0
        count = 0

        # check axis
        if axis == "x":
            L = self.xhi - self.xlo
        elif axis == "y":
            L = self.yhi - self.ylo
        elif axis == "z":
            L = self.zhi - self.zlo
        else:
            print "An axis must be specified. x, y, or z"
            err_flg = 1
             
        if not err_flg:
            print ("Generating density profile for {} molecules."
                   .format(mol_sig))
            vol = ((self.xhi - self.xlo)*
                   (self.yhi - self.ylo)*
                   (self.zhi - self.zlo)/nbin)
            # Loop over molecules 
            count = 0
            for mol in self.mollist:
                # Check first atom in molecule for mol_sig
                if mol.atomlist[0].spec.find(mol_sig) >= 0:
                    # continue with stuff

                    R = self.COM_vector(mol)

                    if axis == "x":
                        q = R[0,0]
                    elif axis == "y":
                        q = R[0,1]
                    elif axis == "z":
                        q = R[0,2]

                    ibin = int(nbin*q/L)

                    bins[ibin] += 1 # increase bin by one                  
                    count += 1                   

            if print_flag:
                f_dens = open(mol_sig+"_density_prof_t_"
                              +str(self.timestep)+".dat", "w")
                f_dens.write("{} (A)      freq.\n".format(axis))  

            header = "{} (A)      freq.        Err".format(axis)


            axis_label = nbin*[0.0]
            for i in range(0, nbin):
                
                bins[i] = bins[i]/vol
                axis_label[i] = (i+0.5)/nbin*L 
                if print_flag: 
                    f_dens.write('{0:<7.3}    {1:<6.4}\n'.format(
                                 (i+0.5)/nbin*L, bins[i]))
            
        if err_flg:
            pass

        if data: 
            data.add_data_set(header, axis_label, bins)  
             
    def generate_density_profile_region(self, axis, nbin, mol_sig, 
                                        qlo, qhi, data, print_flag):
        """ 
        Routine to generate "molecular" number density profile from 
        current configuration.
        axis is specified as string "x", "y", or "z"
        nbin is number of bins for histogram
        mol_sig is the first part of the string that prefixes the atom type
        associated with the atom number. e.g. "MPD" in "MPD_N".
        The density profile will only be created for one particular 
        molecule type. This may be altered later. 
        start at qlo end at qhi 
        """
        
        bins = [0.0]*nbin
        err_flg = 0
        count = 0

        L_x = self.xhi - self.xlo
        L_y = self.yhi - self.ylo
        L_z = self.zhi - self.zlo
        # check axis
        if axis == "x":
            vol = (qhi-qlo)*L_y*L_z/nbin
        elif axis == "y":
            vol = L_x*(qhi-qlo)*L_z/nbin
        elif axis == "z":
            vol = L_x*L_y*(qhi-qlo)/nbin
        else:
            print "An axis must be specified. x, y, or z"
            err_flg = 1
             
        if not err_flg:
            print ("Generating density profile for {} molecules."
                   .format(mol_sig))
            
            # Loop over molecules 
            count = 0
            for mol in self.mollist:
                # Check first atom in molecule for mol_sig
                if mol.atomlist[0].spec.find(mol_sig) >= 0:
                    # continue with stuff

                    R = self.COM_vector(mol)

                    if axis == "x":
                        q = R[0,0]
                    elif axis == "y":
                        q = R[0,1]
                    elif axis == "z":
                        q = R[0,2]

                    ibin = int(nbin*q/L)

                    bins[ibin] += 1 # increase bin by one                  
                    count += 1                   

            if print_flag:
                f_dens = open(mol_sig+"_density_prof_t_"
                              +str(self.timestep)+".dat", "w")
                f_dens.write("{} (A)      freq.\n".format(axis))  

            header = "{} (A)      freq.        Err".format(axis)


            axis_label = nbin*[0.0]
            for i in range(0, nbin):
                
                bins[i] = bins[i]/vol
                axis_label[i] = (i+0.5)/nbin*L 
                if print_flag: 
                    f_dens.write('{0:<7.3}    {1:<6.4}\n'.format(
                                 (i+0.5)/nbin*L, bins[i]))
            
        if err_flg:
            pass

        if data: 
            data.add_data_set(header, axis_label, bins)  
