# Atom class for use when reading lammps dump files
# Contains attributes id, type, molid, position, velocity, mass
from numpy import matrix

class Atom:
    def __init__(self, i, moli, species, x, y, z, xv, yv, zv, m):
        self.atom_id = int(i)
        self.mol_id = int(moli)
        self.spec = species
        self.rx = float(x)
        self.ry = float(y)
        self.rz = float(z)
        self.vx = float(xv)
        self.vy = float(yv)
        self.vz = float(zv)
        self.mass = float(m)
        self.q = 0.0

    def return_r_matrix(self):
        r = matrix([[self.rx, self.ry, self.rz]])
        return r
    
    def return_r_array(self):
        r = numpy.array([self.rx, self.ry, self.rz])
        return r

    def return_v_matrix(self):
        v = matrix([[self.vx, self.vy, self.vz]])
        return v

    def print_atom(self):
            print "{0} {1:6f} {2:6f} {3:6f}".format(
                  self.spec, self.rx, self.ry, self.rz)

    def set_r(self, r):
        """
        Assigns position vector to atom given as numpy matrix.
        """
        self.rx = r[0,0] 
        self.ry = r[0,1] 
        self.rz = r[0,2] 

    def set_v(self, v):
        """
        Assigns velocity vector to atom given as numpy matrix.
        """
        self.vx = v[0,0] 
        self.vy = v[0,1] 
        self.vz = v[0,2] 

    def set_spec(self, s):
        """
        Assigns species to given atom.
        """
        self.spec = s

    def set_charge(self, charge):
        """
        Assigns partial charge value to atom
        """
        self.q = charge

    def set_atom_id(self, aid):
        self.atom_id = aid
 
    def set_mol_id(self, aid):
        self.atom_id = aid

class Molecule:

    def __init__(self, moli):
        self.mol_id = moli
        self.atomlist = []
        self.mol_type = ''

    def add_atom(self, atom):
        self.atomlist.append(atom)       
    
    def set_type(self, label):
        self.mol_type = label
