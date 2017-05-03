## Selection of operations to be performed on vectors - mostly uses numpy 
## matrix class
import numpy
from numpy import matrix
import math

def get_normal(r_ij, r_kl):
    """
    Returns unit vector normal to vector r_ij and r_kl
    """
    n = numpy.cross(r_ij, r_kl)
    n = n/math.sqrt(float(matrix.dot(n, n.T)))

    return n

def get_basis(r_ij, r_kl):
    """
    Returns a list of vectors corresponding to an orthonormal basis
    given two orthogonal vectors r_ij and r_kl
    """

    e = []

    e.append(r_ij/math.sqrt(float(matrix.dot(r_ij, r_ij.T))))
    e.append(r_kl/math.sqrt(float(matrix.dot(r_kl, r_kl.T))))
    e.append(get_normal(r_ij, r_kl))

    return e

def transform_translate(e_prime, r_trans, point_set):
                                  
    # calculate transformation matrix elements L_ij
    e = []
    e.append(matrix([[1.0, 0.0, 0.0]]))
    e.append(matrix([[0.0, 1.0, 0.0]]))
    e.append(matrix([[0.0, 0.0, 1.0]]))
        
    L = numpy.zeros(shape=(3, 3))
        
    for i in range(0, len(e_prime)):
        for j in range(0, len(e)):
            L_ij = float(matrix.dot(e_prime[i], e[j].T))
            L[i][j] = L_ij


    new_set = []

    for r in point_set:
    
        # Rotate about the origin
        r = r*L

        # Translate
        r += r_trans
      
        new_set.append(r)

    return new_set

