#A module to visualize different basis sets.
import numpy as np
import re

#Gaussian basis sets.


#STO-3g

def get_exp_and_coeff(atom_name, basis_str):
    pattern = re.compile("\n" + atom_name + "")

basis_str = ''
with open('sto_3g.dat', 'r') as fin:
 basis_str = fin.read()

x = np.linspace(0,5,200)
y1 = ((2 * 0.1543289673)/np.pi) ** (3./4.) * np.exp(-3.425250914 * x ** 2.)
print (x)
print (y1)
