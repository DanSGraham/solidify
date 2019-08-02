#A module to visualize different basis sets.
import numpy as np
import re
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

#Gaussian basis sets.


#STO-nG
def get_exp_and_coeff(atom_name, basis_str):
    pattern = re.compile(atom_name + "\s*\D\s*\n([^#]*\n)*")
    basis_str = pattern.search(basis_str).group(0)
    split_pattern = re.compile(atom_name + "\s*[A-Z]+ *")
    ao_functions = split_pattern.split(basis_str)[1:]
    labels = split_pattern.findall(basis_str)
    basis_list = []
    for i in range(len(ao_functions)):
        element, ao_name = labels[i].split()
        gauss_primitives = ao_functions[i].split("\n")[1:-1]
        exponent_list = [x.split()[0] for x in gauss_primitives]
        coeff_list = []
        for j in range(len(ao_name)):
            ao_label = str(i+1) + ao_name[j]        
            coeff_list = [x.split()[j+1] for x in gauss_primitives]
            basis_list.append((element, ao_label, exponent_list, coeff_list))

    return basis_list
    
def plot_2d(basis_type, element, ao_name, exponents, coeff):
    plt.title(element + " " + ao_name)
    x = np.linspace(0,5,200)
    total_ao = np.zeros_like(x)
    if basis_type == 'STO-nG':
        if ao_name[1] == 'S':
            for i in range(len(exponents)):
                y = (((2 * float(exponents[i]))/np.pi) ** (3./4.) * np.exp(-float(exponents[i]) * (x ** 2.)))*float(coeff[i])
                plt.plot(x, y)
                total_ao += y
        elif ao_name[1] == 'P':
            for i in range(len(exponents)):
                y = (((128 * float(exponents[i])**5.)/np.pi**3.) ** (1./4.) * x * np.exp(-float(exponents[i]) * (x ** 2.)))*float(coeff[i])
                plt.plot(x, y)
                total_ao += y
        elif ao_name[1] == 'D':
            pass
        else:
            print ("UNKNOWN AO NAME")
    elif basis_type == 'X-YZg':
        pass
    elif basis_type == 'cc':
        pass
    plt.plot(x, total_ao)
    plt.show()
 
basis_str = ''
with open('sto_2g.dat', 'r') as fin:
 basis_str = fin.read()

basis_format = get_exp_and_coeff('C', basis_str)
plot_2d('STO-nG', basis_format[0][0], basis_format[0][1], basis_format[0][2], basis_format[0][3])
plot_2d('STO-nG', basis_format[1][0], basis_format[1][1], basis_format[1][2], basis_format[1][3])
plot_2d('STO-nG', basis_format[2][0], basis_format[2][1], basis_format[2][2], basis_format[2][3])
#for basis in basis_format:
#    x = np.linspace(0,5,200)
#    y1 = (((2 * float(basis[2][0]))/np.pi) ** (3./4.) * np.exp(-float(basis[2][0]) * x ** 2.))*float(basis[3][0])
#    y2 = (((2 * float(basis[2][1]))/np.pi) ** (3./4.) * np.exp(-float(basis[2][1]) * x ** 2.))*float(basis[3][1])
#    y3 = (((2 * float(basis[2][2]))/np.pi) ** (3./4.) * np.exp(-float(basis[2][2]) * x ** 2.))*float(basis[3][2])
#    plt.title(basis[0] + " " + basis[1])
#    plt.plot(x, y1)
#    plt.plot(x, y2)
#    plt.plot(x, y3)
#    plt.plot(x, y1 + y2 + y3)
#    plt.show()
