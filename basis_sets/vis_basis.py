#!/usr/bin/env python3

import sys, argparse, re

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


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("basis_name", help='specify name of basis')
    parser.add_argument("element", help='specify element name')
    parser.add_argument("ao_type", help='specify ao type')
    parser.add_argument("disp_dmn", help='specify display dimensions')
    parser.add_argument("axes_persp", help='specify axis of perspective')
    args = parser.parse_args()
    basis_fname = args.basis_name.replace("-", "_")
    with open(basis_fname + ".dat", 'r') as fin:
        basis_str = fin.read()

    basis_format = get_exp_and_coeff(args.element, basis_str)

if __name__=="__main__":
    main(sys.argv[1:])
