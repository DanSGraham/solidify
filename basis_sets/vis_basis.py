#!/usr/bin/env python3

import sys, argparse, re
import numpy as np
import matplotlib.pyplot as plt

def get_exp_and_coeff(atom_name, basis_str):
    pattern = re.compile(atom_name + "\s*\D\s*\n([^#]*\n)*")
    basis_str = pattern.search(basis_str).group(0)
    split_pattern = re.compile(atom_name + "\s*[A-Z]+ *")
    ao_functions = split_pattern.split(basis_str)[1:]
    labels = split_pattern.findall(basis_str)
    basis_list = []
    s_val = 1
    p_val = 2
    d_val = 3
    f_val = 4
    for i in range(len(ao_functions)):
        element, ao_name = labels[i].split()
        gauss_primitives = ao_functions[i].split("\n")[1:-1]
        exponent_list = [x.split()[0] for x in gauss_primitives]
        coeff_list = []
        for j in range(len(ao_name)):
            if ao_name[j] == 'S':
                ao_label = str(s_val) + ao_name[j].lower()
                s_val += 1
            if ao_name[j] == 'P':
                ao_label = str(p_val) + ao_name[j].lower()
                p_val += 1
            if ao_name[j] == 'D':
                ao_label = str(d_val) + ao_name[j].lower()
                d_val += 1
            if ao_name[j] == 'F':
                ao_label = str(f_val) + ao_name[j].lower()
                f_val += 1
                
            coeff_list = [x.split()[j+1] for x in gauss_primitives]
            basis_list.append((element, ao_label, exponent_list, coeff_list))

    return basis_list

def plot_1d(basis_type, view_axis, element, ao_name, exponents, coeff):
    x = np.linspace(-5,5,400)
    total_ao = np.zeros_like(x)
    threshold = 0.001
    if ao_name[1] == 's':
        fig, ax = plt.subplots(1,1)
        for i in range(len(exponents)):
            x = np.linspace(-5,5,400)
            y = (((2 * float(exponents[i]))/np.pi) ** (3./4.) * np.exp(-float(exponents[i]) * (x ** 2.)))*float(coeff[i])
            total_ao += y
            filter_list = np.array([np.abs(val) > threshold for val in y])
            ax.plot(x[filter_list], y[filter_list], label='Gaussian ' + str(i))
        x = np.linspace(-5,5,400)
        filter_list = np.array([np.abs(val) > threshold for val in total_ao])
        ax.plot(x[filter_list], total_ao[filter_list], label='Total Orbital')
        plt.legend()
        plt.title(element + " " + ao_name + " " + basis_type)
        plt.ylabel('Orbital amplitude')
        plt.xlabel(view_axis)
        ax.yaxis.set_label_coords(0,0.5)
        ax.xaxis.set_label_coords(1,0.54)
        # set the x-spine
        ax.spines['left'].set_position('zero')
        
        # turn off the right spine/ticks
        ax.spines['right'].set_color('none')
        ax.yaxis.tick_left()
        
        # set the y-spine
        ax.spines['bottom'].set_position('zero')
        
        # turn off the top spine/ticks
        ax.spines['top'].set_color('none')
        ax.xaxis.tick_bottom()

        y_tick_labels = ax.get_yticklabels()
        for tick in y_tick_labels:
            tick.set_visible(False)
        plt.show()
    elif ao_name[1] == 'p':
        fig, ax = plt.subplots(1,1)
        for i in range(len(exponents)):
            x = np.linspace(-5,5,400)
            y = (((128 * float(exponents[i])**5.)/np.pi**3.) ** (1./4.) * x * np.exp(-float(exponents[i]) * (x ** 2.)))*float(coeff[i])
            total_ao += y
            filter_list = np.array([np.abs(val) > threshold for val in y])
            ax.plot(x[filter_list], y[filter_list], label='Gaussian ' + str(i))
        x = np.linspace(-5,5,400)
        filter_list = np.array([np.abs(val) > threshold for val in total_ao])
        ax.plot(x[filter_list], total_ao[filter_list], label='Total Orbital')
        plt.legend()
        plt.title(element + " " + ao_name + view_axis + " " + basis_type)
        plt.ylabel('Orbital amplitude')
        plt.xlabel(view_axis)
        ax.yaxis.set_label_coords(0,0.5)
        ax.xaxis.set_label_coords(1,0.54)
        # set the x-spine
        ax.spines['left'].set_position('zero')
        
        # turn off the right spine/ticks
        ax.spines['right'].set_color('none')
        ax.yaxis.tick_left()
        
        # set the y-spine
        ax.spines['bottom'].set_position('zero')
        
        # turn off the top spine/ticks
        ax.spines['top'].set_color('none')
        ax.xaxis.tick_bottom()

        y_tick_labels = ax.get_yticklabels()
        for tick in y_tick_labels:
            tick.set_visible(False)
        plt.show()
    elif ao_name[1] == 'd':
        pass
    elif ao_name[1] == 'f':
        pass
    else:
        print ("AO UNKNOWN")

def limit_scope(threshold, x, y, z_2d):
    x_start = 0
    y_start = 0
    x_rad = False
    y_rad = False
    while x_start < len(x) and not x_rad:
        for val in z_2d[x_start]:
            if val >= threshold:
                x_rad = True
        x_start += 1

    while y_start < len(y) and not y_rad:
        for val in z_2d.T[y_start]:
            if val >= threshold:
                y_rad = True
        y_start += 1

    return x[x_start:-x_start], y[y_start:-y_start], z_2d[x_start:-x_start,y_start:-y_start]
            
    
    
def plot_2d(basis_type, view_axis, element, ao_name, exponents, coeff):
    bounds = 10
    density = 700
    x = np.linspace(-bounds,bounds,density)
    y = np.linspace(-bounds,bounds,density)
    total_ao = np.zeros((len(x), len(y)))
    threshold = 0.01
    if ao_name[1] == 's':
        fig, ax = plt.subplots(1,1)
        for i in range(len(exponents)):
            x = np.linspace(-bounds,bounds,density)
            y = np.linspace(-bounds,bounds,density)
            for m in range(len(x)):
                for n in range(len(y)):
                    z = (((2 * float(exponents[i]))/np.pi) ** (3./4.) * np.exp(-float(exponents[i]) * (x[m] ** 2.))) * np.exp(-float(exponents[i]) * (y[n] ** 2.)) * float(coeff[i])
                    total_ao[m][n] += z
        x = np.linspace(-bounds,bounds,density)
        y = np.linspace(-bounds,bounds,density)
        x, y, total_ao = limit_scope(threshold, x, y, total_ao)
        plt.contourf(x, y, total_ao, label='Total Orbital')
        plt.legend()
        plt.title(element + " " + ao_name + " " + basis_type)
        plt.ylabel('Orbital amplitude')
        plt.xlabel(view_axis)
        #ax.yaxis.set_label_coords(0,0.5)
        #ax.xaxis.set_label_coords(1,0.54)
        # set the x-spine
        #ax.spines['left'].set_position('zero')
        
        # turn off the right spine/ticks
        #ax.spines['right'].set_color('none')
        #ax.yaxis.tick_left()
        
        # set the y-spine
        #ax.spines['bottom'].set_position('zero')
        
        # turn off the top spine/ticks
        #ax.spines['top'].set_color('none')
        #ax.xaxis.tick_bottom()

        #y_tick_labels = ax.get_yticklabels()
        #for tick in y_tick_labels:
        #    tick.set_visible(False)
        plt.show()
    elif ao_name[1] == 'p':
        fig, ax = plt.subplots(1,1)
        for i in range(len(exponents)):
            x = np.linspace(-5,5,400)
            y = (((128 * float(exponents[i])**5.)/np.pi**3.) ** (1./4.) * x * np.exp(-float(exponents[i]) * (x ** 2.)))*float(coeff[i])
            total_ao += y
            filter_list = np.array([np.abs(val) > threshold for val in y])
            ax.plot(x[filter_list], y[filter_list], label='Gaussian ' + str(i))
        x = np.linspace(-5,5,400)
        filter_list = np.array([np.abs(val) > threshold for val in total_ao])
        ax.contourf(x, y, total_ao, label='Total Orbital')
        plt.legend()
        plt.title(element + " " + ao_name + view_axis + " " + basis_type)
        plt.ylabel('Orbital amplitude')
        plt.xlabel(view_axis)
        ax.yaxis.set_label_coords(0,0.5)
        ax.xaxis.set_label_coords(1,0.54)
        # set the x-spine
        ax.spines['left'].set_position('zero')
        
        # turn off the right spine/ticks
        ax.spines['right'].set_color('none')
        ax.yaxis.tick_left()
        
        # set the y-spine
        ax.spines['bottom'].set_position('zero')
        
        # turn off the top spine/ticks
        ax.spines['top'].set_color('none')
        ax.xaxis.tick_bottom()

        y_tick_labels = ax.get_yticklabels()
        for tick in y_tick_labels:
            tick.set_visible(False)
        plt.show()
    elif ao_name[1] == 'd':
        pass
    elif ao_name[1] == 'f':
        pass
    else:
        print ("AO UNKNOWN")

def plot_3d(basis_type, view_axis, element, ao_name, exponents, coeff):
    pass
 


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
    pass

def plot_3d(basis_type, view_axis, element, ao_name, exponents, coeff):
    pass
 


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
    plot_ao = None
    for ao in basis_format:
        if args.ao_type.lower() == ao[1]:
            plot_ao = ao

    if args.disp_dmn == '1d':
        plot_1d(args.basis_name, args.axes_persp, plot_ao[0], plot_ao[1], plot_ao[2], plot_ao[3])
    elif args.disp_dmn == '2d':
        plot_2d(args.basis_name, args.axes_persp, plot_ao[0], plot_ao[1], plot_ao[2], plot_ao[3])
    elif args.disp_dmn == '3d':
        plot_3d(args.basis_name, args.axes_persp, plot_ao[0], plot_ao[1], plot_ao[2], plot_ao[3])
    else:
        print ("AXIS NOT GIVEN")

if __name__=="__main__":
    main(sys.argv[1:])
