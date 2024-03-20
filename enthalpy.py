#!/work/apps/python3/anaconda/2023.1/bin/python
#---Change/remove the above line depending on your environment---


import os
import numpy as np
import subprocess
import math
from pymatgen.io.vasp import Poscar
from pymatgen.analysis.local_env import BrunnerNN_real


def coordination_number_of_max_atomic_number(poscar_path):
    # Read the POSCAR file
    poscar = Poscar.from_file(poscar_path)
    structure = poscar.structure

    # Determine the coordination number using Brunner's algorithm
    nn = BrunnerNN_real(cutoff=10.0)  # Set a different cutoff distance

    # Find the element with the highest atomic number
    max_atomic_number = max(site.specie.number for site in structure)

    # Loop over all atoms in the structure
    for i, site in enumerate(structure):
        # Check if the atom is of the element with the highest atomic number
        if site.specie.number == max_atomic_number:
            info = nn.get_nn_info(structure, i)
            coordination_number = len(info)
            #print(f"The coordination number of {site.species_string} atom at site {i} is {coordination_number}.")
    
    return coordination_number


def grep_tail_awk(file_path, pattern):
    last_match = None
    with open(file_path, 'r') as file:
        for line in file:
            if pattern in line:
                last_match = line.strip().split()[4]
    return last_match

def read_poscar(file_path):
    num_elements = []
    with open(file_path, 'r') as file:
        for line_num, line in enumerate(file, start=1):
            if line_num == 7:
                num_elements = line.strip().split()
    return num_elements

def get_structural(file_path):
    with open(file_path, "r") as f:
        lines = f.readlines()

    # Store lattice vectors
    lattice_vectors = np.zeros((3, 3))
    for i in range(3):
        lattice_vectors[i, :] = list(map(float, lines[2+i].strip().split()))

    # Calculate lattice parameters (a, b, c)
    a = np.sqrt(np.sum(lattice_vectors[0, :] ** 2))
    b = np.sqrt(np.sum(lattice_vectors[1, :] ** 2))
    c = np.sqrt(np.sum(lattice_vectors[2, :] ** 2))

    # Calculate alpha, beta, gamma
    alpha = np.arccos(np.dot(lattice_vectors[1, :], lattice_vectors[2, :]) / (b * c)) * (180 / np.pi)
    beta  = np.arccos(np.dot(lattice_vectors[0, :], lattice_vectors[2, :]) / (a * c)) * (180 / np.pi)
    gamma = np.arccos(np.dot(lattice_vectors[0, :], lattice_vectors[1, :]) / (a * b)) * (180 / np.pi)

    # Calculate volume
    volume = np.dot(lattice_vectors[0, :], np.cross(lattice_vectors[1, :], lattice_vectors[2, :]))

    return lattice_vectors, a, b, c, alpha, beta, gamma, volume

def read_label(file_path):
    element_labels = []
    with open(file_path, 'r') as file:
        for line in file:
            if 'PBE' in line and 'POTCAR' not in line and 'TITEL' not in line:
                fields = line.strip().split()
                if len(fields) >= 2:
                    element_labels.append(fields[1])
    return element_labels

def calculate_compound_energy(poscar_path, outcar_path, element_values):
    compound_energy = 0.0
    num_elements = read_poscar(poscar_path)
    element_labels = read_label(outcar_path)

    for element, quantity in zip(element_labels, num_elements):
        if element in element_values:
            energy = element_values[element]
            quantity = int(quantity)
            compound_energy += energy * quantity
    return compound_energy

def formation_enthalpy(outcar_path, compound_energy, num_elements):
    pattern = 'free  energy'
    total_energy = float(grep_tail_awk(outcar_path, pattern))
    total_num_elements = sum(int(num) for num in num_elements)
    formation_enthalpy_value = (total_energy - compound_energy) / total_num_elements
    return "{:.3f}".format(formation_enthalpy_value)

def sanitize_element(element):
    """Returns the element symbol without suffixes."""
    return element.split('_')[0]

element_values = {
    'H': -3.434, 'He': -0.004, 'Li_sv': -1.731, 'Be': -3.653, 'B': -6.656,
    'C': -9.044, 'N': -8.195, 'O': -4.523, 'F': -1.443, 'Ne': -0.029,
    'Na_pv': -1.196, 'Mg': -1.417, 'Al': -3.66, 'Si': -5.386, 'P': -5.175,
    'S': -3.868, 'Cl': -1.479, 'Ar': -0.006, 'K_sv': -0.987, 'Ca_pv': -1.78,
    'Sc_sv': -6.344, 'Ti': -7.702, 'V': -8.898, 'Cr': -9.463, 'Mn': -8.898,
    'Fe': -8.499, 'Co': -7.078, 'Ni': -5.587, 'Cu': -3.71, 'Zn': -1.157,
    'Ga_d': -2.902, 'Ge_d': -4.522, 'As': -4.593, 'Se': -3.374, 'Br': -1.333,
    'Kr': -0.004, 'Rb_sv': -0.881, 'Sr_sv': -1.549, 'Y_sv': -6.449,
    'Zr_sv': -8.438, 'Nb_sv': -10.017, 'Mo_pv': -10.921, 'Tc_pv': -10.457,
    'Ru': -9.21, 'Rh': -7.319, 'Pd': -5.197, 'Ag': -2.907, 'Cd': -0.861,
    'In_d': -2.609, 'Sn_d': -3.938, 'Sb': -4.155, 'Te': -3.027, 'I': -1.365,
    'Xe': -0.64, 'Cs_sv': -0.744, 'Ba_sv': -1.479, 'La': -4.959, 'Ce_3': -4.564,
    'Pr_3': -4.627, 'Nd_3': -4.697, 'Pm_3': -4.716, 'Sm_3': -4.606,
    'Eu_2': -1.708, 'Eu_3': -4.4758596, 'Gd_3': -4.718, 'Tb_3': -4.731,
    'Dy_3': -4.66, 'Ho_3': -4.573, 'Er_3': -4.582, 'Tm_3': -4.451,
    'Yb_2': -1.125, 'Yb_3': -4.43728973, 'Lu_3': -4.549, 'Hf_pv': -9.902,
    'Ta_pv': -11.941, 'W_pv': -13.13, 'Re': -12.378, 'Os_pv': -11.374,
    'Ir': -8.953, 'Pt': -6.162, 'Au': -3.283, 'Hg': -0.374, 'Tl_d': -2.48, 
    'Pb_d': -3.951, 'Bi_d': -4.199, 'Ac': -4.106, 'Th': -7.237, 'Pa': -9.497, 
    'U': -11.032, 'Np': -12.797,'Pu': -13.95, 'Am': -15.50779025, 'Cm': -17.23058243
}

element_electronegativity = {
    'H': 2.20, 'He': 4.00, 'Li': 0.98, 'Be': 1.57, 'B': 2.04, 'C': 2.55, 'N': 3.04, 'O': 3.44,
    'F': 3.98, 'Ne': 0.00, 'Na': 0.93, 'Mg': 1.31, 'Al': 1.61, 'Si': 1.90, 'P': 2.19, 'S': 2.58,
    'Cl': 3.16, 'Ar': 0.00, 'K': 0.82, 'Ca': 1.00, 'Sc': 1.36, 'Ti': 1.54, 'V': 1.63, 'Cr': 1.66,
    'Mn': 1.55, 'Fe': 1.83, 'Co': 1.88, 'Ni': 1.91, 'Cu': 1.90, 'Zn': 1.65, 'Ga': 1.81, 'Ge': 2.01,
    'As': 2.18, 'Se': 2.55, 'Br': 2.96, 'Kr': 0.00, 'Rb': 0.82, 'Sr': 0.95, 'Y': 1.22, 'Zr': 1.33,
    'Nb': 1.60, 'Mo': 2.16, 'Tc': 1.90, 'Ru': 2.20, 'Rh': 2.28, 'Pd': 2.20, 'Ag': 1.93, 'Cd': 1.69,
    'In': 1.78, 'Sn': 1.96, 'Sb': 2.05, 'Te': 2.10, 'I': 2.66, 'Xe': 2.60, 'Cs': 0.79, 'Ba': 0.89,
    'La': 1.10, 'Ce': 1.12, 'Pr': 1.13, 'Nd': 1.14, 'Pm': 1.13, 'Sm': 1.17, 'Eu': 1.20, 'Gd': 1.20,
    'Tb': 1.10, 'Dy': 1.22, 'Ho': 1.23, 'Er': 1.24, 'Tm': 1.25, 'Yb': 1.10, 'Lu': 1.27, 'Hf': 1.30,
    'Ta': 1.50, 'W': 2.36, 'Re': 1.90, 'Os': 2.20, 'Ir': 2.20, 'Pt': 2.28, 'Au': 2.54, 'Hg': 2.00,
    'Tl': 1.62, 'Pb': 2.33, 'Bi': 2.02, 'Po': 2.00, 'At': 2.20, 'Rn': 0.00, 'Fr': 0.70, 'Ra': 0.90,
    'Ac': 1.10, 'Th': 1.30, 'Pa': 1.50, 'U': 1.38, 'Np': 1.36, 'Pu': 1.28, 'Am': 1.30, 'Cm': 1.30
}

atomic_number = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
    'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19,
    'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28,
    'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37,
    'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46,
    'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55,
    'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64,
    'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73,
    'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82,
    'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91,
    'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96
}


# Check if the file doesn't exist or if it exists but is empty or has less than 5 lines
if (not os.path.exists('CONTCAR') or
    (os.path.exists('CONTCAR') and os.stat('CONTCAR').st_size == 0) or
    (os.path.exists('CONTCAR') and len(open('CONTCAR').readlines()) < 5)):

    # Write the first two lines
    with open('CONTCAR', 'w') as f:
        f.write('first_line\n')
        f.write('second_line\n')

    # Execute the grep and awk command
    grep_cmd = "grep -A3 'direct lattice vectors' OUTCAR | tail -3 | awk '{print $1,$2,$3}'"
    result = subprocess.check_output(grep_cmd, shell=True).decode()

    # Append the result to CONTCAR
    with open('CONTCAR', 'a') as f:
        f.write(result)
 
outcar_path = 'OUTCAR'
poscar_path = 'POSCAR'
contcar_path = 'CONTCAR'
num_elements = read_poscar(poscar_path)

#---to find the greatest common divisor----#
temp_num = [int(item) for item in num_elements]
gcd = temp_num[0] # to initialize 
temp_sum = 0
for item in temp_num:
	temp_sum = temp_sum + item
#print(temp_sum)
for item in temp_num[1:]:
	gcd = math.gcd(gcd,item)
#print(gcd)

#print(num_elements)

element_labels = read_label(outcar_path)
#print(element_labels)

compound_energy = calculate_compound_energy(poscar_path, outcar_path, element_values)
#print("Compound Energy:", compound_energy)

formation_enthalpy_value = formation_enthalpy(outcar_path, compound_energy, num_elements)
#formation_enthalpy_value_kj_mol = ((formation_enthalpy(outcar_path, compound_energy, num_elements)*temp_sum)/gcd)*96.487
#print("Formation Enthalpy:", formation_enthalpy_value)

lattice_vectors, a, b, c, alpha, beta, gamma, volume = get_structural(contcar_path)

#print(element_labels)
#print(num_elements)

denom = 0
for item in num_elements:
	denom = denom + int(item)
#print(denom)

combined_array = [[a,b] for a,b in zip(element_labels,num_elements)]
#print(combined_array)

numerator = 0
for item in combined_array:
	sanitized_item = sanitize_element(item[0])
	numerator = int(item[1])*element_electronegativity[sanitized_item] + numerator

weighted_elec_neg= numerator/denom


coordination_number = coordination_number_of_max_atomic_number(poscar_path)

#for item in element_labels:
#	print(item+' ', end='')
#for item in element_labels:
#    sanitized_item = sanitize_element(item)
#    print(atomic_number[sanitized_item],' ', end='')
#for item in num_elements:
#	print(item+' ', end='')
#for item in element_labels:
#    sanitized_item = sanitize_element(item)
#    print(element_electronegativity[sanitized_item],' ', end='')

print(f"{weighted_elec_neg:.3f}"," ", end='')
print(f"{alpha:.2f}"," ", end='')
print(f"{beta:.2f}"," ", end='')
print(f"{gamma:.2f}"," ", end='')
print(f"{a:.3f}"," " , end='')
print(f"{b:.3f}"," " , end='')
print(f"{c:.3f}"," " , end='')
print(f"{volume:.3f}", " ", end='')
print(coordination_number, " ", end='')
#kj_mol = ((float(formation_enthalpy_value)*temp_sum)/gcd)*96.487
#print(f"{kj_mol:.3f}"," ", end='')
print(formation_enthalpy_value)
#print(type(formation_enthalpy_value))


