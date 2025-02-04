import os
import math

def read_coordinates(file_path):
    coordinates = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith(('ATOM', 'HETATM')):
                parts = line.split()
                atom_id = int(parts[1])
                atom_name = parts[2]
                x = float(parts[5])
                y = float(parts[6])
                z = float(parts[7])
                coordinates.append((atom_id, atom_name, x, y, z))
    return coordinates

def read_connectivity(file_path):
    connectivity = [] 
    with open(file_path, "r") as file: 
        for line in file:  
            if line.startswith('CONECT'):  
                atom_numbers = line.split()[1:]  
                atoms = [int(atom) for atom in atom_numbers]  
                if len(atoms) > 1:  
                    connectivity.append(atoms)  
    return connectivity

def print_connectivity_and_coordinates(coordinates, connectivity):
    print("\nAtom Coordinates:")
    for atom in coordinates:
        print(f"Atom ID: {atom[0]}, Name: {atom[1]}, Coordinates: ({atom[2]:.3f}, {atom[3]:.3f}, {atom[4]:.3f})")
    
    print("\nConnectivity extracted:")
    for bond in connectivity:
        print(bond)

def calculate_bond_length(atom1, atom2):
    x1, y1, z1 = atom1[2], atom1[3], atom1[4]
    x2, y2, z2 = atom2[2], atom2[3], atom2[4]
    return math.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)


def calculate_bond_angle(atom1, atom2, atom3):
    x1, y1, z1 = atom1[2], atom1[3], atom1[4]
    x2, y2, z2 = atom2[2], atom2[3], atom2[4]
    x3, y3, z3 = atom3[2], atom3[3], atom3[4]
    
    vec1 = (x1 - x2, y1 - y2, z1 - z2)
    vec2 = (x3 - x2, y3 - y2, z3 - z2)
    
    dot_product = vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2]
    
    mag1 = math.sqrt(vec1[0]**2 + vec1[1]**2 + vec1[2]**2)
    mag2 = math.sqrt(vec2[0]**2 + vec2[1]**2 + vec2[2]**2)
    
    cos_theta = dot_product / (mag1 * mag2)
    angle = math.acos(cos_theta) * (180.0 / math.pi)
    return angle

def calculate_dihedral_angle(atom1, atom2, atom3, atom4):
    x1, y1, z1 = atom1[2], atom1[3], atom1[4]
    x2, y2, z2 = atom2[2], atom2[3], atom2[4]
    x3, y3, z3 = atom3[2], atom3[3], atom3[4]
    x4, y4, z4 = atom4[2], atom4[3], atom4[4]
    
    vec1 = (x2 - x1, y2 - y1, z2 - z1)  #vec from atom1 to atom2
    vec2 = (x3 - x2, y3 - y2, z3 - z2)  #vec from atom2 to atom3
    vec3 = (x4 - x3, y4 - y3, z4 - z3)  #vec from atom3 to atom4
    
    normal1 = [
        vec1[1] * vec2[2] - vec1[2] * vec2[1],
        vec1[2] * vec2[0] - vec1[0] * vec2[2],
        vec1[0] * vec2[1] - vec1[1] * vec2[0]
    ]
    
    normal2 = [
        vec2[1] * vec3[2] - vec2[2] * vec3[1],
        vec2[2] * vec3[0] - vec2[0] * vec3[2],
        vec2[0] * vec3[1] - vec2[1] * vec3[0]
    ]
    
    mag1 = math.sqrt(normal1[0]**2 + normal1[1]**2 + normal1[2]**2)
    mag2 = math.sqrt(normal2[0]**2 + normal2[1]**2 + normal2[2]**2)
    
    dot_product = normal1[0] * normal2[0] + normal1[1] * normal2[1] + normal1[2] * normal2[2]
    cos_phi = dot_product / (mag1 * mag2)
    dihedral_angle = math.acos(cos_phi) * (180.0 / math.pi)
    
    cross_product = (
        normal1[1] * normal2[2] - normal1[2] * normal2[1],
        normal1[2] * normal2[0] - normal1[0] * normal2[2],
        normal1[0] * normal2[1] - normal1[1] * normal2[0]
    )
    
    if cross_product[2] < 0:
        dihedral_angle = -dihedral_angle
    return dihedral_angle

def write_xyz_file(coordinates, output_file):
    with open(output_file, 'w') as output:  
        output.write(f"{len(coordinates)}\n\n")  
        for atom in coordinates:  
            output.write(f"{atom[1]}    {atom[2]:.5f}    {atom[3]:.5f}    {atom[4]:.5f}\n")  
    print(f"XYZ file created successfully: {output_file}")


def atom_ids(coordinates, atom_id):
    for atom in coordinates:
        if atom[0] == atom_id:  
            return atom
    return None 

def compare_properties(file1_coordinates, atom_id1, file2_coordinates=None, atom_id2=None):
    if file2_coordinates and atom_id2:
        if len(atom_id1) == 2 and len(atom_id2) == 2:
            atom1_1 = atom_ids(file1_coordinates, atom_id1[0])
            atom2_1 = atom_ids(file1_coordinates, atom_id1[1])
            atom1_2 = atom_ids(file2_coordinates, atom_id2[0])
            atom2_2 = atom_ids(file2_coordinates, atom_id2[1])

            if atom1_1 and atom2_1 and atom1_2 and atom2_2:
                bond_length_1 = calculate_bond_length(atom1_1, atom2_1)
                bond_length_2 = calculate_bond_length(atom1_2, atom2_2)
                difference_bond_length = abs(bond_length_1 - bond_length_2)
                print(f"Bond Length for atoms in File 1: {bond_length_1:.2f} Å")
                print(f"Bond Length for atoms in File 2: {bond_length_2:.2f} Å")
                if difference_bond_length > 0.05:
                    print(f"Bond length difference: {difference_bond_length:.2f} Å")
                else:
                    print("No significant bond length difference")
        
        elif len(atom_id1) == 3 and len(atom_id2) == 3:
            atom1_1 = atom_ids(file1_coordinates, atom_id1[0])
            atom2_1 = atom_ids(file1_coordinates, atom_id1[1])
            atom3_1 = atom_ids(file1_coordinates, atom_id1[2])
            atom1_2 = atom_ids(file2_coordinates, atom_id2[0])
            atom2_2 = atom_ids(file2_coordinates, atom_id2[1])
            atom3_2 = atom_ids(file2_coordinates, atom_id2[2])

            if atom1_1 and atom2_1 and atom3_1 and atom1_2 and atom2_2 and atom3_2:
                bond_angle_1 = calculate_bond_angle(atom1_1, atom2_1, atom3_1)
                bond_angle_2 = calculate_bond_angle(atom1_2, atom2_2, atom3_2)
                difference_bond_angle = abs(bond_angle_1 - bond_angle_2)
                print(f"Bond Angle for atoms in File 1: {bond_angle_1:.1f}°")
                print(f"Bond Angle for atoms in File 2: {bond_angle_2:.1f}°")
                if difference_bond_angle > 1:
                    print(f"Bond angle difference : {difference_bond_angle:.1f}°")
                else:
                    print("No significant bond angle difference")
            
        elif len(atom_id1) == 4 and len(atom_id2) == 4:
            atom1_1 = atom_ids(file1_coordinates, atom_id1[0])
            atom2_1 = atom_ids(file1_coordinates, atom_id1[1])
            atom3_1 = atom_ids(file1_coordinates, atom_id1[2])
            atom4_1 = atom_ids(file1_coordinates, atom_id1[3])
            atom1_2 = atom_ids(file2_coordinates, atom_id2[0])
            atom2_2 = atom_ids(file2_coordinates, atom_id2[1])
            atom3_2 = atom_ids(file2_coordinates, atom_id2[2])
            atom4_2 = atom_ids(file2_coordinates, atom_id2[3])

            if atom1_1 and atom2_1 and atom3_1 and atom4_1 and atom1_2 and atom2_2 and atom3_2 and atom4_2:
                dihedral_1 = calculate_dihedral_angle(atom1_1, atom2_1, atom3_1, atom4_1)
                dihedral_2 = calculate_dihedral_angle(atom1_2, atom2_2, atom3_2, atom4_2)
                difference_dihedral_angle = abs(dihedral_1 - dihedral_2)
                print(f"Dihedral Angle for atoms in File 1: {dihedral_1:.0f}°")
                print(f"Dihedral Angle for atoms in File 2: {dihedral_2:.0f}°")
                if difference_dihedral_angle > 5:
                    print(f"Dihedral angle difference : {difference_dihedral_angle: .0f}°")
                else:
                    print("No significant dihedral angle difference")
    else:
        if len(atom_id1) == 2:
            atom1_1 = atom_ids(file1_coordinates, atom_id1[0])
            atom2_1 = atom_ids(file1_coordinates, atom_id1[1])

            if atom1_1 and atom2_1:
                bond_length_1 = calculate_bond_length(atom1_1, atom2_1)
                print(f"Bond Length for atoms in File 1: {bond_length_1:.2f}°")
        elif len(atom_id1) == 3:
            atom1_1 = atom_ids(file1_coordinates, atom_id1[0])
            atom2_1 = atom_ids(file1_coordinates, atom_id1[1])
            atom3_1 = atom_ids(file1_coordinates, atom_id1[2])

            if atom1_1 and atom2_1 and atom3_1:
                bond_angle_1 = calculate_bond_angle(atom1_1, atom2_1, atom3_1)
                print(f"Bond Angle for atoms in File 1: {bond_angle_1:.1f}°")
        elif len(atom_id1) == 4:
            atom1_1 = atom_ids(file1_coordinates, atom_id1[0])
            atom2_1 = atom_ids(file1_coordinates, atom_id1[1])
            atom3_1 = atom_ids(file1_coordinates, atom_id1[2])
            atom4_1 = atom_ids(file1_coordinates, atom_id1[3])

            if atom1_1 and atom2_1 and atom3_1 and atom4_1:
                dihedral_1 = calculate_dihedral_angle(atom1_1, atom2_1, atom3_1, atom4_1)
                print(f"Dihedral Angle for atoms in File 1: {dihedral_1: .0f}°")
def input_part():
    file1_path = input("Enter path to the first PDB file: ").strip()

    if not os.path.exists(file1_path):
        print("The file path is incorrect.")
        return

    file1_coordinates = read_coordinates(file1_path)
    connectivity1 = read_connectivity(file1_path)
    print_connectivity_and_coordinates(file1_coordinates, connectivity1)

    output_xyz1 = file1_path.replace(".pdb", ".xyz")
    write_xyz_file(file1_coordinates, output_xyz1)

    compare_second = input("Do you want to compare with a second file? (yes/no): ").strip()
    
    if compare_second == "yes":
        file2_path = input("Enter path to the second PDB file: ").strip()
        if not os.path.exists(file2_path):
            print("The second file path is incorrect.")
            return
        
        file2_coordinates = read_coordinates(file2_path)
        connectivity2 = read_connectivity(file2_path)
        print_connectivity_and_coordinates(file2_coordinates, connectivity2)

        output_xyz2 = file2_path.replace(".pdb", ".xyz")
        write_xyz_file(file2_coordinates, output_xyz2)

        atom_id1 = list(map(int, input("Enter atom IDs for the first file: ").split(',')))
        atom_id2 = list(map(int, input("Enter atom IDs for the second file: ").split(',')))
        compare_properties(file1_coordinates, atom_id1, file2_coordinates, atom_id2)
    
    elif compare_second == "no":
        atom_id1 = list(map(int, input("Enter atom IDs for the first file: ").split(',')))
        compare_properties(file1_coordinates, atom_id1)

input_part()