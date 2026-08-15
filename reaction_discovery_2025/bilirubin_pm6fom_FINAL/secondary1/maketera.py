import numpy as np
from scipy.spatial.transform import Rotation as R

def read_xyz_file(filepath):
    with open(filepath, 'r') as file:
        num_atoms = int(file.readline().strip())  # Read and store the number of atoms (optional)
        file.readline()  # Skip the second line (comment)

        # Read atom types and coordinates
        coords = []
        for _ in range(num_atoms):
            parts = file.readline().split()
            atom_type = parts[0]
            x, y, z = map(float, parts[1:4])
            coords.append((atom_type, np.array([x, y, z])))

    return coords

def calculate_dihedral(coords, indices):
    """
    Calculate the dihedral angle for the given atom indices.
    """
    p1, p2, p3, p4 = [coords[i][1] for i in indices]

    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3

    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    n1 /= np.linalg.norm(n1)
    n2 /= np.linalg.norm(n2)

    cos_theta = np.dot(n1, n2)
    m1 = np.cross(n1, n2)
    sin_theta = np.dot(m1, b2 / np.linalg.norm(b2))

    angle = np.arctan2(sin_theta, cos_theta)
    angle_degrees = np.degrees(angle)

    return angle_degrees

def calculate_angle(coords, indices):
    """
    Calculate the angle for the given atom indices.
    """
    p1, p2, p3 = [coords[i][1] for i in indices]

    ba = p1 - p2
    bc = p3 - p2

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(np.clip(cosine_angle, -1.0, 1.0))
    angle_degrees = np.degrees(angle)

    return angle_degrees

def calculate_all_dihedrals(coords):
    dihedral_mappings = {
    1: [49,46,34,10],
    2: [10,31,29,23],
    3: [31,29,23,8],
    4: [34,46,49,12],
    5: [12,52,64,66],
    6: [52,64,66,14]
    }

    dihedrals = {}
    for key, indices in dihedral_mappings.items():
        dihedrals[key] = calculate_dihedral(coords, [i for i in indices])
    return dihedrals

def calculate_all_angles(coords):
    angle_mappings = {
#        3: [15, 16, 18],
#        5: [16, 18, 20],
#        7: [18, 20, 22],
#        9: [20, 22, 24]
    }

    angles = {}
    for key, indices in angle_mappings.items():
        angles[key] = calculate_angle(coords, [i for i in indices])
    return angles

def process_input_file(input_filepath, coords, xt):
    initial_dihedrals = calculate_all_dihedrals(coords)
#    initial_dihedrals = {
#        1: 78.29,
#        2: -15.17,
#        3: 170.10,
#        4: -23.13,
#        5: 10.16,
#        6: -163.31
#    }
    initial_angles = calculate_all_angles(coords)

    with open(input_filepath, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line
        numbers = [float(num) for line in lines for num in line.split()]

        final_dihedrals = initial_dihedrals.copy()
        final_angles = initial_angles.copy()
        active_atoms_mapping = {}

        dihedral_mappings = {
        1: [49,46,34,10],
        2: [10,31,29,23],
        3: [31,29,23,8],
        4: [34,46,49,12],
        5: [12,52,64,66],
        6: [52,64,66,14]
        }

        angle_mappings = {
            # 9: [20, 22, 24]
        }

        for i, value in enumerate(numbers, start=1):
            if value != 0.0:
                if i in initial_dihedrals:
                    final_dihedrals[i] += value  # Add the input value
                    active_atoms = [index + 1 for index in dihedral_mappings[i]]
                    active_atoms_mapping[i] = active_atoms
                elif i in initial_angles:
                    final_angles[i] += value  # Add the input value
                    active_atoms = [index + 1 for index in angle_mappings[i]]
                    active_atoms_mapping[i] = active_atoms

    all_active_atoms = sorted(set([atom for atoms in active_atoms_mapping.values() for atom in atoms]))

    for step in range(1, xt + 1):
        output_script = f'constraints_{step}.txt'
        generate_terachem_script(
            coords,
            initial_dihedrals,
            final_dihedrals,
            initial_angles,
            final_angles,
            output_script,
            dihedral_mappings,
            angle_mappings,
            xt,
            step
        )

def generate_terachem_script(coords, initial_dihedrals, final_dihedrals, initial_angles, final_angles, output_script, dihedral_mappings, angle_mappings, xt, step):
    coords_input = "pregeom.xyz" if step == 1 else f"pregeom.xyz"

    terachem_script = f"""
$set
"""

    for key in dihedral_mappings:
        dihedral_value = initial_dihedrals[key] + (final_dihedrals[key] - initial_dihedrals[key]) * step / xt
        dihedral_atoms = " ".join(str(i + 1) for i in dihedral_mappings[key])
        terachem_script += f"dihedral {dihedral_atoms} {dihedral_value:.2f}\n"

    for key in angle_mappings:
        angle_value = initial_angles[key] + (final_angles[key] - initial_angles[key]) * step / xt
        angle_atoms = " ".join(str(i + 1) for i in angle_mappings[key])
        terachem_script += f"angle {angle_atoms} {angle_value:.2f}\n"

#    terachem_script += "$end\n"

    with open(output_script, 'w') as file:
        file.write(terachem_script)

    print(f"TeraChem input file written to {output_script}")

# Example usage
xyz_file_path = 'pregeom.xyz'
input_file_path = 'coord1.txt'
xt = 3  # Number of steps

# Read coordinates
coords = read_xyz_file(xyz_file_path)

# Process the input file and generate multiple ChemShell scripts for each step
process_input_file(input_file_path, coords, xt)
