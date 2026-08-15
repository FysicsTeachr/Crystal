import numpy as np
import math

def read_xyz(filename):
    """
    Reads an XYZ trajectory file and returns a list of geometries (each as a numpy array).
    """
    geometries = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        i = 0
        while i < len(lines):
            num_atoms = int(lines[i].strip())
            i += 2  # Skip the atom count and comment line
            geometry = []
            for _ in range(num_atoms):
                parts = lines[i].strip().split()
                geometry.append([float(x) for x in parts[1:4]])
                i += 1
            geometries.append(np.array(geometry))
    return geometries

def calculate_dihedral(p1, p2, p3, p4):
    """
    Calculates the dihedral angle (in degrees) defined by four points.
    """
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3

    b1xb2 = np.cross(b1, b2)
    b2xb3 = np.cross(b2, b3)
    b1xb2_norm = b1xb2 / np.linalg.norm(b1xb2)
    b2xb3_norm = b2xb3 / np.linalg.norm(b2xb3)

    # Calculate the angle between b1xb2 and b2xb3
    cos_angle = np.dot(b1xb2_norm, b2xb3_norm)
    cos_angle = np.clip(cos_angle, -1.0, 1.0)  # Handle numerical precision issues
    angle = np.arccos(cos_angle)

    # Determine the sign of the angle
    if np.dot(b1xb2_norm, b3) < 0:
        angle = -angle

    return np.degrees(angle)

def compute_dihedrals(geometries, dihedral_indices):
    """
    Computes the specified dihedrals for each geometry.

    Parameters:
    - geometries: List of geometries (numpy arrays).
    - dihedral_indices: List of lists containing indices defining the dihedrals.

    Returns:
    - dihedral_values: List of dictionaries containing dihedrals for each geometry.
    """
    dihedral_values = []
    for geometry in geometries:
        dihedrals = {}
        for idx, atoms in dihedral_indices.items():
            p1, p2, p3, p4 = geometry[atoms[0]], geometry[atoms[1]], geometry[atoms[2]], geometry[atoms[3]]
            dihedrals[idx] = calculate_dihedral(p1, p2, p3, p4)
        dihedral_values.append(dihedrals)
    return dihedral_values

# Main program
if __name__ == "__main__":
    input_file = "yesfile.xyz"
    output_file = "dihedrals.txt"

    # Dihedral definitions
    dihedral_indices = {
        1: [41, 0, 3, 4],
        2: [4, 6, 23, 25],
        3: [6, 23, 25, 26],
        4: [3, 0, 41, 42],
        5: [42, 44, 61, 63],
        6: [44, 61, 63, 64]
    }

    geometries = read_xyz(input_file)
    dihedral_values = compute_dihedrals(geometries, dihedral_indices)

    # Write results to a file
    with open(output_file, 'w') as f:
        for dihedrals in dihedral_values:
            line = " ".join(f"{angle:.2f}" for angle in dihedrals.values())
            f.write(line + "\n")

    print(f"Dihedral angles have been written to {output_file}.")
