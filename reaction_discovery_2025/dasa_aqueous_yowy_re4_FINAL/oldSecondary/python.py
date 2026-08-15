from ase import Atoms
from ase.io import read, write
import numpy as np

def rotate_atoms(atoms, indices, angle, axis, center):
    """Rotate a subset of atoms around a given axis and center."""
    rotation_matrix = rotation_matrix_from_axis_angle(axis, angle)
    for i in indices:
        atoms.positions[i] = np.dot(rotation_matrix, atoms.positions[i] - center) + center

def rotation_matrix_from_axis_angle(axis, angle):
    """Generate a rotation matrix from an axis and an angle."""
    axis = axis / np.linalg.norm(axis)
    a = np.cos(angle / 2.0)
    b, c, d = -axis * np.sin(angle / 2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

# Function to read angles and translation vector from a file
def read_angles_and_translation(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()[1:]  # Skip the first line
        # Flatten the list of lists into a single list of numbers
        numbers = [float(num) for line in lines for num in line.split()]
        angles_degrees = numbers[:11]
        H_translation = numbers[11:]
        return angles_degrees, H_translation

# Read the original XYZ file
file_path = 'pregeom.xyz'  # Change this to the path of your file
atoms = read(file_path)

# Read angles and translation vector from coord.txt
angles_degrees, H_translation = read_angles_and_translation('coord1.txt')
angles = [angle * np.pi / 180 for angle in angles_degrees]  # Convert degrees to radians

# Step 1: Rotate atoms 16, 17, ..., 40 about the axis through 15 parallel to 13-14
a13 = atoms.positions[13]
a15 = atoms.positions[15]
a14 = atoms.positions[14]
axis1 = a14 - a13
rotate_atoms(atoms, list(range(16, 41)), angles[0], axis1, a15)

# Step 2: Rotate atoms 17, 18, 19, ..., 40 about the axis through 15-16
a16 = atoms.positions[16]
axis3 = a16 - a15
rotate_atoms(atoms, list(range(17, 41)), angles[1], axis3, a15)

# Step 3: Rotate atoms 18, 19, ..., 40 about the axis through 16, defined by the cross product of 16-18 and 15-16
a18 = atoms.positions[18]
v1 = a18 - a16
v2 = a15 - a16
axis2 = np.cross(v1, v2)
rotate_atoms(atoms, list(range(18, 41)), angles[2], axis2, a16)
rotate_atoms(atoms, [17], angles[1]/2, axis2, a16)

# Step 4: Rotate atoms 20, 21, ..., 40 about the axis through 16-18
axis5 = a18 - a16
rotate_atoms(atoms, list(range(19, 41)), angles[3], axis5, a16)

# Step 5: Rotate atoms 20, 21, ..., 40 about the axis through 18, defined by the cross product of 18-20 and 16-18
a20 = atoms.positions[20]
v1 = a20 - a18
v2 = a16 - a18
axis4 = np.cross(v1, v2)
rotate_atoms(atoms, list(range(20, 41)), angles[4], axis4, a18)
rotate_atoms(atoms, [19], angles[3]/2, axis4, a18)

# Step 6: Rotate atoms 22, 23, ..., 40 about the axis through 18-20
axis7 = a20 - a18
rotate_atoms(atoms, list(range(21, 41)), angles[5], axis7, a18)

# Step 7: Rotate atoms 22, 23, ..., 40 about the axis through 20, defined by the cross product of 20-22 and 18-20
a22 = atoms.positions[22]
v1 = a22 - a20
v2 = a18 - a20
axis6 = np.cross(v1, v2)
rotate_atoms(atoms, list(range(22, 41)), angles[6], axis6, a20)
rotate_atoms(atoms, [21], angles[5]/2, axis6, a20)

# Step 8: Rotate atoms 24, 25, ..., 40 about the axis through 20-22
axis9 = a22 - a20
rotate_atoms(atoms, list(range(23, 41)), angles[7], axis9, a20)

# Step 9: Rotate atoms 24, 25, ..., 40 about the axis through 22, defined by the cross product of 20-22 and 22-24
a24 = atoms.positions[24]
v1 = a24 - a22
v2 = a20 - a22
axis8 = np.cross(v1, v2)
rotate_atoms(atoms, list(range(24, 41)), angles[8], axis8, a22)
rotate_atoms(atoms, [23], angles[7]/2, axis8, a22)

# Step 10: Rotate atoms 26, 27, ..., 40 about the axis through 22-24
axis11 = a24 - a22
rotate_atoms(atoms, list(range(25, 41)), angles[9], axis11, a22)

# Step 11: Rotate atoms 26, 27, ..., 40 about the axis through 24, defined by the cross product of 22-24 and 24-26
a26 = atoms.positions[26]
v1 = a26 - a24
v2 = a22 - a24
axis10 = np.cross(v1, v2)
rotate_atoms(atoms, list(range(26, 41)), angles[10], axis10, a24)
rotate_atoms(atoms, [25], angles[9]/2, axis10, a24)

# Translate the position of the 41st atom
atoms.positions[41] += H_translation  # Note: index 41 in human terms is index 40 in zero-based indexing

# Save the transformed coordinates to an XYZ file
output_file_path = 'tgeom.xyz'  # Change this to your desired output path
write(output_file_path, atoms)

