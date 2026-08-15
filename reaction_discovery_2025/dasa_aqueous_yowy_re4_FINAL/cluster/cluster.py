import numpy as np

def calculate_distance(atom1, atom2):
    """Calculate the distance between two atoms."""
    return np.linalg.norm(atom1 - atom2)

def refined_calculate_dihedral(p1, p2, p3, p4):
    """
    Calculate the dihedral angle (in degrees) using a refined approach.
    Ensures accurate calculation of angles within [-180, 180].
    """
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3

    b1 = b1 / np.linalg.norm(b1)
    b2 = b2 / np.linalg.norm(b2)
    b3 = b3 / np.linalg.norm(b3)

    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    n1 = n1 / np.linalg.norm(n1)
    n2 = n2 / np.linalg.norm(n2)

    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)
    angle = np.degrees(np.arctan2(y, x))

    return angle

def categorize_dihedrals(dihedrals):
    """Categorize dihedrals as 'cis' or 'trans'."""
    return ['cis' if abs(d) < 90 else 'trans' for d in dihedrals]

def process_geometries(input_file, output_file, lines_per_block=1214, required_lines=44, bond_cutoff=1.6):
    # Dihedral and bond indices
    dihedral_indices = [
        (14, 15, 16, 18),
        (15, 16, 18, 20),
        (16, 18, 20, 22),
        (18, 20, 22, 24)
    ]
    bonded_atoms_indices = [11, 12, 13, 14, 15, 16, 18, 20, 22, 24, 26, 27, 30, 34, 37]

    # Read the input file
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    unique_geometries = []
    geometry_to_original_mapping = {}

    with open(output_file, 'w') as outfile:
        for block_idx in range(0, len(lines), lines_per_block):
            current_block = lines[block_idx:block_idx + required_lines]
            current_block[0] = '42\n'  # Update atom count
            num_atoms = int(current_block[0])
            atom_positions = np.array(
                [list(map(float, l.split()[1:])) for l in current_block[2:]]
            )

            # Calculate dihedrals and bonds
            dihedrals = categorize_dihedrals([
                refined_calculate_dihedral(
                    atom_positions[i1], atom_positions[i2], atom_positions[i3], atom_positions[i4]
                )
                for i1, i2, i3, i4 in dihedral_indices
            ])

            bonds = {index: [] for index in bonded_atoms_indices}
            for idx1 in bonded_atoms_indices:
                for idx2 in range(num_atoms):
                    if idx1 != idx2:
                        dist = calculate_distance(atom_positions[idx1 - 1], atom_positions[idx2 - 1])
                        if dist <= bond_cutoff:
                            bonds[idx1].append(idx2 + 1)

            # Check uniqueness
            is_unique = True
            for geom_idx, geometry in enumerate(unique_geometries, start=1):
                if geometry['dihedrals'] == dihedrals and geometry['bonds'] == bonds:
                    is_unique = False
                    geometry_to_original_mapping[geom_idx] = block_idx // lines_per_block + 1
                    break

            if is_unique:
                unique_geometries.append({'dihedrals': dihedrals, 'bonds': bonds})
                geometry_number = block_idx // lines_per_block + 1
                geometry_to_original_mapping[len(unique_geometries)] = geometry_number
                current_block[1] = current_block[1].strip() + f" Geometry {geometry_number}\n"
                outfile.writelines(current_block)

if __name__ == "__main__":
    input_file = "all_m.xyz"  # Replace with the full path to your input file
    output_file = "unique_geometries.xyz"  # Output file for unique geometries
    process_geometries(input_file, output_file)

