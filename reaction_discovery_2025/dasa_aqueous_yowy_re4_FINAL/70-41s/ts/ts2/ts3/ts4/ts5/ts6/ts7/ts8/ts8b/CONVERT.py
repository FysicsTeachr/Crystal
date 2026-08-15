def rst7_to_xyz(rst7_file, xyz_file, output_file, num_atoms=42):
    """
    Convert the first 42 atoms from an rst7 file to an xyz file format.
    The atom types in the xyz file are used to set the order, and they are capitalized.

    Parameters:
        rst7_file (str): Path to the rst7 file.
        xyz_file (str): Path to the xyz file (to get atom types and order).
        output_file (str): Path to the output xyz file.
        num_atoms (int): Number of atoms to process (default is 42).
    """
    # Read the rst7 file
    with open(rst7_file, 'r') as rst7:
        rst7_lines = rst7.readlines()
    
    # Extract coordinates from rst7 (skipping the first two lines)
    rst7_coords = []
    for line in rst7_lines[2:]:
        # Convert each set of three numbers to floats
        coords = [float(x) for x in line.split()]
        rst7_coords.extend([coords[i:i+3] for i in range(0, len(coords), 3)])
        if len(rst7_coords) >= num_atoms:
            break
    rst7_coords = rst7_coords[:num_atoms]

    # Read the xyz file to get atom types
    with open(xyz_file, 'r') as xyz:
        xyz_lines = xyz.readlines()

    # Extract atom types from xyz
    atom_types = []
    for line in xyz_lines[2:2+num_atoms]:  # Start after the first two header lines
        atom_type = line.split()[0].upper()
        atom_types.append(atom_type)

    # Write the converted coordinates to the output file in xyz format
    with open(output_file, 'w') as out_xyz:
        out_xyz.write(f"{num_atoms}\n")  # Write number of atoms
        out_xyz.write("Converted from rst7\n")  # Write comment/header line
        for atom, coords in zip(atom_types, rst7_coords):
            out_xyz.write(f"{atom}   {coords[0]:.6f}   {coords[1]:.6f}   {coords[2]:.6f}\n")
    
    print(f"Conversion complete. Output written to {output_file}")

# Example usage
rst7_to_xyz('scr.neb1/neb_1.rst7', 'optimized.xyz', 'gas_0.xyz', num_atoms=42)

