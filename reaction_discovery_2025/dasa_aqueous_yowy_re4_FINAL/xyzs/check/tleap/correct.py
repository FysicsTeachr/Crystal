def fix_water_residues(pdb_path, output_path):
    with open(pdb_path, 'r') as file:
        lines = file.readlines()

    fixed_lines = []
    current_residue = None
    hydrogen_count = 0

    for line in lines:
        if line.startswith("HETATM"):
            atom_name = line[12:16].strip()
            residue_name = line[17:20].strip()
            residue_number = line[22:26].strip()

            if residue_name == "HOH":
                if atom_name == "O":
                    # Reset hydrogen counter for a new water molecule
                    hydrogen_count = 0
                    current_residue = residue_number
                    fixed_lines.append(line)
                elif atom_name.startswith("H"):
                    # Name hydrogens as H1, H2
                    hydrogen_count += 1
                    new_atom_name = f"H{hydrogen_count}"
                    line = line[:12] + new_atom_name.ljust(4) + line[16:22] + current_residue.rjust(4) + line[26:]
                    fixed_lines.append(line)
            else:
                fixed_lines.append(line)

    with open(output_path, 'w') as output_file:
        output_file.writelines(fixed_lines)

# Define the input and output file paths
input_pdb_path = '179.pdb'
output_pdb_path = '179_fixed.pdb'

# Run the function to fix water residues
fix_water_residues(input_pdb_path, output_pdb_path)

