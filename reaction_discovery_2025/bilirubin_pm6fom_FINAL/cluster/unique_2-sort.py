def create_unique_geometries_file(xyz_file_path, txt_file_path, output_file_path):
    # Parse all geometries and related information
    with open(xyz_file_path, 'r') as xyz_file:
        xyz_lines = xyz_file.readlines()

    with open(txt_file_path, 'r') as txt_file:
        txt_lines = txt_file.readlines()

    geometry_start_indices = [i for i, line in enumerate(txt_lines) if line.startswith("Geometry")]

    # Helper function to extract a geometry's dihedrals and H-bond identities
    def extract_geometry_info(index):
        geom_start = geometry_start_indices[index]
        dihedral_line = txt_lines[geom_start + 2].split(":")[1].strip()
        dihedral_values = [float(x) for x in dihedral_line.split()]

        hbond_lines = []
        hbond_line_start = geom_start + 3
        while hbond_line_start < len(txt_lines) and "Donor=" in txt_lines[hbond_line_start]:
            hbond = txt_lines[hbond_line_start].strip()
            donor_acceptor = hbond.split(",")[0].strip() + ", " + hbond.split(",")[1].strip()
            hbond_lines.append(donor_acceptor)
            hbond_line_start += 1

        return dihedral_values, hbond_lines

    def extract_energy(geom_index):
        energy_line = xyz_lines[1 + 81 * geom_index].strip()
        if "Energy" in energy_line:
            raw_energy = float(energy_line.split("Energy")[-1].strip())
            converted_energy = raw_energy * 627.509  # Convert Hartree to kcal/mol
            return raw_energy, converted_energy
        return None, None

    # Initialize output and process the geometries
    unique_geometries = []  # Store unique geometries as (energy, dihedrals, hbonds)
    with open(output_file_path, 'w') as output_file:
        for i in range(len(geometry_start_indices)):
            # Extract energy, dihedrals, and hydrogen bonds for geometry `i`
            raw_energy, converted_energy = extract_energy(i)
            dihedrals, hbonds = extract_geometry_info(i)

            if not unique_geometries:
                # Always save the first geometry
                unique_geometries.append((converted_energy, dihedrals, hbonds))
                output_file.writelines(xyz_lines[81 * i : 81 * (i + 1)])
            else:
                # Compare only with the last saved geometry
                last_energy, last_dihedrals, last_hbonds = unique_geometries[-1]
                energy_diff = abs(converted_energy - last_energy)
                dihedral_diffs = [abs(d1 - d2) for d1, d2 in zip(dihedrals, last_dihedrals)]
                max_dihedral_diff = max(dihedral_diffs)
                hbonds_different = set(hbonds) != set(last_hbonds)  # Compare donor-acceptor pairs only

                # Save geometry if it meets uniqueness criteria
                if (
                    energy_diff > 5
                    or max_dihedral_diff > 45
                    or hbonds_different
                ):
                    unique_geometries.append((converted_energy, dihedrals, hbonds))
                    output_file.writelines(xyz_lines[81 * i : 81 * (i + 1)])

    return len(unique_geometries)

# Specify file paths
xyz_file_path = "sorted_all_2.xyz"  # Replace with your XYZ file path
txt_file_path = "2_unique_1-list.txt"  # Replace with your TXT file path
output_file_path = "2_unique_2-sort.xyz"  # Replace with your desired output file path

# Create the output file
num_unique_geometries = create_unique_geometries_file(xyz_file_path, txt_file_path, output_file_path)
print(f"Unique geometries saved to {output_file_path}. Total unique geometries: {num_unique_geometries}")



