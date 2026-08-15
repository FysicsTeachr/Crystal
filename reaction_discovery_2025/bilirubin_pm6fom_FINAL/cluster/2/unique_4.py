import os

# Function to process geometry data
def process_geometry_data(file_path):
    geometries = []
    with open(file_path, 'r') as file:
        current_geometry = {}
        for line in file:
            line = line.strip()
            if line.startswith("Geometry"):
                if current_geometry:
                    geometries.append(current_geometry)
                current_geometry = {
                    "geometry_number": int(line.split()[1][:-1]), 
                    "hydrogen_bonds": [], 
                    "dihedrals": []
                }
            elif "Dihedrals" in line:
                current_geometry["dihedrals"] = tuple(float(x) for x in line.split(":")[1].split())
            elif "Donor=" in line:
                donor_acceptor = line.split(",")[0].split("Donor=")[1]
                acceptor = line.split("Acceptor=")[1].split(",")[0]
                current_geometry["hydrogen_bonds"].append((int(donor_acceptor.split('=')[1]), int(acceptor)))
        if current_geometry:
            geometries.append(current_geometry)
    return geometries

# Function to remove duplicate geometries using pairwise comparisons
def remove_duplicate_geometries(geometries):
    unique_geometries = []
    seen = set()

    for geometry in geometries:
        # Create a hashable representation of each geometry (tuple of dihedrals and sorted hydrogen bonds)
        geom_key = (geometry["dihedrals"], tuple(sorted(geometry["hydrogen_bonds"])))
        
        if geom_key not in seen:
            seen.add(geom_key)
            unique_geometries.append(geometry)

    return unique_geometries

# Function to extract energies from blocks of 81 lines
def extract_energies(file_path, block_size=81, energy_line=2):
    energies = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for i in range(0, len(lines), block_size):
            block = lines[i:i+block_size]
            if len(block) >= energy_line:
                energy_line_content = block[energy_line - 1].strip()
                energy = float(energy_line_content.split()[-1])  # Assuming energy is the last value in the line
                energies.append(energy)
    return energies

# Function to categorize geometries
def categorize_geometries(geometries, energies):
    categories = {'z,z': [], 'e,z': [], 'z,e': [], 'e,e': []}
    
    for geometry in geometries:
        dihedral_3, dihedral_6 = geometry['dihedrals'][2], geometry['dihedrals'][5]
        
        if -45 <= dihedral_3 <= 45 and -45 <= dihedral_6 <= 45:
            categories['z,z'].append(geometry)
        elif (-45 <= abs(dihedral_3) <= 45 and 135 <= abs(dihedral_6) <= 180):
            categories['z,e'].append(geometry)
        elif (135 <= abs(dihedral_3) <= 180 and -45 <= abs(dihedral_6) <= 45):
            categories['e,z'].append(geometry)
        elif (135 <= abs(dihedral_3) <= 180 and 135 <= abs(dihedral_6) <= 180):
            categories['e,e'].append(geometry)

    # Align energies with geometries based on geometry number
    for category, geometries_list in categories.items():
        for geometry in geometries_list:
            geometry_number = geometry["geometry_number"]
            geometry["energy"] = energies[geometry_number] if geometry_number < len(energies) else None

    return categories

# Function to write output to a text file
def write_output(categories, output_file):
    with open(output_file, 'w') as file:
        file.write("Geometry Number\tCategory\tHydrogen Bonds\tDihedral Values\tEnergy\n")
        for category, geometries in categories.items():
            for geometry in geometries[:100]:  # Restrict to first 100 geometries per category
                hydrogen_bonds = ", ".join(map(str, geometry["hydrogen_bonds"]))
                dihedrals = ", ".join(map(str, geometry["dihedrals"]))
                energy = geometry["energy"]
                file.write(f"{geometry['geometry_number']}\t{category}\t{hydrogen_bonds}\t{dihedrals}\t{energy}\n")

# Input file paths
input_geometry_file = "unique_2-info.txt"  # Replace with the first input file path
input_energy_file = "unique_2-sort.xyz"  # Replace with the second input file path
output_file_path = "unique_4-few.xyz"

# Main script execution
if __name__ == "__main__":
    geometries_data = process_geometry_data(input_geometry_file)
    
    # Remove duplicate geometries (pairwise comparison)
    unique_geometries = remove_duplicate_geometries(geometries_data)

    energies = extract_energies(input_energy_file)
    categorized_data = categorize_geometries(unique_geometries, energies)
    
    write_output(categorized_data, output_file_path)
    print(f"Output written to {output_file_path}")

