import numpy as np
import mdtraj as md

# Read atom indices from atom_list.dat
with open('atom_list.dat', 'r') as f:
    indices = np.array(list(map(int, f.readline().split()))) - 1  # zero-indexed

# Load the XYZ trajectory
traj = md.load('interpolated2.xyz', top="dasa-carved_ID.prmtop")

# Load the template .rst7 file
with open('neb_s.rst', 'r') as f:
    template_lines = f.readlines()

# Function to replace coordinates
def replace_coordinates(template, new_coords, frame_indices):
    updated_lines = template[:2]  # Keep header lines as is
    coord_idx = 0
    atom_count = 0
    for i, line in enumerate(template[2:], start=2):
        line_data = line.strip().split()
        new_line_data = []
        # There are 3 coordinates per atom, and each line contains coordinates for 2 atoms
        for j in range(2):  # Two atoms per line
            if atom_count in frame_indices:
                # Replace coordinates for this atom
                new_coords_subset = new_coords[coord_idx]
                coord_idx += 1
                new_line_data.extend(["{:12.7f}".format(coord*10.0) for coord in new_coords_subset])
            else:
                # Keep existing coordinates and ensure they maintain formatting
                new_line_data.extend(["{:12.7f}".format(float(coord)) for coord in line_data[j*3:(j+1)*3]])
            atom_count += 1
        # Ensure spaces between all coordinates
        updated_lines.append("".join(new_line_data) + '\n')
    return updated_lines

# Process each frame
for i, frame in enumerate(traj):
    # Extract coordinates for selected indices
    selected_coords = frame.xyz[0, indices, :]  # Adjusting for the mdtraj xyz shape
    # Replace coordinates in template
    new_rst_content = replace_coordinates(template_lines, selected_coords, indices)
    # Write to a new .rst7 file
    with open(f'neb{i+1}.rst', 'w') as f:
        f.writelines(new_rst_content)

print("All files have been written successfully.")

