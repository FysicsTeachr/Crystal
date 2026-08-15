import numpy as np

###############################################################################
# 1) Parsing the multi-geometry XYZ file
###############################################################################
def parse_xyz_geometries(filename):
    """
    Reads all geometry blocks from a multi-geometry XYZ file.
    Returns a list [(elements, coords), (elements, coords), ...]
       elements: list of atomic symbols, e.g. ['O','O','N','H',...]
       coords  : np.array of shape (n_atoms, 3)
    """
    geometries = []
    with open(filename, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break  # EOF
            line = line.strip()
            # First line should be number_of_atoms
            if not line.isdigit():
                break

            n_atoms = int(line)
            # Next line is a comment line, e.g. "Iteration 62 Energy ..."
            _ = f.readline()

            elems = []
            coords_list = []
            for _ in range(n_atoms):
                row = f.readline()
                if not row:
                    break
                parts = row.split()
                elem = parts[0]
                x, y, z = map(float, parts[1:4])
                elems.append(elem)
                coords_list.append([x, y, z])

            if len(elems) < n_atoms:
                # Possibly hit EOF mid-block
                break

            coords_array = np.array(coords_list, dtype=np.float64)
            geometries.append((elems, coords_array))

    return geometries

###############################################################################
# 2) Basic geometry routines: distance, angle, dihedral
###############################################################################
def distance(p1, p2):
    """Return the Euclidean distance between p1 and p2 (each length-3)."""
    return np.sqrt(np.sum((p2 - p1)**2))

def angle(heavy_coords, h_coords, acceptor_coords):
    """
    Compute the angle (in degrees) at the hydrogen:
    (heavy_atom --- H --- acceptor).
    """
    v1 = heavy_coords - h_coords
    v2 = acceptor_coords - h_coords
    dot_val = np.dot(v1, v2)
    mag1 = np.linalg.norm(v1)
    mag2 = np.linalg.norm(v2)
    cos_theta = dot_val / (mag1 * mag2 + 1e-14)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)  # numerical safety
    theta_rad = np.arccos(cos_theta)
    return np.degrees(theta_rad)

def compute_dihedral(coords, i0, i1, i2, i3):
    """
    Compute the dihedral angle (in degrees) for the atoms (i0,i1,i2,i3),
    using the standard cross-product method.
    Returns an angle in the range (-180, 180].
    """
    p0 = coords[i0]
    p1 = coords[i1]
    p2 = coords[i2]
    p3 = coords[i3]

    # Define the three bond vectors
    b1 = p1 - p0
    b2 = p2 - p1
    b3 = p3 - p2

    # Normal vectors to the planes
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    # Normalize them to unit vectors
    n1_mag = np.linalg.norm(n1) + 1e-14
    n2_mag = np.linalg.norm(n2) + 1e-14
    n1n = n1 / n1_mag
    n2n = n2 / n2_mag

    # cos(phi) = dot(n1n, n2n)
    cos_phi = np.dot(n1n, n2n)
    cos_phi = np.clip(cos_phi, -1.0, 1.0)  # clip just in case
    phi = np.arccos(cos_phi)

    # To get sign, look at cross(n1, n2) dot b2
    sign = np.sign(np.dot(np.cross(n1n, n2n), b2))
    angle_deg = np.degrees(phi) * sign

    # Ensure it's in (-180, 180]
    if angle_deg > 180.0:
        angle_deg -= 360.0
    elif angle_deg <= -180.0:
        angle_deg += 360.0

    return angle_deg

###############################################################################
# 3) Find donors: O/N--H pairs. All O/N are acceptors.
###############################################################################
def find_donors(elements, coords, oh_cutoff=1.1):
    """
    Identify donor pairs (heavy_atom_idx, H_idx) by checking for
    O/N---H distances <= oh_cutoff (1.1 Å).
    """
    donors = []
    # All O/N indices:
    on_indices = [i for i, e in enumerate(elements) if e in ('O', 'N')]
    for h_idx, elem_h in enumerate(elements):
        if elem_h == 'H':
            hpos = coords[h_idx]
            for on_idx in on_indices:
                onpos = coords[on_idx]
                if distance(hpos, onpos) <= oh_cutoff:
                    donors.append((on_idx, h_idx))
    return donors

###############################################################################
# 4) Main
###############################################################################
def main():
#    xyz_file = "2_unique_2-sort.xyz"
#    out_file = "2_unique_2-info.txt"

    xyz_file = "sorted_all_2.xyz"
    out_file = "2_unique_1-list.txt"

    # H-bond geometry criteria
    max_h_acceptor_dist = 2.5   # Å
    min_angle_deg       = 120.0 # degrees

    # Dihedrals of interest (zero-based indices)
    dihedral_indices = [
        [49,46,34,10],
        [10,31,29,23],
        [31,29,23, 8],
        [34,46,49,12],
        [12,52,64,66],
        [52,64,66,14],
    ]

    # Parse the entire file in memory
    all_geometries = parse_xyz_geometries(xyz_file)

    # We build our output lines in memory, then write once for speed
    output_lines = []
    output_lines.append(f"Found {len(all_geometries)} geometries in {xyz_file}\n")

    for g_idx, (elements, coords) in enumerate(all_geometries):
        # 1) Donors & Acceptors
        donors = find_donors(elements, coords, oh_cutoff=1.1)
        acceptors = [i for i, e in enumerate(elements) if e in ('O','N')]

        # 2) Check each donor vs. each acceptor for geometry
        hbonds = []  # list of tuples (heavy_i, H_i, acc_i, distance, angle)
        for (heavy_i, h_i) in donors:
            heavy_pos = coords[heavy_i]
            h_pos = coords[h_i]
            for acc_i in acceptors:
                if acc_i == heavy_i or acc_i == h_i:
                    continue
                a_pos = coords[acc_i]

                dist_ha = distance(h_pos, a_pos)
                if dist_ha <= max_h_acceptor_dist:
                    ang_deg = angle(heavy_pos, h_pos, a_pos)
                    if ang_deg >= min_angle_deg:
                        hbonds.append((heavy_i, h_i, acc_i, dist_ha, ang_deg))

        # 3) Compute the 6 dihedral angles
        dihedral_vals = []
        for (i0, i1, i2, i3) in dihedral_indices:
            # In case any index is out of range, handle gracefully
            if max(i0, i1, i2, i3) < len(elements):
                dval = compute_dihedral(coords, i0, i1, i2, i3)
            else:
                dval = 0.0  # or np.nan, but hopefully all are valid
            dihedral_vals.append(dval)

        # 4) Prepare output
        output_lines.append(f"Geometry {g_idx}:")
        output_lines.append(f"  # of hydrogen bonds = {len(hbonds)}")
        dihedrals_str = " ".join(f"{v:.2f}" for v in dihedral_vals)
        output_lines.append(f"  Dihedrals (degrees): {dihedrals_str}")

        for (hv_i, h_i, ac_i, dist_val, ang_val) in hbonds:
            output_lines.append(
                f"    Donor=(heavy={hv_i}, H={h_i}), Acceptor={ac_i}, "
                f"dist={dist_val:.2f} , angle={ang_val:.1f}"                
#                f"    Donor=(heavy={hv_i}, H={h_i}), Acceptor={ac_i}, "
#                f"dist={dist_val:.2f} Å, angle={ang_val:.1f}°"
            )
        output_lines.append("")  # blank line after each geometry

    # Finally, write the output
    with open(out_file, 'w') as f:
        f.write("\n".join(output_lines))

    print(f"Done. Wrote results to '{out_file}'")

if __name__ == "__main__":
    main()

