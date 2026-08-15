import numpy as np
import mdtraj as md
from mdtraj.core import element, topology

def load_xyz_with_topology(xyz_filename):
    """Load (multi-frame) XYZ and build a minimal mdtraj.Topology."""
    frames_coords = []
    frames_elements = []
    
    with open(xyz_filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    
    i = 0
    while i < len(lines):
        # Attempt to read number of atoms
        try:
            n_atoms = int(lines[i])
        except ValueError:
            # If line isn't an integer (e.g., 'Frame 0'), skip it
            i += 1
            continue
        i += 1  # move past the count line
        
        # The next line is the comment line (which we can ignore or keep)
        # but we do have to skip it
        if i < len(lines):
            comment_line = lines[i]
        i += 1
        
        # Now read n_atoms lines of <Symbol> x y z
        coords = []
        elems  = []
        for _ in range(n_atoms):
            parts = lines[i].split()
            i += 1
            elems.append(parts[0])
            x, y, z = map(float, parts[1:4])
            coords.append([x, y, z])
        
        frames_coords.append(coords)
        frames_elements.append(elems)
    
    # Convert frames_coords to a numpy array of shape (n_frames, n_atoms, 3)
    coords_np = np.array(frames_coords, dtype=np.float32)
    n_frames, n_atoms, _ = coords_np.shape
    
    # Build a minimal topology from the first frame's elements
    top = topology.Topology()
    chain = top.add_chain()
    res = top.add_residue('UNK', chain)
    for symb in frames_elements[0]:
        # Convert element symbol (e.g., 'C') to an mdtraj Element
        # If symbol is not recognized, default to element.carbon or so
        try:
            el = element.get_by_symbol(symb)
        except KeyError:
            el = element.carbon
        top.add_atom(symb, el, res)
    
    # Create mdtraj.Trajectory
    traj = md.Trajectory(xyz=coords_np, topology=top)
    return traj


# ------------------------------------------------------------------------------
# 1) Load the reference geometry using the helper function
# ------------------------------------------------------------------------------
ref = load_xyz_with_topology('32_dft.xyz')
ref.xyz[:] *= 0.1  # Actually convert from Å to nm
print(f"Loaded reference: {ref}")

# ------------------------------------------------------------------------------
# 2) Load the multi-frame trajectory
# ------------------------------------------------------------------------------
traj = load_xyz_with_topology('interpolated_32-41.xyz')
traj.xyz[:] *= 0.1  # Actually convert from Å to nm
print(f"Loaded trajectory with {traj.n_frames} frames, {traj.n_atoms} atoms each.")

# ------------------------------------------------------------------------------
# 3) Read the list of inactive atom indices
#    NOTE: Many chemistry codes use 1-based indexing. mdtraj is 0-based.
# ------------------------------------------------------------------------------
inactive_atom_indices = []
with open('inactive_atoms.txt', 'r') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        # Subtract 1 if your file is 1-based:
        # If your file is already 0-based, remove the -1
        idx = int(line) - 1
        inactive_atom_indices.append(idx)

##print(f"Inactive atom indices: {inactive_atom_indices}")

# ------------------------------------------------------------------------------
# 4) Superpose (align) the trajectory to the reference on the inactive atoms
# ------------------------------------------------------------------------------
traj.superpose(
    reference=ref,
    frame=0,  # align to reference frame 0
    atom_indices=inactive_atom_indices,
    ref_atom_indices=inactive_atom_indices
)

# ------------------------------------------------------------------------------
# 5) Save the aligned trajectory to an XYZ file
# ------------------------------------------------------------------------------
traj.save_xyz('aligned_interpolated_32-41.xyz')
print("Alignment complete.")
## Output: aligned_interpolated_29-51.xyz")

