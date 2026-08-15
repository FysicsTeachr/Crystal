import mdtraj as md
import numpy as np

# 1. Load the PDB directly (no top= needed because the PDB already has topology)
traj = md.load('sorted_all_bili_min.pdb')
n_frames = traj.n_frames
print(f"Loaded {n_frames} frames from all_bili_min.pdb")

# 2. (Optional) compute hydrogen bonds using Baker-Hubbard
hbonds_all = md.baker_hubbard(traj)

# 3. Bucket the h-bonds by frame (frame, donor_idx, acceptor_idx)
hbonds_by_frame = [set() for _ in range(n_frames)]
for (frame_i, d_idx, a_idx) in hbonds_all:
    hbonds_by_frame[frame_i].add((d_idx, a_idx))

# 4. Compute RMSD for each pair of consecutive frames with alignment
indices_to_keep = [0]  # always keep the first frame
for i in range(n_frames - 1):
    # Align frame i+1 to frame i before RMSD calculation
    rmsd_val = md.rmsd(traj[i+1], traj[i], atom_indices=None)  # Use all atoms for alignment
    rmsd_val_nm = rmsd_val[0]  # Extract the single value (in nm)

    # Compare RMSD and H-bonds
    if (rmsd_val_nm > 0.05) or (hbonds_by_frame[i] != hbonds_by_frame[i+1]):
        # if RMSD > 0.05 nm (0.5 Å) or H-bonds differ, keep frame (i+1)
        indices_to_keep.append(i+1)

indices_to_keep = sorted(set(indices_to_keep))
print(f"Selected {len(indices_to_keep)} frames out of {n_frames} total.")

# 5. Save these frames to an XYZ or PDB
subset = traj.slice(indices_to_keep)
subset.save_xyz('h-rmsd.xyz')
# or: subset.save_pdb('unique.pdb')

print("Wrote unique.xyz with the filtered frames.") 
