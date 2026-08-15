import numpy as np
from scipy.spatial.transform import Rotation as R

def read_pun_file(filepath):
    """
    Read the coordinates from a PUN file.
    """
    with open(filepath, 'r') as file:
        lines = file.readlines()[4:]  # Skip the first 4 lines
        coords = []
        for line in lines:
            parts = line.split()
            if len(parts) >= 4:
                try:
                    x, y, z = map(float, parts[1:4])
                    coords.append((parts[0], np.array([x, y, z])))
                except ValueError:
                    continue  # Skip lines that don't contain valid floats
    return coords

def calculate_angle(coords, indices):
    """
    Calculate the angle for the given atom indices.
    """
    p1, p2, p3 = [coords[i][1] for i in indices]

    ba = p1 - p2
    bc = p3 - p2

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(np.clip(cosine_angle, -1.0, 1.0))
    angle_degrees = np.degrees(angle)

    return angle_degrees

def calculate_dihedral(coords, indices):
    """
    Calculate the dihedral angle for the given atom indices.
    """
    p1, p2, p3, p4 = [coords[i][1] for i in indices]
    
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3
    
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    
    n1 /= np.linalg.norm(n1)
    n2 /= np.linalg.norm(n2)
    
    cos_theta = np.dot(n1, n2)
    m1 = np.cross(n1, n2)
    sin_theta = np.dot(m1, b2/np.linalg.norm(b2))
    
    angle = np.arctan2(sin_theta, cos_theta)
    angle_degrees = np.degrees(angle)
    
    return angle_degrees

def calculate_all_dihedrals(coords):
    """
    Calculate all five dihedrals from the atom indices.
    """
    dihedral_mappings = {
        2: [13, 15, 16, 18],
        4: [15, 16, 18, 20],
        6: [16, 18, 20, 22],
        8: [18, 20, 22, 24]
    }
    
    dihedrals = {}
    for key, indices in dihedral_mappings.items():
        dihedrals[key] = calculate_dihedral(coords, [i  for i in indices])  # Adjust indices for zero-based indexing
    return dihedrals

def calculate_all_angles(coords):
    angle_mappings = {
        1: [15, 16, 18],
        3: [16, 18, 20],
        5: [18, 20, 22],
        7: [20, 22, 24]
    }

    angles = {}
    for key, indices in angle_mappings.items():
        angles[key] = calculate_angle(coords, [i for i in indices])
    return angles

def process_input_file(input_filepath, coords, output_script, xt):
    """
    Process the input file to determine changes in dihedrals, and generate the ChemShell script.
    """
    # Calculate the initial dihedrals from the coordinates
    initial_dihedrals = calculate_all_dihedrals(coords)
    initial_angles = calculate_all_angles(coords)
    with open(input_filepath, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line
        numbers = [float(num) for line in lines for num in line.split()]
        
        final_dihedrals = initial_dihedrals.copy()
        final_angles = initial_angles.copy()
        active_atoms_mapping = {}

        dihedral_mappings = {
            2: [13, 15, 16, 18],
            4: [15, 16, 18, 20],
            6: [16, 18, 20, 22],
            8: [18, 20, 22, 24]
        }
        angle_mappings = {
            1: [15, 16, 18],
            3: [16, 18, 20],
            5: [18, 20, 22],
            7: [20, 22, 24]
        }

        for i, value in enumerate(numbers, start=1):
            if value != 0.0:
                if i in initial_dihedrals:
                    final_dihedrals[i] += value  # Add the input value
                    active_atoms = [index + 1 for index in dihedral_mappings[i]]
                    active_atoms_mapping[i] = active_atoms
                elif i in initial_angles:
                    final_angles[i] += value  # Add the input value
                    active_atoms = [index + 1 for index in angle_mappings[i]]
                    active_atoms_mapping[i] = active_atoms                
                # Adjust active atom indices for 1-based indexing and store them
#                active_atoms = [index + 1 for index in dihedral_mappings[i]]  
#                active_atoms_mapping[i] = active_atoms
        
    # Merge all active atoms into one unique list for dl-find
    #all_active_atoms = sorted(set([atom for atoms in active_atoms_mapping.values() for atom in atoms]))
    all_active_atoms = list(range(12, 43))
    # Now pass 'xt' to generate_chemshell_script
    generate_chemshell_script(coords, initial_dihedrals, final_dihedrals, initial_angles,final_angles, output_script, all_active_atoms, dihedral_mappings,angle_mappings, xt)


def generate_chemshell_script(coords, initial_dihedrals, final_dihedrals, initial_angles,final_angles, output_script, all_active_atoms, dihedral_mappings, angle_mappings, xt):
    """
    Generate a ChemShell script with a stepwise dihedral transformation over 'xt' steps,
    following the correct formula for incremental dihedral changes.
    """
    active_atoms_str = " ".join(map(str, all_active_atoms))  # Join all active atoms into a string

    # Initialize the ChemShell script header (unchanged from the original format)
    chemshell_script = f"""
# Generated ChemShell Script

set residues [ pdb_to_res "dasa-carved_ID.pdb" names={{PT1}} ]

energy coords=pregeom.pun theory=dl_poly : [ list \\
                       amber_prmtop_file= dasa-carved_ID.prmtop \\
                       scale14 = [ list [ expr 1 / 1.2 ] 0.5  ] \\
                       use_pairlist = no \\
                       save_dl_poly_files = yes \\
                       list_option=none ]

set atom_charges [ list_amber_atom_charges ]

source qm_tb.txt
source active_tb.txt

set charge 0
set mult 1

source orca-chemsh.tcl
set orcasimpleinput " ! GFN2-xTB TIGHTSCF"
set orcablocks "
%maxcore 2000
%scf MaxIter 500
end
%pal
nprocs 1
end
"
"""

    # Loop over 'xt' steps for the dihedral adjustment
    for step in range(1, xt + 1):
        # Use "pregeom.pun" for the first step, and "optimized_<step-1>.pun" for subsequent steps
        coords_input = "pregeom.pun" if step == 1 else f"optimized_{step-1}.pun"

        chemshell_script += f"""
# Step {step}: Applying 1/{xt}th of the total dihedral change

dl-find coords={coords_input} active_atoms= [ list {active_atoms_str} ] \\
    theory=restraint : [ list \\
    restraints= [ list \\
        [ list angle {angle_mappings[1][0]+1} {angle_mappings[1][1]+1} {angle_mappings[1][2]+1} {initial_angles[1] + (final_angles[1] - initial_angles[1]) * step / xt:.2f} 1.0 ] \\
        [ list angle {angle_mappings[3][0]+1} {angle_mappings[3][1]+1} {angle_mappings[3][2]+1} {initial_angles[3] + (final_angles[3] - initial_angles[3]) * step / xt:.2f} 1.0 ] \\
        [ list angle {angle_mappings[5][0]+1} {angle_mappings[5][1]+1} {angle_mappings[5][2]+1} {initial_angles[5] + (final_angles[5] - initial_angles[5]) * step / xt:.2f} 1.0 ] \\
        [ list angle {angle_mappings[7][0]+1} {angle_mappings[7][1]+1} {angle_mappings[7][2]+1} {initial_angles[7] + (final_angles[7] - initial_angles[7]) * step / xt:.2f} 1.0 ] \\
        [ list torsion {dihedral_mappings[2][0]+1} {dihedral_mappings[2][1]+1} {dihedral_mappings[2][2]+1} {dihedral_mappings[2][3]+1} {initial_dihedrals[2] + (final_dihedrals[2] - initial_dihedrals[2]) * step / xt:.2f} 1.0 ] \\
        [ list torsion {dihedral_mappings[4][0]+1} {dihedral_mappings[4][1]+1} {dihedral_mappings[4][2]+1} {dihedral_mappings[4][3]+1} {initial_dihedrals[4] + (final_dihedrals[4] - initial_dihedrals[4]) * step / xt:.2f} 1.0 ] \\
        [ list torsion {dihedral_mappings[6][0]+1} {dihedral_mappings[6][1]+1} {dihedral_mappings[6][2]+1} {dihedral_mappings[6][3]+1} {initial_dihedrals[6] + (final_dihedrals[6] - initial_dihedrals[6]) * step / xt:.2f} 1.0 ] \\
        [ list torsion {dihedral_mappings[8][0]+1} {dihedral_mappings[8][1]+1} {dihedral_mappings[8][2]+1} {dihedral_mappings[8][3]+1} {initial_dihedrals[8] + (final_dihedrals[8] - initial_dihedrals[8]) * step / xt:.2f} 1.0 ] ] ] \\
    result=adjusted_{step}.pun

dl-find coords=adjusted_{step}.pun \\
       theory=hybrid : [ list \\
         coupling=shift \\
         debug=no \\
         atom_charges= $atom_charges \\
         qm_region= $qmatoms \\
         qm_theory=orca : [ list \\
              executable= orca \\
              charge= $charge \\
              mult= $mult \\
              orcasimpleinput= $orcasimpleinput \\
              orcablocks= $orcablocks ] \\
         mm_theory=dl_poly : [ list \\
              list_option= none \\
              debug=no \\
              mxexcl=500 \\
              exact_srf= yes \\
              use_pairlist= no \\
              cutoff= 1000 \\
              conn= adjusted_{step}.pun \\
              scale14= [ list [ expr 1 / 1.2 ] 0.5  ] \\
              save_dl_poly_files=yes \\
              amber_prmtop_file= dasa-carved_ID.prmtop ] ] \\
       active_atoms= $active_atoms \\
       coordinates= hdlc \\
       residues= $residues \\
       list_option=full \\
       dump=10 \\
       tolerance=0.009 \\
       maxcycle=10000 \\
       maxstep=0.2 \\
       constraints= {{ {{ angle {angle_mappings[1][0]+1} {angle_mappings[1][1]+1} {angle_mappings[1][2]+1} }} \\
                      {{ angle {angle_mappings[3][0]+1} {angle_mappings[3][1]+1} {angle_mappings[3][2]+1} }} \\
                      {{ angle {angle_mappings[5][0]+1} {angle_mappings[5][1]+1} {angle_mappings[5][2]+1} }} \\
                      {{ angle {angle_mappings[7][0]+1} {angle_mappings[7][1]+1} {angle_mappings[7][2]+1} }} \\
                      {{ torsion {dihedral_mappings[2][0]+1} {dihedral_mappings[2][1]+1} {dihedral_mappings[2][2]+1} {dihedral_mappings[2][3]+1} }} \\
                      {{ torsion {dihedral_mappings[4][0]+1} {dihedral_mappings[4][1]+1} {dihedral_mappings[4][2]+1} {dihedral_mappings[4][3]+1} }} \\
                      {{ torsion {dihedral_mappings[6][0]+1} {dihedral_mappings[6][1]+1} {dihedral_mappings[6][2]+1} {dihedral_mappings[6][3]+1} }} \\
                      {{ torsion {dihedral_mappings[8][0]+1} {dihedral_mappings[8][1]+1} {dihedral_mappings[8][2]+1} {dihedral_mappings[8][3]+1} }} }} \\
       result=optimized_{step}.pun
"""


    # Write the final structure after all 'xt' steps
    chemshell_script += f"""
# Write the final structure after all {xt} steps
#write_xyz file=final_structure.xyz coords=optimized_{xt}.pun
#write_pun file=optimized.pun coords=optimized_{xt}.pun
"""

    # Output the ChemShell script
    with open(output_script, 'w') as file:
        file.write(chemshell_script)

    print(f"ChemShell script written to {output_script}")


# Example usage
pun_file_path = 'pregeom.pun'
input_file_path = 'coord1.txt'
output_script = 'move.chm'
xt = 3  # Number of steps

# Read coordinates
coords = read_pun_file(pun_file_path)

# Process the input file and generate the ChemShell script with 'xt' dihedral steps
process_input_file(input_file_path, coords, output_script, xt)

