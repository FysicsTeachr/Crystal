# Define the atom order as per the given example.
atom_order = [
    "C", "H", "H", "H", "C", "H", "H", "H", "C", "O", "O", "O", "O", "C", "C", "C", "C", "H", "C", "O",
    "C", "H", "C", "H", "C", "H", "N", "C", "H", "H", "C", "H", "H", "H", "C", "H", "H", "C", "H", "H",
    "H", "H"
]

def reorder_atoms(block, atom_order):
    reordered_block = [block[0], block[1]]  # The first two lines are kept as is
    for atom, line in zip(atom_order, block[2:]):
        reordered_block.append(f"{atom} {line.strip()}\n")
    return reordered_block

def process_xyz_file(input_file, output_file, atom_order):
    with open(input_file, 'r') as file:
        lines = file.readlines()
    
    blocks = []
    current_block = []

    for line in lines:
        if len(current_block) == 0 or len(current_block) < 44:
            current_block.append(line)
        else:
            blocks.append(current_block)
            current_block = [line]

    if current_block:
        blocks.append(current_block)

    reordered_blocks = []
    for block in blocks:
        reordered_block = reorder_atoms(block, atom_order)
        reordered_blocks.append(reordered_block)

    # Write the reordered content to the output file
    with open(output_file, 'w') as file:
        for block in reordered_blocks:
            file.writelines(block)

# Define input and output filenames
input_filename = 'climbing.xyz'  # Replace with your actual input file path
output_filename = 'climbings.xyz'  # Replace with your desired output file path

# Process the input file and generate the output file
process_xyz_file(input_filename, output_filename, atom_order)

