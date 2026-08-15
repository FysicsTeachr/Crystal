def sort_geometries_by_energy(file_path, block_size=81):
    """
    Sorts geometries in an XYZ file by ascending energy.

    Parameters:
        file_path (str): Path to the XYZ file.
        block_size (int): Number of lines per block (default is 81).

    Returns:
        list: Sorted geometries as (energy, block) tuples.
    """
    sorted_geometries = []

    with open(file_path, 'r') as file:
        while True:
            block = [file.readline() for _ in range(block_size)]
            if not block[0].strip():
                break  # End of file

            try:
                # Extract energy from the second line of the block
                energy_line = block[1].strip()
                energy = float(energy_line.split("Energy")[-1].strip())
            except (IndexError, ValueError) as e:
                print(f"Error parsing energy in block: {block[:2]} - {e}")
                continue

            # Append the block and its energy to the list
            sorted_geometries.append((energy, block))

    # Sort geometries by energy
    sorted_geometries.sort(key=lambda x: x[0])
    return sorted_geometries

# Input file path
file_path = 'all_2.xyz'

# Process the file to sort geometries
sorted_geometries = sort_geometries_by_energy(file_path)

# Save the sorted geometries to a new file
output_file_path = 'sorted_all_2.xyz'
with open(output_file_path, 'w') as sorted_file:
    for _, block in sorted_geometries:
        sorted_file.writelines(block)

print(f"Sorted geometries saved to: {output_file_path}")

