import sys
import parmed as pmd

# Check if the user provided an input filename
if len(sys.argv) != 2:
    print("Usage: python script.py <input_filename_without_extension>")
    sys.exit(1)

# Get the input filename from the command-line arguments
input_name = sys.argv[1]

# Construct input and output filenames
input_file = f"{input_name}.pdb"
output_file = f"{input_name}.rst7"

# Load the structure and save it in the desired format
structure = pmd.load_file(input_file)
structure.save(output_file, format='rst7')
