import sys
import parmed as pmd

def main():
    if len(sys.argv) != 2:
        print("Usage: python convert_pdb_to_rst7.py <input_file.pdb>")
        sys.exit(1)

    # Get input file from command-line arguments
    input_file = sys.argv[1]

    # Check if the input file has the correct extension
    if not input_file.endswith('.pdb'):
        print("Error: The input file must have a .pdb extension.")
        sys.exit(1)

    # Generate the output file name by replacing .pdb with .rst7
    output_file = input_file.replace('.pdb', '.rst7')

    # Load the structure and save in the new format
    structure = pmd.load_file(input_file)
    structure.save(output_file, format='rst7')

    print(f"RST7 file saved to: {output_file}")

if __name__ == "__main__":
    main()

