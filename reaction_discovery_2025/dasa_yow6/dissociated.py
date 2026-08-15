import os

# Define the directory containing the text files
directory = 'dissociated'

# Define the output file name
output_file = 'dissociated.xyz'

# Create or open the output file in write mode
with open(output_file, 'w') as outfile:
    # Iterate over all files in the directory
    for filename in os.listdir(directory):
        # Check if the file is a text file
        if filename.endswith('.xyz'):
            # Get the full path of the text file
            file_path = os.path.join(directory, filename)
            # Open the text file in read mode
            with open(file_path, 'r') as infile:
                # Read the content of the text file
                content = infile.read()
                # Write the content to the output file
                outfile.write(content)
                # Add a newline character to separate files
#                outfile.write('\n')

print(f"All text files in '{directory}' have been combined into '{output_file}'")

