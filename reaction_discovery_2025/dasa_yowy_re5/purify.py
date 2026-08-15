def remove_lines(input_file, output_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()

    with open(output_file, 'w') as file:
        for line in lines:
            if "xi" not in line and "transfer" not in line:
                file.write(line)

# Usage
input_file = 'terachem.o13411698'
output_file = 'output.txt'
remove_lines(input_file, output_file)

