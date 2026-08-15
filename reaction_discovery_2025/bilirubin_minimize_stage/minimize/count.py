import os

# Base directory containing the folders
base_dir = "folders"
gradient_log_count = 0

# Iterate through all subdirectories in the base directory
for root, dirs, files in os.walk(base_dir):
    if "gradient.log" in files:
        gradient_log_path = os.path.join(root, "gradient.log")
        try:
            # Read the last line of the file
            with open(gradient_log_path, "r") as file:
                last_line = file.readlines()[-1].strip()
                # Check if the last line contains "Time elapsed"
                if "Time elapsed" in last_line:
                    gradient_log_count += 1
        except Exception as e:
            print(f"Error reading file {gradient_log_path}: {e}")

print(f"Number of gradient.log files with 'Time elapsed' in the last line: {gradient_log_count}")

