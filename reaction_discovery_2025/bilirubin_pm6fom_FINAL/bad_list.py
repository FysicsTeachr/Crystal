import os
import glob

def check_last_line(base_path):
    count = 0
    files_without_phrase = []
    
    # Loop through potential directories by scanning the base path for any 'secondaryx' pattern
    for dir_name in os.listdir(base_path):
        if dir_name.startswith("secondary") and os.path.isdir(os.path.join(base_path, dir_name)):
            directory = os.path.join(base_path, dir_name, "test/allslurms")
            if os.path.exists(directory):
                # List all .log files in the directory
                log_files = glob.glob(os.path.join(directory, '*.log'))
                for file_path in log_files:
                    count += 1
                    try:
                        with open(file_path, 'r') as file:
                            lines = file.readlines()
                            if "Time elapsed" not in lines[-1].strip():
                                files_without_phrase.append(file_path)
                    except Exception as e:
                        print(f"Error reading file {file_path}: {e}")

    # Output the results
    print(f"Total number of .log files checked: {count}")
    if files_without_phrase:
        print("Files where the last line does not contain 'Time elapsed':")
        for file in files_without_phrase:
            print(file)
    else:
        print("All files contain 'Time elapsed' in the last line.")

# Specify the base path according to the directory you're working from
base_path = "/lustre/work/pan60047/bili-geometric-pm7"
check_last_line(base_path)

