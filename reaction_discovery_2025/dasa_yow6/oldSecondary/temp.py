import sys
import shutil

def append_to_file(temp_file_path, allfiles_path):
    with open(temp_file_path, 'r') as temp_file:
        with open(allfiles_path, 'a') as allfiles:
            shutil.copyfileobj(temp_file, allfiles)

if __name__ == "__main__":
    temp_file_path = sys.argv[1]
    allfiles_path = "../allfiles.txt"
    append_to_file(temp_file_path, allfiles_path)

