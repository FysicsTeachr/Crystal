import numpy as np

def calculate_distance(coord1, coord2):
    return np.sqrt((coord1['x'] - coord2['x'])**2 +
                   (coord1['y'] - coord2['y'])**2 +
                   (coord1['z'] - coord2['z'])**2)

def get_coordinates(line):
    parts = line.split()
    return {
        'x': float(parts[1]),
        'y': float(parts[2]),
        'z': float(parts[3])
    }

def main(file_path):
    pairs_to_check = [
        (24, 26),
        (22, 24),
        (22, 20),
        (20, 18),
        (18, 16),
        (16, 15),
        (15, 14),
        (15, 13),
        (18, 19),
    ]

    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()

        condition_met = False

        # First check the original criteria
        for pair in pairs_to_check:
            index1, index2 = pair[0] + 2, pair[1] + 2
            coord1 = get_coordinates(lines[index1])
            coord2 = get_coordinates(lines[index2])
            distance = calculate_distance(coord1, coord2)
            if distance > 1.7:
                condition_met = True
                break

        # Check the new criteria
        if not condition_met:
            coord_reference = get_coordinates(lines[41 + 2])
            all_far = True
            for i in range(0, 41):
                coord = get_coordinates(lines[i + 2])
                distance = calculate_distance(coord, coord_reference)
                
                # Check if any distance is less than 0.5
                if distance < 0.8:
                    condition_met = True
                    break
                
                if distance <= 1.3:
                    all_far = False

            if all_far:
                condition_met = True

        result = "yes" if condition_met else "no"

        with open('dissociated.txt', 'w') as file:
            file.write(result)
            
    except Exception as e:
        print(f"An error occurred: {e}")

# Running the function with the provided file path
main('tgeom.xyz')

